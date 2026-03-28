"""Query helpers for gene-/model-/lineage-level lookups on DepMap data."""

from pathlib import Path
from typing import Literal

import pandas as pd

from .config import get_logger
from .database.connection import DatabaseManager, get_db_manager

logger = get_logger(__name__)

SummaryGrouping = Literal["lineage", "disease"]
QueryFormat = Literal["table", "csv", "json"]
MutationClass = Literal["any", "likely_lof", "hotspot", "driver"]

_GROUPING_COLUMN_MAP: dict[SummaryGrouping, str] = {
    "lineage": "oncotree_lineage",
    "disease": "oncotree_primary_disease",
}

_MATRIX_TABLES = {
    "dependency": ("gene_effects_wide", {"model_id", "created_at"}),
    "expression": (
        "gene_expression_wide",
        {
            "model_id",
            "sequencing_id",
            "model_condition_id",
            "is_default_entry",
            "is_default_for_mc",
            "created_at",
        },
    ),
}


class GeneQueryService:
    """Query service for gene-centric lookups on wide tables."""

    def __init__(self, db_manager: DatabaseManager | None = None) -> None:
        self.db_manager = db_manager or get_db_manager()

    def get_dependency_summary(
        self,
        gene: str,
        *,
        group_by: SummaryGrouping = "lineage",
        limit: int | None = None,
    ) -> pd.DataFrame:
        """Summarize dependency values for a gene by lineage or disease."""
        gene_column = self._resolve_gene_column("dependency", gene)
        grouping_column = _GROUPING_COLUMN_MAP[group_by]
        limit_clause = f"LIMIT {limit}" if limit is not None else ""

        query = f"""
        SELECT
            m.{grouping_column} AS group_name,
            COUNT(gef."{gene_column}") AS model_count,
            AVG(gef."{gene_column}") AS mean_dependency,
            MEDIAN(gef."{gene_column}") AS median_dependency,
            MIN(gef."{gene_column}") AS min_dependency,
            MAX(gef."{gene_column}") AS max_dependency
        FROM gene_effects_wide gef
        INNER JOIN models m ON m.model_id = gef.model_id
        WHERE m.{grouping_column} IS NOT NULL
          AND gef."{gene_column}" IS NOT NULL
        GROUP BY 1
        ORDER BY mean_dependency ASC, model_count DESC, group_name ASC
        {limit_clause}
        """
        return self.db_manager.fetch_df(query)

    def get_expression_summary(
        self,
        gene: str,
        *,
        group_by: SummaryGrouping = "lineage",
        limit: int | None = None,
    ) -> pd.DataFrame:
        """Summarize expression values for a gene by lineage or disease."""
        gene_column = self._resolve_gene_column("expression", gene)
        grouping_column = _GROUPING_COLUMN_MAP[group_by]
        limit_clause = f"LIMIT {limit}" if limit is not None else ""

        query = f"""
        SELECT
            m.{grouping_column} AS group_name,
            COUNT(gex."{gene_column}") AS model_count,
            AVG(gex."{gene_column}") AS mean_expression,
            MEDIAN(gex."{gene_column}") AS median_expression,
            MIN(gex."{gene_column}") AS min_expression,
            MAX(gex."{gene_column}") AS max_expression
        FROM gene_expression_wide gex
        INNER JOIN models m ON m.model_id = gex.model_id
        WHERE m.{grouping_column} IS NOT NULL
          AND gex."{gene_column}" IS NOT NULL
        GROUP BY 1
        ORDER BY mean_expression DESC, model_count DESC, group_name ASC
        {limit_clause}
        """
        return self.db_manager.fetch_df(query)

    def get_dependency_models(
        self,
        gene: str,
        *,
        lineage: str | None = None,
        disease: str | None = None,
        limit: int | None = None,
    ) -> pd.DataFrame:
        """Return per-model dependency values for a gene."""
        gene_column = self._resolve_gene_column("dependency", gene)
        where_clause, parameters = self._build_model_filters(lineage, disease)
        limit_clause = f"LIMIT {limit}" if limit is not None else ""
        non_null_clause = f'gef."{gene_column}" IS NOT NULL'
        if where_clause:
            where_clause = f"{where_clause} AND {non_null_clause}"
        else:
            where_clause = f"WHERE {non_null_clause}"

        query = f"""
        SELECT
            m.model_id,
            m.cell_line_name,
            m.oncotree_lineage,
            m.oncotree_primary_disease,
            m.oncotree_subtype,
            gef."{gene_column}" AS dependency
        FROM gene_effects_wide gef
        INNER JOIN models m ON m.model_id = gef.model_id
        {where_clause}
        ORDER BY dependency ASC, m.cell_line_name ASC
        {limit_clause}
        """
        return self.db_manager.fetch_df(query, parameters)

    def export_dependency_models_csv(
        self,
        gene: str,
        output_path: Path,
        *,
        lineage: str | None = None,
        disease: str | None = None,
        limit: int | None = None,
    ) -> int:
        """Export per-model dependency values for a gene to CSV."""
        df = self.get_dependency_models(
            gene,
            lineage=lineage,
            disease=disease,
            limit=limit,
        )
        output_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_path, index=False)
        return len(df)

    def get_dependency_by_mutation(
        self,
        target_gene: str,
        *,
        mutation_gene: str,
        mutation_class: MutationClass = "any",
        lineage: str | None = None,
        disease: str | None = None,
    ) -> pd.DataFrame:
        """Compare dependency for mutant vs wild-type cohorts."""
        dependency_column = self._resolve_gene_column("dependency", target_gene)
        where_clause, parameters = self._build_model_filters(lineage, disease)
        mutation_condition = self._mutation_class_condition("mgms", mutation_class)
        if where_clause:
            where_clause = f"{where_clause} AND gef.\"{dependency_column}\" IS NOT NULL"
        else:
            where_clause = f'WHERE gef."{dependency_column}" IS NOT NULL'

        query = f"""
        WITH cohort AS (
            SELECT
                m.model_id,
                m.cell_line_name,
                m.oncotree_lineage,
                m.oncotree_primary_disease,
                gef."{dependency_column}" AS dependency,
                CASE
                    WHEN {mutation_condition} THEN TRUE
                    ELSE FALSE
                END AS is_mutant
            FROM gene_effects_wide gef
            INNER JOIN models m ON m.model_id = gef.model_id
            LEFT JOIN model_gene_mutation_status mgms
                ON mgms.model_id = m.model_id
               AND lower(mgms.gene_symbol) = lower(?)
            {where_clause}
        )
        SELECT
            ? AS target_gene,
            ? AS mutation_gene,
            ? AS mutation_class,
            COUNT(*) FILTER (WHERE is_mutant) AS mutant_model_count,
            COUNT(*) FILTER (WHERE NOT is_mutant) AS wt_model_count,
            AVG(dependency) FILTER (WHERE is_mutant) AS mutant_mean_dependency,
            AVG(dependency) FILTER (WHERE NOT is_mutant) AS wt_mean_dependency,
            MEDIAN(dependency) FILTER (WHERE is_mutant) AS mutant_median_dependency,
            MEDIAN(dependency) FILTER (WHERE NOT is_mutant) AS wt_median_dependency,
            AVG(dependency) FILTER (WHERE is_mutant)
                - AVG(dependency) FILTER (WHERE NOT is_mutant) AS delta_mean
        FROM cohort
        """
        query_parameters = [
            mutation_gene,
            *parameters,
            target_gene,
            mutation_gene,
            mutation_class,
        ]
        return self.db_manager.fetch_df(query, query_parameters)

    def _resolve_gene_column(
        self, data_type: Literal["dependency", "expression"], gene: str
    ) -> str:
        """Resolve a gene symbol to the actual matrix column name."""
        table_name, excluded_columns = _MATRIX_TABLES[data_type]
        placeholders = ", ".join("?" for _ in excluded_columns)
        query = f"""
        SELECT column_name
        FROM information_schema.columns
        WHERE table_name = ?
          AND column_name NOT IN ({placeholders})
          AND lower(column_name) = lower(?)
        ORDER BY column_name
        """
        parameters: list[str] = [table_name, *sorted(excluded_columns), gene]
        result = self.db_manager.execute(query, parameters).fetchall()

        if not result:
            raise ValueError(
                f"Gene '{gene}' was not found in {table_name}."
            )
        if len(result) > 1:
            raise ValueError(
                f"Gene '{gene}' matched multiple columns in {table_name}."
            )

        return str(result[0][0])

    def _build_model_filters(
        self, lineage: str | None, disease: str | None
    ) -> tuple[str, list[str]]:
        """Build SQL filter clause for model-level queries."""
        clauses: list[str] = []
        parameters: list[str] = []

        if lineage:
            clauses.append("lower(m.oncotree_lineage) = lower(?)")
            parameters.append(lineage)
        if disease:
            clauses.append("lower(m.oncotree_primary_disease) = lower(?)")
            parameters.append(disease)

        if not clauses:
            return "", parameters

        where_clause = "WHERE " + " AND ".join(clauses)
        return where_clause, parameters

    def _mutation_class_condition(
        self, alias: str, mutation_class: MutationClass
    ) -> str:
        """Map mutation-class names to model_gene_mutation_status predicates."""
        conditions = {
            "any": f"{alias}.is_mutated = TRUE",
            "likely_lof": f"{alias}.has_likely_lof = TRUE",
            "hotspot": f"{alias}.has_hotspot = TRUE",
            "driver": f"{alias}.has_driver = TRUE",
        }
        return conditions[mutation_class]


class MutationQueryService:
    """Query service for mutation event and prevalence lookups."""

    def __init__(self, db_manager: DatabaseManager | None = None) -> None:
        self.db_manager = db_manager or get_db_manager()

    def get_model_mutations(
        self,
        model_query: str,
        *,
        gene: str | None = None,
        lof_only: bool = False,
        hotspot_only: bool = False,
        driver_only: bool = False,
        limit: int | None = None,
    ) -> pd.DataFrame:
        """Return event-level mutations for a resolved model/cell line."""
        model = self._resolve_model(model_query)
        clauses = ["mu.model_id = ?"]
        parameters: list[str] = [str(model["model_id"])]

        if gene:
            clauses.append("lower(COALESCE(mu.hugo_symbol, mu.hgnc_name)) = lower(?)")
            parameters.append(gene)
        if lof_only:
            clauses.append("mu.likely_lof = TRUE")
        if hotspot_only:
            clauses.append("mu.hotspot = TRUE")
        if driver_only:
            clauses.append("mu.hess_driver = TRUE")

        where_clause = "WHERE " + " AND ".join(clauses)
        limit_clause = f"LIMIT {limit}" if limit is not None else ""

        query = f"""
        SELECT
            mu.model_id,
            ? AS matched_cell_line,
            COALESCE(mu.hugo_symbol, mu.hgnc_name) AS gene,
            mu.protein_change,
            mu.dna_change,
            mu.molecular_consequence,
            mu.vep_impact,
            mu.hotspot,
            mu.likely_lof,
            mu.hess_driver AS driver,
            mu.af,
            mu.dp,
            mu.chrom,
            mu.pos,
            mu.ref,
            mu.alt,
            mu.variant_type
        FROM mutations mu
        {where_clause}
        ORDER BY COALESCE(mu.hugo_symbol, mu.hgnc_name) ASC, mu.pos ASC, mu.protein_change ASC
        {limit_clause}
        """
        return self.db_manager.fetch_df(query, [str(model["cell_line_name"]), *parameters])

    def get_lineage_mutation_frequency(
        self, lineage: str, *, limit: int | None = None
    ) -> pd.DataFrame:
        """Rank most frequently mutated genes in a lineage."""
        limit_clause = f"LIMIT {limit}" if limit is not None else ""
        query = f"""
        WITH lineage_models AS (
            SELECT model_id
            FROM models
            WHERE lower(oncotree_lineage) = lower(?)
        ),
        lineage_totals AS (
            SELECT COUNT(*) AS total_models FROM lineage_models
        )
        SELECT
            mgms.gene_symbol,
            COUNT(*) AS mutated_model_count,
            MAX(lt.total_models) AS total_models,
            COUNT(*)::DOUBLE / NULLIF(MAX(lt.total_models), 0) AS mutation_frequency,
            SUM(CASE WHEN mgms.has_hotspot THEN 1 ELSE 0 END) AS hotspot_model_count,
            SUM(CASE WHEN mgms.has_likely_lof THEN 1 ELSE 0 END) AS likely_lof_model_count,
            SUM(CASE WHEN mgms.has_driver THEN 1 ELSE 0 END) AS driver_model_count
        FROM model_gene_mutation_status mgms
        INNER JOIN lineage_models lm ON lm.model_id = mgms.model_id
        CROSS JOIN lineage_totals lt
        GROUP BY mgms.gene_symbol
        ORDER BY mutated_model_count DESC, mutation_frequency DESC, mgms.gene_symbol ASC
        {limit_clause}
        """
        return self.db_manager.fetch_df(query, [lineage])

    def _resolve_model(self, model_query: str) -> pd.Series:
        """Resolve a model by exact id/name first, then unique partial match."""
        exact_query = """
        SELECT
            model_id,
            cell_line_name,
            stripped_cell_line_name,
            ccle_name,
            oncotree_lineage,
            oncotree_primary_disease,
            CASE
                WHEN lower(model_id) = lower(?) THEN 1
                WHEN lower(cell_line_name) = lower(?) THEN 2
                WHEN lower(COALESCE(stripped_cell_line_name, '')) = lower(?) THEN 3
                WHEN lower(COALESCE(ccle_name, '')) = lower(?) THEN 4
                ELSE 100
            END AS match_rank
        FROM models
        WHERE lower(model_id) = lower(?)
           OR lower(cell_line_name) = lower(?)
           OR lower(COALESCE(stripped_cell_line_name, '')) = lower(?)
           OR lower(COALESCE(ccle_name, '')) = lower(?)
        ORDER BY match_rank, model_id
        """
        exact_params = [model_query, model_query, model_query, model_query] * 2
        exact = self.db_manager.fetch_df(exact_query, exact_params)
        if len(exact) == 1:
            return exact.iloc[0]
        if len(exact) > 1:
            choices = ", ".join(
                f"{row.model_id} ({row.cell_line_name})" for row in exact.itertuples()
            )
            raise ValueError(
                f"Model query '{model_query}' is ambiguous. Matches: {choices}"
            )

        partial_query = """
        SELECT
            model_id,
            cell_line_name,
            stripped_cell_line_name,
            ccle_name,
            oncotree_lineage,
            oncotree_primary_disease
        FROM models
        WHERE lower(model_id) LIKE lower(?)
           OR lower(cell_line_name) LIKE lower(?)
           OR lower(COALESCE(stripped_cell_line_name, '')) LIKE lower(?)
           OR lower(COALESCE(ccle_name, '')) LIKE lower(?)
        ORDER BY cell_line_name, model_id
        """
        like_term = f"%{model_query}%"
        partial = self.db_manager.fetch_df(
            partial_query, [like_term, like_term, like_term, like_term]
        )
        if len(partial) == 1:
            return partial.iloc[0]
        if partial.empty:
            raise ValueError(f"Model '{model_query}' was not found.")

        preview = ", ".join(
            f"{row.model_id} ({row.cell_line_name})"
            for row in partial.head(10).itertuples()
        )
        suffix = "" if len(partial) <= 10 else f" ... and {len(partial) - 10} more"
        raise ValueError(
            f"Model query '{model_query}' matched multiple models: {preview}{suffix}"
        )
