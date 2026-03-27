"""Query helpers for gene-level lookups on wide DepMap matrices."""

from pathlib import Path
from typing import Literal

import pandas as pd

from .config import get_logger
from .database.connection import DatabaseManager, get_db_manager

logger = get_logger(__name__)

SummaryGrouping = Literal["lineage", "disease"]
QueryFormat = Literal["table", "csv", "json"]

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
