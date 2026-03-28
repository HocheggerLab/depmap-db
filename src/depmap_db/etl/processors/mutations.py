"""Processor for OmicsSomaticMutations data."""

from __future__ import annotations

import contextlib
from datetime import datetime
from pathlib import Path

import duckdb
import pandas as pd

from .base import BaseProcessor, ProcessingResult


class MutationsProcessor(BaseProcessor):
    """Processor for OmicsSomaticMutations.csv files.

    Loads mutation event data into the canonical mutations table
    and derives a sparse model_gene_mutation_status table for mutated pairs.
    """

    _TEXT_NULL_MARKERS = ("", "nan", "none", "null")
    _BOOLEAN_TRUE_MARKERS = ("true", "1", "yes", "t")

    def __init__(self, batch_size: int | None = None):
        """Initialize mutations processor."""
        super().__init__("OmicsSomaticMutations", batch_size)

    def get_table_name(self) -> str:
        """Get the target database table name."""
        return "mutations"

    def validate_data(
        self, df: pd.DataFrame
    ) -> tuple[pd.DataFrame, list[str]]:
        """Validate mutation data.

        The real mutation import path is DuckDB-native and does not materialize the
        full dataset in pandas. This method remains for tests/smaller inputs and to
        preserve the BaseProcessor contract.
        """
        warnings: list[str] = []

        if "ModelID" not in df.columns:
            raise ValueError("Missing required columns: ['ModelID']")

        missing_model_ids = df["ModelID"].isna().sum()
        if missing_model_ids > 0:
            warnings.append(
                f"Removing {missing_model_ids} rows with missing ModelID"
            )
            df = df.dropna(subset=["ModelID"])

        return df, warnings

    def transform_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """Transform mutation data.

        For large real-world loads we do this in DuckDB SQL for speed. Returning the
        dataframe unchanged keeps this method available for compatibility.
        """
        return df

    def process_file(
        self, file_path: Path, force_reload: bool = False
    ) -> ProcessingResult:
        """Process mutation file with DuckDB bulk staging and atomic reload."""
        start_time = datetime.now()
        warnings: list[str] = []

        try:
            if not file_path.exists():
                raise FileNotFoundError(f"Data file not found: {file_path}")

            existing_count = self._get_existing_record_count()
            if existing_count > 0 and not force_reload:
                self.logger.info(
                    "Table %s already has %s records",
                    self.get_table_name(),
                    existing_count,
                )
                return ProcessingResult(
                    processor_name=self.__class__.__name__,
                    dataset_name=self.dataset_name,
                    source_file=file_path,
                    records_processed=0,
                    records_inserted=0,
                    records_updated=0,
                    records_skipped=existing_count,
                    processing_time_seconds=0.0,
                    status="skipped",
                )

            source_columns = self._read_source_columns(file_path)
            if "ModelID" not in source_columns:
                raise ValueError("Missing required columns: ['ModelID']")

            conn = self.db_manager.connect()
            total_records = self._count_source_rows(file_path)

            self.logger.info(
                "Bulk staging %s mutation rows from %s",
                f"{total_records:,}",
                file_path.name,
            )
            self._create_stage_table(conn, file_path, source_columns)

            stage_counts = conn.execute(
                """
                SELECT
                    COUNT(*) AS staged_rows,
                    COUNT(DISTINCT model_id) AS distinct_models,
                    COUNT(
                        DISTINCT COALESCE(
                            CAST(hugo_symbol AS VARCHAR),
                            CAST(hgnc_name AS VARCHAR)
                        )
                    ) FILTER (
                        WHERE COALESCE(
                            CAST(hugo_symbol AS VARCHAR),
                            CAST(hgnc_name AS VARCHAR)
                        ) IS NOT NULL
                    ) AS distinct_genes
                FROM mutation_stage
                """
            ).fetchone()
            if stage_counts is None:
                raise RuntimeError("Failed to read staged mutation counts")

            staged_rows = int(stage_counts[0])
            distinct_models = int(stage_counts[1])
            distinct_genes = int(stage_counts[2])
            skipped_rows = total_records - staged_rows
            if skipped_rows > 0:
                warnings.append(
                    "Dropped "
                    f"{skipped_rows:,} source rows during mutation staging "
                    "(missing ModelID and/or exact duplicates after normalization)"
                )

            if force_reload:
                warnings.append(
                    "Force reload enabled: replaced both mutations and "
                    "model_gene_mutation_status from staged data"
                )

            self._reload_tables(conn)

            derived_count_row = conn.execute(
                "SELECT COUNT(*) FROM model_gene_mutation_status"
            ).fetchone()
            if derived_count_row is None:
                raise RuntimeError(
                    "Failed to count model_gene_mutation_status rows"
                )
            derived_count = int(derived_count_row[0])

            self._record_import(file_path, staged_rows, "success")

            processing_time = (datetime.now() - start_time).total_seconds()
            warnings.append(
                "Loaded "
                f"{staged_rows:,} unique mutation events across "
                f"{distinct_models:,} models and {distinct_genes:,} genes"
            )
            warnings.append(
                "Derived "
                f"{int(derived_count):,} model-gene mutation status rows"
            )

            return ProcessingResult(
                processor_name=self.__class__.__name__,
                dataset_name=self.dataset_name,
                source_file=file_path,
                records_processed=total_records,
                records_inserted=staged_rows,
                records_updated=0,
                records_skipped=skipped_rows,
                processing_time_seconds=processing_time,
                status="success",
                warnings=warnings,
            )

        except (OSError, ValueError, RuntimeError) as e:
            processing_time = (datetime.now() - start_time).total_seconds()
            error_msg = f"Processing failed: {e}"
            self.logger.error("%s", error_msg)

            with contextlib.suppress(RuntimeError):
                self._record_import(file_path, 0, "failed", error_msg)

            return ProcessingResult(
                processor_name=self.__class__.__name__,
                dataset_name=self.dataset_name,
                source_file=file_path,
                records_processed=0,
                records_inserted=0,
                records_updated=0,
                records_skipped=0,
                processing_time_seconds=processing_time,
                status="failed",
                error_message=error_msg,
            )

    def get_existing_records(self) -> set[str]:
        """Get set of existing mutation IDs."""
        table_name = self.get_table_name()

        try:
            if not self.db_manager.table_exists(table_name):
                return set()

            result = self.db_manager.execute(
                f"SELECT mutation_id FROM {table_name}"
            )
            existing = {row[0] for row in result.fetchall()}

            self.logger.info(
                "Found %s existing mutation records", len(existing)
            )
            return existing

        except RuntimeError as e:
            self.logger.warning("Could not retrieve existing records: %s", e)
            return set()

    def _read_source_columns(self, file_path: Path) -> set[str]:
        """Read just the CSV header to discover source columns."""
        return set(pd.read_csv(file_path, nrows=0).columns)

    def _count_source_rows(self, file_path: Path) -> int:
        """Count source rows using DuckDB without loading the full file into pandas."""
        result = self.db_manager.execute(
            "SELECT COUNT(*) FROM read_csv_auto(?, header=true, sample_size=-1)",
            [str(file_path)],
        )
        row = result.fetchone()
        if row is None:
            raise RuntimeError("Failed to count mutation source rows")
        return int(row[0])

    def _create_stage_table(
        self,
        conn: duckdb.DuckDBPyConnection,
        file_path: Path,
        source_columns: set[str],
    ) -> None:
        """Create a temporary staged mutation table from the CSV file."""
        transformed_columns = [
            ("model_id", self._normalized_text_expr("ModelID", source_columns)),
            ("chrom", self._normalized_text_expr("Chrom", source_columns)),
            ("pos", self._numeric_expr("Pos", source_columns, "BIGINT")),
            ("ref", self._normalized_text_expr("Ref", source_columns)),
            ("alt", self._normalized_text_expr("Alt", source_columns)),
            (
                "variant_type",
                self._normalized_text_expr("VariantType", source_columns),
            ),
            ("dna_change", self._normalized_text_expr("DNAChange", source_columns)),
            (
                "protein_change",
                self._normalized_text_expr("ProteinChange", source_columns),
            ),
            ("hugo_symbol", self._normalized_text_expr("HugoSymbol", source_columns)),
            (
                "ensembl_gene_id",
                self._normalized_text_expr("EnsemblGeneID", source_columns),
            ),
            (
                "entrez_gene_id",
                self._normalized_text_expr("EntrezGeneID", source_columns),
            ),
            ("hgnc_name", self._normalized_text_expr("HgncName", source_columns)),
            (
                "molecular_consequence",
                self._normalized_text_expr(
                    "MolecularConsequence", source_columns
                ),
            ),
            ("vep_impact", self._normalized_text_expr("VepImpact", source_columns)),
            ("af", self._numeric_expr("AF", source_columns, "DOUBLE")),
            ("dp", self._numeric_expr("DP", source_columns, "INTEGER")),
            (
                "ref_count",
                self._numeric_expr("RefCount", source_columns, "INTEGER"),
            ),
            (
                "alt_count",
                self._numeric_expr("AltCount", source_columns, "INTEGER"),
            ),
            ("gt", self._normalized_text_expr("GT", source_columns)),
            ("sift", self._normalized_text_expr("Sift", source_columns)),
            (
                "polyphen",
                self._normalized_text_expr("Polyphen", source_columns),
            ),
            (
                "gnomade_af",
                self._numeric_expr("GnomadeAF", source_columns, "DOUBLE"),
            ),
            (
                "gnomadg_af",
                self._numeric_expr("GnomadgAF", source_columns, "DOUBLE"),
            ),
            (
                "revel_score",
                self._numeric_expr("RevelScore", source_columns, "DOUBLE"),
            ),
            (
                "likely_lof",
                self._boolean_expr("LikelyLoF", source_columns),
            ),
            ("hotspot", self._boolean_expr("Hotspot", source_columns)),
            (
                "oncogene_high_impact",
                self._boolean_expr("OncogeneHighImpact", source_columns),
            ),
            (
                "tumor_suppressor_high_impact",
                self._boolean_expr(
                    "TumorSuppressorHighImpact", source_columns
                ),
            ),
            ("hess_driver", self._boolean_expr("HessDriver", source_columns)),
        ]

        source_select = ",\n                    ".join(
            f"{expression} AS {alias}"
            for alias, expression in transformed_columns
        )
        hash_columns = [alias for alias, _ in transformed_columns]
        hash_expr = self._mutation_hash_expr(hash_columns)
        stage_columns_sql = ", ".join(hash_columns)

        sql = f"""
        CREATE OR REPLACE TEMP TABLE mutation_stage AS
        WITH source AS (
            SELECT
                {source_select}
            FROM read_csv_auto(
                ?,
                header=true,
                all_varchar=true,
                sample_size=-1
            )
        ),
        deduped AS (
            SELECT DISTINCT {stage_columns_sql}
            FROM source
            WHERE model_id IS NOT NULL
        )
        SELECT
            {hash_expr} AS mutation_id,
            {stage_columns_sql},
            CURRENT_TIMESTAMP AS created_at
        FROM deduped
        """
        conn.execute(sql, [str(file_path)])

    def _reload_tables(self, conn: duckdb.DuckDBPyConnection) -> None:
        """Atomically replace mutation tables from the staged data."""
        try:
            conn.execute("BEGIN")
            conn.execute("DELETE FROM model_gene_mutation_status")
            conn.execute("DELETE FROM mutations")
            conn.execute(
                """
                INSERT INTO mutations
                SELECT * FROM mutation_stage
                """
            )
            conn.execute(
                """
                INSERT INTO model_gene_mutation_status (
                    model_id,
                    gene_symbol,
                    is_mutated,
                    mutation_count,
                    has_likely_lof,
                    has_hotspot,
                    has_oncogene_high_impact,
                    has_tumor_suppressor_high_impact,
                    has_driver,
                    created_at
                )
                SELECT
                    model_id,
                    COALESCE(
                        CAST(hugo_symbol AS VARCHAR),
                        CAST(hgnc_name AS VARCHAR)
                    ) AS gene_symbol,
                    TRUE AS is_mutated,
                    COUNT(*) AS mutation_count,
                    BOOL_OR(likely_lof) AS has_likely_lof,
                    BOOL_OR(hotspot) AS has_hotspot,
                    BOOL_OR(oncogene_high_impact) AS has_oncogene_high_impact,
                    BOOL_OR(tumor_suppressor_high_impact)
                        AS has_tumor_suppressor_high_impact,
                    BOOL_OR(hess_driver) AS has_driver,
                    CURRENT_TIMESTAMP AS created_at
                FROM mutations
                WHERE COALESCE(
                    CAST(hugo_symbol AS VARCHAR),
                    CAST(hgnc_name AS VARCHAR)
                ) IS NOT NULL
                GROUP BY
                    model_id,
                    COALESCE(
                        CAST(hugo_symbol AS VARCHAR),
                        CAST(hgnc_name AS VARCHAR)
                    )
                """
            )
            conn.execute("COMMIT")
        except Exception:
            conn.execute("ROLLBACK")
            raise

    def _normalized_text_expr(
        self, column_name: str, source_columns: set[str]
    ) -> str:
        """Return SQL to normalize string columns to NULL when empty-ish."""
        if column_name not in source_columns:
            return "NULL"

        quoted = self._quote_identifier(column_name)
        null_checks = " OR ".join(
            f"lower(trim(CAST({quoted} AS VARCHAR))) = '{marker}'"
            for marker in self._TEXT_NULL_MARKERS
        )
        return (
            f"CASE WHEN {quoted} IS NULL OR {null_checks} "
            f"THEN NULL ELSE trim(CAST({quoted} AS VARCHAR)) END"
        )

    def _numeric_expr(
        self, column_name: str, source_columns: set[str], sql_type: str
    ) -> str:
        """Return SQL to coerce numeric columns safely."""
        if column_name not in source_columns:
            return f"CAST(NULL AS {sql_type})"

        quoted = self._quote_identifier(column_name)
        return (
            f"TRY_CAST(NULLIF(trim(CAST({quoted} AS VARCHAR)), '') AS {sql_type})"
        )

    def _boolean_expr(
        self, column_name: str, source_columns: set[str]
    ) -> str:
        """Return SQL to coerce boolean-ish columns."""
        if column_name not in source_columns:
            return "FALSE"

        quoted = self._quote_identifier(column_name)
        true_checks = ", ".join(
            f"'{marker}'" for marker in self._BOOLEAN_TRUE_MARKERS
        )
        return (
            f"CASE WHEN {quoted} IS NULL THEN FALSE "
            f"WHEN lower(trim(CAST({quoted} AS VARCHAR))) IN ({true_checks}) "
            f"THEN TRUE ELSE FALSE END"
        )

    def _mutation_hash_expr(self, columns: list[str]) -> str:
        """Build a stable mutation-id hash from normalized staged columns."""
        hash_inputs = ", ".join(
            f"COALESCE(CAST({column} AS VARCHAR), '')" for column in columns
        )
        return f"substr(md5(concat_ws('|', {hash_inputs})), 1, 16)"

    def _quote_identifier(self, identifier: str) -> str:
        """Quote a SQL identifier for DuckDB."""
        escaped = identifier.replace('"', '""')
        return f'"{escaped}"'
