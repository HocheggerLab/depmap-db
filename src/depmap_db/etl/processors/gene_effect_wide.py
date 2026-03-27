"""Wide-format processor for CRISPR gene effect data."""

from pathlib import Path
from typing import Any

import pandas as pd

from .base import BaseProcessor, ProcessingResult


class GeneEffectWideProcessor(BaseProcessor):
    """Processor for CRISPRGeneEffect.csv files - stores in wide format with threshold-based NaN handling."""

    def __init__(self, nan_threshold: int = 50):
        """Initialize wide-format gene effect processor.

        Args:
            nan_threshold: Maximum NaN values allowed per row/column before dropping
        """
        super().__init__("CRISPRGeneEffect", batch_size=None)
        self.nan_threshold = nan_threshold
        self.column_mapping: dict[str, str] = {}

    def get_table_name(self) -> str:
        """Get the target database table name."""
        return "gene_effects_wide"

    def preprocess_and_clean(self, file_path: Path) -> pd.DataFrame:
        """Preprocess the data with threshold-based NaN handling.

        This function:
        1. Loads the data
        2. Removes columns (genes) with > nan_threshold NaNs
        3. Removes rows (models) with > nan_threshold NaNs
        4. Fills remaining NaNs with column median
        """
        self.logger.info("Loading and preprocessing data...")

        # Load data
        df = pd.read_csv(file_path, index_col=0)
        initial_shape = df.shape
        self.logger.info("Initial data shape: %s", initial_shape)

        # Clean gene names in columns (remove Entrez IDs)
        df.columns = df.columns.str.replace(r"\s*\([^)]*\)", "", regex=True)

        # Step 1: Remove genes (columns) with too many NaNs
        nan_per_column = df.isna().sum()
        columns_to_keep = nan_per_column[
            nan_per_column <= self.nan_threshold
        ].index
        df_filtered = df[columns_to_keep]
        genes_removed = initial_shape[1] - len(columns_to_keep)
        self.logger.info(
            "Removed %s genes with >%s NaNs. Kept %s genes.",
            genes_removed,
            self.nan_threshold,
            len(columns_to_keep),
        )

        # Step 2: Remove models (rows) with too many NaNs
        nan_per_row = df_filtered.isna().sum(axis=1)
        rows_to_keep = nan_per_row[nan_per_row <= self.nan_threshold].index
        df_clean = df_filtered.loc[rows_to_keep]
        models_removed = len(df_filtered) - len(rows_to_keep)
        self.logger.info(
            "Removed %s models with >%s NaNs. Kept %s models.",
            models_removed,
            self.nan_threshold,
            len(rows_to_keep),
        )

        # Step 3: Fill remaining NaNs with column median
        initial_nans = df_clean.isna().sum().sum()
        if initial_nans > 0:
            column_medians = df_clean.median()
            df_clean = df_clean.fillna(column_medians)
            self.logger.info(
                "Filled %s remaining NaNs with column medians", initial_nans
            )

        # Verify no NaNs remain
        remaining_nans = df_clean.isna().sum().sum()
        if remaining_nans > 0:
            self.logger.warning(
                "%s NaNs still remain after preprocessing", remaining_nans
            )

        final_shape = df_clean.shape
        self.logger.info(
            "Final data shape: %s (%s models, %s genes removed)",
            final_shape,
            initial_shape[0] - final_shape[0],
            initial_shape[1] - final_shape[1],
        )

        return df_clean

    def create_wide_table_schema(self, df: pd.DataFrame) -> None:
        """Create the wide table with dynamic gene columns."""
        table_name = self.get_table_name()

        # Drop existing table if it exists
        self.db_manager.execute(f"DROP TABLE IF EXISTS {table_name}")

        # Build column definitions
        gene_columns: list[str] = []
        for gene in df.columns:
            # Clean gene name for use as column name
            clean_gene = (
                gene.replace("-", "_").replace(".", "_").replace(" ", "_")
            )
            # Ensure column name doesn't start with number
            if clean_gene[0].isdigit():
                clean_gene = f"gene_{clean_gene}"
            gene_columns.append(f'"{clean_gene}" DOUBLE')

        # Create table SQL
        columns_sql = ",\n                ".join(
            [
                "model_id VARCHAR PRIMARY KEY",
                "created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP",
            ]
            + gene_columns
        )

        create_sql = f"""
        CREATE TABLE {table_name} (
                {columns_sql}
        )
        """

        self.logger.info(
            "Creating wide table with %s gene columns...", len(gene_columns)
        )
        self.db_manager.execute(create_sql)
        self.logger.info("Successfully created %s table", table_name)

        # Store the column mapping for later use
        self.column_mapping = {
            gene: gene.replace("-", "_").replace(".", "_").replace(" ", "_")
            for gene in df.columns
        }
        # Ensure column names don't start with numbers
        for original, clean in self.column_mapping.items():
            if clean[0].isdigit():
                self.column_mapping[original] = f"gene_{clean}"

    def process_file(
        self, file_path: Path, force_reload: bool = False
    ) -> ProcessingResult:
        """Process file and store in wide format."""
        start_time = pd.Timestamp.now()

        try:
            # Preprocess and clean data
            df_clean = self.preprocess_and_clean(file_path)

            # Create the wide table schema
            self.create_wide_table_schema(df_clean)

            # Prepare data for insertion
            df_insert = df_clean.copy()
            df_insert = df_insert.reset_index()  # Make model_id a column
            df_insert = df_insert.rename(columns={"index": "model_id"})

            # Rename gene columns to match database schema
            df_insert = df_insert.rename(columns=self.column_mapping)

            # Add timestamp
            df_insert["created_at"] = pd.Timestamp.now()

            # Ensure column order matches database schema: model_id, created_at, then gene columns
            db_gene_columns = [
                self.column_mapping[gene] for gene in df_clean.columns
            ]
            column_order = ["model_id", "created_at"] + db_gene_columns
            df_insert = df_insert[column_order]

            # Insert data using DuckDB's efficient COPY
            temp_file = Path("/tmp/gene_effects_wide.parquet")
            df_insert.to_parquet(temp_file, index=False)

            self.logger.info("Inserting %s model records...", len(df_insert))
            self.db_manager.execute(f"""
                INSERT INTO {self.get_table_name()}
                SELECT * FROM read_parquet('{temp_file}')
            """)

            # Clean up temp file
            temp_file.unlink()

            processing_time = (pd.Timestamp.now() - start_time).total_seconds()
            total_records = len(df_insert)

            self.logger.info(
                "Successfully loaded %s model records with %s genes in %s seconds",
                f"{total_records:,}",
                f"{len(self.column_mapping):,}",
                f"{processing_time:.1f}",
            )

            # Return processing result
            return ProcessingResult(
                processor_name=self.__class__.__name__,
                dataset_name=self.dataset_name,
                source_file=file_path,
                records_processed=total_records,
                records_inserted=total_records,
                records_updated=0,
                records_skipped=0,
                processing_time_seconds=processing_time,
                status="success",
            )

        except (ValueError, OSError, RuntimeError) as e:
            self.logger.error("Processing failed: %s", e)
            raise

    def validate_data(
        self, df: pd.DataFrame
    ) -> tuple[pd.DataFrame, list[str]]:
        """Minimal validation since we handle cleaning in preprocess."""
        warnings: list[str] = []
        return df, warnings

    def transform_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """Transform handled in process_file."""
        return df

    def get_summary_stats(self) -> dict[str, Any]:
        """Get summary statistics for loaded gene effect data."""
        table_name = self.get_table_name()

        try:
            if not self.db_manager.table_exists(table_name):
                return {"error": "Table does not exist"}

            # Get basic info about the table
            result = self.db_manager.execute(
                f"SELECT COUNT(*) FROM {table_name}"
            )
            model_count = result.fetchone()[0]

            # Get column count (excluding model_id and created_at)
            columns_result = self.db_manager.execute(f"DESCRIBE {table_name}")
            columns = columns_result.fetchall()
            gene_columns = [
                col[0]
                for col in columns
                if col[0] not in ["model_id", "created_at"]
            ]
            gene_count = len(gene_columns)

            # Sample some statistics from a subset of gene columns (first 10)
            sample_genes = gene_columns[:10]
            stats_columns = ", ".join(
                [
                    f"AVG({gene}) as avg_{gene}, MIN({gene}) as min_{gene}, MAX({gene}) as max_{gene}"
                    for gene in sample_genes
                ]
            )

            if stats_columns:
                stats_result = self.db_manager.execute(
                    f"SELECT {stats_columns} FROM {table_name}"
                )
                stats_row = stats_result.fetchone()
                sample_stats: dict[str, dict[str, float]] = {}
                for i, gene in enumerate(sample_genes):
                    sample_stats[gene] = {
                        "avg": float(stats_row[i * 3])
                        if stats_row[i * 3]
                        else 0,
                        "min": float(stats_row[i * 3 + 1])
                        if stats_row[i * 3 + 1]
                        else 0,
                        "max": float(stats_row[i * 3 + 2])
                        if stats_row[i * 3 + 2]
                        else 0,
                    }
            else:
                sample_stats = {}

            return {
                "total_models": model_count,
                "total_genes": gene_count,
                "table_format": "wide",
                "sample_gene_stats": sample_stats,
            }

        except RuntimeError as e:
            return {"error": str(e)}
