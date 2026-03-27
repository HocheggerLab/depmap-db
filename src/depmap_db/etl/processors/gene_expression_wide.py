"""Wide-format processor for gene expression TPM data."""

import re
from pathlib import Path
from typing import Any

import pandas as pd

from .base import BaseProcessor, ProcessingResult


class GeneExpressionWideProcessor(BaseProcessor):
    """Processor for gene expression data in wide format."""

    def __init__(self, nan_threshold: int = 50):
        """Initialize wide-format gene expression processor.

        Args:
            nan_threshold: Maximum NaN values allowed per row/column before dropping
        """
        super().__init__("GeneExpression", batch_size=None)
        self.nan_threshold = nan_threshold
        self.column_mapping: dict[str, str] = {}

    def get_table_name(self) -> str:
        """Get the target database table name."""
        return "gene_expression_wide"

    def preprocess_and_clean(self, file_path: Path) -> pd.DataFrame:
        """Preprocess the data with threshold-based NaN handling.

        This function:
        1. Loads the data
        2. Removes columns (genes) with > nan_threshold NaNs
        3. Removes rows (models) with > nan_threshold NaNs
        4. Fills remaining NaNs with column median
        """
        self.logger.info("Loading and preprocessing gene expression data...")

        # Load data
        df = pd.read_csv(file_path, low_memory=False)
        initial_shape = df.shape
        self.logger.info("Initial data shape: %s", initial_shape)

        # Identify metadata columns
        metadata_cols = [
            "ModelID",
            "SequencingID",
            "ModelConditionID",
            "IsDefaultEntryForModel",
            "IsDefaultEntryForMC",
        ]
        metadata_cols = [col for col in metadata_cols if col in df.columns]

        # Handle duplicate ModelIDs - keep only default entries
        if "ModelID" in df.columns:
            initial_models = df["ModelID"].nunique()

            if "IsDefaultEntryForModel" in df.columns:
                # Filter for default entries only (handle both boolean and string values)
                before_filter = len(df)
                # Check if values are "Yes"/"No" strings or boolean True/False
                if df["IsDefaultEntryForModel"].dtype == "object":
                    df = df[
                        df["IsDefaultEntryForModel"].str.upper() == "YES"
                    ].copy()
                else:
                    df = df[df["IsDefaultEntryForModel"]].copy()
                after_filter = len(df)
                removed = before_filter - after_filter
                if removed > 0:
                    self.logger.info(
                        "Filtered to default entries only: removed %s non-default rows",
                        removed,
                    )

            # Check for any remaining duplicates and keep first occurrence
            duplicates = df.duplicated(subset=["ModelID"], keep="first")
            if duplicates.any():
                num_duplicates = duplicates.sum()
                self.logger.warning(
                    "Found %s duplicate ModelIDs after filtering. Keeping first occurrence.",
                    num_duplicates,
                )
                df = df[~duplicates].copy()

            final_models = len(df)
            self.logger.info(
                "Unique models: %s -> %s", initial_models, final_models
            )

        # Identify gene columns (those with parentheses)
        gene_columns = [
            col
            for col in df.columns
            if col not in metadata_cols and "(" in col
        ]

        if not gene_columns:
            raise ValueError("No gene columns found in data")

        self.logger.info(
            "Found %s metadata columns and %s gene columns",
            len(metadata_cols),
            len(gene_columns),
        )

        # Store metadata for later
        df_metadata = df[metadata_cols].copy()

        # Work with gene columns only for cleaning
        df_genes = df[gene_columns].copy()

        # Clean gene names in columns (remove Entrez IDs)
        gene_name_map: dict[str, str] = {}
        for col in gene_columns:
            clean_name = re.sub(r"\s*\([^)]*\)", "", col).strip()
            gene_name_map[col] = clean_name

        df_genes.columns = pd.Index(
            [gene_name_map[col] for col in df_genes.columns]
        )

        # Step 1: Remove genes (columns) with too many NaNs
        nan_per_column = df_genes.isna().sum()
        columns_to_keep = nan_per_column[
            nan_per_column <= self.nan_threshold
        ].index
        df_genes_filtered = df_genes[columns_to_keep]
        genes_removed = len(gene_columns) - len(columns_to_keep)
        self.logger.info(
            "Removed %s genes with >%s NaNs. Kept %s genes.",
            genes_removed,
            self.nan_threshold,
            len(columns_to_keep),
        )

        # Step 2: Remove models (rows) with too many NaNs
        nan_per_row = df_genes_filtered.isna().sum(axis=1)
        rows_to_keep = nan_per_row[nan_per_row <= self.nan_threshold].index
        df_genes_clean = df_genes_filtered.loc[rows_to_keep]
        df_metadata_clean = df_metadata.loc[rows_to_keep]
        models_removed = len(df_genes_filtered) - len(rows_to_keep)
        self.logger.info(
            "Removed %s models with >%s NaNs. Kept %s models.",
            models_removed,
            self.nan_threshold,
            len(rows_to_keep),
        )

        # Step 3: Fill remaining NaNs with column median
        initial_nans = df_genes_clean.isna().sum().sum()
        if initial_nans > 0:
            column_medians = df_genes_clean.median()
            df_genes_clean = df_genes_clean.fillna(column_medians)
            self.logger.info(
                "Filled %s remaining NaNs with column medians", initial_nans
            )

        # Verify no NaNs remain
        remaining_nans = df_genes_clean.isna().sum().sum()
        if remaining_nans > 0:
            self.logger.warning(
                "%s NaNs still remain after preprocessing", remaining_nans
            )

        # Combine metadata and gene data
        df_clean = pd.concat(
            [
                df_metadata_clean.reset_index(drop=True),
                df_genes_clean.reset_index(drop=True),
            ],
            axis=1,
        )

        final_shape = df_clean.shape
        self.logger.info(
            "Final data shape: %s (%s models, %s columns removed)",
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

        # Identify metadata and gene columns
        metadata_cols = [
            "ModelID",
            "SequencingID",
            "ModelConditionID",
            "IsDefaultEntryForModel",
            "IsDefaultEntryForMC",
        ]
        metadata_cols_in_df = [
            col for col in metadata_cols if col in df.columns
        ]
        gene_columns = [
            col for col in df.columns if col not in metadata_cols_in_df
        ]

        # Build column definitions for gene columns
        gene_column_defs: list[str] = []
        self.column_mapping = {}

        for gene in gene_columns:
            # Clean gene name for use as column name
            clean_gene = (
                gene.replace("-", "_").replace(".", "_").replace(" ", "_")
            )
            # Ensure column name doesn't start with number
            if clean_gene[0].isdigit():
                clean_gene = f"gene_{clean_gene}"

            self.column_mapping[gene] = clean_gene
            gene_column_defs.append(f'"{clean_gene}" DOUBLE')

        # Create table SQL
        columns_sql = ",\n                ".join(
            [
                "model_id VARCHAR PRIMARY KEY",
                "sequencing_id VARCHAR",
                "model_condition_id VARCHAR",
                "is_default_entry BOOLEAN",
                "is_default_for_mc BOOLEAN",
                "created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP",
            ]
            + gene_column_defs
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

            # Rename metadata columns to match database schema
            metadata_mapping = {
                "ModelID": "model_id",
                "SequencingID": "sequencing_id",
                "ModelConditionID": "model_condition_id",
                "IsDefaultEntryForModel": "is_default_entry",
                "IsDefaultEntryForMC": "is_default_for_mc",
            }
            df_insert = df_insert.rename(columns=metadata_mapping)

            # Rename gene columns to match database schema
            df_insert = df_insert.rename(columns=self.column_mapping)

            # Add timestamp
            df_insert["created_at"] = pd.Timestamp.now()

            # Ensure column order matches database schema
            metadata_cols_db = [
                "model_id",
                "sequencing_id",
                "model_condition_id",
                "is_default_entry",
                "is_default_for_mc",
                "created_at",
            ]
            metadata_cols_present = [
                col for col in metadata_cols_db if col in df_insert.columns
            ]
            db_gene_columns = list(self.column_mapping.values())
            column_order = metadata_cols_present + db_gene_columns
            df_insert = df_insert[column_order]

            # Insert data using DuckDB's efficient COPY
            temp_file = Path("/tmp/gene_expression_wide.parquet")
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
        """Get summary statistics for loaded gene expression data."""
        table_name = self.get_table_name()

        try:
            if not self.db_manager.table_exists(table_name):
                return {"error": "Table does not exist"}

            # Get basic info about the table
            result = self.db_manager.execute(
                f"SELECT COUNT(*) FROM {table_name}"
            )
            model_count = result.fetchone()[0]

            # Get column count (excluding metadata columns)
            columns_result = self.db_manager.execute(f"DESCRIBE {table_name}")
            columns = columns_result.fetchall()
            metadata_cols = [
                "model_id",
                "sequencing_id",
                "model_condition_id",
                "is_default_entry",
                "is_default_for_mc",
                "created_at",
            ]
            gene_columns = [
                col[0] for col in columns if col[0] not in metadata_cols
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
