"""Optimized processor for CRISPR gene effect data."""

import gc
from pathlib import Path

import pandas as pd

from .base import BaseProcessor, ProcessingResult


class GeneEffectOptimizedProcessor(BaseProcessor):
    """Optimized processor for CRISPRGeneEffect.csv files."""

    def __init__(self, batch_size: int = 50000, nan_threshold: int = 50):
        """Initialize optimized gene effect processor.

        Args:
            batch_size: Number of records to process at once
            nan_threshold: Maximum NaN values allowed per row/column before dropping
        """
        super().__init__("CRISPRGeneEffect", batch_size)
        self.nan_threshold = nan_threshold
        self.chunk_size = 100  # Process N models at a time to save memory

    def get_table_name(self) -> str:
        """Get the target database table name."""
        return "gene_effects"

    def preprocess_and_clean(self, file_path: Path) -> pd.DataFrame:
        """Preprocess the data efficiently before loading.

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
        df = df[columns_to_keep]
        self.logger.info(
            "Kept %s/%s genes (threshold: %s NaNs)",
            len(columns_to_keep),
            initial_shape[1],
            self.nan_threshold,
        )

        # Step 2: Remove models (rows) with too many NaNs
        nan_per_row = df.isna().sum(axis=1)
        rows_to_keep = nan_per_row[nan_per_row <= self.nan_threshold].index
        df = df.loc[rows_to_keep]
        self.logger.info(
            "Kept %s/%s models (threshold: %s NaNs)",
            len(rows_to_keep),
            initial_shape[0],
            self.nan_threshold,
        )

        # Step 3: Fill remaining NaNs with column median
        column_medians = df.median()
        df = df.fillna(column_medians)
        remaining_nans = df.isna().sum().sum()
        self.logger.info(
            "Filled NaNs with column medians. Remaining NaNs: %s",
            remaining_nans,
        )

        return df

    def process_file_chunked(
        self, file_path: Path, force_reload: bool = False
    ) -> ProcessingResult:
        """Process file in chunks to avoid memory issues."""
        start_time = pd.Timestamp.now()

        try:
            # Preprocess and clean data first
            df_clean = self.preprocess_and_clean(file_path)

            # Clear existing data if force_reload
            if force_reload:
                self.logger.info("Clearing existing data...")
                self.db_manager.execute(f"DELETE FROM {self.get_table_name()}")

            total_records = 0
            model_chunks = [
                df_clean.index[i : i + self.chunk_size]
                for i in range(0, len(df_clean), self.chunk_size)
            ]

            self.logger.info(
                "Processing %s models in %s chunks...",
                len(df_clean),
                len(model_chunks),
            )

            for chunk_idx, model_batch in enumerate(model_chunks):
                # Process a chunk of models
                chunk_df = df_clean.loc[model_batch]

                # Convert to long format for this chunk only
                chunk_long = pd.melt(
                    chunk_df.reset_index(),
                    id_vars=["index"],
                    var_name="gene_id",
                    value_name="gene_effect",
                ).rename(columns={"index": "model_id"})

                # Remove any remaining NaNs (shouldn't be any after preprocessing)
                chunk_long = chunk_long.dropna()

                # Add timestamp
                chunk_long["created_at"] = pd.Timestamp.now()

                # Batch insert using DuckDB's efficient COPY
                if not chunk_long.empty:
                    # Create temporary parquet file for fast loading
                    temp_file = Path(
                        f"/tmp/gene_effects_chunk_{chunk_idx}.parquet"
                    )
                    chunk_long.to_parquet(temp_file, index=False)

                    # Use DuckDB's COPY command for fast insertion
                    self.db_manager.execute(f"""
                        INSERT INTO {self.get_table_name()}
                        SELECT * FROM read_parquet('{temp_file}')
                    """)

                    # Clean up temp file
                    temp_file.unlink()

                    total_records += len(chunk_long)

                    if (chunk_idx + 1) % 10 == 0:
                        self.logger.info(
                            "Processed %s/%s chunks, %s records so far",
                            chunk_idx + 1,
                            len(model_chunks),
                            f"{total_records:,}",
                        )

                # Free memory
                del chunk_long
                gc.collect()

            processing_time = (pd.Timestamp.now() - start_time).total_seconds()
            self.logger.info(
                "Successfully loaded %s records in %s seconds",
                f"{total_records:,}",
                f"{processing_time:.1f}",
            )

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
        """Transform handled in process_file_chunked."""
        return df

    def process_file(
        self, file_path: Path, force_reload: bool = False
    ) -> ProcessingResult:
        """Override to use chunked processing."""
        return self.process_file_chunked(file_path, force_reload)
