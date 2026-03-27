"""Processor for CRISPR gene effect data."""

from typing import Any

import numpy as np
import pandas as pd

from .base import BaseProcessor


class GeneEffectProcessor(BaseProcessor):
    """Processor for CRISPRGeneEffect.csv files."""

    def __init__(self, batch_size: int | None = None):
        """Initialize gene effect processor."""
        super().__init__("CRISPRGeneEffect", batch_size)

    def get_table_name(self) -> str:
        """Get the target database table name."""
        return "gene_effects"

    def validate_data(
        self, df: pd.DataFrame
    ) -> tuple[pd.DataFrame, list[str]]:
        """Validate and clean CRISPR gene effect data.

        Args:
            df: Input DataFrame with models as rows and genes as columns

        Returns:
            Tuple of (cleaned_dataframe, warnings_list)
        """
        warnings: list[str] = []
        original_shape = df.shape

        # Check if first column is model IDs (should start with 'ACH-')
        first_col = df.columns[0]
        if first_col.lower() in ["unnamed: 0", "index"]:
            # First column is likely the index, use it as model_id
            df = df.set_index(df.columns[0])
            df.index.name = "model_id"
        elif not df.index.name and df.index.dtype == "object":
            # Index already contains model IDs
            df.index.name = "model_id"

        # Ensure index contains valid model IDs
        invalid_models: list[str] = []
        for model_id in df.index:
            if not isinstance(model_id, str) or not model_id.startswith(
                "ACH-"
            ):
                invalid_models.append(str(model_id))

        if invalid_models:
            warnings.append(
                f"Found {len(invalid_models)} invalid model IDs: {invalid_models[:5]}..."
            )
            # Remove invalid models
            valid_mask = df.index.str.startswith("ACH-", na=False)
            df = df[valid_mask]

        # Check gene columns - should be HUGO symbols
        gene_columns = df.columns.tolist()

        for gene in gene_columns:
            # Basic validation - gene symbols should be strings without special characters
            if not isinstance(gene, str) or len(gene) == 0:
                pass

        # Handle missing values
        missing_count = df.isna().sum().sum()
        if missing_count > 0:
            warnings.append(
                f"Found {missing_count} missing values, filling with NaN"
            )

        # Check for infinite values
        inf_count = np.isinf(df.select_dtypes(include=[np.number])).sum().sum()
        if inf_count > 0:
            warnings.append(
                f"Found {inf_count} infinite values, replacing with NaN"
            )
            df = df.replace([np.inf, -np.inf], np.nan)

        # Validate data ranges - gene effects should typically be between -5 and +2
        numeric_df = df.select_dtypes(include=[np.number])
        if not numeric_df.empty:
            min_val = numeric_df.min().min()
            max_val = numeric_df.max().max()

            if min_val < -10 or max_val > 5:
                warnings.append(
                    f"Gene effect values outside expected range: min={min_val:.2f}, max={max_val:.2f}"
                )

        self.logger.info(
            "Validation: %s -> %s, %s warnings",
            original_shape,
            df.shape,
            len(warnings),
        )
        return df, warnings

    def transform_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """Transform wide-format gene effect data to long format for database storage.

        Args:
            df: Wide-format DataFrame (models x genes)

        Returns:
            Long-format DataFrame with columns: model_id, gene_id, gene_effect
        """
        self.logger.info("Transforming wide format to long format...")

        # Reset index to make model_id a column
        df_reset = df.reset_index()

        # Melt from wide to long format
        df_long = pd.melt(
            df_reset,
            id_vars=["model_id"],
            var_name="gene_id",
            value_name="gene_effect",
        )

        # Remove rows with NaN gene effects to save space
        # (we can represent missing data by absence of records)
        initial_count = len(df_long)
        df_long = df_long.dropna(subset=["gene_effect"])
        final_count = len(df_long)

        self.logger.info(
            "Removed %s records with missing gene effects",
            initial_count - final_count,
        )

        # Clean up gene IDs - remove Entrez IDs in parentheses (e.g., "A1BG (1)" -> "A1BG")
        df_long["gene_id"] = df_long["gene_id"].str.replace(
            r"\s*\([^)]*\)", "", regex=True
        )

        # Ensure proper data types
        df_long["model_id"] = df_long["model_id"].astype(str)
        df_long["gene_id"] = df_long["gene_id"].astype(str)
        df_long["gene_effect"] = df_long["gene_effect"].astype(float)

        # Add created_at timestamp
        df_long["created_at"] = pd.Timestamp.now()

        # Sort for better insert performance
        df_long = df_long.sort_values(["model_id", "gene_id"])

        self.logger.info("Transformed to %s gene effect records", len(df_long))
        return df_long

    def get_existing_records(self) -> set[tuple[str, str]]:
        """Get set of existing (model_id, gene_id) pairs."""
        table_name = self.get_table_name()

        try:
            if not self.db_manager.table_exists(table_name):
                return set()

            result = self.db_manager.execute(
                f"SELECT model_id, gene_id FROM {table_name}"
            )
            existing = {(row[0], row[1]) for row in result.fetchall()}

            self.logger.info(
                "Found %s existing gene effect records", len(existing)
            )
            return existing

        except RuntimeError as e:
            self.logger.warning("Could not retrieve existing records: %s", e)
            return set()

    def insert_batch(self, df: pd.DataFrame) -> tuple[int, int]:
        """Insert batch of gene effect records, handling duplicates efficiently."""
        if df.empty:
            return 0, 0

        # Filter out records that already exist
        existing_records = self.get_existing_records()

        if existing_records:
            # Create mask for new records
            df["_exists"] = df.apply(
                lambda row: (row["model_id"], row["gene_id"])
                in existing_records,
                axis=1,
            )
            new_records_df = df[~df["_exists"]].drop("_exists", axis=1)
            duplicate_count = len(df) - len(new_records_df)

            if duplicate_count > 0:
                self.logger.debug(
                    "Skipping %s duplicate records in batch",
                    duplicate_count,
                )
        else:
            new_records_df = df

        if new_records_df.empty:
            return 0, 0

        # Use parent class insert logic
        return super().insert_batch(new_records_df)

    def _should_reprocess(self) -> bool:
        """Gene effect data should typically be reprocessed if source is newer."""
        # For now, always allow reprocessing
        return True

    def get_summary_stats(self) -> dict[str, Any]:
        """Get summary statistics for loaded gene effect data."""
        table_name = self.get_table_name()

        try:
            if not self.db_manager.table_exists(table_name):
                return {"error": "Table does not exist"}

            # Get basic counts
            stats_query = """
            SELECT
                COUNT(*) as total_records,
                COUNT(DISTINCT model_id) as unique_models,
                COUNT(DISTINCT gene_id) as unique_genes,
                AVG(gene_effect) as mean_gene_effect,
                MIN(gene_effect) as min_gene_effect,
                MAX(gene_effect) as max_gene_effect,
                COUNT(CASE WHEN gene_effect < -0.5 THEN 1 END) as dependency_count,
                COUNT(CASE WHEN gene_effect < -1.0 THEN 1 END) as strong_dependency_count
            FROM {table_name}
            """

            result = self.db_manager.execute(stats_query)
            row = result.fetchone()

            return {
                "total_records": row[0],
                "unique_models": row[1],
                "unique_genes": row[2],
                "mean_gene_effect": float(row[3]) if row[3] else 0,
                "min_gene_effect": float(row[4]) if row[4] else 0,
                "max_gene_effect": float(row[5]) if row[5] else 0,
                "dependency_count": row[6],
                "strong_dependency_count": row[7],
                "dependency_rate": row[6] / row[0] if row[0] > 0 else 0,
            }

        except RuntimeError as e:
            return {"error": str(e)}
