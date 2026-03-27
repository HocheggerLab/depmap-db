"""Processor for gene expression TPM data (long format)."""

import re

import numpy as np
import pandas as pd

from .base import BaseProcessor


class GeneExpressionProcessor(BaseProcessor):
    """Processor for gene expression data in long format."""

    def __init__(self, batch_size: int | None = None):
        """Initialize gene expression processor."""
        super().__init__("GeneExpression", batch_size)

    def get_table_name(self) -> str:
        """Get the target database table name."""
        return "gene_expression"

    def validate_data(
        self, df: pd.DataFrame
    ) -> tuple[pd.DataFrame, list[str]]:
        """Validate and clean gene expression data.

        Args:
            df: Input DataFrame with models as rows and genes as columns

        Returns:
            Tuple of (cleaned_dataframe, warnings_list)
        """
        warnings: list[str] = []
        original_shape = df.shape

        # Check for required metadata columns
        required_cols = ["ModelID", "SequencingID", "ModelConditionID"]
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            warnings.append(f"Missing required columns: {missing_cols}")

        # Validate ModelID format
        if "ModelID" in df.columns:
            invalid_models = df[
                ~df["ModelID"].str.startswith("ACH-", na=False)
            ]["ModelID"].tolist()
            if invalid_models:
                warnings.append(
                    f"Found {len(invalid_models)} invalid model IDs"
                )
                df = df[df["ModelID"].str.startswith("ACH-", na=False)]

        # Identify gene columns (those with parentheses indicating Entrez IDs)
        gene_columns = [
            col
            for col in df.columns
            if isinstance(col, str) and "(" in col and ")" in col
        ]

        if not gene_columns:
            warnings.append("No gene columns found in data")
            return df, warnings

        # Check expression value ranges (log(TPM+1) should be >= 0)
        for col in gene_columns[:10]:  # Sample first 10 genes
            if col in df.columns:
                min_val = df[col].min()
                max_val = df[col].max()

                if min_val < 0:
                    warnings.append(
                        f"Negative expression values found in {col}: min={min_val:.3f}"
                    )
                if max_val > 20:
                    warnings.append(
                        f"Unusually high expression in {col}: max={max_val:.3f}"
                    )

        # Handle missing values
        missing_count = df[gene_columns].isna().sum().sum()
        if missing_count > 0:
            warnings.append(
                f"Found {missing_count:,} missing expression values"
            )

        # Handle infinite values
        inf_count = 0
        for col in gene_columns:
            if df[col].dtype in [np.float64, np.float32]:
                inf_mask = np.isinf(df[col])
                if inf_mask.any():
                    inf_count += inf_mask.sum()
                    df.loc[inf_mask, col] = np.nan

        if inf_count > 0:
            warnings.append(f"Replaced {inf_count:,} infinite values with NaN")

        self.logger.info(
            "Validation: %s -> %s, %s warnings",
            original_shape,
            df.shape,
            len(warnings),
        )
        return df, warnings

    def transform_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """Transform wide-format expression data to long format.

        Args:
            df: Wide-format DataFrame (models x genes)

        Returns:
            Long-format DataFrame with columns: model_id, gene_id, expression_value, etc.
        """
        self.logger.info("Transforming wide format to long format...")

        # Identify metadata and gene columns
        metadata_cols = [
            "ModelID",
            "SequencingID",
            "ModelConditionID",
            "IsDefaultEntryForModel",
            "IsDefaultEntryForMC",
        ]
        metadata_cols = [col for col in metadata_cols if col in df.columns]

        gene_columns = [
            col
            for col in df.columns
            if col not in metadata_cols and "(" in col
        ]

        # Clean gene names - remove Entrez IDs in parentheses
        gene_name_map: dict[str, str] = {}
        for col in gene_columns:
            clean_name = re.sub(r"\s*\([^)]*\)", "", col).strip()
            gene_name_map[col] = clean_name

        # Melt to long format
        df_long = pd.melt(
            df,
            id_vars=metadata_cols,
            value_vars=gene_columns,
            var_name="gene_raw",
            value_name="expression_value",
        )

        # Map to clean gene names
        df_long["gene_id"] = df_long["gene_raw"].map(gene_name_map)
        df_long = df_long.drop("gene_raw", axis=1)

        # Remove rows with missing expression values
        initial_count = len(df_long)
        df_long = df_long.dropna(subset=["expression_value"])
        final_count = len(df_long)

        self.logger.info(
            "Removed %s records with missing values",
            f"{initial_count - final_count:,}",
        )

        # Rename columns to match schema
        column_mapping = {
            "ModelID": "model_id",
            "SequencingID": "sequencing_id",
            "ModelConditionID": "model_condition_id",
            "IsDefaultEntryForModel": "is_default_entry",
            "IsDefaultEntryForMC": "is_default_for_mc",
        }
        df_long = df_long.rename(columns=column_mapping)

        # Ensure proper data types
        df_long["model_id"] = df_long["model_id"].astype(str)
        df_long["gene_id"] = df_long["gene_id"].astype(str)
        df_long["expression_value"] = df_long["expression_value"].astype(float)

        if "is_default_entry" in df_long.columns:
            df_long["is_default_entry"] = df_long["is_default_entry"].astype(
                bool
            )
        if "is_default_for_mc" in df_long.columns:
            df_long["is_default_for_mc"] = df_long["is_default_for_mc"].astype(
                bool
            )

        # Add timestamp
        df_long["created_at"] = pd.Timestamp.now()

        # Sort for better insert performance
        df_long = df_long.sort_values(["model_id", "gene_id"])

        self.logger.info(
            "Transformed to %s expression records", f"{len(df_long):,}"
        )

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
                "Found %s existing expression records",
                f"{len(existing):,}",
            )
            return existing

        except RuntimeError as e:
            self.logger.warning("Could not retrieve existing records: %s", e)
            return set()
