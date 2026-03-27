"""Processor for Model metadata."""

from typing import Any

import pandas as pd

from .base import BaseProcessor


class ModelProcessor(BaseProcessor):
    """Processor for Model.csv files."""

    def __init__(self, batch_size: int | None = None):
        """Initialize model processor."""
        super().__init__("Model", batch_size)

    def get_table_name(self) -> str:
        """Get the target database table name."""
        return "models"

    def validate_data(
        self, df: pd.DataFrame
    ) -> tuple[pd.DataFrame, list[str]]:
        """Validate and clean model metadata.

        Args:
            df: Input DataFrame with model information

        Returns:
            Tuple of (cleaned_dataframe, warnings_list)
        """
        warnings: list[str] = []
        original_shape = df.shape

        # Check for model_id column (primary key)
        model_id_col = None
        for col in df.columns:
            if col.lower() in ["model_id", "modelid", "depmap_id"]:
                model_id_col = col
                break

        if model_id_col is None:
            # Check if index contains model IDs
            if df.index.name and "model" in df.index.name.lower():
                df = df.reset_index()
                model_id_col = df.columns[0]
            else:
                warnings.append("No model_id column found")
                return df, warnings

        # Standardize model_id column name
        if model_id_col != "model_id":
            df = df.rename(columns={model_id_col: "model_id"})

        # Validate model IDs
        missing_model_ids = df["model_id"].isna().sum()
        if missing_model_ids > 0:
            warnings.append(
                f"Removing {missing_model_ids} rows with missing model IDs"
            )
            df = df.dropna(subset=["model_id"])

        # Check for duplicate model IDs
        duplicates = df["model_id"].duplicated().sum()
        if duplicates > 0:
            warnings.append(
                f"Found {duplicates} duplicate model IDs, keeping first occurrence"
            )
            df = df.drop_duplicates(subset=["model_id"], keep="first")

        # Validate model ID format (should be ACH-XXXXXX)
        if not df.empty:
            valid_format = df["model_id"].str.match(r"^ACH-\d+$", na=False)
            invalid_format = (~valid_format).sum()
            if invalid_format > 0:
                warnings.append(
                    f"Found {invalid_format} model IDs with non-standard format"
                )

        # Standardize column names to match schema
        column_mapping = {
            "celllinename": "cell_line_name",
            "cell_line_name": "cell_line_name",
            "strippedcelllinename": "stripped_cell_line_name",
            "stripped_cell_line_name": "stripped_cell_line_name",
            "depmapmodeltype": "depmap_model_type",
            "depmap_model_type": "depmap_model_type",
            "oncotreelineage": "oncotree_lineage",
            "oncotree_lineage": "oncotree_lineage",
            "oncotreeprimarydisease": "oncotree_primary_disease",
            "oncotree_primary_disease": "oncotree_primary_disease",
            "oncotreesubtype": "oncotree_subtype",
            "oncotree_subtype": "oncotree_subtype",
            "oncotreecode": "oncotree_code",
            "oncotree_code": "oncotree_code",
            "patientid": "patient_id",
            "patient_id": "patient_id",
            "agecategory": "age_category",
            "age_category": "age_category",
            "primaryormetastasis": "primary_or_metastasis",
            "primary_or_metastasis": "primary_or_metastasis",
            "samplecollectionsite": "sample_collection_site",
            "sample_collection_site": "sample_collection_site",
            "sourcetype": "source_type",
            "source_type": "source_type",
            "sourcedetail": "source_detail",
            "source_detail": "source_detail",
            "modeltype": "model_type",
            "model_type": "model_type",
            "tissueorigin": "tissue_origin",
            "tissue_origin": "tissue_origin",
            "growthpattern": "growth_pattern",
            "growth_pattern": "growth_pattern",
            "cclename": "ccle_name",
            "ccle_name": "ccle_name",
            "cosmicid": "cosmic_id",
            "cosmic_id": "cosmic_id",
        }

        # Apply column name mapping (case-insensitive)
        df_columns_lower = {col.lower(): col for col in df.columns}
        for standard_name, schema_name in column_mapping.items():
            if standard_name in df_columns_lower:
                old_col = df_columns_lower[standard_name]
                if old_col != schema_name:
                    df = df.rename(columns={old_col: schema_name})

        # Validate required fields
        if "cell_line_name" in df.columns:
            missing_names = df["cell_line_name"].isna().sum()
            if missing_names > 0:
                warnings.append(
                    f"Found {missing_names} models with missing cell line names"
                )

        # Clean up categorical fields
        categorical_fields = [
            "age_category",
            "sex",
            "primary_or_metastasis",
            "source_type",
            "model_type",
            "tissue_origin",
            "growth_pattern",
        ]

        for field in categorical_fields:
            if field in df.columns:
                # Standardize values
                df[field] = df[field].astype(str).str.strip()
                df[field] = df[field].replace(
                    {"nan": None, "NaN": None, "": None}
                )

        # Validate age if present
        if "age" in df.columns:
            df["age"] = pd.to_numeric(df["age"], errors="coerce")
            invalid_ages = ((df["age"] < 0) | (df["age"] > 120)).sum()
            if invalid_ages > 0:
                warnings.append(
                    f"Found {invalid_ages} models with invalid ages"
                )
                df.loc[(df["age"] < 0) | (df["age"] > 120), "age"] = None

        self.logger.info(
            "Validation: %s -> %s, %s warnings",
            original_shape,
            df.shape,
            len(warnings),
        )
        return df, warnings

    def transform_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """Transform model data to match database schema.

        Args:
            df: Input DataFrame

        Returns:
            Transformed DataFrame ready for database insertion
        """
        self.logger.info("Transforming model metadata...")

        # Define the schema columns we want to keep
        schema_columns = [
            "model_id",
            "patient_id",
            "cell_line_name",
            "stripped_cell_line_name",
            "depmap_model_type",
            "oncotree_lineage",
            "oncotree_primary_disease",
            "oncotree_subtype",
            "oncotree_code",
            "rrid",
            "age",
            "age_category",
            "sex",
            "patient_race",
            "primary_or_metastasis",
            "sample_collection_site",
            "source_type",
            "source_detail",
            "model_type",
            "tissue_origin",
            "growth_pattern",
            "ccle_name",
            "cosmic_id",
        ]

        # Create transformed dataframe
        transformed_df = pd.DataFrame()

        # Copy available columns, fill missing ones with None
        for col in schema_columns:
            if col in df.columns:
                transformed_df[col] = df[col]
            else:
                transformed_df[col] = None

        # Ensure model_id exists (required)
        if (
            "model_id" not in transformed_df.columns
            or transformed_df["model_id"].isna().all()
        ):
            raise ValueError("model_id column is required and missing")

        # Set cell_line_name from model_id if missing
        if (
            "cell_line_name" not in df.columns
            or transformed_df["cell_line_name"].isna().all()
        ):
            # Extract from model_id (ACH-000001 -> something more readable)
            transformed_df["cell_line_name"] = transformed_df["model_id"]

        # Ensure proper data types
        if "age" in transformed_df.columns:
            transformed_df["age"] = transformed_df["age"].astype(
                "Int64"
            )  # Nullable integer

        # Handle string fields - replace empty strings with None
        string_fields = [
            "patient_id",
            "cell_line_name",
            "stripped_cell_line_name",
            "depmap_model_type",
            "oncotree_lineage",
            "oncotree_primary_disease",
            "oncotree_subtype",
            "oncotree_code",
            "rrid",
            "age_category",
            "sex",
            "patient_race",
            "primary_or_metastasis",
            "sample_collection_site",
            "source_type",
            "source_detail",
            "model_type",
            "tissue_origin",
            "growth_pattern",
            "ccle_name",
            "cosmic_id",
        ]

        for field in string_fields:
            if field in transformed_df.columns:
                transformed_df[field] = transformed_df[field].replace("", None)
                transformed_df[field] = transformed_df[field].replace(
                    "nan", None
                )

        # Add timestamps
        transformed_df["created_at"] = pd.Timestamp.now()
        transformed_df["updated_at"] = pd.Timestamp.now()

        # Remove any rows with missing model_id
        before_count = len(transformed_df)
        transformed_df = transformed_df.dropna(subset=["model_id"])
        after_count = len(transformed_df)

        if before_count != after_count:
            self.logger.warning(
                "Removed %s rows with missing model_id",
                before_count - after_count,
            )

        self.logger.info(
            "Transformed to %s model records", len(transformed_df)
        )
        return transformed_df

    def get_existing_records(self) -> set[str]:
        """Get set of existing model IDs."""
        table_name = self.get_table_name()

        try:
            if not self.db_manager.table_exists(table_name):
                return set()

            result = self.db_manager.execute(
                f"SELECT model_id FROM {table_name}"
            )
            existing = {row[0] for row in result.fetchall()}

            self.logger.info("Found %s existing model records", len(existing))
            return existing

        except RuntimeError as e:
            self.logger.warning("Could not retrieve existing records: %s", e)
            return set()

    def get_summary_stats(self) -> dict[str, Any]:
        """Get summary statistics for loaded model data."""
        table_name = self.get_table_name()

        try:
            if not self.db_manager.table_exists(table_name):
                return {"error": "Table does not exist"}

            stats_query = f"""
            SELECT
                COUNT(*) as total_models,
                COUNT(DISTINCT oncotree_lineage) as unique_lineages,
                COUNT(DISTINCT oncotree_primary_disease) as unique_diseases,
                COUNT(CASE WHEN primary_or_metastasis = 'Primary' THEN 1 END) as primary_tumors,
                COUNT(CASE WHEN primary_or_metastasis = 'Metastasis' THEN 1 END) as metastatic_tumors,
                COUNT(CASE WHEN age IS NOT NULL THEN 1 END) as models_with_age,
                AVG(age) as avg_age
            FROM {table_name}
            """

            result = self.db_manager.execute(stats_query)
            row = result.fetchone()

            # Get top lineages
            lineage_query = f"""
            SELECT oncotree_lineage, COUNT(*) as count
            FROM {table_name}
            WHERE oncotree_lineage IS NOT NULL
            GROUP BY oncotree_lineage
            ORDER BY count DESC
            LIMIT 10
            """

            lineage_result = self.db_manager.execute(lineage_query)
            top_lineages = {
                row[0]: row[1] for row in lineage_result.fetchall()
            }

            return {
                "total_models": row[0],
                "unique_lineages": row[1],
                "unique_diseases": row[2],
                "primary_tumors": row[3],
                "metastatic_tumors": row[4],
                "models_with_age": row[5],
                "avg_age": float(row[6]) if row[6] else None,
                "top_lineages": top_lineages,
            }

        except RuntimeError as e:
            return {"error": str(e)}
