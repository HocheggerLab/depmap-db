"""Processor for Gene metadata."""

from typing import Any

import pandas as pd

from .base import BaseProcessor


class GeneProcessor(BaseProcessor):
    """Processor for Gene.csv files."""

    def __init__(self, batch_size: int | None = None):
        """Initialize gene processor."""
        super().__init__("Gene", batch_size)

    def get_table_name(self) -> str:
        """Get the target database table name."""
        return "genes"

    def validate_data(
        self, df: pd.DataFrame
    ) -> tuple[pd.DataFrame, list[str]]:
        """Validate and clean gene metadata.

        Args:
            df: Input DataFrame with gene information

        Returns:
            Tuple of (cleaned_dataframe, warnings_list)
        """
        warnings: list[str] = []
        original_shape = df.shape

        # Check required columns based on typical Gene.csv structure
        expected_columns = ["hugo_symbol", "entrez_id"]
        missing_required: list[str] = []

        # Check for key columns (case-insensitive)
        df_columns_lower = [col.lower() for col in df.columns]

        for expected in expected_columns:
            if expected not in df_columns_lower:
                # Check for common variations
                if expected == "hugo_symbol":
                    if any(
                        col in df_columns_lower
                        for col in ["gene_name", "symbol", "gene_symbol"]
                    ):
                        # Find the actual column name
                        for col in df.columns:
                            if col.lower() in [
                                "gene_name",
                                "symbol",
                                "gene_symbol",
                            ]:
                                df = df.rename(columns={col: "hugo_symbol"})
                                break
                    else:
                        missing_required.append(expected)
                elif expected == "entrez_id":
                    if any(
                        col in df_columns_lower
                        for col in ["entrez", "entrez_gene_id", "gene_id"]
                    ):
                        for col in df.columns:
                            if col.lower() in [
                                "entrez",
                                "entrez_gene_id",
                                "gene_id",
                            ]:
                                df = df.rename(columns={col: "entrez_id"})
                                break
                    else:
                        missing_required.append(expected)

        if missing_required:
            warnings.append(f"Missing required columns: {missing_required}")

        # Create gene_id column if it doesn't exist (use hugo_symbol as primary key)
        if "gene_id" not in df.columns and "hugo_symbol" in df.columns:
            df["gene_id"] = df["hugo_symbol"]
        elif "gene_id" not in df.columns:
            warnings.append(
                "No gene_id column and cannot create from hugo_symbol"
            )

        # Validate hugo_symbol
        if "hugo_symbol" in df.columns:
            # Remove rows with missing gene symbols
            missing_symbols = df["hugo_symbol"].isna().sum()
            if missing_symbols > 0:
                warnings.append(
                    f"Removing {missing_symbols} rows with missing gene symbols"
                )
                df = df.dropna(subset=["hugo_symbol"])

            # Check for duplicated gene symbols
            duplicates = df["hugo_symbol"].duplicated().sum()
            if duplicates > 0:
                warnings.append(
                    f"Found {duplicates} duplicate gene symbols, keeping first occurrence"
                )
                df = df.drop_duplicates(subset=["hugo_symbol"], keep="first")

            # Validate gene symbol format (basic check)
            invalid_symbols = (
                df["hugo_symbol"].str.len() > 50
            )  # Unusually long symbols
            if invalid_symbols.any():
                warnings.append(
                    f"Found {invalid_symbols.sum()} unusually long gene symbols"
                )

        # Validate entrez_id if present
        if "entrez_id" in df.columns:
            # Convert to numeric, handling missing values
            df["entrez_id"] = pd.to_numeric(df["entrez_id"], errors="coerce")

            # Check for missing entrez IDs
            missing_entrez = df["entrez_id"].isna().sum()
            if missing_entrez > 0:
                warnings.append(
                    f"Found {missing_entrez} genes with missing Entrez IDs"
                )

        # Validate ensembl_id if present
        if "ensembl_id" in df.columns and df["ensembl_id"].notna().any():
            # Check Ensembl ID format (should start with ENSG)
            valid_ensembl = df["ensembl_id"].str.startswith("ENSG", na=True)
            invalid_ensembl = (~valid_ensembl & df["ensembl_id"].notna()).sum()
            if invalid_ensembl > 0:
                warnings.append(
                    f"Found {invalid_ensembl} invalid Ensembl ID formats"
                )

        # Handle chromosome information
        if "chromosome" in df.columns:
            # Standardize chromosome names
            df["chromosome"] = df["chromosome"].astype(str).str.upper()
            df["chromosome"] = df["chromosome"].replace({"23": "X", "24": "Y"})

        # Clean up gene_type if present
        if "gene_type" in df.columns:
            df["gene_type"] = df["gene_type"].fillna("unknown")

        self.logger.info(
            "Validation: %s -> %s, %s warnings",
            original_shape,
            df.shape,
            len(warnings),
        )
        return df, warnings

    def transform_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """Transform gene data to match database schema.

        Args:
            df: Input DataFrame

        Returns:
            Transformed DataFrame ready for database insertion
        """
        self.logger.info("Transforming gene metadata...")

        # Ensure required columns exist with proper names
        schema_mapping = {
            "hugo_symbol": "hugo_symbol",
            "entrez_id": "entrez_id",
            "ensembl_id": "ensembl_id",
            "gene_type": "gene_type",
            "chromosome": "chromosome",
            "start_position": "start_position",
            "end_position": "end_position",
            "strand": "strand",
            "description": "description",
            "synonyms": "synonyms",
        }

        # Create transformed dataframe with only the columns we need
        transformed_df = pd.DataFrame()

        # gene_id is required (primary key)
        transformed_df["gene_id"] = df.get("gene_id", df.get("hugo_symbol"))

        # Map existing columns to schema
        for schema_col, source_col in schema_mapping.items():
            if source_col in df.columns:
                transformed_df[schema_col] = df[source_col]
            else:
                # Set defaults for missing columns
                if schema_col == "gene_type":
                    transformed_df[schema_col] = "unknown"
                elif schema_col in [
                    "entrez_id",
                    "start_position",
                    "end_position",
                ]:
                    transformed_df[schema_col] = (
                        None  # Will be handled as NULL in database
                    )
                elif schema_col == "strand":
                    transformed_df[schema_col] = None
                else:
                    transformed_df[schema_col] = None

        # Handle synonyms - convert to string array format if needed
        if (
            "synonyms" in transformed_df.columns
            and transformed_df["synonyms"].notna().any()
        ):
            # If synonyms are pipe-separated, convert to array
            def parse_synonyms(syn: object) -> list[str] | None:
                if pd.isna(syn):
                    return None
                if isinstance(syn, str) and "|" in syn:
                    return syn.split("|")
                return None

            transformed_df["synonyms"] = transformed_df["synonyms"].apply(
                parse_synonyms
            )

        # Ensure proper data types
        if "entrez_id" in transformed_df.columns:
            transformed_df["entrez_id"] = transformed_df["entrez_id"].astype(
                "Int64"
            )  # Nullable integer

        if "start_position" in transformed_df.columns:
            transformed_df["start_position"] = pd.to_numeric(
                transformed_df["start_position"], errors="coerce"
            ).astype("Int64")

        if "end_position" in transformed_df.columns:
            transformed_df["end_position"] = pd.to_numeric(
                transformed_df["end_position"], errors="coerce"
            ).astype("Int64")

        # Add timestamps
        transformed_df["created_at"] = pd.Timestamp.now()
        transformed_df["updated_at"] = pd.Timestamp.now()

        # Remove any rows with missing gene_id (our primary key)
        before_count = len(transformed_df)
        transformed_df = transformed_df.dropna(subset=["gene_id"])
        after_count = len(transformed_df)

        if before_count != after_count:
            self.logger.warning(
                "Removed %s rows with missing gene_id",
                before_count - after_count,
            )

        self.logger.info("Transformed to %s gene records", len(transformed_df))
        return transformed_df

    def get_existing_records(self) -> set[str]:
        """Get set of existing gene IDs."""
        table_name = self.get_table_name()

        try:
            if not self.db_manager.table_exists(table_name):
                return set()

            result = self.db_manager.execute(
                f"SELECT gene_id FROM {table_name}"
            )
            existing = {row[0] for row in result.fetchall()}

            self.logger.info("Found %s existing gene records", len(existing))
            return existing

        except RuntimeError as e:
            self.logger.warning("Could not retrieve existing records: %s", e)
            return set()

    def get_summary_stats(self) -> dict[str, Any]:
        """Get summary statistics for loaded gene data."""
        table_name = self.get_table_name()

        try:
            if not self.db_manager.table_exists(table_name):
                return {"error": "Table does not exist"}

            stats_query = f"""
            SELECT
                COUNT(*) as total_genes,
                COUNT(CASE WHEN entrez_id IS NOT NULL THEN 1 END) as genes_with_entrez,
                COUNT(CASE WHEN ensembl_id IS NOT NULL THEN 1 END) as genes_with_ensembl,
                COUNT(CASE WHEN gene_type IS NOT NULL AND gene_type != 'unknown' THEN 1 END) as genes_with_type,
                COUNT(CASE WHEN chromosome IS NOT NULL THEN 1 END) as genes_with_location,
                COUNT(DISTINCT gene_type) as unique_gene_types
            FROM {table_name}
            """

            result = self.db_manager.execute(stats_query)
            row = result.fetchone()

            return {
                "total_genes": row[0],
                "genes_with_entrez": row[1],
                "genes_with_ensembl": row[2],
                "genes_with_type": row[3],
                "genes_with_location": row[4],
                "unique_gene_types": row[5],
            }

        except RuntimeError as e:
            return {"error": str(e)}
