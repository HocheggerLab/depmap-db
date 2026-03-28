"""Wide-format processor for DepMap harmonized Gygi mass-spec proteomics."""

from __future__ import annotations

from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any

import pandas as pd

from ...utils.constants import DEPMAP_FILES
from ..uniprot_mapping import UniProtMapper
from .base import BaseProcessor, ProcessingResult


class ProteinExpressionMSWideProcessor(BaseProcessor):
    """Processor for the harmonized Gygi CCLE mass-spec matrix.

    This is a phase-1 import path for the DepMap `harmonized_MS_CCLE_Gygi.csv`
    dataset. It stores the wide model-by-protein matrix as-is, without inventing
    additional normalization or imputation semantics beyond the published file.
    """

    DATASET_NAME = "ProteomicsMSGygi"
    TABLE_NAME = "protein_expression_ms_wide"
    BRIDGE_TABLE_NAME = "protein_features"

    def __init__(self) -> None:
        super().__init__(self.DATASET_NAME, batch_size=None)
        self.column_mapping: dict[str, str] = {}
        self.dataset_info = DEPMAP_FILES[self.DATASET_NAME]

    def get_table_name(self) -> str:
        """Get the target database table name."""
        return self.TABLE_NAME

    def validate_data(
        self, df: pd.DataFrame
    ) -> tuple[pd.DataFrame, list[str]]:
        """Validate the source matrix structure."""
        warnings: list[str] = []
        if df.empty:
            raise ValueError("Proteomics MS file is empty")
        if len(df.columns) < 2:
            raise ValueError(
                "Proteomics MS file must contain model and protein columns"
            )

        model_column = df.columns[0]
        missing_model_ids = df[model_column].isna().sum()
        if missing_model_ids > 0:
            warnings.append(
                f"Removing {missing_model_ids} rows with missing model identifiers"
            )
            df = df.dropna(subset=[model_column]).copy()

        duplicate_model_ids = df[model_column].duplicated().sum()
        if duplicate_model_ids > 0:
            warnings.append(
                f"Found {duplicate_model_ids} duplicate model identifiers, keeping first occurrence"
            )
            df = df.drop_duplicates(subset=[model_column], keep="first").copy()

        return df, warnings

    def transform_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """Transform is handled in ``process_file`` for this wide matrix."""
        return df

    def process_file(
        self, file_path: Path, force_reload: bool = False
    ) -> ProcessingResult:
        """Process the Gygi mass-spec file and build the UniProt bridge."""
        start_time = pd.Timestamp.now()

        try:
            df = pd.read_csv(file_path, low_memory=False)
            df, warnings = self.validate_data(df)
            model_column = df.columns[0]
            protein_accessions = [str(column) for column in df.columns[1:]]

            bridge_df = self._build_bridge_dataframe(
                protein_accessions,
                force_refresh=force_reload,
            )
            self._replace_protein_features(bridge_df)

            self.create_wide_table_schema(protein_accessions)
            df_insert = df.rename(
                columns={model_column: "model_id", **self.column_mapping}
            ).copy()
            df_insert["created_at"] = pd.Timestamp.now()
            ordered_columns = [
                "model_id",
                "created_at",
                *self.column_mapping.values(),
            ]
            df_insert = df_insert[ordered_columns]

            with NamedTemporaryFile(
                suffix=".parquet", delete=False
            ) as temp_file:
                temp_path = Path(temp_file.name)
            try:
                df_insert.to_parquet(temp_path, index=False)
                self.db_manager.execute(
                    f"""
                    INSERT INTO {self.TABLE_NAME}
                    SELECT * FROM read_parquet('{temp_path}')
                    """
                )
            finally:
                temp_path.unlink(missing_ok=True)

            mapped_gene_symbols = int(bridge_df["gene_symbol"].notna().sum())
            mapped_local_genes = int(bridge_df["local_gene_id"].notna().sum())
            missing_gene_symbols = len(bridge_df) - mapped_gene_symbols
            missing_local_genes = len(bridge_df) - mapped_local_genes
            warnings.extend(
                [
                    (
                        "Imported harmonized Gygi mass-spec matrix without additional "
                        "normalization or missing-value imputation."
                    ),
                    (
                        "DepMap labels this file as harmonized, but the exact meaning of "
                        "that harmonization is not asserted here beyond the source label."
                    ),
                    (
                        f"UniProt bridge coverage: {mapped_gene_symbols:,}/"
                        f"{len(bridge_df):,} proteins mapped to a gene symbol"
                    ),
                    (
                        f"Local gene-table coverage: {mapped_local_genes:,}/"
                        f"{len(bridge_df):,} proteins mapped to a local genes row"
                    ),
                ]
            )
            if missing_gene_symbols > 0:
                warnings.append(
                    f"{missing_gene_symbols:,} proteins remained without a UniProt gene symbol"
                )
            if missing_local_genes > 0:
                warnings.append(
                    f"{missing_local_genes:,} proteins could not be resolved to the local genes table"
                )

            processing_time = (pd.Timestamp.now() - start_time).total_seconds()
            total_records = len(df_insert)
            self._record_import(file_path, total_records, "success")

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
                warnings=warnings,
            )
        except (OSError, ValueError, RuntimeError) as e:
            processing_time = (pd.Timestamp.now() - start_time).total_seconds()
            error_msg = f"Processing failed: {e}"
            self.logger.error("%s", error_msg)
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

    def create_wide_table_schema(self, protein_accessions: list[str]) -> None:
        """Create the wide table with dynamic protein columns."""
        self.db_manager.execute(f"DROP TABLE IF EXISTS {self.TABLE_NAME}")
        self.column_mapping = {
            accession: self._sanitize_column_name(accession)
            for accession in protein_accessions
        }
        protein_columns_sql = ",\n                ".join(
            f'"{column_name}" DOUBLE'
            for column_name in self.column_mapping.values()
        )
        self.db_manager.execute(
            f"""
            CREATE TABLE {self.TABLE_NAME} (
                model_id VARCHAR PRIMARY KEY,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                {protein_columns_sql},
                FOREIGN KEY (model_id) REFERENCES models(model_id)
            )
            """
        )

    def _build_bridge_dataframe(
        self, protein_accessions: list[str], *, force_refresh: bool
    ) -> pd.DataFrame:
        mapper = UniProtMapper(self.settings.depmap.cache_dir / "uniprot")
        mapping_df = mapper.build_mapping(
            protein_accessions,
            cache_name="harmonized_MS_CCLE_Gygi_uniprot_bridge",
            force_refresh=force_refresh,
        )
        mapping_df["storage_column_name"] = mapping_df[
            "protein_accession"
        ].map(self._sanitize_column_name)
        mapping_df["source_dataset"] = self.DATASET_NAME
        mapping_df["source_filename"] = self.dataset_info.filename
        mapping_df["modality"] = self.dataset_info.modality
        mapping_df["release_label"] = (
            self.dataset_info.release_label_override
            or self.settings.depmap.release_label
        )

        local_genes = self.db_manager.fetch_df(
            """
            SELECT gene_id, hugo_symbol, entrez_id, ensembl_id
            FROM genes
            """
        )
        local_by_symbol = local_genes.copy()
        local_by_symbol["symbol_key"] = (
            local_by_symbol["hugo_symbol"].astype(str).str.lower()
        )
        local_by_symbol = local_by_symbol.drop_duplicates(
            subset=["symbol_key"], keep="first"
        )

        mapping_df["gene_symbol_key"] = (
            mapping_df["gene_symbol"].astype(str).str.lower()
        )
        bridge_df = mapping_df.merge(
            local_by_symbol[
                [
                    "symbol_key",
                    "gene_id",
                    "hugo_symbol",
                    "entrez_id",
                    "ensembl_id",
                ]
            ],
            left_on="gene_symbol_key",
            right_on="symbol_key",
            how="left",
        ).rename(
            columns={
                "gene_id": "local_gene_id",
                "hugo_symbol": "local_hugo_symbol",
                "entrez_id_y": "local_entrez_id",
                "ensembl_id": "local_ensembl_id",
                "entrez_id_x": "entrez_id",
            }
        )

        alias_lookup, uniprot_lookup = self._build_local_gene_lookups(
            local_genes
        )
        self._resolve_local_genes_from_aliases(bridge_df, alias_lookup)
        self._resolve_local_genes_from_uniprot_ids(bridge_df, uniprot_lookup)

        unresolved_mask = (
            bridge_df["local_gene_id"].isna() & bridge_df["entrez_id"].notna()
        )
        if unresolved_mask.any():
            local_by_entrez = local_genes.dropna(
                subset=["entrez_id"]
            ).drop_duplicates(
                subset=["entrez_id"],
                keep="first",
            )
            entrez_matches = bridge_df.loc[
                unresolved_mask, ["protein_accession", "entrez_id"]
            ].merge(
                local_by_entrez[
                    ["entrez_id", "gene_id", "hugo_symbol", "ensembl_id"]
                ],
                on="entrez_id",
                how="left",
            )
            entrez_lookup = entrez_matches.set_index("protein_accession")
            for column_name, source_name in [
                ("local_gene_id", "gene_id"),
                ("local_hugo_symbol", "hugo_symbol"),
                ("local_ensembl_id", "ensembl_id"),
            ]:
                bridge_df.loc[unresolved_mask, column_name] = bridge_df.loc[
                    unresolved_mask, "protein_accession"
                ].map(entrez_lookup[source_name])
            bridge_df.loc[unresolved_mask, "local_entrez_id"] = bridge_df.loc[
                unresolved_mask, "entrez_id"
            ]
            bridge_df.loc[
                unresolved_mask & bridge_df["local_gene_id"].notna(),
                "mapping_method",
            ] = "entrez_fallback"

        bridge_df.loc[
            bridge_df["mapping_status"].eq("mapped")
            & bridge_df["local_gene_id"].notna(),
            "mapping_status",
        ] = "mapped_to_local_gene"
        bridge_df.loc[
            bridge_df["mapping_status"].eq("mapped")
            & bridge_df["local_gene_id"].isna(),
            "mapping_status",
        ] = "mapped_to_uniprot_only"

        bridge_df = bridge_df.drop(
            columns=["gene_symbol_key", "symbol_key"], errors="ignore"
        )
        bridge_df = bridge_df[
            [
                "protein_accession",
                "protein_accession_base",
                "storage_column_name",
                "protein_entry_name",
                "protein_name",
                "gene_symbol",
                "entrez_id",
                "ensembl_transcript_ids",
                "local_gene_id",
                "local_hugo_symbol",
                "local_entrez_id",
                "local_ensembl_id",
                "is_reviewed",
                "mapping_method",
                "mapping_status",
                "source_dataset",
                "source_filename",
                "modality",
                "release_label",
            ]
        ].sort_values("protein_accession")
        return bridge_df

    def _build_local_gene_lookups(
        self, local_genes: pd.DataFrame
    ) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, Any]]]:
        gene_cache_path = self.settings.depmap.cache_dir / "Gene.csv"
        if not gene_cache_path.exists():
            return {}, {}

        gene_metadata = pd.read_csv(
            gene_cache_path,
            low_memory=False,
            usecols=["symbol", "alias_symbol", "prev_symbol", "uniprot_ids"],
        )
        if gene_metadata.empty:
            return {}, {}

        gene_metadata = gene_metadata.rename(columns={"symbol": "hugo_symbol"})
        gene_metadata = gene_metadata.dropna(subset=["hugo_symbol"])
        gene_metadata["hugo_symbol"] = gene_metadata["hugo_symbol"].astype(str)

        local_lookup = local_genes[
            ["gene_id", "hugo_symbol", "entrez_id", "ensembl_id"]
        ].copy()
        local_lookup["hugo_symbol"] = local_lookup["hugo_symbol"].astype(str)
        merged = gene_metadata.merge(
            local_lookup, on="hugo_symbol", how="inner"
        )
        if merged.empty:
            return {}, {}

        alias_rows: list[dict[str, Any]] = []
        uniprot_rows: list[dict[str, Any]] = []
        for _, row in merged.iterrows():
            payload = {
                "local_gene_id": row["gene_id"],
                "local_hugo_symbol": row["hugo_symbol"],
                "local_entrez_id": row["entrez_id"],
                "local_ensembl_id": row["ensembl_id"],
            }
            for column in ["alias_symbol", "prev_symbol"]:
                for token in self._split_lookup_tokens(
                    row.get(column), delimiter="|"
                ):
                    alias_rows.append({"lookup_key": token.lower(), **payload})
            for accession in self._split_lookup_tokens(
                row.get("uniprot_ids"), delimiter="|"
            ):
                uniprot_rows.append({"lookup_key": accession, **payload})

        return self._build_unique_lookup(
            alias_rows
        ), self._build_unique_lookup(uniprot_rows)

    def _resolve_local_genes_from_aliases(
        self, bridge_df: pd.DataFrame, alias_lookup: dict[str, dict[str, Any]]
    ) -> None:
        if not alias_lookup:
            return

        unresolved_indices = bridge_df.index[
            bridge_df["local_gene_id"].isna()
            & bridge_df["gene_symbol"].notna()
        ]
        for idx in unresolved_indices:
            matches = {
                match["local_gene_id"]: match
                for token in self._split_lookup_tokens(
                    bridge_df.at[idx, "gene_symbol"], delimiter=";"
                )
                if (match := alias_lookup.get(token.lower())) is not None
            }
            if len(matches) != 1:
                continue
            self._apply_local_gene_match(
                bridge_df,
                idx,
                next(iter(matches.values())),
                method_suffix="gene_symbol_alias_lookup",
            )

    def _resolve_local_genes_from_uniprot_ids(
        self,
        bridge_df: pd.DataFrame,
        uniprot_lookup: dict[str, dict[str, Any]],
    ) -> None:
        if not uniprot_lookup:
            return

        unresolved_indices = bridge_df.index[bridge_df["local_gene_id"].isna()]
        for idx in unresolved_indices:
            accession = bridge_df.at[idx, "protein_accession_base"]
            if pd.isna(accession):
                continue
            match = uniprot_lookup.get(str(accession))
            if match is None:
                continue
            self._apply_local_gene_match(
                bridge_df,
                idx,
                match,
                method_suffix="uniprot_id_lookup",
            )

    def _apply_local_gene_match(
        self,
        bridge_df: pd.DataFrame,
        idx: Any,
        match: dict[str, Any],
        *,
        method_suffix: str,
    ) -> None:
        bridge_df.at[idx, "local_gene_id"] = match["local_gene_id"]
        bridge_df.at[idx, "local_hugo_symbol"] = match["local_hugo_symbol"]
        bridge_df.at[idx, "local_entrez_id"] = match["local_entrez_id"]
        bridge_df.at[idx, "local_ensembl_id"] = match["local_ensembl_id"]
        bridge_df.at[idx, "mapping_method"] = self._append_mapping_method(
            bridge_df.at[idx, "mapping_method"], method_suffix
        )

    def _build_unique_lookup(
        self, rows: list[dict[str, Any]]
    ) -> dict[str, dict[str, Any]]:
        if not rows:
            return {}

        lookup_df = pd.DataFrame(rows).drop_duplicates()
        unique_keys = (
            lookup_df.groupby("lookup_key")["local_gene_id"]
            .nunique()
            .loc[lambda counts: counts == 1]
            .index
        )
        filtered = lookup_df[lookup_df["lookup_key"].isin(unique_keys)]
        filtered = filtered.drop_duplicates(
            subset=["lookup_key"], keep="first"
        )
        result: dict[str, dict[str, Any]] = filtered.set_index(
            "lookup_key"
        ).to_dict(orient="index")
        return result

    def _split_lookup_tokens(self, value: Any, *, delimiter: str) -> list[str]:
        if value is None or pd.isna(value):
            return []
        return [
            token.strip()
            for token in str(value).split(delimiter)
            if token and token.strip()
        ]

    def _append_mapping_method(self, current: Any, suffix: str) -> str:
        if current is None or pd.isna(current) or str(current).strip() == "":
            return suffix
        current_text = str(current)
        if suffix in current_text.split("+"):
            return current_text
        return f"{current_text}+{suffix}"

    def _replace_protein_features(self, bridge_df: pd.DataFrame) -> None:
        self.db_manager.execute(
            "DELETE FROM protein_features WHERE source_dataset = ?",
            [self.DATASET_NAME],
        )
        bridge_insert = bridge_df.copy()
        bridge_insert["created_at"] = pd.Timestamp.now()
        with NamedTemporaryFile(suffix=".parquet", delete=False) as temp_file:
            temp_path = Path(temp_file.name)
        try:
            bridge_insert.to_parquet(temp_path, index=False)
            self.db_manager.execute(
                f"""
                INSERT INTO {self.BRIDGE_TABLE_NAME}
                SELECT * FROM read_parquet('{temp_path}')
                """
            )
        finally:
            temp_path.unlink(missing_ok=True)

    def _sanitize_column_name(self, accession: str) -> str:
        column_name = accession.replace("-", "_").replace(".", "_")
        if column_name[0].isdigit():
            return f"protein_{column_name}"
        return column_name

    def get_summary_stats(self) -> dict[str, Any]:
        """Get summary statistics for loaded proteomics data."""
        if not self.db_manager.table_exists(self.TABLE_NAME):
            return {"error": "Table does not exist"}

        result = self.db_manager.execute(
            f"SELECT COUNT(*) FROM {self.TABLE_NAME}"
        ).fetchone()
        model_count = int(result[0]) if result else 0
        columns = self.db_manager.execute(
            f"DESCRIBE {self.TABLE_NAME}"
        ).fetchall()
        protein_columns = [
            column[0]
            for column in columns
            if column[0] not in {"model_id", "created_at"}
        ]
        bridge_counts = self.db_manager.execute(
            """
            SELECT
                COUNT(*) AS proteins_total,
                COUNT(*) FILTER (WHERE gene_symbol IS NOT NULL) AS mapped_gene_symbols,
                COUNT(*) FILTER (WHERE local_gene_id IS NOT NULL) AS mapped_local_genes
            FROM protein_features
            WHERE source_dataset = ?
            """,
            [self.DATASET_NAME],
        ).fetchone()
        return {
            "total_models": model_count,
            "total_proteins": len(protein_columns),
            "mapped_gene_symbols": int(bridge_counts[1])
            if bridge_counts
            else 0,
            "mapped_local_genes": int(bridge_counts[2])
            if bridge_counts
            else 0,
            "table_format": "wide",
        }
