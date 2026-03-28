"""Shared helpers for phase-1 PRISM drug-sensitivity support."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from ...config import get_settings


@dataclass
class CompoundTargetResolution:
    """A provenance-aware target mapping result for one parsed target token."""

    broad_id: str
    target_ordinal: int
    target_symbol: str
    source_target_text: str
    source_dataset: str
    source_filename: str
    local_gene_id: str | None
    local_hugo_symbol: str | None
    local_entrez_id: int | None
    mapping_status: str
    mapping_method: str


class PrismCompoundMixin:
    """Helpers for compound metadata upsert and target parsing/resolution.

    Phase 1 keeps this intentionally conservative:
    - parse target strings as source assertions, not truth
    - allow compounds with no resolved local genes
    - only resolve to local genes when the mapping is unique and defensible
    """

    @property
    def _gene_cache_path(self) -> Path:
        settings = get_settings()
        return settings.depmap.cache_dir / "Gene.csv"

    def _split_target_tokens(self, value: Any) -> list[str]:
        if value is None or pd.isna(value):
            return []

        tokens: list[str] = []
        text = str(value)
        for raw_piece in text.replace(";", ",").replace("|", ",").split(","):
            piece = raw_piece.strip()
            if not piece or piece.upper() in {"NA", "N/A", "NONE", "UNKNOWN"}:
                continue
            if piece not in tokens:
                tokens.append(piece)
        return tokens

    def _split_pipe_tokens(self, value: Any) -> list[str]:
        if value is None or pd.isna(value):
            return []
        return [
            token.strip()
            for token in str(value).split("|")
            if token and token.strip()
        ]

    def _build_local_gene_lookup(
        self, local_genes: pd.DataFrame
    ) -> dict[str, dict[str, Any]]:
        lookup_rows: list[dict[str, Any]] = []

        if not local_genes.empty:
            for row in local_genes.itertuples(index=False):
                payload = {
                    "local_gene_id": row.gene_id,
                    "local_hugo_symbol": row.hugo_symbol,
                    "local_entrez_id": row.entrez_id,
                }
                symbol = getattr(row, "hugo_symbol", None)
                if symbol is not None and not pd.isna(symbol):
                    lookup_rows.append(
                        {"lookup_key": str(symbol).strip().lower(), **payload}
                    )

        gene_cache_path = self._gene_cache_path
        if gene_cache_path.exists():
            gene_metadata = pd.read_csv(
                gene_cache_path,
                low_memory=False,
                usecols=["symbol", "alias_symbol", "prev_symbol"],
            )
            gene_metadata = gene_metadata.dropna(subset=["symbol"]).copy()
            gene_metadata["symbol"] = gene_metadata["symbol"].astype(str)

            local_lookup = local_genes[["gene_id", "hugo_symbol", "entrez_id"]].copy()
            local_lookup["hugo_symbol"] = local_lookup["hugo_symbol"].astype(str)
            merged = gene_metadata.merge(
                local_lookup,
                left_on="symbol",
                right_on="hugo_symbol",
                how="inner",
            )

            for _, row in merged.iterrows():
                payload = {
                    "local_gene_id": row["gene_id"],
                    "local_hugo_symbol": row["hugo_symbol"],
                    "local_entrez_id": row["entrez_id"],
                }
                for column in ["alias_symbol", "prev_symbol"]:
                    for token in self._split_pipe_tokens(row.get(column)):
                        lookup_rows.append(
                            {"lookup_key": token.lower(), **payload}
                        )

        if not lookup_rows:
            return {}

        lookup_df = pd.DataFrame(lookup_rows).drop_duplicates()
        unique_keys = (
            lookup_df.groupby("lookup_key")["local_gene_id"]
            .nunique()
            .loc[lambda counts: counts == 1]
            .index
        )
        filtered = lookup_df[lookup_df["lookup_key"].isin(unique_keys)]
        filtered = filtered.drop_duplicates(subset=["lookup_key"], keep="first")
        return filtered.set_index("lookup_key").to_dict(orient="index")

    def _resolve_compound_targets(
        self,
        compounds_df: pd.DataFrame,
        *,
        target_column: str,
        source_dataset: str,
        source_filename: str,
    ) -> pd.DataFrame:
        if compounds_df.empty or target_column not in compounds_df.columns:
            return pd.DataFrame(
                columns=[
                    "broad_id",
                    "target_ordinal",
                    "target_symbol",
                    "source_target_text",
                    "source_dataset",
                    "source_filename",
                    "local_gene_id",
                    "local_hugo_symbol",
                    "local_entrez_id",
                    "mapping_status",
                    "mapping_method",
                ]
            )

        local_genes = self.db_manager.fetch_df(
            """
            SELECT gene_id, hugo_symbol, entrez_id
            FROM genes
            """
        )
        local_lookup = self._build_local_gene_lookup(local_genes)

        rows: list[CompoundTargetResolution] = []
        for compound in compounds_df.itertuples(index=False):
            broad_id = str(getattr(compound, "broad_id"))
            source_target_text = getattr(compound, target_column)
            tokens = self._split_target_tokens(source_target_text)
            for ordinal, token in enumerate(tokens, start=1):
                match = local_lookup.get(token.lower())
                if match is None:
                    rows.append(
                        CompoundTargetResolution(
                            broad_id=broad_id,
                            target_ordinal=ordinal,
                            target_symbol=token,
                            source_target_text="" if pd.isna(source_target_text) else str(source_target_text),
                            source_dataset=source_dataset,
                            source_filename=source_filename,
                            local_gene_id=None,
                            local_hugo_symbol=None,
                            local_entrez_id=None,
                            mapping_status="source_only_unresolved",
                            mapping_method="target_string_parse",
                        )
                    )
                    continue

                rows.append(
                    CompoundTargetResolution(
                        broad_id=broad_id,
                        target_ordinal=ordinal,
                        target_symbol=token,
                        source_target_text="" if pd.isna(source_target_text) else str(source_target_text),
                        source_dataset=source_dataset,
                        source_filename=source_filename,
                        local_gene_id=str(match["local_gene_id"]),
                        local_hugo_symbol=(
                            None
                            if pd.isna(match["local_hugo_symbol"])
                            else str(match["local_hugo_symbol"])
                        ),
                        local_entrez_id=(
                            None
                            if pd.isna(match["local_entrez_id"])
                            else int(match["local_entrez_id"])
                        ),
                        mapping_status="resolved_to_local_gene",
                        mapping_method="target_string_parse+local_gene_lookup",
                    )
                )

        if not rows:
            return pd.DataFrame(
                columns=[
                    "broad_id",
                    "target_ordinal",
                    "target_symbol",
                    "source_target_text",
                    "source_dataset",
                    "source_filename",
                    "local_gene_id",
                    "local_hugo_symbol",
                    "local_entrez_id",
                    "mapping_status",
                    "mapping_method",
                ]
            )

        return pd.DataFrame([row.__dict__ for row in rows]).drop_duplicates()

    def _replace_compounds_table(self, compounds_df: pd.DataFrame) -> None:
        if compounds_df.empty:
            return

        incoming = compounds_df.dropna(subset=["broad_id"]).copy()
        if incoming.empty:
            return

        incoming = incoming.drop_duplicates(subset=["broad_id"], keep="first")
        broad_ids = incoming["broad_id"].astype(str).tolist()
        placeholders = ", ".join("?" for _ in broad_ids)
        existing = self.db_manager.fetch_df(
            f"SELECT * FROM compounds WHERE broad_id IN ({placeholders})",
            broad_ids,
        )

        if not existing.empty:
            existing = existing.drop(columns=["created_at"], errors="ignore")
            existing = existing.set_index("broad_id")
            incoming_indexed = incoming.set_index("broad_id")
            merged = existing.combine_first(incoming_indexed)
            merged.update(incoming_indexed)
            deduped = merged.reset_index()
        else:
            deduped = incoming

        deduped["created_at"] = pd.Timestamp.now()

        temp_name = "prism_compounds_stage"
        conn = self.db_manager.connect()
        self.db_manager.execute(f"DROP VIEW IF EXISTS {temp_name}")
        conn.register(temp_name, deduped)
        try:
            self.db_manager.execute(
                f"DELETE FROM compounds WHERE broad_id IN (SELECT broad_id FROM {temp_name})"
            )
            self.db_manager.execute(
                f"INSERT INTO compounds SELECT * FROM {temp_name}"
            )
        finally:
            conn.unregister(temp_name)

    def _replace_compound_targets(
        self, target_df: pd.DataFrame, source_dataset: str
    ) -> None:
        self.db_manager.execute(
            "DELETE FROM compound_targets WHERE source_dataset = ?",
            [source_dataset],
        )
        if target_df.empty:
            return

        target_insert = target_df.copy()
        target_insert["created_at"] = pd.Timestamp.now()
        temp_name = "prism_compound_targets_stage"
        conn = self.db_manager.connect()
        self.db_manager.execute(f"DROP VIEW IF EXISTS {temp_name}")
        conn.register(temp_name, target_insert)
        try:
            self.db_manager.execute(
                f"INSERT INTO compound_targets SELECT * FROM {temp_name}"
            )
        finally:
            conn.unregister(temp_name)
