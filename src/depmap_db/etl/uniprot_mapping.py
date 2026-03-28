"""Helpers for building a reproducible UniProt-to-gene bridge."""

from __future__ import annotations

import csv
import io
import time
from dataclasses import dataclass
from pathlib import Path

import httpx
import pandas as pd

from ..config import get_logger

logger = get_logger(__name__)


@dataclass(frozen=True)
class UniProtMappingConfig:
    """Configuration for UniProt bridge resolution."""

    batch_size: int = 50
    pause_seconds: float = 0.1
    timeout_seconds: int = 120


class UniProtMapper:
    """Resolve UniProt accessions to basic gene/protein annotations."""

    SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
    FIELD_NAMES = (
        "accession,id,gene_primary,protein_name,xref_ensembl,"
        "xref_geneid,reviewed"
    )

    def __init__(self, cache_dir: Path, config: UniProtMappingConfig | None = None) -> None:
        self.cache_dir = cache_dir
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.config = config or UniProtMappingConfig()

    def build_mapping(
        self,
        accessions: list[str],
        *,
        cache_name: str,
        force_refresh: bool = False,
    ) -> pd.DataFrame:
        """Build or load a cached mapping table for the requested accessions."""
        cache_path = self.cache_dir / f"{cache_name}.tsv"
        if cache_path.exists() and not force_refresh:
            df = pd.read_csv(cache_path, sep="\t")
            cached_accessions = set(df["protein_accession"].astype(str))
            missing = sorted(set(accessions) - cached_accessions)
            if not missing:
                return df

        df = self._fetch_mapping(accessions)
        df.to_csv(cache_path, sep="\t", index=False)
        return df

    def _fetch_mapping(self, accessions: list[str]) -> pd.DataFrame:
        """Fetch UniProt annotations for the requested accessions."""
        unique_accessions = list(dict.fromkeys(accessions))
        exact_rows = self._query_batches(unique_accessions, query_mode="exact")
        exact_df = self._normalize_mapping_rows(exact_rows)
        found_exact = set(exact_df["protein_accession"]) if not exact_df.empty else set()

        unmapped = [acc for acc in unique_accessions if acc not in found_exact]
        isoform_unmapped = [acc for acc in unmapped if "-" in acc]
        fallback_rows: list[dict[str, object]] = []
        if isoform_unmapped:
            base_accessions = list(
                dict.fromkeys(acc.split("-", 1)[0] for acc in isoform_unmapped)
            )
            base_df = self._normalize_mapping_rows(
                self._query_batches(base_accessions, query_mode="base_accession")
            )
            base_lookup = {
                str(row["protein_accession"]): row
                for _, row in base_df.iterrows()
            }
            for accession in isoform_unmapped:
                base_accession = accession.split("-", 1)[0]
                base_row = base_lookup.get(base_accession)
                if base_row is None:
                    continue
                row = base_row.to_dict()
                row["protein_accession"] = accession
                row["protein_accession_base"] = base_accession
                row["mapping_method"] = "base_accession_fallback"
                fallback_rows.append(row)

        combined = pd.concat(
            [exact_df, pd.DataFrame(fallback_rows)],
            ignore_index=True,
            sort=False,
        )
        if combined.empty:
            combined = pd.DataFrame(columns=self._output_columns())

        combined = combined.drop_duplicates(subset=["protein_accession"], keep="first")

        missing_rows = [
            {
                "protein_accession": accession,
                "protein_accession_base": accession.split("-", 1)[0],
                "protein_entry_name": None,
                "protein_name": None,
                "gene_symbol": None,
                "entrez_id": None,
                "ensembl_transcript_ids": None,
                "is_reviewed": None,
                "mapping_method": "unmapped",
                "mapping_status": "unmapped",
            }
            for accession in unique_accessions
            if accession not in set(combined["protein_accession"].astype(str))
        ]
        if missing_rows:
            combined = pd.concat(
                [combined, pd.DataFrame(missing_rows)],
                ignore_index=True,
                sort=False,
            )

        combined = combined[self._output_columns()]
        combined = combined.sort_values("protein_accession").reset_index(drop=True)
        return combined

    def _query_batches(
        self, accessions: list[str], *, query_mode: str
    ) -> list[dict[str, str]]:
        rows: list[dict[str, str]] = []
        if not accessions:
            return rows

        with httpx.Client(timeout=self.config.timeout_seconds) as client:
            for start in range(0, len(accessions), self.config.batch_size):
                batch = accessions[start : start + self.config.batch_size]
                query = " OR ".join(f"accession:{accession}" for accession in batch)
                response = client.get(
                    self.SEARCH_URL,
                    params={
                        "query": query,
                        "format": "tsv",
                        "fields": self.FIELD_NAMES,
                        "size": str(len(batch)),
                    },
                )
                response.raise_for_status()
                reader = csv.DictReader(io.StringIO(response.text), delimiter="\t")
                for row in reader:
                    row["query_mode"] = query_mode
                    rows.append(row)
                time.sleep(self.config.pause_seconds)

        return rows

    def _normalize_mapping_rows(
        self, rows: list[dict[str, str]]
    ) -> pd.DataFrame:
        normalized_rows: list[dict[str, object]] = []
        for row in rows:
            accession = row.get("Entry") or None
            if accession is None:
                continue
            entrez_id = self._parse_first_int(row.get("GeneID"))
            normalized_rows.append(
                {
                    "protein_accession": accession,
                    "protein_accession_base": accession.split("-", 1)[0],
                    "protein_entry_name": self._clean_text(row.get("Entry Name")),
                    "protein_name": self._clean_text(row.get("Protein names")),
                    "gene_symbol": self._clean_text(row.get("Gene Names (primary)")),
                    "entrez_id": entrez_id,
                    "ensembl_transcript_ids": self._normalize_ensembl_list(
                        row.get("Ensembl")
                    ),
                    "is_reviewed": self._parse_reviewed(row.get("Reviewed")),
                    "mapping_method": row.get("query_mode") or "exact",
                    "mapping_status": "mapped",
                }
            )

        if not normalized_rows:
            return pd.DataFrame(columns=self._output_columns())
        return pd.DataFrame(normalized_rows)

    def _normalize_ensembl_list(self, value: str | None) -> str | None:
        cleaned = self._clean_text(value)
        if cleaned is None:
            return None
        tokens = [token.strip() for token in cleaned.split(";") if token.strip()]
        normalized_tokens = [token.split(" ", 1)[0] for token in tokens]
        return ";".join(normalized_tokens) or None

    def _parse_reviewed(self, value: str | None) -> bool | None:
        cleaned = self._clean_text(value)
        if cleaned is None:
            return None
        return cleaned.lower() == "reviewed"

    def _parse_first_int(self, value: str | None) -> int | None:
        cleaned = self._clean_text(value)
        if cleaned is None:
            return None
        first = cleaned.split(";", 1)[0].strip()
        if not first:
            return None
        try:
            return int(first)
        except ValueError:
            return None

    def _clean_text(self, value: str | None) -> str | None:
        if value is None:
            return None
        cleaned = value.strip()
        return cleaned or None

    def _output_columns(self) -> list[str]:
        return [
            "protein_accession",
            "protein_accession_base",
            "protein_entry_name",
            "protein_name",
            "gene_symbol",
            "entrez_id",
            "ensembl_transcript_ids",
            "is_reviewed",
            "mapping_method",
            "mapping_status",
        ]
