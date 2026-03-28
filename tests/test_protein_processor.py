"""Tests for proteomics bridge resolution helpers."""

from pathlib import Path

import pandas as pd

from depmap_db.config import reload_settings
from depmap_db.etl.processors.protein_expression_ms_wide import (
    ProteinExpressionMSWideProcessor,
)


def _make_processor(tmp_path: Path) -> ProteinExpressionMSWideProcessor:
    settings = reload_settings()
    settings.depmap.cache_dir = tmp_path
    processor = ProteinExpressionMSWideProcessor.__new__(ProteinExpressionMSWideProcessor)
    processor.settings = settings
    return processor


def test_local_gene_lookup_uses_unique_aliases_and_uniprot_ids(tmp_path: Path) -> None:
    gene_cache = tmp_path / "Gene.csv"
    gene_cache.write_text(
        "symbol,alias_symbol,prev_symbol,uniprot_ids\n"
        "EPDR1,MERP-1|UCC1,,Q9UM22\n"
        "C3orf52,FLJ23186|TTMP,,Q5BVD1\n"
        "MFSD5,MGC11308|HsMOT2,SLC61A1,Q6N075\n"
        "SERF1A,H4F5|4F5,FAM2A,O75920\n"
        "SERF1B,FAM2B,,O75920\n",
        encoding="utf-8",
    )

    local_genes = pd.DataFrame(
        [
            {
                "gene_id": "EPDR1",
                "hugo_symbol": "EPDR1",
                "entrez_id": 54749,
                "ensembl_id": None,
            },
            {
                "gene_id": "C3orf52",
                "hugo_symbol": "C3orf52",
                "entrez_id": 79669,
                "ensembl_id": None,
            },
            {
                "gene_id": "MFSD5",
                "hugo_symbol": "MFSD5",
                "entrez_id": 84975,
                "ensembl_id": None,
            },
            {
                "gene_id": "SERF1A",
                "hugo_symbol": "SERF1A",
                "entrez_id": 8293,
                "ensembl_id": None,
            },
            {
                "gene_id": "SERF1B",
                "hugo_symbol": "SERF1B",
                "entrez_id": 728492,
                "ensembl_id": None,
            },
        ]
    )

    processor = _make_processor(tmp_path)
    alias_lookup, uniprot_lookup = processor._build_local_gene_lookups(local_genes)

    assert alias_lookup["ucc1"]["local_gene_id"] == "EPDR1"
    assert alias_lookup["ttmp"]["local_gene_id"] == "C3orf52"
    assert alias_lookup["slc61a1"]["local_gene_id"] == "MFSD5"
    assert uniprot_lookup["Q5BVD1"]["local_gene_id"] == "C3orf52"
    assert uniprot_lookup["Q6N075"]["local_gene_id"] == "MFSD5"
    assert "O75920" not in uniprot_lookup


def test_alias_and_uniprot_resolution_skip_ambiguous_matches(tmp_path: Path) -> None:
    gene_cache = tmp_path / "Gene.csv"
    gene_cache.write_text(
        "symbol,alias_symbol,prev_symbol,uniprot_ids\n"
        "EPDR1,MERP-1|UCC1,,Q9UM22\n"
        "C3orf52,FLJ23186|TTMP,,Q5BVD1\n"
        "MFSD5,MGC11308|HsMOT2,SLC61A1,Q6N075\n"
        "SERF1A,H4F5|4F5,FAM2A,O75920\n"
        "SERF1B,FAM2B,,O75920\n",
        encoding="utf-8",
    )

    local_genes = pd.DataFrame(
        [
            {
                "gene_id": "EPDR1",
                "hugo_symbol": "EPDR1",
                "entrez_id": 54749,
                "ensembl_id": None,
            },
            {
                "gene_id": "C3orf52",
                "hugo_symbol": "C3orf52",
                "entrez_id": 79669,
                "ensembl_id": None,
            },
            {
                "gene_id": "MFSD5",
                "hugo_symbol": "MFSD5",
                "entrez_id": 84975,
                "ensembl_id": None,
            },
            {
                "gene_id": "SERF1A",
                "hugo_symbol": "SERF1A",
                "entrez_id": 8293,
                "ensembl_id": None,
            },
            {
                "gene_id": "SERF1B",
                "hugo_symbol": "SERF1B",
                "entrez_id": 728492,
                "ensembl_id": None,
            },
        ]
    )

    bridge_df = pd.DataFrame(
        [
            {
                "protein_accession": "A4D1W8",
                "protein_accession_base": "A4D1W8",
                "gene_symbol": "UCC1",
                "local_gene_id": None,
                "local_hugo_symbol": None,
                "local_entrez_id": None,
                "local_ensembl_id": None,
                "mapping_method": "exact",
            },
            {
                "protein_accession": "Q5BVD1-3",
                "protein_accession_base": "Q5BVD1",
                "gene_symbol": "TTMP",
                "local_gene_id": None,
                "local_hugo_symbol": None,
                "local_entrez_id": None,
                "local_ensembl_id": None,
                "mapping_method": "exact",
            },
            {
                "protein_accession": "Q6N075-2",
                "protein_accession_base": "Q6N075",
                "gene_symbol": "SLC61A1",
                "local_gene_id": None,
                "local_hugo_symbol": None,
                "local_entrez_id": None,
                "local_ensembl_id": None,
                "mapping_method": "exact",
            },
            {
                "protein_accession": "O75920-2",
                "protein_accession_base": "O75920",
                "gene_symbol": "SERF1A; SERF1B",
                "local_gene_id": None,
                "local_hugo_symbol": None,
                "local_entrez_id": None,
                "local_ensembl_id": None,
                "mapping_method": "exact",
            },
        ]
    )

    processor = _make_processor(tmp_path)
    alias_lookup, uniprot_lookup = processor._build_local_gene_lookups(local_genes)
    processor._resolve_local_genes_from_aliases(bridge_df, alias_lookup)
    processor._resolve_local_genes_from_uniprot_ids(bridge_df, uniprot_lookup)

    assert bridge_df.loc[0, "local_gene_id"] == "EPDR1"
    assert bridge_df.loc[0, "mapping_method"] == "exact+gene_symbol_alias_lookup"
    assert bridge_df.loc[1, "local_gene_id"] == "C3orf52"
    assert bridge_df.loc[2, "local_gene_id"] == "MFSD5"
    assert bridge_df.loc[3, "local_gene_id"] is None
