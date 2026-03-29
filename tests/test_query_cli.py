"""Tests for the phase-1 query CLI layer."""

import os
from pathlib import Path

import pandas as pd
import pytest
from click.testing import CliRunner

os.environ.setdefault("LOG_LEVEL", "INFO")
os.environ.setdefault("LOG_FILE_PATH", "logs/app.log")

from depmap_db.cli import cli
from depmap_db.config import reload_settings
from depmap_db.database import create_tables, get_db_manager, reset_db_manager
from depmap_db.query import (
    GeneQueryService,
    MutationQueryService,
    ProteinQueryService,
)


@pytest.fixture()
def populated_database(tmp_path: Path) -> Path:
    """Create a small DuckDB fixture with dependency/expression/mutation data."""
    db_path = tmp_path / "depmap-test.duckdb"

    settings = reload_settings()
    settings.database.path = db_path
    settings.database.memory = False
    reset_db_manager()
    create_tables()

    db = get_db_manager()
    db.execute('ALTER TABLE gene_effects_wide ADD COLUMN "HAPSTR1" DOUBLE')
    db.execute('ALTER TABLE gene_effects_wide ADD COLUMN "MASTL" DOUBLE')
    db.execute('ALTER TABLE gene_expression_wide ADD COLUMN "HAPSTR1" DOUBLE')
    db.execute('ALTER TABLE gene_expression_wide ADD COLUMN "MASTL" DOUBLE')
    db.execute('ALTER TABLE protein_expression_ms_wide ADD COLUMN "protein_A0AV96" DOUBLE')
    db.execute('ALTER TABLE protein_expression_ms_wide ADD COLUMN "protein_Q99729_3" DOUBLE')

    model_rows = [
        (
            "ACH-000001",
            "CellA",
            "CELLA",
            "CELLA_CCLE",
            "Breast",
            "Breast Cancer",
            "SubtypeA",
        ),
        (
            "ACH-000002",
            "CellB",
            "CELLB",
            "CELLB_CCLE",
            "Breast",
            "Breast Cancer",
            "SubtypeB",
        ),
        (
            "ACH-000003",
            "CellC",
            "CELLC",
            "CELLC_CCLE",
            "Lung",
            "Lung Cancer",
            "SubtypeC",
        ),
        (
            "ACH-000004",
            "CellD",
            "CELLD",
            "CELLD_CCLE",
            "Lung",
            "Lung Cancer",
            "SubtypeD",
        ),
    ]
    db.execute_many(
        """
        INSERT INTO models (
            model_id,
            cell_line_name,
            stripped_cell_line_name,
            ccle_name,
            oncotree_lineage,
            oncotree_primary_disease,
            oncotree_subtype
        ) VALUES (?, ?, ?, ?, ?, ?, ?)
        """,
        model_rows,
    )

    gene_rows = [
        ("RBM47", "RBM47", 54502),
        ("HNRNPAB", "HNRNPAB", 3182),
    ]
    db.execute_many(
        "INSERT INTO genes (gene_id, hugo_symbol, entrez_id) VALUES (?, ?, ?)",
        gene_rows,
    )

    dependency_rows = [
        ("ACH-000001", -1.2, -0.8),
        ("ACH-000002", -0.8, -0.3),
        ("ACH-000003", -0.2, -1.5),
        ("ACH-000004", -0.4, -0.7),
    ]
    db.execute_many(
        'INSERT INTO gene_effects_wide (model_id, "HAPSTR1", "MASTL") VALUES (?, ?, ?)',
        dependency_rows,
    )

    expression_rows = [
        ("ACH-000001", 1.0, 10.0),
        ("ACH-000002", 3.0, 12.0),
        ("ACH-000003", 8.0, 5.0),
        ("ACH-000004", 5.0, 6.0),
    ]
    db.execute_many(
        'INSERT INTO gene_expression_wide (model_id, "HAPSTR1", "MASTL") VALUES (?, ?, ?)',
        expression_rows,
    )

    protein_feature_rows = [
        (
            "A0AV96",
            "A0AV96",
            "protein_A0AV96",
            "RBM47_HUMAN",
            "RNA-binding protein 47",
            "RBM47",
            54502,
            "ENST00000295971.12;ENST00000381793.6",
            "RBM47",
            "RBM47",
            54502,
            None,
            True,
            "exact",
            "mapped_to_local_gene",
            "ProteomicsMSGygi",
            "harmonized_MS_CCLE_Gygi.csv",
            "proteomics_ms",
            "Harmonized Public Proteomics 24Q4",
        ),
        (
            "Q99729-3",
            "Q99729",
            "protein_Q99729_3",
            "ROAA_HUMAN",
            "Isoform 3 of Heterogeneous nuclear ribonucleoprotein A/B",
            "HNRNPAB",
            None,
            "ENST00000355836.9;ENST00000506259.5",
            "HNRNPAB",
            "HNRNPAB",
            3182,
            None,
            True,
            "base_accession_fallback",
            "mapped_to_local_gene",
            "ProteomicsMSGygi",
            "harmonized_MS_CCLE_Gygi.csv",
            "proteomics_ms",
            "Harmonized Public Proteomics 24Q4",
        ),
    ]
    db.execute_many(
        """
        INSERT INTO protein_features (
            protein_accession,
            protein_accession_base,
            storage_column_name,
            protein_entry_name,
            protein_name,
            gene_symbol,
            entrez_id,
            ensembl_transcript_ids,
            local_gene_id,
            local_hugo_symbol,
            local_entrez_id,
            local_ensembl_id,
            is_reviewed,
            mapping_method,
            mapping_status,
            source_dataset,
            source_filename,
            modality,
            release_label
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        protein_feature_rows,
    )

    protein_rows = [
        ("ACH-000001", 0.3, 1.2),
        ("ACH-000002", 0.1, 1.0),
        ("ACH-000003", -0.4, -0.5),
        ("ACH-000004", -0.2, -0.7),
    ]
    db.execute_many(
        'INSERT INTO protein_expression_ms_wide (model_id, "protein_A0AV96", "protein_Q99729_3") VALUES (?, ?, ?)',
        protein_rows,
    )

    mutation_rows = [
        (
            "mut-1",
            "ACH-000001",
            "17",
            7577538,
            "C",
            "T",
            "SNP",
            "c.916C>T",
            "p.R306*",
            "TP53",
            "stop_gained",
            "HIGH",
            0.45,
            True,
            False,
            True,
        ),
        (
            "mut-2",
            "ACH-000001",
            "12",
            25398284,
            "C",
            "T",
            "SNP",
            "c.35G>A",
            "p.G12D",
            "KRAS",
            "missense_variant",
            "MODERATE",
            0.31,
            False,
            True,
            True,
        ),
        (
            "mut-3",
            "ACH-000002",
            "17",
            7673803,
            "C",
            "T",
            "SNP",
            "c.818G>A",
            "p.R273H",
            "TP53",
            "missense_variant",
            "MODERATE",
            0.52,
            False,
            True,
            False,
        ),
        (
            "mut-4",
            "ACH-000003",
            "7",
            55249071,
            "G",
            "A",
            "SNP",
            "c.2573T>G",
            "p.L858R",
            "EGFR",
            "missense_variant",
            "MODERATE",
            0.41,
            False,
            True,
            True,
        ),
        (
            "mut-5",
            "ACH-000004",
            "17",
            7579472,
            "G",
            "A",
            "SNP",
            "c.743G>A",
            "p.R248Q",
            "TP53",
            "missense_variant",
            "MODERATE",
            0.33,
            False,
            True,
            True,
        ),
    ]
    db.execute_many(
        """
        INSERT INTO mutations (
            mutation_id,
            model_id,
            chrom,
            pos,
            ref,
            alt,
            variant_type,
            dna_change,
            protein_change,
            hugo_symbol,
            molecular_consequence,
            vep_impact,
            af,
            likely_lof,
            hotspot,
            hess_driver
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        mutation_rows,
    )

    status_rows = [
        ("ACH-000001", "TP53", True, 1, True, False, False, True, True),
        ("ACH-000001", "KRAS", True, 1, False, True, True, False, True),
        ("ACH-000002", "TP53", True, 1, False, True, False, False, False),
        ("ACH-000003", "EGFR", True, 1, False, True, True, False, True),
        ("ACH-000004", "TP53", True, 1, False, True, False, False, True),
    ]
    db.execute_many(
        """
        INSERT INTO model_gene_mutation_status (
            model_id,
            gene_symbol,
            is_mutated,
            mutation_count,
            has_likely_lof,
            has_hotspot,
            has_oncogene_high_impact,
            has_tumor_suppressor_high_impact,
            has_driver
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        status_rows,
    )

    return db_path


def test_dependency_summary_groups_by_lineage(
    populated_database: Path,
) -> None:
    """Dependency summaries should aggregate by lineage."""
    service = GeneQueryService()

    result = service.get_dependency_summary("HAPSTR1", group_by="lineage")

    assert list(result["group_name"]) == ["Breast", "Lung"]
    assert list(result["model_count"]) == [2, 2]
    assert result.loc[0, "mean_dependency"] == pytest.approx(-1.0)
    assert result.loc[1, "mean_dependency"] == pytest.approx(-0.3)


def test_dependency_models_can_filter_by_lineage(
    populated_database: Path,
) -> None:
    """Per-model dependency lookups should support lineage filters."""
    service = GeneQueryService()

    result = service.get_dependency_models("MASTL", lineage="Breast")

    assert list(result["model_id"]) == ["ACH-000001", "ACH-000002"]
    assert list(result["dependency"]) == [-0.8, -0.3]


def test_gene_cli_dependency_models_can_export_csv(
    populated_database: Path, tmp_path: Path
) -> None:
    """The CLI should support first-class CSV export for per-model dependency data."""
    output_path = tmp_path / "mastl_breast.csv"
    runner = CliRunner()

    result = runner.invoke(
        cli,
        [
            "gene",
            "dependency-models",
            "MASTL",
            "--lineage",
            "Breast",
            "--output",
            str(output_path),
        ],
    )

    assert result.exit_code == 0
    assert output_path.exists()

    exported = pd.read_csv(output_path)
    assert list(exported["model_id"]) == ["ACH-000001", "ACH-000002"]
    assert list(exported["dependency"]) == [-0.8, -0.3]
    assert "Exported 2 dependency rows" in result.output


def test_gene_cli_expression_summary_supports_json(
    populated_database: Path,
) -> None:
    """The CLI should expose expression summaries for downstream tooling."""
    runner = CliRunner()

    result = runner.invoke(
        cli,
        [
            "gene",
            "expression-summary",
            "HAPSTR1",
            "--group-by",
            "lineage",
            "--format",
            "json",
            "--limit",
            "5",
        ],
    )

    assert result.exit_code == 0
    assert '"group_name":"Lung"' in result.output
    assert '"mean_expression":6.5' in result.output


def test_gene_query_raises_for_missing_gene(populated_database: Path) -> None:
    """Missing genes should fail cleanly instead of generating broken SQL."""
    service = GeneQueryService()

    with pytest.raises(ValueError, match="NOT_A_GENE"):
        service.get_dependency_summary("NOT_A_GENE")


def test_model_mutations_resolves_cell_line_name_and_filters(
    populated_database: Path,
) -> None:
    """Model mutation lookups should accept cell line names and event filters."""
    service = MutationQueryService()

    result = service.get_model_mutations("cella", hotspot_only=True)

    assert list(result["gene"]) == ["KRAS"]
    assert bool(result.loc[0, "driver"]) is True
    assert result.loc[0, "matched_cell_line"] == "CellA"


def test_lineage_mutation_frequency_ranks_genes(
    populated_database: Path,
) -> None:
    """Lineage prevalence should rank genes by mutated-model count."""
    service = MutationQueryService()

    result = service.get_lineage_mutation_frequency("Breast")

    assert list(result["gene_symbol"]) == ["TP53", "KRAS"]
    assert list(result["mutated_model_count"]) == [2, 1]
    assert list(result["total_models"]) == [2, 2]
    assert result.loc[0, "mutation_frequency"] == pytest.approx(1.0)
    assert result.loc[1, "mutation_frequency"] == pytest.approx(0.5)


def test_dependency_by_mutation_compares_mutant_vs_wt(
    populated_database: Path,
) -> None:
    """Dependency-by-mutation should aggregate mutant and wild-type cohorts."""
    service = GeneQueryService()

    result = service.get_dependency_by_mutation(
        "HAPSTR1",
        mutation_gene="TP53",
        mutation_class="driver",
    )

    row = result.iloc[0]
    assert row["mutant_model_count"] == 2
    assert row["wt_model_count"] == 2
    assert row["mutant_mean_dependency"] == pytest.approx(-0.8)
    assert row["wt_mean_dependency"] == pytest.approx(-0.5)
    assert row["delta_mean"] == pytest.approx(-0.3)


def test_model_mutations_cli_supports_json_output(
    populated_database: Path,
) -> None:
    """The model mutations CLI should expose JSON for downstream scripting."""
    runner = CliRunner()

    result = runner.invoke(
        cli,
        [
            "model",
            "mutations",
            "ACH-000001",
            "--gene",
            "TP53",
            "--format",
            "json",
        ],
    )

    assert result.exit_code == 0
    assert '"gene":"TP53"' in result.output
    assert '"protein_change":"p.R306*"' in result.output


def test_lineage_mutation_frequency_cli_supports_json(
    populated_database: Path,
) -> None:
    """The lineage mutation-frequency CLI should expose JSON."""
    runner = CliRunner()

    result = runner.invoke(
        cli,
        [
            "lineage",
            "mutation-frequency",
            "Lung",
            "--format",
            "json",
            "--limit",
            "5",
        ],
    )

    assert result.exit_code == 0
    assert '"gene_symbol":"EGFR"' in result.output
    assert '"gene_symbol":"TP53"' in result.output


def test_dependency_by_mutation_cli_supports_json(
    populated_database: Path,
) -> None:
    """The gene dependency-by-mutation CLI should expose JSON."""
    runner = CliRunner()

    result = runner.invoke(
        cli,
        [
            "gene",
            "dependency-by-mutation",
            "HAPSTR1",
            "--mutation-gene",
            "TP53",
            "--mutation-class",
            "driver",
            "--format",
            "json",
        ],
    )

    assert result.exit_code == 0
    assert '"mutation_gene":"TP53"' in result.output
    assert '"mutant_model_count":2' in result.output


def test_protein_expression_summary_supports_gene_resolution(
    populated_database: Path,
) -> None:
    """Protein summaries should resolve a unique mapped gene symbol."""
    service = ProteinQueryService()

    result = service.get_expression_summary("RBM47", group_by="lineage")

    assert list(result["group_name"]) == ["Breast", "Lung"]
    assert result.loc[0, "mean_abundance"] == pytest.approx(0.2)
    assert result.loc[1, "mean_abundance"] == pytest.approx(-0.3)


def test_protein_mapping_summary_cli_supports_json(
    populated_database: Path,
) -> None:
    """The protein CLI should expose bridge coverage summaries."""
    runner = CliRunner()

    result = runner.invoke(
        cli,
        ["protein", "mapping-summary", "--format", "json"],
    )

    assert result.exit_code == 0
    assert '"source_dataset":"ProteomicsMSGygi"' in result.output
    assert '"protein_count":2' in result.output


def test_polars_is_exposed_as_a_cli_group() -> None:
    """Polars access should also be reachable from the CLI export surface."""
    runner = CliRunner()

    result = runner.invoke(cli, ["--help"])

    assert result.exit_code == 0
    assert "polars" in result.output


if __name__ == "__main__":
    pytest.main([__file__])
