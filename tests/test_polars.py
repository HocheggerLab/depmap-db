"""Tests for the Polars Python API."""

import os
from pathlib import Path

import polars as pl
import pytest
from click.testing import CliRunner

os.environ.setdefault("LOG_LEVEL", "INFO")
os.environ.setdefault("LOG_FILE_PATH", "logs/app.log")

from depmap_db.cli import cli
from depmap_db.config import reload_settings
from depmap_db.database import create_tables, get_db_manager, reset_db_manager
from depmap_db.polars import (
    export_polars_datasets,
    export_polars_tables,
    get_lazy_datasets,
    get_lazy_tables,
    lazy_dataset,
    lazy_mutations,
    lazy_table,
    prepare_lazy_datasets,
    prepare_lazy_tables,
)


@pytest.fixture()
def polars_database(tmp_path: Path) -> Path:
    """Create a small DuckDB fixture for Polars export tests."""
    db_path = tmp_path / "depmap-polars.duckdb"

    settings = reload_settings()
    settings.database.path = db_path
    settings.database.memory = False
    reset_db_manager()
    create_tables()

    db = get_db_manager()
    db.execute('ALTER TABLE gene_effects_wide ADD COLUMN "TP53" DOUBLE')
    db.execute('ALTER TABLE gene_expression_wide ADD COLUMN "TP53" DOUBLE')
    db.execute(
        'ALTER TABLE protein_expression_ms_wide ADD COLUMN "protein_A0AV96" DOUBLE'
    )
    db.execute(
        'ALTER TABLE drug_response_primary_wide ADD COLUMN "ACH-000001" DOUBLE'
    )
    db.execute(
        'ALTER TABLE drug_response_primary_wide ADD COLUMN "ACH-000002" DOUBLE'
    )

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
        [
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
                "Lung",
                "Lung Cancer",
                "SubtypeB",
            ),
        ],
    )
    db.execute_many(
        'INSERT INTO genes (gene_id, hugo_symbol) VALUES (?, ?)',
        [("1", "TP53"), ("2", "RBM47"), ("3", "CDK4")],
    )
    db.execute_many(
        'INSERT INTO gene_effects_wide (model_id, "TP53") VALUES (?, ?)',
        [("ACH-000001", -1.0), ("ACH-000002", -0.2)],
    )
    db.execute_many(
        'INSERT INTO gene_expression_wide (model_id, "TP53") VALUES (?, ?)',
        [("ACH-000001", 8.0), ("ACH-000002", 5.0)],
    )
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
            local_gene_id,
            local_hugo_symbol,
            local_entrez_id,
            is_reviewed,
            mapping_method,
            mapping_status,
            source_dataset,
            source_filename,
            modality,
            release_label
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        [
            (
                "A0AV96",
                "A0AV96",
                "protein_A0AV96",
                "RBM47_HUMAN",
                "RNA-binding protein 47",
                "RBM47",
                54502,
                "2",
                "RBM47",
                54502,
                True,
                "exact",
                "mapped_to_local_gene",
                "ProteomicsMSGygi",
                "harmonized_MS_CCLE_Gygi.csv",
                "proteomics_ms",
                "Harmonized Public Proteomics 24Q4",
            )
        ],
    )
    db.execute_many(
        'INSERT INTO protein_expression_ms_wide (model_id, "protein_A0AV96") VALUES (?, ?)',
        [("ACH-000001", 0.3), ("ACH-000002", 0.1)],
    )
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
            protein_change,
            hugo_symbol,
            molecular_consequence,
            vep_impact,
            likely_lof,
            hotspot,
            hess_driver
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        [
            (
                "mut-1",
                "ACH-000001",
                "17",
                7579472,
                "G",
                "A",
                "SNP",
                "p.R248Q",
                "TP53",
                "missense_variant",
                "MODERATE",
                False,
                True,
                True,
            ),
        ],
    )
    db.execute_many(
        """
        INSERT INTO model_gene_mutation_status (
            model_id,
            gene_symbol,
            is_mutated,
            mutation_count,
            has_likely_lof,
            has_hotspot,
            has_driver
        ) VALUES (?, ?, ?, ?, ?, ?, ?)
        """,
        [("ACH-000001", "TP53", True, 1, False, True, True)],
    )

    db.execute(
        """
        INSERT INTO compounds (
            broad_id,
            compound_name,
            compound_synonyms,
            moa,
            target_text,
            source_dataset,
            source_filename,
            release_label
        ) VALUES (
            'BRD-K11111111-001-01-1',
            'palbociclib',
            'PALBO',
            'CDK inhibitor',
            'CDK4',
            'PRISMPrimaryRepurposingExtended',
            'Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv',
            'Repurposing Public 24Q2'
        )
        """
    )
    db.execute(
        """
        INSERT INTO drug_screens (
            screen_id,
            screen_kind,
            dataset_name,
            source_filename,
            release_label,
            release_track,
            default_secondary_summary_metric
        ) VALUES (
            'REP.300',
            'primary',
            'PRISMPrimaryRepurposingExtended',
            'Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv',
            'Repurposing Public 24Q2',
            'prism_primary',
            'auc'
        )
        """
    )
    db.execute(
        """
        INSERT INTO drug_screens (
            screen_id,
            screen_kind,
            dataset_name,
            source_filename,
            release_label,
            release_track,
            default_secondary_summary_metric
        ) VALUES (
            'MTS010',
            'secondary',
            'PRISMSecondaryDoseResponseCurveParameters',
            'secondary-screen-dose-response-curve-parameters.csv',
            'Repurposing Public 24Q2',
            'prism_secondary',
            'auc'
        )
        """
    )
    db.execute(
        """
        INSERT INTO compound_targets (
            broad_id,
            target_ordinal,
            target_symbol,
            source_target_text,
            source_dataset,
            source_filename,
            local_gene_id,
            local_hugo_symbol,
            local_entrez_id,
            mapping_status,
            mapping_method
        ) VALUES (
            'BRD-K11111111-001-01-1',
            1,
            'CDK4',
            'CDK4',
            'PRISMPrimaryRepurposingExtended',
            'Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv',
            '3',
            'CDK4',
            1019,
            'resolved_to_local_gene',
            'exact_symbol'
        )
        """
    )
    db.execute_many(
        'INSERT INTO drug_response_primary_wide (broad_id, screen_id, "ACH-000001", "ACH-000002") VALUES (?, ?, ?, ?)',
        [("BRD-K11111111-001-01-1", "REP.300", -0.8, -0.3)],
    )
    db.execute(
        """
        INSERT INTO drug_response_secondary (
            response_id,
            broad_id,
            model_id,
            ccle_name,
            screen_id,
            upper_limit,
            lower_limit,
            slope,
            r2,
            auc,
            ec50,
            ic50,
            compound_name,
            moa,
            target_text,
            source_dataset,
            source_filename
        ) VALUES (
            'BRD-K11111111-001-01-1::ACH-000001::MTS010',
            'BRD-K11111111-001-01-1',
            'ACH-000001',
            'CELLA_CCLE',
            'MTS010',
            1.0,
            -1.0,
            2.0,
            0.9,
            0.81,
            1.2,
            1.5,
            'palbociclib',
            'CDK inhibitor',
            'CDK4',
            'PRISMSecondaryDoseResponseCurveParameters',
            'secondary-screen-dose-response-curve-parameters.csv'
        )
        """
    )

    return db_path


def test_prepare_lazy_tables_returns_lazyframes(
    polars_database: Path, tmp_path: Path
) -> None:
    """Preparing tables should export snapshots and return LazyFrames."""
    output_dir = tmp_path / "polars"

    tables = prepare_lazy_tables(
        output_dir=output_dir,
        tables=["models", "mutations", "model_gene_mutation_status"],
        db_path=polars_database,
    )

    assert isinstance(tables["models"], pl.LazyFrame)
    assert output_dir.joinpath("models.parquet").exists()
    assert tables["models"].select(pl.len()).collect().item() == 2
    assert tables["mutations"].select(pl.len()).collect().item() == 1


def test_export_polars_tables_writes_requested_parquet_files(
    polars_database: Path, tmp_path: Path
) -> None:
    """The Python API should export just the requested snapshot files."""
    output_dir = tmp_path / "lazy"

    exported = export_polars_tables(
        output_dir=output_dir,
        tables=["models", "mutations", "drug_response_secondary"],
        db_path=polars_database,
    )

    assert exported["models"] == output_dir / "models.parquet"
    assert exported["mutations"] == output_dir / "mutations.parquet"
    assert exported["drug_response_secondary"] == (
        output_dir / "drug_response_secondary.parquet"
    )
    assert output_dir.joinpath("models.parquet").exists()
    assert output_dir.joinpath("mutations.parquet").exists()
    assert output_dir.joinpath("drug_response_secondary.parquet").exists()
    assert not output_dir.joinpath("genes.parquet").exists()


def test_get_lazy_tables_reopens_existing_snapshots(
    polars_database: Path, tmp_path: Path
) -> None:
    """Existing snapshots should be reopenable without refresh."""
    output_dir = tmp_path / "polars"
    export_polars_tables(
        output_dir=output_dir,
        tables=["models", "mutations"],
        db_path=polars_database,
    )

    tables = get_lazy_tables(output_dir=output_dir, tables=["models", "mutations"])

    assert tables["models"].select(pl.len()).collect().item() == 2
    assert tables["mutations"].select(pl.len()).collect().item() == 1


def test_get_lazy_tables_requires_existing_snapshots(tmp_path: Path) -> None:
    """Missing snapshots should raise a Python-API-oriented error."""
    with pytest.raises(ValueError, match="prepare_lazy_tables"):
        get_lazy_tables(output_dir=tmp_path / "missing", tables=["models"])



def test_prepare_lazy_datasets_returns_curated_long_frames(
    polars_database: Path, tmp_path: Path
) -> None:
    """Curated datasets should export notebook-friendly long frames."""
    output_dir = tmp_path / "datasets"

    datasets = prepare_lazy_datasets(
        output_dir=output_dir,
        datasets=[
            "mutation_events",
            "proteomics_long",
            "drug_primary_long",
            "drug_secondary_enriched",
        ],
        db_path=polars_database,
    )

    proteomics = datasets["proteomics_long"].collect()
    primary = datasets["drug_primary_long"].collect()
    secondary = datasets["drug_secondary_enriched"].collect()
    mutations = datasets["mutation_events"].collect()

    assert proteomics.height == 2
    assert proteomics.item(0, "protein_accession") == "A0AV96"
    assert set(primary.get_column("model_id").to_list()) == {
        "ACH-000001",
        "ACH-000002",
    }
    assert primary.item(0, "compound_name") == "palbociclib"
    assert secondary.item(0, "default_secondary_summary_metric") == "auc"
    assert mutations.item(0, "cell_line_name") == "CellA"



def test_get_lazy_datasets_reopens_existing_snapshots(
    polars_database: Path, tmp_path: Path
) -> None:
    """Existing curated dataset snapshots should be reopenable."""
    output_dir = tmp_path / "datasets"
    export_polars_datasets(
        output_dir=output_dir,
        datasets=["proteomics_long", "drug_primary_long"],
        db_path=polars_database,
    )

    datasets = get_lazy_datasets(
        output_dir=output_dir,
        datasets=["proteomics_long", "drug_primary_long"],
    )

    assert datasets["proteomics_long"].select(pl.len()).collect().item() == 2
    assert datasets["drug_primary_long"].select(pl.len()).collect().item() == 2



def test_get_lazy_datasets_requires_existing_snapshots(tmp_path: Path) -> None:
    """Missing curated snapshots should raise a dataset-oriented error."""
    with pytest.raises(ValueError, match="prepare_lazy_datasets"):
        get_lazy_datasets(
            output_dir=tmp_path / "missing",
            datasets=["proteomics_long"],
        )



def test_lazy_table_helpers_support_exploratory_work(
    polars_database: Path, tmp_path: Path
) -> None:
    """Single-table helpers should return ready-to-query LazyFrames."""
    output_dir = tmp_path / "polars"

    models = lazy_table("models", output_dir=output_dir, db_path=polars_database)
    mutations = lazy_mutations(output_dir=output_dir, db_path=polars_database)

    joined = mutations.join(
        models.select(["model_id", "cell_line_name"]),
        on="model_id",
        how="left",
    )

    result = joined.collect()
    assert result.height == 1
    assert result.item(0, "cell_line_name") == "CellA"



def test_lazy_dataset_helper_supports_curated_analysis(
    polars_database: Path, tmp_path: Path
) -> None:
    """Curated dataset helpers should behave like the table helpers."""
    output_dir = tmp_path / "datasets"

    proteomics = lazy_dataset(
        "proteomics_long",
        output_dir=output_dir,
        db_path=polars_database,
    )

    result = proteomics.filter(pl.col("cell_line_name") == "CellA").collect()
    assert result.height == 1
    assert result.item(0, "gene_symbol") == "RBM47"



def test_polars_export_cli_supports_tables_and_datasets(
    polars_database: Path, tmp_path: Path
) -> None:
    """The CLI should export both raw tables and curated datasets."""
    output_dir = tmp_path / "cli-polars"
    runner = CliRunner()

    result = runner.invoke(
        cli,
        [
            "polars",
            "export",
            "--table",
            "models",
            "--dataset",
            "proteomics_long",
            "--output-dir",
            str(output_dir),
            "--database-path",
            str(polars_database),
        ],
    )

    assert result.exit_code == 0
    assert output_dir.joinpath("models.parquet").exists()
    assert output_dir.joinpath("proteomics_long.parquet").exists()
    assert "Polars Export Results" in result.output
