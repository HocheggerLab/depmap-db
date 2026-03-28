"""Tests for the Polars Python API."""

import os
from pathlib import Path

import polars as pl
import pytest

os.environ.setdefault("LOG_LEVEL", "INFO")
os.environ.setdefault("LOG_FILE_PATH", "logs/app.log")

from depmap_db.config import reload_settings
from depmap_db.database import create_tables, get_db_manager, reset_db_manager
from depmap_db.polars import (
    export_polars_tables,
    get_lazy_tables,
    lazy_mutations,
    lazy_table,
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

    db.execute_many(
        """
        INSERT INTO models (
            model_id,
            cell_line_name,
            oncotree_lineage,
            oncotree_primary_disease
        ) VALUES (?, ?, ?, ?)
        """,
        [
            ("ACH-000001", "CellA", "Breast", "Breast Cancer"),
            ("ACH-000002", "CellB", "Lung", "Lung Cancer"),
        ],
    )
    db.execute_many(
        'INSERT INTO genes (gene_id, hugo_symbol) VALUES (?, ?)',
        [("1", "TP53"), ("2", "KRAS")],
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
        tables=["models", "mutations"],
        db_path=polars_database,
    )

    assert exported["models"] == output_dir / "models.parquet"
    assert exported["mutations"] == output_dir / "mutations.parquet"
    assert output_dir.joinpath("models.parquet").exists()
    assert output_dir.joinpath("mutations.parquet").exists()
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
