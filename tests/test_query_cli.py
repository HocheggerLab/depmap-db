"""Tests for the phase-1 gene query CLI layer."""

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
from depmap_db.query import GeneQueryService


@pytest.fixture()
def populated_database(tmp_path: Path) -> Path:
    """Create a small DuckDB fixture with wide dependency/expression data."""
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

    model_rows = [
        (
            "ACH-000001",
            "CellA",
            "Breast",
            "Breast Cancer",
            "SubtypeA",
        ),
        (
            "ACH-000002",
            "CellB",
            "Breast",
            "Breast Cancer",
            "SubtypeB",
        ),
        (
            "ACH-000003",
            "CellC",
            "Lung",
            "Lung Cancer",
            "SubtypeC",
        ),
    ]
    db.execute_many(
        """
        INSERT INTO models (
            model_id,
            cell_line_name,
            oncotree_lineage,
            oncotree_primary_disease,
            oncotree_subtype
        ) VALUES (?, ?, ?, ?, ?)
        """,
        model_rows,
    )

    dependency_rows = [
        ("ACH-000001", -1.2, -0.8),
        ("ACH-000002", -0.8, -0.3),
        ("ACH-000003", -0.2, -1.5),
    ]
    db.execute_many(
        'INSERT INTO gene_effects_wide (model_id, "HAPSTR1", "MASTL") VALUES (?, ?, ?)',
        dependency_rows,
    )

    expression_rows = [
        ("ACH-000001", 1.0, 10.0),
        ("ACH-000002", 3.0, 12.0),
        ("ACH-000003", 8.0, 5.0),
    ]
    db.execute_many(
        'INSERT INTO gene_expression_wide (model_id, "HAPSTR1", "MASTL") VALUES (?, ?, ?)',
        expression_rows,
    )

    return db_path


def test_dependency_summary_groups_by_lineage(
    populated_database: Path,
) -> None:
    """Dependency summaries should aggregate by lineage."""
    service = GeneQueryService()

    result = service.get_dependency_summary("HAPSTR1", group_by="lineage")

    assert list(result["group_name"]) == ["Breast", "Lung"]
    assert list(result["model_count"]) == [2, 1]
    assert result.loc[0, "mean_dependency"] == pytest.approx(-1.0)
    assert result.loc[1, "mean_dependency"] == pytest.approx(-0.2)


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
    assert '"mean_expression":8.0' in result.output


def test_gene_query_raises_for_missing_gene(populated_database: Path) -> None:
    """Missing genes should fail cleanly instead of generating broken SQL."""
    service = GeneQueryService()

    with pytest.raises(ValueError, match="NOT_A_GENE"):
        service.get_dependency_summary("NOT_A_GENE")
