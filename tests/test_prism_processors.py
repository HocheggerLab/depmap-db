"""Tests for phase-1 PRISM drug sensitivity support."""

from pathlib import Path

import pandas as pd
import pytest

from depmap_db.config import reload_settings
from depmap_db.database import create_tables, get_db_manager, reset_db_manager
from depmap_db.etl.processors.prism_primary import PrismPrimaryWideProcessor
from depmap_db.etl.processors.prism_secondary import PrismSecondaryProcessor


@pytest.fixture()
def prism_database(tmp_path: Path) -> Path:
    """Create a small DuckDB fixture for PRISM processor tests."""
    db_path = tmp_path / "depmap-prism-test.duckdb"

    settings = reload_settings()
    settings.database.path = db_path
    settings.database.memory = False
    settings.depmap.cache_dir = tmp_path / "cache"
    settings.depmap.cache_dir.mkdir(parents=True, exist_ok=True)

    reset_db_manager()
    create_tables()

    db = get_db_manager()
    db.execute_many(
        "INSERT INTO models (model_id, cell_line_name) VALUES (?, ?)",
        [
            ("ACH-000001", "CellA"),
            ("ACH-000002", "CellB"),
        ],
    )
    db.execute_many(
        "INSERT INTO genes (gene_id, hugo_symbol, entrez_id) VALUES (?, ?, ?)",
        [
            ("CDK4", "CDK4", 1019),
            ("CDK6", "CDK6", 1021),
            ("PARP1", "PARP1", 142),
            ("PARP2", "PARP2", 10038),
        ],
    )

    gene_cache = settings.depmap.cache_dir / "Gene.csv"
    gene_cache.write_text(
        "symbol,alias_symbol,prev_symbol\n"
        "CDK4,,\n"
        "CDK6,,\n"
        "PARP1,,\n"
        "PARP2,,\n",
        encoding="utf-8",
    )

    return db_path


def test_prism_primary_processor_loads_wide_matrix_and_compound_targets(
    prism_database: Path, tmp_path: Path
) -> None:
    settings = reload_settings()
    settings.database.path = prism_database
    settings.database.memory = False
    settings.depmap.cache_dir = tmp_path / "cache"
    reset_db_manager()

    matrix_path = tmp_path / "Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv"
    matrix_path.write_text(
        ",ACH-000001,ACH-000002\n"
        "BRD-K11111111-001-01-1,-0.8,-0.3\n"
        "BRD-K22222222-001-01-1,-0.2,-1.1\n",
        encoding="utf-8",
    )
    compound_list_path = tmp_path / "Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv"
    compound_list_path.write_text(
        "screen,dose,repurposing_target,MOA,IDs,Drug.Name,Synonyms\n"
        "REP.300,2.5,CDK4,CDK inhibitor,BRD-K11111111-001-01-1,palbociclib,PALBO\n"
        "REP.300,2.5,PARP1,PARP inhibitor,BRD-K22222222-001-01-1,olaparib,OLAP\n",
        encoding="utf-8",
    )

    processor = PrismPrimaryWideProcessor()
    result = processor.process_file(matrix_path)

    assert result.status == "success"
    db = get_db_manager()
    assert db.execute("SELECT COUNT(*) FROM drug_response_primary_wide").fetchone()[0] == 2
    assert db.execute("SELECT COUNT(*) FROM compounds").fetchone()[0] == 2
    assert db.execute("SELECT COUNT(*) FROM compound_targets").fetchone()[0] == 2
    assert db.execute("SELECT compound_name FROM compounds WHERE broad_id = 'BRD-K11111111-001-01-1'").fetchone()[0] == "palbociclib"
    assert db.execute(
        "SELECT local_gene_id FROM compound_targets WHERE broad_id = 'BRD-K11111111-001-01-1'"
    ).fetchone()[0] == "CDK4"


def test_prism_secondary_processor_loads_long_table_and_target_bridge(
    prism_database: Path, tmp_path: Path
) -> None:
    settings = reload_settings()
    settings.database.path = prism_database
    settings.database.memory = False
    settings.depmap.cache_dir = tmp_path / "cache"
    reset_db_manager()

    secondary_path = tmp_path / "secondary-screen-dose-response-curve-parameters.csv"
    secondary_df = pd.DataFrame(
        {
            "broad_id": ["BRD-K33333333-001-01-1", "BRD-K44444444-001-01-1"],
            "depmap_id": ["ACH-000001", "ACH-000002"],
            "ccle_name": ["CELLA", "CELLB"],
            "screen_id": ["MTS010", "MTS010"],
            "upper_limit": [1.0, 1.1],
            "lower_limit": [-1.0, -1.2],
            "slope": [2.0, 2.5],
            "r2": [0.9, 0.95],
            "auc": [0.81, 0.27],
            "ec50": [1.2, 0.9],
            "ic50": [1.5, 0.8],
            "passed_str_profiling": [True, True],
            "row_name": ["PR500_ACH-000001", "PR500_ACH-000002"],
            "name": ["ribociclib", "olaparib"],
            "moa": ["CDK inhibitor", "PARP inhibitor"],
            "target": ["CDK4, CDK6", "PARP1, PARP2"],
            "disease.area": ["oncology", "oncology"],
            "indication": ["breast cancer", "ovarian cancer"],
            "smiles": ["SMILES1", "SMILES2"],
            "phase": ["Launched", "Launched"],
        }
    )
    secondary_df.to_csv(secondary_path, index=False)

    processor = PrismSecondaryProcessor()
    result = processor.process_file(secondary_path)

    assert result.status == "success"
    db = get_db_manager()
    assert db.execute("SELECT COUNT(*) FROM drug_response_secondary").fetchone()[0] == 2
    assert db.execute("SELECT COUNT(*) FROM compounds").fetchone()[0] == 2
    assert db.execute(
        "SELECT COUNT(*) FROM compound_targets WHERE broad_id = 'BRD-K33333333-001-01-1'"
    ).fetchone()[0] == 2
    metric = db.execute(
        "SELECT default_secondary_summary_metric FROM drug_screens WHERE screen_id = 'MTS010'"
    ).fetchone()[0]
    assert metric == "auc"
    mapping_status = db.execute(
        "SELECT mapping_status FROM compound_targets WHERE broad_id = 'BRD-K44444444-001-01-1' AND target_symbol = 'PARP1'"
    ).fetchone()[0]
    assert mapping_status == "resolved_to_local_gene"
