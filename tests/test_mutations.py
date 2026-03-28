"""Tests for mutation support."""

import tempfile
from collections.abc import Generator
from pathlib import Path

import pandas as pd
import pytest

from depmap_db.config import reload_settings
from depmap_db.database import (
    create_tables,
    get_current_schema_version,
    get_db_manager,
    reset_db_manager,
)
from depmap_db.etl.processors import MutationsProcessor


@pytest.fixture
def setup_test_db() -> Generator[None]:
    """Set up a clean test database."""
    settings = reload_settings()
    settings.database.memory = True

    reset_db_manager()
    create_tables()

    # Also create a minimal models table entry for FK constraint
    db_manager = get_db_manager()
    db_manager.execute(
        """
        INSERT INTO models (model_id, cell_line_name)
        VALUES ('ACH-000001', 'TEST_CELL_LINE_1'),
               ('ACH-000002', 'TEST_CELL_LINE_2')
    """
    )

    yield

    reset_db_manager()


def test_mutations_schema_exists(setup_test_db: None) -> None:
    """Test that mutations tables are created in schema."""
    schema_version = get_current_schema_version()
    assert schema_version == "1.3.0"

    db_manager = get_db_manager()
    assert db_manager.table_exists("mutations")
    assert db_manager.table_exists("model_gene_mutation_status")


def test_mutations_processor_basic(setup_test_db: None) -> None:
    """Test basic mutation processing."""
    # Create a test mutation file
    with tempfile.TemporaryDirectory() as tmp_dir:
        test_file = Path(tmp_dir) / "test_mutations.csv"

        # Create test data
        test_data = pd.DataFrame(
            {
                "ModelID": ["ACH-000001", "ACH-000001", "ACH-000002"],
                "Chrom": ["1", "1", "17"],
                "Pos": [1000000, 2000000, 7577538],
                "Ref": ["A", "G", "C"],
                "Alt": ["T", "A", "T"],
                "VariantType": ["SNP", "SNP", "SNP"],
                "DNAChange": ["c.100A>T", "c.200G>A", "c.916C>T"],
                "ProteinChange": ["p.K34N", "p.R67H", "p.R306*"],
                "HugoSymbol": ["GENE1", "GENE2", "TP53"],
                "EnsemblGeneID": [
                    "ENSG00000000001",
                    "ENSG00000000002",
                    "ENSG00000141510",
                ],
                "EntrezGeneID": ["1", "2", "7157"],
                "HgncName": ["GENE1", "GENE2", "TP53"],
                "MolecularConsequence": [
                    "missense_variant",
                    "missense_variant",
                    "stop_gained",
                ],
                "VepImpact": ["MODERATE", "MODERATE", "HIGH"],
                "AF": [0.5, 0.6, 0.45],
                "DP": [100, 150, 80],
                "RefCount": [50, 60, 44],
                "AltCount": [50, 90, 36],
                "GT": ["0/1", "0/1", "0/1"],
                "LikelyLoF": [False, False, True],
                "Hotspot": [False, True, False],
                "OncogeneHighImpact": [False, False, False],
                "TumorSuppressorHighImpact": [False, False, True],
                "HessDriver": [False, False, True],
            }
        )

        test_data.to_csv(test_file, index=False)

        # Process the file
        processor = MutationsProcessor()
        result = processor.process_file(test_file)

        assert result.status == "success"
        assert result.records_inserted == 3

        # Check mutations table
        db_manager = get_db_manager()
        mutations_result = db_manager.execute(
            "SELECT COUNT(*) FROM mutations"
        )
        assert mutations_result.fetchone()[0] == 3

        # Check specific mutation fields
        mutation_check = db_manager.execute(
            """
            SELECT hugo_symbol, likely_lof, hotspot,
                   tumor_suppressor_high_impact, hess_driver
            FROM mutations
            WHERE hugo_symbol = 'TP53'
        """
        )
        row = mutation_check.fetchone()
        assert row[0] == "TP53"
        assert row[1] is True  # likely_lof
        assert row[2] is False  # hotspot
        assert row[3] is True  # tumor_suppressor_high_impact
        assert row[4] is True  # hess_driver

        # Check derived model_gene_mutation_status table
        status_result = db_manager.execute(
            "SELECT COUNT(*) FROM model_gene_mutation_status"
        )
        # Should have 3 unique model-gene pairs
        assert status_result.fetchone()[0] == 3

        # Check aggregated flags for TP53 in ACH-000002
        status_check = db_manager.execute(
            """
            SELECT mutation_count, has_likely_lof, has_hotspot,
                   has_tumor_suppressor_high_impact, has_driver
            FROM model_gene_mutation_status
            WHERE model_id = 'ACH-000002' AND gene_symbol = 'TP53'
        """
        )
        row = status_check.fetchone()
        assert row[0] == 1  # mutation_count
        assert row[1] is True  # has_likely_lof
        assert row[2] is False  # has_hotspot
        assert row[3] is True  # has_tumor_suppressor_high_impact
        assert row[4] is True  # has_driver

        # Check GENE2 with hotspot
        status_check2 = db_manager.execute(
            """
            SELECT has_hotspot
            FROM model_gene_mutation_status
            WHERE model_id = 'ACH-000001' AND gene_symbol = 'GENE2'
        """
        )
        row2 = status_check2.fetchone()
        assert row2[0] is True  # has_hotspot


def test_mutations_processor_aggregation(setup_test_db: None) -> None:
    """Test that multiple mutations for same model-gene are aggregated."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        test_file = Path(tmp_dir) / "test_mutations.csv"

        # Create test data with multiple mutations for same model-gene
        test_data = pd.DataFrame(
            {
                "ModelID": ["ACH-000001", "ACH-000001", "ACH-000001"],
                "Chrom": ["17", "17", "17"],
                "Pos": [7577538, 7577539, 7577540],
                "Ref": ["C", "G", "A"],
                "Alt": ["T", "A", "G"],
                "HugoSymbol": ["TP53", "TP53", "TP53"],
                "LikelyLoF": [True, False, False],
                "Hotspot": [False, True, False],
                "HessDriver": [True, False, True],
            }
        )

        test_data.to_csv(test_file, index=False)

        processor = MutationsProcessor()
        result = processor.process_file(test_file)

        assert result.status == "success"
        assert result.records_inserted == 3

        # Check that only one model-gene status record is created
        db_manager = get_db_manager()
        status_result = db_manager.execute(
            """
            SELECT mutation_count, has_likely_lof, has_hotspot, has_driver
            FROM model_gene_mutation_status
            WHERE model_id = 'ACH-000001' AND gene_symbol = 'TP53'
        """
        )
        row = status_result.fetchone()
        assert row[0] == 3  # mutation_count
        assert row[1] is True  # has_likely_lof (any mutation)
        assert row[2] is True  # has_hotspot (any mutation)
        assert row[3] is True  # has_driver (any mutation)


def test_mutations_processor_missing_gene_symbols(setup_test_db: None) -> None:
    """Test that mutations without gene symbols are still loaded but not in status table."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        test_file = Path(tmp_dir) / "test_mutations.csv"

        test_data = pd.DataFrame(
            {
                "ModelID": ["ACH-000001", "ACH-000001"],
                "Chrom": ["1", "2"],
                "Pos": [1000000, 2000000],
                "Ref": ["A", "G"],
                "Alt": ["T", "C"],
                "HugoSymbol": ["GENE1", None],  # Second mutation has no gene
            }
        )

        test_data.to_csv(test_file, index=False)

        processor = MutationsProcessor()
        result = processor.process_file(test_file)

        assert result.status == "success"
        assert result.records_inserted == 2

        db_manager = get_db_manager()

        # Both mutations should be in mutations table
        mutations_result = db_manager.execute(
            "SELECT COUNT(*) FROM mutations"
        )
        assert mutations_result.fetchone()[0] == 2

        # Only one should be in status table (the one with a gene symbol)
        status_result = db_manager.execute(
            "SELECT COUNT(*) FROM model_gene_mutation_status"
        )
        count = status_result.fetchone()[0]

        # Debug: check what's in the status table
        if count != 1:
            debug_result = db_manager.execute(
                "SELECT model_id, gene_symbol FROM model_gene_mutation_status"
            )
            print(f"Status records: {debug_result.fetchall()}")

        assert count == 1


def test_mutations_force_reload_replaces_partial_state(
    setup_test_db: None,
) -> None:
    """Force reload should fully replace both canonical and derived tables."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        first_file = Path(tmp_dir) / "first_mutations.csv"
        second_file = Path(tmp_dir) / "second_mutations.csv"

        first_data = pd.DataFrame(
            {
                "ModelID": ["ACH-000001", "ACH-000001"],
                "Chrom": ["1", "2"],
                "Pos": [1000000, 2000000],
                "Ref": ["A", "G"],
                "Alt": ["T", "C"],
                "HugoSymbol": ["GENE1", "GENE2"],
            }
        )
        second_data = pd.DataFrame(
            {
                "ModelID": ["ACH-000002"],
                "Chrom": ["17"],
                "Pos": [7577538],
                "Ref": ["C"],
                "Alt": ["T"],
                "HugoSymbol": ["TP53"],
                "LikelyLoF": [True],
                "HessDriver": [True],
            }
        )

        first_data.to_csv(first_file, index=False)
        second_data.to_csv(second_file, index=False)

        processor = MutationsProcessor()
        first_result = processor.process_file(first_file)
        assert first_result.status == "success"
        assert first_result.records_inserted == 2

        second_result = processor.process_file(second_file, force_reload=True)
        assert second_result.status == "success"
        assert second_result.records_inserted == 1

        db_manager = get_db_manager()
        assert (
            db_manager.execute("SELECT COUNT(*) FROM mutations").fetchone()[0]
            == 1
        )
        assert (
            db_manager.execute(
                "SELECT COUNT(*) FROM model_gene_mutation_status"
            ).fetchone()[0]
            == 1
        )

        row = db_manager.execute(
            """
            SELECT model_id, gene_symbol, mutation_count, has_likely_lof, has_driver
            FROM model_gene_mutation_status
            """
        ).fetchone()
        assert row == ("ACH-000002", "TP53", 1, True, True)


if __name__ == "__main__":
    pytest.main([__file__])
