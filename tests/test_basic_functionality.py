"""Basic functionality tests for DepMap database infrastructure."""

import tempfile
from pathlib import Path

import pytest

from depmap_db.config import get_settings, reload_settings
from depmap_db.database import (
    create_tables,
    get_current_schema_version,
    get_db_manager,
    reset_db_manager,
)
from depmap_db.downloader import RefreshPlanner
from depmap_db.downloader.file_manager import FileManager
from depmap_db.downloader.releases import ReleaseTracker
from depmap_db.utils.constants import DEPMAP_FILES


def test_settings_configuration() -> None:
    """Test that settings are properly configured."""
    settings = get_settings()

    assert settings.env == "development"
    assert settings.log_level == "INFO"
    assert settings.database.max_memory.endswith("GB")
    assert settings.depmap.max_retries == 3
    assert settings.depmap.release_label == "currentRelease"


def test_dataset_tables_use_wide_matrix_storage() -> None:
    """Core matrix datasets should point at wide canonical tables."""
    assert DEPMAP_FILES["CRISPRGeneEffect"].table_name == "gene_effects_wide"
    assert DEPMAP_FILES["GeneExpression"].table_name == "gene_expression_wide"


def test_database_creation() -> None:
    """Test database schema creation."""
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = Path(tmp_dir) / "depmap-test.duckdb"

        settings = reload_settings()
        settings.database.path = tmp_path
        settings.database.memory = False

        reset_db_manager()
        create_tables()

        schema_version = get_current_schema_version()
        assert schema_version is not None
        assert schema_version == "1.1.0"

        db_manager = get_db_manager()
        assert db_manager.table_exists("schema_version")
        assert db_manager.table_exists("models")
        assert db_manager.table_exists("genes")
        assert db_manager.table_exists("gene_effects_wide")
        assert db_manager.table_exists("gene_expression_wide")
        assert db_manager.table_exists("model_conditions")
        assert db_manager.table_exists("screens")
        assert db_manager.table_exists("data_imports")
        assert not db_manager.table_exists("gene_expression")

        reset_db_manager()


def test_in_memory_database() -> None:
    """Test in-memory database functionality."""
    settings = reload_settings()
    settings.database.memory = True

    reset_db_manager()

    try:
        create_tables()

        schema_version = get_current_schema_version()
        assert schema_version == "1.1.0"

        db_manager = get_db_manager()
        assert db_manager.table_exists("models")

    finally:
        reset_db_manager()


def test_refresh_planner_uses_release_tracking(tmp_path: Path) -> None:
    """Refresh planning should reuse cache when the tracked release matches."""
    tracking_file = tmp_path / "release_state.json"
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()

    settings = reload_settings()
    settings.depmap.release_tracking_file = tracking_file
    settings.depmap.cache_dir = cache_dir
    settings.depmap.release_label = "depmap-24Q4"

    file_manager = FileManager(cache_dir)
    tracker = ReleaseTracker(tracking_file)
    planner = RefreshPlanner(file_manager=file_manager, tracker=tracker)

    initial_plan = planner.build_plan(["Model", "Gene"])
    assert set(initial_plan.datasets_to_download) == {"Model", "Gene"}
    assert initial_plan.cached_datasets == []

    tracker.save(initial_plan.snapshot)
    model_csv = cache_dir / "Model.csv"
    model_csv.write_text("model_id\nACH-000001\n")
    file_manager.register_download(
        "Model",
        model_csv,
        "https://depmap.org/portal/api/download/file/Model.csv",
        force_checksum=False,
    )

    second_plan = planner.build_plan(["Model", "Gene"])
    assert second_plan.cached_datasets == ["Model"]
    assert second_plan.datasets_to_download == ["Gene"]


if __name__ == "__main__":
    pytest.main([__file__])
