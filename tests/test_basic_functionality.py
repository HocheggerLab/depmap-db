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


def test_settings_configuration():
    """Test that settings are properly configured."""
    settings = get_settings()

    assert settings.env == "development"
    assert settings.log_level == "INFO"
    assert settings.database.max_memory == "1GB"
    assert settings.depmap.max_retries == 3


def test_database_creation():
    """Test database schema creation."""
    # Use temporary database for testing
    with tempfile.NamedTemporaryFile(
        suffix=".duckdb", delete=False
    ) as tmp_file:
        tmp_path = Path(tmp_file.name)

    try:
        # Override database path for testing
        settings = reload_settings()
        settings.database.path = tmp_path
        settings.database.memory = False

        # Reset database manager to use new settings
        reset_db_manager()

        # Create tables
        create_tables()

        # Verify schema was created
        schema_version = get_current_schema_version()
        assert schema_version is not None
        assert schema_version == "1.0.0"

        # Verify core tables exist
        db_manager = get_db_manager()
        assert db_manager.table_exists("schema_version")
        assert db_manager.table_exists("models")
        assert db_manager.table_exists("genes")
        assert db_manager.table_exists("gene_effects")
        assert db_manager.table_exists("model_conditions")
        assert db_manager.table_exists("screens")
        assert db_manager.table_exists("data_imports")

    finally:
        # Clean up
        reset_db_manager()
        if tmp_path.exists():
            tmp_path.unlink()


def test_in_memory_database():
    """Test in-memory database functionality."""
    # Configure for in-memory database
    settings = reload_settings()
    settings.database.memory = True

    # Reset database manager
    reset_db_manager()

    try:
        # Create tables
        create_tables()

        # Verify it works
        schema_version = get_current_schema_version()
        assert schema_version == "1.0.0"

        db_manager = get_db_manager()
        assert db_manager.table_exists("models")

    finally:
        reset_db_manager()


if __name__ == "__main__":
    pytest.main([__file__])
