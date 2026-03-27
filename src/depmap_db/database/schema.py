"""Database schema management for DepMap data."""

from dataclasses import dataclass
from enum import Enum

from ..config import get_logger
from .connection import get_db_manager

logger = get_logger(__name__)


class TableType(Enum):
    """Types of tables in the schema."""

    METADATA = "metadata"
    CORE_DATA = "core_data"
    DERIVED = "derived"
    SYSTEM = "system"


@dataclass
class TableSchema:
    """Schema definition for a database table."""

    name: str
    sql: str
    table_type: TableType
    dependencies: list[str]
    description: str


class Schema:
    """Database schema management."""

    CURRENT_VERSION = "1.0.0"

    def __init__(self) -> None:
        self.db_manager = get_db_manager()
        self._tables: dict[str, TableSchema] = {}
        self._define_schema()

    def _define_schema(self) -> None:
        """Define the complete database schema."""

        # System tables
        self._add_table(
            TableSchema(
                name="schema_version",
                table_type=TableType.SYSTEM,
                dependencies=[],
                description="Tracks database schema version",
                sql="""
            CREATE TABLE IF NOT EXISTS schema_version (
                version VARCHAR PRIMARY KEY,
                applied_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                description TEXT
            )
            """,
            )
        )

        # Gene metadata table
        self._add_table(
            TableSchema(
                name="genes",
                table_type=TableType.METADATA,
                dependencies=["schema_version"],
                description="Gene metadata and identifiers",
                sql="""
            CREATE TABLE IF NOT EXISTS genes (
                gene_id VARCHAR PRIMARY KEY,
                hugo_symbol VARCHAR NOT NULL,
                entrez_id INTEGER,
                ensembl_id VARCHAR,
                gene_type VARCHAR,
                chromosome VARCHAR,
                start_position BIGINT,
                end_position BIGINT,
                strand VARCHAR(1),
                description TEXT,
                synonyms TEXT[],
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
            """,
            )
        )

        # Model metadata table
        self._add_table(
            TableSchema(
                name="models",
                table_type=TableType.METADATA,
                dependencies=["schema_version"],
                description="Cancer model/cell line metadata",
                sql="""
            CREATE TABLE IF NOT EXISTS models (
                model_id VARCHAR PRIMARY KEY,
                patient_id VARCHAR,
                cell_line_name VARCHAR NOT NULL,
                stripped_cell_line_name VARCHAR,
                depmap_model_type VARCHAR,
                oncotree_lineage VARCHAR,
                oncotree_primary_disease VARCHAR,
                oncotree_subtype VARCHAR,
                oncotree_code VARCHAR,
                rrid VARCHAR,
                age INTEGER,
                age_category VARCHAR,
                sex VARCHAR,
                patient_race VARCHAR,
                primary_or_metastasis VARCHAR,
                sample_collection_site VARCHAR,
                source_type VARCHAR,
                source_detail TEXT,
                model_type VARCHAR,
                tissue_origin VARCHAR,
                growth_pattern VARCHAR,
                ccle_name VARCHAR,
                cosmic_id VARCHAR,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
            """,
            )
        )

        # Model conditions table
        self._add_table(
            TableSchema(
                name="model_conditions",
                table_type=TableType.METADATA,
                dependencies=["models"],
                description="Conditions under which models were assayed",
                sql="""
            CREATE TABLE IF NOT EXISTS model_conditions (
                model_condition_id VARCHAR PRIMARY KEY,
                model_id VARCHAR NOT NULL,
                parent_model_condition_id VARCHAR,
                data_source VARCHAR,
                cell_format VARCHAR,
                passage_number VARCHAR,
                growth_media VARCHAR,
                plate_coating VARCHAR,
                serum_free_media BOOLEAN,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (model_id) REFERENCES models(model_id)
            )
            """,
            )
        )

        # CRISPR gene effects table (wide format - genes as columns)
        # Note: This table will be created dynamically based on the actual gene columns
        self._add_table(
            TableSchema(
                name="gene_effects_wide",
                table_type=TableType.CORE_DATA,
                dependencies=["models"],
                description="CRISPR gene effect scores in wide format (models x genes matrix)",
                sql="""
            CREATE TABLE IF NOT EXISTS gene_effects_wide (
                model_id VARCHAR PRIMARY KEY,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (model_id) REFERENCES models(model_id)
                -- Gene columns will be added dynamically during data loading
            )
            """,
            )
        )

        # Gene expression long format table
        self._add_table(
            TableSchema(
                name="gene_expression",
                table_type=TableType.CORE_DATA,
                dependencies=["models", "model_conditions"],
                description="Gene expression TPM values in long format",
                sql="""
            CREATE TABLE IF NOT EXISTS gene_expression (
                model_id VARCHAR NOT NULL,
                gene_id VARCHAR NOT NULL,
                expression_value DOUBLE NOT NULL,
                sequencing_id VARCHAR,
                model_condition_id VARCHAR,
                is_default_entry BOOLEAN DEFAULT TRUE,
                is_default_for_mc BOOLEAN DEFAULT TRUE,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                PRIMARY KEY (model_id, gene_id),
                FOREIGN KEY (model_id) REFERENCES models(model_id),
                FOREIGN KEY (model_condition_id) REFERENCES model_conditions(model_condition_id)
            )
            """,
            )
        )

        # Gene expression wide format table
        self._add_table(
            TableSchema(
                name="gene_expression_wide",
                table_type=TableType.CORE_DATA,
                dependencies=["models", "model_conditions"],
                description="Gene expression TPM values in wide format (models x genes matrix)",
                sql="""
            CREATE TABLE IF NOT EXISTS gene_expression_wide (
                model_id VARCHAR PRIMARY KEY,
                sequencing_id VARCHAR,
                model_condition_id VARCHAR,
                is_default_entry BOOLEAN DEFAULT TRUE,
                is_default_for_mc BOOLEAN DEFAULT TRUE,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (model_id) REFERENCES models(model_id),
                FOREIGN KEY (model_condition_id) REFERENCES model_conditions(model_condition_id)
                -- Gene columns will be added dynamically during data loading
            )
            """,
            )
        )

        # Screen metadata table
        self._add_table(
            TableSchema(
                name="screens",
                table_type=TableType.METADATA,
                dependencies=["model_conditions"],
                description="CRISPR screen metadata and QC information",
                sql="""
            CREATE TABLE IF NOT EXISTS screens (
                screen_id VARCHAR PRIMARY KEY,
                model_condition_id VARCHAR NOT NULL,
                model_id VARCHAR NOT NULL,
                screen_type VARCHAR,
                library VARCHAR,
                days INTEGER,
                pdna_batch VARCHAR,
                passes_qc BOOLEAN,
                screen_nnmd DOUBLE,
                screen_roc_auc DOUBLE,
                screen_fpr DOUBLE,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (model_condition_id) REFERENCES model_conditions(model_condition_id),
                FOREIGN KEY (model_id) REFERENCES models(model_id)
            )
            """,
            )
        )

        # Data lineage tracking
        self._add_table(
            TableSchema(
                name="data_imports",
                table_type=TableType.SYSTEM,
                dependencies=["schema_version"],
                description="Track data import operations",
                sql="""
            CREATE TABLE IF NOT EXISTS data_imports (
                import_id VARCHAR PRIMARY KEY,
                import_type VARCHAR NOT NULL,
                source_file VARCHAR NOT NULL,
                file_checksum VARCHAR,
                records_imported INTEGER,
                status VARCHAR NOT NULL,
                started_at TIMESTAMP NOT NULL,
                completed_at TIMESTAMP,
                error_message TEXT,
                metadata JSON
            )
            """,
            )
        )

    def _add_table(self, table_schema: TableSchema) -> None:
        """Add a table schema to the registry."""
        self._tables[table_schema.name] = table_schema

    def get_table_schema(self, table_name: str) -> TableSchema | None:
        """Get schema for a specific table."""
        return self._tables.get(table_name)

    def get_tables_by_type(self, table_type: TableType) -> list[TableSchema]:
        """Get all tables of a specific type."""
        return [t for t in self._tables.values() if t.table_type == table_type]

    def get_creation_order(self) -> list[str]:
        """Get table names in dependency order for creation."""
        ordered: list[str] = []
        remaining = list(self._tables.keys())

        while remaining:
            # Find tables with all dependencies satisfied
            ready = []
            for table_name in remaining:
                table_schema = self._tables[table_name]
                if all(dep in ordered for dep in table_schema.dependencies):
                    ready.append(table_name)

            if not ready:
                raise RuntimeError(
                    f"Circular dependencies detected in tables: {remaining}"
                )

            # Add ready tables to ordered list and remove from remaining
            ordered.extend(ready)
            for table_name in ready:
                remaining.remove(table_name)

        return ordered

    def create_table(self, table_name: str) -> None:
        """Create a single table."""
        if table_name not in self._tables:
            raise ValueError(f"Unknown table: {table_name}")

        table_schema = self._tables[table_name]
        logger.info("Creating table: %s", table_name)

        try:
            self.db_manager.execute(table_schema.sql)
            logger.info("Successfully created table: %s", table_name)
        except (OSError, RuntimeError) as e:
            logger.error("Failed to create table %s: %s", table_name, e)
            raise


def create_tables(tables: list[str] | None = None) -> None:
    """Create database tables in dependency order.

    Args:
        tables: Optional list of specific tables to create. If None, creates all tables.
    """
    schema = Schema()

    if tables is None:
        # Create all tables in dependency order
        table_names = schema.get_creation_order()
    else:
        # Validate requested tables exist
        for table_name in tables:
            if table_name not in schema._tables:
                raise ValueError(f"Unknown table: {table_name}")
        table_names = tables

    logger.info("Creating tables: %s", table_names)

    for table_name in table_names:
        schema.create_table(table_name)

    # Update schema version
    _update_schema_version(schema.CURRENT_VERSION)

    logger.info("Database schema creation completed successfully")


def get_current_schema_version() -> str | None:
    """Get the current database schema version."""
    db_manager = get_db_manager()

    if not db_manager.table_exists("schema_version"):
        return None

    try:
        result = db_manager.execute(
            "SELECT version FROM schema_version ORDER BY applied_at DESC LIMIT 1"
        )
        row = result.fetchone()
        return row[0] if row else None
    except (OSError, RuntimeError) as e:
        logger.error("Failed to get schema version: %s", e)
        return None


def _update_schema_version(
    version: str, description: str | None = None
) -> None:
    """Update the schema version record."""
    db_manager = get_db_manager()

    description = description or f"Schema version {version}"

    db_manager.execute(
        "INSERT INTO schema_version (version, description) VALUES (?, ?)",
        [version, description],
    )

    logger.info("Updated schema version to %s", version)


def drop_all_tables() -> None:
    """Drop all tables (useful for testing)."""
    db_manager = get_db_manager()
    schema = Schema()

    # Drop in reverse dependency order
    table_names = list(reversed(schema.get_creation_order()))

    for table_name in table_names:
        try:
            db_manager.execute(f"DROP TABLE IF EXISTS {table_name}")
            logger.info("Dropped table: %s", table_name)
        except (OSError, RuntimeError) as e:
            logger.error("Failed to drop table %s: %s", table_name, e)

    logger.info("All tables dropped successfully")
