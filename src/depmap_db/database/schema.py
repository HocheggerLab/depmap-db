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

    CURRENT_VERSION = "1.4.0"

    def __init__(self) -> None:
        self.db_manager = get_db_manager()
        self._tables: dict[str, TableSchema] = {}
        self._define_schema()

    def _define_schema(self) -> None:
        """Define the complete database schema."""

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

        self._add_table(
            TableSchema(
                name="protein_features",
                table_type=TableType.METADATA,
                dependencies=["genes"],
                description=(
                    "Protein feature bridge for harmonized CCLE Gygi mass-spec "
                    "proteomics accessions and local gene mappings"
                ),
                sql="""
            CREATE TABLE IF NOT EXISTS protein_features (
                protein_accession VARCHAR PRIMARY KEY,
                protein_accession_base VARCHAR,
                storage_column_name VARCHAR NOT NULL,
                protein_entry_name VARCHAR,
                protein_name TEXT,
                gene_symbol VARCHAR,
                entrez_id INTEGER,
                ensembl_transcript_ids TEXT,
                local_gene_id VARCHAR,
                local_hugo_symbol VARCHAR,
                local_entrez_id INTEGER,
                local_ensembl_id VARCHAR,
                is_reviewed BOOLEAN,
                mapping_method VARCHAR,
                mapping_status VARCHAR,
                source_dataset VARCHAR NOT NULL,
                source_filename VARCHAR,
                modality VARCHAR,
                release_label VARCHAR,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (local_gene_id) REFERENCES genes(gene_id)
            )
            """,
            )
        )

        self._add_table(
            TableSchema(
                name="compounds",
                table_type=TableType.METADATA,
                dependencies=["schema_version"],
                description=(
                    "Compound-level PRISM metadata with explicit provenance and "
                    "room for partial metadata coverage across PRISM releases"
                ),
                sql="""
            CREATE TABLE IF NOT EXISTS compounds (
                broad_id VARCHAR PRIMARY KEY,
                compound_name VARCHAR,
                compound_synonyms TEXT,
                moa TEXT,
                target_text TEXT,
                smiles TEXT,
                phase VARCHAR,
                primary_screen_id VARCHAR,
                primary_dose_um DOUBLE,
                secondary_screen_id VARCHAR,
                secondary_row_name VARCHAR,
                secondary_passed_str_profiling BOOLEAN,
                disease_area VARCHAR,
                indication TEXT,
                source_dataset VARCHAR NOT NULL,
                source_filename VARCHAR,
                release_label VARCHAR,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
            """,
            )
        )

        self._add_table(
            TableSchema(
                name="drug_screens",
                table_type=TableType.METADATA,
                dependencies=["schema_version"],
                description=(
                    "PRISM screen/run metadata used to track release labels, release "
                    "tracks, and the default summary metric for secondary analyses"
                ),
                sql="""
            CREATE TABLE IF NOT EXISTS drug_screens (
                screen_id VARCHAR PRIMARY KEY,
                screen_kind VARCHAR NOT NULL,
                dataset_name VARCHAR NOT NULL,
                source_filename VARCHAR,
                release_label VARCHAR,
                release_track VARCHAR,
                default_secondary_summary_metric VARCHAR,
                notes TEXT,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
            """,
            )
        )

        self._add_table(
            TableSchema(
                name="compound_targets",
                table_type=TableType.METADATA,
                dependencies=["compounds", "genes"],
                description=(
                    "Many-to-many compound target bridge parsed from source target "
                    "strings with optional resolution into the local genes table"
                ),
                sql="""
            CREATE TABLE IF NOT EXISTS compound_targets (
                broad_id VARCHAR NOT NULL,
                target_ordinal INTEGER NOT NULL,
                target_symbol VARCHAR NOT NULL,
                source_target_text TEXT,
                source_dataset VARCHAR NOT NULL,
                source_filename VARCHAR,
                local_gene_id VARCHAR,
                local_hugo_symbol VARCHAR,
                local_entrez_id INTEGER,
                mapping_status VARCHAR NOT NULL,
                mapping_method VARCHAR NOT NULL,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                PRIMARY KEY (broad_id, target_ordinal, source_dataset),
                FOREIGN KEY (broad_id) REFERENCES compounds(broad_id),
                FOREIGN KEY (local_gene_id) REFERENCES genes(gene_id)
            )
            """,
            )
        )

        self._add_table(
            TableSchema(
                name="gene_effects_wide",
                table_type=TableType.CORE_DATA,
                dependencies=["models"],
                description="Canonical CRISPR gene effect matrix (models x genes)",
                sql="""
            CREATE TABLE IF NOT EXISTS gene_effects_wide (
                model_id VARCHAR PRIMARY KEY,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (model_id) REFERENCES models(model_id)
                -- Gene columns are added dynamically during data loading
            )
            """,
            )
        )

        self._add_table(
            TableSchema(
                name="gene_expression_wide",
                table_type=TableType.CORE_DATA,
                dependencies=["models", "model_conditions"],
                description="Canonical gene expression matrix (models x genes)",
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
                -- Gene columns are added dynamically during data loading
            )
            """,
            )
        )

        self._add_table(
            TableSchema(
                name="protein_expression_ms_wide",
                table_type=TableType.CORE_DATA,
                dependencies=["models", "protein_features"],
                description=(
                    "Canonical harmonized Gygi CCLE mass-spec proteomics matrix "
                    "(models x protein accessions)"
                ),
                sql="""
            CREATE TABLE IF NOT EXISTS protein_expression_ms_wide (
                model_id VARCHAR PRIMARY KEY,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (model_id) REFERENCES models(model_id)
                -- Protein accession columns are added dynamically during data loading
            )
            """,
            )
        )

        self._add_table(
            TableSchema(
                name="drug_response_primary_wide",
                table_type=TableType.CORE_DATA,
                dependencies=["compounds", "drug_screens"],
                description=(
                    "Canonical PRISM primary response matrix stored as published in "
                    "compound-by-model wide orientation"
                ),
                sql="""
            CREATE TABLE IF NOT EXISTS drug_response_primary_wide (
                broad_id VARCHAR PRIMARY KEY,
                screen_id VARCHAR,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (broad_id) REFERENCES compounds(broad_id),
                FOREIGN KEY (screen_id) REFERENCES drug_screens(screen_id)
                -- Model columns are added dynamically during data loading
            )
            """,
            )
        )

        self._add_table(
            TableSchema(
                name="drug_response_secondary",
                table_type=TableType.CORE_DATA,
                dependencies=["compounds", "models", "drug_screens"],
                description=(
                    "Canonical PRISM secondary long table with auc/ec50/ic50 and "
                    "dose-response fit terms"
                ),
                sql="""
            CREATE TABLE IF NOT EXISTS drug_response_secondary (
                response_id VARCHAR PRIMARY KEY,
                broad_id VARCHAR NOT NULL,
                model_id VARCHAR NOT NULL,
                ccle_name VARCHAR,
                screen_id VARCHAR,
                upper_limit DOUBLE,
                lower_limit DOUBLE,
                slope DOUBLE,
                r2 DOUBLE,
                auc DOUBLE,
                ec50 DOUBLE,
                ic50 DOUBLE,
                fit_name VARCHAR,
                successful_fit BOOLEAN,
                auc_riemann DOUBLE,
                minimum_dose_um DOUBLE,
                maximum_dose_um DOUBLE,
                source_project_id VARCHAR,
                passed_str_profiling BOOLEAN,
                row_name VARCHAR,
                compound_name VARCHAR,
                moa TEXT,
                target_text TEXT,
                disease_area VARCHAR,
                indication TEXT,
                smiles TEXT,
                phase VARCHAR,
                source_dataset VARCHAR NOT NULL,
                source_filename VARCHAR,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (broad_id) REFERENCES compounds(broad_id),
                FOREIGN KEY (model_id) REFERENCES models(model_id),
                FOREIGN KEY (screen_id) REFERENCES drug_screens(screen_id)
            )
            """,
            )
        )

        self._add_table(
            TableSchema(
                name="drug_response_secondary_dose",
                table_type=TableType.CORE_DATA,
                dependencies=["drug_response_secondary"],
                description=(
                    "Dose-level PRISM-like secondary response summaries collapsed to "
                    "one row per compound-model-screen-dose combination"
                ),
                sql="""
            CREATE TABLE IF NOT EXISTS drug_response_secondary_dose (
                dose_response_id VARCHAR PRIMARY KEY,
                response_id VARCHAR NOT NULL,
                broad_id VARCHAR NOT NULL,
                model_id VARCHAR NOT NULL,
                screen_id VARCHAR NOT NULL,
                dose_um DOUBLE NOT NULL,
                dose_unit VARCHAR,
                median_l2fc DOUBLE,
                median_l2fc_uncorrected DOUBLE,
                num_bio_reps INTEGER,
                pool_id VARCHAR,
                lua VARCHAR,
                cell_set VARCHAR,
                growth_pattern VARCHAR,
                pert_type VARCHAR,
                pert_vehicle VARCHAR,
                pert_plate VARCHAR,
                day INTEGER,
                source_project_id VARCHAR,
                source_dataset VARCHAR NOT NULL,
                source_filename VARCHAR,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (response_id) REFERENCES drug_response_secondary(response_id),
                FOREIGN KEY (broad_id) REFERENCES compounds(broad_id),
                FOREIGN KEY (model_id) REFERENCES models(model_id),
                FOREIGN KEY (screen_id) REFERENCES drug_screens(screen_id)
            )
            """,
            )
        )

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

        self._add_table(
            TableSchema(
                name="mutations",
                table_type=TableType.CORE_DATA,
                dependencies=["models"],
                description="Canonical somatic mutation event table (MAF-like)",
                sql="""
            CREATE TABLE IF NOT EXISTS mutations (
                mutation_id VARCHAR PRIMARY KEY,
                model_id VARCHAR NOT NULL,
                chrom VARCHAR,
                pos BIGINT,
                ref VARCHAR,
                alt VARCHAR,
                variant_type VARCHAR,
                dna_change VARCHAR,
                protein_change VARCHAR,
                hugo_symbol VARCHAR,
                ensembl_gene_id VARCHAR,
                entrez_gene_id VARCHAR,
                hgnc_name VARCHAR,
                molecular_consequence VARCHAR,
                vep_impact VARCHAR,
                af DOUBLE,
                dp INTEGER,
                ref_count INTEGER,
                alt_count INTEGER,
                gt VARCHAR,
                sift VARCHAR,
                polyphen VARCHAR,
                gnomade_af DOUBLE,
                gnomadg_af DOUBLE,
                revel_score DOUBLE,
                likely_lof BOOLEAN,
                hotspot BOOLEAN,
                oncogene_high_impact BOOLEAN,
                tumor_suppressor_high_impact BOOLEAN,
                hess_driver BOOLEAN,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (model_id) REFERENCES models(model_id)
            )
            """,
            )
        )

        self._add_table(
            TableSchema(
                name="model_gene_mutation_status",
                table_type=TableType.DERIVED,
                dependencies=["mutations", "models"],
                description="Sparse model-gene mutation status table (mutated pairs only)",
                sql="""
            CREATE TABLE IF NOT EXISTS model_gene_mutation_status (
                model_id VARCHAR NOT NULL,
                gene_symbol VARCHAR NOT NULL,
                is_mutated BOOLEAN NOT NULL DEFAULT TRUE,
                mutation_count INTEGER NOT NULL DEFAULT 1,
                has_likely_lof BOOLEAN DEFAULT FALSE,
                has_hotspot BOOLEAN DEFAULT FALSE,
                has_oncogene_high_impact BOOLEAN DEFAULT FALSE,
                has_tumor_suppressor_high_impact BOOLEAN DEFAULT FALSE,
                has_driver BOOLEAN DEFAULT FALSE,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                PRIMARY KEY (model_id, gene_symbol),
                FOREIGN KEY (model_id) REFERENCES models(model_id)
            )
            """,
            )
        )

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
            ready = []
            for table_name in remaining:
                table_schema = self._tables[table_name]
                if all(dep in ordered for dep in table_schema.dependencies):
                    ready.append(table_name)

            if not ready:
                raise RuntimeError(
                    f"Circular dependencies detected in tables: {remaining}"
                )

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
    """Create database tables in dependency order."""
    schema = Schema()

    if tables is None:
        table_names = schema.get_creation_order()
    else:
        for table_name in tables:
            if table_name not in schema._tables:
                raise ValueError(f"Unknown table: {table_name}")
        table_names = tables

    logger.info("Creating tables: %s", table_names)

    for table_name in table_names:
        schema.create_table(table_name)

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
        """
        INSERT INTO schema_version (version, description)
        SELECT ?, ?
        WHERE NOT EXISTS (
            SELECT 1 FROM schema_version WHERE version = ?
        )
        """,
        [version, description, version],
    )

    logger.info("Updated schema version to %s", version)


def drop_all_tables() -> None:
    """Drop all tables (useful for testing)."""
    db_manager = get_db_manager()
    schema = Schema()

    table_names = list(reversed(schema.get_creation_order()))

    for table_name in table_names:
        try:
            db_manager.execute(f"DROP TABLE IF EXISTS {table_name}")
            logger.info("Dropped table: %s", table_name)
        except (OSError, RuntimeError) as e:
            logger.error("Failed to drop table %s: %s", table_name, e)

    logger.info("All tables dropped successfully")
