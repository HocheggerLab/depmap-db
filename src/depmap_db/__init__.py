__version__ = "0.1.0"

from .config import set_env_vars

# Initialize environment variables
set_env_vars()

# Public API exports
from .database import (
    DatabaseManager,
    Schema,
    create_tables,
    get_current_schema_version,
    get_db_connection,
    get_db_manager,
)
from .polars import (
    DEFAULT_POLARS_OUTPUT_DIR,
    SUPPORTED_POLARS_TABLES,
    export_polars_tables,
    get_lazy_tables,
    lazy_gene_effects_wide,
    lazy_gene_expression_wide,
    lazy_genes,
    lazy_model_gene_mutation_status,
    lazy_models,
    lazy_mutations,
    lazy_table,
    prepare_lazy_tables,
)



__all__ = [
    "__version__",
    "set_env_vars",
    "DatabaseManager",
    "Schema",
    "create_tables",
    "get_current_schema_version",
    "get_db_connection",
    "get_db_manager",
    "SUPPORTED_POLARS_TABLES",
    "DEFAULT_POLARS_OUTPUT_DIR",
    "export_polars_tables",
    "get_lazy_tables",
    "prepare_lazy_tables",
    "lazy_table",
    "lazy_models",
    "lazy_genes",
    "lazy_gene_effects_wide",
    "lazy_gene_expression_wide",
    "lazy_mutations",
    "lazy_model_gene_mutation_status",
]
