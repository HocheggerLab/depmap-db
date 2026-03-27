__version__ = "0.1.0"

from .config import set_env_vars

# Public API exports
from .database import (
    DatabaseManager,
    Schema,
    create_tables,
    get_current_schema_version,
    get_db_connection,
    get_db_manager,
)

# Initialize environment variables
set_env_vars()

__all__ = [
    "__version__",
    "set_env_vars",
    "DatabaseManager",
    "Schema",
    "create_tables",
    "get_current_schema_version",
    "get_db_connection",
    "get_db_manager",
]
