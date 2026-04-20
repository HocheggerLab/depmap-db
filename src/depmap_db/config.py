import contextlib
import logging
import os
from logging.handlers import RotatingFileHandler
from pathlib import Path

from dotenv import find_dotenv, load_dotenv
from pydantic import BaseModel, Field, field_validator
from pydantic_settings import BaseSettings, SettingsConfigDict

# Define project_root at module level
project_root = Path(__file__).parent.parent.parent.resolve()


def set_env_vars() -> None:
    """
    Load environment variables from configuration files.
    If ENV is not set, defaults to 'development'.
    Tries to load from .env.{ENV} first, then falls back to .env if needed.
    If no files are found, checks for required environment variables.

    Raises:
        OSError: If no configuration exists in files or environment.
    """
    env = os.getenv("ENV", "development").lower()
    package_root = Path(__file__).parent.parent.parent.resolve()

    print(f"[DEBUG set_env_vars] cwd={Path.cwd()}")
    print(f"[DEBUG set_env_vars] package_root={package_root}")
    print(f"[DEBUG set_env_vars] __file__={Path(__file__).resolve()}")

    # Try env-specific file first via directory-tree search
    env_file = find_dotenv(f".env.{env}", usecwd=True)
    print(f"[DEBUG set_env_vars] find_dotenv(.env.{env})={env_file!r}")
    if env_file:
        result = load_dotenv(env_file, override=True)
        print(
            f"[DEBUG set_env_vars] load_dotenv result={result}, LOG_LEVEL={os.getenv('LOG_LEVEL')!r}"
        )
        return

    # Fall back to .env via directory-tree search
    default_file = find_dotenv(".env", usecwd=True)
    print(f"[DEBUG set_env_vars] find_dotenv(.env)={default_file!r}")
    if default_file:
        result = load_dotenv(default_file, override=True)
        print(
            f"[DEBUG set_env_vars] load_dotenv result={result}, LOG_LEVEL={os.getenv('LOG_LEVEL')!r}"
        )
        return

    # Explicit path candidates as last resort
    candidates = dict.fromkeys([package_root, Path.cwd()])
    for root in candidates:
        default_env_path = root / ".env"
        print(
            f"[DEBUG set_env_vars] checking explicit path={default_env_path}, exists={default_env_path.exists()}"
        )
        if default_env_path.exists():
            result = load_dotenv(default_env_path, override=True)
            print(
                f"[DEBUG set_env_vars] load_dotenv result={result}, LOG_LEVEL={os.getenv('LOG_LEVEL')!r}"
            )
            return

    # If no files found, check for required environment variables
    required_vars = [
        "ENV",
        "LOG_LEVEL",
        "LOG_FORMAT",
        "ENABLE_CONSOLE_LOGGING",
        "ENABLE_FILE_LOGGING",
    ]

    if all(os.getenv(var) is not None for var in required_vars):
        # All required variables are present in environment
        return

    # If we get here, no configuration was found
    error_msg = "\n".join(
        [
            "No configuration found!",
            f"Current environment: {env}",
            f"Working directory: {Path.cwd()}",
            "Searched for .env and .env.{ENV} up the directory tree and in working directory.",
            "And checked environment variables for:",
            f"  - {', '.join(required_vars)}",
            "\nPlease create a .env file in the project root or set all required environment variables.",
        ]
    )
    raise OSError(error_msg)


def validate_env_vars() -> None:
    """
    Validate that all required environment variables are set.
    """
    required_vars = ["LOG_LEVEL", "LOG_FILE_PATH"]
    if missing_vars := [var for var in required_vars if not os.getenv(var)]:
        raise OSError(
            f"Missing required environment variables: {', '.join(missing_vars)}"
        )


def configure_log_handler(
    handler: logging.Handler,
    log_level: str,
    formatter: logging.Formatter,
    logger: logging.Logger,
) -> None:
    """Configure a logging handler with the specified settings.

    Args:
        handler: The logging handler to configure
        log_level: The logging level to set
        formatter: The formatter to use for log messages
        logger: The logger to add the handler to
    """
    handler.setLevel(getattr(logging, log_level, logging.DEBUG))
    handler.setFormatter(formatter)
    logger.addHandler(handler)


def get_logger(name: str) -> logging.Logger:
    """
    Get a logger with the given name, ensuring it's properly configured.
    If this is the first call, it will set up the root logger configuration.
    Subsequent calls will return appropriately named loggers that inherit the configuration.

    Args:
        name: The logger name, typically __name__ from the calling module

    Returns:
        logging.Logger: A configured logger instance
    """
    # Handle the case when module is run directly (__main__)
    if name == "__main__":
        # Get the caller's file path
        import inspect

        frame = inspect.stack()[1]
        module_path = Path(frame.filename)
        try:
            # Get relative path from project root to the module
            rel_path = module_path.relative_to(project_root / "src")
            # Convert path to module notation (my_app.submodule.file)
            module_name = str(rel_path.with_suffix("")).replace(os.sep, ".")
            name = module_name
        except ValueError:
            # Fallback if file is not in src directory
            name = module_path.stem

    # Get or create the logger
    logger = logging.getLogger(name)

    # If the root logger isn't configured yet, configure it
    root_logger = logging.getLogger()
    if not root_logger.handlers:
        set_env_vars()
        validate_env_vars()

        # Retrieve logging configurations from environment variables
        LOG_LEVEL = os.getenv("LOG_LEVEL", "WARNING").upper()
        LOG_FORMAT = os.getenv(
            "LOG_FORMAT",
            "%(asctime)s - %(name)s - %(levelname)s - %(filename)s:%(lineno)d - %(message)s",
        )
        ENABLE_CONSOLE_LOGGING = os.getenv(
            "ENABLE_CONSOLE_LOGGING", "False"
        ).lower() in ["true", "1", "yes"]
        ENABLE_FILE_LOGGING = os.getenv(
            "ENABLE_FILE_LOGGING", "False"
        ).lower() in ["true", "1", "yes"]
        LOG_FILE_PATH = os.getenv("LOG_FILE_PATH", "logs/app.log")
        LOG_MAX_BYTES = int(os.getenv("LOG_MAX_BYTES", 1048576))  # 1MB default
        LOG_BACKUP_COUNT = int(os.getenv("LOG_BACKUP_COUNT", 5))

        # Configure the root logger
        root_logger.setLevel(getattr(logging, LOG_LEVEL, logging.DEBUG))

        # Prevent propagation beyond our root logger
        root_logger.propagate = False

        # Formatter
        formatter = logging.Formatter(LOG_FORMAT)

        # Console Handler
        if ENABLE_CONSOLE_LOGGING:
            ch = logging.StreamHandler()
            configure_log_handler(ch, LOG_LEVEL, formatter, root_logger)

        # File Handler
        if ENABLE_FILE_LOGGING:
            log_path = Path(LOG_FILE_PATH)
            if log_dir := log_path.parent:
                log_dir.mkdir(parents=True, exist_ok=True)

            fh = RotatingFileHandler(
                LOG_FILE_PATH,
                maxBytes=LOG_MAX_BYTES,
                backupCount=LOG_BACKUP_COUNT,
            )
            configure_log_handler(fh, LOG_LEVEL, formatter, root_logger)

    return logger


# ========================
# Pydantic Settings Models
# ========================


class DatabaseSettings(BaseModel):
    """Database configuration settings."""

    path: Path = Field(
        default=Path.home() / ".depmap" / "depmap.duckdb",
        description="Path to the DuckDB database file",
    )
    memory: bool = Field(
        default=False, description="Use in-memory database instead of file"
    )
    threads: int | None = Field(
        default=None, description="Number of threads for DuckDB operations"
    )
    max_memory: str = Field(
        default="1GB", description="Maximum memory for DuckDB operations"
    )

    @field_validator("path")
    @classmethod
    def validate_db_path(cls, v: Path) -> Path:
        """Ensure the database directory exists."""
        if not v.parent.exists():
            v.parent.mkdir(parents=True, exist_ok=True)
        return v


class DepMapSettings(BaseModel):
    """DepMap data download configuration."""

    base_url: str = Field(
        default="https://depmap.org/portal/api/download/files",
        description="Manifest CSV endpoint for DepMap downloads",
    )
    datasets_url: str = Field(
        default="https://depmap.org/portal/api/download/datasets",
        description="Dataset catalogue JSON endpoint for DepMap downloads",
    )
    release_url: str = Field(
        default="https://depmap.org/portal/data_page/?tab=currentRelease",
        description="URL to check for current release information",
    )
    release_label: str = Field(
        default="DepMap Public 25Q3",
        description="Configured DepMap release label used for refresh planning",
    )
    release_tracking_file: Path = Field(
        default=Path("data/cache/release_state.json"),
        description="File used to track the last planned/applied release",
    )
    cache_dir: Path = Field(
        default=Path("data/cache"),
        description="Directory for caching downloaded files",
    )
    max_retries: int = Field(
        default=3,
        ge=1,
        le=10,
        description="Maximum number of download retry attempts",
    )
    timeout_seconds: int = Field(
        default=300, ge=30, le=3600, description="Request timeout in seconds"
    )
    batch_size: int = Field(
        default=10000,
        ge=1000,
        le=100000,
        description="Batch size for data processing operations",
    )
    verify_checksums: bool = Field(
        default=True,
        description="Verify file integrity using checksums when available",
    )

    @field_validator("cache_dir")
    @classmethod
    def validate_cache_dir(cls, v: Path) -> Path:
        """Ensure the cache directory exists."""
        v.mkdir(parents=True, exist_ok=True)
        return v

    @field_validator("release_tracking_file")
    @classmethod
    def validate_release_tracking_file(cls, v: Path) -> Path:
        """Ensure the release tracking directory exists."""
        v.parent.mkdir(parents=True, exist_ok=True)
        return v


class AppSettings(BaseSettings):
    """Main application settings combining all configurations."""

    model_config = SettingsConfigDict(
        env_file=".env",
        env_file_encoding="utf-8",
        env_nested_delimiter="__",
        env_prefix="DEPMAP_",
        case_sensitive=False,
    )

    # Environment and basic settings
    env: str = Field(
        default="development", description="Application environment"
    )
    debug: bool = Field(default=False, description="Enable debug mode")

    # Logging settings (keeping existing functionality)
    log_level: str = Field(default="INFO", description="Logging level")
    log_format: str = Field(
        default="%(asctime)s - %(name)s - %(levelname)s - %(filename)s:%(lineno)d - %(message)s",
        description="Log message format",
    )
    enable_console_logging: bool = Field(
        default=True, description="Enable console logging"
    )
    enable_file_logging: bool = Field(
        default=False, description="Enable file logging"
    )
    log_file_path: str = Field(
        default="logs/app.log", description="Log file path"
    )
    log_max_bytes: int = Field(
        default=1048576, description="Maximum log file size in bytes"
    )
    log_backup_count: int = Field(
        default=5, description="Number of backup log files"
    )

    # Nested configuration models
    database: DatabaseSettings = Field(default_factory=DatabaseSettings)
    depmap: DepMapSettings = Field(default_factory=DepMapSettings)

    @field_validator("log_level")
    @classmethod
    def validate_log_level(cls, v: str) -> str:
        """Validate log level is a valid logging level."""
        valid_levels = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
        v_upper = v.upper()
        if v_upper not in valid_levels:
            raise ValueError(
                f"Invalid log level: {v}. Must be one of {valid_levels}"
            )
        return v_upper


# Global settings instance
_settings: AppSettings | None = None


def get_settings() -> AppSettings:
    """Get the global application settings instance."""
    global _settings
    if _settings is None:
        # Load environment variables first using existing logic
        with contextlib.suppress(OSError):
            set_env_vars()

        # Create settings instance
        _settings = AppSettings()

    return _settings


def reload_settings() -> AppSettings:
    """Force reload of settings (useful for testing)."""
    global _settings
    _settings = None
    return get_settings()
