"""Python helpers for working with local DepMap data in Polars.

This module is intentionally a *Python API*, not a CLI surface. DuckDB remains
the source of truth for local DepMap storage, while this module exports
Parquet snapshots that Polars can scan lazily via ``pl.scan_parquet(...)``.
That keeps exploratory notebook workflows ergonomic without pretending the CLI
can hand back live ``polars.LazyFrame`` objects.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Final

import polars as pl

from .config import get_settings
from .database import DatabaseManager

SUPPORTED_POLARS_TABLES: Final[tuple[str, ...]] = (
    "models",
    "genes",
    "gene_effects_wide",
    "gene_expression_wide",
    "mutations",
    "model_gene_mutation_status",
)
DEFAULT_POLARS_OUTPUT_DIR: Final[Path] = Path("data/polars")


def resolve_polars_tables(tables: Iterable[str] | None = None) -> list[str]:
    """Resolve and validate the tables exposed through the Polars API."""
    if tables is None:
        return list(SUPPORTED_POLARS_TABLES)

    requested = list(dict.fromkeys(tables))
    invalid = sorted(set(requested) - set(SUPPORTED_POLARS_TABLES))
    if invalid:
        supported = ", ".join(SUPPORTED_POLARS_TABLES)
        invalid_list = ", ".join(invalid)
        raise ValueError(
            f"Unsupported Polars tables: {invalid_list}. "
            f"Available tables: {supported}"
        )
    return requested


def resolve_database_path(db_path: str | Path | None = None) -> Path:
    """Resolve the file-backed DuckDB database path used by the Polars API."""
    settings = get_settings()

    if db_path is None:
        if settings.database.memory:
            raise ValueError(
                "Polars lazy exports require a file-backed DuckDB database; "
                "in-memory databases are not supported."
            )
        resolved = settings.database.path
    else:
        resolved = Path(db_path)

    if not resolved.exists():
        raise ValueError(f"DuckDB database not found at {resolved}")

    return resolved


def resolve_polars_output_dir(output_dir: str | Path | None = None) -> Path:
    """Resolve the Parquet snapshot directory used for lazy Polars scans."""
    resolved = Path(output_dir) if output_dir is not None else DEFAULT_POLARS_OUTPUT_DIR
    resolved.mkdir(parents=True, exist_ok=True)
    return resolved


def export_polars_tables(
    *,
    output_dir: str | Path | None = None,
    tables: Iterable[str] | None = None,
    overwrite: bool = False,
    db_path: str | Path | None = None,
) -> dict[str, Path]:
    """Export supported DuckDB tables to Parquet snapshots for Polars.

    Existing files are reused unless they are missing, older than the database
    file, or ``overwrite`` is set.
    """
    resolved_tables = resolve_polars_tables(tables)
    resolved_db_path = resolve_database_path(db_path)
    resolved_output_dir = resolve_polars_output_dir(output_dir)
    db_mtime = resolved_db_path.stat().st_mtime

    exported: dict[str, Path] = {}
    with DatabaseManager(db_path=resolved_db_path) as db_manager:
        for table_name in resolved_tables:
            if not db_manager.table_exists(table_name):
                raise ValueError(
                    f"Table '{table_name}' does not exist in {resolved_db_path}"
                )

            output_path = resolved_output_dir / f"{table_name}.parquet"
            needs_refresh = (
                overwrite
                or not output_path.exists()
                or output_path.stat().st_mtime < db_mtime
            )

            if needs_refresh:
                escaped_output = str(output_path).replace("'", "''")
                db_manager.execute(
                    f"COPY (SELECT * FROM {table_name}) "
                    f"TO '{escaped_output}' (FORMAT PARQUET)"
                )

            exported[table_name] = output_path

    return exported


def get_lazy_tables(
    *,
    output_dir: str | Path | None = None,
    tables: Iterable[str] | None = None,
) -> dict[str, pl.LazyFrame]:
    """Open existing Parquet snapshots as Polars LazyFrames."""
    resolved_tables = resolve_polars_tables(tables)
    resolved_output_dir = resolve_polars_output_dir(output_dir)

    lazy_tables: dict[str, pl.LazyFrame] = {}
    for table_name in resolved_tables:
        parquet_path = resolved_output_dir / f"{table_name}.parquet"
        if not parquet_path.exists():
            raise ValueError(
                f"Parquet snapshot missing for '{table_name}' at {parquet_path}. "
                "Run export_polars_tables() or prepare_lazy_tables() first."
            )
        lazy_tables[table_name] = pl.scan_parquet(parquet_path)

    return lazy_tables


def prepare_lazy_tables(
    *,
    output_dir: str | Path | None = None,
    tables: Iterable[str] | None = None,
    overwrite: bool = False,
    db_path: str | Path | None = None,
) -> dict[str, pl.LazyFrame]:
    """Export or refresh Parquet snapshots, then return LazyFrames."""
    export_polars_tables(
        output_dir=output_dir,
        tables=tables,
        overwrite=overwrite,
        db_path=db_path,
    )
    return get_lazy_tables(output_dir=output_dir, tables=tables)


def lazy_table(
    table_name: str,
    *,
    output_dir: str | Path | None = None,
    overwrite: bool = False,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Prepare and return a single supported table as a LazyFrame."""
    return prepare_lazy_tables(
        output_dir=output_dir,
        tables=[table_name],
        overwrite=overwrite,
        db_path=db_path,
    )[table_name]


def lazy_models(
    *,
    output_dir: str | Path | None = None,
    overwrite: bool = False,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Return the ``models`` table as a Polars LazyFrame."""
    return lazy_table(
        "models",
        output_dir=output_dir,
        overwrite=overwrite,
        db_path=db_path,
    )


def lazy_genes(
    *,
    output_dir: str | Path | None = None,
    overwrite: bool = False,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Return the ``genes`` table as a Polars LazyFrame."""
    return lazy_table(
        "genes",
        output_dir=output_dir,
        overwrite=overwrite,
        db_path=db_path,
    )


def lazy_gene_effects_wide(
    *,
    output_dir: str | Path | None = None,
    overwrite: bool = False,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Return the ``gene_effects_wide`` table as a Polars LazyFrame."""
    return lazy_table(
        "gene_effects_wide",
        output_dir=output_dir,
        overwrite=overwrite,
        db_path=db_path,
    )


def lazy_gene_expression_wide(
    *,
    output_dir: str | Path | None = None,
    overwrite: bool = False,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Return the ``gene_expression_wide`` table as a Polars LazyFrame."""
    return lazy_table(
        "gene_expression_wide",
        output_dir=output_dir,
        overwrite=overwrite,
        db_path=db_path,
    )


def lazy_mutations(
    *,
    output_dir: str | Path | None = None,
    overwrite: bool = False,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Return the ``mutations`` table as a Polars LazyFrame."""
    return lazy_table(
        "mutations",
        output_dir=output_dir,
        overwrite=overwrite,
        db_path=db_path,
    )


def lazy_model_gene_mutation_status(
    *,
    output_dir: str | Path | None = None,
    overwrite: bool = False,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Return the ``model_gene_mutation_status`` table as a Polars LazyFrame."""
    return lazy_table(
        "model_gene_mutation_status",
        output_dir=output_dir,
        overwrite=overwrite,
        db_path=db_path,
    )

