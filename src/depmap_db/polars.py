"""Python and CLI-facing helpers for working with local DepMap data in Polars.

DuckDB remains the source of truth for local DepMap storage. This module
provides two ways to get data into Polars:

1. **Direct scans** (``scan_*`` functions): query DuckDB and return Polars
   DataFrames via zero-copy Arrow transfer. No intermediate files — always
   returns current data. Prefer this for interactive / notebook use.

2. **Parquet snapshots** (``lazy_*`` / ``export_*`` functions): materialise
   DuckDB tables to Parquet files that can be opened lazily with
   ``pl.scan_parquet(...)``. Useful for sharing data or working offline.

Both surfaces expose raw tables and curated (analysis-oriented) datasets.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from pathlib import Path
from typing import Final

import duckdb
import polars as pl

from .config import get_settings
from .database import DatabaseManager

SUPPORTED_POLARS_TABLES: Final[tuple[str, ...]] = (
    "models",
    "genes",
    "protein_features",
    "compounds",
    "compound_targets",
    "drug_screens",
    "gene_effects_wide",
    "gene_expression_wide",
    "protein_expression_ms_wide",
    "drug_response_primary_wide",
    "drug_response_secondary",
    "drug_response_secondary_dose",
    "mutations",
    "model_gene_mutation_status",
)
SUPPORTED_POLARS_DATASETS: Final[tuple[str, ...]] = (
    "mutation_events",
    "proteomics_long",
    "drug_primary_long",
    "drug_secondary_enriched",
)
DEFAULT_POLARS_OUTPUT_DIR: Final[Path] = Path("data/polars")

_PROTEOMICS_EXCLUDED_COLUMNS: Final[frozenset[str]] = frozenset(
    {"model_id", "created_at"}
)
_DRUG_PRIMARY_EXCLUDED_COLUMNS: Final[frozenset[str]] = frozenset(
    {"broad_id", "screen_id", "created_at"}
)


def resolve_polars_tables(tables: Iterable[str] | None = None) -> list[str]:
    """Resolve and validate the tables exposed through the Polars API."""
    return _resolve_names(
        requested=tables,
        supported=SUPPORTED_POLARS_TABLES,
        surface_label="tables",
    )


def resolve_polars_datasets(
    datasets: Iterable[str] | None = None,
) -> list[str]:
    """Resolve and validate the curated datasets exposed through the Polars API."""
    return _resolve_names(
        requested=datasets,
        supported=SUPPORTED_POLARS_DATASETS,
        surface_label="datasets",
    )


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
    resolved = (
        Path(output_dir)
        if output_dir is not None
        else DEFAULT_POLARS_OUTPUT_DIR
    )
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
                _copy_query_to_parquet(
                    db_manager,
                    query=f"SELECT * FROM {_quote_identifier(table_name)}",
                    output_path=output_path,
                )

            exported[table_name] = output_path

    return exported


def export_polars_datasets(
    *,
    output_dir: str | Path | None = None,
    datasets: Iterable[str] | None = None,
    overwrite: bool = False,
    db_path: str | Path | None = None,
) -> dict[str, Path]:
    """Export curated Polars-ready datasets to Parquet snapshots."""
    resolved_datasets = resolve_polars_datasets(datasets)
    resolved_db_path = resolve_database_path(db_path)
    resolved_output_dir = resolve_polars_output_dir(output_dir)
    db_mtime = resolved_db_path.stat().st_mtime

    exported: dict[str, Path] = {}
    with DatabaseManager(db_path=resolved_db_path) as db_manager:
        for dataset_name in resolved_datasets:
            output_path = resolved_output_dir / f"{dataset_name}.parquet"
            needs_refresh = (
                overwrite
                or not output_path.exists()
                or output_path.stat().st_mtime < db_mtime
            )

            if needs_refresh:
                query = _build_dataset_query(dataset_name, db_manager)
                _copy_query_to_parquet(
                    db_manager,
                    query=query,
                    output_path=output_path,
                )

            exported[dataset_name] = output_path

    return exported


def get_lazy_tables(
    *,
    output_dir: str | Path | None = None,
    tables: Iterable[str] | None = None,
) -> dict[str, pl.LazyFrame]:
    """Open existing Parquet snapshots for raw tables as LazyFrames."""
    resolved_tables = resolve_polars_tables(tables)
    return _get_lazy_snapshots(
        names=resolved_tables,
        output_dir=output_dir,
        missing_hint="Run export_polars_tables() or prepare_lazy_tables() first.",
    )


def get_lazy_datasets(
    *,
    output_dir: str | Path | None = None,
    datasets: Iterable[str] | None = None,
) -> dict[str, pl.LazyFrame]:
    """Open existing Parquet snapshots for curated datasets as LazyFrames."""
    resolved_datasets = resolve_polars_datasets(datasets)
    return _get_lazy_snapshots(
        names=resolved_datasets,
        output_dir=output_dir,
        missing_hint=(
            "Run export_polars_datasets() or prepare_lazy_datasets() first."
        ),
    )


def prepare_lazy_tables(
    *,
    output_dir: str | Path | None = None,
    tables: Iterable[str] | None = None,
    overwrite: bool = False,
    db_path: str | Path | None = None,
) -> dict[str, pl.LazyFrame]:
    """Export or refresh table snapshots, then return LazyFrames."""
    export_polars_tables(
        output_dir=output_dir,
        tables=tables,
        overwrite=overwrite,
        db_path=db_path,
    )
    return get_lazy_tables(output_dir=output_dir, tables=tables)


def prepare_lazy_datasets(
    *,
    output_dir: str | Path | None = None,
    datasets: Iterable[str] | None = None,
    overwrite: bool = False,
    db_path: str | Path | None = None,
) -> dict[str, pl.LazyFrame]:
    """Export or refresh curated datasets, then return LazyFrames."""
    export_polars_datasets(
        output_dir=output_dir,
        datasets=datasets,
        overwrite=overwrite,
        db_path=db_path,
    )
    return get_lazy_datasets(output_dir=output_dir, datasets=datasets)


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


def lazy_dataset(
    dataset_name: str,
    *,
    output_dir: str | Path | None = None,
    overwrite: bool = False,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Prepare and return one curated dataset as a LazyFrame."""
    return prepare_lazy_datasets(
        output_dir=output_dir,
        datasets=[dataset_name],
        overwrite=overwrite,
        db_path=db_path,
    )[dataset_name]


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


def lazy_protein_expression_ms_wide(
    *,
    output_dir: str | Path | None = None,
    overwrite: bool = False,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Return the ``protein_expression_ms_wide`` table as a LazyFrame."""
    return lazy_table(
        "protein_expression_ms_wide",
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
    """Return the ``model_gene_mutation_status`` table as a LazyFrame."""
    return lazy_table(
        "model_gene_mutation_status",
        output_dir=output_dir,
        overwrite=overwrite,
        db_path=db_path,
    )


def lazy_compounds(
    *,
    output_dir: str | Path | None = None,
    overwrite: bool = False,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Return the ``compounds`` table as a LazyFrame."""
    return lazy_table(
        "compounds",
        output_dir=output_dir,
        overwrite=overwrite,
        db_path=db_path,
    )


def lazy_compound_targets(
    *,
    output_dir: str | Path | None = None,
    overwrite: bool = False,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Return the ``compound_targets`` table as a LazyFrame."""
    return lazy_table(
        "compound_targets",
        output_dir=output_dir,
        overwrite=overwrite,
        db_path=db_path,
    )


def lazy_drug_response_primary_wide(
    *,
    output_dir: str | Path | None = None,
    overwrite: bool = False,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Return the ``drug_response_primary_wide`` table as a LazyFrame."""
    return lazy_table(
        "drug_response_primary_wide",
        output_dir=output_dir,
        overwrite=overwrite,
        db_path=db_path,
    )


def lazy_drug_response_secondary(
    *,
    output_dir: str | Path | None = None,
    overwrite: bool = False,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Return the ``drug_response_secondary`` table as a LazyFrame."""
    return lazy_table(
        "drug_response_secondary",
        output_dir=output_dir,
        overwrite=overwrite,
        db_path=db_path,
    )


def lazy_drug_response_secondary_dose(
    *,
    output_dir: str | Path | None = None,
    overwrite: bool = False,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Return the ``drug_response_secondary_dose`` table as a LazyFrame."""
    return lazy_table(
        "drug_response_secondary_dose",
        output_dir=output_dir,
        overwrite=overwrite,
        db_path=db_path,
    )


def lazy_mutation_events(
    *,
    output_dir: str | Path | None = None,
    overwrite: bool = False,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Return the curated ``mutation_events`` dataset as a LazyFrame."""
    return lazy_dataset(
        "mutation_events",
        output_dir=output_dir,
        overwrite=overwrite,
        db_path=db_path,
    )


def lazy_proteomics_long(
    *,
    output_dir: str | Path | None = None,
    overwrite: bool = False,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Return the curated ``proteomics_long`` dataset as a LazyFrame."""
    return lazy_dataset(
        "proteomics_long",
        output_dir=output_dir,
        overwrite=overwrite,
        db_path=db_path,
    )


def lazy_drug_primary_long(
    *,
    output_dir: str | Path | None = None,
    overwrite: bool = False,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Return the curated ``drug_primary_long`` dataset as a LazyFrame."""
    return lazy_dataset(
        "drug_primary_long",
        output_dir=output_dir,
        overwrite=overwrite,
        db_path=db_path,
    )


def lazy_drug_secondary_enriched(
    *,
    output_dir: str | Path | None = None,
    overwrite: bool = False,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Return the curated ``drug_secondary_enriched`` dataset."""
    return lazy_dataset(
        "drug_secondary_enriched",
        output_dir=output_dir,
        overwrite=overwrite,
        db_path=db_path,
    )


# ---------------------------------------------------------------------------
# Direct DuckDB → Polars scan API (no Parquet intermediate)
# ---------------------------------------------------------------------------


def _open_read_only_connection(
    db_path: str | Path | None = None,
) -> duckdb.DuckDBPyConnection:
    """Open a read-only DuckDB connection for scanning."""
    resolved = resolve_database_path(db_path)
    return duckdb.connect(str(resolved), read_only=True)


def scan_table(
    table_name: str,
    *,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Return a DuckDB table as a Polars LazyFrame via zero-copy Arrow.

    Args:
        table_name: Name of the table (must be in SUPPORTED_POLARS_TABLES).
        db_path: Optional path to the DuckDB database file.

    Returns:
        Polars LazyFrame backed by the query result.
    """
    resolve_polars_tables([table_name])
    conn = _open_read_only_connection(db_path)
    try:
        return (
            conn.sql(f"SELECT * FROM {_quote_identifier(table_name)}")
            .pl()
            .lazy()
        )
    finally:
        conn.close()


def scan_dataset(
    dataset_name: str,
    *,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Return a curated dataset as a Polars LazyFrame via zero-copy Arrow.

    Uses the same SQL reshaping logic as the Parquet export but skips the
    intermediate file — the result always reflects the current database state.

    Args:
        dataset_name: Name of the dataset (must be in SUPPORTED_POLARS_DATASETS).
        db_path: Optional path to the DuckDB database file.

    Returns:
        Polars LazyFrame backed by the query result.
    """
    resolve_polars_datasets([dataset_name])
    resolved_db_path = resolve_database_path(db_path)
    with DatabaseManager(db_path=resolved_db_path) as db_manager:
        query = _build_dataset_query(dataset_name, db_manager)
        conn = db_manager.connect()
        return conn.sql(query).pl().lazy()


def scan_tables(
    tables: Iterable[str] | None = None,
    *,
    db_path: str | Path | None = None,
) -> dict[str, pl.LazyFrame]:
    """Return multiple DuckDB tables as Polars LazyFrames.

    Args:
        tables: Table names to scan. ``None`` returns all supported tables.
        db_path: Optional path to the DuckDB database file.

    Returns:
        Mapping of table name to LazyFrame.
    """
    resolved = resolve_polars_tables(tables)
    return {name: scan_table(name, db_path=db_path) for name in resolved}


def scan_datasets(
    datasets: Iterable[str] | None = None,
    *,
    db_path: str | Path | None = None,
) -> dict[str, pl.LazyFrame]:
    """Return multiple curated datasets as Polars LazyFrames.

    Args:
        datasets: Dataset names to scan. ``None`` returns all supported.
        db_path: Optional path to the DuckDB database file.

    Returns:
        Mapping of dataset name to LazyFrame.
    """
    resolved = resolve_polars_datasets(datasets)
    return {name: scan_dataset(name, db_path=db_path) for name in resolved}


# -- Convenience wrappers (mirror the lazy_* Parquet API) ------------------


def scan_models(*, db_path: str | Path | None = None) -> pl.LazyFrame:
    """Return the ``models`` table directly from DuckDB."""
    return scan_table("models", db_path=db_path)


def scan_genes(*, db_path: str | Path | None = None) -> pl.LazyFrame:
    """Return the ``genes`` table directly from DuckDB."""
    return scan_table("genes", db_path=db_path)


def scan_gene_effects_wide(
    *, db_path: str | Path | None = None
) -> pl.LazyFrame:
    """Return the ``gene_effects_wide`` table directly from DuckDB."""
    return scan_table("gene_effects_wide", db_path=db_path)


def scan_gene_expression_wide(
    *, db_path: str | Path | None = None
) -> pl.LazyFrame:
    """Return the ``gene_expression_wide`` table directly from DuckDB."""
    return scan_table("gene_expression_wide", db_path=db_path)


def scan_protein_expression_ms_wide(
    *,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Return the ``protein_expression_ms_wide`` table directly from DuckDB."""
    return scan_table("protein_expression_ms_wide", db_path=db_path)


def scan_mutations(*, db_path: str | Path | None = None) -> pl.LazyFrame:
    """Return the ``mutations`` table directly from DuckDB."""
    return scan_table("mutations", db_path=db_path)


def scan_model_gene_mutation_status(
    *,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Return ``model_gene_mutation_status`` directly from DuckDB."""
    return scan_table("model_gene_mutation_status", db_path=db_path)


def scan_compounds(*, db_path: str | Path | None = None) -> pl.LazyFrame:
    """Return the ``compounds`` table directly from DuckDB."""
    return scan_table("compounds", db_path=db_path)


def scan_compound_targets(
    *, db_path: str | Path | None = None
) -> pl.LazyFrame:
    """Return the ``compound_targets`` table directly from DuckDB."""
    return scan_table("compound_targets", db_path=db_path)


def scan_drug_response_primary_wide(
    *,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Return ``drug_response_primary_wide`` directly from DuckDB."""
    return scan_table("drug_response_primary_wide", db_path=db_path)


def scan_drug_response_secondary(
    *,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Return ``drug_response_secondary`` directly from DuckDB."""
    return scan_table("drug_response_secondary", db_path=db_path)


def scan_drug_response_secondary_dose(
    *,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Return ``drug_response_secondary_dose`` directly from DuckDB."""
    return scan_table("drug_response_secondary_dose", db_path=db_path)


def scan_mutation_events(*, db_path: str | Path | None = None) -> pl.LazyFrame:
    """Return the curated ``mutation_events`` dataset directly from DuckDB."""
    return scan_dataset("mutation_events", db_path=db_path)


def scan_proteomics_long(*, db_path: str | Path | None = None) -> pl.LazyFrame:
    """Return the curated ``proteomics_long`` dataset directly from DuckDB."""
    return scan_dataset("proteomics_long", db_path=db_path)


def scan_drug_primary_long(
    *, db_path: str | Path | None = None
) -> pl.LazyFrame:
    """Return the curated ``drug_primary_long`` dataset directly from DuckDB."""
    return scan_dataset("drug_primary_long", db_path=db_path)


def scan_drug_secondary_enriched(
    *,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Return curated ``drug_secondary_enriched`` directly from DuckDB."""
    return scan_dataset("drug_secondary_enriched", db_path=db_path)


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _resolve_names(
    *,
    requested: Iterable[str] | None,
    supported: tuple[str, ...],
    surface_label: str,
) -> list[str]:
    """Resolve and validate requested names against a supported registry."""
    if requested is None:
        return list(supported)

    resolved = list(dict.fromkeys(requested))
    invalid = sorted(set(resolved) - set(supported))
    if invalid:
        supported_list = ", ".join(supported)
        invalid_list = ", ".join(invalid)
        raise ValueError(
            f"Unsupported Polars {surface_label}: {invalid_list}. "
            f"Available {surface_label}: {supported_list}"
        )
    return resolved


def _copy_query_to_parquet(
    db_manager: DatabaseManager,
    *,
    query: str,
    output_path: Path,
) -> None:
    """Copy the result of a SQL query to a Parquet file."""
    escaped_output = str(output_path).replace("'", "''")
    db_manager.execute(
        f"COPY ({query}) TO '{escaped_output}' (FORMAT PARQUET)"
    )


def _get_lazy_snapshots(
    *,
    names: Iterable[str],
    output_dir: str | Path | None,
    missing_hint: str,
) -> dict[str, pl.LazyFrame]:
    """Open named Parquet snapshots as Polars LazyFrames."""
    resolved_output_dir = resolve_polars_output_dir(output_dir)

    lazy_frames: dict[str, pl.LazyFrame] = {}
    for name in names:
        parquet_path = resolved_output_dir / f"{name}.parquet"
        if not parquet_path.exists():
            raise ValueError(
                f"Parquet snapshot missing for '{name}' at {parquet_path}. "
                f"{missing_hint}"
            )
        lazy_frames[name] = pl.scan_parquet(parquet_path)

    return lazy_frames


def _build_dataset_query(
    dataset_name: str, db_manager: DatabaseManager
) -> str:
    """Build the SQL query used to export one curated dataset."""
    builders: dict[str, Callable[[DatabaseManager], str]] = {
        "mutation_events": _build_mutation_events_query,
        "proteomics_long": _build_proteomics_long_query,
        "drug_primary_long": _build_drug_primary_long_query,
        "drug_secondary_enriched": _build_drug_secondary_enriched_query,
    }
    return builders[dataset_name](db_manager)


def _build_mutation_events_query(_: DatabaseManager) -> str:
    """Build the enriched mutation event dataset query."""
    return """
    SELECT
        mu.mutation_id,
        mu.model_id,
        m.cell_line_name,
        m.stripped_cell_line_name,
        m.ccle_name,
        m.oncotree_lineage,
        m.oncotree_primary_disease,
        m.oncotree_subtype,
        COALESCE(mu.hugo_symbol, mu.hgnc_name) AS gene_symbol,
        mu.hugo_symbol,
        mu.hgnc_name,
        mu.ensembl_gene_id,
        mu.entrez_gene_id,
        mu.chrom,
        mu.pos,
        mu.ref,
        mu.alt,
        mu.variant_type,
        mu.dna_change,
        mu.protein_change,
        mu.molecular_consequence,
        mu.vep_impact,
        mu.af,
        mu.dp,
        mu.ref_count,
        mu.alt_count,
        mu.gt,
        mu.sift,
        mu.polyphen,
        mu.gnomade_af,
        mu.gnomadg_af,
        mu.revel_score,
        mu.likely_lof,
        mu.hotspot,
        mu.oncogene_high_impact,
        mu.tumor_suppressor_high_impact,
        mu.hess_driver,
        mu.created_at
    FROM mutations mu
    LEFT JOIN models m ON m.model_id = mu.model_id
    """


def _build_proteomics_long_query(db_manager: DatabaseManager) -> str:
    """Build a long, model-protein proteomics dataset query."""
    if not db_manager.table_exists("protein_expression_ms_wide"):
        raise ValueError("Table 'protein_expression_ms_wide' does not exist")
    if not db_manager.table_exists("models"):
        raise ValueError("Table 'models' does not exist")
    if not db_manager.table_exists("protein_features"):
        raise ValueError("Table 'protein_features' does not exist")

    protein_columns = _get_dynamic_columns(
        db_manager,
        table_name="protein_expression_ms_wide",
        excluded_columns=_PROTEOMICS_EXCLUDED_COLUMNS,
    )
    unpivot_columns = ", ".join(
        _quote_identifier(name) for name in protein_columns
    )

    return f"""
    WITH long_matrix AS (
        SELECT
            model_id,
            storage_column_name,
            abundance
        FROM protein_expression_ms_wide
        UNPIVOT (abundance FOR storage_column_name IN ({unpivot_columns}))
    )
    SELECT
        lm.model_id,
        m.cell_line_name,
        m.stripped_cell_line_name,
        m.ccle_name,
        m.oncotree_lineage,
        m.oncotree_primary_disease,
        m.oncotree_subtype,
        pf.protein_accession,
        pf.protein_accession_base,
        pf.storage_column_name,
        pf.protein_entry_name,
        pf.protein_name,
        pf.gene_symbol,
        pf.entrez_id,
        pf.local_gene_id,
        pf.local_hugo_symbol,
        pf.local_entrez_id,
        pf.local_ensembl_id,
        pf.mapping_method,
        pf.mapping_status,
        pf.release_label,
        lm.abundance
    FROM long_matrix lm
    INNER JOIN models m ON m.model_id = lm.model_id
    LEFT JOIN protein_features pf
        ON pf.storage_column_name = lm.storage_column_name
    """


def _build_drug_primary_long_query(db_manager: DatabaseManager) -> str:
    """Build a long, model-compound primary PRISM response dataset query."""
    required_tables = (
        "drug_response_primary_wide",
        "models",
        "compounds",
        "drug_screens",
        "compound_targets",
    )
    for table_name in required_tables:
        if not db_manager.table_exists(table_name):
            raise ValueError(f"Table '{table_name}' does not exist")

    model_columns = _get_dynamic_columns(
        db_manager,
        table_name="drug_response_primary_wide",
        excluded_columns=_DRUG_PRIMARY_EXCLUDED_COLUMNS,
    )
    unpivot_columns = ", ".join(
        _quote_identifier(name) for name in model_columns
    )

    return f"""
    WITH target_summary AS (
        SELECT
            broad_id,
            string_agg(DISTINCT target_symbol, ', ') AS target_symbols,
            string_agg(
                DISTINCT COALESCE(local_hugo_symbol, target_symbol),
                ', '
            ) AS local_target_symbols,
            COUNT(*) AS target_count
        FROM compound_targets
        GROUP BY broad_id
    ),
    long_matrix AS (
        SELECT
            broad_id,
            screen_id,
            model_id,
            primary_response
        FROM drug_response_primary_wide
        UNPIVOT (primary_response FOR model_id IN ({unpivot_columns}))
    )
    SELECT
        lm.broad_id,
        c.compound_name,
        c.compound_synonyms,
        c.moa,
        c.target_text,
        ts.target_symbols,
        ts.local_target_symbols,
        ts.target_count,
        lm.screen_id,
        ds.screen_kind,
        ds.dataset_name,
        ds.release_label,
        ds.release_track,
        lm.model_id,
        m.cell_line_name,
        m.stripped_cell_line_name,
        m.ccle_name,
        m.oncotree_lineage,
        m.oncotree_primary_disease,
        m.oncotree_subtype,
        lm.primary_response
    FROM long_matrix lm
    LEFT JOIN compounds c ON c.broad_id = lm.broad_id
    LEFT JOIN target_summary ts ON ts.broad_id = lm.broad_id
    LEFT JOIN drug_screens ds ON ds.screen_id = lm.screen_id
    LEFT JOIN models m ON m.model_id = lm.model_id
    """


def _build_drug_secondary_enriched_query(db_manager: DatabaseManager) -> str:
    """Build an enriched secondary PRISM response dataset query."""
    required_tables = (
        "drug_response_secondary",
        "models",
        "compounds",
        "drug_screens",
        "compound_targets",
    )
    for table_name in required_tables:
        if not db_manager.table_exists(table_name):
            raise ValueError(f"Table '{table_name}' does not exist")

    return """
    WITH target_summary AS (
        SELECT
            broad_id,
            string_agg(DISTINCT target_symbol, ', ') AS target_symbols,
            string_agg(
                DISTINCT COALESCE(local_hugo_symbol, target_symbol),
                ', '
            ) AS local_target_symbols,
            COUNT(*) AS target_count
        FROM compound_targets
        GROUP BY broad_id
    )
    SELECT
        drs.response_id,
        drs.broad_id,
        COALESCE(drs.compound_name, c.compound_name) AS compound_name,
        c.compound_synonyms,
        COALESCE(drs.moa, c.moa) AS moa,
        COALESCE(drs.target_text, c.target_text) AS target_text,
        ts.target_symbols,
        ts.local_target_symbols,
        ts.target_count,
        drs.model_id,
        m.cell_line_name,
        m.stripped_cell_line_name,
        m.ccle_name,
        m.oncotree_lineage,
        m.oncotree_primary_disease,
        m.oncotree_subtype,
        drs.screen_id,
        ds.screen_kind,
        ds.dataset_name,
        ds.release_label,
        ds.release_track,
        ds.default_secondary_summary_metric,
        drs.upper_limit,
        drs.lower_limit,
        drs.slope,
        drs.r2,
        drs.auc,
        drs.ec50,
        drs.ic50,
        drs.fit_name,
        drs.successful_fit,
        drs.auc_riemann,
        drs.minimum_dose_um,
        drs.maximum_dose_um,
        drs.source_project_id,
        drs.passed_str_profiling,
        drs.row_name,
        drs.disease_area,
        drs.indication,
        drs.smiles,
        drs.phase,
        drs.source_dataset,
        drs.source_filename,
        drs.created_at
    FROM drug_response_secondary drs
    LEFT JOIN compounds c ON c.broad_id = drs.broad_id
    LEFT JOIN target_summary ts ON ts.broad_id = drs.broad_id
    LEFT JOIN drug_screens ds ON ds.screen_id = drs.screen_id
    LEFT JOIN models m ON m.model_id = drs.model_id
    """


def _get_dynamic_columns(
    db_manager: DatabaseManager,
    *,
    table_name: str,
    excluded_columns: frozenset[str],
) -> list[str]:
    """Return non-key dynamic columns for a wide matrix table."""
    column_names = [row[1] for row in db_manager.get_table_info(table_name)]
    dynamic_columns = [
        str(column_name)
        for column_name in column_names
        if str(column_name) not in excluded_columns
    ]
    if not dynamic_columns:
        excluded = ", ".join(sorted(excluded_columns))
        raise ValueError(
            f"Table '{table_name}' does not contain any Polars-exportable "
            f"value columns after excluding: {excluded}"
        )
    return dynamic_columns


def _quote_identifier(identifier: str) -> str:
    """Quote a SQL identifier for DuckDB."""
    escaped = identifier.replace('"', '""')
    return f'"{escaped}"'
