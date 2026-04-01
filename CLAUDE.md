# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

**depmap-db** — local DuckDB-based system for downloading, storing, querying, and analyzing DepMap (Cancer Dependency Map) datasets. Provides CLI (`depmap-db`), Python API, and Polars-based analysis surfaces.

## Commands

```bash
# Install / sync
uv sync

# Run CLI
depmap-db init               # Create schema
depmap-db status             # Show DB stats
depmap-db download           # Download datasets
depmap-db load-folder --folder PATH

# Tests
uv run pytest                           # All tests
uv run pytest tests/test_mutations.py -v  # Single file

# Lint / format / type check
uv run ruff check .
uv run ruff format .
uv run mypy src
```

## Architecture

### Layered design

1. **CLI** (`cli.py`) — Click + Rich. Entry point: `depmap-db` command.
2. **Query layer** (`query.py`) — GeneQueryService, ProteinQueryService, MutationQueryService for analysis queries.
3. **ETL pipeline** (`etl/processors/`) — BaseProcessor abstract class with `load_file → validate_data → transform_data → insert_batch` contract. One processor per dataset type.
4. **Download manager** (`downloader/`) — DepMapClient (HTTP/async), FileManager (caching/checksums), ReleaseTracker (incremental refresh via `release_state.json`).
5. **Database** (`database/`) — ConnectionManager (DuckDB, singleton via `get_db_manager()`), SchemaManager (25+ tables, version-tracked).
6. **Polars integration** (`polars.py`) — Raw table exports to Parquet + curated long-format datasets (`mutation_events`, `proteomics_long`, `drug_primary_long`, `drug_secondary_enriched`).

### Data model decisions

- **Wide format is canonical** for CRISPR (`gene_effects_wide`), expression (`gene_expression_wide`), proteomics (`protein_expression_ms_wide`), and PRISM primary (`drug_response_primary_wide`) — matches DepMap's published orientation.
- **Long format for mutations** — event-based (MAF-like), variable events per model-gene pair. Derived `model_gene_mutation_status` provides sparse indicator for joins.
- **PRISM uses both** — primary is wide as published; secondary dose-response is long (`drug_response_secondary`). Separate release tracks from core DepMap.
- **Release tracking is per-dataset-type** — core DepMap, proteomics, and PRISM each have independent release labels.

### Adding a new dataset type

1. Create processor in `etl/processors/` extending `BaseProcessor`
2. Add table schema in `database/schema.py`
3. Register processor in `cli.py` (`_get_processor` mapping)
4. Add dataset info in `utils/constants.py` (`DEPMAP_FILES`)
5. Add tests

## Configuration

Pydantic v2 settings (`config.py`) with env var support: `DEPMAP__DATABASE__PATH`, `DEPMAP__DEPMAP__BATCH_SIZE`, etc. Priority: env vars > `.env.{ENV}` > `.env` > defaults.

## Conventions

- **Python 3.13+**, `uv` for deps (never pip)
- **ruff** for lint+format (line-length 79), **mypy** strict mode
- **Conventional commits** via commitizen (`cz_conventional_commits`)
- Tests use in-memory DuckDB (`:memory:`) fixtures
- Pre-commit hooks enforced
