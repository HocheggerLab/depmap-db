# depmap-db

Local DuckDB tooling for downloading, storing, and exporting DepMap datasets.

## Current design

The repository now treats the two large matrix datasets the same way:

- **CRISPR gene effect** → canonical table: `gene_effects_wide`
- **Gene expression** → canonical table: `gene_expression_wide`

This is a deliberate wide-matrix choice.

### Why wide is canonical here

- DepMap publishes both datasets as model-by-gene matrices already.
- The export layer and downstream analysis in this repo are wide-oriented.
- Supporting a partial long-path for expression but a wide-path for dependency made the code harder to reason about.
- A single canonical storage model is simpler than maintaining half-integrated wide/long toggles.

Long-format projections can still be added later as **derived views or exports** if a concrete use-case needs them, but the database now has one primary storage model for both matrix datasets.

## Install

```bash
git clone <repo>
cd depmap-db
uv venv
source .venv/bin/activate
uv sync
```

## CLI

### Initialize the database

```bash
depmap-db init
```

### Download phase-1 datasets

```bash
depmap-db download
```

### Download and load datasets

```bash
depmap-db download --datasets Model --datasets Gene --load-data
```

### Load a local folder of CSVs

```bash
depmap-db load-folder --folder /path/to/depmap/csvs
```

### Plan or apply a refresh

`refresh` now has a first clean abstraction point for repeatable updates:

- a configured **release label** (`DEPMAP_DEPMAP__RELEASE_LABEL`)
- a persisted **release tracking file**
- a refresh planner that decides which datasets are already cached for the tracked release

Plan only:

```bash
depmap-db refresh --datasets Model --datasets Gene
```

Apply the refresh plan:

```bash
depmap-db refresh --datasets Model --datasets Gene --apply --load-data
```

## Practical update strategy

A pragmatic recurring update loop should work like this:

1. **Detect or configure the target release**
   - simplest robust path: configure a release label manually per DepMap release
   - future improvement: fetch the DepMap release page and parse a release identifier into that same abstraction
2. **Map datasets to filenames**
   - keep this in a single dataset registry (`DEPMAP_FILES`)
3. **Build a refresh plan**
   - compare the target release label + filename mapping against the last applied snapshot
   - reuse valid cached files when nothing changed
4. **Download missing assets**
   - register cached files with metadata and checksums
5. **Load in a repeatable order**
   - metadata first (`Model`, `Gene`, `ModelCondition`), then matrix datasets
6. **Record applied state**
   - write the applied release snapshot so the next run can be incremental

### What should stay configurable

- release label
- cache directory
- release tracking file
- requested datasets
- whether refresh also loads data into DuckDB
- whether an apply should force-reload tables

## Validation

Run the standard checks:

```bash
uv run pytest
uv run ruff check .
uv run mypy src
```
