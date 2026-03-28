# depmap-db

Local DuckDB tooling for downloading, storing, and exporting DepMap datasets.

## Current design

The repository treats datasets according to their natural structure:

- **CRISPR gene effect** → canonical table: `gene_effects_wide` (wide matrix)
- **Gene expression** → canonical table: `gene_expression_wide` (wide matrix)
- **Somatic mutations** → canonical table: `mutations` (long event table) + derived `model_gene_mutation_status` (sparse)
- **Proteomics (Gygi MS)** → canonical table: `protein_expression_ms_wide` (wide matrix) + `protein_features` UniProt→gene bridge

### Why wide is canonical for dependency/expression/proteomics

- DepMap publishes these modalities as model-by-feature matrices already.
- The export layer and downstream analysis in this repo are wide-oriented.
- Supporting a partial long-path for expression but a wide-path for dependency made the code harder to reason about.
- A single canonical storage model is simpler than maintaining half-integrated wide/long toggles.

For proteomics, the initial implementation targets the **DepMap Harmonized MS CCLE Gygi** matrix only (`harmonized_MS_CCLE_Gygi.csv`). The corresponding release label currently resolves to **Harmonized Public Proteomics 24Q4**, which may move on a different release track from the core DepMap release. This repo records that release label explicitly rather than assuming it matches the main DepMap quarter.

The file is stored as-published. DepMap labels it as *harmonized*, but this repo does **not** claim a more specific normalization meaning unless DepMap documents that clearly in the source release.

Long-format projections can still be added later as **derived views or exports** if a concrete use-case needs them, but the database now has one primary storage model for both matrix datasets.

### Why long is canonical for mutations

- Mutations are naturally event-based (MAF-like) with variable numbers of events per model-gene pair.
- The `mutations` table preserves all source annotation columns for filtering and interpretation.
- The derived `model_gene_mutation_status` table provides a sparse model-gene mutation indicator for joining to wide dependency/expression tables.

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

### Gene-level query helpers

The first query-oriented CLI layer is intentionally small and composable:

```bash
# Which lineages are most dependent on HAPSTR1?
depmap-db gene dependency-summary HAPSTR1 --group-by lineage --limit 10

# Same idea, grouped by primary disease instead
depmap-db gene dependency-summary HAPSTR1 --group-by disease --limit 10

# Export per-model dependency values for MASTL in breast models
depmap-db gene dependency-models MASTL --lineage Breast --output mastl_breast.csv

# Summarise expression by lineage for a gene
depmap-db gene expression-summary HAPSTR1 --group-by lineage --limit 10

# Inspect Gygi MS bridge coverage
# shows how many UniProt accessions map to gene symbols / local genes
# and reports the proteomics release label tracked for the dataset

depmap-db protein mapping-summary

# Search the Gygi protein bridge by accession, gene symbol, or protein label

depmap-db protein search RBM47

# Summarise Gygi MS abundance by lineage for one protein accession

depmap-db protein expression-summary A0AV96 --group-by lineage --limit 10
```

`dependency-models` can also print to stdout in table, CSV, or JSON format via `--format`.

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

## Mutation query surfaces

The first mutation-aware query layer is now available:

```bash
# Inspect event-level mutations for a model / cell line
# accepts model_id, cell line name, stripped name, or unique partial match
# filters can be combined
# driver currently means the Hess driver flag
# output supports table/csv/json

depmap-db model mutations DLD-1 --gene KRAS --hotspot-only --format json

# Rank the most frequently mutated genes in a lineage
# includes total models, mutation frequency, and hotspot/LoF/driver model counts

depmap-db lineage mutation-frequency Lung --limit 10

# Compare dependency between mutant and WT cohorts
# useful for first-pass stratified hypothesis checks

depmap-db gene dependency-by-mutation KRAS --mutation-gene TP53 --mutation-class driver --lineage Lung
```

Notes:

- `model mutations` returns event-level rows from the canonical `mutations` table.
- `lineage mutation-frequency` is based on the derived `model_gene_mutation_status` table.
- `gene dependency-by-mutation` joins `model_gene_mutation_status` to `gene_effects_wide` and reports mutant/WT counts plus mean/median dependency and `delta_mean`.

## Polars Python API

Polars support is exposed as a **Python API for exploratory analysis**, not as a CLI surface.

### Why it works this way

A shell command cannot return a live `polars.LazyFrame`, and the cleanest lazy workflow today is still:

- DuckDB as the local source of truth
- export selected DuckDB tables to Parquet snapshots
- reopen those snapshots with `pl.scan_parquet(...)`

That is what `depmap_db.polars` does for you.

### Supported tables

- `models`
- `genes`
- `protein_features`
- `gene_effects_wide`
- `gene_expression_wide`
- `protein_expression_ms_wide`
- `mutations`
- `model_gene_mutation_status`

### Recommended usage pattern

In a notebook or REPL, call `prepare_lazy_tables(...)` for the tables you want to explore. This will refresh Parquet snapshots when needed and return Polars `LazyFrame`s ready for composition.

```python
import polars as pl
from depmap_db.polars import prepare_lazy_tables

tables = prepare_lazy_tables(
    tables=["models", "mutations", "model_gene_mutation_status"],
)

models = tables["models"]
mutations = tables["mutations"]
mutation_status = tables["model_gene_mutation_status"]

tp53_events = (
    mutations
    .filter(pl.col("hugo_symbol") == "TP53")
    .join(
        models.select(["model_id", "cell_line_name", "oncotree_lineage"]),
        on="model_id",
        how="left",
    )
)

print(tp53_events.limit(10).collect())
```

### Reuse existing snapshots

If the Parquet snapshots already exist and you just want to reopen them, use `get_lazy_tables(...)`:

```python
from depmap_db.polars import get_lazy_tables

tables = get_lazy_tables(tables=["models", "mutations"])
```

### Single-table helpers

For quick exploratory work, the module also exposes convenience helpers for individual tables:

```python
from depmap_db.polars import lazy_models, lazy_mutations

models = lazy_models()
mutations = lazy_mutations()
```

### Snapshot location

By default, snapshots are written to `data/polars/`. Override that with `output_dir=...` if you want notebook-specific or temporary outputs. You can also pass `db_path=...` to target a specific DuckDB file.

## Validation

Run the standard checks:

```bash
uv run pytest
uv run ruff check .
uv run mypy src
```
