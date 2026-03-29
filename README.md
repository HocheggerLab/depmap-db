# depmap-db

Local DuckDB tooling for downloading, storing, and exporting DepMap datasets.

## Current design

The repository treats datasets according to their natural structure:

- **CRISPR gene effect** → canonical table: `gene_effects_wide` (wide matrix)
- **Gene expression** → canonical table: `gene_expression_wide` (wide matrix)
- **Somatic mutations** → canonical table: `mutations` (long event table) + derived `model_gene_mutation_status` (sparse)
- **Proteomics (Gygi MS)** → canonical table: `protein_expression_ms_wide` (wide matrix) + `protein_features` UniProt→gene bridge
- **PRISM primary repurposing** → canonical table: `drug_response_primary_wide` (compound × model wide matrix, stored as published) + `compounds` + `compound_targets`
- **PRISM secondary repurposing** → canonical table: `drug_response_secondary` (long dose-response table with AUC/EC50/IC50 + fit terms) + `compounds` + `compound_targets` + `drug_screens`

### Why wide is canonical for dependency/expression/proteomics

- DepMap publishes these modalities as model-by-feature matrices already.
- The export layer and downstream analysis in this repo are wide-oriented.
- Supporting a partial long-path for expression but a wide-path for dependency made the code harder to reason about.
- A single canonical storage model is simpler than maintaining half-integrated wide/long toggles.

For proteomics, the initial implementation targets the **DepMap Harmonized MS CCLE Gygi** matrix only (`harmonized_MS_CCLE_Gygi.csv`). The corresponding release label currently resolves to **Harmonized Public Proteomics 24Q4**, which may move on a different release track from the core DepMap release. This repo records that release label explicitly rather than assuming it matches the main DepMap quarter.

The file is stored as-published. DepMap labels it as *harmonized*, but this repo does **not** claim a more specific normalization meaning unless DepMap documents that clearly in the source release.

Long-format projections can still be added later as **derived views or exports** if a concrete use-case needs them, but the database now has one primary storage model for both matrix datasets.

### Why PRISM phase-1 uses both wide and long tables

PRISM is different enough from CRISPR/expression/proteomics that phase 1 uses the published structure of each release rather than forcing everything into one drug-response abstraction immediately.

- **Primary repurposing matrix (`Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv`)** is stored in `drug_response_primary_wide` in the same **compound-by-model** orientation DepMap publishes.
- **Secondary dose-response parameters (`secondary-screen-dose-response-curve-parameters.csv`)** are stored canonically in `drug_response_secondary` as a long table.
- **Compound metadata** lives in `compounds`.
- **Compound→gene links** live in `compound_targets` and are explicitly provenance-aware. A row in `compound_targets` means the source metadata named that target token; it does **not** claim biochemical certainty or exclusivity.
- **Screen metadata / release semantics** live in `drug_screens`, which lets PRISM sit on its own release tracks (`prism_primary`, `prism_secondary`) rather than pretending it always matches the core DepMap quarterly release.

For secondary-screen analyses, this repo keeps multiple response measures but recommends **AUC** as the default summary metric in phase 1 because it is usually more stable than leaning on a single fitted concentration summary alone. EC50, IC50, slope, upper/lower limits, and fit quality are still preserved.

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

### Load PRISM phase-1 datasets

```bash
# download the primary matrix + secondary dose-response table
# release labels are tracked explicitly per PRISM dataset

depmap-db download \
  --datasets PRISMPrimaryRepurposingExtended \
  --datasets PRISMSecondaryDoseResponseCurveParameters

# load a local folder that contains the PRISM files
# if the primary compound list sits next to the primary matrix,
# it will also be used to populate compound metadata + source target strings

depmap-db load-folder --folder /path/to/depmap/prism/files
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

## Polars API and export surface

Polars support now has two layers:

1. **raw table snapshots** for direct access to canonical DuckDB tables
2. **curated datasets** for notebook-friendly long/enriched analysis surfaces

The core idea is still the same:

- DuckDB is the local source of truth
- selected tables or curated dataset queries are exported to Parquet snapshots
- Polars reopens those snapshots lazily via `pl.scan_parquet(...)`

A CLI cannot return a live `polars.LazyFrame`, but it *can* materialize the Parquet snapshots you want.

### CLI: export Polars snapshots

```bash
# raw tables

depmap-db polars export \
  --table models \
  --table mutations

# curated long-format datasets for the main downstream modalities

depmap-db polars export \
  --dataset proteomics_long \
  --dataset mutation_events \
  --dataset drug_primary_long \
  --dataset drug_secondary_enriched
```

### Raw tables exposed to Polars

- `models`
- `genes`
- `protein_features`
- `compounds`
- `compound_targets`
- `drug_screens`
- `gene_effects_wide`
- `gene_expression_wide`
- `protein_expression_ms_wide`
- `drug_response_primary_wide`
- `drug_response_secondary`
- `drug_response_secondary_dose`
- `mutations`
- `model_gene_mutation_status`

### Curated Polars datasets

These are the recommended notebook-facing entry points when you want the main Helfrid-relevant modalities in a comparable way.

- `mutation_events` → mutation event table enriched with model metadata
- `proteomics_long` → Gygi MS matrix reshaped to one row per model × protein with protein bridge metadata
- `drug_primary_long` → PRISM primary matrix reshaped to one row per compound × model with compound/target metadata
- `drug_secondary_enriched` → PRISM secondary dose-response summaries enriched with model/compound/screen metadata

### Recommended Python usage

For Helfrid-style downstream analysis, the intended starting point is:

- `proteomics_long` for model × protein abundance
- `mutation_events` for event-level mutation context
- `drug_primary_long` for broad primary PRISM sensitivity screens
- `drug_secondary_enriched` for fitted secondary PRISM response summaries

Use `prepare_lazy_datasets(...)` when you want analysis-ready long datasets:

```python
import polars as pl
from depmap_db.polars import prepare_lazy_datasets

frames = prepare_lazy_datasets(
    datasets=[
        "proteomics_long",
        "mutation_events",
        "drug_primary_long",
        "drug_secondary_enriched",
    ]
)

proteomics = frames["proteomics_long"]
mutations = frames["mutation_events"]
drug_primary = frames["drug_primary_long"]
drug_secondary = frames["drug_secondary_enriched"]

rbm47_tp53 = (
    proteomics
    .filter(pl.col("gene_symbol") == "RBM47")
    .join(
        mutations
        .filter(pl.col("gene_symbol") == "TP53")
        .select(["model_id", "gene_symbol", "protein_change", "hotspot"]),
        on="model_id",
        how="left",
    )
)

print(rbm47_tp53.limit(10).collect())
```

Use `prepare_lazy_tables(...)` when you want the canonical stored tables directly:

```python
from depmap_db.polars import prepare_lazy_tables

tables = prepare_lazy_tables(
    tables=["models", "mutations", "model_gene_mutation_status"],
)
```

### Reuse existing snapshots

If the Parquet snapshots already exist and you just want to reopen them:

```python
from depmap_db.polars import get_lazy_datasets, get_lazy_tables

frames = get_lazy_datasets(datasets=["proteomics_long", "drug_primary_long"])
tables = get_lazy_tables(tables=["models", "mutations"])
```

### Single-object helpers

For quick exploratory work, the module also exposes convenience helpers for both tables and curated datasets:

```python
from depmap_db.polars import lazy_models, lazy_mutations, lazy_proteomics_long

models = lazy_models()
mutations = lazy_mutations()
proteomics = lazy_proteomics_long()
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
