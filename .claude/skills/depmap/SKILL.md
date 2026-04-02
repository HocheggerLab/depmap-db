---
name: depmap
description: >
  Analyse DepMap cancer dependency data. Use when the user asks about gene
  dependencies, expression, mutations, drug sensitivity, or cell line
  comparisons across cancer types. Supports DuckDB queries, Polars analysis,
  and publication-ready plotting.
allowed-tools: Bash(uv run python *), Read, Grep, Glob
argument-hint: "[question about DepMap data]"
---

# DepMap Data Analysis Skill

You are an expert analyst for the Cancer Dependency Map (DepMap) dataset stored
in a local DuckDB database. Answer the user's question: **$ARGUMENTS**

## Strategy — pick the right tool for the job

### 1. Simple lookups and counts → DuckDB SQL

For questions like "how many models have X?", "list all Y", or filtering a
single table, use DuckDB directly:

```python
import duckdb

con = duckdb.connect("/Users/hh65/.depmap/depmap.duckdb", read_only=True)
result = con.sql("SELECT ... FROM table_name WHERE ...").fetchdf()
con.close()
```

Key tables:
- `models` — cell line metadata (`model_id`, `cell_line_name`, `oncotree_lineage`, `oncotree_primary_disease`, `oncotree_subtype`)
- `genes` — gene identifiers (`hugo_symbol`, `entrez_id`, `ensembl_gene_id`)
- `mutations` — somatic variants (MAF-like, one row per mutation call)
- `model_gene_mutation_status` — sparse binary indicator (absent = WT)
- `compounds` — drug metadata (`broad_id`, `compound_name`, `moa`, `target_text`)
- `compound_targets` — compound-to-gene target bridge table
- `drug_screens` — screen metadata

### 2. Multi-table analysis → Polars scan API

For joins, reshaping, or multi-step analysis, use the Polars scan functions.
These return `pl.LazyFrame` via zero-copy Arrow from DuckDB.

```python
from depmap_db import scan_models, scan_gene_effects_wide, scan_mutations
import polars as pl
```

**Raw table scans** (return DuckDB tables directly):

| Function | Description |
|---|---|
| `scan_models()` | Cell line metadata |
| `scan_genes()` | Gene identifiers |
| `scan_gene_effects_wide()` | CRISPR dependency (wide: rows=models, cols=genes) |
| `scan_gene_expression_wide()` | RNA expression log2(TPM+1) (wide) |
| `scan_protein_expression_ms_wide()` | Proteomics (wide) |
| `scan_mutations()` | Mutation events (long) |
| `scan_model_gene_mutation_status()` | Sparse mutation indicator |
| `scan_compounds()` | Drug metadata |
| `scan_compound_targets()` | Compound-target mappings |
| `scan_drug_response_primary_wide()` | PRISM primary viability (wide) |
| `scan_drug_response_secondary()` | PRISM secondary dose-response fits |
| `scan_drug_response_secondary_dose()` | Per-dose viability |

**Curated datasets** (pre-joined, analysis-ready long format):

| Function | Description |
|---|---|
| `scan_mutation_events()` | Mutations + model metadata + resolved gene_symbol |
| `scan_proteomics_long()` | Unpivoted proteomics + model + protein annotations |
| `scan_drug_primary_long()` | Primary PRISM + compound + target + model metadata |
| `scan_drug_secondary_enriched()` | Secondary fits fully annotated |

**Important**: Wide table gene columns use bare symbols (e.g. `"MASTL"`).
To safely resolve a gene column:

```python
names = data.collect_schema().names()
col = gene if gene in names else next(c for c in names if c.startswith(f"{gene} ("))
```

### 3. Common visualisations → plotting library

For standard analysis patterns, use the built-in plotting functions which
produce publication-ready figures with the lab style sheet.

```python
from depmap_db.plots import lineage_analysis, mutation_analysis
```

**Dependency or expression by lineage** (box + strip, sorted by median):

```python
from depmap_db.plots import lineage_analysis
fig = lineage_analysis("MASTL", assay="dependency")   # or assay="expression"
fig.savefig("mastl_dependency_by_lineage.png", dpi=300)
```

**Mutation-stratified comparison** (WT vs mutant box + strip with Mann-Whitney U):

```python
from depmap_db.plots import mutation_analysis
# Tumour suppressor: filters by likely_lof OR tumor_suppressor_high_impact
fig = mutation_analysis("MASTL", "RB1", assay="dependency", gene_type="suppressor")

# Oncogene: filters by hotspot OR oncogene_high_impact OR hess_driver
fig = mutation_analysis("MAP2K1", "KRAS", assay="dependency", gene_type="oncogene")
```

Raises `ValueError` if fewer than 20 mutant models are found (configurable via `min_models`).

**Step-by-step access** (when you need the data before plotting):

```python
from depmap_db.plots import analyse_by_mutation, plot_mutation_comparison

df = analyse_by_mutation("MASTL", "TP53", assay="dependency", gene_type="suppressor")
# inspect df, filter further, then plot:
fig, ax = plot_mutation_comparison(df, "MASTL", "TP53", assay="dependency")
```

## Mutation impact filters

When filtering mutations, use the pre-computed boolean flags rather than
parsing `molecular_consequence`:

| Flag | Use for | Captures |
|---|---|---|
| `likely_lof` | Tumour suppressors | Nonsense, frameshift, splice-site |
| `tumor_suppressor_high_impact` | Tumour suppressors | DepMap curated high-impact |
| `hotspot` | Oncogenes | Recurrent mutation positions |
| `oncogene_high_impact` | Oncogenes | DepMap curated high-impact |
| `hess_driver` | Either | Hess et al. driver classification |

## Data model reminders

- **Wide format is canonical** for CRISPR, expression, proteomics, and PRISM primary
- **Long format** for mutations (event-based) and PRISM secondary (dose-response)
- Gene effect scores: negative = dependency, ~0 = no effect, -1 ~ median essential
- Expression values: log2(TPM+1)
- Drug viability: lower = more sensitive
- Join everything via `model_id`; genes via `hugo_symbol` or column name

## Output guidelines

- Write Python scripts and run with `uv run python script.py`
- Save figures to `reports/` directory as PNG (300 DPI)
- When showing data, print concise summaries — not raw DataFrames with thousands of rows
- Always report sample sizes (n) when making comparisons
- For statistical tests, report the test used, p-value, and effect size where appropriate
