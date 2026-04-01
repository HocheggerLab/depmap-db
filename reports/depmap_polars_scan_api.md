# depmap-db Polars Scan API

> Direct DuckDB-to-Polars access via zero-copy Arrow transfer. All functions return `pl.LazyFrame` — call `.collect()` to materialise.

## Core Functions

| Function | Returns |
| --- | --- |
| `scan_table(name)` | Any single supported table as a LazyFrame |
| `scan_tables(names=None)` | Dict of name → LazyFrame (all tables if None) |
| `scan_dataset(name)` | Any single curated dataset as a LazyFrame |
| `scan_datasets(names=None)` | Dict of name → LazyFrame (all datasets if None) |

## Raw Table Scans

| Function | Table | Format | Description |
| --- | --- | --- | --- |
| `scan_models()` | `models` | Dimension | Cell line / cancer model metadata |
| `scan_genes()` | `genes` | Dimension | Gene identifiers (Hugo, Entrez, Ensembl, HGNC) |
| `scan_gene_effects_wide()` | `gene_effects_wide` | Wide matrix | CRISPR gene effect scores (Chronos/DEMETER2) |
| `scan_gene_expression_wide()` | `gene_expression_wide` | Wide matrix | RNA-seq expression, log2(TPM+1) |
| `scan_protein_expression_ms_wide()` | `protein_expression_ms_wide` | Wide matrix | Mass-spec proteomics (Gygi dataset) |
| `scan_mutations()` | `mutations` | Long (MAF-like) | Somatic mutation events |
| `scan_model_gene_mutation_status()` | `model_gene_mutation_status` | Sparse indicator | Binary mutated/WT per model-gene pair |
| `scan_compounds()` | `compounds` | Dimension | PRISM compound metadata |
| `scan_compound_targets()` | `compound_targets` | Bridge | Compound-target mappings (many-to-many) |
| `scan_drug_response_primary_wide()` | `drug_response_primary_wide` | Wide matrix | PRISM primary screen viability |
| `scan_drug_response_secondary()` | `drug_response_secondary` | Long | PRISM secondary dose-response fits |
| `scan_drug_response_secondary_dose()` | `drug_response_secondary_dose` | Long | Per-dose viability measurements |

Additional tables available via `scan_table()`: `drug_screens`, `protein_features`.

## Curated Dataset Scans

Pre-joined, analysis-ready long-format datasets with model and gene annotations.

| Function | Dataset | Description |
| --- | --- | --- |
| `scan_mutation_events()` | `mutation_events` | Mutations + model metadata (lineage, disease) + resolved gene_symbol |
| `scan_proteomics_long()` | `proteomics_long` | Unpivoted protein expression + model metadata + protein annotations |
| `scan_drug_primary_long()` | `drug_primary_long` | Unpivoted primary PRISM + compound metadata + targets + model metadata |
| `scan_drug_secondary_enriched()` | `drug_secondary_enriched` | Secondary fits + compound metadata + targets + model metadata |
