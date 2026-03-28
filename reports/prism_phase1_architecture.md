# PRISM phase-1 architecture note

## Goal

Add practical first-pass drug sensitivity support focused on DepMap PRISM without pretending to solve all pharmacology forever.

## Scope chosen for phase 1

Included now:

- **PRISM primary repurposing 24Q2**
  - `Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv`
  - `Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv`
- **PRISM secondary repurposing 20Q2**
  - `secondary-screen-dose-response-curve-parameters.csv`

Explicitly *not* solved yet:

- a grand unified pharmacology ontology
- target truth adjudication across external knowledge bases
- cross-dataset harmonisation across every drug screen in DepMap
- first-class query CLI for PRISM analysis beyond storage/ETL

## Release handling

PRISM is treated as its own release family rather than piggy-backing on the core DepMap release label.

Dataset registry additions:

- `PRISMPrimaryRepurposingExtended`
  - release track: `prism_primary`
  - release label override: `PRISM Primary Repurposing DepMap Public 24Q2`
- `PRISMPrimaryRepurposingCompoundList`
  - release track: `prism_primary`
  - release label override: `PRISM Primary Repurposing DepMap Public 24Q2`
- `PRISMSecondaryDoseResponseCurveParameters`
  - release track: `prism_secondary`
  - release label override: `PRISM Secondary Repurposing 20Q2`

This follows the same pattern already used for Gygi proteomics: each dataset can carry an explicit release label / track when it diverges from the main DepMap quarterly release.

## Schema

### `compounds`

Phase-1 compound metadata table keyed by `broad_id`.

Stores:

- source-facing compound name / synonyms
- MOA text
- source target text
- SMILES / phase / indication fields when present
- primary and secondary screen-specific metadata columns where useful
- dataset + source filename + release label provenance

### `compound_targets`

Many-to-many compound→target bridge.

Important semantics:

- one row per parsed target token from the source metadata
- keeps the original source target string
- may or may not resolve to a local `genes` row
- mapping is conservative and nullable
- a resolved row means “this target token could be uniquely matched to a local gene record”, **not** “this is confirmed mechanism-of-action truth”

### `drug_screens`

Minimal screen/run metadata table for PRISM.

Used for:

- release label / release track tracking
- screen kind (`primary` / `secondary`)
- recording that **AUC** is the default secondary-screen summary metric in this repo
- attaching notes about storage / semantics

### `drug_response_primary_wide`

Stores the PRISM primary matrix **as published** in wide compound-by-model form.

Why keep it wide?

- the source file is already wide and that orientation is meaningful
- phase 1 should preserve the source faithfully
- it avoids an expensive and less obviously useful pivot into a huge long table before we know the main analysis patterns

### `drug_response_secondary`

Canonical long table for secondary dose-response fits.

Stores:

- `auc`, `ec50`, `ic50`
- fit terms: `upper_limit`, `lower_limit`, `slope`, `r2`
- model / compound / screen identifiers
- source metadata columns (`moa`, `target_text`, `disease_area`, etc.)

## ETL design

### Primary processor

`PrismPrimaryWideProcessor`

Responsibilities:

- load / validate the wide matrix
- create a dynamic wide table with one column per model
- create / update a `drug_screens` row for the primary screen
- load the companion compound list when present
- populate `compounds`
- parse `repurposing_target` into `compound_targets`

Fallback behavior:

- if the companion compound list is missing, still load the wide matrix and create minimal `compounds` rows from the broad IDs

### Secondary processor

`PrismSecondaryProcessor`

Responsibilities:

- load / validate the long dose-response parameter table
- populate `drug_screens`
- derive / upsert `compounds`
- parse `target` into `compound_targets`
- load the canonical long table into `drug_response_secondary`

### Target resolution strategy

Current resolution path is intentionally conservative:

1. parse comma-separated target tokens from the source string
2. try exact match against local `genes.hugo_symbol`
3. if available, also allow unique alias / previous-symbol resolution via cached `Gene.csv`
4. if a token cannot be uniquely resolved, keep it in `compound_targets` with `mapping_status = 'source_only_unresolved'`

No external target database is consulted in phase 1.

## Response semantics

### Primary matrix

Stored exactly as a primary response matrix. No attempt is made to reinterpret or normalise the published values beyond preserving them in DuckDB.

### Secondary matrix

The full fitted response record is retained.

Default analytical summary for phase 1: **AUC**.

Reason:

- more robust as a first-pass summary than over-relying on a single fitted concentration metric
- still keeps EC50 / IC50 available for downstream work
- consistent with the idea that the canonical table should preserve the full fit, while the repo can recommend one default summary metric for common analyses

## What remains for next pass

- PRISM query helpers (compound-centric / target-centric summaries)
- optional exports / derived views for long projections of the primary matrix
- better metadata enrichment using PRISM cell-line / treatment sidecar tables
- more nuanced target-token parsing (multi-target phrases, pathway-level terms, uncertainty categories)
- tests against small real PRISM fixtures or sampled release snippets
- possible support for incremental per-dataset upserts rather than full-table replacement during reloads
