# Gygi proteomics unmapped / partially resolved protein audit

Post-fix audit for `ProteomicsMSGygi` on branch `feature/depmap-proteomics-gygi-ms`.

## Coverage change

- Total protein features: 12,558
- UniProt gene-symbol coverage: 12,539 / 12,558 (unchanged; 19 still lack a UniProt primary gene symbol in the bridge search fields)
- Local gene-table coverage before fix: 12,533 / 12,558
- Local gene-table coverage after fix: 12,536 / 12,558
- Net improvement: **+3 local gene links**

## What changed

Implemented a safe local-gene resolution step in the proteomics bridge that uses raw `Gene.csv` metadata to:

1. resolve **unique HGNC alias / previous-symbol** matches (for example `UCC1 -> EPDR1`, `TTMP -> C3orf52`, `SLC61A1 -> MFSD5`), and
2. resolve **unique UniProt-ID -> gene** matches from `Gene.csv` when the accession belongs to exactly one local gene.

Ambiguous shared-accession cases are intentionally left unresolved.

## Classification summary

- `ambiguous_shared_uniprot_accession`: 3
- `deleted_uniprot_accession`: 1
- `fixed_alias_lookup`: 3
- `legacy_or_provisional_symbol_absent_from_local_genes`: 1
- `missing_uniprot_gene_symbol_entrez_resolved`: 4
- `missing_uniprot_gene_symbol_unresolved`: 13
- `nonstandard_retroelement_absent_from_local_genes`: 4

## Detailed accession audit

| accession | category | UniProt mapping result | local gene-link result | recommended handling |
| --- | --- | --- | --- | --- |
| `A4D1W8` | `fixed_alias_lookup` | TrEMBL A4D1W8_HUMAN; genes include UCC1 and EPDR1. | EPDR1 (EPDR1) (`exact+gene_symbol_alias_lookup` / `mapped_to_local_gene`) | Keep alias-based local lookup; resolved safely to EPDR1. |
| `B4DKC2` | `missing_uniprot_gene_symbol_entrez_resolved` | TrEMBL B4DKC2_HUMAN; no primary gene symbol; Entrez 22902. | RUFY3 (RUFY3) (`entrez_fallback` / `mapped_to_local_gene`) | Leave as Entrez fallback to RUFY3. |
| `B4DN88` | `missing_uniprot_gene_symbol_unresolved` | TrEMBL B4DN88_HUMAN; no primary gene symbol; no Entrez. | unresolved (`exact` / `mapped_to_uniprot_only`) | Document and leave unresolved unless UniProt/HGNC metadata improves. |
| `B4DXA9` | `missing_uniprot_gene_symbol_entrez_resolved` | TrEMBL B4DXA9_HUMAN; no primary gene symbol; Entrez 57541. | ZNF398 (ZNF398) (`entrez_fallback` / `mapped_to_local_gene`) | Leave as Entrez fallback to ZNF398. |
| `B4DXK4` | `missing_uniprot_gene_symbol_entrez_resolved` | TrEMBL B4DXK4_HUMAN; no primary gene symbol; Entrez 140807. | KRT72 (KRT72) (`entrez_fallback` / `mapped_to_local_gene`) | Leave as Entrez fallback to KRT72. |
| `B4E1Z4` | `missing_uniprot_gene_symbol_unresolved` | TrEMBL B4E1Z4_HUMAN; no primary gene symbol; no Entrez. | unresolved (`exact` / `mapped_to_uniprot_only`) | Document and leave unresolved; avoid name-based guessing. |
| `B7Z1Y9` | `missing_uniprot_gene_symbol_entrez_resolved` | TrEMBL B7Z1Y9_HUMAN; no primary gene symbol; Entrez 54933. | RHBDL2 (RHBDL2) (`entrez_fallback` / `mapped_to_local_gene`) | Leave as Entrez fallback to RHBDL2. |
| `E9PLD3` | `missing_uniprot_gene_symbol_unresolved` | TrEMBL E9PLD3_HUMAN; no primary gene symbol; no Entrez. | unresolved (`exact` / `mapped_to_uniprot_only`) | Document and leave unresolved. |
| `E9PSI1` | `missing_uniprot_gene_symbol_unresolved` | TrEMBL E9PSI1_HUMAN; no primary gene symbol; no Entrez. | unresolved (`exact` / `mapped_to_uniprot_only`) | Document and leave unresolved. |
| `F8W031` | `missing_uniprot_gene_symbol_unresolved` | TrEMBL F8W031_HUMAN; no primary gene symbol; no Entrez. | unresolved (`exact` / `mapped_to_uniprot_only`) | Document and leave unresolved. |
| `H0YHG0` | `missing_uniprot_gene_symbol_unresolved` | TrEMBL H0YHG0_HUMAN; no primary gene symbol; no Entrez. | unresolved (`exact` / `mapped_to_uniprot_only`) | Document and leave unresolved. |
| `H3BMH7` | `missing_uniprot_gene_symbol_unresolved` | TrEMBL H3BMH7_HUMAN; no primary gene symbol; no Entrez. | unresolved (`exact` / `mapped_to_uniprot_only`) | Document and leave unresolved. |
| `H3BQF6` | `missing_uniprot_gene_symbol_unresolved` | TrEMBL H3BQF6_HUMAN; no primary gene symbol; no Entrez. | unresolved (`exact` / `mapped_to_uniprot_only`) | Document and leave unresolved. |
| `H3BV83` | `missing_uniprot_gene_symbol_unresolved` | TrEMBL H3BV83_HUMAN; no primary gene symbol; no Entrez. | unresolved (`exact` / `mapped_to_uniprot_only`) | Document and leave unresolved. |
| `H7C469` | `missing_uniprot_gene_symbol_unresolved` | TrEMBL H7C469_HUMAN; no primary gene symbol; no Entrez. | unresolved (`exact` / `mapped_to_uniprot_only`) | Document and leave unresolved. |
| `I3L0E3` | `missing_uniprot_gene_symbol_unresolved` | TrEMBL I3L0E3_HUMAN; no primary gene symbol in search fields; no Entrez. | unresolved (`exact` / `mapped_to_uniprot_only`) | Document and leave unresolved. |
| `J3KQV8` | `deleted_uniprot_accession` | Inactive UniProt accession; deleted from sequence source. | unresolved (`exact` / `mapped_to_uniprot_only`) | Keep unresolved and optionally flag as obsolete/deleted in downstream reporting. |
| `O00370` | `nonstandard_retroelement_absent_from_local_genes` | Swiss-Prot LORF2_HUMAN; retrotransposable element ORF2; no primary gene symbol. | unresolved (`exact` / `mapped_to_uniprot_only`) | Treat as non-standard retroelement / control-like feature unless local gene model expands. |
| `O75920-2` | `ambiguous_shared_uniprot_accession` | Swiss-Prot SERF1_HUMAN; gene field includes SERF1A and SERF1B. | unresolved (`exact` / `mapped_to_uniprot_only`) | Leave unresolved unless a dataset-specific rule can disambiguate isoform-to-locus. |
| `P60507` | `nonstandard_retroelement_absent_from_local_genes` | Swiss-Prot EFC1_HUMAN; gene symbol ERVFC1; Entrez 105373297. | unresolved (`exact` / `mapped_to_uniprot_only`) | Treat as non-standard retroelement feature; do not force-map into HGNC genes. |
| `Q16385-2` | `ambiguous_shared_uniprot_accession` | Swiss-Prot SSX2_HUMAN; gene field includes SSX2 and SSX2B. | unresolved (`exact` / `mapped_to_uniprot_only`) | Leave unresolved unless a dataset-specific rule can disambiguate isoform-to-locus. |
| `Q5BVD1-3` | `fixed_alias_lookup` | Swiss-Prot TTMP_HUMAN; genes include TTMP and C3orf52. | C3orf52 (C3orf52) (`exact+gene_symbol_alias_lookup` / `mapped_to_local_gene`) | Keep alias-based local lookup; resolved safely to C3orf52. |
| `Q5HYB6` | `legacy_or_provisional_symbol_absent_from_local_genes` | TrEMBL Q5HYB6_HUMAN; gene symbol DKFZp686J1372; no Entrez. | unresolved (`exact` / `mapped_to_uniprot_only`) | Leave unresolved until a stable HGNC mapping exists in source gene metadata. |
| `Q6N075-2` | `fixed_alias_lookup` | Swiss-Prot MFSD5_HUMAN; genes include SLC61A1 and MFSD5. | MFSD5 (MFSD5) (`exact+gene_symbol_alias_lookup` / `mapped_to_local_gene`) | Keep alias/previous-symbol lookup; resolved safely to MFSD5. |
| `Q6ZSR9` | `missing_uniprot_gene_symbol_unresolved` | Swiss-Prot YJ005_HUMAN; no primary gene symbol; caution on sequence. | unresolved (`exact` / `mapped_to_uniprot_only`) | Document and leave unresolved. |
| `Q902F9` | `nonstandard_retroelement_absent_from_local_genes` | Swiss-Prot EN113_HUMAN; gene symbol HERVK_113. | unresolved (`exact` / `mapped_to_uniprot_only`) | Treat as non-standard retroelement feature; do not force-map into HGNC genes. |
| `Q96JG8-4` | `ambiguous_shared_uniprot_accession` | Swiss-Prot MAGD4_HUMAN; gene field includes MAGED4 and MAGED4B. | unresolved (`exact` / `mapped_to_uniprot_only`) | Leave unresolved unless a dataset-specific rule can disambiguate isoform-to-locus. |
| `Q9UF83` | `missing_uniprot_gene_symbol_unresolved` | Swiss-Prot YM012_HUMAN; no primary gene symbol; sequence caution present. | unresolved (`exact` / `mapped_to_uniprot_only`) | Document and leave unresolved. |
| `Q9UN81` | `nonstandard_retroelement_absent_from_local_genes` | Swiss-Prot LORF1_HUMAN; gene symbol L1RE1. | unresolved (`exact` / `mapped_to_uniprot_only`) | Treat as non-standard retroelement feature; do not force-map into HGNC genes. |


## Notes

- The four `entrez_fallback` cases (`B4DKC2`, `B4DXA9`, `B4DXK4`, `B7Z1Y9`) were already linkable to local genes because UniProt returned an Entrez cross-reference even though it did not expose a primary gene symbol in the search response.
- Three cases were improved safely by the new alias/previous-symbol logic: `A4D1W8 -> EPDR1`, `Q5BVD1-3 -> C3orf52`, and `Q6N075-2 -> MFSD5`.
- Three reviewed/shared-accession cases remain intentionally unresolved because a single UniProt accession maps to multiple local genes: `O75920-2`, `Q16385-2`, and `Q96JG8-4`.
- Four reviewed retroelement / endogenous retrovirus proteins (`O00370`, `P60507`, `Q902F9`, `Q9UN81`) are not represented in the current local genes table and should not be force-mapped to conventional HGNC genes.
- `J3KQV8` is an inactive/deleted UniProt accession and should remain unresolved.

Companion CSV: `reports/proteomics_gygi_unmapped_audit.csv`
