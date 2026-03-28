"""Constants and configuration data for DepMap datasets."""

from __future__ import annotations

from typing import NamedTuple


class DatasetInfo(NamedTuple):
    """Information about a DepMap dataset file."""

    filename: str
    description: str
    table_name: str
    priority: int
    modality: str = "other"
    release_track: str = "depmap_core"
    release_label_override: str | None = None


DEPMAP_FILES: dict[str, DatasetInfo] = {
    "CRISPRGeneEffect": DatasetInfo(
        filename="CRISPRGeneEffect.csv",
        description="Gene effect estimates for all models, integrated using Chronos",
        table_name="gene_effects_wide",
        priority=1,
        modality="crispr",
    ),
    "Model": DatasetInfo(
        filename="Model.csv",
        description="Metadata describing all cancer models/cell lines",
        table_name="models",
        priority=1,
        modality="metadata",
    ),
    "Gene": DatasetInfo(
        filename="Gene.csv",
        description="Gene metadata and identifiers",
        table_name="genes",
        priority=1,
        modality="metadata",
    ),
    "ModelCondition": DatasetInfo(
        filename="ModelCondition.csv",
        description="Conditions under which models were assayed",
        table_name="model_conditions",
        priority=1,
        modality="metadata",
    ),
    "CRISPRScreenMap": DatasetInfo(
        filename="CRISPRScreenMap.csv",
        description="Map from ModelID to ScreenIDs",
        table_name="screen_mappings",
        priority=2,
        modality="crispr",
    ),
    "AchillesScreenQCReport": DatasetInfo(
        filename="AchillesScreenQCReport.csv",
        description="Screen-level quality control metrics",
        table_name="screens",
        priority=2,
        modality="crispr",
    ),
    "CRISPRGeneEffectUncorrected": DatasetInfo(
        filename="CRISPRGeneEffectUncorrected.csv",
        description="Gene effect estimates without copy number correction",
        table_name="gene_effects_uncorrected",
        priority=2,
        modality="crispr",
    ),
    "CRISPRGeneDependency": DatasetInfo(
        filename="CRISPRGeneDependency.csv",
        description="Gene dependency probability estimates",
        table_name="gene_dependencies",
        priority=2,
        modality="crispr",
    ),
    "AchillesCommonEssentialControls": DatasetInfo(
        filename="AchillesCommonEssentialControls.csv",
        description="List of common essential control genes",
        table_name="essential_controls",
        priority=2,
        modality="crispr",
    ),
    "AchillesNonessentialControls": DatasetInfo(
        filename="AchillesNonessentialControls.csv",
        description="List of nonessential control genes",
        table_name="nonessential_controls",
        priority=2,
        modality="crispr",
    ),
    "GeneExpression": DatasetInfo(
        filename="OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
        description="Gene expression data (TPM, log-transformed)",
        table_name="gene_expression_wide",
        priority=1,
        modality="rna",
    ),
    "ProteomicsMSGygi": DatasetInfo(
        filename="harmonized_MS_CCLE_Gygi.csv",
        description=(
            "DepMap harmonized Gygi CCLE mass-spectrometry proteomics matrix. "
            "This proteomics release track can differ from the core DepMap release."
        ),
        table_name="protein_expression_ms_wide",
        priority=2,
        modality="proteomics_ms",
        release_track="proteomics",
        release_label_override="Harmonized Public Proteomics 24Q4",
    ),
    "PRISMPrimaryRepurposingExtended": DatasetInfo(
        filename="Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv",
        description=(
            "PRISM primary repurposing response matrix stored as published in "
            "compound-by-model wide form."
        ),
        table_name="drug_response_primary_wide",
        priority=2,
        modality="drug_sensitivity",
        release_track="prism_primary",
        release_label_override="PRISM Primary Repurposing DepMap Public 24Q2",
    ),
    "PRISMPrimaryRepurposingCompoundList": DatasetInfo(
        filename="Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv",
        description=(
            "Companion compound metadata for the PRISM primary repurposing matrix, "
            "including source target strings and MOA labels."
        ),
        table_name="compounds",
        priority=2,
        modality="drug_sensitivity_metadata",
        release_track="prism_primary",
        release_label_override="PRISM Primary Repurposing DepMap Public 24Q2",
    ),
    "PRISMSecondaryDoseResponseCurveParameters": DatasetInfo(
        filename="secondary-screen-dose-response-curve-parameters.csv",
        description=(
            "PRISM secondary repurposing dose-response curve parameters with auc, "
            "ec50, ic50, fit terms, and compound metadata columns."
        ),
        table_name="drug_response_secondary",
        priority=2,
        modality="drug_sensitivity",
        release_track="prism_secondary",
        release_label_override="PRISM Secondary Repurposing 20Q2",
    ),
    "OmicsSomaticMutations": DatasetInfo(
        filename="OmicsSomaticMutations.csv",
        description="Somatic mutations in cell lines",
        table_name="mutations",
        priority=3,
        modality="mutation",
    ),
    "OmicsCNGeneWGS": DatasetInfo(
        filename="OmicsCNGeneWGS.csv",
        description="Gene-level copy number data from WGS",
        table_name="copy_numbers",
        priority=3,
        modality="copy_number",
    ),
}

PHASE_1_DATASETS = [
    name for name, info in DEPMAP_FILES.items() if info.priority == 1
]
PHASE_2_DATASETS = [
    name for name, info in DEPMAP_FILES.items() if info.priority == 2
]
PHASE_3_DATASETS = [
    name for name, info in DEPMAP_FILES.items() if info.priority == 3
]

SUPPORTED_DATASETS: list[str] = list(DEPMAP_FILES.keys())
DEFAULT_DATASETS = PHASE_1_DATASETS

DEPMAP_BASE_URL = "https://depmap.org/portal/api/download/files"
DEPMAP_RELEASE_INFO_URL = (
    "https://depmap.org/portal/data_page/?tab=currentRelease"
)

SUPPORTED_FILE_EXTENSIONS = [".csv", ".txt", ".maf"]
MAX_FILE_SIZE_GB = 50

DEFAULT_BATCH_SIZE = 10000
DEFAULT_CHUNK_SIZE = 1000000


def get_dataset_release_label(dataset_name: str, default_release_label: str) -> str:
    """Return the release label that should be used for a dataset."""
    dataset = DEPMAP_FILES[dataset_name]
    return dataset.release_label_override or default_release_label
