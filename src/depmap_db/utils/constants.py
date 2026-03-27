"""Constants and configuration data for DepMap datasets."""

from typing import NamedTuple


class DatasetInfo(NamedTuple):
    """Information about a DepMap dataset file."""

    filename: str
    description: str
    table_name: str
    priority: int


DEPMAP_FILES: dict[str, DatasetInfo] = {
    "CRISPRGeneEffect": DatasetInfo(
        filename="CRISPRGeneEffect.csv",
        description="Gene effect estimates for all models, integrated using Chronos",
        table_name="gene_effects_wide",
        priority=1,
    ),
    "Model": DatasetInfo(
        filename="Model.csv",
        description="Metadata describing all cancer models/cell lines",
        table_name="models",
        priority=1,
    ),
    "Gene": DatasetInfo(
        filename="Gene.csv",
        description="Gene metadata and identifiers",
        table_name="genes",
        priority=1,
    ),
    "ModelCondition": DatasetInfo(
        filename="ModelCondition.csv",
        description="Conditions under which models were assayed",
        table_name="model_conditions",
        priority=1,
    ),
    "CRISPRScreenMap": DatasetInfo(
        filename="CRISPRScreenMap.csv",
        description="Map from ModelID to ScreenIDs",
        table_name="screen_mappings",
        priority=2,
    ),
    "AchillesScreenQCReport": DatasetInfo(
        filename="AchillesScreenQCReport.csv",
        description="Screen-level quality control metrics",
        table_name="screens",
        priority=2,
    ),
    "CRISPRGeneEffectUncorrected": DatasetInfo(
        filename="CRISPRGeneEffectUncorrected.csv",
        description="Gene effect estimates without copy number correction",
        table_name="gene_effects_uncorrected",
        priority=2,
    ),
    "CRISPRGeneDependency": DatasetInfo(
        filename="CRISPRGeneDependency.csv",
        description="Gene dependency probability estimates",
        table_name="gene_dependencies",
        priority=2,
    ),
    "AchillesCommonEssentialControls": DatasetInfo(
        filename="AchillesCommonEssentialControls.csv",
        description="List of common essential control genes",
        table_name="essential_controls",
        priority=2,
    ),
    "AchillesNonessentialControls": DatasetInfo(
        filename="AchillesNonessentialControls.csv",
        description="List of nonessential control genes",
        table_name="nonessential_controls",
        priority=2,
    ),
    "GeneExpression": DatasetInfo(
        filename="OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
        description="Gene expression data (TPM, log-transformed)",
        table_name="gene_expression_wide",
        priority=2,
    ),
    "OmicsSomaticMutations": DatasetInfo(
        filename="OmicsSomaticMutations.csv",
        description="Somatic mutations in cell lines",
        table_name="mutations",
        priority=3,
    ),
    "OmicsCNGeneWGS": DatasetInfo(
        filename="OmicsCNGeneWGS.csv",
        description="Gene-level copy number data from WGS",
        table_name="copy_numbers",
        priority=3,
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

DEPMAP_BASE_URL = "https://depmap.org/portal/api/download/file/"
DEPMAP_RELEASE_INFO_URL = (
    "https://depmap.org/portal/data_page/?tab=currentRelease"
)

SUPPORTED_FILE_EXTENSIONS = [".csv", ".txt", ".maf"]
MAX_FILE_SIZE_GB = 50

DEFAULT_BATCH_SIZE = 10000
DEFAULT_CHUNK_SIZE = 1000000
