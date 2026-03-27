"""Constants and configuration data for DepMap datasets."""

from typing import NamedTuple


class DatasetInfo(NamedTuple):
    """Information about a DepMap dataset file."""

    filename: str
    description: str
    table_name: str
    priority: int  # Lower numbers indicate higher priority for initial implementation


# Mapping of dataset names to file information
DEPMAP_FILES: dict[str, DatasetInfo] = {
    # Core gene effect data (Phase 1 priority)
    "CRISPRGeneEffect": DatasetInfo(
        filename="CRISPRGeneEffect.csv",
        description="Gene effect estimates for all models, integrated using Chronos",
        table_name="gene_effects",
        priority=1,
    ),
    # Essential metadata (Phase 1 priority)
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
    # Screen data (Phase 2 priority)
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
    # Additional gene effect data (Phase 2 priority)
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
    # Control gene lists (Phase 2 priority)
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
    # Expression data (Phase 2 priority - now implemented)
    "GeneExpression": DatasetInfo(
        filename="OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
        description="Gene expression data (TPM, log-transformed)",
        table_name="gene_expression",
        priority=2,
    ),
    # Mutation data (Phase 3 priority)
    "OmicsSomaticMutations": DatasetInfo(
        filename="OmicsSomaticMutations.csv",
        description="Somatic mutations in cell lines",
        table_name="mutations",
        priority=3,
    ),
    # Copy number data (Phase 3 priority)
    "OmicsCNGeneWGS": DatasetInfo(
        filename="OmicsCNGeneWGS.csv",
        description="Gene-level copy number data from WGS",
        table_name="copy_numbers",
        priority=3,
    ),
}

# Datasets grouped by implementation phase
PHASE_1_DATASETS = [
    name for name, info in DEPMAP_FILES.items() if info.priority == 1
]

PHASE_2_DATASETS = [
    name for name, info in DEPMAP_FILES.items() if info.priority == 2
]

PHASE_3_DATASETS = [
    name for name, info in DEPMAP_FILES.items() if info.priority == 3
]

# All supported dataset names
SUPPORTED_DATASETS: list[str] = list(DEPMAP_FILES.keys())

# Default datasets for initial download
DEFAULT_DATASETS = PHASE_1_DATASETS

# DepMap API endpoints and URLs
DEPMAP_BASE_URL = "https://depmap.org/portal/api/download/file/"
DEPMAP_RELEASE_INFO_URL = (
    "https://depmap.org/portal/data_page/?tab=currentRelease"
)

# File validation settings
SUPPORTED_FILE_EXTENSIONS = [".csv", ".txt", ".maf"]
MAX_FILE_SIZE_GB = 50  # Maximum expected file size

# Database settings
DEFAULT_BATCH_SIZE = 10000
DEFAULT_CHUNK_SIZE = 1000000  # For large file processing
