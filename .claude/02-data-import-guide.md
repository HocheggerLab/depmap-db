# DepMap Data Import Guide

## Overview

This guide explains how to import DepMap CSV files into your local DuckDB database. It covers the complete workflow from downloading data to querying imported datasets, with special focus on implementing gene expression data import.

## Quick Start

### 1. Initialize Database

```bash
# Create the database schema
depmap-db init

# Verify initialization
depmap-db status
```

### 2. Download and Load Data

```bash
# Option A: Download from DepMap portal and load automatically
depmap-db download --load-data

# Option B: Load from local folder
depmap-db load-folder --folder /path/to/depmap_data
```

### 3. Verify Import

```bash
# Check loaded data
depmap-db status

# Run a test query
depmap-db sql --sql "SELECT COUNT(*) FROM models"
```

## Data Types and Structure

### Currently Supported Datasets

#### 1. **Model Metadata** (`Model.csv`)

**Purpose**: Cancer cell line/model information

**Structure**:
```
ModelID, PatientID, CellLineName, OncotreeLineage, OncotreePrimaryDisease,
Age, Sex, PrimaryOrMetastasis, ...
```

**Import Process**:
- Processor: `ModelProcessor`
- Target Table: `models`
- Primary Key: `model_id`
- Records: ~1,900 cancer models

**Example**:
```bash
depmap-db load-folder --folder ./data --pattern "*Model*.csv"
```

#### 2. **Gene Metadata** (`Gene.csv`)

**Purpose**: Gene identifiers and annotations

**Structure**:
```
GeneID, HugoSymbol, EntrezID, EnsemblID, GeneType, Chromosome,
StartPosition, EndPosition, Strand, Description
```

**Import Process**:
- Processor: `GeneProcessor`
- Target Table: `genes`
- Primary Key: `gene_id`
- Records: ~19,000 human genes

#### 3. **CRISPR Gene Effects** (`CRISPRGeneEffect.csv`)

**Purpose**: Gene dependency scores from CRISPR knockout screens

**Structure**: Wide format matrix
```
ModelID, GENE1 (12345), GENE2 (67890), ...
ACH-000001, -0.234, 0.123, ...
```

**Values**:
- Negative scores indicate dependency (gene is essential)
- Scores < -1.0 indicate strong dependency
- Positive scores indicate non-essential genes

**Import Options**:

**Wide Format** (default, recommended for most use cases):
```bash
depmap-db load-folder --folder ./data --wide-format
```
- One row per model, genes as columns
- Fast for model-centric queries
- Table: `gene_effects_wide`
- Efficient storage with DuckDB columnar compression

**Long Format** (for gene-centric analyses):
```bash
depmap-db load-folder --folder ./data --long-format
```
- One row per (model, gene) pair
- Fast for gene-centric queries
- Table: `gene_effects`
- Larger file size but better for aggregations

### Gene Expression Data (Implementation Needed)

#### File: `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`

**Purpose**: RNA-seq gene expression levels

**Structure**:
```
SequencingID, ModelID, IsDefaultEntryForModel, ModelConditionID,
IsDefaultEntryForMC, GENE1 (12345), GENE2 (67890), ...
```

**Metadata Columns**:
- `SequencingID`: Unique identifier for sequencing run
- `ModelID`: Links to models table (ACH-XXXXXX format)
- `IsDefaultEntryForModel`: Boolean flag for primary measurement
- `ModelConditionID`: Links to experimental conditions
- `IsDefaultEntryForMC`: Boolean flag for default condition

**Gene Columns**: ~18,000 genes in format `"SYMBOL (EntrezID)"`

**Values**:
- log(TPM + 1) transformed expression values
- Range: typically 0-15 (0 = not expressed, >10 = highly expressed)
- TPM = Transcripts Per Million (normalized read count)

## Import Implementation

### For Currently Supported Data

All currently supported data types work out-of-the-box:

```bash
# Step 1: Initialize database
depmap-db init

# Step 2: Load all data from folder
depmap-db load-folder --folder /Users/hh65/depmap_data

# The CLI automatically detects:
# - Model.csv → ModelProcessor → models table
# - Gene.csv → GeneProcessor → genes table
# - CRISPRGeneEffect.csv → GeneEffectWideProcessor → gene_effects_wide table
```

### For Gene Expression Data (Implementation Guide)

Gene expression data requires implementing two new processors. Here's the complete implementation strategy:

#### Step 1: Update Database Schema

Add to `src/depmap_db/database/schema.py`:

```python
# Add after gene_effects_wide table definition (line 171)

# Gene expression long format
self._add_table(
    TableSchema(
        name="gene_expression",
        table_type=TableType.CORE_DATA,
        dependencies=["models", "model_conditions"],
        description="Gene expression TPM values in long format",
        sql="""
        CREATE TABLE IF NOT EXISTS gene_expression (
            model_id VARCHAR NOT NULL,
            gene_id VARCHAR NOT NULL,
            expression_value DOUBLE NOT NULL,
            sequencing_id VARCHAR,
            model_condition_id VARCHAR,
            is_default_entry BOOLEAN DEFAULT TRUE,
            is_default_for_mc BOOLEAN DEFAULT TRUE,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            PRIMARY KEY (model_id, gene_id),
            FOREIGN KEY (model_id) REFERENCES models(model_id),
            FOREIGN KEY (model_condition_id) REFERENCES model_conditions(model_condition_id)
        )
        """,
    )
)

# Gene expression wide format
self._add_table(
    TableSchema(
        name="gene_expression_wide",
        table_type=TableType.CORE_DATA,
        dependencies=["models", "model_conditions"],
        description="Gene expression TPM values in wide format (models × genes matrix)",
        sql="""
        CREATE TABLE IF NOT EXISTS gene_expression_wide (
            model_id VARCHAR PRIMARY KEY,
            sequencing_id VARCHAR,
            model_condition_id VARCHAR,
            is_default_entry BOOLEAN DEFAULT TRUE,
            is_default_for_mc BOOLEAN DEFAULT TRUE,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            FOREIGN KEY (model_id) REFERENCES models(model_id),
            FOREIGN KEY (model_condition_id) REFERENCES model_conditions(model_condition_id)
            -- Gene columns will be added dynamically during data loading
        )
        """,
    )
)
```

#### Step 2: Create Gene Expression Processor (Long Format)

Create `src/depmap_db/etl/processors/gene_expression.py`:

```python
"""Processor for gene expression TPM data (long format)."""

import re
from typing import List

import numpy as np
import pandas as pd

from .base import BaseProcessor


class GeneExpressionProcessor(BaseProcessor):
    """Processor for gene expression data in long format."""

    def __init__(self, batch_size: int = None):
        """Initialize gene expression processor."""
        super().__init__("GeneExpression", batch_size)

    def get_table_name(self) -> str:
        """Get the target database table name."""
        return "gene_expression"

    def validate_data(self, df: pd.DataFrame) -> tuple[pd.DataFrame, List[str]]:
        """Validate and clean gene expression data.

        Args:
            df: Input DataFrame with models as rows and genes as columns

        Returns:
            Tuple of (cleaned_dataframe, warnings_list)
        """
        warnings = []
        original_shape = df.shape

        # Check for required metadata columns
        required_cols = ['ModelID', 'SequencingID', 'ModelConditionID']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            warnings.append(f"Missing required columns: {missing_cols}")

        # Validate ModelID format
        if 'ModelID' in df.columns:
            invalid_models = df[~df['ModelID'].str.startswith('ACH-', na=False)]['ModelID'].tolist()
            if invalid_models:
                warnings.append(f"Found {len(invalid_models)} invalid model IDs")
                df = df[df['ModelID'].str.startswith('ACH-', na=False)]

        # Identify gene columns (those with parentheses indicating Entrez IDs)
        gene_columns = [col for col in df.columns
                       if isinstance(col, str) and '(' in col and ')' in col]

        if not gene_columns:
            warnings.append("No gene columns found in data")
            return df, warnings

        # Check expression value ranges (log(TPM+1) should be >= 0)
        for col in gene_columns[:10]:  # Sample first 10 genes
            if col in df.columns:
                min_val = df[col].min()
                max_val = df[col].max()

                if min_val < 0:
                    warnings.append(f"Negative expression values found in {col}: min={min_val:.3f}")
                if max_val > 20:
                    warnings.append(f"Unusually high expression in {col}: max={max_val:.3f}")

        # Handle missing values
        missing_count = df[gene_columns].isna().sum().sum()
        if missing_count > 0:
            warnings.append(f"Found {missing_count:,} missing expression values")

        # Handle infinite values
        inf_count = 0
        for col in gene_columns:
            if df[col].dtype in [np.float64, np.float32]:
                inf_mask = np.isinf(df[col])
                if inf_mask.any():
                    inf_count += inf_mask.sum()
                    df.loc[inf_mask, col] = np.nan

        if inf_count > 0:
            warnings.append(f"Replaced {inf_count:,} infinite values with NaN")

        self.logger.info(f"Validation: {original_shape} -> {df.shape}, {len(warnings)} warnings")
        return df, warnings

    def transform_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """Transform wide-format expression data to long format.

        Args:
            df: Wide-format DataFrame (models × genes)

        Returns:
            Long-format DataFrame with columns: model_id, gene_id, expression_value, etc.
        """
        self.logger.info("Transforming wide format to long format...")

        # Identify metadata and gene columns
        metadata_cols = ['ModelID', 'SequencingID', 'ModelConditionID',
                        'IsDefaultEntryForModel', 'IsDefaultEntryForMC']
        metadata_cols = [col for col in metadata_cols if col in df.columns]

        gene_columns = [col for col in df.columns
                       if col not in metadata_cols and '(' in col]

        # Clean gene names - remove Entrez IDs in parentheses
        gene_name_map = {}
        for col in gene_columns:
            clean_name = re.sub(r'\s*\([^)]*\)', '', col).strip()
            gene_name_map[col] = clean_name

        # Melt to long format
        df_long = pd.melt(
            df,
            id_vars=metadata_cols,
            value_vars=gene_columns,
            var_name='gene_raw',
            value_name='expression_value'
        )

        # Map to clean gene names
        df_long['gene_id'] = df_long['gene_raw'].map(gene_name_map)
        df_long = df_long.drop('gene_raw', axis=1)

        # Remove rows with missing expression values
        initial_count = len(df_long)
        df_long = df_long.dropna(subset=['expression_value'])
        final_count = len(df_long)

        self.logger.info(f"Removed {initial_count - final_count:,} records with missing values")

        # Rename columns to match schema
        column_mapping = {
            'ModelID': 'model_id',
            'SequencingID': 'sequencing_id',
            'ModelConditionID': 'model_condition_id',
            'IsDefaultEntryForModel': 'is_default_entry',
            'IsDefaultEntryForMC': 'is_default_for_mc'
        }
        df_long = df_long.rename(columns=column_mapping)

        # Ensure proper data types
        df_long['model_id'] = df_long['model_id'].astype(str)
        df_long['gene_id'] = df_long['gene_id'].astype(str)
        df_long['expression_value'] = df_long['expression_value'].astype(float)

        if 'is_default_entry' in df_long.columns:
            df_long['is_default_entry'] = df_long['is_default_entry'].astype(bool)
        if 'is_default_for_mc' in df_long.columns:
            df_long['is_default_for_mc'] = df_long['is_default_for_mc'].astype(bool)

        # Add timestamp
        df_long['created_at'] = pd.Timestamp.now()

        # Sort for better insert performance
        df_long = df_long.sort_values(['model_id', 'gene_id'])

        self.logger.info(f"Transformed to {len(df_long):,} expression records")

        return df_long

    def get_existing_records(self) -> set:
        """Get set of existing (model_id, gene_id) pairs."""
        table_name = self.get_table_name()

        try:
            if not self.db_manager.table_exists(table_name):
                return set()

            result = self.db_manager.execute(
                f"SELECT model_id, gene_id FROM {table_name}"
            )
            existing = {(row[0], row[1]) for row in result.fetchall()}

            self.logger.info(f"Found {len(existing):,} existing expression records")
            return existing

        except Exception as e:
            self.logger.warning(f"Could not retrieve existing records: {e}")
            return set()
```

#### Step 3: Create Gene Expression Wide Processor

Create `src/depmap_db/etl/processors/gene_expression_wide.py`:

```python
"""Processor for gene expression TPM data (wide format)."""

import re
from typing import List

import pandas as pd

from .base import BaseProcessor


class GeneExpressionWideProcessor(BaseProcessor):
    """Processor for gene expression data in wide format."""

    def __init__(self, batch_size: int = None):
        """Initialize gene expression wide processor."""
        super().__init__("GeneExpression", batch_size)

    def get_table_name(self) -> str:
        """Get the target database table name."""
        return "gene_expression_wide"

    def validate_data(self, df: pd.DataFrame) -> tuple[pd.DataFrame, List[str]]:
        """Validate gene expression data (same as long format processor)."""
        # Reuse validation from GeneExpressionProcessor
        from .gene_expression import GeneExpressionProcessor
        temp_processor = GeneExpressionProcessor()
        return temp_processor.validate_data(df)

    def transform_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """Transform data for wide format storage.

        Args:
            df: Wide-format DataFrame (already in correct format)

        Returns:
            DataFrame with cleaned column names and proper types
        """
        self.logger.info("Preparing wide format data...")

        # Identify metadata columns
        metadata_cols = ['ModelID', 'SequencingID', 'ModelConditionID',
                        'IsDefaultEntryForModel', 'IsDefaultEntryForMC']

        # Clean gene column names - remove Entrez IDs
        new_columns = {}
        for col in df.columns:
            if col in metadata_cols:
                # Rename metadata columns
                clean_col = col
                if col == 'ModelID':
                    clean_col = 'model_id'
                elif col == 'SequencingID':
                    clean_col = 'sequencing_id'
                elif col == 'ModelConditionID':
                    clean_col = 'model_condition_id'
                elif col == 'IsDefaultEntryForModel':
                    clean_col = 'is_default_entry'
                elif col == 'IsDefaultEntryForMC':
                    clean_col = 'is_default_for_mc'
                new_columns[col] = clean_col
            elif '(' in col:  # Gene column
                # Remove Entrez ID in parentheses
                clean_col = re.sub(r'\s*\([^)]*\)', '', col).strip()
                new_columns[col] = clean_col

        df = df.rename(columns=new_columns)

        # Add timestamp
        df['created_at'] = pd.Timestamp.now()

        self.logger.info(f"Prepared {len(df)} models with {len(df.columns)-6} genes")

        return df

    def insert_batch(self, df: pd.DataFrame) -> tuple[int, int]:
        """Insert batch using DuckDB's efficient bulk insert."""
        if df.empty:
            return 0, 0

        table_name = self.get_table_name()

        try:
            # Check if table exists and has the correct schema
            if not self.db_manager.table_exists(table_name):
                # Create table from DataFrame structure
                metadata_cols = ['model_id', 'sequencing_id', 'model_condition_id',
                               'is_default_entry', 'is_default_for_mc', 'created_at']
                gene_cols = [col for col in df.columns if col not in metadata_cols]

                # Create table with metadata columns
                create_sql = f"""
                CREATE TABLE {table_name} (
                    model_id VARCHAR PRIMARY KEY,
                    sequencing_id VARCHAR,
                    model_condition_id VARCHAR,
                    is_default_entry BOOLEAN,
                    is_default_for_mc BOOLEAN,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
                """
                self.db_manager.execute(create_sql)

                # Add gene columns dynamically
                for gene in gene_cols:
                    safe_gene = gene.replace("'", "''")
                    self.db_manager.execute(
                        f"ALTER TABLE {table_name} ADD COLUMN \"{safe_gene}\" DOUBLE"
                    )

            # Insert using DuckDB's DataFrame support
            self.db_manager.execute(f"""
                INSERT OR REPLACE INTO {table_name}
                SELECT * FROM df
            """)

            return len(df), 0

        except Exception as e:
            self.logger.error(f"Failed to insert batch: {e}")
            raise
```

#### Step 4: Register Processors in CLI

Update `src/depmap_db/cli.py`:

```python
# Add to imports (line 17-22)
from .etl.processors import (
    GeneEffectProcessor,
    GeneEffectWideProcessor,
    GeneExpressionProcessor,        # Add this
    GeneExpressionWideProcessor,    # Add this
    GeneProcessor,
    ModelProcessor,
)

# Update _get_processor function (line 32-43)
def _get_processor(dataset_name: str, wide_format: bool = True):
    """Get the appropriate processor for a dataset."""
    processors = {
        "CRISPRGeneEffect": GeneEffectWideProcessor if wide_format else GeneEffectProcessor,
        "GeneExpression": GeneExpressionWideProcessor if wide_format else GeneExpressionProcessor,  # Add this
        "Gene": GeneProcessor,
        "Model": ModelProcessor,
    }

    processor_class = processors.get(dataset_name)
    return processor_class() if processor_class else None

# Update _detect_dataset_from_filename function (line 46-68)
def _detect_dataset_from_filename(filename: str) -> Optional[str]:
    """Detect dataset type from filename."""
    filename_lower = filename.lower()

    if "crisprgeneeffect" in filename_lower or "crispr_gene_effect" in filename_lower:
        return "CRISPRGeneEffect"
    elif "omicsexpression" in filename_lower or "expression" in filename_lower:  # Add this
        return "GeneExpression"                                                   # Add this
    elif ("model" in filename_lower and "data" in filename_lower) or filename_lower.startswith("model"):
        return "Model"
    elif filename_lower.startswith("gene") and filename_lower.endswith(".csv"):
        return "Gene"
    elif "modelcondition" in filename_lower or "model_condition" in filename_lower:
        return "ModelCondition"

    return None
```

#### Step 5: Update Constants

Update `src/depmap_db/utils/constants.py`:

```python
# Add to SUPPORTED_DATASETS list
SUPPORTED_DATASETS = [
    "Model",
    "Gene",
    "CRISPRGeneEffect",
    "GeneExpression",  # Add this
]

# Optionally add to PHASE_1_DATASETS if you want it auto-downloaded
PHASE_1_DATASETS = [
    "Model",
    "Gene",
    "CRISPRGeneEffect",
    # "GeneExpression",  # Uncomment when ready for automatic downloads
]
```

### Usage After Implementation

Once implemented, importing gene expression data becomes as simple as:

```bash
# Initialize database (if not done already)
depmap-db init

# Load gene expression data
depmap-db load-folder --folder /Users/hh65/depmap_data \
    --pattern "*OmicsExpression*.csv" \
    --wide-format

# Verify import
depmap-db status
depmap-db sql --sql "SELECT COUNT(*) FROM gene_expression_wide"
```

## Data Import Best Practices

### 1. Import Order

Always import in this order to satisfy foreign key constraints:

```bash
# 1. Core metadata
depmap-db load-folder --folder ./data --pattern "*Model*.csv"
depmap-db load-folder --folder ./data --pattern "*Gene*.csv"

# 2. Assay data
depmap-db load-folder --folder ./data --pattern "*CRISPRGeneEffect*.csv"
depmap-db load-folder --folder ./data --pattern "*OmicsExpression*.csv"
```

Or import all at once (CLI handles dependencies):
```bash
depmap-db load-folder --folder ./data
```

### 2. Performance Optimization

For large datasets:

```bash
# Increase batch size for faster imports (uses more memory)
export DEPMAP_ETL__BATCH_SIZE=50000
depmap-db load-folder --folder ./data

# Use more memory for DuckDB
export DEPMAP_DATABASE__MAX_MEMORY=16GB
depmap-db load-folder --folder ./data
```

### 3. Error Handling

```bash
# Force reload if import fails partway
depmap-db load-folder --folder ./data --force

# Check for errors in specific file
depmap-db load-folder --folder ./data --pattern "*problem_file*.csv" --debug
```

### 4. Data Validation

After import, run validation queries:

```sql
-- Check for orphaned records
SELECT COUNT(*) FROM gene_expression_wide WHERE model_id NOT IN (SELECT model_id FROM models);

-- Check data completeness
SELECT
    'models' as table_name, COUNT(*) as count FROM models
UNION ALL
SELECT 'genes', COUNT(*) FROM genes
UNION ALL
SELECT 'gene_effects_wide', COUNT(*) FROM gene_effects_wide
UNION ALL
SELECT 'gene_expression_wide', COUNT(*) FROM gene_expression_wide;

-- Sample data check
SELECT * FROM gene_expression_wide LIMIT 5;
```

## Querying Imported Data

### Basic Queries

```sql
-- Get expression for a specific gene across all models
SELECT m.cell_line_name, m.oncotree_lineage, ge.expression_value
FROM gene_expression ge
JOIN models m ON ge.model_id = m.model_id
WHERE ge.gene_id = 'TP53'
ORDER BY ge.expression_value DESC;

-- Compare CRISPR dependency and expression
SELECT
    m.cell_line_name,
    gef.gene_effect as crispr_score,
    gex.expression_value as expression_tpm
FROM models m
JOIN gene_effects gef ON m.model_id = gef.model_id
JOIN gene_expression gex ON m.model_id = gex.model_id AND gef.gene_id = gex.gene_id
WHERE gef.gene_id = 'MYC'
  AND gef.gene_effect < -0.5;
```

### Advanced Queries

```sql
-- Find genes with high expression but low dependency (overexpressed non-essential)
WITH gene_stats AS (
    SELECT
        ge.gene_id,
        AVG(ge.gene_effect) as avg_dependency,
        AVG(gx.expression_value) as avg_expression
    FROM gene_effects ge
    JOIN gene_expression gx ON ge.model_id = gx.model_id AND ge.gene_id = gx.gene_id
    GROUP BY ge.gene_id
)
SELECT gene_id, avg_dependency, avg_expression
FROM gene_stats
WHERE avg_expression > 5  -- High expression (log scale)
  AND avg_dependency > 0   -- Not essential
ORDER BY avg_expression DESC
LIMIT 50;

-- Tissue-specific gene expression
SELECT
    m.oncotree_lineage,
    AVG(gx.expression_value) as mean_expression,
    STDDEV(gx.expression_value) as stddev_expression
FROM gene_expression gx
JOIN models m ON gx.model_id = m.model_id
WHERE gx.gene_id = 'ALB'  -- Albumin (liver-specific)
GROUP BY m.oncotree_lineage
ORDER BY mean_expression DESC;
```

## Troubleshooting

### Common Issues

#### Issue: "Table already exists"

**Solution**: Use `--force` flag to reload
```bash
depmap-db load-folder --folder ./data --force
```

#### Issue: "Foreign key constraint failed"

**Cause**: Trying to import assay data before model metadata

**Solution**: Import models first
```bash
depmap-db load-folder --folder ./data --pattern "*Model*.csv"
depmap-db load-folder --folder ./data --pattern "*CRISPRGeneEffect*.csv"
```

#### Issue: "Out of memory"

**Solution**: Reduce batch size or increase memory limit
```bash
export DEPMAP_ETL__BATCH_SIZE=5000
export DEPMAP_DATABASE__MAX_MEMORY=16GB
depmap-db load-folder --folder ./data
```

#### Issue: "Cannot detect dataset type"

**Cause**: Filename doesn't match expected patterns

**Solution**: Rename file or specify processor manually (requires code modification)

### Data Integrity Checks

```sql
-- Check for null values in key columns
SELECT COUNT(*) FROM gene_expression WHERE model_id IS NULL;
SELECT COUNT(*) FROM gene_expression WHERE gene_id IS NULL;
SELECT COUNT(*) FROM gene_expression WHERE expression_value IS NULL;

-- Check for duplicate entries
SELECT model_id, gene_id, COUNT(*) as count
FROM gene_expression
GROUP BY model_id, gene_id
HAVING COUNT(*) > 1;

-- Verify value ranges
SELECT
    MIN(expression_value) as min_expr,
    MAX(expression_value) as max_expr,
    AVG(expression_value) as avg_expr
FROM gene_expression;
```

## Summary

### Current Workflow (Supported Data)

```bash
# 1. Initialize
depmap-db init

# 2. Load data
depmap-db load-folder --folder /Users/hh65/depmap_data

# 3. Query
depmap-db sql --sql "SELECT * FROM models LIMIT 10"
```

### Future Workflow (With Gene Expression)

```bash
# 1. Initialize
depmap-db init

# 2. Load all data including expression
depmap-db load-folder --folder /Users/hh65/depmap_data --wide-format

# 3. Query integrated data
depmap-db sql --sql "
    SELECT m.cell_line_name,
           gef.TP53 as tp53_dependency,
           gex.TP53 as tp53_expression
    FROM models m
    JOIN gene_effects_wide gef ON m.model_id = gef.model_id
    JOIN gene_expression_wide gex ON m.model_id = gex.model_id
    LIMIT 100
"
```

The implementation guide above provides everything needed to add gene expression support following the established patterns in the codebase.
