"""Export functions for generating CSV files from the database."""

from pathlib import Path
from typing import Any

from ..config import get_logger
from ..database.connection import get_db_manager

logger = get_logger(__name__)


def export_gene_expression_csv(
    output_path: Path,
    metadata_columns: list[str] | None = None,
    models_filter: str | None = None,
) -> int:
    """Export gene expression data with model metadata to CSV.

    This function creates a CSV file in wide format (models x genes) with
    metadata columns matching the DepMap portal format.

    Args:
        output_path: Path where the CSV file should be saved
        metadata_columns: Optional list of metadata columns from models table to include.
                        If None, uses default columns that match gene effects export.
        models_filter: Optional SQL WHERE clause to filter models (e.g., "oncotree_lineage = 'Lung'")

    Returns:
        Number of models exported

    Example:
        >>> from pathlib import Path
        >>> from depmap_db.export import export_gene_expression_csv
        >>>
        >>> # Export all models
        >>> export_gene_expression_csv(Path("gene_expression_export.csv"))
        >>>
        >>> # Export only lung cancer models
        >>> export_gene_expression_csv(
        ...     Path("lung_expression.csv"),
        ...     models_filter="oncotree_lineage = 'Lung'"
        ... )
    """
    db_manager = get_db_manager()

    # Default metadata columns that match gene effects export format
    if metadata_columns is None:
        metadata_columns = [
            "model_id",
            "cell_line_name",
            "stripped_cell_line_name",
            "oncotree_lineage",
            "oncotree_primary_disease",
            "oncotree_subtype",
            "oncotree_code",
            "ccle_name",
            "depmap_model_type",
        ]

    logger.info("Exporting gene expression data to %s", output_path)
    logger.info("Including %s metadata columns", len(metadata_columns))

    # Build the query
    # Get gene columns from gene_expression_wide (excluding metadata)
    gene_expression_cols_query = """
    SELECT column_name
    FROM information_schema.columns
    WHERE table_name = 'gene_expression_wide'
      AND column_name NOT IN ('model_id', 'sequencing_id', 'model_condition_id',
                               'is_default_entry', 'is_default_for_mc', 'created_at')
    ORDER BY column_name
    """

    gene_cols_result = db_manager.execute(gene_expression_cols_query)
    gene_columns = [row[0] for row in gene_cols_result.fetchall()]

    logger.info("Found %s gene columns", len(gene_columns))

    # Build SELECT clause for metadata
    metadata_select = ", ".join([f"m.{col}" for col in metadata_columns])

    # Build SELECT clause for gene expression values
    gene_select = ", ".join([f'gex."{col}"' for col in gene_columns])

    # Build WHERE clause
    where_clause = f"WHERE {models_filter}" if models_filter else ""

    # Full query
    export_query = f"""
    SELECT
        {metadata_select},
        {gene_select}
    FROM models m
    INNER JOIN gene_expression_wide gex ON m.model_id = gex.model_id
    {where_clause}
    ORDER BY m.cell_line_name
    """

    logger.info("Executing export query using optimized COPY command...")

    # Use DuckDB's COPY TO command for memory-efficient export
    # This writes directly to CSV without loading everything into memory
    output_path.parent.mkdir(parents=True, exist_ok=True)

    copy_query = f"""
    COPY (
        {export_query}
    ) TO '{output_path}' (HEADER, DELIMITER ',')
    """

    db_manager.execute(copy_query)

    # Count the exported records
    count_query = f"""
    SELECT COUNT(*) as count FROM (
        {export_query}
    )
    """
    result = db_manager.execute(count_query)
    num_records: int = result.fetchone()[0]

    file_size = output_path.stat().st_size
    logger.info(
        "Successfully exported %s models to %s (%s bytes)",
        num_records,
        output_path,
        file_size,
    )

    return num_records


def export_gene_effects_csv(
    output_path: Path,
    metadata_columns: list[str] | None = None,
    models_filter: str | None = None,
) -> int:
    """Export CRISPR gene effect data with model metadata to CSV.

    This function creates a CSV file in wide format (models x genes) with
    metadata columns matching the DepMap portal format.

    Args:
        output_path: Path where the CSV file should be saved
        metadata_columns: Optional list of metadata columns from models table to include.
                        If None, uses default columns.
        models_filter: Optional SQL WHERE clause to filter models (e.g., "oncotree_lineage = 'Lung'")

    Returns:
        Number of models exported

    Example:
        >>> from pathlib import Path
        >>> from depmap_db.export import export_gene_effects_csv
        >>>
        >>> # Export all models
        >>> export_gene_effects_csv(Path("gene_effects_export.csv"))
        >>>
        >>> # Export only breast cancer models
        >>> export_gene_effects_csv(
        ...     Path("breast_crispr.csv"),
        ...     models_filter="oncotree_lineage = 'Breast'"
        ... )
    """
    db_manager = get_db_manager()

    # Default metadata columns
    if metadata_columns is None:
        metadata_columns = [
            "model_id",
            "cell_line_name",
            "stripped_cell_line_name",
            "oncotree_lineage",
            "oncotree_primary_disease",
            "oncotree_subtype",
            "oncotree_code",
            "ccle_name",
            "depmap_model_type",
        ]

    logger.info("Exporting CRISPR gene effect data to %s", output_path)
    logger.info("Including %s metadata columns", len(metadata_columns))

    # Get gene columns from gene_effects_wide (excluding metadata)
    gene_effects_cols_query = """
    SELECT column_name
    FROM information_schema.columns
    WHERE table_name = 'gene_effects_wide'
      AND column_name NOT IN ('model_id', 'created_at')
    ORDER BY column_name
    """

    gene_cols_result = db_manager.execute(gene_effects_cols_query)
    gene_columns = [row[0] for row in gene_cols_result.fetchall()]

    logger.info("Found %s gene columns", len(gene_columns))

    # Build SELECT clause for metadata
    metadata_select = ", ".join([f"m.{col}" for col in metadata_columns])

    # Build SELECT clause for gene effect values
    gene_select = ", ".join([f'gef."{col}"' for col in gene_columns])

    # Build WHERE clause
    where_clause = f"WHERE {models_filter}" if models_filter else ""

    # Full query
    export_query = f"""
    SELECT
        {metadata_select},
        {gene_select}
    FROM models m
    INNER JOIN gene_effects_wide gef ON m.model_id = gef.model_id
    {where_clause}
    ORDER BY m.cell_line_name
    """

    logger.info("Executing export query using optimized COPY command...")

    # Use DuckDB's COPY TO command for memory-efficient export
    # This writes directly to CSV without loading everything into memory
    output_path.parent.mkdir(parents=True, exist_ok=True)

    copy_query = f"""
    COPY (
        {export_query}
    ) TO '{output_path}' (HEADER, DELIMITER ',')
    """

    db_manager.execute(copy_query)

    # Count the exported records
    count_query = f"""
    SELECT COUNT(*) as count FROM (
        {export_query}
    )
    """
    result = db_manager.execute(count_query)
    num_records: int = result.fetchone()[0]

    file_size = output_path.stat().st_size
    logger.info(
        "Successfully exported %s models to %s (%s bytes)",
        num_records,
        output_path,
        file_size,
    )

    return num_records


def export_integrated_data_csv(
    output_path: Path,
    metadata_columns: list[str] | None = None,
    genes_to_include: list[str] | None = None,
    models_filter: str | None = None,
) -> int:
    """Export integrated gene expression and CRISPR gene effect data.

    This function creates a CSV with both gene expression and gene dependency
    data side-by-side for the same models and genes.

    Args:
        output_path: Path where the CSV file should be saved
        metadata_columns: Optional list of metadata columns from models table
        genes_to_include: Optional list of specific genes to include. If None, includes all genes.
        models_filter: Optional SQL WHERE clause to filter models

    Returns:
        Number of models exported

    Example:
        >>> from pathlib import Path
        >>> from depmap_db.export import export_integrated_data_csv
        >>>
        >>> # Export TP53 and MYC for all models
        >>> export_integrated_data_csv(
        ...     Path("tp53_myc_integrated.csv"),
        ...     genes_to_include=["TP53", "MYC"]
        ... )
    """
    db_manager = get_db_manager()

    # Default metadata columns
    if metadata_columns is None:
        metadata_columns = [
            "model_id",
            "cell_line_name",
            "oncotree_lineage",
            "oncotree_primary_disease",
            "oncotree_subtype",
        ]

    logger.info("Exporting integrated data to %s", output_path)

    # Determine which genes to include
    if genes_to_include:
        genes_list = genes_to_include
        logger.info("Including %s specified genes", len(genes_list))
    else:
        # Get genes that are in both tables
        common_genes_query = """
        SELECT gex_cols.column_name
        FROM information_schema.columns gex_cols
        WHERE gex_cols.table_name = 'gene_expression_wide'
          AND gex_cols.column_name NOT IN ('model_id', 'sequencing_id', 'model_condition_id',
                                           'is_default_entry', 'is_default_for_mc', 'created_at')
          AND EXISTS (
              SELECT 1
              FROM information_schema.columns gef_cols
              WHERE gef_cols.table_name = 'gene_effects_wide'
                AND gef_cols.column_name = gex_cols.column_name
          )
        ORDER BY gex_cols.column_name
        """
        result = db_manager.execute(common_genes_query)
        genes_list = [row[0] for row in result.fetchall()]
        logger.info("Found %s genes in both datasets", len(genes_list))

    if not genes_list:
        logger.error("No genes to export")
        return 0

    # Build SELECT clauses
    metadata_select = ", ".join([f"m.{col}" for col in metadata_columns])

    # For each gene, select both expression and effect
    gene_selects = []
    for gene in genes_list:
        gene_selects.append(f'gex."{gene}" as {gene}_expression')
        gene_selects.append(f'gef."{gene}" as {gene}_dependency')

    gene_select = ", ".join(gene_selects)

    # Build WHERE clause
    where_clause = f"WHERE {models_filter}" if models_filter else ""

    # Full query with both joins
    export_query = f"""
    SELECT
        {metadata_select},
        {gene_select}
    FROM models m
    LEFT JOIN gene_expression_wide gex ON m.model_id = gex.model_id
    LEFT JOIN gene_effects_wide gef ON m.model_id = gef.model_id
    {where_clause}
    ORDER BY m.cell_line_name
    """

    logger.info("Executing export query using optimized COPY command...")

    # Use DuckDB's COPY TO command for memory-efficient export
    # This writes directly to CSV without loading everything into memory
    output_path.parent.mkdir(parents=True, exist_ok=True)

    copy_query = f"""
    COPY (
        {export_query}
    ) TO '{output_path}' (HEADER, DELIMITER ',')
    """

    db_manager.execute(copy_query)

    # Count the exported records
    count_query = f"""
    SELECT COUNT(*) as count FROM (
        {export_query}
    )
    """
    result = db_manager.execute(count_query)
    num_records: int = result.fetchone()[0]

    file_size = output_path.stat().st_size
    logger.info(
        "Successfully exported %s models with %s genes (expression + dependency) to %s (%s bytes)",
        num_records,
        len(genes_list),
        output_path,
        file_size,
    )

    return num_records


def export_matched_expression_csv(
    output_path: Path,
    metadata_columns: list[str] | None = None,
    models_filter: str | None = None,
    report_missing: bool = True,
) -> dict[str, Any]:
    """Export gene expression data matched to CRISPR gene effect data.

    This function exports only models and genes that are present in the CRISPR
    gene_effects_wide table. If data exists in CRISPR but not in expression,
    it will still be exported but a report will be generated.

    Args:
        output_path: Path where the CSV file should be saved
        metadata_columns: Optional list of metadata columns from models table to include.
                        If None, uses default columns matching gene effects export.
        models_filter: Optional SQL WHERE clause to filter models (e.g., "oncotree_lineage = 'Lung'")
        report_missing: If True, generates a detailed report of missing data

    Returns:
        Dictionary with export statistics and missing data report

    Example:
        >>> from pathlib import Path
        >>> from depmap_db.export import export_matched_expression_csv
        >>>
        >>> # Export expression data matched to CRISPR
        >>> result = export_matched_expression_csv(Path("matched_expression.csv"))
        >>> print(f"Exported {result['models_exported']} models")
        >>> print(f"Missing models: {result['models_in_crispr_not_expression']}")
    """
    db_manager = get_db_manager()

    # Default metadata columns matching gene effects
    if metadata_columns is None:
        metadata_columns = [
            "model_id",
            "cell_line_name",
            "stripped_cell_line_name",
            "oncotree_lineage",
            "oncotree_primary_disease",
            "oncotree_subtype",
            "oncotree_code",
            "ccle_name",
            "depmap_model_type",
        ]

    logger.info("Analyzing CRISPR and expression data for matching...")

    # Get models and genes from CRISPR data
    crispr_models_query = "SELECT DISTINCT model_id FROM gene_effects_wide"
    if models_filter:
        crispr_models_query = f"""
        SELECT DISTINCT gef.model_id
        FROM gene_effects_wide gef
        INNER JOIN models m ON gef.model_id = m.model_id
        WHERE {models_filter}
        """

    crispr_models_result = db_manager.execute(crispr_models_query)
    crispr_models = {row[0] for row in crispr_models_result.fetchall()}

    logger.info("Found %s models in CRISPR data", len(crispr_models))

    # Get genes from CRISPR
    crispr_genes_query = """
    SELECT column_name
    FROM information_schema.columns
    WHERE table_name = 'gene_effects_wide'
      AND column_name NOT IN ('model_id', 'created_at')
    ORDER BY column_name
    """
    crispr_genes_result = db_manager.execute(crispr_genes_query)
    crispr_genes = {row[0] for row in crispr_genes_result.fetchall()}

    logger.info("Found %s genes in CRISPR data", len(crispr_genes))

    # Get models and genes from expression data
    expression_models_query = (
        "SELECT DISTINCT model_id FROM gene_expression_wide"
    )
    if models_filter:
        expression_models_query = f"""
        SELECT DISTINCT gex.model_id
        FROM gene_expression_wide gex
        INNER JOIN models m ON gex.model_id = m.model_id
        WHERE {models_filter}
        """

    expression_models_result = db_manager.execute(expression_models_query)
    expression_models = {row[0] for row in expression_models_result.fetchall()}

    logger.info("Found %s models in expression data", len(expression_models))

    # Get genes from expression
    expression_genes_query = """
    SELECT column_name
    FROM information_schema.columns
    WHERE table_name = 'gene_expression_wide'
      AND column_name NOT IN ('model_id', 'sequencing_id', 'model_condition_id',
                               'is_default_entry', 'is_default_for_mc', 'created_at')
    ORDER BY column_name
    """
    expression_genes_result = db_manager.execute(expression_genes_query)
    expression_genes = {row[0] for row in expression_genes_result.fetchall()}

    logger.info("Found %s genes in expression data", len(expression_genes))

    # Calculate intersections and differences
    models_intersection = crispr_models & expression_models
    genes_intersection = crispr_genes & expression_genes

    models_in_crispr_not_expression = crispr_models - expression_models
    genes_in_crispr_not_expression = crispr_genes - expression_genes
    models_in_expression_not_crispr = expression_models - crispr_models
    genes_in_expression_not_crispr = expression_genes - crispr_genes

    # Create report dictionary
    report: dict[str, Any] = {
        "total_crispr_models": len(crispr_models),
        "total_crispr_genes": len(crispr_genes),
        "total_expression_models": len(expression_models),
        "total_expression_genes": len(expression_genes),
        "models_in_both": len(models_intersection),
        "genes_in_both": len(genes_intersection),
        "models_in_crispr_not_expression": len(
            models_in_crispr_not_expression
        ),
        "genes_in_crispr_not_expression": len(genes_in_crispr_not_expression),
        "models_in_expression_not_crispr": len(
            models_in_expression_not_crispr
        ),
        "genes_in_expression_not_crispr": len(genes_in_expression_not_crispr),
    }

    if report_missing:
        logger.info("%s", "\n" + "=" * 60)
        logger.info("DATA MATCHING REPORT")
        logger.info("%s", "=" * 60)
        logger.info(
            "CRISPR data: %s models, %s genes",
            report["total_crispr_models"],
            report["total_crispr_genes"],
        )
        logger.info(
            "Expression data: %s models, %s genes",
            report["total_expression_models"],
            report["total_expression_genes"],
        )
        logger.info("Models in both datasets: %s", report["models_in_both"])
        logger.info("Genes in both datasets: %s", report["genes_in_both"])

        if models_in_crispr_not_expression:
            logger.warning(
                "%s models in CRISPR but NOT in expression data",
                len(models_in_crispr_not_expression),
            )
            logger.warning(
                "   Missing models: %s",
                ", ".join(sorted(models_in_crispr_not_expression)[:20]),
            )

        if genes_in_crispr_not_expression:
            logger.warning(
                "%s genes in CRISPR but NOT in expression data",
                len(genes_in_crispr_not_expression),
            )
            logger.warning(
                "   Missing genes: %s",
                ", ".join(sorted(genes_in_crispr_not_expression)[:20]),
            )

        if models_in_expression_not_crispr:
            logger.info(
                "%s models in expression but not in CRISPR (will be excluded from export)",
                len(models_in_expression_not_crispr),
            )

        if genes_in_expression_not_crispr:
            logger.info(
                "%s genes in expression but not in CRISPR (will be excluded from export)",
                len(genes_in_expression_not_crispr),
            )

        logger.info("%s", "=" * 60)

    # Check if we have data to export
    if not models_intersection:
        logger.error(
            "No overlapping models found between CRISPR and expression data!"
        )
        report["models_exported"] = 0
        return report

    if not genes_intersection:
        logger.error(
            "No overlapping genes found between CRISPR and expression data!"
        )
        report["models_exported"] = 0
        return report

    # Build export query for matching models and genes only
    logger.info(
        "Exporting %s models with %s genes...",
        len(models_intersection),
        len(genes_intersection),
    )

    # Build SELECT clause for metadata
    metadata_select = ", ".join([f"m.{col}" for col in metadata_columns])

    # Build SELECT clause for gene expression values (only genes in CRISPR)
    gene_select = ", ".join(
        [f'gex."{gene}"' for gene in sorted(genes_intersection)]
    )

    # Build WHERE clause for models in CRISPR
    models_list = "', '".join(models_intersection)
    crispr_models_filter = f"gex.model_id IN ('{models_list}')"

    if models_filter:
        where_clause = f"WHERE {crispr_models_filter} AND ({models_filter})"
    else:
        where_clause = f"WHERE {crispr_models_filter}"

    # Full query
    export_query = f"""
    SELECT
        {metadata_select},
        {gene_select}
    FROM models m
    INNER JOIN gene_expression_wide gex ON m.model_id = gex.model_id
    {where_clause}
    ORDER BY m.cell_line_name
    """

    logger.info("Executing export query using optimized COPY command...")

    # Use DuckDB's COPY TO command for memory-efficient export
    output_path.parent.mkdir(parents=True, exist_ok=True)

    copy_query = f"""
    COPY (
        {export_query}
    ) TO '{output_path}' (HEADER, DELIMITER ',')
    """

    db_manager.execute(copy_query)

    # Count the exported records
    count_query = f"""
    SELECT COUNT(*) as count FROM (
        {export_query}
    )
    """
    result = db_manager.execute(count_query)
    num_records: int = result.fetchone()[0]

    file_size = output_path.stat().st_size
    logger.info(
        "Successfully exported %s models to %s (%s bytes)",
        num_records,
        output_path,
        file_size,
    )

    report["models_exported"] = num_records
    report["genes_exported"] = len(genes_intersection)
    report["output_file"] = str(output_path)
    report["file_size_bytes"] = file_size

    return report


def export_joint_crispr_expression_csv(
    expression_output_path: Path,
    crispr_output_path: Path,
    metadata_columns: list[str] | None = None,
    models_filter: str | None = None,
    report_missing: bool = True,
) -> dict[str, Any]:
    """Export perfectly matched CRISPR and expression data for cross-correlation analysis.

    This function ensures that both exports contain exactly the same models and genes,
    making them suitable for direct correlation analysis. Only the intersection of
    models and genes present in both datasets will be exported.

    Args:
        expression_output_path: Path for expression CSV export
        crispr_output_path: Path for CRISPR gene effects CSV export
        metadata_columns: Optional list of metadata columns from models table to include.
                        If None, uses default columns.
        models_filter: Optional SQL WHERE clause to filter models (e.g., "oncotree_lineage = 'Lung'")
        report_missing: If True, generates a detailed report of missing data

    Returns:
        Dictionary with export statistics and missing data report

    Example:
        >>> from pathlib import Path
        >>> from depmap_db.export import export_joint_crispr_expression_csv
        >>>
        >>> # Export perfectly matched datasets
        >>> result = export_joint_crispr_expression_csv(
        ...     expression_output_path=Path("expression_matched.csv"),
        ...     crispr_output_path=Path("crispr_matched.csv")
        ... )
        >>> print(f"Exported {result['models_exported']} models x {result['genes_exported']} genes")
    """
    db_manager = get_db_manager()

    # Default metadata columns
    if metadata_columns is None:
        metadata_columns = [
            "model_id",
            "cell_line_name",
            "stripped_cell_line_name",
            "oncotree_lineage",
            "oncotree_primary_disease",
            "oncotree_subtype",
            "oncotree_code",
            "ccle_name",
            "depmap_model_type",
        ]

    logger.info("Analyzing CRISPR and expression data for perfect matching...")

    # Get models and genes from CRISPR data
    crispr_models_query = "SELECT DISTINCT model_id FROM gene_effects_wide"
    if models_filter:
        crispr_models_query = f"""
        SELECT DISTINCT gef.model_id
        FROM gene_effects_wide gef
        INNER JOIN models m ON gef.model_id = m.model_id
        WHERE {models_filter}
        """

    crispr_models_result = db_manager.execute(crispr_models_query)
    crispr_models = {row[0] for row in crispr_models_result.fetchall()}

    # Get genes from CRISPR
    crispr_genes_query = """
    SELECT column_name
    FROM information_schema.columns
    WHERE table_name = 'gene_effects_wide'
      AND column_name NOT IN ('model_id', 'created_at')
    ORDER BY column_name
    """
    crispr_genes_result = db_manager.execute(crispr_genes_query)
    crispr_genes = {row[0] for row in crispr_genes_result.fetchall()}

    # Get models and genes from expression data
    expression_models_query = (
        "SELECT DISTINCT model_id FROM gene_expression_wide"
    )
    if models_filter:
        expression_models_query = f"""
        SELECT DISTINCT gex.model_id
        FROM gene_expression_wide gex
        INNER JOIN models m ON gex.model_id = m.model_id
        WHERE {models_filter}
        """

    expression_models_result = db_manager.execute(expression_models_query)
    expression_models = {row[0] for row in expression_models_result.fetchall()}

    # Get genes from expression
    expression_genes_query = """
    SELECT column_name
    FROM information_schema.columns
    WHERE table_name = 'gene_expression_wide'
      AND column_name NOT IN ('model_id', 'sequencing_id', 'model_condition_id',
                               'is_default_entry', 'is_default_for_mc', 'created_at')
    ORDER BY column_name
    """
    expression_genes_result = db_manager.execute(expression_genes_query)
    expression_genes = {row[0] for row in expression_genes_result.fetchall()}

    # Calculate intersections
    models_intersection = crispr_models & expression_models
    genes_intersection = crispr_genes & expression_genes

    models_in_crispr_not_expression = crispr_models - expression_models
    genes_in_crispr_not_expression = crispr_genes - expression_genes
    models_in_expression_not_crispr = expression_models - crispr_models
    genes_in_expression_not_crispr = expression_genes - crispr_genes

    # Create report dictionary
    report: dict[str, Any] = {
        "total_crispr_models": len(crispr_models),
        "total_crispr_genes": len(crispr_genes),
        "total_expression_models": len(expression_models),
        "total_expression_genes": len(expression_genes),
        "models_in_both": len(models_intersection),
        "genes_in_both": len(genes_intersection),
        "models_in_crispr_not_expression": len(
            models_in_crispr_not_expression
        ),
        "genes_in_crispr_not_expression": len(genes_in_crispr_not_expression),
        "models_in_expression_not_crispr": len(
            models_in_expression_not_crispr
        ),
        "genes_in_expression_not_crispr": len(genes_in_expression_not_crispr),
    }

    if report_missing:
        logger.info("%s", "\n" + "=" * 60)
        logger.info("JOINT EXPORT - PERFECT MATCHING REPORT")
        logger.info("%s", "=" * 60)
        logger.info(
            "CRISPR data: %s models, %s genes",
            report["total_crispr_models"],
            report["total_crispr_genes"],
        )
        logger.info(
            "Expression data: %s models, %s genes",
            report["total_expression_models"],
            report["total_expression_genes"],
        )
        logger.info(
            "Perfect intersection: %s models x %s genes",
            report["models_in_both"],
            report["genes_in_both"],
        )
        logger.info(
            "  Both exports will contain exactly the same models and genes"
        )

        if models_in_crispr_not_expression:
            logger.warning(
                "Excluding %s CRISPR-only models",
                len(models_in_crispr_not_expression),
            )
            if len(models_in_crispr_not_expression) <= 10:
                logger.warning(
                    "   Models: %s",
                    ", ".join(sorted(models_in_crispr_not_expression)),
                )

        if genes_in_crispr_not_expression:
            logger.warning(
                "Excluding %s CRISPR-only genes",
                len(genes_in_crispr_not_expression),
            )
            if len(genes_in_crispr_not_expression) <= 10:
                logger.warning(
                    "   Genes: %s",
                    ", ".join(sorted(genes_in_crispr_not_expression)),
                )

        if models_in_expression_not_crispr:
            logger.warning(
                "Excluding %s expression-only models",
                len(models_in_expression_not_crispr),
            )

        if genes_in_expression_not_crispr:
            logger.warning(
                "Excluding %s expression-only genes",
                len(genes_in_expression_not_crispr),
            )

        logger.info("%s", "=" * 60)

    # Check if we have data to export
    if not models_intersection:
        logger.error(
            "No overlapping models found between CRISPR and expression data!"
        )
        report["models_exported"] = 0
        report["genes_exported"] = 0
        return report

    if not genes_intersection:
        logger.error(
            "No overlapping genes found between CRISPR and expression data!"
        )
        report["models_exported"] = 0
        report["genes_exported"] = 0
        return report

    logger.info(
        "Exporting perfectly matched datasets: %s models x %s genes...",
        len(models_intersection),
        len(genes_intersection),
    )

    # Build SELECT clause for metadata
    metadata_select = ", ".join([f"m.{col}" for col in metadata_columns])

    # Build SELECT clause for genes (only genes in intersection)
    genes_sorted = sorted(genes_intersection)

    # Build WHERE clause for models in intersection
    models_list = "', '".join(models_intersection)
    matched_models_filter = f"m.model_id IN ('{models_list}')"

    if models_filter:
        where_clause = f"WHERE {matched_models_filter} AND ({models_filter})"
    else:
        where_clause = f"WHERE {matched_models_filter}"

    # === Export Expression Data ===
    logger.info("Exporting expression data...")

    expression_gene_select = ", ".join(
        [f'gex."{gene}"' for gene in genes_sorted]
    )

    expression_query = f"""
    SELECT
        {metadata_select},
        {expression_gene_select}
    FROM models m
    INNER JOIN gene_expression_wide gex ON m.model_id = gex.model_id
    {where_clause}
    ORDER BY m.cell_line_name
    """

    expression_output_path.parent.mkdir(parents=True, exist_ok=True)

    expression_copy_query = f"""
    COPY (
        {expression_query}
    ) TO '{expression_output_path}' (HEADER, DELIMITER ',')
    """

    db_manager.execute(expression_copy_query)
    expression_file_size = expression_output_path.stat().st_size

    # === Export CRISPR Data ===
    logger.info("Exporting CRISPR gene effects data...")

    crispr_gene_select = ", ".join([f'gef."{gene}"' for gene in genes_sorted])

    crispr_query = f"""
    SELECT
        {metadata_select},
        {crispr_gene_select}
    FROM models m
    INNER JOIN gene_effects_wide gef ON m.model_id = gef.model_id
    {where_clause}
    ORDER BY m.cell_line_name
    """

    crispr_output_path.parent.mkdir(parents=True, exist_ok=True)

    crispr_copy_query = f"""
    COPY (
        {crispr_query}
    ) TO '{crispr_output_path}' (HEADER, DELIMITER ',')
    """

    db_manager.execute(crispr_copy_query)
    crispr_file_size = crispr_output_path.stat().st_size

    # Count the exported records (should be same for both)
    count_query = f"""
    SELECT COUNT(*) as count FROM (
        {expression_query}
    )
    """
    result = db_manager.execute(count_query)
    num_records: int = result.fetchone()[0]

    logger.info("Successfully exported perfectly matched datasets:")
    logger.info(
        "  Expression: %s models x %s genes -> %s (%s bytes)",
        num_records,
        len(genes_intersection),
        expression_output_path,
        expression_file_size,
    )
    logger.info(
        "  CRISPR:     %s models x %s genes -> %s (%s bytes)",
        num_records,
        len(genes_intersection),
        crispr_output_path,
        crispr_file_size,
    )

    report["models_exported"] = num_records
    report["genes_exported"] = len(genes_intersection)
    report["expression_output_file"] = str(expression_output_path)
    report["crispr_output_file"] = str(crispr_output_path)
    report["expression_file_size_bytes"] = expression_file_size
    report["crispr_file_size_bytes"] = crispr_file_size

    return report
