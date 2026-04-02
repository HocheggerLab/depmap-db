"""Focused tests for the plotting API."""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import pytest

os.environ.setdefault("LOG_LEVEL", "INFO")
os.environ.setdefault("LOG_FILE_PATH", "logs/app.log")
matplotlib.use("Agg")

from depmap_db.config import reload_settings
from depmap_db.database import create_tables, get_db_manager, reset_db_manager
from depmap_db.plots import (
    analyse_biomarker_volcano,
    analyse_comutation,
    analyse_drug_by_lineage,
    analyse_drug_by_mutation,
    analyse_drug_dependency,
    analyse_expr_vs_dep,
    analyse_gene_gene,
    analyse_gene_heatmap,
    analyse_mutation_waterfall,
    analyse_selectivity,
    analyse_top_correlations,
    biomarker_volcano,
    drug_by_lineage,
    drug_by_mutation,
    drug_dependency_scatter,
    expr_vs_dep,
    gene_gene_scatter,
    gene_heatmap,
    mutation_waterfall,
    plot_biomarker_volcano,
    plot_comutation,
    plot_drug_by_lineage,
    plot_drug_by_mutation,
    plot_drug_dependency,
    plot_expr_vs_dep,
    plot_gene_gene,
    plot_gene_heatmap,
    plot_mutation_waterfall,
    plot_selectivity,
    plot_top_correlations,
    selectivity_plot,
    top_correlations,
)


@pytest.fixture()
def plot_database(tmp_path: Path) -> Path:
    """Create a compact DuckDB fixture that covers plot inputs."""
    db_path = tmp_path / "depmap-plots.duckdb"

    settings = reload_settings()
    settings.database.path = db_path
    settings.database.memory = False
    reset_db_manager()
    create_tables()

    db = get_db_manager()
    for gene in ["TP53", "MDM2", "EGFR", "CDK4"]:
        db.execute(f'ALTER TABLE gene_effects_wide ADD COLUMN "{gene}" DOUBLE')
        db.execute(
            f'ALTER TABLE gene_expression_wide ADD COLUMN "{gene}" DOUBLE'
        )

    db.execute_many(
        """
        INSERT INTO models (
            model_id,
            cell_line_name,
            stripped_cell_line_name,
            ccle_name,
            oncotree_lineage,
            oncotree_primary_disease,
            oncotree_subtype
        ) VALUES (?, ?, ?, ?, ?, ?, ?)
        """,
        [
            (
                "ACH-000001",
                "CellA",
                "CELLA",
                "CELLA_CCLE",
                "Breast",
                "Breast Cancer",
                "SubtypeA",
            ),
            (
                "ACH-000002",
                "CellB",
                "CELLB",
                "CELLB_CCLE",
                "Breast",
                "Breast Cancer",
                "SubtypeA",
            ),
            (
                "ACH-000003",
                "CellC",
                "CELLC",
                "CELLC_CCLE",
                "Lung",
                "Lung Cancer",
                "SubtypeB",
            ),
            (
                "ACH-000004",
                "CellD",
                "CELLD",
                "CELLD_CCLE",
                "Lung",
                "Lung Cancer",
                "SubtypeB",
            ),
        ],
    )

    db.execute_many(
        'INSERT INTO gene_effects_wide (model_id, "TP53", "MDM2", "EGFR", "CDK4") VALUES (?, ?, ?, ?, ?)',
        [
            ("ACH-000001", -1.8, -1.6, -0.2, -1.4),
            ("ACH-000002", -1.2, -1.0, -0.4, -1.0),
            ("ACH-000003", -0.5, -0.4, -1.4, -0.6),
            ("ACH-000004", -0.2, -0.1, -1.8, -0.3),
        ],
    )
    db.execute_many(
        'INSERT INTO gene_expression_wide (model_id, "TP53", "MDM2", "EGFR", "CDK4") VALUES (?, ?, ?, ?, ?)',
        [
            ("ACH-000001", 9.5, 9.0, 2.0, 8.0),
            ("ACH-000002", 8.7, 8.3, 2.5, 7.6),
            ("ACH-000003", 4.1, 4.2, 8.8, 4.5),
            ("ACH-000004", 3.0, 3.2, 9.4, 3.8),
        ],
    )

    db.execute_many(
        """
        INSERT INTO mutations (
            mutation_id,
            model_id,
            chrom,
            pos,
            ref,
            alt,
            variant_type,
            protein_change,
            hugo_symbol,
            molecular_consequence,
            vep_impact,
            likely_lof,
            hotspot,
            hess_driver,
            tumor_suppressor_high_impact,
            oncogene_high_impact
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        [
            (
                "mut-1",
                "ACH-000001",
                "17",
                7579472,
                "G",
                "A",
                "SNP",
                "p.R248Q",
                "TP53",
                "missense_variant",
                "MODERATE",
                True,
                True,
                True,
                True,
                False,
            ),
            (
                "mut-2",
                "ACH-000002",
                "17",
                7574000,
                "C",
                "T",
                "DEL",
                "p.E285*",
                "TP53",
                "stop_gained",
                "HIGH",
                True,
                False,
                False,
                True,
                False,
            ),
            (
                "mut-3",
                "ACH-000003",
                "7",
                55249071,
                "G",
                "A",
                "SNP",
                "p.L858R",
                "EGFR",
                "missense_variant",
                "MODERATE",
                False,
                True,
                True,
                False,
                True,
            ),
            (
                "mut-4",
                "ACH-000004",
                "7",
                55242465,
                "C",
                "T",
                "SNP",
                "p.T790M",
                "EGFR",
                "missense_variant",
                "MODERATE",
                False,
                True,
                True,
                False,
                True,
            ),
        ],
    )

    db.execute(
        """
        INSERT INTO compounds (
            broad_id,
            compound_name,
            moa,
            target_text,
            source_dataset,
            source_filename,
            release_label
        ) VALUES (
            'BRD-K11111111-001-01-1',
            'palbociclib',
            'CDK inhibitor',
            'CDK4',
            'PRISMSecondaryDoseResponseCurveParameters',
            'secondary.csv',
            'Repurposing Public 24Q2'
        )
        """
    )
    db.execute(
        """
        INSERT INTO drug_screens (
            screen_id,
            screen_kind,
            dataset_name,
            source_filename,
            release_label,
            release_track,
            default_secondary_summary_metric
        ) VALUES (
            'MTS010',
            'secondary',
            'PRISMSecondaryDoseResponseCurveParameters',
            'secondary.csv',
            'Repurposing Public 24Q2',
            'prism_secondary',
            'auc'
        )
        """
    )

    db.execute_many(
        """
        INSERT INTO drug_response_secondary (
            response_id,
            broad_id,
            model_id,
            ccle_name,
            screen_id,
            auc,
            ec50,
            ic50,
            compound_name,
            moa,
            target_text,
            source_dataset,
            source_filename
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        [
            (
                "resp-1",
                "BRD-K11111111-001-01-1",
                "ACH-000001",
                "CELLA_CCLE",
                "MTS010",
                0.25,
                0.8,
                1.0,
                "palbociclib",
                "CDK inhibitor",
                "CDK4",
                "PRISMSecondaryDoseResponseCurveParameters",
                "secondary.csv",
            ),
            (
                "resp-2",
                "BRD-K11111111-001-01-1",
                "ACH-000002",
                "CELLB_CCLE",
                "MTS010",
                0.35,
                1.0,
                1.2,
                "palbociclib",
                "CDK inhibitor",
                "CDK4",
                "PRISMSecondaryDoseResponseCurveParameters",
                "secondary.csv",
            ),
            (
                "resp-3",
                "BRD-K11111111-001-01-1",
                "ACH-000003",
                "CELLC_CCLE",
                "MTS010",
                0.72,
                3.0,
                4.0,
                "palbociclib",
                "CDK inhibitor",
                "CDK4",
                "PRISMSecondaryDoseResponseCurveParameters",
                "secondary.csv",
            ),
            (
                "resp-4",
                "BRD-K11111111-001-01-1",
                "ACH-000004",
                "CELLD_CCLE",
                "MTS010",
                0.88,
                4.2,
                5.1,
                "palbociclib",
                "CDK inhibitor",
                "CDK4",
                "PRISMSecondaryDoseResponseCurveParameters",
                "secondary.csv",
            ),
        ],
    )

    reset_db_manager()
    return db_path


def test_plot_analyses_cover_all_requested_families(plot_database: Path) -> None:
    """Analysis helpers should produce tidy frames for each plot family."""
    gene_gene = analyse_gene_gene(
        "TP53",
        "MDM2",
        color_by="mutation",
        mutation_gene="TP53",
        db_path=plot_database,
    )
    assert gene_gene.height == 4
    assert set(gene_gene.columns) == {"model_id", "TP53", "MDM2", "mutation_status"}

    top = analyse_top_correlations("TP53", db_path=plot_database, n=3)
    assert top.height == 3
    assert "gene" in top.columns

    expr_dep = analyse_expr_vs_dep(
        "TP53",
        color_by="mutation",
        mutation_gene="TP53",
        db_path=plot_database,
    )
    assert set(expr_dep.columns) == {"model_id", "expression", "dependency", "mutation_status"}

    volcano = analyse_biomarker_volcano("TP53", db_path=plot_database)
    assert volcano.height >= 3
    assert {"gene", "effect_size", "p", "neg_log10_p"}.issubset(volcano.columns)

    drug_lineage = analyse_drug_by_lineage(
        "palbociclib", db_path=plot_database
    )
    assert set(drug_lineage.get_column("oncotree_lineage")) == {"Breast", "Lung"}

    drug_mut = analyse_drug_by_mutation(
        "palbociclib", "TP53", db_path=plot_database
    )
    assert set(drug_mut.get_column("mutation_status")) == {"WT", "mutant"}

    drug_dep = analyse_drug_dependency(
        "palbociclib", "CDK4", db_path=plot_database
    )
    assert set(drug_dep.columns) == {"model_id", "auc", "dependency"}

    waterfall = analyse_mutation_waterfall(
        "TP53", "TP53", db_path=plot_database
    )
    assert {"model_id", "value", "rank", "mutation_status", "variant_type"}.issubset(
        waterfall.columns
    )

    comutation = analyse_comutation(
        "TP53",
        "TP53",
        "EGFR",
        gene_type_b="oncogene",
        db_path=plot_database,
    )
    assert set(comutation.get_column("group")) == {
        "TP53 only",
        "EGFR only",
    }

    selectivity = analyse_selectivity("TP53", db_path=plot_database)
    assert selectivity.height == 4
    assert "rank" in selectivity.columns

    heatmap = analyse_gene_heatmap(
        ["TP53", "MDM2", "EGFR"],
        assay="expression",
        db_path=plot_database,
    )
    assert heatmap.columns == ["model_id", "TP53", "MDM2", "EGFR"]


def test_plot_functions_return_figures_or_clustergrid(plot_database: Path) -> None:
    """Plotting helpers should render without interactive backends."""
    gene_gene_df = analyse_gene_gene(
        "TP53",
        "MDM2",
        color_by="mutation",
        mutation_gene="TP53",
        db_path=plot_database,
    )
    top_df = analyse_top_correlations("TP53", db_path=plot_database, n=3)
    expr_dep_df = analyse_expr_vs_dep(
        "TP53",
        color_by="mutation",
        mutation_gene="TP53",
        db_path=plot_database,
    )
    volcano_df = analyse_biomarker_volcano("TP53", db_path=plot_database)
    drug_lineage_df = analyse_drug_by_lineage("palbociclib", db_path=plot_database)
    drug_mut_df = analyse_drug_by_mutation(
        "palbociclib", "TP53", db_path=plot_database
    )
    drug_dep_df = analyse_drug_dependency(
        "palbociclib", "CDK4", db_path=plot_database
    )
    waterfall_df = analyse_mutation_waterfall(
        "TP53", "TP53", db_path=plot_database
    )
    comutation_df = analyse_comutation(
        "TP53",
        "TP53",
        "EGFR",
        gene_type_b="oncogene",
        db_path=plot_database,
    )
    selectivity_df = analyse_selectivity("TP53", db_path=plot_database)
    heatmap_df = analyse_gene_heatmap(
        ["TP53", "MDM2", "EGFR"],
        assay="expression",
        db_path=plot_database,
    )

    figures = [
        plot_gene_gene(gene_gene_df, "TP53", "MDM2", color_by="mutation")[0],
        gene_gene_scatter(
            "TP53",
            "MDM2",
            color_by="mutation",
            mutation_gene="TP53",
            db_path=plot_database,
        ),
        plot_top_correlations(top_df, "TP53")[0],
        top_correlations("TP53", db_path=plot_database, n=3),
        plot_expr_vs_dep(expr_dep_df, "TP53", color_by="mutation")[0],
        expr_vs_dep(
            "TP53",
            color_by="mutation",
            mutation_gene="TP53",
            db_path=plot_database,
        ),
        plot_biomarker_volcano(volcano_df, "TP53")[0],
        biomarker_volcano("TP53", db_path=plot_database),
        plot_drug_by_lineage(
            drug_lineage_df,
            "palbociclib",
            min_models=1,
        )[0],
        drug_by_lineage("palbociclib", min_models=1, db_path=plot_database),
        plot_drug_by_mutation(drug_mut_df, "palbociclib", "TP53")[0],
        drug_by_mutation("palbociclib", "TP53", db_path=plot_database),
        plot_drug_dependency(drug_dep_df, "palbociclib", "CDK4")[0],
        drug_dependency_scatter(
            "palbociclib", "CDK4", db_path=plot_database
        ),
        plot_mutation_waterfall(waterfall_df, "TP53", "TP53")[0],
        mutation_waterfall("TP53", "TP53", db_path=plot_database),
        plot_comutation(comutation_df, "TP53", "TP53", "EGFR")[0],
        plot_selectivity(selectivity_df, "TP53")[0],
        selectivity_plot("TP53", db_path=plot_database),
    ]

    for figure in figures:
        assert hasattr(figure, "axes")
        assert figure.axes
        plt.close(figure)

    cluster_grid = plot_gene_heatmap(
        heatmap_df,
        ["TP53", "MDM2", "EGFR"],
        assay="expression",
        cluster_rows=False,
        cluster_cols=False,
    )
    assert hasattr(cluster_grid, "fig")
    assert cluster_grid.fig.axes
    plt.close(cluster_grid.fig)

    heatmap_figure = gene_heatmap(
        ["TP53", "MDM2", "EGFR"],
        assay="expression",
        cluster_rows=False,
        cluster_cols=False,
        db_path=plot_database,
    )
    assert heatmap_figure.axes
    plt.close(heatmap_figure)
