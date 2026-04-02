"""Mutation waterfall and co-mutation enrichment plots.

Covers requested plot families:
8) Mutation waterfall / oncoplot strip: dependency ranked with mutation annotation.
9) Co-mutation enrichment: grouped box plot across 4 mutation groups.
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import polars as pl
import seaborn as sns

from depmap_db import scan_mutations

from ._helpers import (
    extract_gene_values,
    get_mutant_model_ids,
    label_mutation_status,
    make_axes,
    use_style,
)
from .colors import COLOR

# ---------------------------------------------------------------------------
# 8. Mutation waterfall
# ---------------------------------------------------------------------------


def analyse_mutation_waterfall(
    target_gene: str,
    mutation_gene: str,
    assay: Literal["dependency", "expression"] = "dependency",
    *,
    gene_type: Literal["suppressor", "oncogene"] = "suppressor",
    db_path: str | Path | None = None,
) -> pl.DataFrame:
    """Build a ranked waterfall DataFrame annotated with mutation status.

    Args:
        target_gene: Gene whose assay values are plotted.
        mutation_gene: Gene whose mutation status annotates the bars.
        assay: Assay to pull values from.
        gene_type: Mutation filter mode.
        db_path: Optional DuckDB path override.

    Returns:
        DataFrame with ``model_id``, ``value``, ``rank``,
        ``mutation_status``, and optionally ``variant_type`` (if available).
    """
    vals = (
        extract_gene_values(assay, target_gene, db_path=db_path)
        .rename({target_gene: "value"})
        .collect()
        .drop_nulls("value")
    )
    mutant_ids = get_mutant_model_ids(
        mutation_gene, gene_type, db_path=db_path
    )
    df = label_mutation_status(vals, mutant_ids)

    # attempt to attach variant_type for mutant models
    try:
        vt = (
            scan_mutations(db_path=db_path)
            .filter(pl.col("hugo_symbol") == mutation_gene)
            .select("model_id", "variant_type")
            .unique(subset=["model_id"])
            .collect()
        )
        df = df.join(vt, on="model_id", how="left")
    except (pl.exceptions.ColumnNotFoundError, pl.exceptions.SchemaError):
        df = df.with_columns(pl.lit(None).alias("variant_type"))

    return (
        df.sort("value")
        .with_row_index("rank")
        .with_columns(pl.col("rank").cast(pl.Int64))
    )


def plot_mutation_waterfall(
    df: pl.DataFrame,
    target_gene: str,
    mutation_gene: str,
    assay: Literal["dependency", "expression"] = "dependency",
    *,
    figsize: tuple[float, float] | None = None,
    ax: plt.Axes | None = None,
) -> tuple[plt.Figure, plt.Axes]:
    """Waterfall bar chart of ranked assay values with mutation annotation.

    Mutant cell lines are coloured distinctly; WT is grey.

    Args:
        df: DataFrame from :func:`analyse_mutation_waterfall`.
        target_gene: For axis label.
        mutation_gene: For title.
        assay: Controls y-axis label.
        figsize: Override figure size.
        ax: Optional pre-existing axes.

    Returns:
        Tuple of (Figure, Axes).
    """
    use_style()
    n = df.height
    fig, ax = make_axes(figsize, ax, default_figsize=(4.0, 2.0))
    pdf = df.to_pandas()

    colors = [
        COLOR.PINK.value if s == "mutant" else COLOR.GREY.value
        for s in pdf["mutation_status"]
    ]
    ax.bar(
        pdf["rank"], pdf["value"],
        width=1.0, color=colors, edgecolor="none",
    )

    if assay == "dependency":
        ax.axhline(-1, ls="--", lw=0.5, color=COLOR.DARKGREY.value, zorder=0)
        ax.set_ylabel(f"{target_gene} gene effect")
    else:
        ax.set_ylabel(f"{target_gene} log2(TPM+1)")

    n_mut = int((pdf["mutation_status"] == "mutant").sum())
    ax.set_xlabel(f"Cell lines (n={n})")
    ax.set_title(
        f"{target_gene} — ranked by {assay}, {mutation_gene} mutants "
        f"highlighted (n={n_mut})",
        fontsize=6,
    )
    ax.set_xticks([])

    fig.tight_layout()
    return fig, ax


def mutation_waterfall(
    target_gene: str,
    mutation_gene: str,
    assay: Literal["dependency", "expression"] = "dependency",
    *,
    gene_type: Literal["suppressor", "oncogene"] = "suppressor",
    db_path: str | Path | None = None,
    figsize: tuple[float, float] | None = None,
) -> plt.Figure:
    """Convenience: analyse and plot mutation waterfall in one step."""
    df = analyse_mutation_waterfall(
        target_gene, mutation_gene, assay,
        gene_type=gene_type, db_path=db_path,
    )
    fig, _ = plot_mutation_waterfall(
        df, target_gene, mutation_gene, assay, figsize=figsize,
    )
    return fig


# ---------------------------------------------------------------------------
# 9. Co-mutation enrichment
# ---------------------------------------------------------------------------


def analyse_comutation(
    target_gene: str,
    gene_a: str,
    gene_b: str,
    assay: Literal["dependency", "expression"] = "dependency",
    *,
    gene_type_a: Literal["suppressor", "oncogene"] = "suppressor",
    gene_type_b: Literal["suppressor", "oncogene"] = "suppressor",
    db_path: str | Path | None = None,
) -> pl.DataFrame:
    """Build DataFrame splitting cell lines into 4 co-mutation groups.

    Groups: ``both``, ``A_only`` (*gene_a* mutant only),
    ``B_only`` (*gene_b* mutant only), ``neither``.

    Args:
        target_gene: Gene whose assay values are measured.
        gene_a: First mutation gene.
        gene_b: Second mutation gene.
        assay: Which assay.
        gene_type_a: Mutation filter for gene_a.
        gene_type_b: Mutation filter for gene_b.
        db_path: Optional DuckDB path override.

    Returns:
        DataFrame with ``model_id``, ``value``, ``group``.
    """
    vals = (
        extract_gene_values(assay, target_gene, db_path=db_path)
        .rename({target_gene: "value"})
        .collect()
        .drop_nulls("value")
    )
    mut_a = get_mutant_model_ids(gene_a, gene_type_a, db_path=db_path)
    mut_b = get_mutant_model_ids(gene_b, gene_type_b, db_path=db_path)

    list_a = mut_a.to_list()
    list_b = mut_b.to_list()
    return vals.with_columns(
        pl.when(
            pl.col("model_id").is_in(list_a)
            & pl.col("model_id").is_in(list_b)
        )
        .then(pl.lit("both"))
        .when(pl.col("model_id").is_in(list_a))
        .then(pl.lit(f"{gene_a} only"))
        .when(pl.col("model_id").is_in(list_b))
        .then(pl.lit(f"{gene_b} only"))
        .otherwise(pl.lit("neither"))
        .alias("group")
    )


def plot_comutation(
    df: pl.DataFrame,
    target_gene: str,
    gene_a: str,
    gene_b: str,
    assay: Literal["dependency", "expression"] = "dependency",
    *,
    figsize: tuple[float, float] | None = None,
    ax: plt.Axes | None = None,
) -> tuple[plt.Figure, plt.Axes]:
    """Grouped box plot of co-mutation enrichment.

    Args:
        df: DataFrame from :func:`analyse_comutation`.
        target_gene: For axis label.
        gene_a: First mutation gene.
        gene_b: Second mutation gene.
        assay: Controls y-axis label.
        figsize: Override figure size.
        ax: Optional pre-existing axes.

    Returns:
        Tuple of (Figure, Axes).
    """
    use_style()
    fig, ax = make_axes(figsize, ax, default_figsize=(3.0, 2.5))
    pdf = df.to_pandas()

    order = ["neither", f"{gene_a} only", f"{gene_b} only", "both"]
    palette = {
        "neither": COLOR.GREY.value,
        f"{gene_a} only": COLOR.BLUE.value,
        f"{gene_b} only": COLOR.LIGHT_BLUE.value,
        "both": COLOR.PINK.value,
    }

    sns.boxplot(
        data=pdf,
        x="group",
        y="value",
        order=order,
        palette=palette,
        fliersize=0,
        linewidth=0.5,
        width=0.6,
        ax=ax,
    )
    sns.stripplot(
        data=pdf,
        x="group",
        y="value",
        order=order,
        color=COLOR.DARKGREY.value,
        size=1.5,
        alpha=0.4,
        jitter=True,
        ax=ax,
    )

    if assay == "dependency":
        ax.axhline(-1, ls="--", lw=0.5, color=COLOR.DARKGREY.value, zorder=0)
        ax.set_ylabel(f"{target_gene} gene effect")
    else:
        ax.set_ylabel(f"{target_gene} log2(TPM+1)")

    ax.set_xlabel("")
    plt.setp(
        ax.get_xticklabels(), rotation=30, ha="right", rotation_mode="anchor"
    )
    ax.set_title(
        f"{target_gene} — {gene_a}/{gene_b} co-mutation", fontsize=7
    )

    fig.tight_layout()
    return fig, ax


def comutation_analysis(
    target_gene: str,
    gene_a: str,
    gene_b: str,
    assay: Literal["dependency", "expression"] = "dependency",
    *,
    gene_type_a: Literal["suppressor", "oncogene"] = "suppressor",
    gene_type_b: Literal["suppressor", "oncogene"] = "suppressor",
    db_path: str | Path | None = None,
    figsize: tuple[float, float] | None = None,
) -> plt.Figure:
    """Convenience: analyse and plot co-mutation enrichment in one step."""
    df = analyse_comutation(
        target_gene, gene_a, gene_b, assay,
        gene_type_a=gene_type_a, gene_type_b=gene_type_b, db_path=db_path,
    )
    fig, _ = plot_comutation(
        df, target_gene, gene_a, gene_b, assay, figsize=figsize,
    )
    return fig
