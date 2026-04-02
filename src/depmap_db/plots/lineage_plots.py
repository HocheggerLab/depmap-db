"""Lineage-level plotting for DepMap gene dependency and expression."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import polars as pl
import seaborn as sns

from depmap_db import (
    scan_gene_effects_wide,
    scan_gene_expression_wide,
    scan_models,
)

from .colors import COLOR

_STYLE = Path(__file__).parent / "hhlab_style01.mplstyle"


def analyse_by_lineage(
    gene: str,
    assay: Literal["dependency", "expression"] = "dependency",
    db_path: str | Path | None = None,
) -> pl.DataFrame:
    """Get CRISPR gene effect or expression across lineages."""
    models = scan_models(db_path=db_path).select(
        "model_id", "oncotree_lineage"
    )
    data = (
        scan_gene_effects_wide(db_path=db_path)
        if assay == "dependency"
        else scan_gene_expression_wide(db_path=db_path)
    )
    # match exact name or "MASTL (84930)" format
    names = data.collect_schema().names()
    if gene in names:
        col = gene
    else:
        col = next(c for c in names if c.startswith(f"{gene} ("))
    return models.join(
        data.select("model_id", pl.col(col)), on="model_id"
    ).rename({col: gene})


def plot_lineage(
    df: pl.DataFrame,
    gene: str,
    assay: Literal["dependency", "expression"] = "dependency",
    *,
    min_models: int = 10,
    figsize: tuple[float, float] | None = None,
    ax: plt.Axes | None = None,
) -> tuple[plt.Figure, plt.Axes]:
    """Box + strip plot of gene dependency or expression by oncotree lineage.

    Args:
        df: DataFrame with columns ``oncotree_lineage`` and ``gene``.
        gene: Gene symbol (must match a column in *df*).
        assay: Controls axis label and reference line.
        min_models: Exclude lineages with fewer than this many models.
        figsize: Override default figure size.
        ax: Optional pre-existing axes to draw on.

    Returns:
        Tuple of (Figure, Axes).
    """
    plt.style.use(str(_STYLE))

    # collect if lazy, filter nulls and small lineages
    if isinstance(df, pl.LazyFrame):
        df = df.collect()
    filtered = df.drop_nulls(gene).filter(
        pl.col("oncotree_lineage").count().over("oncotree_lineage")
        >= min_models
    )

    # sort lineages by median value
    order = (
        filtered.group_by("oncotree_lineage")
        .agg(pl.col(gene).median().alias("median"))
        .sort("median")
        .get_column("oncotree_lineage")
        .to_list()
    )

    # convert to pandas at the seaborn boundary
    pdf = filtered.to_pandas()

    if ax is None:
        fig, ax = plt.subplots(
            figsize=figsize or (max(4.0, len(order) * 0.25), 2.5)
        )
    else:
        fig = ax.get_figure()

    sns.boxplot(
        data=pdf,
        x="oncotree_lineage",
        y=gene,
        order=order,
        color=COLOR.GREY.value,
        fliersize=0,
        linewidth=0.5,
        width=0.6,
        ax=ax,
    )
    sns.stripplot(
        data=pdf,
        x="oncotree_lineage",
        y=gene,
        order=order,
        color=COLOR.BLUE.value,
        size=1.5,
        alpha=0.5,
        jitter=True,
        ax=ax,
    )

    # reference line
    if assay == "dependency":
        ax.axhline(-1, ls="--", lw=0.5, color=COLOR.DARKGREY.value, zorder=0)
        ax.set_ylabel(f"{gene} gene effect")
    else:
        ax.set_ylabel(f"{gene} log2(TPM+1)")

    ax.set_xlabel("")
    plt.setp(
        ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor"
    )
    ax.set_title(f"{gene} — {assay} by lineage", fontsize=7)

    fig.tight_layout()
    return fig, ax


def lineage_analysis(
    gene: str,
    assay: Literal["dependency", "expression"] = "dependency",
    db_path: str | Path | None = None,
    min_models: int = 3,
    figsize: tuple[float, float] | None = None,
) -> plt.Figure:
    """Convenience function to run lineage analysis and plotting in one step."""
    df = analyse_by_lineage(gene=gene, assay=assay, db_path=db_path)
    fig, _ = plot_lineage(
        df=df, gene=gene, assay=assay, min_models=min_models, figsize=figsize
    )
    return fig
