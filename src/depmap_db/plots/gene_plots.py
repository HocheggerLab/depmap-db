"""Gene-gene scatter and top-correlated-genes plots.

Covers requested plot families:
1) Gene-gene scatter (dependency or expression) with correlation annotation.
2) Top correlated genes bar/lollipop chart ranked by correlation coefficient.
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import seaborn as sns
from scipy import stats

from ._helpers import (
    annotate_correlation,
    assay_ylabel,
    extract_gene_values,
    get_mutant_model_ids,
    join_models,
    label_mutation_status,
    make_axes,
    pearson_spearman,
    resolve_gene_column,
    scan_assay,
    use_style,
)
from .colors import COLOR

# ---------------------------------------------------------------------------
# 1. Gene–gene scatter
# ---------------------------------------------------------------------------


def analyse_gene_gene(
    gene_x: str,
    gene_y: str,
    assay: Literal["dependency", "expression"] = "dependency",
    *,
    color_by: Literal["lineage", "mutation", "none"] = "none",
    mutation_gene: str | None = None,
    gene_type: Literal["suppressor", "oncogene"] = "suppressor",
    db_path: str | Path | None = None,
) -> pl.DataFrame:
    """Build a DataFrame joining two gene values across cell lines.

    Args:
        gene_x: First gene (x-axis).
        gene_y: Second gene (y-axis).
        assay: Which assay to pull from.
        color_by: ``"lineage"`` attaches ``oncotree_lineage`` and
            ``"mutation"`` attaches ``mutation_status``.
        mutation_gene: Gene whose mutation status colours points.
        gene_type: Mutation filter mode when *color_by* is ``"mutation"``.
        db_path: Optional DuckDB path override.

    Returns:
        DataFrame with ``model_id``, *gene_x*, *gene_y*, and optionally
        ``oncotree_lineage`` or ``mutation_status``.
    """
    vx = extract_gene_values(assay, gene_x, db_path=db_path)
    vy = extract_gene_values(assay, gene_y, db_path=db_path)
    merged = vx.join(vy, on="model_id")

    if color_by == "lineage":
        merged = join_models(merged, db_path=db_path)

    df = merged.collect().drop_nulls([gene_x, gene_y])

    if color_by == "mutation" and mutation_gene is not None:
        mutant_ids = get_mutant_model_ids(
            mutation_gene, gene_type, db_path=db_path
        )
        df = label_mutation_status(df, mutant_ids)

    return df


def plot_gene_gene(
    df: pl.DataFrame,
    gene_x: str,
    gene_y: str,
    assay: Literal["dependency", "expression"] = "dependency",
    *,
    color_by: Literal["lineage", "mutation", "none"] = "none",
    method: Literal["pearson", "spearman"] = "pearson",
    figsize: tuple[float, float] | None = None,
    ax: plt.Axes | None = None,
) -> tuple[plt.Figure, plt.Axes]:
    """Scatter plot of two genes with regression line and r annotation.

    Args:
        df: DataFrame from :func:`analyse_gene_gene`.
        gene_x: Column for x-axis.
        gene_y: Column for y-axis.
        assay: Controls axis labels.
        color_by: ``"lineage"`` or ``"mutation"`` colors points by the
            corresponding grouping column.
        method: Correlation method shown in annotation.
        figsize: Override figure size.
        ax: Optional pre-existing axes.

    Returns:
        Tuple of (Figure, Axes).
    """
    use_style()
    fig, ax = make_axes(figsize, ax)
    pdf = df.to_pandas()

    hue = None
    palette = None
    if color_by == "lineage":
        hue = "oncotree_lineage"
    elif color_by == "mutation":
        hue = "mutation_status"
        palette = {"WT": COLOR.GREY.value, "mutant": COLOR.PINK.value}

    sns.scatterplot(
        data=pdf,
        x=gene_x,
        y=gene_y,
        hue=hue,
        palette=palette,
        s=4,
        alpha=0.6,
        linewidth=0,
        ax=ax,
        legend=hue is not None,
    )

    # regression line
    x, y = pdf[gene_x].values, pdf[gene_y].values
    m, b = np.polyfit(x, y, 1)
    xs = np.linspace(x.min(), x.max(), 100)
    ax.plot(xs, m * xs + b, lw=0.7, color=COLOR.PINK.value, zorder=5)

    # correlation
    cs = pearson_spearman(x, y)
    r_key = f"{method}_r"
    p_key = f"{method}_p"
    annotate_correlation(ax, cs[r_key], cs[p_key], label=method[0])

    ax.set_xlabel(assay_ylabel(gene_x, assay))
    ax.set_ylabel(assay_ylabel(gene_y, assay))
    ax.set_title(f"{gene_x} vs {gene_y}", fontsize=7)

    if hue is not None and ax.get_legend() is not None:
        ax.legend(fontsize=4, markerscale=0.5, loc="best", frameon=True)

    fig.tight_layout()
    return fig, ax


def gene_gene_scatter(
    gene_x: str,
    gene_y: str,
    assay: Literal["dependency", "expression"] = "dependency",
    *,
    color_by: Literal["lineage", "mutation", "none"] = "none",
    mutation_gene: str | None = None,
    gene_type: Literal["suppressor", "oncogene"] = "suppressor",
    method: Literal["pearson", "spearman"] = "pearson",
    db_path: str | Path | None = None,
    figsize: tuple[float, float] | None = None,
) -> plt.Figure:
    """Convenience: analyse and plot gene-gene scatter in one step."""
    df = analyse_gene_gene(
        gene_x,
        gene_y,
        assay,
        color_by=color_by,
        mutation_gene=mutation_gene,
        gene_type=gene_type,
        db_path=db_path,
    )
    fig, _ = plot_gene_gene(
        df,
        gene_x,
        gene_y,
        assay,
        color_by=color_by,
        method=method,
        figsize=figsize,
    )
    return fig


# ---------------------------------------------------------------------------
# 2. Top correlated genes
# ---------------------------------------------------------------------------


def analyse_top_correlations(
    query_gene: str,
    assay: Literal["dependency", "expression"] = "dependency",
    *,
    n: int = 20,
    method: Literal["pearson", "spearman"] = "pearson",
    db_path: str | Path | None = None,
) -> pl.DataFrame:
    """Compute correlations between *query_gene* and all other genes.

    Args:
        query_gene: Hugo symbol of the query gene.
        assay: Assay to use.
        n: Number of top (absolute) correlations to return.
        method: ``"pearson"`` or ``"spearman"``.
        db_path: Optional DuckDB path override.

    Returns:
        DataFrame with columns ``gene``, ``r``, ``p`` sorted by
        descending absolute ``r``, limited to *n* rows.
    """
    data = scan_assay(assay, db_path=db_path).collect()
    names = data.columns
    query_col = resolve_gene_column(names, query_gene)
    query_vals = data.get_column(query_col).to_numpy()

    gene_cols = [
        c for c in names
        if c != "model_id" and c != query_col
    ]

    records: list[dict[str, object]] = []
    corr_fn = stats.pearsonr if method == "pearson" else stats.spearmanr
    for col in gene_cols:
        vals = data.get_column(col).to_numpy()
        mask = np.isfinite(query_vals) & np.isfinite(vals)
        if mask.sum() < 3:
            continue
        r, p = corr_fn(query_vals[mask], vals[mask])
        # strip entrez suffix for display
        gene_label = col.split(" (")[0] if " (" in col else col
        records.append({"gene": gene_label, "r": float(r), "p": float(p)})

    result = pl.DataFrame(records)
    return (
        result
        .with_columns(pl.col("r").abs().alias("abs_r"))
        .sort("abs_r", descending=True)
        .head(n)
        .drop("abs_r")
    )


def plot_top_correlations(
    df: pl.DataFrame,
    query_gene: str,
    *,
    method: Literal["pearson", "spearman"] = "pearson",
    figsize: tuple[float, float] | None = None,
    ax: plt.Axes | None = None,
) -> tuple[plt.Figure, plt.Axes]:
    """Horizontal lollipop chart of top correlated genes.

    Args:
        df: DataFrame from :func:`analyse_top_correlations`.
        query_gene: Used for the title.
        method: Label for axis.
        figsize: Override figure size.
        ax: Optional pre-existing axes.

    Returns:
        Tuple of (Figure, Axes).
    """
    use_style()
    n = df.height
    fig, ax = make_axes(figsize, ax, default_figsize=(3.0, max(1.5, n * 0.15)))

    genes = df.get_column("gene").to_list()[::-1]
    r_vals = df.get_column("r").to_list()[::-1]
    colors = [
        COLOR.PINK.value if v < 0 else COLOR.BLUE.value for v in r_vals
    ]

    y_pos = range(len(genes))
    ax.barh(y_pos, r_vals, height=0.6, color=colors, edgecolor="none")
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(genes)
    ax.axvline(0, lw=0.5, color=COLOR.DARKGREY.value)
    ax.set_xlabel(f"{method} r")
    ax.set_title(f"Top correlates of {query_gene}", fontsize=7)

    fig.tight_layout()
    return fig, ax


def top_correlations(
    query_gene: str,
    assay: Literal["dependency", "expression"] = "dependency",
    *,
    n: int = 20,
    method: Literal["pearson", "spearman"] = "pearson",
    db_path: str | Path | None = None,
    figsize: tuple[float, float] | None = None,
) -> plt.Figure:
    """Convenience: analyse and plot top correlations in one step."""
    df = analyse_top_correlations(
        query_gene, assay, n=n, method=method, db_path=db_path,
    )
    fig, _ = plot_top_correlations(
        df, query_gene, method=method, figsize=figsize,
    )
    return fig
