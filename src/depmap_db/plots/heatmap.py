"""Selectivity waterfall and multi-gene heatmap plots.

Covers requested plot families:
10) Selectivity plot / ranked dependency waterfall.
11) Multi-gene heatmap of dependency or expression for a gene set.
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import polars as pl
import seaborn as sns

from ._helpers import (
    join_models,
    make_axes,
    resolve_gene_column,
    scan_assay,
    use_style,
)
from .colors import COLOR

# ---------------------------------------------------------------------------
# 10. Selectivity / ranked dependency waterfall
# ---------------------------------------------------------------------------


def analyse_selectivity(
    gene: str,
    assay: Literal["dependency", "expression"] = "dependency",
    *,
    db_path: str | Path | None = None,
) -> pl.DataFrame:
    """Rank all cell lines by assay value for *gene*.

    Args:
        gene: Hugo symbol.
        assay: Which assay.
        db_path: Optional DuckDB path override.

    Returns:
        DataFrame with ``model_id``, ``oncotree_lineage``, ``value``,
        ``rank`` sorted from lowest to highest value.
    """
    data = scan_assay(assay, db_path=db_path)
    col = resolve_gene_column(data.collect_schema().names(), gene)
    values = data.select("model_id", pl.col(col)).rename({col: "value"})
    merged = join_models(values, db_path=db_path)

    return (
        merged.collect()
        .drop_nulls("value")
        .sort("value")
        .with_row_index("rank")
        .with_columns(pl.col("rank").cast(pl.Int64))
    )


def plot_selectivity(
    df: pl.DataFrame,
    gene: str,
    assay: Literal["dependency", "expression"] = "dependency",
    *,
    highlight_lineage: str | None = None,
    figsize: tuple[float, float] | None = None,
    ax: plt.Axes | None = None,
) -> tuple[plt.Figure, plt.Axes]:
    """Waterfall bar chart ranking all cell lines for a gene.

    Args:
        df: DataFrame from :func:`analyse_selectivity`.
        gene: Gene symbol for labels.
        assay: Controls y-axis label.
        highlight_lineage: If given, highlight this lineage in colour.
        figsize: Override figure size.
        ax: Optional pre-existing axes.

    Returns:
        Tuple of (Figure, Axes).
    """
    use_style()
    fig, ax = make_axes(figsize, ax, default_figsize=(4.0, 2.0))
    pdf = df.to_pandas()

    if highlight_lineage:
        colors = [
            COLOR.PINK.value if lin == highlight_lineage else COLOR.GREY.value
            for lin in pdf["oncotree_lineage"]
        ]
    else:
        colors = COLOR.BLUE.value

    ax.bar(
        pdf["rank"], pdf["value"],
        width=1.0, color=colors, edgecolor="none",
    )

    if assay == "dependency":
        ax.axhline(-1, ls="--", lw=0.5, color=COLOR.DARKGREY.value, zorder=0)
        ax.set_ylabel(f"{gene} gene effect")
    else:
        ax.set_ylabel(f"{gene} log2(TPM+1)")

    ax.set_xlabel(f"Cell lines (n={df.height})")
    title = f"{gene} — {assay} selectivity"
    if highlight_lineage:
        title += f" ({highlight_lineage} highlighted)"
    ax.set_title(title, fontsize=7)
    ax.set_xticks([])

    fig.tight_layout()
    return fig, ax


def selectivity_plot(
    gene: str,
    assay: Literal["dependency", "expression"] = "dependency",
    *,
    highlight_lineage: str | None = None,
    db_path: str | Path | None = None,
    figsize: tuple[float, float] | None = None,
) -> plt.Figure:
    """Convenience: analyse and plot selectivity in one step."""
    df = analyse_selectivity(gene, assay, db_path=db_path)
    fig, _ = plot_selectivity(
        df, gene, assay,
        highlight_lineage=highlight_lineage, figsize=figsize,
    )
    return fig


# ---------------------------------------------------------------------------
# 11. Multi-gene heatmap
# ---------------------------------------------------------------------------


def analyse_gene_heatmap(
    genes: list[str],
    assay: Literal["dependency", "expression"] = "dependency",
    *,
    db_path: str | Path | None = None,
) -> pl.DataFrame:
    """Build a matrix of assay values for a set of genes across cell lines.

    Args:
        genes: List of Hugo symbols.
        assay: Which assay.
        db_path: Optional DuckDB path override.

    Returns:
        DataFrame with ``model_id`` and one column per gene (using the
        short Hugo symbol as column name).
    """
    data = scan_assay(assay, db_path=db_path)
    names = data.collect_schema().names()

    select_cols = ["model_id"]
    rename_map: dict[str, str] = {}
    for gene in genes:
        col = resolve_gene_column(names, gene)
        select_cols.append(col)
        if col != gene:
            rename_map[col] = gene

    df = data.select(select_cols).collect()
    if rename_map:
        df = df.rename(rename_map)
    return df


def plot_gene_heatmap(
    df: pl.DataFrame,
    genes: list[str],
    assay: Literal["dependency", "expression"] = "dependency",
    *,
    cluster_rows: bool = True,
    cluster_cols: bool = True,
    figsize: tuple[float, float] | None = None,
) -> sns.matrix.ClusterGrid:
    """Clustered heatmap of a gene set across cell lines.

    Args:
        df: DataFrame from :func:`analyse_gene_heatmap`.
        genes: Gene columns to include.
        assay: Controls the colour bar label.
        cluster_rows: Cluster cell lines.
        cluster_cols: Cluster genes.
        figsize: Override figure size.

    Returns:
        seaborn ClusterGrid object (access ``.fig`` for the Figure).
    """
    use_style()

    matrix = df.select(genes).to_pandas()
    matrix.index = df.get_column("model_id").to_list()
    matrix = matrix.dropna(how="all")

    if figsize is None:
        figsize = (max(2.5, len(genes) * 0.3), max(3.0, len(matrix) * 0.01))

    cbar_label = "gene effect" if assay == "dependency" else "log2(TPM+1)"

    cg = sns.clustermap(
        matrix,
        cmap="RdBu_r" if assay == "dependency" else "viridis",
        center=0 if assay == "dependency" else None,
        row_cluster=cluster_rows,
        col_cluster=cluster_cols,
        figsize=figsize,
        linewidths=0,
        xticklabels=True,
        yticklabels=False,
        cbar_kws={"label": cbar_label},
    )
    cg.ax_heatmap.set_xlabel("")
    cg.ax_heatmap.set_ylabel("")
    plt.setp(
        cg.ax_heatmap.get_xticklabels(),
        rotation=45, ha="right", rotation_mode="anchor", fontsize=5,
    )

    return cg


def gene_heatmap(
    genes: list[str],
    assay: Literal["dependency", "expression"] = "dependency",
    *,
    cluster_rows: bool = True,
    cluster_cols: bool = True,
    db_path: str | Path | None = None,
    figsize: tuple[float, float] | None = None,
) -> plt.Figure:
    """Convenience: analyse and plot multi-gene heatmap in one step."""
    df = analyse_gene_heatmap(genes, assay, db_path=db_path)
    cg = plot_gene_heatmap(
        df, genes, assay,
        cluster_rows=cluster_rows, cluster_cols=cluster_cols,
        figsize=figsize,
    )
    return cg.fig
