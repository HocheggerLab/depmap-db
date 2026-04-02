"""Expression vs dependency scatter and biomarker volcano.

Covers requested plot families:
3) Expression vs dependency scatter for a single gene.
4) Biomarker volcano: target gene dependency vs all genes expression.
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import seaborn as sns
from scipy import stats

from depmap_db import (
    scan_gene_effects_wide,
    scan_gene_expression_wide,
)

from ._helpers import (
    annotate_correlation,
    extract_gene_values,
    get_mutant_model_ids,
    join_models,
    label_mutation_status,
    make_axes,
    pearson_spearman,
    resolve_gene_column,
    use_style,
)
from .colors import COLOR

# ---------------------------------------------------------------------------
# 3. Expression vs dependency scatter
# ---------------------------------------------------------------------------


def analyse_expr_vs_dep(
    gene: str,
    *,
    color_by: Literal["lineage", "mutation", "none"] = "none",
    mutation_gene: str | None = None,
    gene_type: Literal["suppressor", "oncogene"] = "suppressor",
    db_path: str | Path | None = None,
) -> pl.DataFrame:
    """Build expression-vs-dependency DataFrame for a single gene.

    Args:
        gene: Hugo symbol.
        color_by: ``"lineage"`` or ``"mutation"`` (requires *mutation_gene*).
        mutation_gene: Gene whose mutation status colours points.
        gene_type: Mutation filter mode when *color_by* is ``"mutation"``.
        db_path: Optional DuckDB path override.

    Returns:
        DataFrame with ``model_id``, ``expression``, ``dependency`` and
        optionally ``oncotree_lineage`` or ``mutation_status``.
    """
    expr = extract_gene_values("expression", gene, db_path=db_path).rename(
        {gene: "expression"}
    )
    dep = extract_gene_values("dependency", gene, db_path=db_path).rename(
        {gene: "dependency"}
    )
    merged = expr.join(dep, on="model_id")

    if color_by == "lineage":
        merged = join_models(merged, db_path=db_path)

    df = merged.collect().drop_nulls(["expression", "dependency"])

    if color_by == "mutation" and mutation_gene is not None:
        mutant_ids = get_mutant_model_ids(
            mutation_gene, gene_type, db_path=db_path
        )
        df = label_mutation_status(df, mutant_ids)

    return df


def plot_expr_vs_dep(
    df: pl.DataFrame,
    gene: str,
    *,
    color_by: Literal["lineage", "mutation", "none"] = "none",
    method: Literal["pearson", "spearman"] = "pearson",
    figsize: tuple[float, float] | None = None,
    ax: plt.Axes | None = None,
) -> tuple[plt.Figure, plt.Axes]:
    """Scatter plot of expression vs dependency for one gene.

    Args:
        df: DataFrame from :func:`analyse_expr_vs_dep`.
        gene: Gene symbol (for labels/title).
        color_by: Colour scheme matching what was used in analysis.
        method: Correlation method for annotation.
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
        x="expression",
        y="dependency",
        hue=hue,
        palette=palette,
        s=4,
        alpha=0.6,
        linewidth=0,
        ax=ax,
    )

    x, y = pdf["expression"].values, pdf["dependency"].values
    cs = pearson_spearman(x, y)
    r_key = f"{method}_r"
    p_key = f"{method}_p"
    annotate_correlation(ax, cs[r_key], cs[p_key], label=method[0])

    ax.set_xlabel(f"{gene} log2(TPM+1)")
    ax.set_ylabel(f"{gene} gene effect")
    ax.set_title(f"{gene} — expression vs dependency", fontsize=7)
    ax.axhline(-1, ls="--", lw=0.5, color=COLOR.DARKGREY.value, zorder=0)

    if hue and ax.get_legend() is not None:
        ax.legend(fontsize=4, markerscale=0.5, loc="best", frameon=True)

    fig.tight_layout()
    return fig, ax


def expr_vs_dep(
    gene: str,
    *,
    color_by: Literal["lineage", "mutation", "none"] = "none",
    mutation_gene: str | None = None,
    gene_type: Literal["suppressor", "oncogene"] = "suppressor",
    method: Literal["pearson", "spearman"] = "pearson",
    db_path: str | Path | None = None,
    figsize: tuple[float, float] | None = None,
) -> plt.Figure:
    """Convenience: analyse and plot expression vs dependency in one step."""
    df = analyse_expr_vs_dep(
        gene,
        color_by=color_by,
        mutation_gene=mutation_gene,
        gene_type=gene_type,
        db_path=db_path,
    )
    fig, _ = plot_expr_vs_dep(
        df, gene, color_by=color_by, method=method, figsize=figsize,
    )
    return fig


# ---------------------------------------------------------------------------
# 4. Biomarker volcano
# ---------------------------------------------------------------------------


def analyse_biomarker_volcano(
    target_gene: str,
    *,
    db_path: str | Path | None = None,
) -> pl.DataFrame:
    """Correlate *target_gene* dependency with all genes' expression.

    Returns a DataFrame with columns ``gene``, ``effect_size`` (Pearson r),
    ``p``, and ``neg_log10_p``.

    Note:
        This scans the full expression matrix and can be slow for the first
        call on a large database.
    """
    dep_data = scan_gene_effects_wide(db_path=db_path).collect()
    dep_col = resolve_gene_column(dep_data.columns, target_gene)
    dep_vals = dep_data.select("model_id", pl.col(dep_col)).rename(
        {dep_col: "dep"}
    )

    expr_data = scan_gene_expression_wide(db_path=db_path).collect()
    # align on common models
    merged = dep_vals.join(
        expr_data, on="model_id", how="inner"
    ).drop_nulls("dep")

    dep_arr = merged.get_column("dep").to_numpy()
    expr_cols = [
        c for c in merged.columns
        if c not in {"model_id", "dep"}
        and merged.schema[c].is_numeric()
    ]

    records: list[dict[str, object]] = []
    for col in expr_cols:
        vals = merged.get_column(col).to_numpy().astype(float)
        mask = np.isfinite(dep_arr) & np.isfinite(vals)
        if mask.sum() < 3:
            continue
        r, p = stats.pearsonr(dep_arr[mask], vals[mask])
        gene_label = col.split(" (")[0] if " (" in col else col
        nlp = -np.log10(max(p, 1e-300))
        records.append({
            "gene": gene_label,
            "effect_size": float(r),
            "p": float(p),
            "neg_log10_p": float(nlp),
        })

    return pl.DataFrame(records).sort("neg_log10_p", descending=True)


def plot_biomarker_volcano(
    df: pl.DataFrame,
    target_gene: str,
    *,
    effect_threshold: float = 0.2,
    p_threshold: float = 1e-5,
    figsize: tuple[float, float] | None = None,
    ax: plt.Axes | None = None,
) -> tuple[plt.Figure, plt.Axes]:
    """Volcano plot: effect size vs -log10(p) for biomarker discovery.

    Args:
        df: DataFrame from :func:`analyse_biomarker_volcano`.
        target_gene: Used for the title.
        effect_threshold: Absolute r threshold for coloring.
        p_threshold: P-value threshold for coloring.
        figsize: Override figure size.
        ax: Optional pre-existing axes.

    Returns:
        Tuple of (Figure, Axes).
    """
    use_style()
    fig, ax = make_axes(figsize, ax, default_figsize=(3.5, 2.5))
    pdf = df.to_pandas()

    nlp_thresh = -np.log10(p_threshold)
    sig = (pdf["effect_size"].abs() >= effect_threshold) & (
        pdf["neg_log10_p"] >= nlp_thresh
    )

    ax.scatter(
        pdf.loc[~sig, "effect_size"],
        pdf.loc[~sig, "neg_log10_p"],
        s=2,
        alpha=0.3,
        color=COLOR.GREY.value,
        linewidth=0,
    )
    ax.scatter(
        pdf.loc[sig, "effect_size"],
        pdf.loc[sig, "neg_log10_p"],
        s=3,
        alpha=0.7,
        color=COLOR.PINK.value,
        linewidth=0,
    )

    ax.axhline(nlp_thresh, ls="--", lw=0.4, color=COLOR.DARKGREY.value)
    ax.axvline(-effect_threshold, ls="--", lw=0.4, color=COLOR.DARKGREY.value)
    ax.axvline(effect_threshold, ls="--", lw=0.4, color=COLOR.DARKGREY.value)

    ax.set_xlabel("Pearson r (expression vs dependency)")
    ax.set_ylabel("-log10(p)")
    ax.set_title(
        f"Biomarkers of {target_gene} dependency", fontsize=7
    )

    fig.tight_layout()
    return fig, ax


def biomarker_volcano(
    target_gene: str,
    *,
    effect_threshold: float = 0.2,
    p_threshold: float = 1e-5,
    db_path: str | Path | None = None,
    figsize: tuple[float, float] | None = None,
) -> plt.Figure:
    """Convenience: analyse and plot biomarker volcano in one step."""
    df = analyse_biomarker_volcano(target_gene, db_path=db_path)
    fig, _ = plot_biomarker_volcano(
        df, target_gene,
        effect_threshold=effect_threshold,
        p_threshold=p_threshold,
        figsize=figsize,
    )
    return fig
