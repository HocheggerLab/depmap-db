"""Drug sensitivity plots: by lineage, by mutation, and drug-dependency correlation.

Covers requested plot families:
5) Drug sensitivity by lineage.
6) Drug sensitivity by mutation status.
7) Drug-dependency correlation scatter.

Schema notes:
    - PRISM primary is wide (model_ids as columns). The ``drug_primary_long``
      curated dataset provides a tidy long-format view.
    - PRISM secondary (``drug_response_secondary``) is already long with
      ``auc``, ``ec50``, ``ic50`` columns.
    - Drug identification uses ``broad_id`` + ``compound_name``.
    - If the ``drug_secondary_enriched`` curated dataset is unavailable, these
      functions fall back to ``drug_response_secondary`` directly.
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
    scan_drug_response_secondary,
    scan_models,
)

from ._helpers import (
    annotate_correlation,
    extract_gene_values,
    get_mutant_model_ids,
    label_mutation_status,
    make_axes,
    pearson_spearman,
    use_style,
)
from .colors import COLOR

# ---------------------------------------------------------------------------
# 5. Drug sensitivity by lineage
# ---------------------------------------------------------------------------


def analyse_drug_by_lineage(
    compound_name: str,
    *,
    metric: Literal["auc", "ic50", "ec50"] = "auc",
    db_path: str | Path | None = None,
) -> pl.DataFrame:
    """Build a DataFrame of drug sensitivity values by lineage.

    Uses the ``drug_response_secondary`` table which has per-model
    dose-response summary metrics.

    Args:
        compound_name: Name of the compound (case-insensitive match).
        metric: Which summary metric to use.
        db_path: Optional DuckDB path override.

    Returns:
        DataFrame with ``model_id``, ``oncotree_lineage``, *metric*.
    """
    drug = (
        scan_drug_response_secondary(db_path=db_path)
        .filter(pl.col("compound_name").str.to_lowercase()
                == compound_name.lower())
        .select("model_id", pl.col(metric))
    )
    models = scan_models(db_path=db_path).select(
        "model_id", "oncotree_lineage"
    )
    return (
        models.join(drug, on="model_id")
        .collect()
        .drop_nulls([metric, "oncotree_lineage"])
    )


def plot_drug_by_lineage(
    df: pl.DataFrame,
    compound_name: str,
    *,
    metric: Literal["auc", "ic50", "ec50"] = "auc",
    min_models: int = 3,
    figsize: tuple[float, float] | None = None,
    ax: plt.Axes | None = None,
) -> tuple[plt.Figure, plt.Axes]:
    """Box + strip plot of drug sensitivity by lineage.

    Args:
        df: DataFrame from :func:`analyse_drug_by_lineage`.
        compound_name: For the title.
        metric: Column to plot.
        min_models: Exclude lineages below this threshold.
        figsize: Override figure size.
        ax: Optional pre-existing axes.

    Returns:
        Tuple of (Figure, Axes).
    """
    use_style()

    filtered = df.filter(
        pl.col("oncotree_lineage").count().over("oncotree_lineage")
        >= min_models
    )
    order = (
        filtered.group_by("oncotree_lineage")
        .agg(pl.col(metric).median().alias("median"))
        .sort("median")
        .get_column("oncotree_lineage")
        .to_list()
    )

    fig, ax = make_axes(
        figsize, ax,
        default_figsize=(max(4.0, len(order) * 0.25), 2.5),
    )
    pdf = filtered.to_pandas()

    sns.boxplot(
        data=pdf,
        x="oncotree_lineage",
        y=metric,
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
        y=metric,
        order=order,
        color=COLOR.TURQUOISE.value,
        size=1.5,
        alpha=0.5,
        jitter=True,
        ax=ax,
    )

    ax.set_xlabel("")
    ax.set_ylabel(metric.upper())
    plt.setp(
        ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor"
    )
    ax.set_title(f"{compound_name} — {metric} by lineage", fontsize=7)

    fig.tight_layout()
    return fig, ax


def drug_by_lineage(
    compound_name: str,
    *,
    metric: Literal["auc", "ic50", "ec50"] = "auc",
    min_models: int = 3,
    db_path: str | Path | None = None,
    figsize: tuple[float, float] | None = None,
) -> plt.Figure:
    """Convenience: analyse and plot drug sensitivity by lineage."""
    df = analyse_drug_by_lineage(
        compound_name, metric=metric, db_path=db_path,
    )
    fig, _ = plot_drug_by_lineage(
        df, compound_name, metric=metric,
        min_models=min_models, figsize=figsize,
    )
    return fig


# ---------------------------------------------------------------------------
# 6. Drug sensitivity by mutation status
# ---------------------------------------------------------------------------


def analyse_drug_by_mutation(
    compound_name: str,
    mutation_gene: str,
    *,
    metric: Literal["auc", "ic50", "ec50"] = "auc",
    gene_type: Literal["suppressor", "oncogene"] = "suppressor",
    db_path: str | Path | None = None,
) -> pl.DataFrame:
    """Build DataFrame of drug sensitivity split by mutation status.

    Args:
        compound_name: Compound name (case-insensitive).
        mutation_gene: Gene whose mutation status defines groups.
        metric: Summary metric.
        gene_type: Mutation filter mode.
        db_path: Optional DuckDB path override.

    Returns:
        DataFrame with ``model_id``, *metric*, ``mutation_status``.
    """
    drug = (
        scan_drug_response_secondary(db_path=db_path)
        .filter(pl.col("compound_name").str.to_lowercase()
                == compound_name.lower())
        .select("model_id", pl.col(metric))
        .collect()
        .drop_nulls(metric)
    )
    mutant_ids = get_mutant_model_ids(
        mutation_gene, gene_type, db_path=db_path
    )
    return label_mutation_status(drug, mutant_ids)


def plot_drug_by_mutation(
    df: pl.DataFrame,
    compound_name: str,
    mutation_gene: str,
    *,
    metric: Literal["auc", "ic50", "ec50"] = "auc",
    figsize: tuple[float, float] | None = None,
    ax: plt.Axes | None = None,
) -> tuple[plt.Figure, plt.Axes]:
    """Box + strip plot of drug sensitivity by mutation status.

    Args:
        df: DataFrame from :func:`analyse_drug_by_mutation`.
        compound_name: For the title.
        mutation_gene: For the title.
        metric: Column to plot.
        figsize: Override figure size.
        ax: Optional pre-existing axes.

    Returns:
        Tuple of (Figure, Axes).
    """
    use_style()
    fig, ax = make_axes(figsize, ax, default_figsize=(2.0, 2.5))
    pdf = df.to_pandas()

    order = ["WT", "mutant"]
    palette = {"WT": COLOR.GREY.value, "mutant": COLOR.PINK.value}

    sns.boxplot(
        data=pdf,
        x="mutation_status",
        y=metric,
        order=order,
        palette=palette,
        fliersize=0,
        linewidth=0.5,
        width=0.6,
        ax=ax,
    )
    sns.stripplot(
        data=pdf,
        x="mutation_status",
        y=metric,
        order=order,
        color=COLOR.DARKGREY.value,
        size=1.5,
        alpha=0.4,
        jitter=True,
        ax=ax,
    )

    # stats annotation
    wt = pdf.loc[pdf["mutation_status"] == "WT", metric]
    mut = pdf.loc[pdf["mutation_status"] == "mutant", metric]
    if len(wt) > 0 and len(mut) > 0:
        _, pval = stats.mannwhitneyu(wt, mut, alternative="two-sided")
        ax.set_title(
            f"{compound_name} — {mutation_gene} mut vs WT\n"
            f"(n={len(mut)} vs {len(wt)}, p={pval:.2e})",
            fontsize=7,
        )
    else:
        ax.set_title(
            f"{compound_name} — {mutation_gene} mut vs WT", fontsize=7
        )

    ax.set_xlabel("")
    ax.set_ylabel(metric.upper())
    fig.tight_layout()
    return fig, ax


def drug_by_mutation(
    compound_name: str,
    mutation_gene: str,
    *,
    metric: Literal["auc", "ic50", "ec50"] = "auc",
    gene_type: Literal["suppressor", "oncogene"] = "suppressor",
    db_path: str | Path | None = None,
    figsize: tuple[float, float] | None = None,
) -> plt.Figure:
    """Convenience: analyse and plot drug sensitivity by mutation."""
    df = analyse_drug_by_mutation(
        compound_name, mutation_gene,
        metric=metric, gene_type=gene_type, db_path=db_path,
    )
    fig, _ = plot_drug_by_mutation(
        df, compound_name, mutation_gene,
        metric=metric, figsize=figsize,
    )
    return fig


# ---------------------------------------------------------------------------
# 7. Drug-dependency correlation scatter
# ---------------------------------------------------------------------------


def analyse_drug_dependency(
    compound_name: str,
    target_gene: str,
    *,
    metric: Literal["auc", "ic50", "ec50"] = "auc",
    db_path: str | Path | None = None,
) -> pl.DataFrame:
    """Build DataFrame pairing drug sensitivity with gene dependency.

    Args:
        compound_name: Compound name (case-insensitive).
        target_gene: Gene whose dependency to correlate.
        metric: Drug response metric.
        db_path: Optional DuckDB path override.

    Returns:
        DataFrame with ``model_id``, *metric*, ``dependency``.
    """
    drug = (
        scan_drug_response_secondary(db_path=db_path)
        .filter(pl.col("compound_name").str.to_lowercase()
                == compound_name.lower())
        .select("model_id", pl.col(metric))
    )
    dep = extract_gene_values(
        "dependency", target_gene, db_path=db_path
    ).rename({target_gene: "dependency"})

    return (
        drug.join(dep, on="model_id")
        .collect()
        .drop_nulls([metric, "dependency"])
    )


def plot_drug_dependency(
    df: pl.DataFrame,
    compound_name: str,
    target_gene: str,
    *,
    metric: Literal["auc", "ic50", "ec50"] = "auc",
    method: Literal["pearson", "spearman"] = "pearson",
    figsize: tuple[float, float] | None = None,
    ax: plt.Axes | None = None,
) -> tuple[plt.Figure, plt.Axes]:
    """Scatter: drug sensitivity vs gene dependency with correlation.

    Args:
        df: DataFrame from :func:`analyse_drug_dependency`.
        compound_name: For the title.
        target_gene: For axis label.
        metric: Drug response metric column.
        method: Correlation method for annotation.
        figsize: Override figure size.
        ax: Optional pre-existing axes.

    Returns:
        Tuple of (Figure, Axes).
    """
    use_style()
    fig, ax = make_axes(figsize, ax)
    pdf = df.to_pandas()

    ax.scatter(
        pdf[metric],
        pdf["dependency"],
        s=4,
        alpha=0.5,
        color=COLOR.BLUE.value,
        linewidth=0,
    )

    x, y = pdf[metric].values, pdf["dependency"].values
    m, b = np.polyfit(x, y, 1)
    xs = np.linspace(x.min(), x.max(), 100)
    ax.plot(xs, m * xs + b, lw=0.7, color=COLOR.PINK.value, zorder=5)

    cs = pearson_spearman(x, y)
    r_key = f"{method}_r"
    p_key = f"{method}_p"
    annotate_correlation(ax, cs[r_key], cs[p_key], label=method[0])

    ax.set_xlabel(f"{compound_name} {metric.upper()}")
    ax.set_ylabel(f"{target_gene} gene effect")
    ax.set_title(
        f"{compound_name} vs {target_gene} dependency", fontsize=7
    )
    ax.axhline(-1, ls="--", lw=0.5, color=COLOR.DARKGREY.value, zorder=0)

    fig.tight_layout()
    return fig, ax


def drug_dependency_scatter(
    compound_name: str,
    target_gene: str,
    *,
    metric: Literal["auc", "ic50", "ec50"] = "auc",
    method: Literal["pearson", "spearman"] = "pearson",
    db_path: str | Path | None = None,
    figsize: tuple[float, float] | None = None,
) -> plt.Figure:
    """Convenience: analyse and plot drug-dependency correlation."""
    df = analyse_drug_dependency(
        compound_name, target_gene, metric=metric, db_path=db_path,
    )
    fig, _ = plot_drug_dependency(
        df, compound_name, target_gene,
        metric=metric, method=method, figsize=figsize,
    )
    return fig
