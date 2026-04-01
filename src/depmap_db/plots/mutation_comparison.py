"""Compare gene dependency or expression by mutation status of another gene."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import polars as pl
import seaborn as sns
from scipy import stats

from depmap_db import (
    scan_gene_effects_wide,
    scan_gene_expression_wide,
    scan_models,
    scan_mutations,
)

from .colors import COLOR

_STYLE = Path(__file__).parent / "hhlab_style01.mplstyle"


def get_mutant_models(
    mutation_gene: str,
    gene_type: Literal["suppressor", "oncogene"] = "suppressor",
) -> pl.LazyFrame:
    """Return model_ids carrying damaging mutations in *mutation_gene*.

    Args:
        mutation_gene: Hugo symbol of the mutated gene to filter on.
        gene_type: Determines which impact flags to use.
            ``"suppressor"`` — likely_lof OR tumor_suppressor_high_impact.
            ``"oncogene"``   — hotspot OR oncogene_high_impact OR hess_driver.

    Returns:
        LazyFrame with unique ``model_id`` values for mutant models.
    """
    muts = scan_mutations().filter(pl.col("hugo_symbol") == mutation_gene)

    if gene_type == "suppressor":
        muts = muts.filter(
            pl.col("likely_lof") | pl.col("tumor_suppressor_high_impact")
        )
    else:
        muts = muts.filter(
            pl.col("hotspot")
            | pl.col("oncogene_high_impact")
            | pl.col("hess_driver")
        )

    return muts.select("model_id").unique()


def analyse_by_mutation(
    target_gene: str,
    mutation_gene: str,
    assay: Literal["dependency", "expression"] = "dependency",
    gene_type: Literal["suppressor", "oncogene"] = "suppressor",
    min_models: int = 20,
) -> pl.DataFrame:
    """Build a DataFrame comparing *target_gene* values in mutant vs WT models.

    Args:
        target_gene: Gene whose dependency/expression to measure.
        mutation_gene: Gene whose mutation status defines the groups.
        assay: Which assay to pull values from.
        gene_type: Mutation filter mode (see :func:`get_mutant_models`).
        min_models: Minimum mutant models required; raises if below.

    Returns:
        DataFrame with columns: ``model_id``, ``oncotree_lineage``,
        *target_gene*, ``mutation_status`` (``"mutant"`` / ``"WT"``).

    Raises:
        ValueError: If fewer than *min_models* mutant models are found.
    """
    # assay data
    data = (
        scan_gene_effects_wide()
        if assay == "dependency"
        else scan_gene_expression_wide()
    )
    names = data.collect_schema().names()
    col = (
        target_gene
        if target_gene in names
        else next(c for c in names if c.startswith(f"{target_gene} ("))
    )

    models = scan_models().select("model_id", "oncotree_lineage")
    values = data.select("model_id", pl.col(col)).rename({col: target_gene})

    # mutant model ids
    mutant_ids_df = get_mutant_models(
        mutation_gene, gene_type=gene_type
    ).collect()
    n_mutant = mutant_ids_df.height

    if n_mutant < min_models:
        raise ValueError(
            f"Only {n_mutant} models with damaging {mutation_gene} mutations "
            f"({gene_type} filter). Minimum is {min_models}."
        )

    mutant_series = mutant_ids_df.get_column("model_id")

    # join and label
    df = (
        models.join(values, on="model_id")
        .collect()
        .drop_nulls(target_gene)
        .with_columns(
            pl.when(pl.col("model_id").is_in(mutant_series))
            .then(pl.lit("mutant"))
            .otherwise(pl.lit("WT"))
            .alias("mutation_status")
        )
    )

    return df


def plot_mutation_comparison(
    df: pl.DataFrame,
    target_gene: str,
    mutation_gene: str,
    assay: Literal["dependency", "expression"] = "dependency",
    *,
    figsize: tuple[float, float] | None = None,
    ax: plt.Axes | None = None,
) -> tuple[plt.Figure, plt.Axes]:
    """Box + strip plot comparing mutant vs WT for a target gene.

    Args:
        df: DataFrame from :func:`analyse_by_mutation`.
        target_gene: Gene symbol (column in *df*).
        mutation_gene: Used for plot title.
        assay: Controls axis label and reference line.
        figsize: Override default figure size.
        ax: Optional pre-existing axes.

    Returns:
        Tuple of (Figure, Axes).
    """
    plt.style.use(str(_STYLE))

    pdf = df.to_pandas()
    order = ["WT", "mutant"]
    palette = {
        "WT": COLOR.GREY.value,
        "mutant": COLOR.PINK.value,
    }

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize or (2.0, 2.5))
    else:
        fig = ax.get_figure()

    sns.boxplot(
        data=pdf,
        x="mutation_status",
        y=target_gene,
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
        y=target_gene,
        order=order,
        color=COLOR.DARKGREY.value,
        size=1.5,
        alpha=0.4,
        jitter=True,
        ax=ax,
    )

    # stats annotation
    wt_vals = pdf.loc[pdf["mutation_status"] == "WT", target_gene]
    mut_vals = pdf.loc[pdf["mutation_status"] == "mutant", target_gene]
    _, p_value = stats.mannwhitneyu(wt_vals, mut_vals, alternative="two-sided")
    n_wt, n_mut = len(wt_vals), len(mut_vals)

    ax.set_title(
        f"{target_gene} — {mutation_gene} mutant vs WT\n"
        f"(n={n_mut} vs {n_wt}, p={p_value:.2e})",
        fontsize=7,
    )

    # reference line
    if assay == "dependency":
        ax.axhline(-1, ls="--", lw=0.5, color=COLOR.DARKGREY.value, zorder=0)
        ax.set_ylabel(f"{target_gene} gene effect")
    else:
        ax.set_ylabel(f"{target_gene} log2(TPM+1)")

    ax.set_xlabel("")
    fig.tight_layout()
    return fig, ax


def mutation_analysis(
    target_gene: str,
    mutation_gene: str,
    assay: Literal["dependency", "expression"] = "dependency",
    gene_type: Literal["suppressor", "oncogene"] = "suppressor",
    min_models: int = 20,
    figsize: tuple[float, float] | None = None,
) -> plt.Figure:
    """Convenience function: analyse and plot in one step."""
    df = analyse_by_mutation(
        target_gene=target_gene,
        mutation_gene=mutation_gene,
        assay=assay,
        gene_type=gene_type,
        min_models=min_models,
    )
    fig, _ = plot_mutation_comparison(
        df=df,
        target_gene=target_gene,
        mutation_gene=mutation_gene,
        assay=assay,
        figsize=figsize,
    )
    return fig
