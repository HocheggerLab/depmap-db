"""Shared helpers for the plots sub-package."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
from scipy import stats

from depmap_db import (
    scan_gene_effects_wide,
    scan_gene_expression_wide,
    scan_models,
    scan_mutations,
)

_STYLE = Path(__file__).parent / "hhlab_style01.mplstyle"


def use_style() -> None:
    """Activate the hhlab matplotlib stylesheet."""
    plt.style.use(str(_STYLE))


def resolve_gene_column(
    schema_names: list[str],
    gene: str,
) -> str:
    """Find the exact column name for *gene* in a wide table.

    Matches ``"TP53"`` or the ``"TP53 (7157)"`` Entrez-suffixed form.

    Raises:
        StopIteration: If *gene* is not present in *schema_names*.
    """
    if gene in schema_names:
        return gene
    return next(c for c in schema_names if c.startswith(f"{gene} ("))


def scan_assay(
    assay: Literal["dependency", "expression"],
    *,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Return the LazyFrame for the requested assay."""
    if assay == "dependency":
        return scan_gene_effects_wide(db_path=db_path)
    return scan_gene_expression_wide(db_path=db_path)


def extract_gene_values(
    assay: Literal["dependency", "expression"],
    gene: str,
    *,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Return ``model_id`` + *gene* value column from an assay table."""
    data = scan_assay(assay, db_path=db_path)
    col = resolve_gene_column(data.collect_schema().names(), gene)
    return data.select("model_id", pl.col(col)).rename({col: gene})


def join_models(
    values: pl.LazyFrame,
    *,
    db_path: str | Path | None = None,
) -> pl.LazyFrame:
    """Join *values* with the models table for lineage metadata."""
    models = scan_models(db_path=db_path).select(
        "model_id", "oncotree_lineage"
    )
    return models.join(values, on="model_id")


def get_mutant_model_ids(
    mutation_gene: str,
    gene_type: Literal["suppressor", "oncogene"] = "suppressor",
    *,
    db_path: str | Path | None = None,
) -> pl.Series:
    """Return a Series of model_ids carrying damaging mutations."""
    muts = scan_mutations(db_path=db_path).filter(
        pl.col("hugo_symbol") == mutation_gene
    )
    if gene_type == "suppressor":
        muts = muts.filter(
            pl.col("likely_lof")
            | pl.col("tumor_suppressor_high_impact")
        )
    else:
        muts = muts.filter(
            pl.col("hotspot")
            | pl.col("oncogene_high_impact")
            | pl.col("hess_driver")
        )
    return muts.select("model_id").unique().collect().get_column("model_id")


def label_mutation_status(
    df: pl.DataFrame,
    mutant_ids: pl.Series,
) -> pl.DataFrame:
    """Add a ``mutation_status`` column (``"mutant"`` / ``"WT"``)."""
    id_list = mutant_ids.to_list()
    return df.with_columns(
        pl.when(pl.col("model_id").is_in(id_list))
        .then(pl.lit("mutant"))
        .otherwise(pl.lit("WT"))
        .alias("mutation_status")
    )


def pearson_spearman(
    x: np.ndarray,
    y: np.ndarray,
) -> dict[str, float]:
    """Compute Pearson and Spearman correlations with p-values."""
    pr, pp = stats.pearsonr(x, y)
    sr, sp = stats.spearmanr(x, y)
    return {
        "pearson_r": float(pr),
        "pearson_p": float(pp),
        "spearman_r": float(sr),
        "spearman_p": float(sp),
    }


def annotate_correlation(
    ax: plt.Axes,
    r: float,
    p: float,
    label: str = "r",
) -> None:
    """Add a correlation annotation in the top-left corner of *ax*."""
    ax.text(
        0.03,
        0.97,
        f"{label} = {r:.3f}  (p = {p:.2e})",
        transform=ax.transAxes,
        fontsize=5,
        va="top",
        ha="left",
    )


def make_axes(
    figsize: tuple[float, float] | None = None,
    ax: plt.Axes | None = None,
    default_figsize: tuple[float, float] = (2.5, 2.5),
) -> tuple[plt.Figure, plt.Axes]:
    """Return (Figure, Axes), creating a new figure if *ax* is None."""
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize or default_figsize)
    else:
        fig = ax.get_figure()
    return fig, ax


def assay_ylabel(gene: str, assay: Literal["dependency", "expression"]) -> str:
    """Return a standard y-axis label for the given assay."""
    if assay == "dependency":
        return f"{gene} gene effect"
    return f"{gene} log2(TPM+1)"
