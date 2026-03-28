"""Helpers for pairwise DepMap expression analyses in notebooks."""

from __future__ import annotations

from difflib import get_close_matches
from pathlib import Path
from typing import Final

import polars as pl

from depmap_db.polars import prepare_lazy_tables

GENE_ALIASES: Final[dict[str, str]] = {
    "CDK1": "CDK1",
    "Greatwall": "MASTL",
    "MASTL": "MASTL",
    "PPP2R2A": "PPP2R2A",
    "B55alpha": "PPP2R2A",
}

HIGHLIGHT_TARGETS: Final[tuple[str, ...]] = ("HELA", "RPE-1")
MODEL_NAME_COLUMNS: Final[tuple[str, ...]] = (
    "cell_line_name",
    "stripped_cell_line_name",
    "ccle_name",
)


def _normalise_label(value: str) -> str:
    return "".join(char for char in value.upper() if char.isalnum())


def build_expression_dataset(
    *,
    db_path: str | Path | None = None,
    output_dir: str | Path | None = None,
    overwrite: bool = False,
) -> pl.DataFrame:
    """Load the model metadata and selected expression columns via Polars."""
    genes = ["CDK1", "MASTL", "PPP2R2A"]
    tables = prepare_lazy_tables(
        tables=["models", "gene_expression_wide"],
        db_path=db_path,
        output_dir=output_dir,
        overwrite=overwrite,
    )

    models = tables["models"].select(
        [
            "model_id",
            "cell_line_name",
            "stripped_cell_line_name",
            "ccle_name",
            "depmap_model_type",
            "oncotree_lineage",
        ]
    )
    expression = tables["gene_expression_wide"].select(["model_id", *genes])

    return (
        models.join(expression, on="model_id", how="inner")
        .sort("cell_line_name")
        .collect()
    )


def pairwise_frame(data: pl.DataFrame, x_gene: str, y_gene: str) -> pl.DataFrame:
    """Prepare one pairwise comparison with TCGA-style colour metric."""
    return (
        data.select(
            [
                "model_id",
                "cell_line_name",
                "stripped_cell_line_name",
                "ccle_name",
                "depmap_model_type",
                "oncotree_lineage",
                pl.col(x_gene).alias("x"),
                pl.col(y_gene).alias("y"),
            ]
        )
        .filter(pl.col("x").is_not_null() & pl.col("y").is_not_null())
        .with_columns(
            pl.when((pl.col("x") + pl.col("y")) == 0)
            .then(None)
            .otherwise((pl.col("x") - pl.col("y")) / (pl.col("x") + pl.col("y")))
            .alias("diff")
        )
    )


def pearson_r(data: pl.DataFrame, x_col: str = "x", y_col: str = "y") -> float:
    """Compute Pearson's r for an eager Polars DataFrame."""
    return float(data.select(pl.corr(x_col, y_col)).item())


def exact_model_matches(data: pl.DataFrame, target: str) -> pl.DataFrame:
    """Return exact-ish matches after punctuation-insensitive normalisation."""
    normalised_target = _normalise_label(target)
    match_expr = pl.lit(False)
    for column in MODEL_NAME_COLUMNS:
        match_expr = match_expr | (
            pl.col(column)
            .fill_null("")
            .map_elements(_normalise_label, return_dtype=pl.String)
            == normalised_target
        )

    return (
        data.filter(match_expr)
        .select(
            [
                "model_id",
                "cell_line_name",
                "stripped_cell_line_name",
                "ccle_name",
                "depmap_model_type",
                "oncotree_lineage",
                "CDK1",
                "MASTL",
                "PPP2R2A",
            ]
        )
        .unique()
        .sort("cell_line_name")
    )


def nearest_model_candidates(
    data: pl.DataFrame,
    target: str,
    *,
    limit: int = 5,
) -> pl.DataFrame:
    """Return nearest available model-name matches for a missing target."""
    records = list(
        data.select(
            [
                "model_id",
                "cell_line_name",
                "stripped_cell_line_name",
                "ccle_name",
                "depmap_model_type",
                "oncotree_lineage",
                "CDK1",
                "MASTL",
                "PPP2R2A",
            ]
        ).iter_rows(named=True)
    )
    target_token = _normalise_label(target)

    token_matches: list[dict[str, object]] = []
    for row in records:
        for column in MODEL_NAME_COLUMNS:
            raw_name = row[column]
            if raw_name is None:
                continue
            candidate_token = _normalise_label(str(raw_name))
            if target_token in candidate_token or candidate_token in target_token:
                token_matches.append(row)
                break

    if token_matches:
        return pl.DataFrame(token_matches).unique().sort("cell_line_name").head(limit)

    names: dict[str, dict[str, object]] = {}
    for row in records:
        for column in MODEL_NAME_COLUMNS:
            raw_name = row[column]
            if raw_name is None:
                continue
            names[str(raw_name)] = row

    close_names = get_close_matches(target, list(names), n=limit, cutoff=0.4)
    rows = [names[name] for name in close_names]
    if not rows:
        return pl.DataFrame(
            schema={
                "model_id": pl.String,
                "cell_line_name": pl.String,
                "stripped_cell_line_name": pl.String,
                "ccle_name": pl.String,
                "depmap_model_type": pl.String,
                "oncotree_lineage": pl.String,
                "CDK1": pl.Float64,
                "MASTL": pl.Float64,
                "PPP2R2A": pl.Float64,
            }
        )
    return pl.DataFrame(rows).unique().sort("cell_line_name")


def build_highlight_table(data: pl.DataFrame) -> pl.DataFrame:
    """Combine exact highlight hits into a tidy table for reporting."""
    frames: list[pl.DataFrame] = []
    for target in HIGHLIGHT_TARGETS:
        matches = exact_model_matches(data, target)
        if matches.height == 0:
            continue
        frames.append(matches.with_columns(pl.lit(target).alias("highlight_target")))

    if not frames:
        return pl.DataFrame(
            schema={
                "highlight_target": pl.String,
                "model_id": pl.String,
                "cell_line_name": pl.String,
                "stripped_cell_line_name": pl.String,
                "ccle_name": pl.String,
                "depmap_model_type": pl.String,
                "oncotree_lineage": pl.String,
                "CDK1": pl.Float64,
                "MASTL": pl.Float64,
                "PPP2R2A": pl.Float64,
            }
        )

    return pl.concat(frames).select(
        [
            "highlight_target",
            "model_id",
            "cell_line_name",
            "stripped_cell_line_name",
            "ccle_name",
            "depmap_model_type",
            "oncotree_lineage",
            "CDK1",
            "MASTL",
            "PPP2R2A",
        ]
    )
