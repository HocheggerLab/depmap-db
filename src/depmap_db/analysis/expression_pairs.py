"""Helpers for pairwise DepMap expression analyses in notebooks."""

from __future__ import annotations

from dataclasses import dataclass
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


@dataclass(frozen=True)
class HighlightTarget:
    """One explicitly requested model highlight for notebook reporting."""

    label: str
    plot_label: str
    aliases: tuple[str, ...]


HIGHLIGHT_TARGETS: Final[tuple[HighlightTarget, ...]] = (
    HighlightTarget(
        label="HeLa",
        plot_label="HeLa",
        aliases=("ACH-001086", "HeLa", "HELA_CERVIX"),
    ),
    HighlightTarget(
        label="RPE1-ss111",
        plot_label="ss111",
        aliases=("ACH-002466", "RPE1-ss111", "RPE1SS111_ENGINEERED"),
    ),
    HighlightTarget(
        label="RPE1-ss119",
        plot_label="ss119",
        aliases=("ACH-002465", "RPE1-ss119", "RPE1SS119_ENGINEERED"),
    ),
    HighlightTarget(
        label="RPE1-ss77",
        plot_label="ss77",
        aliases=("ACH-002463", "RPE1-ss77", "RPE1SS77_ENGINEERED"),
    ),
    HighlightTarget(
        label="RPE1-ss51",
        plot_label="ss51",
        aliases=("ACH-002467", "RPE1-ss51", "RPE1SS51_ENGINEERED"),
    ),
    HighlightTarget(
        label="RPE1-ss48",
        plot_label="ss48",
        aliases=("ACH-002462", "RPE1-ss48", "RPE1SS48_ENGINEERED"),
    ),
    HighlightTarget(
        label="RPE1-ss6",
        plot_label="ss6",
        aliases=("ACH-002464", "RPE1-ss6", "RPE1SS6_ENGINEERED"),
    ),
)

MODEL_NAME_COLUMNS: Final[tuple[str, ...]] = (
    "cell_line_name",
    "stripped_cell_line_name",
    "ccle_name",
)
MODEL_MATCH_COLUMNS: Final[tuple[str, ...]] = ("model_id", *MODEL_NAME_COLUMNS)


def _normalise_label(value: str) -> str:
    return "".join(char for char in value.upper() if char.isalnum())


def _target_aliases(target: str | HighlightTarget) -> tuple[str, ...]:
    if isinstance(target, HighlightTarget):
        return target.aliases
    return (target,)


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


def exact_model_matches(
    data: pl.DataFrame,
    target: str | HighlightTarget,
) -> pl.DataFrame:
    """Return exact-ish matches after punctuation-insensitive normalisation."""
    normalised_targets = {_normalise_label(alias) for alias in _target_aliases(target)}
    match_expr = pl.lit(False)
    for column in MODEL_MATCH_COLUMNS:
        normalised_column = pl.col(column).fill_null("").map_elements(
            _normalise_label,
            return_dtype=pl.String,
        )
        for normalised_target in normalised_targets:
            match_expr = match_expr | (normalised_column == normalised_target)

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
        for column in MODEL_MATCH_COLUMNS:
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
        for column in MODEL_MATCH_COLUMNS:
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


def build_highlight_table(
    data: pl.DataFrame,
    targets: tuple[HighlightTarget, ...] = HIGHLIGHT_TARGETS,
) -> pl.DataFrame:
    """Combine exact highlight hits into a tidy table for reporting."""
    frames: list[pl.DataFrame] = []
    for target in targets:
        matches = exact_model_matches(data, target)
        if matches.height == 0:
            continue
        frames.append(matches.with_columns(pl.lit(target.label).alias("highlight_target")))

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
