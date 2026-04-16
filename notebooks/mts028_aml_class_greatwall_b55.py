import marimo

__generated_with = "0.20.4"
app = marimo.App(width="wide")


@app.cell
def _():
    import os
    import sys
    from pathlib import Path

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import polars as pl

    PROJECT_ROOT = Path(__file__).resolve().parents[1]
    SRC_ROOT = PROJECT_ROOT / "src"
    if str(SRC_ROOT) not in sys.path:
        sys.path.insert(0, str(SRC_ROOT))

    from depmap_db.polars import prepare_lazy_datasets, prepare_lazy_tables

    plt.style.use("seaborn-v0_8-whitegrid")

    DEFAULT_DB_PATH = PROJECT_ROOT / "data" / "depmap.duckdb"
    DB_PATH = Path(os.environ.get("DEPMAP_DB_PATH", DEFAULT_DB_PATH))
    POLARS_DIR = PROJECT_ROOT / "data" / "polars"

    AML_DISEASE = "Acute Myeloid Leukemia"
    TARGET_GENES = [
        "MASTL",
        "ENSA",
        "ARPP19",
        "PPP2R2A",
        "PPP2R2B",
        "PPP2R2C",
        "PPP2R2D",
    ]
    MODALITY_LABELS = {
        "dependency": "CRISPR gene effect",
        "expression": "RNA expression",
        "proteomics": "Gygi MS proteomics",
    }
    return (
        AML_DISEASE,
        DB_PATH,
        MODALITY_LABELS,
        POLARS_DIR,
        TARGET_GENES,
        mo,
        np,
        pd,
        pl,
        plt,
        prepare_lazy_datasets,
        prepare_lazy_tables,
    )


@app.cell
def _(mo):
    mo.md(r"""
    # AML as a class in the Greatwall/B55 axis

    This notebook asks a narrower, more class-level question than the first AML responder notebook:

    **Does AML stand out relative to other tumour groups in the MASTL/ENSA/ARPP19/PPP2R2 axis?**

    It uses four DepMap surfaces where available:

    - **dependency** (`gene_effects_wide`)
    - **expression** (`gene_expression_wide`)
    - **proteomics** (`proteomics_long`, derived from Gygi MS)
    - **mutation** (`model_gene_mutation_status` and `mutation_events`)

    The notebook stays deliberately descriptive and hypothesis-generating:

    - AML is defined as `oncotree_primary_disease == "Acute Myeloid Leukemia"`.
    - Non-AML comparators are shown both as a pooled **non-AML** group and as **non-AML lineages**.
    - Effects are summarized with counts, medians, and AML-minus-non-AML deltas.
    - Sparse modalities are flagged rather than stretched into strong claims.
    """)
    return


@app.cell
def _(DB_PATH, mo):
    min_lineage_n = mo.ui.slider(
        start=5,
        stop=40,
        step=1,
        value=15,
        label="Minimum lineage size for lineage-level summaries",
    )
    show_top_lineages = mo.ui.slider(
        start=5,
        stop=20,
        step=1,
        value=10,
        label="Number of lineage rows to show around AML in rankings",
    )
    mo.vstack(
        [
            mo.md(f"**Database path:** `{DB_PATH}`"),
            mo.hstack([min_lineage_n, show_top_lineages], justify="start"),
        ]
    )
    return min_lineage_n, show_top_lineages


@app.cell
def _(AML_DISEASE, MODALITY_LABELS, TARGET_GENES, np, pd, plt):
    def add_aml_flags(df: pd.DataFrame) -> pd.DataFrame:
        out = df.copy()
        out["is_aml"] = out["oncotree_primary_disease"].eq(AML_DISEASE)
        out["aml_status"] = np.where(out["is_aml"], "AML", "Non-AML")
        out["analysis_group"] = np.where(
            out["is_aml"],
            "AML",
            out["oncotree_lineage"].fillna("Unknown"),
        )
        return out

    def melt_wide_targets(
        wide_df: pd.DataFrame,
        *,
        modality: str,
        meta_cols: list[str],
    ) -> pd.DataFrame:
        feature_cols = [gene for gene in TARGET_GENES if gene in wide_df.columns]
        if not feature_cols:
            return pd.DataFrame(
                columns=meta_cols + ["feature", "value", "modality"]
            )
        long_df = wide_df.loc[:, meta_cols + feature_cols].melt(
            id_vars=meta_cols,
            value_vars=feature_cols,
            var_name="feature",
            value_name="value",
        )
        long_df["modality"] = modality
        return long_df.dropna(subset=["value"]).reset_index(drop=True)

    def summarize_binary_comparison(long_df: pd.DataFrame) -> pd.DataFrame:
        if long_df.empty:
            return pd.DataFrame()
        grouped = (
            long_df.groupby(["modality", "feature", "aml_status"], as_index=False)
            .agg(
                n_models=("model_id", "nunique"),
                mean_value=("value", "mean"),
                median_value=("value", "median"),
                q1_value=("value", lambda x: x.quantile(0.25)),
                q3_value=("value", lambda x: x.quantile(0.75)),
            )
        )
        wide = grouped.pivot(
            index=["modality", "feature"],
            columns="aml_status",
            values=["n_models", "mean_value", "median_value", "q1_value", "q3_value"],
        )
        wide.columns = [f"{left}_{right.lower().replace('-', '_')}" for left, right in wide.columns]
        out = wide.reset_index()
        for column in [
            "n_models_aml",
            "n_models_non_aml",
            "mean_value_aml",
            "mean_value_non_aml",
            "median_value_aml",
            "median_value_non_aml",
        ]:
            if column not in out.columns:
                out[column] = np.nan
        out["aml_minus_non_aml_median"] = (
            out["median_value_aml"] - out["median_value_non_aml"]
        )
        out["aml_minus_non_aml_mean"] = (
            out["mean_value_aml"] - out["mean_value_non_aml"]
        )
        out["modality_label"] = out["modality"].map(MODALITY_LABELS)
        return out.sort_values(["modality", "feature"]).reset_index(drop=True)

    def summarize_lineage_medians(
        long_df: pd.DataFrame,
        *,
        min_n: int,
    ) -> pd.DataFrame:
        if long_df.empty:
            return pd.DataFrame()
        grouped = (
            long_df.groupby(["modality", "feature", "analysis_group"], as_index=False)
            .agg(
                n_models=("model_id", "nunique"),
                median_value=("value", "median"),
                mean_value=("value", "mean"),
            )
        )
        grouped = grouped.loc[grouped["n_models"] >= min_n].copy()
        if grouped.empty:
            return grouped

        rank_frames: list[pd.DataFrame] = []
        for (modality, feature), sub in grouped.groupby(["modality", "feature"]):
            ranking = sub.sort_values("median_value", ascending=False).reset_index(drop=True)
            ranking["median_rank_desc"] = np.arange(1, len(ranking) + 1)
            ranking["rank_out_of"] = len(ranking)
            ranking["modality"] = modality
            ranking["feature"] = feature
            rank_frames.append(ranking)
        return pd.concat(rank_frames, ignore_index=True)

    def summarize_feature_coverage(long_df: pd.DataFrame) -> pd.DataFrame:
        if long_df.empty:
            return pd.DataFrame()
        return (
            long_df.groupby(["modality", "feature", "aml_status"], as_index=False)
            .agg(n_models=("model_id", "nunique"))
            .pivot(index=["modality", "feature"], columns="aml_status", values="n_models")
            .reset_index()
            .rename(columns={"AML": "aml_models", "Non-AML": "non_aml_models"})
            .fillna(0)
            .sort_values(["modality", "feature"])
            .reset_index(drop=True)
        )

    def summarize_mutation_axis(
        mutation_status_df: pd.DataFrame,
        mutation_events_df: pd.DataFrame,
    ) -> pd.DataFrame:
        status_summary = (
            mutation_status_df.groupby(["aml_status", "gene_symbol"], as_index=False)
            .agg(
                mutated_models=("model_id", "nunique"),
                mutation_status_rows=("model_id", "size"),
                hotspot_models=("has_hotspot", "sum"),
                lof_models=("has_likely_lof", "sum"),
                driver_models=("has_driver", "sum"),
            )
        ) if not mutation_status_df.empty else pd.DataFrame(
            columns=["aml_status", "gene_symbol", "mutated_models", "mutation_status_rows", "hotspot_models", "lof_models", "driver_models"]
        )

        event_summary = (
            mutation_events_df.groupby(["aml_status", "gene_symbol"], as_index=False)
            .agg(
                mutation_events=("mutation_id", "nunique"),
                event_models=("model_id", "nunique"),
                hotspot_events=("hotspot", "sum"),
                lof_events=("likely_lof", "sum"),
                driver_events=("hess_driver", "sum"),
            )
        ) if not mutation_events_df.empty else pd.DataFrame(
            columns=["aml_status", "gene_symbol", "mutation_events", "event_models", "hotspot_events", "lof_events", "driver_events"]
        )

        template = pd.MultiIndex.from_product(
            [["AML", "Non-AML"], TARGET_GENES],
            names=["aml_status", "gene_symbol"],
        ).to_frame(index=False)
        out = template.merge(status_summary, on=["aml_status", "gene_symbol"], how="left")
        out = out.merge(event_summary, on=["aml_status", "gene_symbol"], how="left")
        numeric_cols = [
            "mutated_models",
            "mutation_status_rows",
            "hotspot_models",
            "lof_models",
            "driver_models",
            "mutation_events",
            "event_models",
            "hotspot_events",
            "lof_events",
            "driver_events",
        ]
        out[numeric_cols] = out[numeric_cols].fillna(0).astype(int)
        return out

    def plot_binary_boxgrid(long_df: pd.DataFrame, modality: str):
        plot_df = long_df.loc[long_df["modality"].eq(modality)].copy()
        features = [gene for gene in TARGET_GENES if gene in plot_df["feature"].unique()]
        if not features:
            return None

        fig, axes = plt.subplots(
            nrows=len(features),
            ncols=1,
            figsize=(8.5, max(2.6 * len(features), 4.5)),
            sharex=False,
        )
        if len(features) == 1:
            axes = [axes]

        color_map = {"AML": "#b2182b", "Non-AML": "#4c78a8"}
        for ax, feature in zip(axes, features, strict=False):
            sub = plot_df.loc[plot_df["feature"].eq(feature)].copy()
            groups = [
                sub.loc[sub["aml_status"].eq("AML"), "value"].dropna().to_numpy(),
                sub.loc[sub["aml_status"].eq("Non-AML"), "value"].dropna().to_numpy(),
            ]
            bp = ax.boxplot(groups, labels=["AML", "Non-AML"], patch_artist=True)
            for patch, label in zip(bp["boxes"], ["AML", "Non-AML"], strict=False):
                patch.set_facecolor(color_map[label])
                patch.set_alpha(0.7)
            ax.set_title(feature, loc="left", fontsize=11)
            ax.set_ylabel(MODALITY_LABELS[modality])
        fig.suptitle(f"AML vs pooled non-AML: {MODALITY_LABELS[modality]}", y=1.0)
        fig.tight_layout()
        return fig

    def plot_lineage_ranking(
        lineage_summary_df: pd.DataFrame,
        *,
        modality: str,
        feature: str,
        top_n: int,
    ):
        sub = lineage_summary_df.loc[
            lineage_summary_df["modality"].eq(modality)
            & lineage_summary_df["feature"].eq(feature)
        ].copy()
        if sub.empty:
            return None
        sub = sub.sort_values("median_value", ascending=False).reset_index(drop=True)
        aml_row = sub.loc[sub["analysis_group"].eq("AML")]
        if aml_row.empty:
            plot_df = sub.head(top_n)
        else:
            aml_rank = int(aml_row["median_rank_desc"].iloc[0]) - 1
            low = max(0, aml_rank - top_n // 2)
            high = min(len(sub), low + top_n)
            plot_df = sub.iloc[low:high].copy()
        colors = ["#b2182b" if group == "AML" else "#808080" for group in plot_df["analysis_group"]]
        fig, ax = plt.subplots(figsize=(9, max(4.5, 0.45 * len(plot_df))))
        ax.barh(plot_df["analysis_group"], plot_df["median_value"], color=colors, alpha=0.9)
        for i, row in enumerate(plot_df.itertuples(index=False)):
            ax.text(
                float(row.median_value),
                i,
                f"  n={row.n_models}, rank={row.median_rank_desc}/{row.rank_out_of}",
                va="center",
                fontsize=8.5,
            )
        ax.invert_yaxis()
        ax.set_xlabel(f"Median {MODALITY_LABELS[modality]}")
        ax.set_ylabel("")
        ax.set_title(f"{feature}: AML positioned among lineage medians")
        fig.tight_layout()
        return fig

    def plot_delta_heatmap(binary_summary_df: pd.DataFrame):
        if binary_summary_df.empty:
            return None
        heat_df = (
            binary_summary_df.pivot(
                index="modality",
                columns="feature",
                values="aml_minus_non_aml_median",
            )
            .reindex(index=[m for m in ["dependency", "expression", "proteomics"] if m in binary_summary_df["modality"].unique()])
        )
        if heat_df.empty:
            return None
        matrix = heat_df.to_numpy(dtype=float)
        fig, ax = plt.subplots(figsize=(1.25 * len(heat_df.columns) + 2, 3.8))
        im = ax.imshow(matrix, cmap="coolwarm", aspect="auto")
        ax.set_xticks(np.arange(len(heat_df.columns)), heat_df.columns)
        ax.set_yticks(np.arange(len(heat_df.index)), [MODALITY_LABELS[idx] for idx in heat_df.index])
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                value = matrix[i, j]
                label = "NA" if pd.isna(value) else f"{value:.2f}"
                ax.text(j, i, label, ha="center", va="center", fontsize=9)
        ax.set_title("AML minus pooled non-AML median by modality and feature")
        fig.colorbar(im, ax=ax, shrink=0.8)
        fig.tight_layout()
        return fig

    return (
        add_aml_flags,
        melt_wide_targets,
        plot_binary_boxgrid,
        plot_delta_heatmap,
        plot_lineage_ranking,
        summarize_binary_comparison,
        summarize_feature_coverage,
        summarize_lineage_medians,
        summarize_mutation_axis,
    )


@app.cell
def _(
    DB_PATH,
    POLARS_DIR,
    TARGET_GENES,
    add_aml_flags,
    pd,
    pl,
    prepare_lazy_datasets,
    prepare_lazy_tables,
):
    model_cols = [
        "model_id",
        "cell_line_name",
        "stripped_cell_line_name",
        "ccle_name",
        "oncotree_lineage",
        "oncotree_primary_disease",
        "oncotree_subtype",
    ]

    tables = prepare_lazy_tables(
        output_dir=POLARS_DIR,
        db_path=DB_PATH,
        tables=[
            "models",
            "gene_effects_wide",
            "gene_expression_wide",
            "model_gene_mutation_status",
        ],
    )

    expression_keep = [gene for gene in TARGET_GENES if gene in tables["gene_expression_wide"].collect_schema().names()]
    dependency_keep = [gene for gene in TARGET_GENES if gene in tables["gene_effects_wide"].collect_schema().names()]

    expression_df = (
        tables["gene_expression_wide"]
        .select(["model_id", *expression_keep])
        .join(tables["models"].select(model_cols), on="model_id", how="left")
        .collect()
        .to_pandas()
    )
    expression_df = add_aml_flags(expression_df)

    dependency_df = (
        tables["gene_effects_wide"]
        .select(["model_id", *dependency_keep])
        .join(tables["models"].select(model_cols), on="model_id", how="left")
        .collect()
        .to_pandas()
    )
    dependency_df = add_aml_flags(dependency_df)

    mutation_status_df = (
        tables["model_gene_mutation_status"]
        .filter(pl.col("gene_symbol").is_in(TARGET_GENES))
        .join(tables["models"].select(model_cols), on="model_id", how="left")
        .collect()
        .to_pandas()
    )
    mutation_status_df = add_aml_flags(mutation_status_df)

    mutation_events_lazy = prepare_lazy_datasets(
        output_dir=POLARS_DIR,
        db_path=DB_PATH,
        datasets=["mutation_events"],
    )["mutation_events"]
    mutation_events_df = (
        mutation_events_lazy
        .filter(pl.col("gene_symbol").is_in(TARGET_GENES))
        .collect()
        .to_pandas()
    )
    mutation_events_df = add_aml_flags(mutation_events_df)

    proteomics_issue = None
    proteomics_df = pd.DataFrame()
    try:
        proteomics_lazy = prepare_lazy_datasets(
            output_dir=POLARS_DIR,
            db_path=DB_PATH,
            datasets=["proteomics_long"],
        )["proteomics_long"]
        proteomics_df = (
            proteomics_lazy
            .filter(pl.col("gene_symbol").is_in(TARGET_GENES))
            .select(model_cols + ["gene_symbol", "protein_accession", "abundance"])
            .collect()
            .to_pandas()
        )
        proteomics_df = proteomics_df.rename(columns={"gene_symbol": "feature", "abundance": "value"})
        proteomics_df = add_aml_flags(proteomics_df)
    except (OSError, RuntimeError, ValueError) as exc:  # pragma: no cover - notebook fallback path
        proteomics_issue = f"{type(exc).__name__}: {exc}"

    model_summary = (
        add_aml_flags(
            tables["models"].select(model_cols).collect().to_pandas()
        )
        .groupby("aml_status", as_index=False)
        .agg(n_models=("model_id", "nunique"))
    )
    return (
        dependency_df,
        dependency_keep,
        expression_df,
        expression_keep,
        model_cols,
        model_summary,
        mutation_events_df,
        mutation_status_df,
        proteomics_df,
        proteomics_issue,
    )


@app.cell
def _(
    DB_PATH,
    dependency_keep,
    expression_keep,
    mo,
    model_summary,
    proteomics_df,
    proteomics_issue,
):
    proteomics_genes = sorted(proteomics_df["feature"].dropna().unique().tolist()) if not proteomics_df.empty else []
    mo.vstack(
        [
            mo.md("## Data coverage and schema-facing feature availability"),
            mo.md(
                f"- DB found at `{DB_PATH}`\n"
                f"- Expression columns present for target genes: **{', '.join(expression_keep)}**\n"
                f"- Dependency columns present for target genes: **{', '.join(dependency_keep)}**\n"
                f"- Proteomics genes present in `proteomics_long`: **{', '.join(proteomics_genes) if proteomics_genes else 'none'}**"
            ),
            model_summary,
            mo.md(
                "**Immediate caveat:** dependency does not expose `PPP2R2D` in the current wide table, and proteomics coverage is more limited than RNA."
                if proteomics_issue is None
                else f"**Proteomics caveat:** `proteomics_long` could not be loaded here ({proteomics_issue})."
            ),
        ]
    )
    return


@app.cell
def _(
    dependency_df,
    expression_df,
    melt_wide_targets,
    model_cols,
    pd,
    proteomics_df,
):
    expression_long = melt_wide_targets(
        expression_df,
        modality="expression",
        meta_cols=model_cols + ["is_aml", "aml_status", "analysis_group"],
    )
    dependency_long = melt_wide_targets(
        dependency_df,
        modality="dependency",
        meta_cols=model_cols + ["is_aml", "aml_status", "analysis_group"],
    )
    modality_long = pd.concat(
        [dependency_long, expression_long, proteomics_df],
        ignore_index=True,
        sort=False,
    )
    return (modality_long,)


@app.cell
def _(
    min_lineage_n,
    modality_long,
    mutation_events_df,
    mutation_status_df,
    summarize_binary_comparison,
    summarize_feature_coverage,
    summarize_lineage_medians,
    summarize_mutation_axis,
):
    binary_summary = summarize_binary_comparison(modality_long)
    coverage_summary = summarize_feature_coverage(modality_long)
    lineage_summary = summarize_lineage_medians(
        modality_long,
        min_n=min_lineage_n.value,
    )
    mutation_axis_summary = summarize_mutation_axis(
        mutation_status_df,
        mutation_events_df,
    )
    return (
        binary_summary,
        coverage_summary,
        lineage_summary,
        mutation_axis_summary,
    )


@app.cell
def _(binary_summary, coverage_summary, mo, plot_delta_heatmap):
    delta_plot = plot_delta_heatmap(binary_summary)
    mo.md("## AML vs pooled non-AML across the target axis")
    mo.vstack(
        [
            mo.md(
                "This is the cleanest first pass: for each modality and axis feature, compare AML to pooled non-AML using observed models only."
            ),
            coverage_summary,
            binary_summary[
                [
                    "modality_label",
                    "feature",
                    "n_models_aml",
                    "n_models_non_aml",
                    "median_value_aml",
                    "median_value_non_aml",
                    "aml_minus_non_aml_median",
                    "mean_value_aml",
                    "mean_value_non_aml",
                ]
            ],
            delta_plot,
        ]
    )
    return


@app.cell
def _(MODALITY_LABELS, mo, modality_long):
    available_modalities = list(dict.fromkeys(modality_long["modality"].tolist()))
    modality_picker = mo.ui.dropdown(
        options={modality: MODALITY_LABELS[modality] for modality in available_modalities},
        value=available_modalities[0],
        label="Modality for AML vs non-AML distribution plots",
    )
    mo.vstack([modality_picker])
    return available_modalities, modality_picker


@app.cell
def _(mo, modality_long, modality_picker, plot_binary_boxgrid):
    binary_plot = plot_binary_boxgrid(modality_long, modality_picker.value)
    mo.vstack(
        [
            mo.md("## Distribution view: AML vs pooled non-AML"),
            mo.md(
                "Use this to see whether AML looks shifted relative to the rest of DepMap for a whole modality, while keeping the actual spread visible."
            ),
            binary_plot,
        ]
    )
    return


@app.cell
def _(TARGET_GENES, available_modalities, mo):
    feature_picker = mo.ui.dropdown(
        options=TARGET_GENES,
        value="MASTL",
        label="Feature for lineage ranking",
    )
    lineage_modality_picker = mo.ui.dropdown(
        options=available_modalities,
        value=available_modalities[0],
        label="Modality for lineage ranking",
    )
    mo.hstack([lineage_modality_picker, feature_picker], justify="start")
    return feature_picker, lineage_modality_picker


@app.cell
def _(
    feature_picker,
    lineage_modality_picker,
    lineage_summary,
    mo,
    plot_lineage_ranking,
    show_top_lineages,
):
    selected_lineage_summary = lineage_summary.loc[
        lineage_summary["modality"].eq(lineage_modality_picker.value)
        & lineage_summary["feature"].eq(feature_picker.value)
    ].copy()
    lineage_plot = plot_lineage_ranking(
        lineage_summary,
        modality=lineage_modality_picker.value,
        feature=feature_picker.value,
        top_n=show_top_lineages.value,
    )
    mo.vstack(
        [
            mo.md("## AML compared with other lineages"),
            mo.md(
                "AML is treated as its own analysis group; all non-AML models are grouped by lineage. The rank column shows where AML sits among lineage medians for the selected feature and modality."
            ),
            selected_lineage_summary,
            lineage_plot,
        ]
    )
    return


@app.cell
def _(mo, mutation_axis_summary):
    aml_mut_count = int(
        mutation_axis_summary.loc[
            mutation_axis_summary["aml_status"].eq("AML"), "mutated_models"
        ].sum()
    )
    mutation_comment = (
        "AML carries essentially no mutations in the focal axis genes in this DepMap slice, so mutation is not an obvious within-axis explanation for any AML class signal here."
        if aml_mut_count == 0
        else "A small number of AML models do carry axis-gene mutations, but counts are still sparse enough that they should be treated cautiously."
    )
    mo.vstack(
        [
            mo.md("## Mutation summary for axis genes"),
            mo.md(mutation_comment),
            mutation_axis_summary,
        ]
    )
    return


@app.cell
def _(binary_summary, lineage_summary, mo, mutation_axis_summary, pd):
    takeaway_rows = []
    for modality in ["dependency", "expression", "proteomics"]:
        sub = binary_summary.loc[binary_summary["modality"].eq(modality)].copy()
        if sub.empty:
            continue
        top_pos = sub.sort_values("aml_minus_non_aml_median", ascending=False).head(2)
        top_neg = sub.sort_values("aml_minus_non_aml_median", ascending=True).head(2)
        for direction, frame in [("AML higher", top_pos), ("AML lower", top_neg)]:
            for row in frame.itertuples(index=False):
                takeaway_rows.append(
                    {
                        "modality": modality,
                        "direction": direction,
                        "feature": row.feature,
                        "aml_minus_non_aml_median": row.aml_minus_non_aml_median,
                        "aml_n": row.n_models_aml,
                        "non_aml_n": row.n_models_non_aml,
                    }
                )
    takeaway_df = pd.DataFrame(takeaway_rows)

    aml_ranks = (
        lineage_summary.loc[lineage_summary["analysis_group"].eq("AML")]
        .loc[:, ["modality", "feature", "median_rank_desc", "rank_out_of", "median_value", "n_models"]]
        .sort_values(["modality", "median_rank_desc", "feature"])
        .reset_index(drop=True)
    )

    aml_axis_mutated_models = int(
        mutation_axis_summary.loc[
            mutation_axis_summary["aml_status"].eq("AML"), "mutated_models"
        ].sum()
    )

    mo.vstack(
        [
            mo.md("## Working interpretation"),
            mo.md(
                f"""
                A few patterns are worth carrying forward, but none should be overclaimed from this notebook alone.

                - **Expression and proteomics can point in different directions from dependency.** That is expected in this axis and is part of the biological question, not a failure state.
                - **Coverage is uneven across modalities.** RNA is broadest, dependency is narrower, and proteomics is the sparsest.
                - **Axis-gene mutation is not a strong AML-level driver in this slice.** Total AML mutated-model count across the focal genes here: **{aml_axis_mutated_models}**.
                - **AML rank is feature-specific rather than universally extreme.** Check the lineage-rank table before telling a single global story about “the axis”.
                """
            ),
            mo.md("### Largest AML-minus-non-AML median shifts by modality"),
            takeaway_df,
            mo.md("### AML ranks among lineage medians"),
            aml_ranks,
        ]
    )
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Practical caveats

    1. **This is a class-level descriptive notebook, not a responder model.**
    2. **`PPP2R2D` is absent from the current dependency wide table**, so that part of the axis is incomplete in CRISPR space.
    3. **Proteomics is sparse and partially mapped at the isoform/accession level.** Absence of a protein here can mean assay coverage/mapping limits rather than true biological absence.
    4. **Mutation scarcity is itself informative** for these genes in AML, but it also means mutation-based comparison is low-power.
    5. If AML looks unusual for one feature/modality, the clean next step is to validate that specific contrast elsewhere rather than promoting a whole-axis claim immediately.
    """)
    return


if __name__ == "__main__":
    app.run()
