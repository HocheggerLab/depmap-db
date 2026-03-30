import marimo

__generated_with = "0.21.1"
app = marimo.App(width="wide")


@app.cell
def _():
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

    from depmap_db.polars import prepare_lazy_tables

    plt.style.use("seaborn-v0_8-whitegrid")
    return PROJECT_ROOT, mo, np, pd, pl, plt, prepare_lazy_tables


@app.cell
def _(mo):
    mo.md(r"""
    # Phase 1 AML analysis for the C-604 MTS028 screen

    This notebook is a first-pass look at the internal **MTS028 c-604** screen.

    It focuses on only three questions:

    1. **Is AML more sensitive than other lineages?**
    2. **Is AML internally heterogeneous?**
    3. **Are AML dose-response curves shaped differently from non-AML curves?**

    Notes:
    - The fitted screen summaries come from `drug_response_secondary`.
    - The dose-level collapsed measurements come from `drug_response_secondary_dose`.
    - **AML** is defined here as `oncotree_primary_disease == "Acute Myeloid Leukemia"`.
    - Lower fitted **AUC** means greater sensitivity.
    - For dose-level plots, growth fraction is approximated as `2 ** median_l2fc`, so values near `1` indicate little effect and lower values indicate stronger inhibition.
    """)
    return


@app.cell
def _(PROJECT_ROOT, mo):
    SCREEN_ID = "MTS028_BIAS_CORRECTED"
    AML_DISEASE = "Acute Myeloid Leukemia"
    DB_PATH = PROJECT_ROOT / "data" / "depmap.duckdb"
    POLARS_DIR = PROJECT_ROOT / "data" / "polars"

    min_lineage_n = mo.ui.slider(
        start=5,
        stop=30,
        step=1,
        value=8,
        label="Minimum group size for lineage comparison",
    )
    mo.hstack([min_lineage_n], justify="start")
    return AML_DISEASE, DB_PATH, POLARS_DIR, SCREEN_ID, min_lineage_n


@app.cell
def _(AML_DISEASE, np, pd, plt):
    def add_analysis_labels(df: pd.DataFrame) -> pd.DataFrame:
        out = df.copy()
        out["is_aml"] = out["oncotree_primary_disease"].eq(AML_DISEASE)
        out["analysis_group"] = np.where(
            out["is_aml"],
            "AML",
            out["oncotree_lineage"].fillna("Unknown"),
        )
        out["aml_status"] = np.where(out["is_aml"], "AML", "Non-AML")
        return out

    def summarize_lineages(df: pd.DataFrame, min_n: int = 8) -> pd.DataFrame:
        summary = (
            df.dropna(subset=["auc"])
            .groupby("analysis_group", as_index=False)
            .agg(
                n=("model_id", "size"),
                median_auc=("auc", "median"),
                mean_auc=("auc", "mean"),
                q1_auc=("auc", lambda x: x.quantile(0.25)),
                q3_auc=("auc", lambda x: x.quantile(0.75)),
                median_ic50=("ic50", "median"),
            )
            .query("n >= @min_n")
            .sort_values(["median_auc", "n"], ascending=[True, False])
            .reset_index(drop=True)
        )
        return summary

    def pick_quantile_example(
        df: pd.DataFrame, *, is_aml: bool, quantile: float, label: str
    ) -> pd.Series:
        subset = df.loc[df["is_aml"].eq(is_aml) & df["auc"].notna()].copy()
        target = subset["auc"].quantile(quantile)
        row = subset.iloc[(subset["auc"] - target).abs().argsort()].iloc[0].copy()
        row["example_label"] = label
        return row

    def select_representative_models(df: pd.DataFrame) -> pd.DataFrame:
        rows = [
            pick_quantile_example(
                df, is_aml=True, quantile=0.25, label="AML — sensitive"
            ),
            pick_quantile_example(
                df, is_aml=True, quantile=0.50, label="AML — typical"
            ),
            pick_quantile_example(
                df, is_aml=False, quantile=0.25, label="Non-AML — sensitive"
            ),
            pick_quantile_example(
                df, is_aml=False, quantile=0.50, label="Non-AML — typical"
            ),
        ]
        cols = [
            "example_label",
            "model_id",
            "cell_line_name",
            "analysis_group",
            "oncotree_primary_disease",
            "auc",
            "ic50",
        ]
        return pd.DataFrame(rows)[cols]

    def build_curve_features(dose_df: pd.DataFrame) -> pd.DataFrame:
        feature_rows = []
        for response_id, group in dose_df.groupby("response_id", sort=False):
            group = group.sort_values("dose_um").copy()
            log_dose = np.log10(group["dose_um"].to_numpy())
            growth = group["growth_fraction"].to_numpy()
            diffs = np.diff(growth)
            log_steps = np.diff(log_dose)
            slopes = np.divide(
                diffs,
                log_steps,
                out=np.full_like(diffs, np.nan),
                where=log_steps != 0,
            )
            feature_rows.append(
                {
                    "response_id": response_id,
                    "model_id": group["model_id"].iloc[0],
                    "cell_line_name": group["cell_line_name"].iloc[0],
                    "oncotree_lineage": group["oncotree_lineage"].iloc[0],
                    "oncotree_primary_disease": group[
                        "oncotree_primary_disease"
                    ].iloc[0],
                    "is_aml": bool(group["is_aml"].iloc[0]),
                    "aml_status": group["aml_status"].iloc[0],
                    "low_dose_growth": float(growth[0]),
                    "high_dose_growth": float(growth[-1]),
                    "best_growth": float(np.nanmin(growth)),
                    "dynamic_range_growth": float(growth[0] - np.nanmin(growth)),
                    "delta_end_start_growth": float(growth[-1] - growth[0]),
                    "monotonic_drop_fraction": float(np.mean(diffs <= 0)),
                    "steepest_drop_per_log10_uM": float(np.nanmin(slopes)),
                    "logdose_auc_growth": float(
                        np.trapezoid(growth, x=log_dose)
                        / (log_dose[-1] - log_dose[0])
                    ),
                }
            )
        return pd.DataFrame(feature_rows)

    def plot_lineage_comparison(summary: pd.DataFrame):
        fig, ax = plt.subplots(figsize=(10, max(5, 0.38 * len(summary))))
        plot_df = summary.sort_values("median_auc", ascending=True).copy()
        colors = ["#b2182b" if g == "AML" else "#4c78a8" for g in plot_df["analysis_group"]]
        ax.barh(plot_df["analysis_group"], plot_df["median_auc"], color=colors, alpha=0.9)
        ax.errorbar(
            plot_df["median_auc"],
            plot_df["analysis_group"],
            xerr=[
                plot_df["median_auc"] - plot_df["q1_auc"],
                plot_df["q3_auc"] - plot_df["median_auc"],
            ],
            fmt="none",
            ecolor="black",
            elinewidth=1,
            capsize=2,
        )
        for i, row in enumerate(plot_df.itertuples(index=False)):
            ax.text(row.median_auc + 0.005, i, f"n={row.n}", va="center", fontsize=9)
        ax.set_xlabel("Median fitted AUC (lower = more sensitive)")
        ax.set_ylabel("")
        ax.set_title("Lineage comparison with AML separated from other Myeloid models")
        ax.set_xlim(0.5, 1.03)
        return fig

    def plot_aml_waterfall(aml_df: pd.DataFrame):
        plot_df = aml_df.dropna(subset=["auc"]).sort_values("auc").reset_index(drop=True)
        fig, ax = plt.subplots(figsize=(11, 4.5))
        ax.bar(np.arange(len(plot_df)), plot_df["auc"], color="#b2182b", alpha=0.9)
        ax.axhline(plot_df["auc"].median(), color="black", linestyle="--", linewidth=1)
        ax.set_ylabel("Fitted AUC (lower = more sensitive)")
        ax.set_xlabel("AML models, sorted by AUC")
        ax.set_title("AML waterfall shows broad heterogeneity in C-604 response")
        ax.set_xticks([])
        return fig

    def plot_group_curves(group_curve_summary: pd.DataFrame):
        fig, ax = plt.subplots(figsize=(8, 5))
        palette = {"AML": "#b2182b", "Non-AML": "#4c78a8"}
        for label, sub in group_curve_summary.groupby("aml_status"):
            sub = sub.sort_values("dose_um")
            ax.plot(
                sub["dose_um"],
                sub["median_growth"],
                marker="o",
                linewidth=2,
                label=label,
                color=palette[label],
            )
            ax.fill_between(
                sub["dose_um"],
                sub["q1_growth"],
                sub["q3_growth"],
                alpha=0.18,
                color=palette[label],
            )
        ax.set_xscale("log")
        ax.set_xlabel("Dose (uM, log scale)")
        ax.set_ylabel("Growth fraction (2 ** median_l2fc)")
        ax.set_title("AML group curves shift lower at higher doses")
        ax.legend(frameon=True)
        return fig

    def plot_representative_curves(curve_df: pd.DataFrame):
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
        panel_map = {
            "Sensitive examples": ["AML — sensitive", "Non-AML — sensitive"],
            "Typical examples": ["AML — typical", "Non-AML — typical"],
        }
        palette = {
            "AML — sensitive": "#b2182b",
            "AML — typical": "#d6604d",
            "Non-AML — sensitive": "#2166ac",
            "Non-AML — typical": "#4393c3",
        }
        for ax, (title, labels) in zip(axes, panel_map.items()):
            for label in labels:
                sub = curve_df.loc[curve_df["example_label"].eq(label)].sort_values("dose_um")
                ax.plot(
                    sub["dose_um"],
                    sub["growth_fraction"],
                    marker="o",
                    linewidth=2,
                    label=f"{label}: {sub['cell_line_name'].iloc[0]}",
                    color=palette[label],
                )
            ax.set_xscale("log")
            ax.set_title(title)
            ax.set_xlabel("Dose (uM, log scale)")
            ax.legend(fontsize=8, frameon=True)
        axes[0].set_ylabel("Growth fraction")
        fig.suptitle("Representative AML and non-AML curves", y=1.02)
        fig.tight_layout()
        return fig

    def plot_curve_features(features_df: pd.DataFrame):
        feature_cols = [
            ("high_dose_growth", "High-dose growth"),
            ("dynamic_range_growth", "Dynamic range"),
            ("steepest_drop_per_log10_uM", "Steepest drop / log10 uM"),
            ("logdose_auc_growth", "Dose-level AUC of growth"),
        ]
        fig, axes = plt.subplots(2, 2, figsize=(10, 7))
        palette = {"AML": "#b2182b", "Non-AML": "#4c78a8"}
        for ax, (col, title) in zip(axes.flat, feature_cols):
            data = [
                features_df.loc[features_df["aml_status"].eq("AML"), col].dropna(),
                features_df.loc[features_df["aml_status"].eq("Non-AML"), col].dropna(),
            ]
            bp = ax.boxplot(data, tick_labels=["AML", "Non-AML"], patch_artist=True)
            for patch, label in zip(bp["boxes"], ["AML", "Non-AML"]):
                patch.set_facecolor(palette[label])
                patch.set_alpha(0.75)
            ax.set_title(title)
        fig.suptitle("First-pass curve-shape features from dose-level data", y=1.02)
        fig.tight_layout()
        return fig

    return (
        add_analysis_labels,
        build_curve_features,
        plot_aml_waterfall,
        plot_curve_features,
        plot_group_curves,
        plot_lineage_comparison,
        plot_representative_curves,
        select_representative_models,
        summarize_lineages,
    )


@app.cell
def _(
    DB_PATH,
    POLARS_DIR,
    SCREEN_ID,
    add_analysis_labels,
    pl,
    prepare_lazy_tables,
):
    tables = prepare_lazy_tables(
        output_dir=POLARS_DIR,
        db_path=DB_PATH,
        tables=["models", "drug_response_secondary", "drug_response_secondary_dose"],
    )

    responses = (
        tables["drug_response_secondary"]
        .filter(pl.col("screen_id") == SCREEN_ID)
        .join(
            tables["models"].select(
                [
                    "model_id",
                    "cell_line_name",
                    "oncotree_lineage",
                    "oncotree_primary_disease",
                ]
            ),
            on="model_id",
            how="left",
        )
        .collect()
        .to_pandas()
    )
    responses = add_analysis_labels(responses)

    doses = (
        tables["drug_response_secondary_dose"]
        .filter(pl.col("screen_id") == SCREEN_ID)
        .join(
            tables["models"].select(
                [
                    "model_id",
                    "cell_line_name",
                    "oncotree_lineage",
                    "oncotree_primary_disease",
                ]
            ),
            on="model_id",
            how="left",
        )
        .collect()
        .to_pandas()
    )
    doses = add_analysis_labels(doses)
    doses["growth_fraction"] = 2.0 ** doses["median_l2fc"]
    return doses, responses


@app.cell
def _(AML_DISEASE, doses, mo, responses):
    fitted_n = responses["auc"].notna().sum()
    aml_n = responses.loc[responses["oncotree_primary_disease"].eq(AML_DISEASE), "model_id"].nunique()
    non_aml_n = responses.loc[~responses["oncotree_primary_disease"].eq(AML_DISEASE), "model_id"].nunique()
    dose_levels = sorted(doses["dose_um"].dropna().unique().tolist())

    mo.md(
        f"""
        ## Data coverage

        - Screen: **`{responses['screen_id'].dropna().iloc[0]}`**
        - Total response rows: **{len(responses):,}**
        - Response rows with fitted AUC: **{fitted_n:,}**
        - AML models in this screen: **{aml_n}**
        - Non-AML models in this screen: **{non_aml_n}**
        - Dose levels (uM): **{', '.join(f'{x:.6g}' for x in dose_levels)}**
        """
    )
    return


@app.cell
def _(min_lineage_n, responses, summarize_lineages):
    lineage_summary = summarize_lineages(responses, min_n=min_lineage_n.value)
    aml_vs_rest = (
        responses.dropna(subset=["auc"])
        .groupby("aml_status", as_index=False)
        .agg(
            n=("model_id", "size"),
            median_auc=("auc", "median"),
            mean_auc=("auc", "mean"),
            median_ic50=("ic50", "median"),
        )
        .sort_values("median_auc")
    )
    return aml_vs_rest, lineage_summary


@app.cell
def _(aml_vs_rest, lineage_summary, mo, plot_lineage_comparison):
    mo.md("## 1. Is AML more sensitive than other lineages?")
    lineage_plot = plot_lineage_comparison(lineage_summary)
    mo.vstack(
        [
            mo.md(
                "AML is broken out explicitly rather than being absorbed into the broader `Myeloid` lineage."
            ),
            aml_vs_rest,
            lineage_summary,
            lineage_plot,
        ]
    )
    return


@app.cell
def _(AML_DISEASE, mo, plot_aml_waterfall, responses):
    aml_responses = responses.loc[
        responses["oncotree_primary_disease"].eq(AML_DISEASE)
    ].copy()
    aml_waterfall = plot_aml_waterfall(aml_responses)
    aml_summary = (
        aml_responses.dropna(subset=["auc"])
        .sort_values("auc")
        .loc[:, ["cell_line_name", "model_id", "auc", "ic50"]]
        .reset_index(drop=True)
    )
    mo.md("## 2. Is AML internally heterogeneous?")
    mo.vstack(
        [
            mo.md(
                "The waterfall plot makes it easy to see whether AML behaves as a tight cluster or spans a wide sensitivity range."
            ),
            aml_waterfall,
            aml_summary,
        ]
    )
    return


@app.cell
def _(build_curve_features, doses):
    group_curve_summary = (
        doses.groupby(["aml_status", "dose_um"], as_index=False)
        .agg(
            median_growth=("growth_fraction", "median"),
            q1_growth=("growth_fraction", lambda x: x.quantile(0.25)),
            q3_growth=("growth_fraction", lambda x: x.quantile(0.75)),
        )
        .sort_values(["aml_status", "dose_um"])
    )
    curve_features = build_curve_features(doses)
    feature_summary = (
        curve_features.groupby("aml_status", as_index=False)
        .agg(
            n=("model_id", "size"),
            median_high_dose_growth=("high_dose_growth", "median"),
            median_dynamic_range=("dynamic_range_growth", "median"),
            median_steepest_drop=("steepest_drop_per_log10_uM", "median"),
            median_logdose_auc=("logdose_auc_growth", "median"),
        )
    )
    return curve_features, feature_summary, group_curve_summary


@app.cell
def _(
    doses,
    group_curve_summary,
    plot_group_curves,
    plot_representative_curves,
    responses,
    select_representative_models,
):
    representative_models = select_representative_models(responses)
    representative_curves = doses.merge(
        representative_models[["example_label", "model_id"]],
        on="model_id",
        how="inner",
    )
    group_curve_plot = plot_group_curves(group_curve_summary)
    representative_curve_plot = plot_representative_curves(representative_curves)
    return group_curve_plot, representative_curve_plot, representative_models


@app.cell
def _(
    curve_features,
    feature_summary,
    group_curve_plot,
    mo,
    plot_curve_features,
    representative_curve_plot,
    representative_models,
):
    feature_plot = plot_curve_features(curve_features)
    mo.md("## 3. Are AML dose-response curves shaped differently from non-AML curves?")
    mo.vstack(
        [
            mo.md(
                "This section combines group-level dose curves, a few representative examples, and simple shape features derived directly from the eight dose points."
            ),
            group_curve_plot,
            representative_models,
            representative_curve_plot,
            feature_summary,
            feature_plot,
        ]
    )
    return


@app.cell
def _(AML_DISEASE, curve_features, mo, responses):
    aml_auc = responses.loc[
        responses["oncotree_primary_disease"].eq(AML_DISEASE) & responses["auc"].notna(),
        "auc",
    ]
    non_aml_auc = responses.loc[
        ~responses["oncotree_primary_disease"].eq(AML_DISEASE) & responses["auc"].notna(),
        "auc",
    ]

    aml_feature_row = (
        curve_features.groupby("aml_status", as_index=False)
        .agg(
            median_high_dose_growth=("high_dose_growth", "median"),
            median_dynamic_range=("dynamic_range_growth", "median"),
            median_steepest_drop=("steepest_drop_per_log10_uM", "median"),
            median_logdose_auc=("logdose_auc_growth", "median"),
        )
        .set_index("aml_status")
    )

    mo.md(
        f"""
        ## Provisional takeaways

        1. **AML looks more sensitive than most comparison groups in this screen.**  
           AML median fitted AUC is **{aml_auc.median():.3f}** versus **{non_aml_auc.median():.3f}** for non-AML models, and AML sits near the most sensitive end of the lineage ranking when separated from the rest of Myeloid.

        2. **AML is clearly heterogeneous.**  
           The AML waterfall spans from very strong responders to near-insensitive models, so this does not look like a uniformly AML-wide effect.

        3. **AML curves look deeper and somewhat steeper at the high-dose end.**  
           Relative to non-AML, AML shows:
           - lower median high-dose growth (**{aml_feature_row.loc['AML', 'median_high_dose_growth']:.3f}** vs **{aml_feature_row.loc['Non-AML', 'median_high_dose_growth']:.3f}**)
           - larger median dynamic range (**{aml_feature_row.loc['AML', 'median_dynamic_range']:.3f}** vs **{aml_feature_row.loc['Non-AML', 'median_dynamic_range']:.3f}**)
           - more negative steepest dose step (**{aml_feature_row.loc['AML', 'median_steepest_drop']:.3f}** vs **{aml_feature_row.loc['Non-AML', 'median_steepest_drop']:.3f}**)
           - lower dose-level growth AUC (**{aml_feature_row.loc['AML', 'median_logdose_auc']:.3f}** vs **{aml_feature_row.loc['Non-AML', 'median_logdose_auc']:.3f}**)

        These are descriptive phase-1 results, not a final biomarker model.
        Good next steps would be AML subtype stratification, integration with genomic features, and replicate-aware sensitivity modeling.
        """
    )
    return


if __name__ == "__main__":
    app.run()
