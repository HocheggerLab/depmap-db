import marimo

__generated_with = "0.20.4"
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
    # C-604 / MASTL esophagus–stomach notebook — 2 µM endpoint framing

    This notebook now treats the screen as a **2 µM endpoint assay**.

    Why that is the defensible framing:

    - The final part of the dose series jumps from **0.6667 to 2.0 µM**.
    - Many apparent response transitions sit somewhere inside that gap.
    - Because the transition region is undersampled, the notebook should **not** present cliff location, transition dose, or curve-shape taxonomy as primary inference.

    So the main questions become:

    1. **How does the esophagus/stomach lineage behave at 2 µM relative to other groups?**
    2. **Within the lineage, is there a robust endpoint split between ESCC and esophagogastric adenocarcinoma?**
    3. **Does any gastric / esophagogastric state split remain visible at 2 µM, and how strongly can it be claimed?**
    4. **Do comparator summaries still work when reduced to endpoint language?**

    Conventions:

    - Primary phenotype: **growth fraction = `2 ** median_l2fc` at 2 µM**.
    - Lower values indicate stronger inhibition.
    - `auc` is shown only as secondary context.
    - Dose-series plots remain descriptive background only.
    """)
    return


@app.cell
def _(PROJECT_ROOT, mo):
    from pathlib import Path

    SCREEN_ID = "MTS028_BIAS_CORRECTED"
    TARGET_LINEAGE = "Esophagus/Stomach"
    ESCC_DISEASE = "Esophageal Squamous Cell Carcinoma"
    EGA_DISEASE = "Esophagogastric Adenocarcinoma"
    TOP_DOSE_UM = 2.0
    DB_PATH = Path.home() / ".depmap" / "depmap.duckdb"
    POLARS_DIR = PROJECT_ROOT / "data" / "polars"

    min_group_n = mo.ui.slider(
        start=5,
        stop=30,
        step=1,
        value=8,
        label="Minimum group size for endpoint comparator ranking",
    )
    mo.hstack([min_group_n], justify="start")
    return (
        DB_PATH,
        EGA_DISEASE,
        ESCC_DISEASE,
        POLARS_DIR,
        SCREEN_ID,
        TARGET_LINEAGE,
        TOP_DOSE_UM,
        min_group_n,
    )


@app.cell
def _(EGA_DISEASE, ESCC_DISEASE, TARGET_LINEAGE, TOP_DOSE_UM, np, pd, plt):
    def classify_upper_gi_state(df: pd.DataFrame) -> pd.DataFrame:
        out = df.copy()
        subtype = out["oncotree_subtype"].fillna("").str.lower()
        disease = out["oncotree_primary_disease"].fillna("")

        out["is_target_lineage"] = out["oncotree_lineage"].eq(TARGET_LINEAGE)
        out["upper_gi_group"] = np.select(
            [
                disease.eq(ESCC_DISEASE),
                disease.eq(EGA_DISEASE),
            ],
            ["ESCC", "Esophagogastric adenocarcinoma"],
            default=out["oncotree_lineage"].fillna("Unknown"),
        )
        out["upper_gi_state"] = np.select(
            [
                disease.eq(ESCC_DISEASE),
                subtype.str.contains("stomach", regex=False),
                subtype.str.contains("esophag", regex=False)
                | subtype.str.contains("gastroesophageal junction", regex=False),
            ],
            [
                "ESCC",
                "Gastric-state adenocarcinoma",
                "Esophageal-state adenocarcinoma",
            ],
            default="Other adenocarcinoma / mixed",
        )
        return out

    def summarize_endpoint(
        df: pd.DataFrame,
        *,
        group_col: str,
        response_df: pd.DataFrame | None = None,
        min_n: int = 1,
    ) -> pd.DataFrame:
        summary = (
            df.groupby(group_col, as_index=False)
            .agg(
                n_models=("model_id", "nunique"),
                median_growth_fraction=("growth_fraction", "median"),
                mean_growth_fraction=("growth_fraction", "mean"),
                q1_growth_fraction=("growth_fraction", lambda x: x.quantile(0.25)),
                q3_growth_fraction=("growth_fraction", lambda x: x.quantile(0.75)),
            )
            .query("n_models >= @min_n")
            .sort_values(["median_growth_fraction", "n_models"], ascending=[True, False])
            .reset_index(drop=True)
        )
        if response_df is not None:
            auc_summary = (
                response_df.dropna(subset=["auc"])
                .groupby(group_col, as_index=False)
                .agg(median_auc=("auc", "median"))
            )
            summary = summary.merge(auc_summary, on=group_col, how="left")
        return summary

    def build_lineage_comparator(top_dose_df: pd.DataFrame, response_df: pd.DataFrame, *, min_n: int = 5) -> pd.DataFrame:
        comp = top_dose_df.copy()
        comp["comparator_group"] = np.select(
            [
                comp["upper_gi_state"].eq("ESCC"),
                comp["upper_gi_state"].eq("Gastric-state adenocarcinoma"),
                comp["upper_gi_state"].eq("Esophageal-state adenocarcinoma"),
                comp["upper_gi_group"].eq("Esophagogastric adenocarcinoma"),
            ],
            [
                "ESCC",
                "Gastric-state adenocarcinoma",
                "Esophageal-state adenocarcinoma",
                "Esophagogastric adenocarcinoma",
            ],
            default=comp["oncotree_lineage"].fillna("Unknown"),
        )
        return summarize_endpoint(
            comp,
            group_col="comparator_group",
            response_df=response_df,
            min_n=min_n,
        )

    def plot_comparator_ranking(summary: pd.DataFrame):
        plot_df = summary.sort_values("median_growth_fraction", ascending=False)
        fig, ax = plt.subplots(figsize=(10, max(5, 0.38 * len(plot_df))))
        highlight = {
            "ESCC": "#b2182b",
            "Gastric-state adenocarcinoma": "#ef8a62",
            "Esophagogastric adenocarcinoma": "#fddbc7",
            "Esophageal-state adenocarcinoma": "#2166ac",
        }
        colors = [highlight.get(g, "#4c78a8") for g in plot_df.iloc[:, 0]]
        ax.barh(plot_df.iloc[:, 0], plot_df["median_growth_fraction"], color=colors, alpha=0.92)
        ax.errorbar(
            plot_df["median_growth_fraction"],
            plot_df.iloc[:, 0],
            xerr=[
                plot_df["median_growth_fraction"] - plot_df["q1_growth_fraction"],
                plot_df["q3_growth_fraction"] - plot_df["median_growth_fraction"],
            ],
            fmt="none",
            ecolor="black",
            elinewidth=1,
            capsize=2,
        )
        for i, row in enumerate(plot_df.itertuples(index=False)):
            ax.text(float(row.median_growth_fraction) + 0.01, i, f"n={row.n_models}", va="center", fontsize=9)
        ax.set_xlabel(f"Median growth fraction at {TOP_DOSE_UM:g} µM (lower = more sensitive)")
        ax.set_ylabel("")
        ax.set_title("Endpoint comparator ranking")
        return fig

    def plot_upper_gi_waterfall(df: pd.DataFrame):
        plot_df = df.sort_values("growth_fraction").reset_index(drop=True)
        palette = {
            "ESCC": "#b2182b",
            "Gastric-state adenocarcinoma": "#ef8a62",
            "Esophageal-state adenocarcinoma": "#2166ac",
            "Other adenocarcinoma / mixed": "#bdbdbd",
        }
        fig, ax = plt.subplots(figsize=(12, 4.8))
        ax.bar(
            np.arange(len(plot_df)),
            plot_df["growth_fraction"],
            color=[palette.get(g, "#bdbdbd") for g in plot_df["upper_gi_state"]],
            alpha=0.95,
        )
        ax.axhline(plot_df["growth_fraction"].median(), color="black", linestyle="--", linewidth=1)
        ax.set_xlabel("Esophagus/stomach models, sorted by 2 µM growth fraction")
        ax.set_ylabel(f"Growth fraction at {TOP_DOSE_UM:g} µM")
        ax.set_title("Upper-GI endpoint waterfall")
        ax.set_xticks([])
        return fig

    def plot_group_dose_curves(dose_df: pd.DataFrame, *, group_col: str, title: str):
        summary = (
            dose_df.groupby([group_col, "dose_um"], as_index=False)
            .agg(
                median_growth_fraction=("growth_fraction", "median"),
                q1_growth_fraction=("growth_fraction", lambda x: x.quantile(0.25)),
                q3_growth_fraction=("growth_fraction", lambda x: x.quantile(0.75)),
            )
            .sort_values([group_col, "dose_um"])
        )
        fig, ax = plt.subplots(figsize=(8.5, 5))
        for label, sub in summary.groupby(group_col):
            ax.plot(
                sub["dose_um"],
                sub["median_growth_fraction"],
                marker="o",
                linewidth=2,
                label=label,
            )
            ax.fill_between(
                sub["dose_um"],
                sub["q1_growth_fraction"],
                sub["q3_growth_fraction"],
                alpha=0.15,
            )
        ax.set_xscale("log")
        ax.set_xlabel("Dose (µM, log scale)")
        ax.set_ylabel("Growth fraction")
        ax.set_title(title)
        ax.legend(frameon=True)
        return fig

    def plot_endpoint_vs_auc(df: pd.DataFrame, *, title: str):
        plot_df = df.dropna(subset=["auc"]).copy()
        fig, ax = plt.subplots(figsize=(6.5, 5))
        ax.scatter(plot_df["growth_fraction"], plot_df["auc"], color="#7f0000", alpha=0.8, s=55)
        for row in plot_df.nsmallest(5, "growth_fraction").itertuples(index=False):
            ax.annotate(
                row.cell_line_name,
                (row.growth_fraction, row.auc),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=8,
            )
        ax.set_xlabel(f"Growth fraction at {TOP_DOSE_UM:g} µM")
        ax.set_ylabel("Fitted AUC (secondary context)")
        ax.set_title(title)
        return fig

    return (
        build_lineage_comparator,
        classify_upper_gi_state,
        plot_comparator_ranking,
        plot_endpoint_vs_auc,
        plot_group_dose_curves,
        plot_upper_gi_waterfall,
        summarize_endpoint,
    )


@app.cell
def _(
    DB_PATH,
    POLARS_DIR,
    SCREEN_ID,
    TOP_DOSE_UM,
    classify_upper_gi_state,
    pl,
    prepare_lazy_tables,
):
    tables = prepare_lazy_tables(
        output_dir=POLARS_DIR,
        db_path=DB_PATH,
        tables=["models", "drug_response_secondary", "drug_response_secondary_dose"],
    )

    model_cols = [
        "model_id",
        "cell_line_name",
        "oncotree_lineage",
        "oncotree_primary_disease",
        "oncotree_subtype",
    ]

    responses = (
        tables["drug_response_secondary"]
        .filter(pl.col("screen_id") == SCREEN_ID)
        .join(tables["models"].select(model_cols), on="model_id", how="left")
        .collect()
        .to_pandas()
    )
    responses = classify_upper_gi_state(responses)

    doses = (
        tables["drug_response_secondary_dose"]
        .filter(pl.col("screen_id") == SCREEN_ID)
        .join(tables["models"].select(model_cols), on="model_id", how="left")
        .with_columns((2 ** pl.col("median_l2fc")).alias("growth_fraction"))
        .collect()
        .to_pandas()
    )
    doses = classify_upper_gi_state(doses)

    top_dose = doses.loc[doses["dose_um"].eq(TOP_DOSE_UM)].copy()
    screen_dose_levels = sorted(doses["dose_um"].dropna().unique().tolist())
    return doses, responses, screen_dose_levels, top_dose


@app.cell
def _(DB_PATH, TARGET_LINEAGE, TOP_DOSE_UM, mo, responses, screen_dose_levels, top_dose):
    cohort_n = top_dose.loc[top_dose["oncotree_lineage"].eq(TARGET_LINEAGE), "model_id"].nunique()
    fitted_n = responses["auc"].notna().sum()

    mo.md(
        f"""
        ## Data coverage

        - Screen: **`{responses['screen_id'].dropna().iloc[0]}`**
        - DuckDB source: **`{DB_PATH}`**
        - Total response rows in screen: **{len(responses):,}**
        - Models with fitted AUC available: **{fitted_n:,}**
        - Esophagus/stomach models with 2 µM endpoint data: **{cohort_n}**
        - Primary endpoint used here: **{TOP_DOSE_UM:g} µM growth fraction**
        - Dose levels: **{", ".join(f"{x:.6g}" for x in screen_dose_levels)}**

        The key interpretive constraint is the **0.6667 → 2.0 µM gap**, so downstream claims are endpoint-first by design.
        """
    )
    return


@app.cell
def _(TARGET_LINEAGE, doses, mo, plot_group_dose_curves):
    lineage_doses = doses.loc[doses["oncotree_lineage"].eq(TARGET_LINEAGE)].copy()
    descriptive_groups = lineage_doses.loc[
        lineage_doses["upper_gi_group"].isin(["ESCC", "Esophagogastric adenocarcinoma"])
    ].copy()

    mo.md("## Dose-series context (descriptive only)")
    mo.vstack(
        [
            mo.md(
                "The curves below are retained only as context. The notebook does **not** infer a reliable cliff location or transition dose inside the final wide step."
            ),
            plot_group_dose_curves(
                descriptive_groups,
                group_col="upper_gi_group",
                title="ESCC versus esophagogastric adenocarcinoma — dose-series context only",
            ),
        ]
    )
    return


@app.cell
def _(TARGET_LINEAGE, top_dose):
    upper_gi_top_dose = top_dose.loc[top_dose["oncotree_lineage"].eq(TARGET_LINEAGE)].copy()
    return (upper_gi_top_dose,)


@app.cell
def _(responses, summarize_endpoint, upper_gi_top_dose):
    disease_summary = summarize_endpoint(
        upper_gi_top_dose,
        group_col="upper_gi_group",
        response_df=responses,
        min_n=1,
    )
    state_summary = summarize_endpoint(
        upper_gi_top_dose,
        group_col="upper_gi_state",
        response_df=responses,
        min_n=1,
    )
    return disease_summary, state_summary


@app.cell
def _(disease_summary, mo, plot_upper_gi_waterfall, upper_gi_top_dose):
    mo.md("## 1. Upper-GI endpoint heterogeneity")
    mo.vstack(
        [
            mo.md(
                "The waterfall below is ordered by **2 µM growth fraction**. This is the main response phenotype for the cohort."
            ),
            disease_summary,
            plot_upper_gi_waterfall(upper_gi_top_dose),
            upper_gi_top_dose[
                [
                    "cell_line_name",
                    "model_id",
                    "oncotree_primary_disease",
                    "oncotree_subtype",
                    "upper_gi_state",
                    "growth_fraction",
                    "auc",
                    "ic50",
                ]
            ],
        ]
    )
    return


@app.cell
def _(mo, plot_group_dose_curves, state_summary, upper_gi_top_dose):
    state_dose_subset = upper_gi_top_dose.loc[
        upper_gi_top_dose["upper_gi_state"].isin(
            ["ESCC", "Gastric-state adenocarcinoma", "Esophageal-state adenocarcinoma"]
        )
    ].copy()

    mo.md("## 2. Does a gastric / esophagogastric state split survive at 2 µM?")
    mo.vstack(
        [
            mo.md(
                "Yes, but the strength of the claim depends on which split is meant. The **robust endpoint split** is **ESCC vs esophagogastric adenocarcinoma**. Within adenocarcinoma, a gastric-state versus esophageal-state direction is still visible, but the esophageal-state subgroup is very small here and stays descriptive only."
            ),
            state_summary,
            plot_group_dose_curves(
                state_dose_subset,
                group_col="upper_gi_state",
                title="Upper-GI state summaries across dose (descriptive only)",
            ),
        ]
    )
    return


@app.cell
def _(
    build_lineage_comparator,
    min_group_n,
    mo,
    plot_comparator_ranking,
    responses,
    top_dose,
):
    comparator_summary = build_lineage_comparator(
        top_dose,
        responses,
        min_n=min_group_n.value,
    )

    mo.md("## 3. Endpoint comparator logic")
    mo.vstack(
        [
            mo.md(
                "The comparator table still works under endpoint framing. In other words, the notebook can place ESCC and gastric-state adenocarcinoma into the broader 2 µM sensitivity landscape without relying on curve classes."
            ),
            comparator_summary,
            plot_comparator_ranking(comparator_summary),
        ]
    )
    return comparator_summary


@app.cell
def _(mo, plot_endpoint_vs_auc, upper_gi_top_dose):
    mo.md("## 4. AUC as secondary context only")
    mo.vstack(
        [
            mo.md(
                "AUC is kept only as a supporting summary. It broadly tracks the endpoint ordering, but the notebook avoids treating it as the primary biology language because the high-dose transition region is too sparsely sampled."
            ),
            plot_endpoint_vs_auc(
                upper_gi_top_dose,
                title="Upper-GI 2 µM endpoint versus fitted AUC",
            ),
        ]
    )
    return


@app.cell
def _(comparator_summary, disease_summary, mo, state_summary):
    disease_row_escc = disease_summary.loc[
        disease_summary["upper_gi_group"].eq("ESCC")
    ].iloc[0]
    disease_row_ega = disease_summary.loc[
        disease_summary["upper_gi_group"].eq("Esophagogastric adenocarcinoma")
    ].iloc[0]
    gastric_state_row = state_summary.loc[
        state_summary["upper_gi_state"].eq("Gastric-state adenocarcinoma")
    ].iloc[0]
    esophageal_state_row = state_summary.loc[
        state_summary["upper_gi_state"].eq("Esophageal-state adenocarcinoma")
    ].iloc[0]

    escc_rank = (
        comparator_summary.reset_index(drop=True)
        .assign(rank=lambda x: x.index + 1)
        .loc[lambda x: x["comparator_group"].eq("ESCC"), "rank"]
        .iloc[0]
    )
    gastric_rank = (
        comparator_summary.reset_index(drop=True)
        .assign(rank=lambda x: x.index + 1)
        .loc[
            lambda x: x["comparator_group"].eq("Gastric-state adenocarcinoma"),
            "rank",
        ]
        .iloc[0]
    )

    mo.md(
        f"""
        ## What is retained versus weakened after reframing

        1. **Retained: ESCC is more sensitive than esophagogastric adenocarcinoma at 2 µM.**  
           ESCC median growth fraction is **{disease_row_escc['median_growth_fraction']:.3f}** (**n={int(disease_row_escc['n_models'])}**) versus **{disease_row_ega['median_growth_fraction']:.3f}** for esophagogastric adenocarcinoma (**n={int(disease_row_ega['n_models'])}**).

        2. **Retained, but with narrower wording: endpoint comparator logic still works.**  
           In the broader ranking, **ESCC** sits at rank **{int(escc_rank)}** and **gastric-state adenocarcinoma** at rank **{int(gastric_rank)}** among comparator groups that meet the minimum size threshold.

        3. **Retained only descriptively: a gastric / esophagogastric state split is still visible at 2 µM.**  
           Gastric-state adenocarcinoma has median growth fraction **{gastric_state_row['median_growth_fraction']:.3f}** (**n={int(gastric_state_row['n_models'])}**), while esophageal-state adenocarcinoma is **{esophageal_state_row['median_growth_fraction']:.3f}** (**n={int(esophageal_state_row['n_models'])}**).  
           The **direction** survives, but the esophageal-state subgroup is too small here for a strong inferential claim.

        4. **Explicitly weakened from older curve-first language:**  
           This notebook does **not** claim a reliable transition dose, cliff position, or mechanistic curve class within upper-GI models. The endpoint assay supports endpoint statements; it does not support overconfident interpolation inside the final dose gap.

        In short: the scientifically defensible biology here is **endpoint sensitivity at 2 µM**, especially the **ESCC vs adenocarcinoma split**. Any finer state decomposition inside adenocarcinoma is kept **visible but clearly demoted to descriptive status** when sample size is thin.
        """
    )
    return


if __name__ == "__main__":
    app.run()
