import marimo

__generated_with = "0.20.4"
app = marimo.App(width="wide")


@app.cell
def _():
    import math
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
    return (
        PROJECT_ROOT,
        math,
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
    # C-604 / MASTL AML notebook — endpoint-framed phase 1

    This notebook intentionally treats **2 µM growth fraction** as the primary phenotype.

    Why this framing changed:

    - The screen has eight doses, but the largest informative gap is between **0.6667 and 2.0 µM**.
    - That is exactly the range where many response transitions appear to occur.
    - Because that transition region is undersampled, **curve-shape taxonomy, transition-dose inference, and AUC-first language are too confident for the design**.

    So the notebook now asks:

    1. **Is AML sensitive at the 2 µM endpoint relative to other groups?**
    2. **Is AML heterogeneous at 2 µM?**
    3. **What candidate molecular features separate strong versus weak AML 2 µM responders?**

    Conventions:

    - Primary phenotype: **growth fraction = `2 ** median_l2fc` at 2 µM**.
      Lower values mean stronger inhibition.
    - Fitted `auc`, `ic50`, and `ec50` remain visible as **secondary context only**.
    - Dose-series plots are retained as descriptive background, not as mechanistic evidence about curve class.
    """)
    return


@app.cell
def _(PROJECT_ROOT, mo):
    from pathlib import Path

    SCREEN_ID = "MTS028_BIAS_CORRECTED"
    AML_DISEASE = "Acute Myeloid Leukemia"
    TOP_DOSE_UM = 2.0
    DB_PATH = Path.home() / ".depmap" / "depmap.duckdb"
    POLARS_DIR = PROJECT_ROOT / "data" / "polars"

    min_lineage_n = mo.ui.slider(
        start=5,
        stop=30,
        step=1,
        value=8,
        label="Minimum group size for lineage ranking",
    )
    responder_tail_pct = mo.ui.slider(
        start=25,
        stop=40,
        step=5,
        value=33,
        label="Responder tail size (% of AML models per strong/weak tail)",
    )
    mo.hstack([min_lineage_n, responder_tail_pct], justify="start")
    return (
        AML_DISEASE,
        DB_PATH,
        POLARS_DIR,
        SCREEN_ID,
        TOP_DOSE_UM,
        min_lineage_n,
        responder_tail_pct,
    )


@app.cell
def _(AML_DISEASE, TOP_DOSE_UM, math, np, pd, plt):
    RESPONSE_META_COLS = {
        "model_id",
        "cell_line_name",
        "oncotree_lineage",
        "oncotree_primary_disease",
        "oncotree_subtype",
        "response_id",
        "screen_id",
        "auc",
        "ic50",
        "ec50",
        "dose_um",
        "median_l2fc",
        "growth_fraction",
        "analysis_group",
        "aml_status",
        "is_aml",
        "responder_group",
    }

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

    def summarize_lineages_at_top_dose(
        top_dose_df: pd.DataFrame,
        response_df: pd.DataFrame,
        *,
        min_n: int = 8,
    ) -> pd.DataFrame:
        dose_summary = (
            top_dose_df.groupby("analysis_group", as_index=False)
            .agg(
                n_models=("model_id", "nunique"),
                median_growth_fraction=("growth_fraction", "median"),
                mean_growth_fraction=("growth_fraction", "mean"),
                q1_growth_fraction=(
                    "growth_fraction",
                    lambda x: x.quantile(0.25),
                ),
                q3_growth_fraction=(
                    "growth_fraction",
                    lambda x: x.quantile(0.75),
                ),
            )
        )
        auc_summary = (
            response_df.dropna(subset=["auc"])
            .groupby("analysis_group", as_index=False)
            .agg(
                median_auc=("auc", "median"),
                mean_auc=("auc", "mean"),
            )
        )
        return (
            dose_summary.merge(auc_summary, on="analysis_group", how="left")
            .query("n_models >= @min_n")
            .sort_values(
                ["median_growth_fraction", "n_models"],
                ascending=[True, False],
            )
            .reset_index(drop=True)
        )

    def assign_aml_responder_groups(
        aml_top_dose_df: pd.DataFrame,
        *,
        tail_fraction: float,
    ) -> tuple[pd.DataFrame, dict[str, float]]:
        lower_q = float(
            aml_top_dose_df["growth_fraction"].quantile(tail_fraction)
        )
        upper_q = float(
            aml_top_dose_df["growth_fraction"].quantile(1.0 - tail_fraction)
        )

        out = aml_top_dose_df.copy()
        out["responder_group"] = "Middle"
        out.loc[out["growth_fraction"] <= lower_q, "responder_group"] = (
            "Strong responder"
        )
        out.loc[out["growth_fraction"] >= upper_q, "responder_group"] = (
            "Weak responder"
        )
        return out, {
            "strong_max_growth_fraction": lower_q,
            "weak_min_growth_fraction": upper_q,
        }

    def summarize_responder_groups(
        aml_grouped_df: pd.DataFrame,
        aml_response_df: pd.DataFrame,
    ) -> pd.DataFrame:
        summary = (
            aml_grouped_df.groupby("responder_group", as_index=False)
            .agg(
                n_models=("model_id", "nunique"),
                median_growth_fraction=("growth_fraction", "median"),
                min_growth_fraction=("growth_fraction", "min"),
                max_growth_fraction=("growth_fraction", "max"),
            )
            .sort_values("median_growth_fraction")
            .reset_index(drop=True)
        )
        auc_summary = (
            aml_grouped_df.merge(
                aml_response_df.loc[:, ["model_id", "auc"]],
                on="model_id",
                how="left",
            )
            .groupby("responder_group", as_index=False)
            .agg(median_auc=("auc", "median"))
        )
        return summary.merge(auc_summary, on="responder_group", how="left")

    def build_aml_model_summary(
        aml_top_dose_df: pd.DataFrame,
        aml_response_df: pd.DataFrame,
    ) -> pd.DataFrame:
        return (
            aml_top_dose_df.merge(
                aml_response_df.loc[:, ["model_id", "auc", "ic50", "ec50"]],
                on="model_id",
                how="left",
            )
            .sort_values("growth_fraction")
            .reset_index(drop=True)
        )

    def scan_wide_feature_table(
        wide_df: pd.DataFrame,
        *,
        group_assignments: pd.DataFrame,
        meta_cols: list[str],
        min_present_per_group: int = 4,
        min_abs_mean_diff: float = 0.5,
        top_n: int = 12,
        direction: str = "weak_minus_strong",
    ) -> pd.DataFrame:
        work = wide_df.copy()
        work = work.merge(
            group_assignments.loc[:, ["model_id", "responder_group"]],
            on="model_id",
            how="inner",
        )
        work = work.loc[
            work["responder_group"].isin(["Strong responder", "Weak responder"])
        ].copy()

        feature_cols = [
            col
            for col in work.columns
            if col not in set(meta_cols) | RESPONSE_META_COLS
        ]

        rows: list[dict[str, float | int | str]] = []
        for col in feature_cols:
            strong = pd.to_numeric(
                work.loc[work["responder_group"].eq("Strong responder"), col],
                errors="coerce",
            ).dropna()
            weak = pd.to_numeric(
                work.loc[work["responder_group"].eq("Weak responder"), col],
                errors="coerce",
            ).dropna()
            if len(strong) < min_present_per_group or len(weak) < min_present_per_group:
                continue

            mean_strong = float(strong.mean())
            mean_weak = float(weak.mean())
            mean_diff = mean_weak - mean_strong
            if abs(mean_diff) < min_abs_mean_diff:
                continue

            var_strong = float(strong.var(ddof=1)) if len(strong) > 1 else float("nan")
            var_weak = float(weak.var(ddof=1)) if len(weak) > 1 else float("nan")
            pooled_sd = float("nan")
            if not math.isnan(var_strong) and not math.isnan(var_weak):
                pooled_num = (len(strong) - 1) * var_strong + (len(weak) - 1) * var_weak
                pooled_den = max(len(strong) + len(weak) - 2, 1)
                if pooled_num > 0:
                    pooled_sd = math.sqrt(pooled_num / pooled_den)

            effect_size = (
                mean_diff / pooled_sd
                if not math.isnan(pooled_sd) and pooled_sd > 0
                else np.nan
            )

            rows.append(
                {
                    "feature_id": col,
                    "n_strong": int(len(strong)),
                    "n_weak": int(len(weak)),
                    "median_strong": float(strong.median()),
                    "median_weak": float(weak.median()),
                    "mean_diff_weak_minus_strong": mean_diff,
                    "effect_size": effect_size,
                }
            )

        if not rows:
            return pd.DataFrame(
                columns=[
                    "feature_id",
                    "n_strong",
                    "n_weak",
                    "median_strong",
                    "median_weak",
                    "mean_diff_weak_minus_strong",
                    "effect_size",
                ]
            )

        out = pd.DataFrame(rows)
        ascending = direction != "weak_minus_strong"
        return (
            out.sort_values(
                ["mean_diff_weak_minus_strong", "effect_size"],
                ascending=[ascending, ascending],
            )
            .head(top_n)
            .reset_index(drop=True)
        )

    def scan_mutation_status(
        mutation_status_df: pd.DataFrame,
        *,
        group_assignments: pd.DataFrame,
        top_n: int = 12,
        direction: str = "weak_minus_strong",
    ) -> pd.DataFrame:
        grouped = group_assignments.loc[
            group_assignments["responder_group"].isin(["Strong responder", "Weak responder"]),
            ["model_id", "responder_group"],
        ].copy()
        if grouped.empty:
            return pd.DataFrame()

        n_strong_total = int(grouped["responder_group"].eq("Strong responder").sum())
        n_weak_total = int(grouped["responder_group"].eq("Weak responder").sum())

        merged = mutation_status_df.merge(grouped, on="model_id", how="inner")
        if merged.empty:
            return pd.DataFrame()

        counts = (
            merged.groupby(["gene_symbol", "responder_group"])["model_id"]
            .nunique()
            .unstack(fill_value=0)
            .reset_index()
        )
        for column in ["Strong responder", "Weak responder"]:
            if column not in counts.columns:
                counts[column] = 0

        counts["strong_freq"] = counts["Strong responder"] / max(n_strong_total, 1)
        counts["weak_freq"] = counts["Weak responder"] / max(n_weak_total, 1)
        counts["weak_minus_strong_freq"] = counts["weak_freq"] - counts["strong_freq"]
        counts = counts.loc[
            (counts["Strong responder"] >= 3) | (counts["Weak responder"] >= 3)
        ].copy()
        ascending = direction != "weak_minus_strong"
        return (
            counts.sort_values("weak_minus_strong_freq", ascending=ascending)
            .head(top_n)
            .reset_index(drop=True)
        )

    def scan_proteomics_long(
        proteomics_df: pd.DataFrame,
        *,
        group_assignments: pd.DataFrame,
        top_n: int = 12,
        direction: str = "weak_minus_strong",
    ) -> pd.DataFrame:
        grouped = group_assignments.loc[
            group_assignments["responder_group"].isin(["Strong responder", "Weak responder"]),
            ["model_id", "responder_group"],
        ].copy()
        merged = proteomics_df.merge(grouped, on="model_id", how="inner")
        if merged.empty:
            return pd.DataFrame()

        summary = (
            merged.groupby(["gene_symbol", "responder_group"], dropna=False)
            .agg(
                n_obs=("abundance", "size"),
                median_abundance=("abundance", "median"),
                mean_abundance=("abundance", "mean"),
            )
            .reset_index()
        )
        wide = summary.pivot(
            index="gene_symbol",
            columns="responder_group",
            values=["n_obs", "median_abundance", "mean_abundance"],
        )
        wide.columns = [f"{a}_{b}" for a, b in wide.columns]
        wide = wide.reset_index()
        required_cols = [
            "n_obs_Strong responder",
            "n_obs_Weak responder",
            "median_abundance_Strong responder",
            "median_abundance_Weak responder",
            "mean_abundance_Strong responder",
            "mean_abundance_Weak responder",
        ]
        for col in required_cols:
            if col not in wide.columns:
                wide[col] = np.nan
        wide = wide.loc[
            (wide["n_obs_Strong responder"] >= 4)
            & (wide["n_obs_Weak responder"] >= 4)
        ].copy()
        wide["mean_diff_weak_minus_strong"] = (
            wide["mean_abundance_Weak responder"] - wide["mean_abundance_Strong responder"]
        )
        ascending = direction != "weak_minus_strong"
        return (
            wide.sort_values("mean_diff_weak_minus_strong", ascending=ascending)
            .head(top_n)
            .reset_index(drop=True)
        )

    def plot_lineage_ranking(lineage_summary: pd.DataFrame):
        plot_df = lineage_summary.sort_values("median_growth_fraction", ascending=False)
        fig, ax = plt.subplots(figsize=(10, max(5, 0.38 * len(plot_df))))
        colors = [
            "#b2182b" if group == "AML" else "#4c78a8"
            for group in plot_df["analysis_group"]
        ]
        ax.barh(
            plot_df["analysis_group"],
            plot_df["median_growth_fraction"],
            color=colors,
            alpha=0.92,
        )
        ax.errorbar(
            plot_df["median_growth_fraction"],
            plot_df["analysis_group"],
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
            ax.text(
                float(row.median_growth_fraction) + 0.01,
                i,
                f"n={row.n_models}",
                va="center",
                fontsize=9,
            )
        ax.set_xlabel(
            f"Median growth fraction at {TOP_DOSE_UM:g} µM (lower = more sensitive)"
        )
        ax.set_ylabel("")
        ax.set_title("Lineage comparison using the 2 µM endpoint")
        return fig

    def plot_aml_waterfall(aml_grouped_df: pd.DataFrame):
        plot_df = aml_grouped_df.sort_values("growth_fraction").reset_index(drop=True)
        color_map = {
            "Strong responder": "#b2182b",
            "Middle": "#bdbdbd",
            "Weak responder": "#2166ac",
        }
        fig, ax = plt.subplots(figsize=(11, 4.8))
        ax.bar(
            np.arange(len(plot_df)),
            plot_df["growth_fraction"],
            color=[color_map[g] for g in plot_df["responder_group"]],
            alpha=0.95,
        )
        ax.axhline(
            plot_df["growth_fraction"].median(),
            color="black",
            linestyle="--",
            linewidth=1,
        )
        ax.set_xlabel("AML models, sorted by 2 µM growth fraction")
        ax.set_ylabel(f"Growth fraction at {TOP_DOSE_UM:g} µM")
        ax.set_title("AML spans a broad response range at the 2 µM endpoint")
        ax.set_xticks([])
        return fig

    def plot_aml_distribution(aml_grouped_df: pd.DataFrame):
        plot_df = aml_grouped_df.copy()
        color_map = {
            "Strong responder": "#b2182b",
            "Middle": "#969696",
            "Weak responder": "#2166ac",
        }
        x_positions = {
            "Strong responder": 0,
            "Middle": 1,
            "Weak responder": 2,
        }
        rng = np.random.default_rng(7)

        fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
        axes[0].hist(plot_df["growth_fraction"], bins=10, color="#7f7f7f", edgecolor="white")
        axes[0].axvline(
            plot_df["growth_fraction"].median(),
            color="black",
            linestyle="--",
            linewidth=1,
        )
        axes[0].set_xlabel(f"Growth fraction at {TOP_DOSE_UM:g} µM")
        axes[0].set_ylabel("AML model count")
        axes[0].set_title("AML distribution at 2 µM")

        for group, sub in plot_df.groupby("responder_group"):
            xs = np.full(len(sub), x_positions[group]) + rng.uniform(-0.11, 0.11, len(sub))
            axes[1].scatter(
                xs,
                sub["growth_fraction"],
                s=55,
                alpha=0.88,
                color=color_map[group],
                label=group,
            )
            axes[1].plot(
                [x_positions[group] - 0.2, x_positions[group] + 0.2],
                [sub["growth_fraction"].median(), sub["growth_fraction"].median()],
                color="black",
                linewidth=2,
            )
        axes[1].set_xticks([0, 1, 2], ["Strong", "Middle", "Weak"])
        axes[1].set_ylabel(f"Growth fraction at {TOP_DOSE_UM:g} µM")
        axes[1].set_title("Strong/weak tails used for endpoint feature scans")
        return fig

    def plot_top_dose_vs_auc(aml_model_summary: pd.DataFrame):
        plot_df = aml_model_summary.dropna(subset=["auc"]).copy()
        fig, ax = plt.subplots(figsize=(6.5, 5))
        ax.scatter(
            plot_df["growth_fraction"],
            plot_df["auc"],
            color="#7f0000",
            alpha=0.8,
            s=55,
        )
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
        ax.set_title("AUC broadly tracks the endpoint ordering, but is not the main readout")
        return fig

    def plot_group_dose_curves(dose_df: pd.DataFrame):
        summary = (
            dose_df.groupby(["aml_status", "dose_um"], as_index=False)
            .agg(
                median_growth_fraction=("growth_fraction", "median"),
                q1_growth_fraction=("growth_fraction", lambda x: x.quantile(0.25)),
                q3_growth_fraction=("growth_fraction", lambda x: x.quantile(0.75)),
            )
            .sort_values(["aml_status", "dose_um"])
        )
        fig, ax = plt.subplots(figsize=(8, 5))
        palette = {"AML": "#b2182b", "Non-AML": "#4c78a8"}
        for label, sub in summary.groupby("aml_status"):
            ax.plot(
                sub["dose_um"],
                sub["median_growth_fraction"],
                marker="o",
                linewidth=2,
                label=label,
                color=palette[label],
            )
            ax.fill_between(
                sub["dose_um"],
                sub["q1_growth_fraction"],
                sub["q3_growth_fraction"],
                alpha=0.18,
                color=palette[label],
            )
        ax.set_xscale("log")
        ax.set_xlabel("Dose (µM, log scale)")
        ax.set_ylabel("Growth fraction")
        ax.set_title("Dose-series context only: AML trends lower at the high end")
        ax.legend(frameon=True)
        return fig

    return (
        add_analysis_labels,
        assign_aml_responder_groups,
        build_aml_model_summary,
        plot_aml_distribution,
        plot_aml_waterfall,
        plot_group_dose_curves,
        plot_lineage_ranking,
        plot_top_dose_vs_auc,
        scan_mutation_status,
        scan_proteomics_long,
        scan_wide_feature_table,
        summarize_lineages_at_top_dose,
        summarize_responder_groups,
    )


@app.cell
def _(
    AML_DISEASE,
    DB_PATH,
    POLARS_DIR,
    SCREEN_ID,
    TOP_DOSE_UM,
    add_analysis_labels,
    pd,
    pl,
    prepare_lazy_datasets,
    prepare_lazy_tables,
):
    required_tables = [
        "models",
        "drug_response_secondary",
        "drug_response_secondary_dose",
        "gene_expression_wide",
        "gene_effects_wide",
        "model_gene_mutation_status",
    ]

    tables = prepare_lazy_tables(
        output_dir=POLARS_DIR,
        db_path=DB_PATH,
        tables=required_tables,
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
    responses = add_analysis_labels(responses)

    doses = (
        tables["drug_response_secondary_dose"]
        .filter(pl.col("screen_id") == SCREEN_ID)
        .join(tables["models"].select(model_cols), on="model_id", how="left")
        .with_columns((2 ** pl.col("median_l2fc")).alias("growth_fraction"))
        .collect()
        .to_pandas()
    )
    doses = add_analysis_labels(doses)

    screen_dose_levels = sorted(doses["dose_um"].dropna().unique().tolist())
    top_dose = doses.loc[doses["dose_um"].eq(TOP_DOSE_UM)].copy()

    aml_model_ids = (
        top_dose.loc[top_dose["oncotree_primary_disease"].eq(AML_DISEASE), "model_id"]
        .drop_duplicates()
        .tolist()
    )

    expression_df = (
        tables["gene_expression_wide"]
        .filter(pl.col("model_id").is_in(aml_model_ids))
        .collect()
        .to_pandas()
    )
    dependency_df = (
        tables["gene_effects_wide"]
        .filter(pl.col("model_id").is_in(aml_model_ids))
        .collect()
        .to_pandas()
    )
    mutation_status_df = (
        tables["model_gene_mutation_status"]
        .filter(pl.col("model_id").is_in(aml_model_ids))
        .collect()
        .to_pandas()
    )

    dataset_status = [
        {
            "surface": "gene_expression_wide (AML subset)",
            "status": "available",
            "rows": len(expression_df),
        },
        {
            "surface": "gene_effects_wide (AML subset)",
            "status": "available",
            "rows": len(dependency_df),
        },
        {
            "surface": "model_gene_mutation_status (AML subset)",
            "status": "available",
            "rows": len(mutation_status_df),
        },
    ]

    proteomics_df = None
    try:
        proteomics_lazy = prepare_lazy_datasets(
            output_dir=POLARS_DIR,
            db_path=DB_PATH,
            datasets=["proteomics_long"],
        )["proteomics_long"]
        proteomics_df = (
            proteomics_lazy.filter(pl.col("model_id").is_in(aml_model_ids))
            .collect()
            .to_pandas()
        )
        dataset_status.append(
            {
                "surface": "proteomics_long (AML subset)",
                "status": "available",
                "rows": len(proteomics_df),
            }
        )
    except Exception as exc:  # pragma: no cover - notebook fallback
        dataset_status.append(
            {
                "surface": "proteomics_long",
                "status": f"unavailable: {type(exc).__name__}",
                "rows": 0,
            }
        )

    dataset_status_df = pd.DataFrame(dataset_status)
    return (
        dataset_status_df,
        dependency_df,
        doses,
        expression_df,
        mutation_status_df,
        proteomics_df,
        responses,
        screen_dose_levels,
        top_dose,
    )


@app.cell
def _(
    AML_DISEASE,
    DB_PATH,
    TOP_DOSE_UM,
    dataset_status_df,
    mo,
    responses,
    screen_dose_levels,
):
    aml_models = responses.loc[
        responses["oncotree_primary_disease"].eq(AML_DISEASE), "model_id"
    ].nunique()
    fitted_n = responses["auc"].notna().sum()

    mo.md(
        f"""
        ## Data coverage

        - Screen: **`{responses['screen_id'].dropna().iloc[0]}`**
        - DuckDB source: **`{DB_PATH}`**
        - Response rows in screen: **{len(responses):,}**
        - Models with fitted AUC available: **{fitted_n:,}**
        - AML models in screen: **{aml_models}**
        - Primary endpoint used here: **{TOP_DOSE_UM:g} µM growth fraction**
        - All dose levels in the screen: **{", ".join(f"{x:.6g}" for x in screen_dose_levels)}**

        The large spacing between **0.6667 and 2.0 µM** is why the notebook avoids curve-class inference.
        Supporting molecular surfaces for the AML feature scan are listed below.
        """
    )
    dataset_status_df
    return


@app.cell
def _(doses, mo, plot_group_dose_curves):
    mo.md("## Dose-series context (descriptive only)")
    mo.vstack(
        [
            mo.md(
                "The dose curves below are useful context, but the endpoint notebook does **not** treat the apparent high-dose transition shape as reliably localised because the final dose step is too wide."
            ),
            plot_group_dose_curves(doses),
        ]
    )
    return


@app.cell
def _(min_lineage_n, responses, summarize_lineages_at_top_dose, top_dose):
    lineage_summary = summarize_lineages_at_top_dose(
        top_dose,
        responses,
        min_n=min_lineage_n.value,
    )
    aml_vs_rest = (
        top_dose.groupby("aml_status", as_index=False)
        .agg(
            n_models=("model_id", "nunique"),
            median_growth_fraction=("growth_fraction", "median"),
            q1_growth_fraction=("growth_fraction", lambda x: x.quantile(0.25)),
            q3_growth_fraction=("growth_fraction", lambda x: x.quantile(0.75)),
        )
        .sort_values("median_growth_fraction")
        .reset_index(drop=True)
    )
    return aml_vs_rest, lineage_summary


@app.cell
def _(aml_vs_rest, lineage_summary, mo, plot_lineage_ranking):
    lineage_plot = plot_lineage_ranking(lineage_summary)
    mo.md("## 1. Is AML sensitive at the 2 µM endpoint?")
    mo.vstack(
        [
            mo.md(
                "This ranking uses **median 2 µM growth fraction** as the primary statistic. Lower values indicate stronger inhibition. AUC is reported only as a cross-check, not as the language of the main claim."
            ),
            aml_vs_rest,
            lineage_summary,
            lineage_plot,
        ]
    )
    return


@app.cell
def _(
    AML_DISEASE,
    assign_aml_responder_groups,
    build_aml_model_summary,
    mo,
    plot_aml_distribution,
    plot_aml_waterfall,
    plot_top_dose_vs_auc,
    responder_tail_pct,
    responses,
    summarize_responder_groups,
    top_dose,
):
    aml_top_dose = top_dose.loc[top_dose["oncotree_primary_disease"].eq(AML_DISEASE)].copy()
    aml_responses = responses.loc[
        responses["oncotree_primary_disease"].eq(AML_DISEASE)
    ].copy()

    tail_fraction = responder_tail_pct.value / 100.0
    aml_grouped, responder_cutoffs = assign_aml_responder_groups(
        aml_top_dose,
        tail_fraction=tail_fraction,
    )
    responder_summary = summarize_responder_groups(aml_grouped, aml_responses)
    aml_model_summary = build_aml_model_summary(aml_grouped, aml_responses)

    mo.md("## 2. Is AML heterogeneous at 2 µM?")
    mo.vstack(
        [
            mo.md(
                f"The strong and weak AML tails below are a **descriptive endpoint split** using the bottom and top **{responder_tail_pct.value}%** of the AML 2 µM distribution. This is suitable for hypothesis generation, not for claiming a final clinical threshold."
            ),
            mo.md(
                f"- Strong responder cutoff: growth fraction ≤ **{responder_cutoffs['strong_max_growth_fraction']:.3f}**\n"
                f"- Weak responder cutoff: growth fraction ≥ **{responder_cutoffs['weak_min_growth_fraction']:.3f}**"
            ),
            responder_summary,
            plot_aml_waterfall(aml_grouped),
            plot_aml_distribution(aml_grouped),
            mo.md(
                "AUC broadly tracks the endpoint ordering, but that agreement is treated as reassurance rather than as the main biological framing."
            ),
            plot_top_dose_vs_auc(aml_model_summary),
            aml_model_summary[
                [
                    "cell_line_name",
                    "model_id",
                    "oncotree_subtype",
                    "responder_group",
                    "growth_fraction",
                    "auc",
                    "ic50",
                ]
            ],
        ]
    )
    return aml_grouped, aml_model_summary, responder_summary


@app.cell
def _(
    aml_grouped,
    dependency_df,
    expression_df,
    mutation_status_df,
    proteomics_df,
    scan_mutation_status,
    scan_proteomics_long,
    scan_wide_feature_table,
):
    strong_weak_assignments = aml_grouped.loc[
        aml_grouped["responder_group"].isin(["Strong responder", "Weak responder"])
    ].copy()

    expression_meta_cols = [
        "model_id",
        "sequencing_id",
        "model_condition_id",
        "is_default_entry",
        "is_default_for_mc",
        "created_at",
    ]
    dependency_meta_cols = ["model_id", "created_at"]

    expression_higher_in_weak = scan_wide_feature_table(
        expression_df,
        group_assignments=strong_weak_assignments,
        meta_cols=expression_meta_cols,
        min_present_per_group=4,
        min_abs_mean_diff=2.0,
        top_n=12,
        direction="weak_minus_strong",
    )
    expression_higher_in_strong = scan_wide_feature_table(
        expression_df,
        group_assignments=strong_weak_assignments,
        meta_cols=expression_meta_cols,
        min_present_per_group=4,
        min_abs_mean_diff=2.0,
        top_n=12,
        direction="strong_minus_weak",
    )

    dependency_higher_in_weak = scan_wide_feature_table(
        dependency_df,
        group_assignments=strong_weak_assignments,
        meta_cols=dependency_meta_cols,
        min_present_per_group=4,
        min_abs_mean_diff=0.3,
        top_n=12,
        direction="weak_minus_strong",
    )
    dependency_higher_in_strong = scan_wide_feature_table(
        dependency_df,
        group_assignments=strong_weak_assignments,
        meta_cols=dependency_meta_cols,
        min_present_per_group=4,
        min_abs_mean_diff=0.3,
        top_n=12,
        direction="strong_minus_weak",
    )

    mutation_enriched_in_weak = scan_mutation_status(
        mutation_status_df,
        group_assignments=strong_weak_assignments,
        top_n=12,
        direction="weak_minus_strong",
    )
    mutation_enriched_in_strong = scan_mutation_status(
        mutation_status_df,
        group_assignments=strong_weak_assignments,
        top_n=12,
        direction="strong_minus_weak",
    )

    if proteomics_df is not None:
        proteomics_higher_in_weak = scan_proteomics_long(
            proteomics_df,
            group_assignments=strong_weak_assignments,
            top_n=12,
            direction="weak_minus_strong",
        )
        proteomics_higher_in_strong = scan_proteomics_long(
            proteomics_df,
            group_assignments=strong_weak_assignments,
            top_n=12,
            direction="strong_minus_weak",
        )
    else:
        proteomics_higher_in_weak = None
        proteomics_higher_in_strong = None
    return (
        dependency_higher_in_strong,
        dependency_higher_in_weak,
        expression_higher_in_strong,
        expression_higher_in_weak,
        mutation_enriched_in_strong,
        mutation_enriched_in_weak,
        proteomics_higher_in_strong,
        proteomics_higher_in_weak,
        strong_weak_assignments,
    )


@app.cell
def _(
    dependency_higher_in_strong,
    dependency_higher_in_weak,
    expression_higher_in_strong,
    expression_higher_in_weak,
    mo,
    mutation_enriched_in_strong,
    mutation_enriched_in_weak,
    proteomics_higher_in_strong,
    proteomics_higher_in_weak,
    strong_weak_assignments,
):
    mo.md("## 3. Endpoint-linked AML feature scan")

    blocks = [
        mo.md(
            f"These are **descriptive 2 µM responder-tail comparisons** across **{len(strong_weak_assignments)}** AML models. They are useful for prioritising follow-up biology, but the notebook does **not** present them as validated biomarkers."
        ),
        mo.md("### Expression features higher in weaker AML responders"),
        expression_higher_in_weak,
        mo.md("### Expression features higher in stronger AML responders"),
        expression_higher_in_strong,
        mo.md("### Dependency features shifted in weaker AML responders"),
        dependency_higher_in_weak,
        mo.md("### Dependency features shifted in stronger AML responders"),
        dependency_higher_in_strong,
        mo.md("### Mutations enriched in weaker AML responders"),
        mutation_enriched_in_weak,
        mo.md("### Mutations enriched in stronger AML responders"),
        mutation_enriched_in_strong,
    ]

    if proteomics_higher_in_weak is not None:
        blocks.extend(
            [
                mo.md("### Proteins higher in weaker AML responders"),
                proteomics_higher_in_weak,
                mo.md("### Proteins higher in stronger AML responders"),
                proteomics_higher_in_strong,
            ]
        )
    else:
        blocks.append(
            mo.md(
                "### Proteomics\n`proteomics_long` was unavailable here, so proteomics is omitted rather than inferred indirectly."
            )
        )

    mo.vstack(blocks)
    return


@app.cell
def _(lineage_summary, mo, mutation_enriched_in_weak, responder_summary):
    aml_lineage_row = lineage_summary.loc[
        lineage_summary["analysis_group"].eq("AML")
    ].iloc[0]

    best_mutation_text = "none under the current filters"
    if len(mutation_enriched_in_weak) > 0:
        top_gene = mutation_enriched_in_weak.iloc[0]
        best_mutation_text = (
            f"{top_gene['gene_symbol']} "
            f"({top_gene['Weak responder']}/{int(responder_summary.loc[responder_summary['responder_group'].eq('Weak responder'), 'n_models'].iloc[0])} weak vs "
            f"{top_gene['Strong responder']}/{int(responder_summary.loc[responder_summary['responder_group'].eq('Strong responder'), 'n_models'].iloc[0])} strong)"
        )

    median_aml_growth = aml_lineage_row["median_growth_fraction"]
    median_aml_auc = aml_lineage_row["median_auc"]
    median_weak = responder_summary.loc[
        responder_summary["responder_group"].eq("Weak responder"),
        "median_growth_fraction",
    ].iloc[0]
    median_strong = responder_summary.loc[
        responder_summary["responder_group"].eq("Strong responder"),
        "median_growth_fraction",
    ].iloc[0]

    mo.md(
        f"""
        ## What survives the endpoint reframing

        1. **Retained: AML is one of the more sensitive groups at 2 µM.**  
           AML median 2 µM growth fraction is **{median_aml_growth:.3f}**. Median AUC is **{median_aml_auc:.3f}**, but that is secondary context rather than the main lens.

        2. **Retained: AML is clearly heterogeneous at 2 µM.**  
           The strong and weak AML tails remain well separated (**{median_strong:.3f}** vs **{median_weak:.3f}** median growth fraction), so the endpoint assay still supports real within-AML spread.

        3. **Retained but still preliminary: endpoint-linked feature tables can nominate candidates.**  
           Under the current filters, the clearest mutation-side weak-responder signal is **{best_mutation_text}**.

        4. **Explicitly weakened versus older curve-first framing:**  
           This notebook does **not** claim a reliable transition dose, cliff location, or AML-specific curve class. The gap from **0.6667 to 2.0 µM** is too large for that level of inference.

        In short: the **2 µM endpoint signal in AML is robust enough to keep**, but the **fine-grained dose-response interpretation is intentionally demoted**.
        """
    )
    return


if __name__ == "__main__":
    app.run()
