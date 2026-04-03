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

    PROJECT_ROOT = Path(__file__).resolve().parents[1]
    SRC_ROOT = PROJECT_ROOT / "src"
    if str(SRC_ROOT) not in sys.path:
        sys.path.insert(0, str(SRC_ROOT))

    from depmap_db.analysis.c604_mastl_overlap import (
        POSITIVE_FDR_CUTOFF,
        TOP_DOSE_UM,
        run_analysis,
    )

    plt.style.use("seaborn-v0_8-whitegrid")
    return (
        POSITIVE_FDR_CUTOFF,
        PROJECT_ROOT,
        TOP_DOSE_UM,
        mo,
        np,
        pd,
        plt,
        run_analysis,
    )


@app.cell
def _(mo):
    mo.md(r"""
    # C-604 endpoint vs MASTL dependency network overlap

    This notebook asks whether **C-604 sensitivity** and **MASTL / Greatwall dependency** converge on a coherent dependency network.

    ## How the two ranked lists are built

    1. **C-604 2 µM endpoint vs dependency strength**
       - Start from the **2 µM endpoint** for the MTS028 screen.
       - Convert `median_l2fc` to **growth fraction** via `2 ** median_l2fc`.
       - Define **C-604 sensitivity** as `-growth_fraction`, so **larger values mean stronger drug sensitivity**.
       - Convert DepMap gene effect to **dependency strength** via `-gene_effect`, so **larger values mean stronger dependency**.
       - For each gene, compute the Pearson correlation between **C-604 sensitivity** and **that gene's dependency strength** across models with both measurements.
       - Positive correlation therefore means: **models that are more C-604-sensitive also tend to be more dependent on that gene**.

    2. **MASTL dependency profile vs all other dependency strengths**
       - Use the same dependency-strength transform: `-gene_effect`.
       - Treat **MASTL dependency strength** as the anchor profile.
       - For each other gene, compute the Pearson correlation between **MASTL dependency strength** and that gene's dependency strength.
       - Positive correlation means: **models more dependent on MASTL also tend to be more dependent on that gene**.

    ## Pathway analysis choice

    The local setup does not ship a GSEA package, so the notebook uses a **defensible preranked enrichment equivalent**:

    - rank genes by the signed correlation coefficient,
    - test whether pathway genes sit preferentially near the **positive/top** end of the ranking,
    - use a **one-sided Mann–Whitney rank enrichment** statistic plus Benjamini–Hochberg FDR.

    Pathway libraries are pulled from Enrichr text exports and cached locally on first run.
    """)
    return


@app.cell
def _(PROJECT_ROOT, run_analysis):
    outputs = run_analysis(project_root=PROJECT_ROOT)
    return (outputs,)


@app.cell
def _(outputs, pd):
    analysis_summary = outputs["analysis_summary"]
    c604_endpoint = outputs["c604_endpoint"]
    top100_c604 = outputs["top100_c604"]
    top100_mastl = outputs["top100_mastl"]
    overlap_table = outputs["overlap_table"]
    gene_overlap_summary = outputs["gene_overlap_summary"]
    c604_pathways = outputs["c604_pathways"]
    mastl_pathways = outputs["mastl_pathways"]
    pathway_concordance_summary = outputs["pathway_concordance_summary"]
    shared_pathways = outputs["shared_pathways"]

    lineage_summary = (
        c604_endpoint.groupby("oncotree_lineage", as_index=False)
        .agg(
            n_models=("model_id", "nunique"),
            median_growth_fraction=("growth_fraction", "median"),
        )
        .sort_values(["n_models", "median_growth_fraction"], ascending=[False, True])
        .reset_index(drop=True)
    )
    return (
        analysis_summary,
        c604_endpoint,
        c604_pathways,
        gene_overlap_summary,
        lineage_summary,
        mastl_pathways,
        overlap_table,
        pathway_concordance_summary,
        shared_pathways,
        top100_c604,
        top100_mastl,
    )


@app.cell
def _(TOP_DOSE_UM, analysis_summary, gene_overlap_summary, mo):
    summary = analysis_summary.iloc[0]
    overlap = gene_overlap_summary.iloc[0]
    mo.md(
        f"""
        ## Dataset summary

        - **C-604 endpoint models in screen:** {int(summary['c604_screen_models'])}
        - **C-604 models with matched dependency profiles:** {int(summary['c604_models_with_dependency'])}
        - **Models in the MASTL dependency network:** {int(summary['mastl_models_with_dependency'])}
        - **Common gene universe for ranking overlap:** {int(summary['gene_universe'])}
        - **Top C-604-associated dependency:** **{summary['top_c604_gene']}** (r = {summary['top_c604_r']:.3f})
        - **Top MASTL co-dependency:** **{summary['top_mastl_gene']}** (r = {summary['top_mastl_r']:.3f})
        - **Shared genes among the top 100 lists:** **{int(overlap['shared_top_genes'])}** / 100
        - **Top-100 overlap hypergeometric p-value:** **{overlap['shared_top_gene_p_value']:.3g}**

        Reminder: this is all framed around the **{TOP_DOSE_UM:g} µM endpoint**. Lower growth fraction means more inhibition, but the correlation ranking uses **`-growth_fraction`** so that **larger = more sensitive**.
        """
    )
    return


@app.cell
def _(lineage_summary, mo):
    mo.vstack(
        [
            mo.md("### C-604 endpoint lineage coverage"),
            lineage_summary.head(15),
        ]
    )
    return


@app.cell
def _(overlap_table, pd, top100_c604, top100_mastl):
    scatter_df = top100_c604.loc[:, ["gene", "pearson_r", "rank_positive"]].merge(
        top100_mastl.loc[:, ["gene", "pearson_r", "rank_positive"]],
        on="gene",
        suffixes=("_c604", "_mastl"),
        how="outer",
    )
    scatter_df["shared_top100"] = scatter_df["gene"].isin(overlap_table["gene"])
    return (scatter_df,)


@app.cell
def _(np, overlap_table, plt, scatter_df):
    plot_df = scatter_df.dropna(subset=["pearson_r_c604", "pearson_r_mastl"]).copy()
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(
        plot_df["pearson_r_c604"],
        plot_df["pearson_r_mastl"],
        color="#bdbdbd",
        s=36,
        alpha=0.75,
    )
    shared = plot_df.loc[plot_df["shared_top100"]].copy()
    ax.scatter(
        shared["pearson_r_c604"],
        shared["pearson_r_mastl"],
        color="#b2182b",
        s=58,
        alpha=0.95,
    )
    for row in shared.itertuples(index=False):
        ax.text(row.pearson_r_c604 + 0.002, row.pearson_r_mastl + 0.002, row.gene, fontsize=9)
    ax.axvline(0, color="black", linewidth=1, linestyle="--")
    ax.axhline(0, color="black", linewidth=1, linestyle="--")
    ax.set_xlabel("C-604 sensitivity vs dependency strength (Pearson r)")
    ax.set_ylabel("MASTL dependency vs dependency strength (Pearson r)")
    ax.set_title("Shared genes between the two top-100 positive tails")
    fig
    return (fig,)


@app.cell
def _(gene_overlap_summary, mo, overlap_table):
    overlap = gene_overlap_summary.iloc[0]
    mo.vstack(
        [
            mo.md(
                f"### Gene-level overlap\n\nShared top-100 genes: **{int(overlap['shared_top_genes'])}** with hypergeometric **p = {overlap['shared_top_gene_p_value']:.3g}**"
            ),
            overlap_table,
        ]
    )
    return


@app.cell
def _(mo, top100_c604, top100_mastl):
    mo.vstack(
        [
            mo.md("### Top 100 positive dependency correlations for the C-604 endpoint"),
            top100_c604,
            mo.md("### Top 100 positive co-dependencies for the MASTL dependency profile"),
            top100_mastl,
        ]
    )
    return


@app.cell
def _(POSITIVE_FDR_CUTOFF, c604_pathways, mastl_pathways, mo, pathway_concordance_summary):
    mo.md(
        f"""
        ## Pathway-level overlap and concordance

        Positive pathway signal is called at **FDR < {POSITIVE_FDR_CUTOFF:.2f}**.

        The table below summarizes two things for each library:

        - how many pathways are significant in each ranking,
        - how many are shared,
        - and the global **Spearman concordance** of pathway scores across all tested terms.
        """
    )
    pathway_concordance_summary
    return


@app.cell
def _(c604_pathways, mastl_pathways, np, plt):
    def top_positive(df, library, n=12):
        return (
            df.loc[df["library"].eq(library)]
            .sort_values(["fdr", "mean_r"], ascending=[True, False])
            .head(n)
            .iloc[::-1]
        )

    fig, axes = plt.subplots(1, 2, figsize=(14, 8), sharex=False)
    c604_top = top_positive(c604_pathways, "Reactome_2022")
    mastl_top = top_positive(mastl_pathways, "Reactome_2022")

    axes[0].barh(c604_top["term"], c604_top["mean_r"], color="#2166ac")
    axes[0].set_title("Top Reactome pathways: C-604 endpoint")
    axes[0].set_xlabel("Mean pathway gene correlation")

    axes[1].barh(mastl_top["term"], mastl_top["mean_r"], color="#b2182b")
    axes[1].set_title("Top Reactome pathways: MASTL dependency")
    axes[1].set_xlabel("Mean pathway gene correlation")

    fig.tight_layout()
    fig
    return (fig,)


@app.cell
def _(c604_pathways, mastl_pathways, plt):
    def top_positive(df, library, n=10):
        return (
            df.loc[df["library"].eq(library)]
            .sort_values(["fdr", "mean_r"], ascending=[True, False])
            .head(n)
            .iloc[::-1]
        )

    fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharex=False)
    c604_top = top_positive(c604_pathways, "MSigDB_Hallmark_2020")
    mastl_top = top_positive(mastl_pathways, "MSigDB_Hallmark_2020")

    axes[0].barh(c604_top["term"], c604_top["mean_r"], color="#4393c3")
    axes[0].set_title("Top Hallmark pathways: C-604 endpoint")
    axes[0].set_xlabel("Mean pathway gene correlation")

    axes[1].barh(mastl_top["term"], mastl_top["mean_r"], color="#d6604d")
    axes[1].set_title("Top Hallmark pathways: MASTL dependency")
    axes[1].set_xlabel("Mean pathway gene correlation")

    fig.tight_layout()
    fig
    return (fig,)


@app.cell
def _(POSITIVE_FDR_CUTOFF, mo, shared_pathways):
    shared_display = shared_pathways.loc[
        :,
        [
            "library",
            "term",
            "mean_r_c604",
            "fdr_c604",
            "mean_r_mastl",
            "fdr_mastl",
        ],
    ].head(40)
    mo.vstack(
        [
            mo.md(
                f"### Shared positive pathways (both rankings at FDR < {POSITIVE_FDR_CUTOFF:.2f})"
            ),
            shared_display,
        ]
    )
    return


@app.cell
def _(mo, pathway_concordance_summary, shared_pathways, top100_c604, top100_mastl):
    reactome = pathway_concordance_summary.loc[
        pathway_concordance_summary["library"].eq("Reactome_2022")
    ].iloc[0]
    hallmark = pathway_concordance_summary.loc[
        pathway_concordance_summary["library"].eq("MSigDB_Hallmark_2020")
    ].iloc[0]

    c604_examples = ", ".join(top100_c604.head(8)["gene"].tolist())
    mastl_examples = ", ".join(top100_mastl.head(8)["gene"].tolist())

    mo.md(
        f"""
        ## Interpretation

        **1. Gene-level convergence is real but modest.**

        The two top-100 positive tails are not identical. That is actually the honest result here. The overlap is small rather than huge, which argues against claiming a single crisp shared gene module. Still, the overlap is not random, and the shared genes include repair / chromosome-end related candidates such as **MRE11**, **TEN1**, and **CDIN1**.

        **2. The C-604 endpoint ranking is headed by MASTL itself.**

        That is the cleanest internal consistency check in the whole analysis: the strongest positive dependency correlation with C-604 sensitivity is **MASTL**. In plain language, models that are more inhibited by C-604 at 2 µM also tend to show stronger genetic dependency on MASTL.

        Leading C-604 hits: {c604_examples}.

        **3. The MASTL dependency network is more canonical cell-cycle / genome-maintenance biology.**

        The MASTL co-dependency list is enriched for centromere, checkpoint, telomere, and DNA repair biology rather than just a vague stress signature.

        Leading MASTL hits: {mastl_examples}.

        **4. Pathway-level convergence is much stronger than gene-level overlap.**

        - **Hallmark:** weaker but still directionally concordant overall (Spearman rho = {hallmark['spearman_rho']:.3f}).
        - **Reactome:** strong positive concordance across pathway scores (Spearman rho = {reactome['spearman_rho']:.3f}), with **{int(reactome['shared_sig_terms'])} shared significant pathways** and a hypergeometric p-value of **{reactome['shared_sig_p_value']:.3g}**.

        The shared Reactome signal is the important part biologically: it includes **homologous recombination / D-loop processing**, **NHEJ**, **telomere packaging / telomere maintenance**, **centromere/CENPA deposition**, **rRNA / translation**, and in the C-604 ranking especially a striking **mitochondrial translation / respiratory chain** component.

        **5. Practical readout for Helfrid.**

        If the question is whether **C-604 sensitivity and MASTL dependency point into the same biological neighborhood**, the answer is **yes at pathway/network scale, but only partially at single-gene scale**.

        ## Caveats

        - The C-604 network is built only on models with both **MTS028 endpoint data** and **DepMap dependency data**, so its model universe is smaller than the MASTL-only network.
        - Correlations are modest in magnitude; this is not a near-duplicate profile.
        - The pathway analysis is a **rank-based enrichment equivalent**, not a classical permutation GSEA implementation.
        - The **2 µM endpoint framing** is deliberate: it avoids over-claiming dose-transition features that the assay spacing cannot support.
        """
    )
    return


if __name__ == "__main__":
    app.run()
