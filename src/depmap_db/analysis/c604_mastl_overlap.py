from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import polars as pl
import requests
from scipy import stats

from depmap_db.polars import prepare_lazy_tables

SCREEN_ID = "MTS028_BIAS_CORRECTED"
COMPOUND_NAME = "c-604"
TOP_DOSE_UM = 2.0
MIN_MODELS_PER_GENE = 8
MIN_PATHWAY_SIZE = 10
MAX_PATHWAY_SIZE = 500
POSITIVE_FDR_CUTOFF = 0.1
PATHWAY_LIBRARIES = ("MSigDB_Hallmark_2020", "Reactome_2022")


def benjamini_hochberg(p_values: pd.Series | np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    n = len(p)
    order = np.argsort(np.nan_to_num(p, nan=np.inf))
    ranked = p[order]
    q = np.empty(n, dtype=float)
    running = 1.0
    for i in range(n - 1, -1, -1):
        value = ranked[i] * n / (i + 1)
        running = min(running, value)
        q[i] = running
    out = np.empty(n, dtype=float)
    out[order] = np.clip(q, 0.0, 1.0)
    return out


def _fetch_gene_set_library(library_name: str, cache_dir: Path) -> dict[str, set[str]]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / f"{library_name}.txt"
    if not cache_path.exists():
        url = (
            "https://maayanlab.cloud/Enrichr/geneSetLibrary"
            f"?mode=text&libraryName={library_name}"
        )
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        cache_path.write_text(response.text)

    gene_sets: dict[str, set[str]] = {}
    for line in cache_path.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        term = parts[0].strip()
        genes = {gene.strip() for gene in parts[2:] if gene.strip()}
        if term and genes:
            gene_sets[term] = genes
    return gene_sets


def _load_inputs(
    *,
    project_root: Path,
    db_path: Path,
    polars_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    tables = prepare_lazy_tables(
        output_dir=polars_dir,
        db_path=db_path,
        tables=[
            "models",
            "drug_response_secondary",
            "drug_response_secondary_dose",
            "gene_effects_wide",
        ],
    )

    model_cols = [
        "model_id",
        "cell_line_name",
        "oncotree_lineage",
        "oncotree_primary_disease",
        "oncotree_subtype",
    ]

    c604_endpoint = (
        tables["drug_response_secondary_dose"]
        .join(
            tables["drug_response_secondary"].select(
                ["response_id", "compound_name", "target_text"]
            ),
            on="response_id",
            how="left",
        )
        .filter(
            (pl.col("screen_id") == SCREEN_ID)
            & (pl.col("compound_name").str.to_lowercase() == COMPOUND_NAME)
            & (pl.col("dose_um") == TOP_DOSE_UM)
        )
        .join(tables["models"].select(model_cols), on="model_id", how="left")
        .with_columns(
            (2 ** pl.col("median_l2fc")).alias("growth_fraction"),
            (-(2 ** pl.col("median_l2fc")).alias("c604_sensitivity")),
        )
        .select(
            [
                "model_id",
                "cell_line_name",
                "oncotree_lineage",
                "oncotree_primary_disease",
                "oncotree_subtype",
                "median_l2fc",
                "growth_fraction",
                "c604_sensitivity",
            ]
        )
        .collect()
        .to_pandas()
    )

    dependency = tables["gene_effects_wide"].collect().to_pandas()
    numeric_gene_cols = [
        col
        for col in dependency.columns
        if col != "model_id" and pd.api.types.is_numeric_dtype(dependency[col])
    ]
    dependency.loc[:, numeric_gene_cols] = -dependency.loc[:, numeric_gene_cols]
    return c604_endpoint, dependency[["model_id", *numeric_gene_cols]].copy()


def _correlate_profile(
    phenotype: pd.Series,
    feature_matrix: pd.DataFrame,
    *,
    exclude: set[str] | None = None,
    min_n: int = MIN_MODELS_PER_GENE,
    profile_name: str,
) -> pd.DataFrame:
    exclude = exclude or set()
    feature_cols = [
        col for col in feature_matrix.columns if col != "model_id" and col not in exclude
    ]
    X = feature_matrix.loc[:, feature_cols]
    r_values = X.corrwith(phenotype)

    rows: list[dict[str, Any]] = []
    phenotype_notna = phenotype.notna()
    for gene in feature_cols:
        gene_values = X[gene]
        mask = phenotype_notna & gene_values.notna()
        n_models = int(mask.sum())
        if n_models < min_n:
            continue
        r = float(r_values[gene])
        if math.isnan(r):
            continue
        t_stat = r * math.sqrt((n_models - 2) / max(1e-12, 1.0 - r * r))
        p_value = 2 * stats.t.sf(abs(t_stat), df=n_models - 2)
        rows.append(
            {
                "gene": gene,
                "profile": profile_name,
                "n_models": n_models,
                "pearson_r": r,
                "p_value": float(p_value),
            }
        )

    correlations = pd.DataFrame(rows).sort_values(
        "pearson_r", ascending=False
    ).reset_index(drop=True)
    correlations["fdr"] = benjamini_hochberg(correlations["p_value"].to_numpy())
    correlations["rank_positive"] = np.arange(1, len(correlations) + 1)
    return correlations


def _rank_pathways(
    ranked_genes: pd.DataFrame,
    gene_sets: dict[str, set[str]],
    *,
    min_size: int = MIN_PATHWAY_SIZE,
    max_size: int = MAX_PATHWAY_SIZE,
    library_name: str,
    profile_name: str,
) -> pd.DataFrame:
    ordered = ranked_genes.loc[:, ["gene", "pearson_r"]].reset_index(drop=True)
    gene_to_index = {gene: i for i, gene in enumerate(ordered["gene"])}
    values = ordered["pearson_r"].to_numpy()

    rows: list[dict[str, Any]] = []
    for term, genes in gene_sets.items():
        hit_indices = [gene_to_index[gene] for gene in genes if gene in gene_to_index]
        overlap_size = len(hit_indices)
        if overlap_size < min_size or overlap_size > max_size:
            continue

        in_set = values[hit_indices]
        out_mask = np.ones(len(values), dtype=bool)
        out_mask[hit_indices] = False
        out_set = values[out_mask]

        mw = stats.mannwhitneyu(in_set, out_set, alternative="greater")
        auc = mw.statistic / (len(in_set) * len(out_set))
        rows.append(
            {
                "library": library_name,
                "profile": profile_name,
                "term": term,
                "overlap_size": overlap_size,
                "mean_r": float(np.mean(in_set)),
                "median_r": float(np.median(in_set)),
                "auc": float(auc),
                "p_value": float(mw.pvalue),
            }
        )

    enriched = pd.DataFrame(rows).sort_values(
        ["p_value", "mean_r"], ascending=[True, False]
    )
    enriched["fdr"] = benjamini_hochberg(enriched["p_value"].to_numpy())
    enriched["rank_positive"] = np.arange(1, len(enriched) + 1)
    return enriched


def _hypergeom_overlap(
    *,
    universe_size: int,
    set_a_size: int,
    set_b_size: int,
    overlap_size: int,
) -> float:
    if overlap_size <= 0:
        return 1.0
    distribution = stats.hypergeom(M=universe_size, n=set_a_size, N=set_b_size)
    return float(distribution.sf(overlap_size - 1))


def _pathway_concordance(
    pathway_a: pd.DataFrame,
    pathway_b: pd.DataFrame,
    *,
    fdr_cutoff: float = POSITIVE_FDR_CUTOFF,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    merged = pathway_a.merge(
        pathway_b,
        on=["library", "term"],
        suffixes=("_c604", "_mastl"),
    )
    shared_sig = merged.loc[
        (merged["fdr_c604"] < fdr_cutoff) & (merged["fdr_mastl"] < fdr_cutoff)
    ].copy()
    shared_sig = shared_sig.sort_values(
        ["mean_r_c604", "mean_r_mastl"], ascending=False
    ).reset_index(drop=True)

    summary_rows: list[dict[str, Any]] = []
    for library, sub in merged.groupby("library", sort=True):
        sig_a = set(sub.loc[sub["fdr_c604"] < fdr_cutoff, "term"])
        sig_b = set(sub.loc[sub["fdr_mastl"] < fdr_cutoff, "term"])
        overlap = sig_a & sig_b
        spearman = stats.spearmanr(sub["mean_r_c604"], sub["mean_r_mastl"])
        summary_rows.append(
            {
                "library": library,
                "n_terms": int(len(sub)),
                "c604_sig_terms": int(len(sig_a)),
                "mastl_sig_terms": int(len(sig_b)),
                "shared_sig_terms": int(len(overlap)),
                "shared_sig_p_value": _hypergeom_overlap(
                    universe_size=len(sub),
                    set_a_size=len(sig_a),
                    set_b_size=len(sig_b),
                    overlap_size=len(overlap),
                ),
                "spearman_rho": float(spearman.statistic),
                "spearman_p_value": float(spearman.pvalue),
            }
        )

    summary = pd.DataFrame(summary_rows).sort_values("library").reset_index(drop=True)
    return summary, shared_sig


def run_analysis(
    *,
    project_root: str | Path,
    db_path: str | Path | None = None,
    polars_dir: str | Path | None = None,
    pathway_cache_dir: str | Path | None = None,
    report_dir: str | Path | None = None,
) -> dict[str, Any]:
    project_root = Path(project_root)
    db_path = Path(db_path) if db_path is not None else Path.home() / ".depmap" / "depmap.duckdb"
    polars_dir = Path(polars_dir) if polars_dir is not None else project_root / "data" / "polars"
    pathway_cache_dir = (
        Path(pathway_cache_dir)
        if pathway_cache_dir is not None
        else project_root / "data" / "pathways"
    )
    report_dir = (
        Path(report_dir)
        if report_dir is not None
        else project_root / "reports" / "mts028_c604_mastl_overlap"
    )
    report_dir.mkdir(parents=True, exist_ok=True)

    c604_endpoint, dependency_strength = _load_inputs(
        project_root=project_root,
        db_path=db_path,
        polars_dir=polars_dir,
    )

    c604_dependency = c604_endpoint.merge(dependency_strength, on="model_id", how="inner")
    c604_phenotype = c604_dependency["c604_sensitivity"]
    c604_gene_correlations = _correlate_profile(
        c604_phenotype,
        c604_dependency.drop(columns=[
            "cell_line_name",
            "oncotree_lineage",
            "oncotree_primary_disease",
            "oncotree_subtype",
            "median_l2fc",
            "growth_fraction",
            "c604_sensitivity",
        ]),
        profile_name="C-604 2 uM endpoint vs dependency strength",
    )

    mastl_gene_correlations = _correlate_profile(
        dependency_strength["MASTL"],
        dependency_strength,
        exclude={"MASTL"},
        profile_name="MASTL dependency strength vs dependency strength",
    )

    libraries = {
        library_name: _fetch_gene_set_library(library_name, pathway_cache_dir)
        for library_name in PATHWAY_LIBRARIES
    }
    c604_pathways = pd.concat(
        [
            _rank_pathways(
                c604_gene_correlations,
                gene_sets,
                library_name=library_name,
                profile_name="C-604 2 uM endpoint",
            )
            for library_name, gene_sets in libraries.items()
        ],
        ignore_index=True,
    )
    mastl_pathways = pd.concat(
        [
            _rank_pathways(
                mastl_gene_correlations,
                gene_sets,
                library_name=library_name,
                profile_name="MASTL dependency profile",
            )
            for library_name, gene_sets in libraries.items()
        ],
        ignore_index=True,
    )

    top100_c604 = c604_gene_correlations.head(100).copy()
    top100_mastl = mastl_gene_correlations.head(100).copy()
    overlap_genes = sorted(set(top100_c604["gene"]) & set(top100_mastl["gene"]))
    overlap_table = (
        top100_c604.loc[top100_c604["gene"].isin(overlap_genes)]
        .merge(
            top100_mastl.loc[top100_mastl["gene"].isin(overlap_genes)],
            on="gene",
            suffixes=("_c604", "_mastl"),
        )
        .sort_values(["rank_positive_c604", "rank_positive_mastl"])
        .reset_index(drop=True)
    )

    gene_overlap_summary = pd.DataFrame(
        [
            {
                "gene_universe": int(
                    len(set(c604_gene_correlations["gene"]) & set(mastl_gene_correlations["gene"]))
                ),
                "c604_top_n": 100,
                "mastl_top_n": 100,
                "shared_top_genes": int(len(overlap_genes)),
                "shared_top_gene_p_value": _hypergeom_overlap(
                    universe_size=len(
                        set(c604_gene_correlations["gene"]) & set(mastl_gene_correlations["gene"])
                    ),
                    set_a_size=100,
                    set_b_size=100,
                    overlap_size=len(overlap_genes),
                ),
            }
        ]
    )

    pathway_concordance_summary, shared_pathways = _pathway_concordance(
        c604_pathways,
        mastl_pathways,
    )

    analysis_summary = pd.DataFrame(
        [
            {
                "c604_models_with_dependency": int(c604_dependency["model_id"].nunique()),
                "c604_screen_models": int(c604_endpoint["model_id"].nunique()),
                "mastl_models_with_dependency": int(dependency_strength["model_id"].nunique()),
                "gene_universe": int(
                    len(set(c604_gene_correlations["gene"]) & set(mastl_gene_correlations["gene"]))
                ),
                "top_c604_gene": str(top100_c604.iloc[0]["gene"]),
                "top_c604_r": float(top100_c604.iloc[0]["pearson_r"]),
                "top_mastl_gene": str(top100_mastl.iloc[0]["gene"]),
                "top_mastl_r": float(top100_mastl.iloc[0]["pearson_r"]),
            }
        ]
    )

    outputs = {
        "analysis_summary": analysis_summary,
        "c604_endpoint": c604_endpoint,
        "top100_c604": top100_c604,
        "top100_mastl": top100_mastl,
        "c604_gene_correlations": c604_gene_correlations,
        "mastl_gene_correlations": mastl_gene_correlations,
        "overlap_table": overlap_table,
        "gene_overlap_summary": gene_overlap_summary,
        "c604_pathways": c604_pathways,
        "mastl_pathways": mastl_pathways,
        "pathway_concordance_summary": pathway_concordance_summary,
        "shared_pathways": shared_pathways,
    }

    for name, frame in outputs.items():
        if isinstance(frame, pd.DataFrame):
            frame.to_csv(report_dir / f"{name}.csv", index=False)

    return outputs
