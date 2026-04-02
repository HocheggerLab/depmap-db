"""DepMap plotting functions."""

from .colors import COLOR
from .drug_plots import (
    analyse_drug_by_lineage,
    analyse_drug_by_mutation,
    analyse_drug_dependency,
    drug_by_lineage,
    drug_by_mutation,
    drug_dependency_scatter,
    plot_drug_by_lineage,
    plot_drug_by_mutation,
    plot_drug_dependency,
)
from .expression_dependency import (
    analyse_biomarker_volcano,
    analyse_expr_vs_dep,
    biomarker_volcano,
    expr_vs_dep,
    plot_biomarker_volcano,
    plot_expr_vs_dep,
)
from .gene_plots import (
    analyse_gene_gene,
    analyse_top_correlations,
    gene_gene_scatter,
    plot_gene_gene,
    plot_top_correlations,
    top_correlations,
)
from .heatmap import (
    analyse_gene_heatmap,
    analyse_selectivity,
    gene_heatmap,
    plot_gene_heatmap,
    plot_selectivity,
    selectivity_plot,
)
from .lineage_plots import lineage_analysis, plot_lineage
from .mutation_comparison import (
    analyse_by_mutation,
    mutation_analysis,
    plot_mutation_comparison,
)
from .mutation_plots import (
    analyse_comutation,
    analyse_mutation_waterfall,
    comutation_analysis,
    mutation_waterfall,
    plot_comutation,
    plot_mutation_waterfall,
)

__all__ = [
    # colors
    "COLOR",
    # lineage (existing)
    "lineage_analysis",
    "plot_lineage",
    # mutation comparison (existing)
    "analyse_by_mutation",
    "mutation_analysis",
    "plot_mutation_comparison",
    # gene-gene scatter (1)
    "analyse_gene_gene",
    "gene_gene_scatter",
    "plot_gene_gene",
    # top correlations (2)
    "analyse_top_correlations",
    "plot_top_correlations",
    "top_correlations",
    # expression vs dependency (3)
    "analyse_expr_vs_dep",
    "expr_vs_dep",
    "plot_expr_vs_dep",
    # biomarker volcano (4)
    "analyse_biomarker_volcano",
    "biomarker_volcano",
    "plot_biomarker_volcano",
    # drug by lineage (5)
    "analyse_drug_by_lineage",
    "drug_by_lineage",
    "plot_drug_by_lineage",
    # drug by mutation (6)
    "analyse_drug_by_mutation",
    "drug_by_mutation",
    "plot_drug_by_mutation",
    # drug-dependency correlation (7)
    "analyse_drug_dependency",
    "drug_dependency_scatter",
    "plot_drug_dependency",
    # mutation waterfall (8)
    "analyse_mutation_waterfall",
    "mutation_waterfall",
    "plot_mutation_waterfall",
    # co-mutation enrichment (9)
    "analyse_comutation",
    "comutation_analysis",
    "plot_comutation",
    # selectivity waterfall (10)
    "analyse_selectivity",
    "plot_selectivity",
    "selectivity_plot",
    # multi-gene heatmap (11)
    "analyse_gene_heatmap",
    "gene_heatmap",
    "plot_gene_heatmap",
]
