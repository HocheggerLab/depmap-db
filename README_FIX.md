# README_FIX.md

Paste the section below over the current `## Figure API pattern` block in `README.md`.

---

## Figure API pattern

Most plotting helpers in `depmap_db.plots` follow a three-layer design:

1. `analyse_*()` → returns a DataFrame used for plotting
2. `plot_*()` → plots a supplied DataFrame
3. a convenience wrapper (for example `expr_vs_dep()` or `lineage_analysis()`) → does both in one step and returns a Matplotlib figure

### Choosing your workflow

Depending on whether you want convenience or performance, you should choose one of two patterns:

#### 1. The "Convenience" Pattern (Quick & Easy)
Best for generating a single plot quickly. It handles the database connection, querying, and plotting in one call.

```python
import depmap_db.plots as plots

fig = plots.lineage_analysis("MASTL", assay="dependency", db_path=DB)
fig.savefig("lineage.png")
```

#### 2. The "Efficient" Pattern (Best for multiple plots)
Best if you are building a multi-panel figure. You load the data into memory **once**, then pass it to multiple plotting functions. This is much faster for complex workflows.

```python
from depmap_db import scan_gene_effects_wide, scan_gene_expression_wide
from depmap_db.plots import plot_lineage, plot_expr_vs_dep

# 1. Load the data into memory once
dep_df = scan_gene_effects_wide(db_path=DB).collect()
expr_df = scan_gene_expression_wide(db_path=DB).collect()

# 2. Pass the pre-loaded DataFrames to plotting functions
fig1, ax1 = plot_lineage(dep_df, gene="MASTL", assay="dependency")
fig2, ax2 = plot_expr_vs_dep(expr_df, gene="MASTL")
```

That means you can either:

- do everything in one line, or
- split analysis and plotting when you want to inspect or modify the intermediate data

For some plot families, the convenience function is what is exported at the top level; if you want the lower-level analysis helper, import it from the specific module inside `depmap_db.plots`.
