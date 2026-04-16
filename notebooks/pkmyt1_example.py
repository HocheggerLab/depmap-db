"""Notebook-style PKMYT1 example for depmap-db.

Run with:
    uv run python notebooks/pkmyt1_example.py

This generates example figure assets for the README.
"""

from __future__ import annotations

from pathlib import Path

import depmap_db.plots as plots

DB_PATH = Path.home() / ".depmap" / "depmap.duckdb"
OUTPUT_DIR = Path(__file__).resolve().parents[1] / "reports" / "readme_assets"
GENE = "PKMYT1"


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    fig = plots.lineage_analysis(GENE, assay="dependency", db_path=DB_PATH)
    fig.savefig(
        OUTPUT_DIR / "pkmyt1_lineage_dependency.png",
        dpi=300,
        bbox_inches="tight",
    )

    fig = plots.expr_vs_dep(GENE, db_path=DB_PATH)
    fig.savefig(
        OUTPUT_DIR / "pkmyt1_expr_vs_dep.png",
        dpi=300,
        bbox_inches="tight",
    )

    print(f"Saved example PKMYT1 figures to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
