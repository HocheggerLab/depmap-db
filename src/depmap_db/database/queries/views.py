"""Predefined views for common DepMap queries."""

from dataclasses import dataclass

from ...config import get_logger
from ..connection import get_db_manager

logger = get_logger(__name__)


@dataclass
class ViewDefinition:
    """Definition for a database view."""

    name: str
    sql: str
    description: str
    dependencies: list[str]


class ViewManager:
    """Manages database views for common queries."""

    def __init__(self) -> None:
        self.db_manager = get_db_manager()
        self._views: dict[str, ViewDefinition] = {}
        self._define_views()

    def _define_views(self) -> None:
        """Define all predefined views."""

        # Top dependencies across all models
        self._add_view(
            ViewDefinition(
                name="top_dependencies",
                dependencies=["gene_effects", "genes"],
                description="Genes with strongest average dependency across all models",
                sql="""
            CREATE OR REPLACE VIEW top_dependencies AS
            SELECT
                g.hugo_symbol,
                g.gene_id,
                COUNT(*) as model_count,
                AVG(ge.gene_effect) as avg_gene_effect,
                STDDEV(ge.gene_effect) as stddev_gene_effect,
                MIN(ge.gene_effect) as min_gene_effect,
                MAX(ge.gene_effect) as max_gene_effect,
                PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY ge.gene_effect) as q1_gene_effect,
                PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY ge.gene_effect) as median_gene_effect,
                PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY ge.gene_effect) as q3_gene_effect
            FROM gene_effects ge
            JOIN genes g ON ge.gene_id = g.gene_id
            WHERE ge.gene_effect < -0.5  -- Focus on dependencies
            GROUP BY g.hugo_symbol, g.gene_id
            HAVING COUNT(*) >= 10  -- At least 10 models
            ORDER BY avg_gene_effect ASC
            """,
            )
        )

        # Model vulnerabilities - strongest dependencies per model
        self._add_view(
            ViewDefinition(
                name="model_vulnerabilities",
                dependencies=["gene_effects", "genes", "models"],
                description="Strongest genetic dependencies for each cancer model",
                sql="""
            CREATE OR REPLACE VIEW model_vulnerabilities AS
            SELECT
                m.model_id,
                m.cell_line_name,
                m.oncotree_lineage,
                m.oncotree_primary_disease,
                g.hugo_symbol,
                g.gene_id,
                ge.gene_effect,
                RANK() OVER (PARTITION BY m.model_id ORDER BY ge.gene_effect ASC) as dependency_rank
            FROM gene_effects ge
            JOIN models m ON ge.model_id = m.model_id
            JOIN genes g ON ge.gene_id = g.gene_id
            WHERE ge.gene_effect < -1.0  -- Strong dependencies only
            QUALIFY dependency_rank <= 20  -- Top 20 dependencies per model
            """,
            )
        )

        # Gene selectivity - context-specific vs pan-essential
        self._add_view(
            ViewDefinition(
                name="gene_selectivity",
                dependencies=["gene_effects", "genes", "models"],
                description="Analysis of gene dependency selectivity across cancer contexts",
                sql="""
            CREATE OR REPLACE VIEW gene_selectivity AS
            WITH gene_stats AS (
                SELECT
                    ge.gene_id,
                    g.hugo_symbol,
                    COUNT(*) as total_models,
                    COUNT(CASE WHEN ge.gene_effect < -0.5 THEN 1 END) as dependent_models,
                    COUNT(CASE WHEN ge.gene_effect < -1.0 THEN 1 END) as strongly_dependent_models,
                    AVG(ge.gene_effect) as avg_gene_effect,
                    STDDEV(ge.gene_effect) as stddev_gene_effect
                FROM gene_effects ge
                JOIN genes g ON ge.gene_id = g.gene_id
                GROUP BY ge.gene_id, g.hugo_symbol
            ),
            lineage_stats AS (
                SELECT
                    ge.gene_id,
                    m.oncotree_lineage,
                    COUNT(*) as lineage_models,
                    COUNT(CASE WHEN ge.gene_effect < -0.5 THEN 1 END) as lineage_dependent_models,
                    AVG(ge.gene_effect) as lineage_avg_effect
                FROM gene_effects ge
                JOIN models m ON ge.model_id = m.model_id
                WHERE m.oncotree_lineage IS NOT NULL
                GROUP BY ge.gene_id, m.oncotree_lineage
            )
            SELECT
                gs.gene_id,
                gs.hugo_symbol,
                gs.total_models,
                gs.dependent_models,
                gs.strongly_dependent_models,
                gs.dependent_models::FLOAT / gs.total_models as dependency_frequency,
                gs.avg_gene_effect,
                gs.stddev_gene_effect,
                COUNT(DISTINCT ls.oncotree_lineage) as lineages_tested,
                COUNT(CASE WHEN ls.lineage_dependent_models::FLOAT / ls.lineage_models > 0.8 THEN 1 END) as lineages_pan_essential,
                CASE
                    WHEN gs.dependent_models::FLOAT / gs.total_models > 0.8 THEN 'Pan-essential'
                    WHEN gs.dependent_models::FLOAT / gs.total_models < 0.2 THEN 'Non-essential'
                    ELSE 'Context-specific'
                END as essentiality_class
            FROM gene_stats gs
            LEFT JOIN lineage_stats ls ON gs.gene_id = ls.gene_id
            WHERE gs.total_models >= 10
            GROUP BY gs.gene_id, gs.hugo_symbol, gs.total_models, gs.dependent_models,
                     gs.strongly_dependent_models, gs.avg_gene_effect, gs.stddev_gene_effect
            ORDER BY dependency_frequency DESC
            """,
            )
        )

        # Quality metrics summary
        self._add_view(
            ViewDefinition(
                name="quality_metrics",
                dependencies=["models", "screens", "model_conditions"],
                description="Summary of screen and model quality metrics",
                sql="""
            CREATE OR REPLACE VIEW quality_metrics AS
            SELECT
                m.model_id,
                m.cell_line_name,
                m.oncotree_lineage,
                m.oncotree_primary_disease,
                COUNT(s.screen_id) as screen_count,
                COUNT(CASE WHEN s.passes_qc THEN 1 END) as passing_screens,
                AVG(s.screen_nnmd) as avg_nnmd,
                AVG(s.screen_roc_auc) as avg_roc_auc,
                AVG(s.screen_fpr) as avg_fpr,
                STRING_AGG(DISTINCT s.library, ', ') as libraries_used,
                MIN(s.days) as min_days,
                MAX(s.days) as max_days
            FROM models m
            LEFT JOIN screens s ON m.model_id = s.model_id
            GROUP BY m.model_id, m.cell_line_name, m.oncotree_lineage, m.oncotree_primary_disease
            ORDER BY passing_screens DESC, avg_roc_auc DESC
            """,
            )
        )

        # Lineage summary statistics
        self._add_view(
            ViewDefinition(
                name="lineage_summary",
                dependencies=["models", "gene_effects"],
                description="Summary statistics by cancer lineage",
                sql="""
            CREATE OR REPLACE VIEW lineage_summary AS
            SELECT
                m.oncotree_lineage,
                COUNT(DISTINCT m.model_id) as model_count,
                COUNT(DISTINCT m.oncotree_primary_disease) as disease_count,
                COUNT(DISTINCT ge.gene_id) as genes_tested,
                AVG(ge.gene_effect) as avg_gene_effect,
                COUNT(CASE WHEN ge.gene_effect < -0.5 THEN 1 END) as dependency_events,
                COUNT(CASE WHEN ge.gene_effect < -1.0 THEN 1 END) as strong_dependency_events
            FROM models m
            LEFT JOIN gene_effects ge ON m.model_id = ge.model_id
            WHERE m.oncotree_lineage IS NOT NULL
            GROUP BY m.oncotree_lineage
            HAVING COUNT(DISTINCT m.model_id) >= 5  -- At least 5 models
            ORDER BY model_count DESC
            """,
            )
        )

    def _add_view(self, view_definition: ViewDefinition) -> None:
        """Add a view definition to the registry."""
        self._views[view_definition.name] = view_definition

    def get_view_definition(self, view_name: str) -> ViewDefinition:
        """Get definition for a specific view."""
        if view_name not in self._views:
            raise ValueError(f"Unknown view: {view_name}")
        return self._views[view_name]

    def list_views(self) -> list[str]:
        """Get list of available view names."""
        return list(self._views.keys())

    def create_view(self, view_name: str) -> None:
        """Create a specific view."""
        if view_name not in self._views:
            raise ValueError(f"Unknown view: {view_name}")

        view_def = self._views[view_name]
        logger.info("Creating view: %s", view_name)

        try:
            self.db_manager.execute(view_def.sql)
            logger.info("Successfully created view: %s", view_name)
        except (OSError, RuntimeError) as e:
            logger.error("Failed to create view %s: %s", view_name, e)
            raise

    def create_all_views(self) -> None:
        """Create all predefined views."""
        logger.info("Creating all views: %s", list(self._views.keys()))

        for view_name in self._views:
            self.create_view(view_name)

        logger.info("All views created successfully")

    def drop_view(self, view_name: str) -> None:
        """Drop a specific view."""
        try:
            self.db_manager.execute(f"DROP VIEW IF EXISTS {view_name}")
            logger.info("Dropped view: %s", view_name)
        except (OSError, RuntimeError) as e:
            logger.error("Failed to drop view %s: %s", view_name, e)
            raise


def create_all_views() -> None:
    """Create all predefined views."""
    view_manager = ViewManager()
    view_manager.create_all_views()


def get_available_views() -> list[str]:
    """Get list of available predefined views."""
    view_manager = ViewManager()
    return view_manager.list_views()


def create_view(view_name: str) -> None:
    """Create a specific predefined view."""
    view_manager = ViewManager()
    view_manager.create_view(view_name)
