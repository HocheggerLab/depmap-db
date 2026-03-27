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
        logger.info("No predefined views configured for wide-only schema")

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
