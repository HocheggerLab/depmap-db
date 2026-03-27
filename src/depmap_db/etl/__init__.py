"""ETL (Extract, Transform, Load) pipeline for DepMap data."""

from .processors import BaseProcessor, ProcessingResult

__all__ = ["BaseProcessor", "ProcessingResult"]
