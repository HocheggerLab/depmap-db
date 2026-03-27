"""Data processors for different DepMap file types."""

from .base import BaseProcessor, ProcessingResult
from .gene import GeneProcessor
from .gene_effect import GeneEffectProcessor
from .gene_effect_wide import GeneEffectWideProcessor
from .gene_expression import GeneExpressionProcessor
from .gene_expression_wide import GeneExpressionWideProcessor
from .model import ModelProcessor

__all__ = [
    "BaseProcessor",
    "ProcessingResult",
    "GeneEffectProcessor",
    "GeneEffectWideProcessor",
    "GeneExpressionProcessor",
    "GeneExpressionWideProcessor",
    "GeneProcessor",
    "ModelProcessor",
]
