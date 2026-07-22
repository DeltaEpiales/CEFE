import os
import sys

# Add the directory containing the compiled extension to the path
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

try:
    from .cefe_core import (
        BoundaryType,
        GridPoint,
        CausalDiamondGrid,
        EntropicFieldEngine,
        Metric,
        MinkowskiMetric,
        SchwarzschildMetric
    )
except ImportError as e:
    raise ImportError(f"Failed to load C++ extension: {e}. Please build the project first.")

__all__ = [
    "BoundaryType",
    "GridPoint",
    "CausalDiamondGrid",
    "EntropicFieldEngine",
    "Metric",
    "MinkowskiMetric",
    "SchwarzschildMetric"
]
