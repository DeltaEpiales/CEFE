"""
cefe_py -- Causal Entropic Field Engine
========================================

A high-performance QFT/QM simulation engine built on C++ with Python bindings.
Computes entanglement entropy, field dynamics, and quantum wave evolution on
causal diamond lattices with support for curved spacetime backgrounds.

Three simulation engines:
  - EntropicFieldEngine:  Entanglement entropy via Srednicki covariance matrices
  - LatticeQFTEngine:     Real scalar lambda-phi^4 field theory (Velocity Verlet)
  - QuantumWaveEngine:    Schrodinger equation via Crank-Nicolson

Five spacetime metrics:
  - MinkowskiMetric:      Flat spacetime
  - SchwarzschildMetric:  Static spherically-symmetric black hole
  - AntiDeSitterMetric:   Anti-de Sitter space (negative cosmological constant)
  - FLRWMetric:           Friedmann-Lemaitre-Robertson-Walker cosmology
  - KerrMetric:           Rotating black hole

Quick Start::

    import cefe_py as ce

    # Build a causal diamond grid in flat spacetime
    grid = ce.CausalDiamondGrid(radius=3.0, spacing=0.5)

    # Compute entanglement entropy of a subregion
    engine = ce.EntropicFieldEngine(grid)
    engine.set_mass(1.0)
    engine.evolve_to_slice(0.0)
    S = engine.compute_entanglement_entropy(1.0)
    print(f"Entanglement entropy: {S:.4f}")
"""

__version__ = "0.2.0"
__author__ = "Ryan Osacra"

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
        SchwarzschildMetric,
        AntiDeSitterMetric,
        FLRWMetric,
        KerrMetric,
        DynamicTensorEngine,
        QuantumWaveEngine,
        LatticeQFTEngine
    )
except ImportError as e:
    raise ImportError(
        f"Failed to load C++ extension: {e}. "
        "Please build the project first with: python -m pip install -e ."
    )


def launch_gui():
    """Launch the interactive CEFE 3D high-fidelity scientific simulator GUI.

    Requires PyQt6 and pyqtgraph. Install with:
        pip install PyQt6 pyqtgraph PyOpenGL
    """
    from .gui import main as gui_main
    gui_main()


def run_cli():
    """Run the CEFE command-line simulator and renderer.

    Usage examples:
        cefe-cli --mode qft --steps 200 --render 3d
        cefe-cli --mode heg --steps 500 --render none
        cefe-cli --mode qm --potential double_slit --render 3d
    """
    from .cli import main as cli_main
    cli_main()


__all__ = [
    # Core types
    "BoundaryType",
    "GridPoint",
    "CausalDiamondGrid",
    # Simulation engines
    "EntropicFieldEngine",
    "LatticeQFTEngine",
    "QuantumWaveEngine",
    "DynamicTensorEngine",  # Experimental
    # Spacetime metrics
    "Metric",
    "MinkowskiMetric",
    "SchwarzschildMetric",
    "AntiDeSitterMetric",
    "FLRWMetric",
    "KerrMetric",
    # Launchers
    "launch_gui",
    "run_cli",
]
