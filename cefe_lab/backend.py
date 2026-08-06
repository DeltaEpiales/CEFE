"""
cefe_lab.backend -- engine / NumPy backend selection.

The compiled CEFE engine (cefe_core >= 0.2.0) exposes the Srednicki
covariance matrices of a causal-diamond slice directly.  When it is not
available (extension not built for the running interpreter), cefe_lab uses
an identical NumPy construction.  All Gaussian-state math (symplectic
spectra, entropies, modular Hamiltonians) is NumPy in both cases -- the same
formulas the engine evaluates internally.
"""

import numpy as np

try:
    import cefe_py as _ce
    HAVE_ENGINE = hasattr(_ce.EntropicFieldEngine, "get_covariance_C")
except Exception:
    _ce = None
    HAVE_ENGINE = False

NAME = "cefe_core (compiled engine)" if HAVE_ENGINE else "NumPy (built-in)"


def engine_vacuum(radius, spacing, mass, temperature=0.0):
    """Vacuum/thermal covariance state from the compiled engine.

    Returns (coords, C, P) for the t=0 causal-diamond slice in Minkowski
    space.  Raises RuntimeError if the engine is unavailable.
    """
    if not HAVE_ENGINE:
        raise RuntimeError("cefe_core is not importable in this interpreter")
    grid = _ce.CausalDiamondGrid(radius, spacing, _ce.MinkowskiMetric())
    eng = _ce.EntropicFieldEngine(grid)
    eng.set_mass(mass)
    eng.evolve_to_slice(0.0, temperature)
    C = np.asarray(eng.get_covariance_C(), dtype=float)
    P = np.asarray(eng.get_covariance_P(), dtype=float)
    idx = list(eng.get_slice_indices())
    pts = grid.get_points()
    coords = np.array([[p.x, p.y, p.z] for p in (pts[i] for i in idx)])
    return coords, C, P
