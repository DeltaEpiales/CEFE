#!/usr/bin/env python
"""
engine_backend.py -- unified backend for the entropic-gravity verification suite.

Preference order:
  1. the compiled CEFE engine (cefe_core >= v0.2.0), which exposes the
     Srednicki covariance matrices directly via
         EntropicFieldEngine.get_covariance_C / get_covariance_P
     and the slice-to-grid map via get_slice_indices.
     With the engine present, the paper's numbers come straight from
     cefe_core.
  2. the line-identical NumPy port (entropic_gravity_verification.py),
     used automatically when the extension is not built for the running
     interpreter.

Both paths return (coords, C, P) for the t=0 causal-diamond slice in
Minkowski space. All entropy / modular-Hamiltonian math layered on top is
plain NumPy and identical in both cases -- the same formulas the engine
evaluates internally in compute_entanglement_entropy_indices.
"""

import numpy as np

try:
    import cefe_py as ce
    HAVE_ENGINE = hasattr(ce.EntropicFieldEngine, "get_covariance_C")
except Exception:
    ce = None
    HAVE_ENGINE = False

from entropic_gravity_verification import (
    build_t0_slice, build_laplacian, covariance_matrices,
)

BACKEND = "cefe_core (compiled engine)" if HAVE_ENGINE else "NumPy port (line-identical)"


def vacuum_slice(R_d, a, m, temperature=0.0, return_engine=False):
    """(coords, C, P) of the scalar vacuum/thermal state on the t=0 slice.

    coords : (N,3) slice-point coordinates
    C, P   : (N,N) Srednicki covariance matrices  C = 1/2 W^-1/2, P = 1/2 W^+1/2
             (both with the coth(omega/2T) thermal factor when temperature > 0)
    """
    if HAVE_ENGINE:
        grid = ce.CausalDiamondGrid(R_d, a, ce.MinkowskiMetric())
        eng = ce.EntropicFieldEngine(grid)
        eng.set_mass(m)
        eng.evolve_to_slice(0.0, temperature)
        C = np.asarray(eng.get_covariance_C(), dtype=float)
        P = np.asarray(eng.get_covariance_P(), dtype=float)
        idx = list(eng.get_slice_indices())
        pts = grid.get_points()
        coords = np.array([[pts[i].x, pts[i].y, pts[i].z] for i in idx])
        return (coords, C, P, eng) if return_engine else (coords, C, P)
    pts = build_t0_slice(R_d, a)
    coords = np.array([[p[3], p[4], p[5]] for p in pts])
    C, P = covariance_matrices(build_laplacian(pts, a), m, temperature=temperature)
    return (coords, C, P, None) if return_engine else (coords, C, P)


def engine_crosscheck(R_d=4.0, a=0.8, m=1.0, radii=None):
    """Three-way consistency check (engine build only):
       (i)   engine's internal compute_entanglement_entropy(R)
       (ii)  symplectic entropy from the engine's exported C, P matrices
       (iii) the NumPy port, same lattice parameters
    Prints the maximum relative deviations."""
    if not HAVE_ENGINE:
        return
    from entropic_gravity_verification import entanglement_entropy
    if radii is None:
        radii = np.linspace(1.2, 3.2, 5)
    coords, C, P, eng = vacuum_slice(R_d, a, m, return_engine=True)
    S_mat = np.array([entanglement_entropy(C, P, coords, r) for r in radii])
    S_int = np.array([eng.compute_entanglement_entropy(r) for r in radii])
    pts_p = build_t0_slice(R_d, a)
    coords_p = np.array([[p[3], p[4], p[5]] for p in pts_p])
    C_p, P_p = covariance_matrices(build_laplacian(pts_p, a), m)
    S_port = np.array([entanglement_entropy(C_p, P_p, coords_p, r) for r in radii])
    d_int = np.max(np.abs(S_mat - S_int) / np.maximum(S_int, 1e-12))
    d_port = np.max(np.abs(S_mat - S_port) / np.maximum(S_port, 1e-12))
    print(f"    cross-check (a={a}, m={m}): engine internal vs exported matrices: "
          f"max rel dev {d_int:.2e}")
    print(f"    cross-check (a={a}, m={m}): engine vs NumPy port:               "
          f"max rel dev {d_port:.2e}")
