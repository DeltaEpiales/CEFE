"""
cefe_lab -- theory verification laboratory for the CEFE engine.

A composable API for testing quantum / QFT / emergent-gravity theories
against exact lattice computations.  Write a script, build a state, measure
an observable, compare against your theory's prediction.

Backend selection is automatic:
  * if the compiled CEFE engine (cefe_core >= 0.2.0) is importable, vacuum
    states are produced by the C++ engine (EntropicFieldEngine covariance
    export);
  * otherwise an identical NumPy construction is used (same Srednicki
    covariance-matrix method, same graph-Laplacian discretization).

Example
-------
>>> import cefe_lab as lab
>>> vac = lab.VacuumState.disk(radius=4.0, spacing=0.2, mass=1.0)
>>> S = [vac.entropy(lab.ball_region(vac.coords, r)) for r in (1.0, 1.5, 2.0)]
"""

from .lattice import (
    disk_lattice, ball_lattice, graph_laplacian,
    schwarzschild_proper_distance, Lattice,
)
from .gaussian import (
    VacuumState, symplectic_eigenvalues, entropy_from_nu,
    ball_region, annulus_region, mutual_information,
)
from .screens import (
    entropy_profile, verlinde_force, power_law_fit, g_eff_per_species,
)
from .tests import TheoryTest, TestResult
from . import backend

__version__ = "0.1.0"
__all__ = [
    "disk_lattice", "ball_lattice", "graph_laplacian",
    "schwarzschild_proper_distance", "Lattice",
    "VacuumState", "symplectic_eigenvalues", "entropy_from_nu",
    "ball_region", "annulus_region", "mutual_information",
    "entropy_profile", "verlinde_force", "power_law_fit", "g_eff_per_species",
    "TheoryTest", "TestResult", "backend",
]
