"""
cefe_lab.gaussian -- Gaussian (vacuum/thermal) states and their entropies.

A VacuumState bundles the t=0-slice coordinates with the Srednicki
covariance pair
    C = 1/2 W^{-1/2} ,  P = 1/2 W^{+1/2} ,   W = L + diag(m_i^2)
(with the coth(omega/2T) thermal factor at T > 0).  A position-dependent
mass profile m_i models localized "matter": a heavy region is a test mass.

All entropy math uses the numerically stable Williamson form:
symplectic eigenvalues are sqrt(eig(C^{1/2} P C^{1/2})).
"""

import numpy as np

from . import backend
from .lattice import Lattice, graph_laplacian


def _covariances_from_W(W, temperature=0.0):
    w2, U = np.linalg.eigh(W)
    w2 = np.clip(w2, 1e-30, None)
    w = np.sqrt(w2)
    if temperature > 0.0:
        f = 1.0 / np.tanh(w / (2.0 * temperature))   # coth
    else:
        f = np.ones_like(w)
    C = (U * (0.5 * f / w)[None, :]) @ U.T
    P = (U * (0.5 * f * w)[None, :]) @ U.T
    return C, P


class VacuumState:
    """Gaussian state of a free scalar on a lattice slice."""
    def __init__(self, coords, C, P, mass, temperature=0.0, source=""):
        self.coords = coords
        self.C = C
        self.P = P
        self.mass = mass
        self.temperature = temperature
        self.source = source

    # -- constructors ------------------------------------------------------
    @classmethod
    def from_lattice(cls, lat: Lattice, mass=1.0, temperature=0.0,
                     mass_profile=None, weight=None):
        """NumPy construction.  mass_profile(coords)->(N,) overrides `mass`
        point-by-point (localized heavy regions = test masses)."""
        L = graph_laplacian(lat, weight=weight)
        m2 = np.full(len(lat), float(mass) ** 2)
        if mass_profile is not None:
            m2 = np.asarray(mass_profile(lat.coords), dtype=float) ** 2
        W = L + np.diag(m2)
        C, P = _covariances_from_W(W, temperature)
        return cls(lat.coords, C, P, mass, temperature, source="numpy")

    @classmethod
    def disk(cls, radius, spacing, mass=1.0, temperature=0.0, **kw):
        from .lattice import disk_lattice
        return cls.from_lattice(disk_lattice(radius, spacing), mass,
                                temperature, **kw)

    @classmethod
    def ball(cls, radius, spacing, mass=1.0, temperature=0.0, **kw):
        from .lattice import ball_lattice
        return cls.from_lattice(ball_lattice(radius, spacing), mass,
                                temperature, **kw)

    @classmethod
    def engine(cls, radius, spacing, mass=1.0, temperature=0.0):
        """From the compiled CEFE engine (3D causal-diamond slice)."""
        coords, C, P = backend.engine_vacuum(radius, spacing, mass, temperature)
        return cls(coords, C, P, mass, temperature, source="cefe_core")

    # -- observables -------------------------------------------------------
    def entropy(self, mask):
        """Entanglement entropy of the subregion selected by boolean mask."""
        nu = symplectic_eigenvalues(self.C[np.ix_(mask, mask)],
                                    self.P[np.ix_(mask, mask)])
        return float(entropy_from_nu(nu))

    def mutual_information(self, maskA, maskB):
        return mutual_information(self.C, self.P, maskA, maskB)


def symplectic_eigenvalues(CA, PA):
    """Symplectic eigenvalues of the covariance pair (Williamson, stable)."""
    ev, U = np.linalg.eigh(CA)
    ev = np.clip(ev, 1e-30, None)
    Chalf = (U * np.sqrt(ev)[None, :]) @ U.T
    B = Chalf @ PA @ Chalf
    return np.sqrt(np.clip(np.linalg.eigvalsh(B), 0.25, None))


def entropy_from_nu(nu):
    """S = sum_k (nu+1/2) ln(nu+1/2) - (nu-1/2) ln(nu-1/2), 0 log 0 = 0."""
    nu = np.asarray(nu, dtype=float)
    x1 = nu + 0.5
    x2 = np.maximum(0.0, nu - 0.5)
    S = np.sum(x1 * np.log(x1))
    pos = x2 > 0
    S -= np.sum(x2[pos] * np.log(x2[pos]))
    return float(S)


# -- region helpers ---------------------------------------------------------

def ball_region(coords, radius, center=None):
    d = coords if center is None else coords - np.asarray(center)
    return np.linalg.norm(d, axis=1) <= radius


def annulus_region(coords, r_in, r_out, center=None):
    d = coords if center is None else coords - np.asarray(center)
    r = np.linalg.norm(d, axis=1)
    return (r > r_in) & (r <= r_out)


def mutual_information(C, P, maskA, maskB):
    """I(A:B) = S_A + S_B - S_AB  (UV-finite for disjoint closures)."""
    SA = entropy_from_nu(symplectic_eigenvalues(C[np.ix_(maskA, maskA)],
                                                P[np.ix_(maskA, maskA)]))
    SB = entropy_from_nu(symplectic_eigenvalues(C[np.ix_(maskB, maskB)],
                                                P[np.ix_(maskB, maskB)]))
    AB = maskA | maskB
    SAB = entropy_from_nu(symplectic_eigenvalues(C[np.ix_(AB, AB)],
                                                 P[np.ix_(AB, AB)]))
    return float(SA + SB - SAB)
