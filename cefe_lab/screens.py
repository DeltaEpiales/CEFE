"""
cefe_lab.screens -- holographic-screen utilities: entropy profiles,
the Verlinde loop, and force-law systematics.
"""

import numpy as np


def entropy_profile(state, radii, region_fn=None):
    """S(R) for a sequence of ball radii (or a custom region function)."""
    from .gaussian import ball_region
    fn = region_fn or ball_region
    return np.array([state.entropy(fn(state.coords, r)) for r in radii])


def power_law_fit(x, y):
    """log-log fit y ~ x^beta: returns (prefactor, beta, R^2)."""
    lx, ly = np.log(np.asarray(x, float)), np.log(np.asarray(y, float))
    beta, loga = np.polyfit(lx, ly, 1)
    pred = loga + beta * lx
    ss_res = np.sum((ly - pred) ** 2)
    ss_tot = np.sum((ly - ly.mean()) ** 2)
    return float(np.exp(loga)), float(beta), float(1 - ss_res / ss_tot)


def verlinde_force(S, M_enclosed=10.0, m_test=1.0):
    """The Verlinde loop on a measured entropy profile:
       N = 4S bits,  T = 2M/N (equipartition),  F = T * dS/dx,
       with the universal gradient dS/dx = 2 pi m (Axiom 2)."""
    S = np.asarray(S, float)
    N = 4.0 * S
    T = 2.0 * M_enclosed / N
    return T * (2.0 * np.pi * m_test)


def g_eff_per_species(S, a, radius):
    """G_eff / a^2 = 1/(4 c2)  with  S = c2 A / a^2,  A = 4 pi R^2."""
    c2 = S * a * a / (4.0 * np.pi * radius ** 2)
    return 1.0 / (4.0 * c2), c2
