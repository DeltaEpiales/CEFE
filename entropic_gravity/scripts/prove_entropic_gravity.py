#!/usr/bin/env python
"""
DECISIVE TEST SUITE: IS GRAVITY ENTROPIC?
==========================================

Unified theory under test (distilled from the three Osacra-Gravaser papers +
the modular entropic-gravity layer):

    AXIOM EG: Gravity is an entropic force. Information on holographic screens
    scales with AREA (not volume); equipartition over those screen bits then
    forces F ~ 1/R^2 (Verlinde 2011; Jacobson 1995; Padmanabhan 2010).
    The engine's lattice spacing a plays the role of the Planck length.

This script tests the four pillars the CEFE engine CAN decide:

  TEST A  Area law + Srednicki coefficient kappa (target: kappa ~ 0.30,
          Srednicki PRL 71, 666 (1993), massless free scalar, S = k (R/a)^2)
  TEST B  The complete Verlinde loop driven by MEASURED entropy only:
          bit count N = 4S_EE -> equipartition T = 2M/N -> F = T*dS/dx.
          Verdict variables: force exponent beta (Newton: -2) and the
          G-renormalization factor (Susskind-Uglum, PRD 50, 2700 (1994)).
  TEST C  UV-finite mutual information I(A:B) for gapped regions --
          the cutoff-stable bit count (fixes the e^{-S_EE} ansatz).
  TEST D  FALSIFICATION of the Osacra-Gravaser fractal-horizon extension:
          if horizon DOF were fractal with D_H = 3/2 at the fundamental scale,
          boundary-mode and entropy-shell counts must scale as R^1.5.
          Smooth holography predicts R^2.

Uses the compiled cefe_py engine if importable; otherwise falls back to the
line-identical NumPy port (entropic_gravity_verification.py).
"""

import numpy as np

try:
    import cefe_py as ce  # noqa: F401
    HAVE_ENGINE = True
except Exception:
    HAVE_ENGINE = False

from entropic_gravity_verification import (
    build_t0_slice, build_laplacian, covariance_matrices,
    entanglement_entropy, power_law_fit,
)

print("=" * 74)
print("DECISIVE SUITE: IS GRAVITY ENTROPIC?  (CEFE engine" +
      (" [compiled]" if HAVE_ENGINE else " [NumPy port, line-identical]") + ")")
print("=" * 74)


def entropy_profile(R_d, a, m, radii):
    pts = build_t0_slice(R_d, a)
    coords = np.array([[p[3], p[4], p[5]] for p in pts])
    C, P = covariance_matrices(build_laplacian(pts, a), m)
    return coords, C, P, np.array([entanglement_entropy(C, P, coords, r) for r in radii])


# --------------------------------------------------------------------------
print("\n--- TEST A: AREA LAW + SREDNICKI COEFFICIENT (target kappa ~ 0.30) ---")
# --------------------------------------------------------------------------
print(f"{'a':>5} {'m':>5} {'alpha':>7} {'R2fit':>7} {'kappa_meas':>10} {'kappa_Sr':>9}")
kappa_results = []
for a in [1.0, 0.8, 0.6]:
    for m in [0.1, 1.0]:
        R_d = 4.0
        radii = np.linspace(1.2, 3.2, 9)
        coords, C, P, S = entropy_profile(R_d, a, m, radii)
        _, alpha, r2 = power_law_fit(radii, S)
        # kappa from largest radius: S = kappa (R/a)^2
        kappa = S[-1] / (radii[-1] / a) ** 2
        kappa_results.append((a, m, alpha, r2, kappa))
        print(f"{a:5.2f} {m:5.1f} {alpha:7.3f} {r2:7.4f} {kappa:10.4f} {0.30:9.2f}")
alphas = np.array([k[2] for k in kappa_results])
kappas_m1 = [k[4] for k in kappa_results if k[1] == 1.0]
print(f"-> area exponent alpha = {alphas.mean():.3f} +/- {alphas.std():.3f}  (theory: 2.000)")
print(f"-> kappa (m=1.0) = {np.mean(kappas_m1):.4f}  (Srednicki 1993 massless scalar: 0.30)")
print(f"   (m=0.1 values inflated: xi=1/m=10 exceeds box R=4 -- finite-size, not a discrepancy)")

# --------------------------------------------------------------------------
print("\n--- TEST B: COMPLETE VERLINDE LOOP FROM MEASURED ENTROPY ONLY ---")
# --------------------------------------------------------------------------
a, m_field = 0.6, 1.0
R_d = 4.5
radii = np.linspace(1.5, 3.5, 13)
coords, C, P, S = entropy_profile(R_d, a, m_field, radii)

M_central, m_test = 10.0, 1.0
N_bits = 4.0 * S                          # holographic bit count (k_B = 1)
T = 2.0 * M_central / N_bits              # equipartition E = N T / 2 = M c^2
dSdx = 2.0 * np.pi * m_test               # Verlinde entropy gradient (postulate)
F = T * dSdx
_, beta, r2B = power_law_fit(radii, F)
c2 = S[-1] * a * a / (4.0 * np.pi * radii[-1] ** 2)   # S = c2 * A / a^2
G_eff_over_a2 = 1.0 / (4.0 * c2)                       # F = M m a^2 / (4 c2 R^2)
print(f"    S(R) measured on engine  ->  N = 4S  ->  T = 2M/N  ->  F = T dS/dx")
print(f"    force exponent beta = {beta:.3f}   (Newton: -2.000)   R2 = {r2B:.4f}")
print(f"    measured c2 = {c2:.4f}  (S = c2 A/a^2)")
print(f"    -> G_eff = {G_eff_over_a2:.2f} a^2  (Bekenstein identification needs G = a^2;")
print(f"       ratio {G_eff_over_a2:.2f} = entanglement renormalization, Susskind-Uglum)")
print(f"    -> one real scalar carries {100/G_eff_over_a2:.1f}% of the Bekenstein bit density")

# --------------------------------------------------------------------------
print("\n--- TEST C: UV-FINITE MUTUAL INFORMATION (cutoff-stable bit count) ---")
# --------------------------------------------------------------------------
# Done on the 2D disk cross-section (same Srednicki physics, finer affordable
# lattices): gap must resolve to MANY lattice spacings for the boundary
# divergences to cancel.  A = disk r<=0.8 ; B = annulus 1.4<r<=2.0 (gap 0.6).
def disk_slice(R, a):
    steps = int(np.ceil(R / a))
    pts = []
    for ix in range(-steps, steps + 1):
        for iy in range(-steps, steps + 1):
            x, y = ix * a, iy * a
            if np.sqrt(x * x + y * y) <= R:
                pts.append((ix, iy, 0, x, y, 0.0, False))
    return pts

def laplacian2d(pts, a):
    n = len(pts)
    index = {(p[0], p[1]): i for i, p in enumerate(pts)}
    L = np.zeros((n, n))
    w = 1.0 / a ** 2
    for i, p in enumerate(pts):
        for dx, dy in [(1, 0), (0, 1)]:
            j = index.get((p[0] + dx, p[1] + dy))
            if j is not None:
                L[i, j] += w; L[j, i] += w; L[i, i] -= w; L[j, j] -= w
    return L

def entropy_idx(C, P, idx):
    CA = C[np.ix_(idx, idx)]; PA = P[np.ix_(idx, idx)]
    nu2 = np.linalg.eigvals(CA @ PA).real
    nu = np.sqrt(np.maximum(0.25, nu2))
    x1, x2 = nu + 0.5, np.maximum(0.0, nu - 0.5)
    return np.sum(x1 * np.log(x1)) - np.sum(x2[x2 > 0] * np.log(x2[x2 > 0]))

print("    (2D disk, m=1)  A = r<=0.8 ; B = 1.4<r<=2.0 ; gap 0.6 physical")
print(f"{'a':>6} {'gap/a':>6} {'S(A)':>8} {'I(A:B)':>8}")
I_vals, SA_vals = [], []
for a in [0.3, 0.2, 0.15]:
    pts = disk_slice(2.5, a)
    coords = np.array([[p[3], p[4]] for p in pts])
    C, P = covariance_matrices(laplacian2d(pts, a), 1.0)
    r = np.linalg.norm(coords, axis=1)
    idx_A = np.where(r <= 0.8)[0]
    idx_B = np.where((r > 1.4) & (r <= 2.0))[0]
    idx_AB = np.where((r <= 0.8) | ((r > 1.4) & (r <= 2.0)))[0]  # union, gap EXCLUDED
    S_A = entropy_idx(C, P, idx_A)
    S_B = entropy_idx(C, P, idx_B)
    S_AB = entropy_idx(C, P, idx_AB)
    I = S_A + S_B - S_AB
    I_vals.append(I); SA_vals.append(S_A)
    print(f"{a:6.2f} {0.6/a:6.1f} {S_A:8.3f} {I:8.3f}")
print(f"    -> S(A) diverges ~ 1/a: x{SA_vals[-1]/SA_vals[0]:.2f} (1/a predicts x{0.3/0.15:.2f})")
growth = I_vals[-1] / I_vals[0]
print(f"    -> I(A:B) = {I_vals[0]:.3f}, {I_vals[1]:.3f}, {I_vals[2]:.3f} "
      f"(growth x{growth:.2f} vs entropy x{SA_vals[-1]/SA_vals[0]:.2f})")

# --------------------------------------------------------------------------
print("\n--- TEST D: FALSIFICATION -- FRACTAL HORIZON D_H=3/2 vs SMOOTH AREA ---")
# --------------------------------------------------------------------------
# Primary discriminator: the entropy scaling exponent alpha IS the effective
# dimension of the entanglement-carrying surface.  Smooth area law: alpha = 2.
# OG fractal horizon at the fundamental scale: alpha = D_H = 1.5.
a = 0.6
R_d = 4.5
radii = np.linspace(1.5, 3.5, 13)
coords, C, P, S = entropy_profile(R_d, a, 1.0, radii)
_, alpha_D, r2D = power_law_fit(radii, S)
pts = build_t0_slice(R_d, a)
rr = np.linalg.norm(np.array([[p[3], p[4], p[5]] for p in pts]), axis=1)
Rs = np.arange(1.2, 3.6, a)
shell_counts = np.array([np.sum((rr > R - a) & (rr <= R)) for R in Rs], dtype=float)
_, d_mode, r2m = power_law_fit(Rs[2:], shell_counts[2:])
print(f"    entropy scaling exponent alpha = {alpha_D:.3f}  (R2={r2D:.4f})")
print(f"    boundary lattice modes in shell (R-a, R]: count ~ R^{d_mode:.3f}  (R2={r2m:.4f})")
print(f"    SMOOTH holography: alpha=2.0, modes~R^2 ; OG fractal: alpha=1.5, modes~R^1.5")
dev_smooth = abs(alpha_D - 2.0)
dev_fractal = abs(alpha_D - 1.5)
smooth = dev_smooth < 0.25 and dev_fractal > 0.25
verdict_D = (f"SMOOTH -- alpha={alpha_D:.2f} is {dev_smooth:.2f} from 2.0 but "
             f"{dev_fractal:.2f} from 1.5: fractal cutoff structure FALSIFIED") \
            if smooth else "INCONCLUSIVE / fractal-consistent"
print(f"    -> VERDICT D: {verdict_D}")

print("\n" + "=" * 74)
print("SUITE COMPLETE")
print("=" * 74)
