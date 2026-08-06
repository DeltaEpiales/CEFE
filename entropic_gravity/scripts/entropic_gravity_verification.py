#!/usr/bin/env python
"""
Independent verification harness for the CEFE engine + entropic-gravity chain.

The compiled cefe_core.cp311 .pyd cannot be loaded by the Pythons available on
this machine (3.12 / 3.14), and no C++ toolchain is present to rebuild it.
This script therefore ports the EXACT algorithms from the C++ sources
(src/geometry/CausalDiamond.cpp, src/core/EntropicFieldEngine.cpp) to NumPy and
re-runs the physics checks, then pushes the entropic-gravity test one level
deeper than tests/test_entropic_gravity.py does:

  - The existing test ASSUMES T = M/R^2 (so F ~ 1/R^2 by construction).
  - Here the screen bit-count N(R) is instead derived from the ENGINE-MEASURED
    entanglement entropy S_EE(R) via the holographic identification N = 4 S/k_B,
    so the inverse-square force must EMERGE from the measured area law.

Physics chain (Verlinde, JHEP 04 (2011) 029, arXiv:1001.0785):
    N  = A c^3 / (G hbar)  = 4 S_screen / k_B
    E  = (1/2) N k_B T = M c^2          ->  T = 2 M c^2 / (N k_B)
    dS = 2 pi k_B m c / hbar dx         (constant in R)
    F  = T dS/dx = 4 pi M m c^3 / (N hbar)  ~  1 / S_EE  ~  1 / R^2
"""

import numpy as np

# ----------------------------------------------------------------------------
# Exact NumPy port of CausalDiamondGrid (Minkowski, t = 0 slice)
# ----------------------------------------------------------------------------

def build_t0_slice(R, ds):
    """Integer-lattice points inside |r| <= R at t=0, with boundary flag.
    Mirrors CausalDiamondGrid::generate_grid + MinkowskiMetric."""
    steps = int(np.ceil(R / ds))
    pts = []
    for ix in range(-steps, steps + 1):
        for iy in range(-steps, steps + 1):
            for iz in range(-steps, steps + 1):
                x, y, z = ix * ds, iy * ds, iz * ds
                r = np.sqrt(x * x + y * y + z * z)
                if r <= R:
                    is_b = abs(r - R) <= ds / 2.0
                    pts.append((ix, iy, iz, x, y, z, is_b))
    return pts


def build_laplacian(pts, ds):
    """Graph Laplacian with 1/proper_ds^2 weights (Minkowski: proper_ds = ds).
    Mirrors CausalDiamondGrid::build_spatial_laplacian (positive-offset
    neighbour search, negative diagonal)."""
    n = len(pts)
    index = {(p[0], p[1], p[2]): i for i, p in enumerate(pts)}
    L = np.zeros((n, n))
    inv_ds2 = 1.0 / (ds * ds)
    for i, p in enumerate(pts):
        for d in range(3):
            nb = [p[0], p[1], p[2]]
            nb[d] += 1
            j = index.get(tuple(nb))
            if j is not None:
                L[i, j] += inv_ds2
                L[j, i] += inv_ds2
                L[i, i] -= inv_ds2
                L[j, j] -= inv_ds2
    return L


# ----------------------------------------------------------------------------
# Exact NumPy port of EntropicFieldEngine::evolve_to_slice
#                       + compute_entanglement_entropy(_indices)
# ----------------------------------------------------------------------------

def covariance_matrices(L, mass, temperature=0.0):
    """Ground-state (or thermal) covariance matrices of H = 1/2 pi^2 + 1/2 phi W phi.
    C = (1/2) W^{-1/2} coth(beta W^{1/2}/2),  P = (1/2) W^{+1/2} coth(...)."""
    W = -L + mass * mass * np.eye(L.shape[0])
    evals, evecs = np.linalg.eigh(W)
    w = np.sqrt(np.maximum(evals, 1e-30))
    factor = np.ones_like(w)
    if temperature > 0.0:
        beta = 1.0 / temperature
        factor = 1.0 / np.tanh(beta * w / 2.0)
    c = 0.5 * (1.0 / w) * factor
    p = 0.5 * w * factor
    C = (evecs * c) @ evecs.T
    P = (evecs * p) @ evecs.T
    return C, P


def entanglement_entropy(C, P, coords, subregion_r):
    """Srednicki (1993) / Bombelli et al. (1986) covariance-matrix entropy."""
    r = np.linalg.norm(coords, axis=1)
    A = np.where(r <= subregion_r)[0]
    if len(A) == 0:
        return 0.0
    CA = C[np.ix_(A, A)]
    PA = P[np.ix_(A, A)]
    nu2 = np.linalg.eigvals(CA @ PA).real
    nu = np.sqrt(np.maximum(0.25, nu2))     # enforce nu >= 1/2 (uncertainty)
    x1 = nu + 0.5
    x2 = np.maximum(0.0, nu - 0.5)
    S = np.sum(x1 * np.log(x1))
    S -= np.sum(x2[x2 > 0] * np.log(x2[x2 > 0]))
    return S


def power_law_fit(x, y):
    lx, ly = np.log(x), np.log(y)
    b, log_a = np.polyfit(lx, ly, 1)
    pred = log_a + b * lx
    r2 = 1 - np.sum((ly - pred) ** 2) / np.sum((ly - ly.mean()) ** 2)
    return np.exp(log_a), b, r2


# ----------------------------------------------------------------------------
# Verification suite
# ----------------------------------------------------------------------------

def main():
    print("=" * 72)
    print("CEFE ENGINE REPLICATION + ENTROPIC GRAVITY VERIFICATION (NumPy port)")
    print("=" * 72)

    # ---- Test 1: Area law (Srednicki 1993: S ~ R^{d-1} = R^2 in 3+1 D) ----
    print("\n[1] AREA LAW  S(R) ~ R^alpha   (expect alpha ~ 2)")
    R_d, ds, m = 5.0, 0.8, 1.0
    pts = build_t0_slice(R_d, ds)
    coords = np.array([[p[3], p[4], p[5]] for p in pts])
    L = build_laplacian(pts, ds)
    C, P = covariance_matrices(L, m)
    radii = np.array([1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0])
    S = np.array([entanglement_entropy(C, P, coords, r) for r in radii])
    a, alpha, r2 = power_law_fit(radii, S)
    for rr, ss in zip(radii, S):
        print(f"    R={rr:4.1f}   S={ss:9.4f}")
    print(f"    -> alpha = {alpha:.3f},  R^2 fit = {r2:.4f}   (README reports 2.12 / 0.995)")

    # ---- Test 2: UV scaling S ~ 1/eps^2 (Bombelli/Srednicki divergence) ----
    print("\n[2] UV SCALING  S vs 1/eps^2 at fixed geometry (expect ~linear)")
    inv_eps2, S_uv = [], []
    for eps in [1.2, 1.0, 0.8, 0.7, 0.6]:
        pts = build_t0_slice(4.0, eps)
        coords = np.array([[p[3], p[4], p[5]] for p in pts])
        C, P = covariance_matrices(build_laplacian(pts, eps), 1.0)
        S_uv.append(entanglement_entropy(C, P, coords, 1.5))
        inv_eps2.append(1.0 / eps ** 2)
        print(f"    eps={eps:4.1f}   1/eps^2={inv_eps2[-1]:6.3f}   S={S_uv[-1]:8.4f}")
    slope, intercept = np.polyfit(inv_eps2, S_uv, 1)
    pred = slope * np.array(inv_eps2) + intercept
    r2uv = 1 - np.sum((np.array(S_uv) - pred) ** 2) / np.sum((np.array(S_uv) - np.mean(S_uv)) ** 2)
    print(f"    -> linear fit R^2 = {r2uv:.4f}   (README reports 0.947)")

    # ---- Test 3: Mass gap (entanglement suppressed by mass) ----
    print("\n[3] MASS GAP  S(m) monotonically decreasing")
    S_m = []
    for mm in [0.1, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0]:
        pts = build_t0_slice(4.0, 0.8)
        coords = np.array([[p[3], p[4], p[5]] for p in pts])
        C, P = covariance_matrices(build_laplacian(pts, 0.8), mm)
        S_m.append(entanglement_entropy(C, P, coords, 1.5))
        print(f"    m={mm:5.1f}   S={S_m[-1]:8.4f}")
    print(f"    -> monotonic decrease: {all(S_m[i] >= S_m[i+1] for i in range(len(S_m)-1))}")

    # ---- Test 4: Thermal states S(T) >= S(0) ----
    print("\n[4] THERMAL  S(T) >= S(0), monotonically increasing")
    pts = build_t0_slice(3.0, 0.8)
    coords = np.array([[p[3], p[4], p[5]] for p in pts])
    Lt = build_laplacian(pts, 0.8)
    S_T = []
    for T in [0.0, 0.5, 1.0, 2.0, 5.0, 10.0]:
        C, P = covariance_matrices(Lt, 1.0, temperature=T)
        S_T.append(entanglement_entropy(C, P, coords, 1.0))
        print(f"    T={T:5.1f}   S={S_T[-1]:8.4f}")
    print(f"    -> monotonic increase: {all(S_T[i] <= S_T[i+1] for i in range(len(S_T)-1))}")

    # ---- Test 5: ENTROPIC GRAVITY -- force derived from measured S_EE ----
    print("\n[5] ENTROPIC GRAVITY: F(R) from engine-measured entropy (not assumed)")
    print("    Verlinde chain: N=4S/kB, T=2Mc^2/(NkB), F=T*dS/dx => F ~ 1/S_EE")
    R_d, ds, m = 5.0, 0.8, 1.0
    pts = build_t0_slice(R_d, ds)
    coords = np.array([[p[3], p[4], p[5]] for p in pts])
    C, P = covariance_matrices(build_laplacian(pts, ds), m)
    Rg = np.linspace(1.5, 4.0, 11)
    Sg = np.array([entanglement_entropy(C, P, coords, r) for r in Rg])
    # N(R) = 4 S_EE(R)/k_B  (holographic identification, k_B = 1 in nat. units)
    N_bits = 4.0 * Sg
    M_central, m_test = 10.0, 1.0
    dSdx = 2.0 * np.pi * m_test          # Verlinde eq. (3.2), c = hbar = k_B = 1
    T_screen = 2.0 * M_central / N_bits  # equipartition E = 1/2 N T = M
    F = T_screen * dSdx                  # entropic force
    aF, betaF, r2F = power_law_fit(Rg, F)
    print(f"    {'R':>5} {'S_EE':>9} {'T(R)':>9} {'F(R)':>9}")
    for i in range(len(Rg)):
        print(f"    {Rg[i]:5.2f} {Sg[i]:9.4f} {T_screen[i]:9.4f} {F[i]:9.4f}")
    print(f"    -> F ~ R^({betaF:.3f})   (Newton: -2.000)   R^2 fit = {r2F:.4f}")
    # Exact ratio check at the endpoints
    print(f"    -> F(R1)/F(R2) = {F[0]/F[-1]:.3f}  vs  (R2/R1)^2 = {(Rg[-1]/Rg[0])**2:.3f}")

    # ---- Test 6: Entropic atom spectrum sanity (tests/theory/entropic_atom) --
    print("\n[6] ENTROPIC ATOM: hydrogen-like spectrum of H = -1/2 lap - alpha/r")
    print("    (resolvable regime alpha=0.8, box R=6, ds=0.75; Dirichlet at wall)")
    alpha_c, R_a, ds_a = 0.8, 6.0, 0.75
    pts = build_t0_slice(R_a, ds_a)
    coords = np.array([[p[3], p[4], p[5]] for p in pts])
    La = build_laplacian(pts, ds_a)
    rr = np.linalg.norm(coords, axis=1)
    V = -alpha_c / (rr + 0.1 * ds_a)
    H = -0.5 * La + np.diag(V)
    evals = np.linalg.eigvalsh(H)[:6]
    E1_cont = -0.5 * alpha_c ** 2   # continuum infinite-space ground state
    print(f"    N_dof = {len(pts)}")
    print(f"    lowest eigenvalues: {np.array2string(evals, precision=4)}")
    print(f"    continuum E_1 = {E1_cont:.4f};  lattice E_1 = {evals[0]:.4f}")
    if evals[0] < 0 and evals[1] < 0:
        ratio = evals[0] / evals[1]
        print(f"    E_1/E_2 = {ratio:.3f}  (hydrogen: 4.000 in infinite space;")
        print(f"               Dirichlet box pushes excited states up -> ratio < 4 expected)")

    print("\n" + "=" * 72)
    print("VERIFICATION COMPLETE")
    print("=" * 72)


if __name__ == "__main__":
    main()
