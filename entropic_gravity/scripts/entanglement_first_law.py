#!/usr/bin/env python
"""
ENTANGLEMENT FIRST LAW + ARROW OF TIME  (CEFE engine / NumPy port)
==================================================================

Two tests at the heart of the unified entropic-gravity theory:

PART 1 -- FIRST LAW OF ENTANGLEMENT:  dS_A = d<K_A>
    For any small perturbation of a state, the change in subregion entropy
    equals the change in expectation of the modular Hamiltonian K_A
    (relative-entropy positivity; the identity underlying Jacobson's
    "Entanglement Equilibrium and the Einstein Equation", PRL 116, 201101
    (2016)). K_A for a Gaussian state is built via Williamson decomposition
    of the vacuum covariance pair (C_A, P_A):
        symplectic eigenvalues nu_k  ->  modular energies
        eps_k = ln[(nu_k + 1/2)/(nu_k - 1/2)].
    Perturbations: (i) mass quench m -> m + dm, (ii) thermal excitation T > 0.

PART 2 -- ARROW OF TIME: fine- vs coarse-grained entropy under unitary
    evolution. The entropic-time axiom (Osacra-Gravaser papers 2-3;
    Barontini et al. arXiv:2509.07745) requires unitary dynamics to keep
    fine-grained entropy constant while coarse-grained entropy grows.
    Schrodinger evolution (Crank-Nicolson, as in QuantumWaveEngine) of a
    Gaussian packet on the t=0 causal-diamond slice:
        purity Tr(rho^2) must stay exactly 1;
        coarse-grained position Shannon entropy must grow.
"""

import numpy as np
from entropic_gravity_verification import (
    build_t0_slice, build_laplacian, covariance_matrices, entanglement_entropy,
)

print("=" * 74)
print("PART 1: FIRST LAW OF ENTANGLEMENT  dS_A = d<K_A>")
print("=" * 74)

# ---- vacuum state on t=0 slice --------------------------------------------
a, R_d, m0, r_A = 0.8, 4.0, 1.0, 1.6
pts = build_t0_slice(R_d, a)
coords = np.array([[p[3], p[4], p[5]] for p in pts])
L = build_laplacian(pts, a)
C, P = covariance_matrices(L, m0)

r = np.linalg.norm(coords, axis=1)
A = np.where(r <= r_A)[0]
CA = C[np.ix_(A, A)]
PA = P[np.ix_(A, A)]
nA = len(A)
print(f"subregion A: {nA} dof (r <= {r_A}), full slice: {len(pts)} dof")

# ---- Williamson decomposition of (CA, PA) ---------------------------------
evalsC, UC = np.linalg.eigh(CA)
Chalf = UC @ np.diag(np.sqrt(np.maximum(evalsC, 1e-30))) @ UC.T
Chalfinv = UC @ np.diag(1.0 / np.sqrt(np.maximum(evalsC, 1e-30))) @ UC.T

B = Chalf @ PA @ Chalf
evalsB, O = np.linalg.eigh(B)
nu = np.sqrt(np.maximum(evalsB, 0.250001))          # symplectic eigenvalues (clamped)
nu_k = nu

# normal-coordinate maps:  q = S phi ,  p~ = S^{-T} pi
nu12 = np.sqrt(nu)
S_map = np.diag(nu12) @ O.T @ Chalfinv               # S
SinvT = np.diag(1.0 / nu12) @ O.T @ Chalf            # S^{-T}
# modular Hamiltonian of a Gaussian mode: K_k = eps_k (n_k + 1/2) = (eps_k/2)(q_k^2 + p~_k^2)
eps = np.log((nu + 0.5) / (nu - 0.5))                # modular energies
eps_k = eps

def entropy_from_nu(nu_vals):
    x1 = nu_vals + 0.5
    x2 = np.maximum(0.0, nu_vals - 0.5)
    out = np.sum(x1 * np.log(x1))
    out -= np.sum(x2[x2 > 0] * np.log(x2[x2 > 0]))
    return out

def symplectic_spectrum(Cx, Px):
    """Symplectic eigenvalues of a covariance pair (for the perturbed entropy)."""
    ev, Ux = np.linalg.eigh(Cx)
    Ch = Ux @ np.diag(np.sqrt(np.maximum(ev, 1e-30))) @ Ux.T
    Bx = Ch @ Px @ Ch
    return np.sqrt(np.maximum(np.linalg.eigvalsh(Bx), 0.25))

def K_expectation(Cx, Px):
    """<K_A> of a state with covariances (Cx, Px) on A, in the VACUUM modular basis.
    K_A = sum_k (eps_k/2)(q_k^2 + p~_k^2); vacuum value = sum_k eps_k nu_k."""
    Qq = S_map @ Cx @ S_map.T            # <q q^T>
    Qp = SinvT @ Px @ SinvT.T            # <p~ p~^T>
    q2 = np.diag(Qq)
    p2 = np.diag(Qp)
    return np.sum(eps_k / 2.0 * (q2 + p2))

S_vac = entropy_from_nu(nu)
K_vac = K_expectation(CA, PA)
print(f"S_A(vac) = {S_vac:.6f}   <K_A>(vac) = {K_vac:.6f}")

# ---- perturbation (i): mass quench -----------------------------------------
print("\n(i) mass quench  m -> m + dm")
print(f"{'dm':>7} {'dS_A':>11} {'d<K_A>':>11} {'|diff|/dS':>11}")
for dm in [0.02, 0.05, 0.10, 0.20]:
    C2, P2 = covariance_matrices(L, m0 + dm)
    CA2, PA2 = C2[np.ix_(A, A)], P2[np.ix_(A, A)]
    S_new = entropy_from_nu(symplectic_spectrum(CA2, PA2))
    K_new = K_expectation(CA2, PA2)
    dS, dK = S_new - S_vac, K_new - K_vac
    rel = abs(dS - dK) / abs(dS) if abs(dS) > 1e-14 else 0.0
    print(f"{dm:7.2f} {dS:11.6f} {dK:11.6f} {rel:11.2e}")

# ---- perturbation (ii): thermal --------------------------------------------
print("\n(ii) thermal excitation  T > 0")
print(f"{'T':>7} {'dS_A':>11} {'d<K_A>':>11} {'|diff|/dS':>11}")
for T in [0.02, 0.05, 0.10, 0.20]:
    C2, P2 = covariance_matrices(L, m0, temperature=T)
    CA2, PA2 = C2[np.ix_(A, A)], P2[np.ix_(A, A)]
    S_new = entropy_from_nu(symplectic_spectrum(CA2, PA2))
    K_new = K_expectation(CA2, PA2)
    dS, dK = S_new - S_vac, K_new - K_vac
    rel = abs(dS - dK) / abs(dS) if abs(dS) > 1e-14 else 0.0
    print(f"{T:7.2f} {dS:11.6f} {dK:11.6f} {rel:11.2e}")

print("\n" + "=" * 74)
print("PART 2: ARROW OF TIME -- fine-grained constant, coarse-grained grows")
print("=" * 74)

# 2D cross-section (z=0) of the causal diamond t=0 slice
R2d, a2 = 2.5, 0.25
pts2 = build_t0_slice(R2d, a2)
sel = [i for i, p in enumerate(pts2) if abs(p[5]) < 1e-12]
coords2 = np.array([[pts2[i][3], pts2[i][4]] for i in sel])
ix2 = np.round(coords2[:, 0] / a2).astype(int)
iy2 = np.round(coords2[:, 1] / a2).astype(int)
index2 = {(ix2[i], iy2[i]): i for i in range(len(coords2))}
n2 = len(coords2)
L2 = np.zeros((n2, n2))
inv_a2 = 1.0 / a2 ** 2
for i in range(n2):
    for dx, dy in [(1, 0), (0, 1)]:
        j = index2.get((ix2[i] + dx, iy2[i] + dy))
        if j is not None:
            L2[i, j] += inv_a2; L2[j, i] += inv_a2
            L2[i, i] -= inv_a2; L2[j, j] -= inv_a2
m_qm = 1.0
H = -0.5 / m_qm * L2
print(f"2D slice dof = {n2}")

# initial Gaussian packet, momentum px
x0, sig, px = -1.2, 0.4, 3.0
dx = coords2[:, 0] - x0
psi = np.exp(-(dx**2 + coords2[:, 1]**2) / (4 * sig**2)) * np.exp(1j * px * dx)
psi /= np.linalg.norm(psi)

dt = 0.02
steps = 220
Idt = np.eye(n2, dtype=complex)
CN_lhs = Idt + 1j * dt / 2 * H
CN_rhs = Idt - 1j * dt / 2 * H
lhs_inv = np.linalg.inv(CN_lhs)

bin_w = 4 * a2
bx = np.floor(coords2[:, 0] / bin_w).astype(int)
by = np.floor(coords2[:, 1] / bin_w).astype(int)
bins = {}
for i in range(n2):
    bins.setdefault((bx[i], by[i]), []).append(i)

def coarse_entropy(psi_v):
    p = np.abs(psi_v) ** 2
    S = 0.0
    for cell in bins.values():
        pc = p[cell].sum()
        if pc > 1e-15:
            S -= pc * np.log(pc)
    return S

print(f"{'t':>6} {'purity':>9} {'S_coarse':>9}")
S_cg = [coarse_entropy(psi)]
purities = [float(np.vdot(psi, psi).real ** 2)]
times = [0.0]
for n in range(1, steps + 1):
    psi = lhs_inv @ (CN_rhs @ psi)
    if n % 20 == 0:
        t = n * dt
        pur = float(np.vdot(psi, psi).real ** 2)
        sc = coarse_entropy(psi)
        purities.append(pur); S_cg.append(sc); times.append(t)
        print(f"{t:6.2f} {pur:9.6f} {sc:9.4f}")

print(f"\npurity drift over run: {abs(purities[-1] - 1.0):.2e}  (unitary: exact 0)")
print(f"coarse-grained entropy: {S_cg[0]:.4f} -> {max(S_cg):.4f}  (grew: {max(S_cg) > S_cg[0]})")
print("=" * 74)
print("DONE")
print("=" * 74)
