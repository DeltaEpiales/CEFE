#!/usr/bin/env python
"""
Entropic Gravity — 3D Schwarzschild visualization with the quantum substrate.

A known gravitational source (Schwarzschild black hole, Rs = 0.5) rendered in
3D, driven entirely by engine-measured quantum mechanics:

  * the holographic screen (wireframe sphere) carries N = 4*S_EE(R) bits,
  * screen temperature T = 2M/N (equipartition) colors the bits,
  * a test particle rides the emergent force F = T*dS/dx,
  * UNDERNEATH: the quantum mechanical functions measured on the same lattice:
      - S_EE(R)  Schwarzschild vs Minkowski (area law on curved background)
      - C(r) = <phi(x0) phi(x)> two-point vacuum correlator (the substrate)
      - F(R) emergent vs Newtonian 1/R^2 reference

The Schwarzschild lattice uses the engine's own proper-distance weighting
(tortoise-coordinate causal structure ported from Metric.hpp).

Output: entropic_gravity_3d.gif
"""

from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
from daimon_runtime import setup_plot
setup_plot()

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

# ---------------------------------------------------------------- metrics
RS = 0.5

def proper_distance_schw(p1, p2, steps=10):
    """Port of SchwarzschildMetric::get_spatial_distance (midpoint Simpson)."""
    p1 = np.asarray(p1); p2 = np.asarray(p2)
    d = (p2 - p1) / steps
    dl2 = float(d @ d) * steps * steps / steps**2  # |p2-p1|^2 / steps^2
    dist = 0.0
    for i in range(steps):
        c = p1 + (i + 0.5) * d
        cr = np.linalg.norm(c)
        if cr <= RS:
            return float(np.linalg.norm(p2 - p1))
        r0 = np.linalg.norm(p1 + i * d)
        r1 = np.linalg.norm(p1 + (i + 1) * d)
        dr = r1 - r0
        dist += np.sqrt(dl2 + (RS / (cr - RS)) * dr * dr)
    return dist

def build_t0_slice(R, ds, Rs=None):
    steps = int(np.ceil(R / ds))
    pts = []
    for ix in range(-steps, steps + 1):
        for iy in range(-steps, steps + 1):
            for iz in range(-steps, steps + 1):
                x, y, z = ix * ds, iy * ds, iz * ds
                r = np.sqrt(x * x + y * y + z * z)
                if r <= R and (Rs is None or r > Rs):
                    pts.append((ix, iy, iz, x, y, z))
    return pts

def build_laplacian(pts, ds, dist_fn=None):
    n = len(pts)
    index = {(p[0], p[1], p[2]): i for i, p in enumerate(pts)}
    L = np.zeros((n, n))
    for i, p in enumerate(pts):
        for d3 in range(3):
            nb = [p[0], p[1], p[2]]; nb[d3] += 1
            j = index.get(tuple(nb))
            if j is None:
                continue
            if dist_fn is None:
                w = 1.0 / (ds * ds)
            else:
                pd = dist_fn(p[3:6], pts[j][3:6])
                w = 1.0 / (pd * pd)
            L[i, j] += w; L[j, i] += w
            L[i, i] -= w; L[j, j] -= w
    return L

def covariance(L, mass):
    W = -L + mass * mass * np.eye(L.shape[0])
    evals, evecs = np.linalg.eigh(W)
    w = np.sqrt(np.maximum(evals, 1e-30))
    C = (evecs * (0.5 / w)) @ evecs.T
    P = (evecs * (0.5 * w)) @ evecs.T
    return C, P

def entropy_idx(C, P, idx):
    CA = C[np.ix_(idx, idx)]; PA = P[np.ix_(idx, idx)]
    ev, U = np.linalg.eigh(CA)
    Ch = U @ np.diag(np.sqrt(np.maximum(ev, 1e-30))) @ U.T
    B = Ch @ PA @ Ch
    nu = np.sqrt(np.maximum(np.linalg.eigvalsh(B), 0.25))
    x1 = nu + 0.5
    x2 = np.maximum(0.0, nu - 0.5)
    return np.sum(x1 * np.log(x1)) - np.sum(x2[x2 > 0] * np.log(x2[x2 > 0]))

# ---------------------------------------------------------------- compute
a, R_d, m_field = 0.7, 4.5, 1.0
radii = np.linspace(1.0, 3.5, 12)

print("[1/2] Schwarzschild background (proper-distance weighted lattice)...")
pts_s = build_t0_slice(R_d, a, Rs=RS)
coords_s = np.array([[p[3], p[4], p[5]] for p in pts_s])
r_s = np.linalg.norm(coords_s, axis=1)
C_s, P_s = covariance(build_laplacian(pts_s, a, proper_distance_schw), m_field)
S_schw = np.array([entropy_idx(C_s, P_s, np.where(r_s <= R)[0]) for R in radii])
print(f"      {len(pts_s)} dof; S(1.0)={S_schw[0]:.3f}  S(3.5)={S_schw[-1]:.3f}")

print("[2/2] Minkowski reference...")
pts_m = build_t0_slice(R_d, a)
coords_m = np.array([[p[3], p[4], p[5]] for p in pts_m])
r_m = np.linalg.norm(coords_m, axis=1)
C_m, P_m = covariance(build_laplacian(pts_m, a), m_field)
S_mink = np.array([entropy_idx(C_m, P_m, np.where(r_m <= R)[0]) for R in radii])
print(f"      {len(pts_m)} dof; S(1.0)={S_mink[0]:.3f}  S(3.5)={S_mink[-1]:.3f}")

# two-point correlator from a reference point on +x axis at r ~ 1.0
def correlator(C, coords, r_vals, r_ref=1.0):
    i_ref = np.argmin(np.abs(np.linalg.norm(coords, axis=1) - r_ref) + np.abs(coords[:, 1]) + np.abs(coords[:, 2]))
    rr = np.linalg.norm(coords, axis=1)
    out = []
    for rv in r_vals:
        sel = np.where(np.abs(rr - rv) < 0.6 * a)[0]
        out.append(np.mean(np.abs(C[i_ref, sel])) if len(sel) else np.nan)
    return np.array(out)

r_corr = np.linspace(1.0, 3.5, 26)
corr_s = correlator(C_s, coords_s, r_corr)
corr_m = correlator(C_m, coords_m, r_corr)

# entropic chain from MEASURED Schwarzschild entropy
M_central, m_test = 10.0, 1.0
S_use = np.maximum.accumulate(S_schw)
N_bits = 4.0 * S_use
T_scr = 2.0 * M_central / N_bits
F_em = T_scr * (2.0 * np.pi * m_test)
F_newton = M_central * m_test / radii**2

# ---------------------------------------------------------------- figure
plt.style.use("dark_background")
fig = plt.figure(figsize=(16, 9), facecolor="#0a0a0f")
fig.suptitle("Entropic Gravity of a Schwarzschild Black Hole — quantum substrate underneath, emergent force on top",
             color="white", fontsize=15, fontweight="bold")
CY, RD, GN, OR, GY = "#00e5ff", "#ff0055", "#00ff88", "#ff9900", "#888899"

ax3d = fig.add_axes([0.01, 0.03, 0.50, 0.88], projection="3d", facecolor="#0a0a0f")
gs = fig.add_gridspec(3, 1, left=0.56, right=0.97, top=0.90, bottom=0.08, hspace=0.55)
axS = fig.add_subplot(gs[0]); axC = fig.add_subplot(gs[1]); axF = fig.add_subplot(gs[2])
for ax in (axS, axC, axF):
    ax.set_facecolor("#111118"); ax.tick_params(colors="#666666", labelsize=8)
    ax.grid(True, ls=":", alpha=0.3, color="#666666")

# --- 3D static content
eq = coords_s[np.abs(coords_s[:, 2]) < a / 2]
sub = eq[np.random.default_rng(1).choice(len(eq), min(700, len(eq)), replace=False)]
ax3d.scatter(sub[:, 0], sub[:, 1], sub[:, 2], s=3, color="#3a3a4d", alpha=0.6)

u, v = np.mgrid[0:np.pi:12j, 0:2 * np.pi:24j]
ax3d.plot_surface(RS * np.sin(u) * np.cos(v), RS * np.sin(u) * np.sin(v), RS * np.cos(u),
                  color=(0.02, 0.02, 0.04), shade=False, zorder=10)
th = np.linspace(0, 2 * np.pi, 100)  # photon sphere ring at 1.5 Rs
ax3d.plot(1.5 * RS * np.cos(th), 1.5 * RS * np.sin(th), 0 * th, color=OR, lw=1.2, alpha=0.8)

uu, vv = np.mgrid[0:np.pi:9j, 0:2 * np.pi:13j]
def sphere_pts(R):
    return (R * np.sin(uu) * np.cos(vv), R * np.sin(uu) * np.sin(vv), R * np.cos(uu))
wire = ax3d.plot_wireframe(*sphere_pts(radii[-1]), color=CY, alpha=0.45, lw=0.7)
bits3d = ax3d.scatter([], [], [], s=24, c=[], cmap="inferno", vmin=T_scr.min(), vmax=T_scr.max())
particle, = ax3d.plot([], [], [], "o", color=RD, ms=9, mec="white", mew=1)
txt = ax3d.text2D(0.03, 0.95, "", transform=ax3d.transAxes, color="white",
                  fontsize=10, family="monospace", va="top",
                  bbox=dict(facecolor="#1a1a1d", edgecolor="#333333"))
ax3d.set(xlim=(-3.8, 3.8), ylim=(-3.8, 3.8), zlim=(-3.8, 3.8))
ax3d.set_box_aspect((1, 1, 1))
ax3d.set_axis_off()
ax3d.set_title("holographic screen + bit field (schematic 3D; panels: engine-measured)",
               color=CY, fontsize=11, pad=-2)

# --- 2D panels
axS.plot(radii, S_mink, "--", color=GY, lw=1.4, label="Minkowski (flat)")
axS.plot(radii, S_schw, "-", color=CY, lw=2.2, label="Schwarzschild (curved)")
ptS, = axS.plot([], [], "o", color="white", ms=6)
axS.set_title(r"quantum substrate I: entanglement entropy $S_{EE}(R)$", color=OR, fontsize=10)
axS.legend(fontsize=8, loc="upper left", facecolor="#1a1a1d", edgecolor="#333", labelcolor="#cccccc")

axC.plot(r_corr, corr_m, "--", color=GY, lw=1.4, label="Minkowski")
axC.plot(r_corr, corr_s, "-", color=CY, lw=2.2, label="Schwarzschild")
axC.set_yscale("log")
axC.set_title(r"quantum substrate II: vacuum correlator $\langle\phi(x_0)\phi(x)\rangle$", color=OR, fontsize=10)
axC.legend(fontsize=8, loc="upper right", facecolor="#1a1a1d", edgecolor="#333", labelcolor="#cccccc")

axF.plot(radii, F_newton, "-", color=GN, lw=1.6, alpha=0.5, label=r"Newton $1/R^2$")
axF.plot(radii, F_em, "-", color=RD, lw=1.2, alpha=0.3)
lineF, = axF.plot([], [], "-", color=RD, lw=2.6, label=r"emergent $F=T\nabla S$")
ptF, = axF.plot([], [], "o", color="white", ms=6)
axF.set_ylim(0, F_em.max() * 1.05)
axF.set_title("emergent force from measured bits", color=OR, fontsize=10)
axF.legend(fontsize=8, loc="upper right", facecolor="#1a1a1d", edgecolor="#333", labelcolor="#cccccc")
axF.set_xlabel("R", color="#aaaaaa", fontsize=9)

# ---------------------------------------------------------------- animation
n_frames = 60
sweep = np.concatenate([np.linspace(11, 0, 30).astype(int), np.linspace(0, 11, 30).astype(int)])

def fib_sphere(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5**0.5) * i
    return np.c_[np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)]

def update(f):
    global wire
    i = sweep[f]
    R = radii[i]
    wire.remove()
    wire = ax3d.plot_wireframe(*sphere_pts(R), color=CY, alpha=0.45, lw=0.7)
    nb = int(10 + 60 * N_bits[i] / N_bits.max())
    bpos = R * fib_sphere(nb)
    bits3d._offsets3d = (bpos[:, 0], bpos[:, 1], bpos[:, 2])
    bits3d.set_array(np.full(nb, T_scr[i]))
    ang = f * 0.35
    particle.set_data([R * np.cos(ang)], [R * np.sin(ang)])
    particle.set_3d_properties([0.0])
    txt.set_text(f"R = {R:.2f}\nN = {N_bits[i]:6.2f} bits\nT = {T_scr[i]:5.2f}\nF = {F_em[i]:6.2f}")
    ax3d.view_init(elev=24, azim=-55 + f * 0.9)
    ptS.set_data([R], [S_schw[i]])
    lineF.set_data(radii[:i + 1], F_em[:i + 1])
    ptF.set_data([R], [F_em[i]])
    return wire, bits3d, particle, txt, ptS, lineF, ptF

ani = animation.FuncAnimation(fig, update, frames=n_frames, blit=False, interval=120)
out = "entropic_gravity_3d.gif"
print(f"Rendering {n_frames} frames -> {out}")
ani.save(out, writer="pillow", fps=10, dpi=85)
plt.close(fig)
print("Saved", out)
