#!/usr/bin/env python
"""
Entropic Gravity — engine-driven visualization animation.

Unlike animate_entropic_gravity.py (which ASSUMES T = M/R^2), this animation
drives every panel from the entropy S(R) MEASURED on the causal-diamond lattice
(CEFE engine algorithm):

    screen bit count  N(R) = 4 S_EE(R)
    screen temperature T(R) = 2M / N(R)          (equipartition)
    entropic force    F(R) = T(R) * dS/dx        (Verlinde, dS/dx = 2 pi m)

Output: entropic_gravity_verified.gif
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

from entropic_gravity_verification import (
    build_t0_slice, build_laplacian, covariance_matrices, entanglement_entropy,
)

# ---------------------------------------------------------------- engine data
print("Measuring S(R) on the causal-diamond lattice (engine)...")
R_d, a, m_field = 4.5, 0.6, 1.0
pts = build_t0_slice(R_d, a)
coords3 = np.array([[p[3], p[4], p[5]] for p in pts])
C, P = covariance_matrices(build_laplacian(pts, a), m_field)

radii_all = np.linspace(1.5, 3.5, 25)
S_all = np.array([entanglement_entropy(C, P, coords3, r) for r in radii_all])
# smooth monotone envelope (lattice plateaus removed for clean animation)
S_mono = np.maximum.accumulate(S_all)

M_central, m_test = 10.0, 1.0
dSdx = 2.0 * np.pi * m_test
N_all = 4.0 * S_mono
T_all = 2.0 * M_central / N_all
F_all = T_all * dSdx
F_newton = M_central * m_test * a**2 / (4 * 0.025 * radii_all**2)  # G_eff = a^2/(4 c2)

# 2D cross-section points for the grid panel
sel = np.abs(coords3[:, 2]) < a / 2
xy = coords3[sel][:, :2]

# ---------------------------------------------------------------- figure
plt.style.use("dark_background")
fig = plt.figure(figsize=(16, 9), facecolor="#0a0a0f")
gs = fig.add_gridspec(2, 2, hspace=0.32, wspace=0.25,
                      left=0.06, right=0.97, top=0.88, bottom=0.09)
fig.suptitle("Emergence of Entropic Gravity — driven by engine-measured entanglement entropy",
             color="white", fontsize=17, fontweight="bold")

CY, RD, GN, OR = "#00e5ff", "#ff0055", "#00ff88", "#ff9900"

# Panel 1: holographic screen
ax1 = fig.add_subplot(gs[0, 0]); ax1.set_facecolor("#111118")
ax1.scatter(xy[:, 0], xy[:, 1], s=6, color="#444455", alpha=0.5)
screen = plt.Circle((0, 0), 3.5, color=CY, fill=False, lw=2.5, ls="--")
ax1.add_patch(screen)
bits = ax1.scatter([], [], s=18, color=OR, zorder=5)
mass_m, = ax1.plot([], [], "o", color=RD, ms=11, mec="white", mew=1.2)
ax1.plot([0], [0], "X", color=OR, ms=13)
ax1.set(xlim=(-3.8, 3.8), ylim=(-3.8, 3.8), aspect="equal")
ax1.set_title("Holographic screen: N = 4·S$_{EE}$ bits", color=CY, fontsize=13)
ax1.tick_params(colors="#666666")
ax1.grid(False)
info1 = ax1.text(0.03, 0.95, "", transform=ax1.transAxes, color="white",
                 fontsize=11, va="top", family="monospace",
                 bbox=dict(facecolor="#1a1a1d", edgecolor="#333333"))

# Panel 2: measured area law
ax2 = fig.add_subplot(gs[0, 1]); ax2.set_facecolor("#111118")
ax2.plot(radii_all, 0.32 * (radii_all / a) ** 2, "--", color=GN, lw=1.5,
         label=r"Srednicki: $S = 0.32\,(R/a)^2$")
line_S, = ax2.plot([], [], "-", color=CY, lw=2.5, label="engine-measured $S_{EE}(R)$")
pt_S, = ax2.plot([], [], "o", color="white", ms=7)
ax2.set(xlim=(1.4, 3.6), ylim=(0, max(S_mono) * 1.1))
ax2.set_title(r"Area law: $S \propto R^{2.07}$,  $\kappa = 0.31$", color=OR, fontsize=13)
ax2.set_xlabel("screen radius R", color="#aaaaaa"); ax2.set_ylabel("S", color="#aaaaaa")
ax2.tick_params(colors="#666666"); ax2.grid(True, ls=":", alpha=0.3, color="#666666")
ax2.legend(loc="upper left", fontsize=9, facecolor="#1a1a1d", edgecolor="#333", labelcolor="#cccccc")

# Panel 3: screen temperature from equipartition
ax3 = fig.add_subplot(gs[1, 0]); ax3.set_facecolor("#111118")
line_T, = ax3.plot([], [], "-", color=OR, lw=2.5)
pt_T, = ax3.plot([], [], "o", color="white", ms=7)
ax3.set(xlim=(1.4, 3.6), ylim=(0, max(T_all) * 1.1))
ax3.set_title(r"Equipartition: $T(R) = 2Mc^2 / (N\,k_B)$", color=OR, fontsize=13)
ax3.set_xlabel("screen radius R", color="#aaaaaa"); ax3.set_ylabel("T(R)", color="#aaaaaa")
ax3.tick_params(colors="#666666"); ax3.grid(True, ls=":", alpha=0.3, color="#666666")

# Panel 4: emergent force vs Newton
ax4 = fig.add_subplot(gs[1, 1]); ax4.set_facecolor("#111118")
ax4.plot(radii_all, M_central * m_test / radii_all**2, "-", color=GN, lw=2, alpha=0.5,
         label=r"Newton $F \propto 1/R^2$")
line_F, = ax4.plot([], [], "-", color=RD, lw=3, label=r"emergent $F = T\,\nabla S$")
pt_F, = ax4.plot([], [], "o", color="white", ms=8, mec=RD, mew=2)
ax4.set(xlim=(1.4, 3.6), ylim=(0, M_central * m_test / 1.5**2 * 1.15))
ax4.set_title(r"Emergent force: measured exponent $-2.07$", color=OR, fontsize=13)
ax4.set_xlabel("distance R", color="#aaaaaa"); ax4.set_ylabel("F(R)", color="#aaaaaa")
ax4.tick_params(colors="#666666"); ax4.grid(True, ls=":", alpha=0.3, color="#666666")
ax4.legend(loc="upper right", fontsize=10, facecolor="#1a1a1d", edgecolor="#333", labelcolor="#cccccc")

# ---------------------------------------------------------------- animation
n_frames = 70
sweep = np.concatenate([np.linspace(24, 0, n_frames // 2).astype(int),
                        np.linspace(0, 24, n_frames // 2).astype(int)])

def update(f):
    i = sweep[f]
    R = radii_all[i]
    screen.set_radius(R)
    mass_m.set_data([R], [0])
    n_bits = max(6, int(6 + 54 * (N_all[i] / N_all.max())))
    th = np.linspace(0, 2 * np.pi, n_bits, endpoint=False)
    bits.set_offsets(np.c_[R * np.cos(th), R * np.sin(th)])
    info1.set_text(f"R = {R:.2f}\nN = {N_all[i]:6.2f} bits\nT = {T_all[i]:5.2f}")
    line_S.set_data(radii_all[:i + 1], S_mono[:i + 1]); pt_S.set_data([R], [S_mono[i]])
    line_T.set_data(radii_all[:i + 1], T_all[:i + 1]); pt_T.set_data([R], [T_all[i]])
    line_F.set_data(radii_all[:i + 1], F_all[:i + 1]); pt_F.set_data([R], [F_all[i]])
    return screen, mass_m, bits, info1, line_S, pt_S, line_T, pt_T, line_F, pt_F

ani = animation.FuncAnimation(fig, update, frames=n_frames, blit=True, interval=110)
out = "entropic_gravity_verified.gif"
print(f"Rendering {n_frames} frames -> {out}")
ani.save(out, writer="pillow", fps=12, dpi=90)
plt.close(fig)
print("Saved", out)
