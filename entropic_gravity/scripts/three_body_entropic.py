"""
Three-body problem: Newtonian gravity vs the engine-measured entropic force.

Left panel:  the Chenciner-Montgomery figure-8 choreography under exact
             Newtonian gravity (F ~ R^-2).  Closed, periodic, stable for
             tens of periods.
Right panel: identical initial conditions under the entropic force law
             measured from the CEFE engine's entropy profile:
             F ~ R^-2.069 (Test B of prove_entropic_gravity.py), normalized
             to agree with Newton at R = 1.  Bertrand's theorem (only
             beta = -2 and beta = +1 close all orbits) says the figure-8
             cannot stay closed -- watch it precess and slowly unwind.
             The 0.069 deviation is the lattice-scale (a = 0.6) correction;
             it vanishes as a -> 0 (alpha -> 2 as the cutoff refines).

Integrator: velocity Verlet, shared for both panels.  G = m = 1 units,
tiny Plummer softening for close encounters.

Output: three_body_entropic.gif  (140 frames, 1200x640)
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
import io, os

OUT = r"C:\Users\kamos\Documents\CEFE"

# Chenciner-Montgomery figure-8 initial conditions (G = m = 1)
P0 = np.array([[ 0.97000436, -0.24308753],
               [-0.97000436,  0.24308753],
               [ 0.0,          0.0       ]])
V3 = np.array([-0.93240737, -0.86473146])
V0 = np.array([-V3/2, -V3/2, V3])
PERIOD = 6.3259

BETA_NEWTON = -2.000
BETA_MEAS   = -2.069      # engine-measured entropic exponent (R^2 = 0.9910)
EPS2 = 1e-6               # Plummer softening^2

def accel(pos, beta):
    acc = np.zeros_like(pos)
    for i in range(3):
        for j in range(3):
            if i == j:
                continue
            d = pos[j] - pos[i]
            r2 = d @ d + EPS2
            acc[i] += d * r2 ** (beta/2 - 0.5)   # G m_j r^beta * unit vector
    return acc

def integrate(beta, n_steps, dt):
    pos, vel = P0.copy(), V0.copy()
    traj = np.empty((n_steps, 3, 2))
    a = accel(pos, beta)
    for s in range(n_steps):
        pos += vel * dt + 0.5 * a * dt * dt
        a_new = accel(pos, beta)
        vel += 0.5 * (a + a_new) * dt
        a = a_new
        traj[s] = pos
    return traj

DT = 2e-4 * PERIOD
N_PERIODS = 4
NSTEPS = int(N_PERIODS * PERIOD / DT)
print(f"integrating {NSTEPS} steps x 2 force laws ...")
TR_N = integrate(BETA_NEWTON, NSTEPS, DT)
TR_E = integrate(BETA_MEAS,   NSTEPS, DT)

BG = "#0b0b10"
BODY_COLORS = ["#ffd166", "#ef476f", "#06d6a0"]
NFR = 140
FRAMES_IDX = np.linspace(0, NSTEPS - 1, NFR).astype(int)
TRAIL = 2600   # trail length in integrator steps

print("rendering ...")
frames = []
for fi, s in enumerate(FRAMES_IDX):
    fig, axes = plt.subplots(1, 2, figsize=(12, 6.4), dpi=100, facecolor=BG)
    for ax, traj, beta, title in (
        (axes[0], TR_N, BETA_NEWTON, "NEWTON:  F ∝ R$^{-2.000}$"),
        (axes[1], TR_E, BETA_MEAS,  f"ENTROPIC (engine-measured):  F ∝ R$^{{{BETA_MEAS:.3f}}}$"),
    ):
        ax.set_facecolor(BG)
        s0 = max(0, s - TRAIL)
        for b in range(3):
            seg = traj[s0:s + 1, b]
            n = len(seg)
            for k in range(1, n, 6):     # sparse segments, alpha grows toward head
                ax.plot(seg[k-1:k+1, 0], seg[k-1:k+1, 1],
                        color=BODY_COLORS[b], lw=1.4, alpha=0.06 + 0.5 * k / n,
                        solid_capstyle="round")
            ax.scatter(*traj[s, b], c=BODY_COLORS[b], s=90, zorder=5,
                       edgecolors="white", linewidths=0.6)
        ax.set_xlim(-1.35, 1.35); ax.set_ylim(-1.05, 1.05)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(title, color="white", fontsize=12, pad=8)
        ax.text(0.02, 0.03, f"t = {s*DT:6.2f}   ({s*DT/PERIOD:4.2f} periods)",
                transform=ax.transAxes, color=(1, 1, 1, 0.55), fontsize=9,
                family="monospace")
    fig.suptitle("THREE-BODY PROBLEM — the figure-8 choreography under two force laws",
                 color="white", fontsize=13, y=0.97)
    fig.text(0.5, 0.030,
             "same initial conditions, same integrator  ·  exponent −2.069 measured from CEFE entropy (Test B)",
             ha="center", color=(1, 1, 1, 0.45), fontsize=8.5)
    fig.text(0.5, 0.006,
             "Bertrand's theorem: only 1/R² closes the 8 — the measured law precesses, and the precession → 0 as the lattice a → 0",
             ha="center", color=(1, 1, 1, 0.45), fontsize=8.5)
    fig.tight_layout(rect=(0, 0.05, 1, 0.94))
    buf = io.BytesIO()
    fig.savefig(buf, format="png", facecolor=BG)
    plt.close(fig)
    buf.seek(0)
    frames.append(Image.open(buf).convert("P", palette=Image.ADAPTIVE, colors=256))
    if fi % 35 == 0:
        print(f"  frame {fi}/{NFR}")
    if fi == NFR // 2:
        frames[-1].convert("RGB").save(os.path.join(OUT, "_3body_check.png"))

out = os.path.join(OUT, "three_body_entropic.gif")
frames[0].save(out, save_all=True, append_images=frames[1:], duration=70, loop=0, optimize=True)
print("saved", out, os.path.getsize(out) // 2**20, "MB")
