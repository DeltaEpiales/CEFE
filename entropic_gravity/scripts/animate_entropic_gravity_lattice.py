"""
Entropic gravity in the 'atomic_energy_model' visual idiom:
dark scene, the lattice itself rendered as discrete dots colored by the
live quantum field value.

Physics (same engine construction as the verification suite):
  - Cartesian lattice of sites in the spherical shell  Rs < r < Rmax
    around a Schwarzschild source (Rs = 0.5).
  - Radial links weighted by proper distance  l_p = int dr/sqrt(1-Rs/r)
    (midpoint Simpson) -- the Metric.hpp construction ported to NumPy.
  - Graph Laplacian K = L + m0^2, diagonalized once: K = U w^2 U^T.
  - Exact KG vacuum time evolution:
        phi(t) = U diag(1/sqrt(2 w)) [ cos(w t) xi1 + sin(w t) xi2 ],
    xi1, xi2 iid N(0,1)  -> every frame is a bona fide snapshot of the
    fluctuating quantum vacuum on curved space.
  - Dot color/size = |phi(x,t)|  (magma), the hole in the middle IS the
    black hole; holographic bits tile the screen at R_screen; a test
    particle rides the orbit the entropic force produces.

Outputs:
  entropic_gravity_lattice_3d.gif    (full 3D, rotating camera)
  entropic_gravity_lattice_slice.gif (z=0 slice, like the 2D reference)
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from PIL import Image
import io, os, time

OUT = r"C:\Users\kamos\Documents\CEFE"

# ---------------- lattice + Schwarzschild metric ----------------
RS    = 0.5      # Schwarzschild radius
RMAX  = 3.4
H     = 0.40     # lattice spacing
MASS2 = 0.05     # small IR mass gap (regulates zero mode)
R_SCREEN = 1.8   # holographic screen radius
M_SRC = 10.0     # source mass (engine units, Test-B normalization)

def proper_radial(r1, r2, n=8):
    """Proper length of a radial link: int sqrt(g_rr) dr, Simpson."""
    xs = np.linspace(r1, r2, n + 1)
    fx = 1.0 / np.sqrt(1.0 - RS / xs)
    h_ = (r2 - r1) / n
    return h_ / 3.0 * (fx[0] + fx[-1] + 4 * fx[1:-1:2].sum() + 2 * fx[2:-1:2].sum())

print("building lattice ...")
t0 = time.time()
grid = {}
pts = []
n_side = int(np.ceil(RMAX / H)) + 1
for ix in range(-n_side, n_side + 1):
    for iy in range(-n_side, n_side + 1):
        for iz in range(-n_side, n_side + 1):
            x, y, z = ix * H, iy * H, iz * H
            r = np.sqrt(x * x + y * y + z * z)
            if RS * 1.02 < r <= RMAX:
                grid[(ix, iy, iz)] = len(pts)
                pts.append((x, y, z))
P = np.array(pts)
N = len(P)
RHO = np.linalg.norm(P, axis=1)
print(f"  {N} sites  ({time.time()-t0:.1f}s)")

# neighbor links (6-connected), weights
rows, cols, wts = [], [], []
DEG = np.zeros(N)
for (ix, iy, iz), i in grid.items():
    for d in ((1,0,0),(0,1,0),(0,0,1)):
        j = grid.get((ix+d[0], iy+d[1], iz+d[2]))
        if j is None:
            continue
        ri, rj = RHO[i], RHO[j]
        if abs(ri - rj) > 1e-9:                    # radial link -> proper weight
            lp = proper_radial(min(ri, rj), max(ri, rj))
            w = (H / lp) / H**2
        else:                                      # tangential link
            w = 1.0 / H**2
        rows += [i, j]; cols += [j, i]; wts += [-w, -w]
        DEG[i] += w; DEG[j] += w
L = np.zeros((N, N))
L[rows, cols] = wts
L[np.arange(N), np.arange(N)] = DEG
K = L + MASS2 * np.eye(N)

print("diagonalizing ...")
t0 = time.time()
w2, U = np.linalg.eigh(K)
w2 = np.clip(w2, 1e-12, None)
OM = np.sqrt(w2)
print(f"  done ({time.time()-t0:.1f}s),  omega in [{OM[0]:.3f}, {OM[-1]:.3f}]")

# vacuum covariance factor:  phi(t) = A [cos(wt) xi1 + sin(wt) xi2]
A = U / np.sqrt(2.0 * OM)[None, :]
rng = np.random.default_rng(7)
XI1 = rng.standard_normal(N)
XI2 = rng.standard_normal(N)

def field(t):
    return A @ (np.cos(OM * t) * XI1 + np.sin(OM * t) * XI2)

# fixed global color scale from the 97th percentile over a few probes
probe = np.concatenate([np.abs(field(t)) for t in (0.0, 2.1, 4.2, 6.3)])
VMAX = np.percentile(probe, 97.5)
print(f"color scale vmax = {VMAX:.3f}")

# holographic screen bits (Fibonacci sphere)
def fib_sphere(n, R):
    k = np.arange(n) + 0.5
    th = np.arccos(1 - 2 * k / n)
    ph = np.pi * (1 + 5**0.5) * k
    return np.c_[R*np.sin(th)*np.cos(ph), R*np.sin(th)*np.sin(ph), R*np.cos(th)]
BITS = fib_sphere(240, R_SCREEN)

# horizon sphere mesh (low res, black)
u, v = np.mgrid[0:2*np.pi:26j, 0:np.pi:16j]
HS = np.c_[np.cos(u.ravel())*np.sin(v.ravel()),
           np.sin(u.ravel())*np.sin(v.ravel()),
           np.cos(v.ravel())] * RS
HX = HS[:,0].reshape(u.shape); HY = HS[:,1].reshape(u.shape); HZ = HS[:,2].reshape(u.shape)

BG = "#0b0b10"

def style_3d(ax, azim, elev=22):
    ax.set_facecolor(BG)
    for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
        pane.fill = False; pane.set_edgecolor((1,1,1,0.12))
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis._axinfo["grid"]["color"] = (1,1,1,0.10)
        axis._axinfo["tick"]["color"] = (1,1,1,0.25)
        axis.set_tick_params(colors=(1,1,1,0.35), labelsize=7)
        axis.label.set_color((1,1,1,0.5))
    LIM = RMAX
    ax.set_xlim(-LIM, LIM); ax.set_ylim(-LIM, LIM); ax.set_zlim(-LIM, LIM)
    ax.set_box_aspect((1,1,1))
    ax.view_init(elev=elev, azim=azim)
    ax.set_xlabel("x", fontsize=8); ax.set_ylabel("y", fontsize=8); ax.set_zlabel("z", fontsize=8)

# ---------------- frame renderers ----------------
NFR = 120
T_MAX = 9.0
TIMES = np.linspace(0, T_MAX, NFR)

def draw_common_3d(ax, phi, iframe, azim):
    val = np.clip(np.abs(phi) / VMAX, 0, 1)
    col = val**0.55                       # gamma lift: mid-tones visible
    sz  = val**1.5
    ax.scatter(P[:,0], P[:,1], P[:,2],
               c=col, cmap="magma", vmin=0, vmax=1,
               s=1.5 + 34*sz, alpha=0.9, linewidths=0, depthshade=False)
    # black hole (black sphere + faint rim so the void reads)
    ax.plot_surface(HX, HY, HZ, color="black", shade=False, alpha=1.0, zorder=5)
    th = np.linspace(0, 2*np.pi, 120)
    ax.plot(RS*np.cos(th), RS*np.sin(th), 0*th, color=(0.6,0.4,0.9,0.35), lw=0.9)
    # photon sphere ring
    ax.plot(1.5*RS*np.cos(th), 1.5*RS*np.sin(th), 0*th,
            color="#ff7f0e", lw=1.0, alpha=0.8)
    # holographic bits
    ax.scatter(BITS[:,0], BITS[:,1], BITS[:,2], c="#59c2ff", s=6,
               alpha=0.8, linewidths=0, depthshade=False)
    # test particle on circular orbit at R_SCREEN (equatorial)
    thp = 2*np.pi * 2.2 * iframe / NFR
    ax.scatter([R_SCREEN*np.cos(thp)], [R_SCREEN*np.sin(thp)], [0],
               c="#ff2d2d", s=55, depthshade=False, zorder=10)
    style_3d(ax, azim)

def render_3d():
    print("rendering 3D gif ...")
    frames = []
    for i, t in enumerate(TIMES):
        phi = field(t)
        fig = plt.figure(figsize=(10, 10), dpi=100, facecolor=BG)
        ax = fig.add_subplot(111, projection="3d", facecolor=BG)
        draw_common_3d(ax, phi, i, azim=18 + 45 * i / NFR)
        fig.suptitle("Entropic Gravity — Scalar Vacuum on Schwarzschild Lattice (3D)",
                     color="white", fontsize=13, y=0.96)
        ax.set_title(f"t = {t:4.1f}   |   dots = |φ(x,t)| of the quantum vacuum   |   "
                     f"blue = holographic bits at R = {R_SCREEN}",
                     color=(1,1,1,0.55), fontsize=9, pad=2)
        fig.tight_layout(rect=(0, 0, 0, 0.94))
        buf = io.BytesIO()
        fig.savefig(buf, format="png", facecolor=BG)
        plt.close(fig)
        buf.seek(0)
        frames.append(Image.open(buf).convert("P", palette=Image.ADAPTIVE, colors=256))
        if i % 20 == 0: print(f"  frame {i}/{NFR}")
        if i == NFR // 2:
            frames[-1].convert("RGB").save(os.path.join(OUT, "_lat3d_check.png"))
    out = os.path.join(OUT, "entropic_gravity_lattice_3d.gif")
    frames[0].save(out, save_all=True, append_images=frames[1:],
                   duration=80, loop=0, optimize=True)
    print("  saved", out, os.path.getsize(out)//2**20, "MB")

def render_slice():
    print("rendering slice gif ...")
    sel = np.abs(P[:,2]) < H/2 * 1.01
    Q = P[sel]
    idx = np.where(sel)[0]
    frames = []
    for i, t in enumerate(TIMES):
        phi = field(t)[idx]
        val = np.clip(np.abs(phi) / VMAX, 0, 1)
        col = val**0.55
        fig, ax = plt.subplots(figsize=(8, 8), dpi=100, facecolor=BG)
        ax.set_facecolor(BG)
        ax.scatter(Q[:,0], Q[:,1], c=col, cmap="magma", vmin=0, vmax=1,
                   s=18 + 55*val**1.5, alpha=0.95, linewidths=0)
        # horizon + photon sphere + screen + particle
        ax.add_patch(plt.Circle((0,0), RS, color="black", zorder=5))
        ax.add_patch(plt.Circle((0,0), RS, fill=False, ec=(1,1,1,0.25), lw=0.8, zorder=6))
        th = np.linspace(0, 2*np.pi, 200)
        ax.plot(1.5*RS*np.cos(th), 1.5*RS*np.sin(th), color="#ff7f0e", lw=1.0, alpha=0.8)
        ax.plot(R_SCREEN*np.cos(th), R_SCREEN*np.sin(th), color="#59c2ff",
                lw=0.9, alpha=0.6, ls="--")
        # holographic bits on the screen circle
        thb = np.linspace(0, 2*np.pi, 60, endpoint=False)
        ax.scatter(R_SCREEN*np.cos(thb), R_SCREEN*np.sin(thb),
                   c="#59c2ff", s=10, alpha=0.85, zorder=8, linewidths=0)
        thp = 2*np.pi * 2.2 * i / NFR
        ax.scatter([R_SCREEN*np.cos(thp)], [R_SCREEN*np.sin(thp)],
                   c="#ff2d2d", s=70, zorder=10)
        ax.set_xlim(-RMAX, RMAX); ax.set_ylim(-RMAX, RMAX)
        ax.set_aspect("equal")
        ax.tick_params(colors=(1,1,1,0.45), labelsize=8)
        for s in ax.spines.values(): s.set_color((1,1,1,0.2))
        ax.set_xlabel("x", color=(1,1,1,0.6)); ax.set_ylabel("y", color=(1,1,1,0.6))
        ax.set_title(f"Entropic Gravity — Quantum Vacuum Fluctuations (z = 0 slice)   t = {t:4.1f}",
                     color="white", fontsize=11, pad=10)
        fig.tight_layout()
        buf = io.BytesIO()
        fig.savefig(buf, format="png", facecolor=BG)
        plt.close(fig)
        buf.seek(0)
        frames.append(Image.open(buf).convert("P", palette=Image.ADAPTIVE, colors=256))
        if i % 20 == 0: print(f"  frame {i}/{NFR}")
        if i == NFR // 2:
            frames[-1].convert("RGB").save(os.path.join(OUT, "_latslice_check.png"))
    out = os.path.join(OUT, "entropic_gravity_lattice_slice.gif")
    frames[0].save(out, save_all=True, append_images=frames[1:],
                   duration=80, loop=0, optimize=True)
    print("  saved", out, os.path.getsize(out)//2**20, "MB")

render_3d()
render_slice()
print("done.")
