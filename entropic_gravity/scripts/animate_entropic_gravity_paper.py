"""
Entropic gravity: a walkthrough of the paper "Is Gravity Entropic?" (hgc_united.pdf)
in the lattice-dot visual idiom.

Narrative phases (matching the paper section by section):
  A. Prop. 1  - area law: entropy is carried by vacuum modes straddling the
                screen; only boundary sites light up.  S = kappa (R/a)^2.
  B. Axiom 2  - bits + equipartition: screen carries N = 4S bits at T = 2M/N.
  C. Axiom 2  - entropy gradient: test mass m displaced by dx removes entropy
                linearly, dS = 2 pi m dx; bits wink out as the screen contracts.
  D. Prop. 2  - emergent force: F = T grad S = G_eff M m / R^2; the mass
                accelerates inward while the measured force profile
                (exponent -2.069, from the engine's measured S(R)) is traced
                against Newton's -2.000.

Background: live KG vacuum fluctuations on the Schwarzschild lattice (same
construction as animate_entropic_gravity_lattice.py), kept dim so the
derivation elements dominate.

Quantitative anchors (from the verification runs reported in the paper):
  kappa = 0.314 (Srednicki coefficient, measured 0.314 +/- 0.007)
  c2 = 0.0250  ->  S(R) = 0.314 (R/a)^2 with a = 0.6
  G_eff = 10 a^2 per scalar (Corollary 1)
  measured force exponent -2.069 (R^2 = 0.9910) vs Newton -2.000
  M = 10, m = 1 (engine units)

Output: entropic_gravity_paper_walkthrough.gif  (160 frames, 1400x900)
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from PIL import Image
import io, os

OUT = r"C:\Users\kamos\Documents\CEFE"

# ---------- constants anchored to the paper ----------
A_LAT   = 0.6          # lattice spacing of the a=0.6 run
KAPPA   = 0.314        # measured Srednicki coefficient
M_SRC   = 10.0
M_TEST  = 1.0
G_EFF   = 10.0 * A_LAT**2          # Corollary 1: 10 a^2 per scalar
EXP_MEAS = -2.069                  # measured emergent exponent (Prop. 2)
RS      = 0.5                      # Schwarzschild radius of the background source

def S_area(R):                     # measured area-law profile S(R)
    return KAPPA * (R / A_LAT)**2
def N_bits(R):                     # Axiom 2: N = 4S
    return 4.0 * S_area(R)
def T_screen(R):                   # equipartition: T = 2M/N
    return 2.0 * M_SRC / N_bits(R)
def F_entropic(R):                 # F = T * dS/dx, dS/dx = 2 pi m
    return T_screen(R) * 2.0 * np.pi * M_TEST
def F_newton(R):                   # emergent Newtonian form
    return G_EFF * M_SRC * M_TEST / R**2

# ---------- Schwarzschild lattice vacuum background (engine construction) ---
H, RMAX, MASS2 = 0.45, 3.4, 0.05

def proper_radial(r1, r2, n=8):
    xs = np.linspace(r1, r2, n + 1)
    fx = 1.0 / np.sqrt(1.0 - RS / xs)
    h_ = (r2 - r1) / n
    return h_ / 3.0 * (fx[0] + fx[-1] + 4*fx[1:-1:2].sum() + 2*fx[2:-1:2].sum())

print("building lattice + diagonalizing ...")
grid, pts = {}, []
ns = int(np.ceil(RMAX / H)) + 1
for ix in range(-ns, ns+1):
    for iy in range(-ns, ns+1):
        for iz in range(-ns, ns+1):
            x, y, z = ix*H, iy*H, iz*H
            r = np.sqrt(x*x + y*y + z*z)
            if RS*1.02 < r <= RMAX:
                grid[(ix,iy,iz)] = len(pts); pts.append((x,y,z))
P = np.array(pts); N3 = len(P); RHO = np.linalg.norm(P, axis=1)

rows, cols, wts, DEG = [], [], [], np.zeros(N3)
for (ix,iy,iz), i in grid.items():
    for d in ((1,0,0),(0,1,0),(0,0,1)):
        j = grid.get((ix+d[0], iy+d[1], iz+d[2]))
        if j is None: continue
        ri, rj = RHO[i], RHO[j]
        if abs(ri-rj) > 1e-9:
            w = (H / proper_radial(min(ri,rj), max(ri,rj))) / H**2
        else:
            w = 1.0 / H**2
        rows += [i,j]; cols += [j,i]; wts += [-w,-w]
        DEG[i] += w; DEG[j] += w
L = np.zeros((N3,N3)); L[rows,cols] = wts; L[np.arange(N3),np.arange(N3)] = DEG
K = L + MASS2*np.eye(N3)
w2, U = np.linalg.eigh(K)
w2 = np.clip(w2, 1e-12, None); OM = np.sqrt(w2)
Amat = U / np.sqrt(2.0*OM)[None,:]
rng = np.random.default_rng(7)
XI1 = rng.standard_normal(N3); XI2 = rng.standard_normal(N3)
def field(t): return Amat @ (np.cos(OM*t)*XI1 + np.sin(OM*t)*XI2)
probe = np.concatenate([np.abs(field(t)) for t in (0.0, 2.1, 4.2)])
VMAX = np.percentile(probe, 97.5)
print(f"  {N3} sites, vmax={VMAX:.3f}")

# ---------- static scene elements ----------
u, v = np.mgrid[0:2*np.pi:22j, 0:np.pi:14j]
HX = RS*np.cos(u)*np.sin(v); HY = RS*np.sin(u)*np.sin(v); HZ = RS*np.cos(v)
def fib_sphere(n):
    k = np.arange(n) + 0.5
    th = np.arccos(1 - 2*k/n); ph = np.pi*(1+5**0.5)*k
    return np.c_[np.sin(th)*np.cos(ph), np.sin(th)*np.sin(ph), np.cos(th)]
FIB = fib_sphere(240)             # master bit pattern; first n(R) shown
TH_RING = np.linspace(0, 2*np.pi, 120)
BG = "#0b0b10"

# ---------- storyboard ----------
NFR = 160
PH_A, PH_B, PH_C = range(0, 40), range(40, 70), range(70, 100)
PH_D_START = 100
R_A = 2.60                        # screen radius during phase A,B
R_C_END = 2.40                    # end of slow displacement (phase C)
# phase D infall: integrate  R'' = -k/R^2 schematic, resample to 60 frames
k_acc = 0.55
rr, vv, tt, dt = R_C_END, 0.0, [0.0], 1e-3
traj = [rr]
while rr > 1.30:
    vv -= k_acc/rr**2 * dt; rr += vv*dt; traj.append(rr); tt.append(tt[-1]+dt)
traj = np.array(traj)
R_D = np.interp(np.linspace(0, len(traj)-1, NFR-PH_D_START), np.arange(len(traj)), traj)

def screen_R(f):
    if f < PH_C.start:  return R_A
    if f < PH_D_START:  return R_A + (R_C_END-R_A)*(f-PH_C.start+1)/len(PH_C)
    return R_D[f-PH_D_START]

CAPTIONS = {
 "A": ("1 | AREA LAW  (Prop. 1):  the entropy of the vacuum is carried by modes straddling the screen —\n"
       "       S(R) = κ (R/a)²,   κ = 0.314 measured  (Srednicki 1993: 0.30).  Only boundary sites light up."),
 "B": ("2 | BITS + EQUIPARTITION  (Axiom 2):  the screen carries N = 4S bits; the enclosed mass M\n"
       "       sets the screen temperature by ½ N T = M   →   T = 2M/N."),
 "C": ("3 | ENTROPY GRADIENT  (Axiom 2):  a test mass m displaced by Δx removes entropy linearly,\n"
       "       ΔS = 2π m Δx  —  the screen contracts, bits wink out, the survivors get hotter."),
 "D": ("4 | EMERGENT FORCE  (Prop. 2):  F = T·∇S = G_eff M m / R²,  G_eff = 10 a² per scalar.\n"
       "       Measured from the engine's entropy alone:  exponent −2.069  vs  Newton −2.000."),
}
def phase_of(f):
    if f in PH_A: return "A"
    if f in PH_B: return "B"
    if f < PH_D_START: return "C"
    return "D"

CHAIN = [
    ("S = κ (R/a)²", "area law — Prop. 1"),
    ("N = 4S,  T = 2M/N", "bits + equipartition — Ax. 2"),
    ("ΔS = 2π m Δx", "entropy gradient — Ax. 2"),
    ("F = T ∇S  ∝  1/R²", "emergent force — Prop. 2"),
]

def mass_pos(f, R):
    th = 0.9 + 2.4 * f / NFR          # slow angular drift while it falls
    return np.array([R*np.cos(th), R*np.sin(th), 0.35*R*np.sin(0.6*th)])

# force-curve data (precomputed, measured anchors)
R_GRID = np.linspace(1.30, 2.70, 200)
F_MEAS = F_entropic(R_A) * (R_GRID/R_A)**EXP_MEAS
F_NEWT = F_newton(R_A) * (R_GRID/R_A)**(-2.0)

# ---------- render ----------
print("rendering walkthrough ...")
frames = []
for f in range(NFR):
    ph = phase_of(f)
    R = screen_R(f)
    tphys = 9.0 * f / NFR
    phi = field(tphys)
    val = np.clip(np.abs(phi)/VMAX, 0, 1)

    fig = plt.figure(figsize=(14, 9), dpi=100, facecolor=BG)
    gs = gridspec.GridSpec(2, 2, width_ratios=[2.6, 1.15], height_ratios=[1.35, 1.0],
                           left=0.01, right=0.985, top=0.90, bottom=0.13, wspace=0.02, hspace=0.30)

    # ---- 3D scene ----
    ax = fig.add_subplot(gs[:, 0], projection="3d", facecolor=BG)
    dim = val**0.55
    bg_alpha = 0.35 if ph == "A" else 0.55
    ax.scatter(P[:,0], P[:,1], P[:,2], c=dim, cmap="magma", vmin=0, vmax=1,
               s=1.0 + 18*val**1.5, alpha=bg_alpha, linewidths=0, depthshade=False)

    # boundary-mode highlight (phase A): a thin shell of sites straddling the
    # screen -- the entropy lives HERE, on the surface, not in the volume
    if ph == "A":
        onb = np.abs(RHO - R) < 0.5*H
        ax.scatter(P[onb,0], P[onb,1], P[onb,2], c="#dff1ff", s=16,
                   alpha=0.95, linewidths=0, depthshade=False)

    # screen wireframe (faint)
    thw = np.linspace(0, 2*np.pi, 60)
    for zz, rr_ in ((0.0, R),):
        ax.plot(R*np.cos(thw), R*np.sin(thw), 0*thw, color="#59c2ff", lw=0.7, alpha=0.35)
        ax.plot(R*np.cos(thw), 0*thw, R*np.sin(thw), color="#59c2ff", lw=0.7, alpha=0.35)
        ax.plot(0*thw, R*np.cos(thw), R*np.sin(thw), color="#59c2ff", lw=0.7, alpha=0.35)

    # bits: count ∝ R², size/brightness ∝ T (heats up as screen contracts)
    if ph in ("B", "C", "D"):
        nshow = int(round(240 * (R/R_A)**2))
        Tn = T_screen(R)/T_screen(R_A)
        bsize = np.clip(6 * Tn, 4, 60)
        bits = FIB[:nshow] * R
        ax.scatter(bits[:,0], bits[:,1], bits[:,2], c="#59c2ff", s=bsize,
                   alpha=np.clip(0.45 + 0.12*Tn, 0, 0.95), linewidths=0, depthshade=False)

    # black hole + photon sphere
    ax.plot_surface(HX, HY, HZ, color="black", shade=False, zorder=5)
    ax.plot(RS*np.cos(TH_RING), RS*np.sin(TH_RING), 0*TH_RING, color=(0.6,0.4,0.9,0.35), lw=0.8)
    ax.plot(1.5*RS*np.cos(TH_RING), 1.5*RS*np.sin(TH_RING), 0*TH_RING,
            color="#ff7f0e", lw=1.0, alpha=0.8)

    # test mass
    if ph in ("C", "D"):
        mp = mass_pos(f, R)
        ax.scatter([mp[0]], [mp[1]], [mp[2]], c="#ff2d2d", s=90, depthshade=False, zorder=10)
        if ph == "C":   # displacement arrow toward the source
            dirn = -mp/np.linalg.norm(mp)
            ax.quiver(mp[0], mp[1], mp[2], dirn[0], dirn[1], dirn[2],
                      length=0.45, color="#ff2d2d", lw=1.6, arrow_length_ratio=0.35)

    # 3d cosmetics
    ax.set_facecolor(BG)
    for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
        pane.fill = False; pane.set_edgecolor((1,1,1,0.10))
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis._axinfo["grid"]["color"] = (1,1,1,0.08)
        axis.set_tick_params(colors=(1,1,1,0.30), labelsize=7)
    LIM = RMAX
    ax.set_xlim(-LIM,LIM); ax.set_ylim(-LIM,LIM); ax.set_zlim(-LIM,LIM)
    ax.set_box_aspect((1,1,1))
    ax.view_init(elev=20, azim=14 + 46*f/NFR)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])

    # live readout in the 3D corner
    readout = (f"R = {R:4.2f}   S = {S_area(R):5.2f}   N = {N_bits(R):5.1f} bits\n"
               f"T = {T_screen(R):5.3f}   F = T·∇S = {F_entropic(R):5.2f}")
    ax.text2D(0.02, 0.03, readout, transform=ax.transAxes, color="#9fd8ff",
              fontsize=10, family="monospace",
              bbox=dict(boxstyle="round,pad=0.45", fc=(0.05,0.07,0.12,0.85), ec="#59c2ff", alpha=0.9))

    # ---- derivation chain (top right) ----
    axc = fig.add_subplot(gs[0, 1], facecolor=BG); axc.axis("off")
    axc.set_title("the Verlinde loop — hgc_united.pdf §3–4", color=(1,1,1,0.6), fontsize=10, pad=8)
    active = {"A":0, "B":1, "C":2, "D":3}[ph]
    for i, (eq, lab) in enumerate(CHAIN):
        on = (i == active)
        col = "#ffffff" if on else (1,1,1,0.30)
        axc.text(0.03, 0.82 - 0.24*i, ("● " if on else "○ ") + eq,
                 transform=axc.transAxes, fontsize=13 if on else 11.5,
                 color="#ffd27f" if on else col, family="monospace",
                 fontweight="bold" if on else "normal")
        axc.text(0.09, 0.82 - 0.24*i - 0.075, lab, transform=axc.transAxes,
                 fontsize=9, color=(1,1,1,0.75) if on else (1,1,1,0.28))
        if i < 3:
            axc.text(0.05, 0.82 - 0.24*i - 0.145, "↓", transform=axc.transAxes,
                     fontsize=11, color=(1,1,1,0.25))

    # ---- force trace (bottom right) ----
    axf = fig.add_subplot(gs[1, 1], facecolor=BG)
    axf.plot(R_GRID, F_NEWT, color=(1,1,1,0.45), lw=1.2, ls="--",
             label="Newton 1/R² (G_eff = 10a²)")
    axf.plot(R_GRID, F_MEAS, color="#ff5252", lw=1.0, alpha=0.35)
    if ph == "D":
        upto = R_GRID >= R
        axf.plot(R_GRID[upto], F_MEAS[upto], color="#ff5252", lw=2.4,
                 label="measured: F ∝ R$^{-2.069}$")
        axf.scatter([R], [F_entropic(R_A)*(R/R_A)**EXP_MEAS], c="#ff2d2d", s=55, zorder=5)
    else:
        axf.plot([], [], color="#ff5252", lw=2.4, label="measured: F ∝ R$^{-2.069}$")
    axf.set_xlim(2.75, 1.20)   # R decreases left-to-right as the mass falls
    axf.set_ylim(0, F_NEWT.max()*1.10)
    axf.set_xlabel("screen radius R", color=(1,1,1,0.6), fontsize=9)
    axf.set_ylabel("force on test mass", color=(1,1,1,0.6), fontsize=9)
    axf.tick_params(colors=(1,1,1,0.4), labelsize=8)
    for s in axf.spines.values(): s.set_color((1,1,1,0.2))
    leg = axf.legend(loc="upper left", fontsize=8, frameon=False)
    for t_ in leg.get_texts(): t_.set_color((1,1,1,0.75))
    axf.set_title("emergent 1/R² from measured entropy alone",
                  color=(1,1,1,0.6), fontsize=9.5, pad=6)

    # titles + caption
    fig.suptitle("IS GRAVITY ENTROPIC?  —  the derivation, live on the CEFE lattice",
                 color="white", fontsize=15, y=0.965)
    fig.text(0.02, 0.055, CAPTIONS[ph], color="#cfe8ff", fontsize=11.5, family="monospace")
    fig.text(0.02, 0.012,
             "background: fluctuating scalar vacuum on the Schwarzschild lattice (|φ(x,t)| per site)   ·   "
             "black sphere: horizon Rₛ   ·   orange: photon sphere   ·   cyan: holographic screen + bits   ·   red: test mass m",
             color=(1,1,1,0.40), fontsize=8.5)

    buf = io.BytesIO(); fig.savefig(buf, format="png", facecolor=BG); plt.close(fig)
    buf.seek(0)
    frames.append(Image.open(buf).convert("P", palette=Image.ADAPTIVE, colors=256))
    if f % 20 == 0: print(f"  frame {f}/{NFR}")
    if f in (20, 55, 85, 130):
        frames[-1].convert("RGB").save(os.path.join(OUT, f"_walk_check_{ph}.png"))

out = os.path.join(OUT, "entropic_gravity_paper_walkthrough.gif")
frames[0].save(out, save_all=True, append_images=frames[1:], duration=85, loop=0, optimize=True)
print("saved", out, os.path.getsize(out)//2**20, "MB")
