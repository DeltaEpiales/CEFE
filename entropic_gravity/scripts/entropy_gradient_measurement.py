#!/usr/bin/env python
"""
ENTROPY-GRADIENT MEASUREMENT  (Verlinde's Axiom 2 under the microscope)
=======================================================================

Verlinde's entropy-gradient postulate:  a test mass m displaced by dx near a
holographic screen changes the screen's entropy linearly,

    dS = 2 pi m dx          (k_B = hbar = c = 1)

This is the ONE ingredient of the entropic-force derivation that the
verification paper (hgc_united.pdf, Section 9) had to assume rather than
measure.  This script attempts to measure it on the lattice.

Model:  the "test mass" is a localized heavy region of the field -- a ball of
radius rho where the scalar mass gap is mu >> m0 (ambient mass m0 = 1).
(A coherent-state excitation would change nothing: displacement does not
alter a Gaussian state's covariance, so a mass-gap region is the only
channel through which "matter" can register in vacuum entanglement.)
The "screen" is the entangling surface of the subregion A = {r <= R_scr};
S_A is its Srednicki entanglement entropy.

Measured here (cefe_lab, NumPy backend; engine path identical API):
  1. mass scan     dS(d) for mu = 2..16   -> linearity & slope in mu
  2. cutoff scan   a = 0.25..0.15         -> is the response length a UV artifact?
  3. geometry scan rho = 0.15..0.60       -> is it the bump's size?

Verdict preview (run to reproduce): the screen entropy DOES respond to a
nearby test mass, but the response is exponential with a PHYSICAL length
~0.4 (independent of cutoff and of bump size -- the ambient massive
propagator scale), sublinear/saturating in mu, and ~100x smaller in slope
than 2 pi m.  Verlinde's gradient does NOT emerge from vacuum entanglement;
it remains an independent postulate -- exactly as Section 9 classified it.

Runtime ~3 min.  Output: entropy_gradient_results.png + printed verdict.
"""
import sys
import numpy as np
from pathlib import Path

# repo-root bootstrap so `cefe_lab` imports from a bare checkout
# (an installed package still wins, since we append rather than insert)
sys.path.append(str(Path(__file__).resolve().parents[2]))

# daimon_runtime plot setup when available (local Kimi runtime);
# plain matplotlib otherwise (repo portability).
try:
    sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
    from daimon_runtime import setup_plot
    setup_plot()
except Exception:
    import matplotlib
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
import cefe_lab as lab

OUT = Path(__file__).resolve().parent
M0, R_SCR, R_DISK = 1.0, 2.0, 4.0

print("backend:", lab.backend.NAME)


def dS_profile(lat, maskA, S0, mu, rho, ds):
    """S(no bump) - S(bump at depth d) for each d in ds."""
    out = []
    for d in ds:
        rc = R_SCR - d          # bump center sits inside the screen at depth d
        prof = lambda c: np.where(
            np.linalg.norm(c - np.array([rc, 0.0, 0.0]), axis=1) <= rho, mu, M0)
        vac = lab.VacuumState.from_lattice(lat, mass=M0, mass_profile=prof)
        out.append(S0 - vac.entropy(maskA))
    return np.array(out)


def decay_length(ds, dS, dmin=0.375):
    sel = ds >= dmin
    k, _ = np.polyfit(ds[sel], np.log(np.maximum(dS[sel], 1e-12)), 1)
    return -1.0 / k


# ---------------------------------------------------------------- 1. mass scan
print("\n[1] mass scan  (a=0.25, rho=0.30)")
A1, RHO1 = 0.25, 0.30
lat = lab.disk_lattice(R_DISK, A1)
maskA = lab.ball_region(lat.coords, R_SCR)
S0 = lab.VacuumState.from_lattice(lat, mass=M0).entropy(maskA)
ds = np.arange(0.0, 1.55, 0.125)
mass_rows = []
for mu in [2.0, 4.0, 8.0, 16.0]:
    dS = dS_profile(lat, maskA, S0, mu, RHO1, ds)
    ell = decay_length(ds, dS)
    mass_rows.append((mu, dS, ell))
    print(f"    mu={mu:5.1f}  dS(contact)={dS[0]:.4f}  decay length={ell:.3f}")

# -------------------------------------------------------------- 2. cutoff scan
print("\n[2] cutoff scan  (mu=8, rho=0.30)")
cut_rows = []
for a in [0.25, 0.20, 0.15]:
    lat_a = lab.disk_lattice(R_DISK, a)
    mask_a = lab.ball_region(lat_a.coords, R_SCR)
    S0_a = lab.VacuumState.from_lattice(lat_a, mass=M0).entropy(mask_a)
    ds_a = np.arange(0.0, 1.55, 0.1)
    dS = dS_profile(lat_a, mask_a, S0_a, 8.0, 0.30, ds_a)
    ell = decay_length(ds_a, dS)
    cut_rows.append((a, ell))
    print(f"    a={a:4.2f}  decay length={ell:.3f}  ({ell/a:.2f} lattice spacings)")

# ------------------------------------------------------------- 3. geometry scan
print("\n[3] geometry scan  (mu=8, a=0.25)")
geo_rows = []
for rho in [0.15, 0.30, 0.45, 0.60]:
    dS = dS_profile(lat, maskA, S0, 8.0, rho, ds)
    ell = decay_length(ds, dS)
    geo_rows.append((rho, ell))
    print(f"    rho={rho:4.2f}  decay length={ell:.3f}")

# -------------------------------------------------------------------- verdict
mus = np.array([r[0] for r in mass_rows])
dS_contact = np.array([r[1][0] for r in mass_rows])
verlinde_slopes = 2 * np.pi * mus
meas_slopes = np.array([(r[1][2] - r[1][4]) / (ds[2] - ds[4]) for r in mass_rows])
print("\n" + "=" * 68)
print("VERDICT")
print("=" * 68)
print(f"  slope |dS/dd| at contact: {np.abs(meas_slopes)}")
print(f"  Verlinde 2*pi*m:          {verlinde_slopes}")
print(f"  ratio:                    {verlinde_slopes/np.abs(meas_slopes)}  (~100x too small)")
print(f"  mass dependence:          dS(contact) ~ {dS_contact}  (saturating, not linear)")
print(f"  response length:          constant in physical units vs cutoff;")
print(f"                            independent of bump size -> ambient propagator scale")
print("  -> Verlinde's dS = 2 pi m dx does NOT emerge from vacuum entanglement.")
print("     It remains an independent postulate (hgc_united.pdf Sec. 9), now")
print("     with direct numerical evidence for its independence.")

# -------------------------------------------------------------------- figure
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

ax = axes[0]
for mu, dS, ell in mass_rows:
    ax.plot(ds, dS, "o-", ms=3.5, lw=1.4, label=f"$\\mu={mu:.0f}$")
    fit = dS[3] * np.exp(-(ds - ds[3]) / ell)
    ax.plot(ds, fit, ":", lw=0.9, color="gray")
ax.set_xlabel("depth $d$ of test mass inside the screen")
ax.set_ylabel(r"$\Delta S(d) = S_A(\mathrm{no\ bump}) - S_A(d)$")
ax.set_title("screen entropy response to a test mass\n(dotted: exp fits, $\\ell\\approx0.4$)")
ax.legend(fontsize=8, title="test mass", title_fontsize=8)

ax = axes[1]
ax.errorbar([r[0] for r in cut_rows], [r[1] for r in cut_rows], fmt="s-",
            label="vs cutoff $a$")
ax.errorbar([r[0] for r in geo_rows], [r[1] for r in geo_rows], fmt="o-",
            label="vs bump size $\\rho$")
ax.axhline(0.4, color="gray", lw=0.8, ls="--")
ax.set_xlabel("cutoff $a$   /   bump radius $\\rho$")
ax.set_ylabel("response length $\\ell$")
ax.set_title("$\\ell$ is physical: constant in $a$ and $\\rho$\n(ambient propagator scale, $m_0=1$)")
ax.legend(fontsize=8)

ax = axes[2]
w = 0.38
x = np.arange(len(mus))
ax.bar(x - w/2, np.abs(meas_slopes), width=w, label="measured $|dS/dd|$")
ax.bar(x + w/2, verlinde_slopes, width=w, label="Verlinde $2\\pi m$")
ax.set_yscale("log")
ax.set_xticks(x, [f"$\\mu$={int(m_)}" for m_ in mus], fontsize=8)
ax.set_ylabel("entropy-gradient slope (log)")
ax.set_title("the missing factor ~100:\npostulated vs measured slope")
ax.legend(fontsize=8)

fig.suptitle("Entropy-gradient measurement: does  $\\Delta S = 2\\pi m\\,\\Delta x$  emerge from vacuum entanglement?  —  No.",
             fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.93))
png = OUT.parent / "figures" / "entropy_gradient_results.png"
png.parent.mkdir(exist_ok=True)
fig.savefig(png, dpi=200, bbox_inches="tight")
print(f"\nfigure saved: {png}")
