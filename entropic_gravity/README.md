# Is Gravity Entropic? — Theory & Verification Package

This branch adds the entropic-gravity research program on top of the CEFE engine
(`main`): the unified-theory writeup, the verification paper, every script that
produced its numbers, and the animations.

**Paper:** `entropic_gravity_paper.tex` — *"Is Gravity Entropic? Numerical
Verification of the Entropic Gravity Hypothesis on Causal Diamond Lattices"*
(compiled PDF: `osacra_gravaser/reference papers/hgc_united.pdf`).

## Headline results (all measured, nothing fitted)

| # | Test | Theory | Measured |
|---|------|--------|----------|
| A | Area-law exponent | 2.000 | 2.09 ± 0.04 |
| A | Srednicki coefficient κ | 0.30 | 0.314 ± 0.007 |
| B | Emergent force exponent (measured S only) | −2.000 | −2.069 |
| B | G renormalization per scalar | O(few) (Susskind–Uglum) | 10.0 a² |
| C | Mutual information vs cutoff | finite | 0.097 / 0.130 / 0.107 |
| D | Fractal horizon D_H = 3/2 | — | **falsified** (α = 2.07) |
| E | Entanglement first law δS = δ⟨K⟩ | O(δ²) remainder | 3.7 % @ δm = 0.02 |
| F | Arrow of time: purity / coarse entropy | 0 drift / grows | 0 exactly / 1.48 → 3.28 |

## Reproduce

```bash
python scripts/prove_entropic_gravity.py     # tests A–D
python scripts/entanglement_first_law.py     # tests E–F
```

With engine **v0.2.0+** built (`pip install .`, Python 3.10/3.11), Tests A, B, D
and the first-law suite run directly on the engine's own Srednicki covariance
matrices via the `get_covariance_C` / `get_covariance_P` bindings, and both
scripts print an internal consistency cross-check (engine's
`compute_entanglement_entropy` vs exported matrices vs the NumPy port; the
first-law script also runs a 3D `QuantumWaveEngine` unitarity check). If the
extension is not built for your Python, everything falls back to a
line-identical NumPy port of the same construction — the paper's headline
numbers above are the port values; engine values agree within lattice-detail
tolerances. Test C (2D disk mutual-information study) and the 2D
arrow-of-time run stay on the port by design: the engine's causal-diamond
slices are intrinsically 3+1D. Requirements: `numpy`, `matplotlib`, `pillow`.

## Animations (`figures/`)

- `entropic_gravity_paper_walkthrough.gif` — the paper's argument in four
  narrated phases: area law → bits + equipartition → entropy gradient →
  emergent 1/R² force, live on the Schwarzschild lattice.
- `entropic_gravity_lattice_3d.gif` / `..._slice.gif` — the fluctuating scalar
  vacuum |φ(x,t)| on the Schwarzschild lattice (exact KG mode evolution),
  with horizon, photon sphere, holographic bits, and test mass.
- `entropic_gravity_3d.gif`, `entropic_gravity_verified.gif` — earlier
  schematic + measured-panels versions.

Regenerate with `python scripts/animate_entropic_gravity_<name>.py`.

## Documents

- `UNIFIED_THEORY.md` — how the three Osacra–Gravaser preprints merge into one
  axiomatic theory, and which parts are testable.
- `Osacra_Gravaser_theory_review.md` — line-by-line review of the theory
  papers, including the known open issues (photonic dispersion normalization,
  muon g−2 calibration, e^(−S_EE) ansatz stability).
- `CEFE_analysis_report.md` — engine-wide analysis and literature cross-check.
- `osacra_gravaser/` — the theory's phenomenology package (ghost-condensate
  matter coupling, backreaction) and the source preprints.
