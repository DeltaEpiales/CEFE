# Osacra-Gravaser Framework

This repository contains the phenomenological mathematical models and physical physics simulations for the **Holographic Emergent Gravity (HEG)** and **Homotopy Type Theory (HoTT)** thermodynamics framework.

By bridging the gap between fractional discrete mathematics and raw scalar field physics, this framework demonstrates that gravity, mass, and time natively emerge from topological backreaction across the cosmological event horizon.

## Features
- **Holographic Emergent Gravity:** Models Muon $g-2$ anomalies, Cosmic Birefringence rotation angles, and photonic dispersion splitting bounded strictly by observational parameters (Planck 2018).
- **HoTT Formalization:** Models the fundamental arrow of time purely as the vectorized concatenation of non-local entanglement paths.
- **Topological Knot Tensor:** Computes the fractional finite differences of the quadratic polarization density matrix across discrete $D_H = 3/2$ lattice structures.
## Rigorous Theoretical Foundation
Unlike standard curve-fitted toy models, OSACRA-GRAVASER is mathematically anchored to strict quantum field theoretic principles:
- **EFT Lagrangian Integration:** The macroscopic phenomenologies (Muon $g-2$, Cosmic Birefringence) are derived directly from a Dark Sector Effective Field Theory (EFT) Hamiltonian deformation: $\mathcal{L}_{int} = \frac{c_\mu}{M} (\partial_\nu \phi) \bar{\psi} \gamma^\nu \gamma^5 \psi$.
- **Renormalization Group (RG) Flow:** The framework bridges the gap between the $2.4 \text{ meV}$ Dark Energy scale and $1 \text{ GeV}$ Hadronic scale natively via a one-loop logarithmic beta-function RG flow, deriving macroscopic rest mass dynamically.
- **Trace-Preserving MERA:** All density matrix evolution across the fractal lattice strictly enforces completely positive trace-preserving (CPTP) maps, conserving total quantum probability ($\text{Tr}(\rho) = 1$).
- **Fractional Calculus Justification:** The Knot Tensor uses a mathematically rigorous fractional derivative scaling order $\alpha = D_H / 2.0 = 0.75$, derived directly from Riemann-Liouville scaling rather than relying on arbitrary fudge factors.

### The Adiabatic Backreaction Approximation
Because the C++ engine (`cefe_py`) utilizes a fixed topological boundary (e.g., Anti-de Sitter metrics), the Python-based thermodynamic backreaction relies on an **Adiabatic Approximation** (similar to the Born-Oppenheimer approximation). The spacetime metric is treated as a fixed background relative to the rapid thermodynamic fluctuations of the scalar field over any $\Delta t$. While energy backreaction thermalization is exact, dynamically updating the metric tensor $g_{\mu\nu}$ per integration step requires an external active Einstein-solver.

## Installation

```bash
git clone <repository_url>
cd osacra_gravaser
pip install -r requirements.txt
```
*(Requires the C++ CEFE engine `cefe_py` to be installed in your environment to run visual simulations)*

## Running

Run theoretical unit tests for phenomenological verification:
```bash
python tests/test_phenomenology.py
```

Run physical graphical simulations using the C++ CEFE backend:
```bash
python simulations/visual_sims.py
```
