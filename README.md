# CEFE — Causal Entropic Field Engine

[![Python 3.8+](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Build Status](https://img.shields.io/badge/build-passing-brightgreen.svg)]()

A high-performance **quantum field theory and quantum mechanics simulation engine** built on C++ with Python bindings. CEFE computes entanglement entropy, real-time field dynamics, and quantum wave evolution on **causal diamond lattices** with support for curved spacetime backgrounds (Schwarzschild, Anti-de Sitter, Kerr, FLRW).

> **"Something physically accurate that doesn't take a supercomputer to run."**

## Key Features

- **Entanglement Entropy** — Exact Srednicki covariance matrix method for computing von Neumann entanglement entropy of arbitrary subregions
- **Area Law Verification** — Confirms S ~ r² scaling (Bekenstein-Hawking area law) on causal diamond lattices
- **Real Scalar λφ⁴ Theory** — Interacting QFT with Velocity Verlet symplectic integration and adaptive sub-stepping
- **Quantum Mechanics** — Schrödinger equation solved via Crank-Nicolson (unitary, probability-conserving)
- **Curved Spacetime** — Schwarzschild, Anti-de Sitter, Kerr, and FLRW metric backgrounds
- **Interactive 3D GUI** — PyQt6/OpenGL workbench for real-time simulation visualization
- **Publication-Ready Outputs** — matplotlib plots and GIF animations suitable for papers

## Quick Start

```bash
# Clone and install
git clone https://github.com/DeltaEpiales/CEFE.git
cd CEFE
pip install -e .

# Run the physics validation suite (all 6 experiments)
python run_final_experiments.py

# Run the complete demo suite (generates plots + animations)
python demos/run_all_demos.py

# Launch the interactive 3D GUI
python run_gui.py
```

## Python API

```python
import cefe_py as ce

# 1. Compute entanglement entropy on a causal diamond
grid = ce.CausalDiamondGrid(radius=5.0, spacing=0.8)
engine = ce.EntropicFieldEngine(grid)
engine.set_mass(1.0)
engine.evolve_to_slice(0.0)
S = engine.compute_entanglement_entropy(2.0)
print(f"Entanglement entropy: {S:.4f}")

# 2. Simulate interacting QFT (lambda-phi^4)
grid = ce.CausalDiamondGrid(radius=2.0, spacing=0.1)
qft = ce.LatticeQFTEngine(grid)
qft.set_mass(1.0)
qft.set_coupling(5.0)
qft.initialize_state()
qft.set_opposing_packets(1.0, 0.2, 2.0, 5.0)
qft.build_laplacian()
for _ in range(1000):
    qft.step_forward(0.005)
print(f"Energy: {qft.get_total_energy():.4f}")

# 3. Curved spacetime (Schwarzschild black hole)
metric = ce.SchwarzschildMetric(Rs=0.5)
grid = ce.CausalDiamondGrid(radius=4.0, spacing=0.8, metric=metric)
engine = ce.EntropicFieldEngine(grid)
engine.set_mass(1.0)
engine.evolve_to_slice(0.0)
S = engine.compute_entanglement_entropy(1.5)
print(f"Schwarzschild entropy: {S:.4f}")

# 4. High-level experiment API
from cefe_py.experiments import run_area_law_test
result = run_area_law_test()
print(f"Area law exponent: {result['alpha']:.3f} (expect ~2.0)")
```

## Three Simulation Engines

| Engine | Physics | Method |
|--------|---------|--------|
| `EntropicFieldEngine` | Entanglement entropy (Srednicki 1993) | Covariance matrix eigendecomposition |
| `LatticeQFTEngine` | Real scalar λφ⁴ field theory | Velocity Verlet symplectic integrator |
| `QuantumWaveEngine` | Schrödinger equation | Crank-Nicolson implicit scheme |

## Five Spacetime Metrics

| Metric | Description |
|--------|-------------|
| `MinkowskiMetric` | Flat spacetime (default) |
| `SchwarzschildMetric(Rs)` | Static black hole with horizon radius Rs |
| `AntiDeSitterMetric(L)` | AdS space with curvature radius L |
| `FLRWMetric` | Expanding universe (matter-dominated) |
| `KerrMetric(M, a)` | Rotating black hole |

## Physics Validation

The engine has been validated against foundational QFT results:

| Experiment | Result | Status |
|-----------|--------|--------|
| **Area Law** (S ~ r²) | α = 2.12, R² = 0.995 | ✅ Confirmed |
| **UV Scaling** (S ~ 1/ε²) | R² = 0.947 | ✅ Confirmed |
| **Curved Spacetime** | Schwarzschild↓ / AdS↑ | ✅ Confirmed |
| **Mass Gap** | S monotonically decreasing | ✅ Confirmed |
| **Thermal States** | S(T) ≥ S(0) | ✅ Confirmed |
| **Energy Conservation** | |δE/E| < 1% over 5000 steps | ✅ Confirmed |
| **Unit Tests** | 11/11 passing | ✅ All passed |

## Project Structure

```
CEFE/
├── cefe_py/               # Python package
│   ├── __init__.py        # Package API
│   ├── experiments.py     # High-level experiment runners
│   ├── gui.py             # PyQt6 interactive workbench
│   ├── cli.py             # Command-line simulator
│   └── render_engine.py   # 3D OpenGL renderer
├── include/cefe/          # C++ headers
│   ├── core/              # EntropicFieldEngine, LatticeQFTEngine, QuantumWaveEngine
│   └── geometry/          # CausalDiamond, Metric (Minkowski/Schwarzschild/AdS/Kerr/FLRW)
├── src/                   # C++ implementation
├── demos/                 # Demonstration scripts
│   ├── demo_qft_scattering.py
│   ├── demo_double_slit.py
│   ├── demo_entanglement_entropy.py
│   └── run_all_demos.py
├── tests/                 # pytest suite
├── run_final_experiments.py  # Definitive 6-experiment validation
├── run_gui.py             # GUI launcher
└── run_cli.py             # CLI launcher
```

## Command-Line Interface

```bash
# QFT simulation with 3D rendering
python run_cli.py --mode qft --steps 200 --render 3d

# Holographic gravity (headless)
cefe-cli --mode heg --steps 500 --render none

# QM double-slit with rendering
cefe-cli --mode qm --potential double_slit --render 3d
```

## Theoretical Foundation

CEFE implements the entanglement entropy calculation from:

- **Srednicki (1993)**: "Entropy and area" — establishes the area law for free scalar fields
- **Bombelli, Koul, Lee, Sorkin (1986)**: "Quantum source of entropy for black holes"
- **Wissner-Gross & Freer (2013)**: "Causal entropic forces"

The entanglement entropy is computed as:
$$S = \sum_i \left[(\nu_i + \tfrac{1}{2})\ln(\nu_i + \tfrac{1}{2}) - (\nu_i - \tfrac{1}{2})\ln(\nu_i - \tfrac{1}{2})\right]$$

where ν_i are the symplectic eigenvalues of the reduced covariance matrix C_A · P_A.

## Requirements

- Python ≥ 3.8
- CMake ≥ 3.14
- C++17 compiler (MSVC 2022, GCC 9+, Clang 10+)
- NumPy, Matplotlib
- **Optional**: PyQt6, pyqtgraph, PyOpenGL (for GUI)

## Citation

```bibtex
@software{cefe2026,
  author = {Osacra, Ryan},
  title = {CEFE: Causal Entropic Field Engine},
  year = {2026},
  url = {https://github.com/DeltaEpiales/CEFE}
}
```

## License

MIT License — see [LICENSE](LICENSE) for details.
