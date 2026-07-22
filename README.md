# Causal Entropic Field Engine (CEFE)

CEFE is a high-performance C++ physics engine (with Python bindings) designed to model the causal entropic field framework, bringing together concepts from causal diamonds, holographic entanglement entropy, and quantum field mechanics.

## Features

- **Advanced Spacetime Geometries:** Simulate scalar fields on Minkowski, Anti-de Sitter (AdS), FLRW (Cosmological), and Kerr (Rotating Black Hole) metrics.
- **Thermal Mixed States:** Initialize the field in a finite-temperature thermal state $T>0$ to model heat baths and Hawking radiation approximations.
- **Interacting Quantum Fields:** Supports $\lambda \phi^4$ interaction coupling with dynamic time-evolution integrators.
- **Arbitrary Entangling Surfaces:** Compute the entanglement entropy of any arbitrary subregion (e.g., tori, disconnected regions, half-spaces) by passing custom node indices.
- **Scalable Sparse Eigensolvers:** Bypasses dense $O(N^3)$ matrix inversion bottlenecks by utilizing the `Spectra` library to compute only the lowest-energy modes via iterative shift-and-invert techniques.
- **High Performance:** Core algorithms are heavily parallelized in C++ using `Eigen` and `OpenMP`.

## Installation

### Prerequisites
- Python 3.7+
- CMake 3.14+
- A C++17 compatible compiler (e.g., GCC, Clang, MSVC)

### Build from Source

Clone the repository and install it directly via `pip`:

```bash
git clone <repository_url>
cd CEFE
pip install .
```

This will automatically fetch dependencies (Eigen, pybind11, Spectra) via CMake and compile the C++ extensions.

## Usage

Here is a simple example showing how to initialize the engine in Python using the new advanced features:

```python
import cefe_py as ce

# 1. Initialize an Anti-de Sitter (AdS) geometry
ads_metric = ce.AntiDeSitterMetric(L=10.0)
grid = ce.CausalDiamondGrid(radius=2.0, spacing=1.0, metric=ads_metric)

# 2. Setup the engine with Interacting Fields (lambda phi^4)
engine = ce.EntropicFieldEngine(grid)
engine.set_mass(1.0)
engine.set_interaction_coupling(lambda_val=0.1)

# 3. Use the Scalable Sparse Eigensolver to compute a Thermal State (T=5.0)
# This avoids the O(N^3) memory crash by only solving for the lowest 10 modes
engine.evolve_to_slice_sparse(target_t=0.0, num_modes=10, temperature=5.0)

# 4. Compute Entanglement Entropy of an Arbitrary Subregion
# Example: computing entropy of a disconnected checkerboard pattern
indices = [i for i in range(engine.get_state_size()) if i % 2 == 0]
entropy = engine.compute_entanglement_entropy_indices(indices)
print(f"Entanglement Entropy: {entropy}")
```

## Testing

Run tests with `pytest`:
```bash
pip install pytest
pytest tests/
python test_simulation.py
```
