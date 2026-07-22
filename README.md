# Causal Entropic Field Engine (CEFE)

CEFE is an engine that models the causal entropic field framework, bringing together concepts from causal diamonds, entropic bounds, and field mechanics.

## Features

- **Causal Diamond Grids:** Calculates the geometry of causal diamonds and Bekenstein bounds.
- **Entropic Field Mechanics:** Computes degrees of freedom and state sizes within localized fields.
- **High Performance:** Core algorithms are implemented in C++ (with Eigen) for maximum performance.
- **Python Bindings:** Provides a seamless interface to use the engine directly from Python using `pybind11`.

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

This will automatically invoke CMake and compile the C++ extensions.

## Usage

Here is a simple example showing how to initialize the engine in Python:

```python
import cefe_core

# Initialize the causal diamond and compute Bekenstein bounds
diamond = cefe_core.CausalDiamond(2.0, 1.0)
print(f"Degrees of Freedom: {diamond.get_degrees_of_freedom()}")
print(f"Bekenstein Bound: {diamond.get_bekenstein_bound()}")

# Run the entropic field engine
engine = cefe_core.EntropicFieldEngine()
dof = diamond.get_degrees_of_freedom()
engine.initialize_state(dof)

print(f"Computed Entropy: {engine.compute_entropy()}")
```

## Testing

Run tests with `pytest`:
```bash
pip install pytest
pytest tests/
```
