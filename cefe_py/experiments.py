"""
cefe_py.experiments -- High-level experiment runner functions
=============================================================

Provides convenient Python functions for running the core physics
experiments programmatically. Each function returns structured results
suitable for analysis in notebooks or scripts.

Example::

    from cefe_py.experiments import run_area_law_test
    result = run_area_law_test()
    print(f"Area law exponent: {result['alpha']:.3f}")
"""

import numpy as np
import cefe_py as ce


def _power_law_fit(x, y):
    """Fit y = a * x^alpha via log-log linear regression."""
    lx = np.log(np.array(x, dtype=float))
    ly = np.log(np.array(y, dtype=float))
    n = len(lx)
    sx, sy = lx.sum(), ly.sum()
    sxx, sxy = (lx * lx).sum(), (lx * ly).sum()
    alpha = (n * sxy - sx * sy) / (n * sxx - sx * sx)
    log_a = (sy - alpha * sx) / n
    ss_res = ((ly - log_a - alpha * lx) ** 2).sum()
    ss_tot = ((ly - ly.mean()) ** 2).sum()
    return np.exp(log_a), alpha, 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0


def _linear_fit(x, y):
    """Fit y = a*x + b via OLS."""
    x, y = np.array(x, dtype=float), np.array(y, dtype=float)
    n = len(x)
    sx, sy = x.sum(), y.sum()
    sxx, sxy = (x * x).sum(), (x * y).sum()
    a = (n * sxy - sx * sy) / (n * sxx - sx * sx)
    b = (sy - a * sx) / n
    ss_res = ((y - a * x - b) ** 2).sum()
    ss_tot = ((y - y.mean()) ** 2).sum()
    return a, b, 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0


def run_area_law_test(radius=5.0, spacing=0.8, mass=1.0,
                      radii=None, metric=None):
    """Test entanglement area law: S(r) ~ r^alpha.

    Args:
        radius: Causal diamond radius
        spacing: Lattice spacing
        mass: Scalar field mass
        radii: List of subregion radii to test (default: 8 log-spaced values)
        metric: Optional spacetime metric (default: Minkowski)

    Returns:
        dict with keys: radii, entropies, alpha, a, r_squared
    """
    if radii is None:
        radii = [1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0]

    if metric is not None:
        grid = ce.CausalDiamondGrid(radius=radius, spacing=spacing, metric=metric)
    else:
        grid = ce.CausalDiamondGrid(radius=radius, spacing=spacing)

    engine = ce.EntropicFieldEngine(grid)
    engine.set_mass(mass)
    engine.evolve_to_slice(0.0)

    entropies = [engine.compute_entanglement_entropy(r) for r in radii]
    a, alpha, r_sq = _power_law_fit(radii, entropies)

    return {
        "radii": radii,
        "entropies": entropies,
        "alpha": alpha,
        "a": a,
        "r_squared": r_sq,
        "dof": engine.get_state_size(),
    }


def run_continuum_limit_test(radius=4.0, subregion_r=1.5, mass=1.0,
                             spacings=None):
    """Test UV scaling: S ~ 1/eps^2 at fixed geometry.

    Args:
        radius: Causal diamond radius
        subregion_r: Fixed subregion radius
        mass: Scalar field mass
        spacings: List of lattice spacings (default: [1.2, 1.0, 0.8, 0.7, 0.6])

    Returns:
        dict with keys: spacings, inv_eps2, entropies, slope, intercept, r_squared
    """
    if spacings is None:
        spacings = [1.2, 1.0, 0.8, 0.7, 0.6]

    entropies = []
    inv_eps2 = []
    for eps in spacings:
        grid = ce.CausalDiamondGrid(radius=radius, spacing=eps)
        engine = ce.EntropicFieldEngine(grid)
        engine.set_mass(mass)
        engine.evolve_to_slice(0.0)
        S = engine.compute_entanglement_entropy(subregion_r)
        entropies.append(S)
        inv_eps2.append(1.0 / (eps * eps))

    a, b, r_sq = _linear_fit(inv_eps2, entropies)

    return {
        "spacings": spacings,
        "inv_eps2": inv_eps2,
        "entropies": entropies,
        "slope": a,
        "intercept": b,
        "r_squared": r_sq,
    }


def run_mass_dependence_test(radius=4.0, spacing=0.8, subregion_r=1.5,
                             masses=None):
    """Test mass gap suppression of entanglement.

    Args:
        radius: Causal diamond radius
        spacing: Lattice spacing
        subregion_r: Subregion radius
        masses: List of masses to test

    Returns:
        dict with keys: masses, entropies, monotonically_decreasing
    """
    if masses is None:
        masses = [0.1, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0]

    entropies = []
    for m in masses:
        grid = ce.CausalDiamondGrid(radius=radius, spacing=spacing)
        engine = ce.EntropicFieldEngine(grid)
        engine.set_mass(m)
        engine.evolve_to_slice(0.0)
        S = engine.compute_entanglement_entropy(subregion_r)
        entropies.append(S)

    monotonic = all(entropies[i] >= entropies[i + 1] for i in range(len(entropies) - 1))

    return {
        "masses": masses,
        "entropies": entropies,
        "monotonically_decreasing": monotonic,
    }


def run_thermal_test(radius=3.0, spacing=0.8, mass=1.0, subregion_r=1.0,
                     temperatures=None):
    """Test thermal state entropy scaling: S(T) >= S(0).

    Args:
        radius: Causal diamond radius
        spacing: Lattice spacing
        mass: Scalar field mass
        subregion_r: Subregion radius
        temperatures: List of temperatures

    Returns:
        dict with keys: temperatures, entropies, monotonically_increasing, above_vacuum
    """
    if temperatures is None:
        temperatures = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0]

    entropies = []
    for T in temperatures:
        grid = ce.CausalDiamondGrid(radius=radius, spacing=spacing)
        engine = ce.EntropicFieldEngine(grid)
        engine.set_mass(mass)
        engine.evolve_to_slice(0.0, temperature=T)
        S = engine.compute_entanglement_entropy(subregion_r)
        entropies.append(S)

    monotonic = all(entropies[i] <= entropies[i + 1] for i in range(len(entropies) - 1))
    above_vacuum = all(S >= entropies[0] - 1e-10 for S in entropies)

    return {
        "temperatures": temperatures,
        "entropies": entropies,
        "monotonically_increasing": monotonic,
        "above_vacuum": above_vacuum,
    }


def run_energy_conservation_test(radius=3.0, spacing=0.5, mass=1.0,
                                 coupling=0.5, steps=5000, dt=0.01):
    """Test symplectic energy conservation over long evolution.

    Args:
        radius: Causal diamond radius
        spacing: Lattice spacing
        mass: Scalar field mass
        coupling: Interaction coupling lambda
        steps: Number of time steps
        dt: Time step size

    Returns:
        dict with keys: times, energies, E0, max_drift, final_drift
    """
    ads = ce.AntiDeSitterMetric(L=10.0)
    grid = ce.CausalDiamondGrid(radius=radius, spacing=spacing, metric=ads)
    engine = ce.EntropicFieldEngine(grid)
    engine.set_mass(mass)
    engine.set_interaction_coupling(coupling)
    engine.initialize_field_state(target_t=0.0, temperature=5.0)

    E0 = engine.get_field_energy()
    times = [0.0]
    energies = [E0]

    for i in range(1, steps + 1):
        engine.step_forward(dt)
        if i % (steps // 20) == 0 or i == steps:
            E = engine.get_field_energy()
            times.append(i * dt)
            energies.append(E)

    max_drift = max(abs(E - E0) / abs(E0) for E in energies) if abs(E0) > 0 else 0
    final_drift = abs(energies[-1] - E0) / abs(E0) if abs(E0) > 0 else 0

    return {
        "times": times,
        "energies": energies,
        "E0": E0,
        "max_drift": max_drift,
        "final_drift": final_drift,
        "steps": steps,
        "dt": dt,
    }
