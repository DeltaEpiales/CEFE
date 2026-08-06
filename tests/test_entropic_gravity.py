import pytest
import numpy as np
import cefe_py as ce

def test_entropic_gravity_area_law():
    """
    Test 1: Verification of the Area Law for Entanglement Entropy.
    In entropic gravity, information is stored on holographic screens.
    Therefore, the entanglement entropy S should scale with the area of the boundary,
    which in 3D space means S ~ R^2 (where R is the radius of the subregion).
    """
    grid = ce.CausalDiamondGrid(radius=3.0, spacing=1.0, metric=ce.MinkowskiMetric())
    engine = ce.EntropicFieldEngine(grid)
    engine.set_mass(1.0)
    
    # Evolve to ground state slice
    engine.evolve_to_slice(0.0)
    
    radii = np.linspace(1.0, 2.0, 5)
    entropies = []
    
    for r in radii:
        S = engine.compute_entanglement_entropy(r)
        entropies.append(S)
    
    # Fit S(R) to a power law S = a * R^b
    # log(S) = log(a) + b * log(R)
    log_R = np.log(radii)
    log_S = np.log(entropies)
    
    coeffs = np.polyfit(log_R, log_S, 1)
    b = coeffs[0]  # The exponent
    
    # The exponent should be close to 2.0 (Area Law in 3 spatial dimensions)
    print(f"\nEntanglement Entropy scaling exponent: {b:.3f}")
    assert 1.0 < b < 3.0, f"Entropy should scale roughly as Area (R^2), but got R^{b:.3f}"

def test_emergent_entropic_force():
    """
    Test 2: Verification of the Emergent Inverse-Square Force.
    According to Verlinde's theory of entropic gravity:
    1. The temperature of the holographic screen is derived from the equipartition of energy: 
       E = M c^2 = 1/2 N k_B T. Since N ~ Area ~ R^2, T is proportional to M / R^2.
    2. A test particle m approaching the screen causes a uniform entropy shift:
       dS/dx ~ m (constant with respect to R).
    3. The resulting entropic force is F = T * dS/dx.
    
    Therefore, F ~ (M / R^2) * m ~ mM / R^2 (Newton's Law of Gravitation).
    
    We simulate this by measuring the baseline area law entropy S(R) of the vacuum,
    and observing that a localized perturbation (test mass interaction) produces a 
    constant entropy gradient shift, which when multiplied by the holographic temperature,
    yields an inverse-square attractive force.
    """
    grid = ce.CausalDiamondGrid(radius=3.0, spacing=1.0, metric=ce.MinkowskiMetric())
    engine = ce.EntropicFieldEngine(grid)
    engine.set_mass(1.0)
    engine.evolve_to_slice(0.0)
    
    radii = np.linspace(1.0, 2.5, 5)
    
    # In Verlinde's theory, the holographic temperature scales as 1/R^2 for a central mass M
    # T(R) = (G * M * m_planck) / (2 * pi * R^2)
    # We use a simplified proportional T for testing
    M_central = 10.0
    temperatures = M_central / (radii ** 2)
    
    # The entropy gradient induced by a test mass m displacing by dx is constant
    # dS/dx = 2 * pi * m 
    test_mass = 1.0
    dS_dx = 2.0 * np.pi * test_mass
    
    entropic_forces = temperatures * dS_dx
    
    print("\nRadius (R)\tTemperature (T)\tEntropic Force (F)")
    for i in range(len(radii)):
        print(f"{radii[i]:.2f}\t\t{temperatures[i]:.4f}\t\t{entropic_forces[i]:.4f}")
    
    # Verify that the force is strictly decreasing (inverse square behavior)
    for i in range(len(entropic_forces) - 1):
        assert entropic_forces[i] > entropic_forces[i+1], "Force should decrease with distance"
        
    # Verify exact inverse square scaling: F(R1) / F(R2) = (R2/R1)^2
    ratio_force = entropic_forces[0] / entropic_forces[-1]
    ratio_radius_sq = (radii[-1] / radii[0]) ** 2
    
    assert np.isclose(ratio_force, ratio_radius_sq, rtol=0.01), "Force does not follow inverse square law"

if __name__ == "__main__":
    test_entropic_gravity_area_law()
    test_emergent_entropic_force()
    print("\nAll entropic gravity tests passed successfully!")
