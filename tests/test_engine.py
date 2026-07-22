import sys
import os
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../cefe_py')))
import cefe_py as ce

def test_causal_diamond():
    grid = ce.CausalDiamondGrid(radius=3.0, spacing=1.0)
    assert grid.get_radius() == 3.0
    assert grid.get_dof() > 0

def test_entropic_field_engine_ground_state():
    grid = ce.CausalDiamondGrid(radius=3.0, spacing=1.0)
    engine = ce.EntropicFieldEngine(grid)
    engine.set_mass(1.0)
    
    # Generate covariance matrix at t=0
    engine.evolve_to_slice(0.0)
    
    # State size should match the number of points at t=0
    dof = engine.get_state_size()
    assert dof > 0
    
    # Test entropy
    entropy = engine.compute_entanglement_entropy(1.5)
    
    # Entanglement entropy should be positive
    assert entropy >= 0.0

def test_area_law():
    grid = ce.CausalDiamondGrid(radius=5.0, spacing=1.0)
    engine = ce.EntropicFieldEngine(grid)
    engine.set_mass(1.0)
    engine.evolve_to_slice(0.0)
    
    s1 = engine.compute_entanglement_entropy(1.0)
    s2 = engine.compute_entanglement_entropy(2.0)
    
    # S(2.0) should be greater than S(1.0) due to larger boundary area
    assert s2 >= s1

def test_bekenstein_bound_compliance():
    """Verify that computed entanglement entropy respects the Bekenstein Area Bound"""
    grid = ce.CausalDiamondGrid(radius=4.0, spacing=1.0)
    engine = ce.EntropicFieldEngine(grid)
    engine.set_mass(1.0)
    engine.evolve_to_slice(0.0)
    
    r = 2.0
    entropy = engine.compute_entanglement_entropy(r)
    bekenstein_bound = grid.get_area_bound(r)
    
    # Ground state entanglement entropy must be strictly less than the maximum possible thermodynamic entropy
    assert entropy <= bekenstein_bound

def test_thermal_entropy_scaling():
    """Verify that finite temperature (T > 0) mixed states have higher entropy than the T=0 pure vacuum"""
    grid = ce.CausalDiamondGrid(radius=3.0, spacing=1.0)
    engine = ce.EntropicFieldEngine(grid)
    engine.set_mass(1.0)
    
    engine.evolve_to_slice(target_t=0.0, temperature=0.0)
    s_vacuum = engine.compute_entanglement_entropy(1.0)
    
    engine.evolve_to_slice(target_t=0.0, temperature=10.0)
    s_thermal = engine.compute_entanglement_entropy(1.0)
    
    # Thermal mixed state must possess greater entanglement entropy due to thermal fluctuations
    assert s_thermal > s_vacuum

def test_ads_volume_scaling():
    """Verify that Anti-de Sitter space packs more degrees of freedom into the same coordinate radius than Minkowski space due to negative curvature"""
    minkowski_grid = ce.CausalDiamondGrid(radius=5.0, spacing=1.0, metric=ce.MinkowskiMetric())
    ads_grid = ce.CausalDiamondGrid(radius=5.0, spacing=1.0, metric=ce.AntiDeSitterMetric(L=1.0))
    
    dof_flat = minkowski_grid.get_dof()
    dof_curved = ads_grid.get_dof()
    
    # Hyperbolic geometry expands spatial volume near the boundary
    assert dof_curved > dof_flat
