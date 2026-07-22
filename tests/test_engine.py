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
