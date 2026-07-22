# Installed via pip, no sys.path hack needed.

import cefe_py as ce

def test_causal_diamond():
    print("Testing CausalDiamondGrid with Minkowski...")
    grid = ce.CausalDiamondGrid(radius=2.2, spacing=1.5)
    
    print(f"Radius: {grid.get_radius()}")
    print(f"Spacing: {grid.get_spacing()}")
    print(f"Degrees of Freedom: {grid.get_dof()}")
    
    bound = grid.get_area_bound(1.0)
    print(f"Bekenstein Bound for r=1.0: {bound}")

def test_entropic_engine():
    print("Testing EntropicFieldEngine (Minkowski Ground State)...")
    grid = ce.CausalDiamondGrid(radius=2.2, spacing=1.5)
    
    engine = ce.EntropicFieldEngine(grid)
    engine.set_mass(1.0)
    
    # Static entanglement entropy calculation
    engine.evolve_to_slice(0.0)
    print(f"3D Slice State Size: {engine.get_state_size()}")
    
    entropy = engine.compute_entanglement_entropy(1.0)
    print(f"Computed Entropy: {entropy}")

def test_full_framework():
    print("\n--- Testing Full Framework Features ---")
    print("1. Schwarzschild Metric (Rs = 0.5)")
    schwarzschild = ce.SchwarzschildMetric(0.5)
    curved_grid = ce.CausalDiamondGrid(radius=2.2, spacing=1.0, metric=schwarzschild)
    print(f"Curved Degrees of Freedom: {curved_grid.get_dof()}")

    print("2. Dynamic Time Evolution")
    engine = ce.EntropicFieldEngine(curved_grid)
    engine.set_mass(1.0)
    
    # Initialize the field on slice t=0
    engine.initialize_field_state(0.0)
    energy_initial = engine.get_field_energy()
    print(f"Initial Field Energy: {energy_initial}")
    
    # Step forward in time (dt = 0.01) for 100 steps
    for _ in range(100):
        engine.step_forward(0.01)
        
    energy_final = engine.get_field_energy()
    print(f"Final Field Energy after t=1.0: {energy_final}")
    
if __name__ == "__main__":
    test_causal_diamond()
    test_entropic_engine()
    test_full_framework()
    print("All tests passed!")
