import cefe_py as ce
import time

def run_experiment_1():
    print("=" * 50)
    print("EXPERIMENT 1: The Area Law of Entanglement Entropy")
    print("=" * 50)
    
    # Create a large grid
    grid = ce.CausalDiamondGrid(radius=3.5, spacing=0.4)
    engine = ce.EntropicFieldEngine(grid)
    engine.set_mass(1.0)
    
    print(f"Grid generated with {grid.get_dof()} degrees of freedom.")
    
    # Evolve to vacuum slice
    engine.evolve_to_slice(0.0)
    
    radii = [0.4, 0.8, 1.2, 1.6, 2.0, 2.4]
    
    print(f"{'Radius (r)':<15} | {'Area (~r^2)':<15} | {'Volume (~r^3)':<15} | {'Entropy (S)':<15}")
    print("-" * 65)
    for r in radii:
        S = engine.compute_entanglement_entropy(r)
        area_proxy = r**2
        vol_proxy = r**3
        print(f"{r:<15.1f} | {area_proxy:<15.2f} | {vol_proxy:<15.2f} | {S:<15.4f}")
        
    print("\nIf S scales with r^2, the Area Law is validated.\n")

def run_experiment_2():
    print("=" * 50)
    print("EXPERIMENT 2: Symplectic Energy Conservation")
    print("=" * 50)
    
    grid = ce.CausalDiamondGrid(radius=3.0, spacing=0.5)
    engine = ce.EntropicFieldEngine(grid)
    engine.set_mass(1.0)
    engine.initialize_field_state(0.0)
    
    steps = 5000
    dt = 0.01
    
    print(f"Simulating {steps} steps with dt={dt}...")
    
    energy_history = []
    
    # Run the Leapfrog evolution
    for i in range(steps + 1):
        if i % 1000 == 0:
            energy_history.append((i, engine.get_field_energy()))
        engine.step_forward(dt)
        
    print(f"{'Step':<10} | {'Field Energy':<20}")
    print("-" * 33)
    for step, E in energy_history:
        print(f"{step:<10} | {E:<20.8f}")
        
    print("\nIf energy oscillates but does not continuously diverge, conservation is validated.\n")

def run_experiment_3():
    print("=" * 50)
    print("EXPERIMENT 3: Flat vs. Curved Spacetime Entanglement")
    print("=" * 50)
    
    # Fixed parameters
    radius = 3.5
    spacing = 0.4
    subregion = 1.6
    
    # 1. Minkowski
    print("Setting up Minkowski (Flat) Spacetime...")
    grid_flat = ce.CausalDiamondGrid(radius=radius, spacing=spacing)
    engine_flat = ce.EntropicFieldEngine(grid_flat)
    engine_flat.evolve_to_slice(0.0)
    S_flat = engine_flat.compute_entanglement_entropy(subregion)
    
    # 2. Schwarzschild Rs = 0.5
    print("Setting up Schwarzschild (Rs = 0.5)...")
    metric_s1 = ce.SchwarzschildMetric(0.5)
    grid_s1 = ce.CausalDiamondGrid(radius=radius, spacing=spacing, metric=metric_s1)
    engine_s1 = ce.EntropicFieldEngine(grid_s1)
    engine_s1.evolve_to_slice(0.0)
    S_s1 = engine_s1.compute_entanglement_entropy(subregion)
    
    # 3. Schwarzschild Rs = 1.0
    print("Setting up Schwarzschild (Rs = 1.0)...")
    metric_s2 = ce.SchwarzschildMetric(1.0)
    grid_s2 = ce.CausalDiamondGrid(radius=radius, spacing=spacing, metric=metric_s2)
    engine_s2 = ce.EntropicFieldEngine(grid_s2)
    engine_s2.evolve_to_slice(0.0)
    S_s2 = engine_s2.compute_entanglement_entropy(subregion)
    
    print("\nResults for subregion r =", subregion)
    print(f"{'Spacetime Geometry':<25} | {'Degrees of Freedom':<20} | {'Entropy (S)':<15}")
    print("-" * 65)
    print(f"{'Minkowski (Rs=0)':<25} | {grid_flat.get_dof():<20} | {S_flat:<15.4f}")
    print(f"{'Schwarzschild (Rs=0.5)':<25} | {grid_s1.get_dof():<20} | {S_s1:<15.4f}")
    print(f"{'Schwarzschild (Rs=1.0)':<25} | {grid_s2.get_dof():<20} | {S_s2:<15.4f}")
    
    print("\nIf entropy and DOF change significantly due to horizon dilation, metric coupling is validated.\n")


if __name__ == "__main__":
    t0 = time.time()
    run_experiment_1()
    run_experiment_2()
    run_experiment_3()
    print(f"All experiments completed in {time.time() - t0:.2f} seconds.")
