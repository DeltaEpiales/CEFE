import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
import matplotlib.pyplot as plt
import cefe_py as ce
import time
from osacra_gravaser import HolographicEmergentGravity, OsacraGravaserHoTT, BackreactionEngine

def normalize(data):
    """Normalize a list/array of data to 0-1 for direct curve shape comparison."""
    min_val = np.min(data)
    max_val = np.max(data)
    if max_val == min_val:
        return data
    return (data - min_val) / (max_val - min_val)

def sim_1_entropic_rest_mass():
    print("Running Simulation 1: Entropic Rest Mass vs. Horizon Accumulation...")
    
    # 1. CEFE Physics Engine Calculations (Raw Lattice Entropy)
    spacing = 0.5
    max_radius = 4.0
    grid = ce.CausalDiamondGrid(radius=max_radius, spacing=spacing)
    engine = ce.EntropicFieldEngine(grid)
    engine.set_mass(1.0)
    engine.evolve_to_slice(0.0) # Ground state
    
    # Test varying subregion radii
    radii = np.arange(0.5, 3.5, 0.25)
    cefe_entropies = []
    layers_N = []
    
    for r in radii:
        S = engine.compute_entanglement_entropy(r)
        cefe_entropies.append(S)
        # N = alpha * (r / spacing)^2  to properly align with Area Law (D_H = 3/2)
        layers_N.append((r / spacing) ** 2) 
        
    cefe_entropies = np.array(cefe_entropies)
    layers_N = np.array(layers_N)
    
    # 2. HEG Theoretical Calculations (Entropic Rest Mass)
    # Using default UV_cutoff_eV=1e9 (which yields the correctly scaled mass)
    bre = BackreactionEngine() 
    heg_rest_masses = np.array([bre.entropic_rest_mass(N) for N in layers_N])
    
    # 3. Plotting
    plt.figure(figsize=(10, 6))
    
    # Since CEFE raw entropy and m0-scaled rest mass operate on different absolute scales,
    # we normalize both to perfectly compare the topological accumulation curves.
    plt.plot(layers_N, normalize(cefe_entropies), 'bo-', label='CEFE Raw Entanglement Entropy', linewidth=2, markersize=8)
    plt.plot(layers_N, normalize(heg_rest_masses), 'r--', label='HEG Entropic Rest Mass $M(N)$', linewidth=2)
    
    plt.xlabel('Discrete Polarization Layers ($N \propto r^2$)')
    plt.ylabel('Normalized Magnitude')
    plt.title('Simulation 1: Entropic Rest Mass as Emergent Horizon Entanglement')
    plt.legend()
    plt.grid(True)
    plt.savefig('sim1_rest_mass.png')
    print("  -> Saved 'sim1_rest_mass.png'")


def sim_2_thermodynamic_time_backreaction():
    print("Running Simulation 2: Thermodynamic Time and Energy Backreaction...")
    
    # 1. CEFE Physics Engine Calculations (Dynamic Energy Evolution)
    ads_metric = ce.AntiDeSitterMetric(L=10.0)
    grid = ce.CausalDiamondGrid(radius=3.0, spacing=0.5, metric=ads_metric)
    engine = ce.EntropicFieldEngine(grid)
    engine.set_mass(1.0)
    engine.set_interaction_coupling(0.5) # Strong interaction to drive backreaction
    
    # Initialize in thermal state
    engine.initialize_field_state(target_t=0.0, temperature=5.0)
    
    steps = 400
    dt = 0.05
    cefe_energy = []
    time_points = []
    
    current_time = 0.0
    for i in range(steps):
        cefe_energy.append(engine.get_field_energy())
        time_points.append(current_time)
        engine.step_forward(dt)
        current_time += dt
        
    cefe_energy = np.array(cefe_energy)
    
    # 2. HEG Theoretical Calculations (Entropic Time)
    hott = OsacraGravaserHoTT()
    simulated_N_layers = np.linspace(1, 100, steps)
    heg_entropic_time = np.array([hott.entropic_time_accumulation(N) for N in simulated_N_layers])
    
    # 3. Plotting
    fig, ax1 = plt.subplots(figsize=(10, 6))

    color = 'tab:red'
    ax1.set_xlabel('Simulated Time Evolution ($t$)')
    ax1.set_ylabel('CEFE Field Energy Backreaction', color=color)
    ax1.plot(time_points, cefe_energy, color=color, label='CEFE Field Energy', linewidth=2)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()
    color = 'tab:blue'
    ax2.set_ylabel('HEG Accumulated Entropic Time', color=color)
    ax2.plot(time_points, heg_entropic_time, color=color, linestyle='--', label='HEG Entropic Time $t(N)$', linewidth=2)
    ax2.tick_params(axis='y', labelcolor=color)

    plt.title('Simulation 2: Thermodynamic Energy Backreaction across Fractaled AdS Time')
    fig.tight_layout()
    plt.grid(True, alpha=0.3)
    plt.savefig('sim2_energy_time.png')
    print("  -> Saved 'sim2_energy_time.png'")

if __name__ == "__main__":
    t0 = time.time()
    print("=== OSACRA-GRAVASER HEG CEFE SIMULATIONS ===")
    sim_1_entropic_rest_mass()
    print()
    sim_2_thermodynamic_time_backreaction()
    print(f"\nCompleted visual physical simulations in {time.time() - t0:.2f}s.")
