import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
import cefe_py as ce

def analyze_frames():
    print("--- CEFE Engine Frame-by-Frame Analysis ---")
    
    radius = 3.0
    spacing = 0.5
    ads_metric = ce.AntiDeSitterMetric(L=10.0)
    grid = ce.CausalDiamondGrid(radius=radius, spacing=spacing, metric=ads_metric)
    
    engine = ce.EntropicFieldEngine(grid)
    engine.set_mass(1.0)
    engine.set_interaction_coupling(0.5)
    
    engine.initialize_field_state(target_t=0.0, temperature=5.0)

    frames = 60
    dt = 0.05
    current_time = 0.0
    
    print(f"{'Frame':<8} | {'Time (t)':<10} | {'Field Energy (E)':<18} | {'dE/dt':<15} | {'Entropy (S)':<15}")
    print("-" * 75)
    
    prev_energy = engine.get_field_energy()
    
    for i in range(frames):
        energy = engine.get_field_energy()
        entropy = engine.compute_entanglement_entropy(1.5)
        
        dE_dt = (energy - prev_energy) / dt if i > 0 else 0.0
        
        print(f"{i:<8} | {current_time:<10.2f} | {energy:<18.5f} | {dE_dt:<15.5f} | {entropy:<15.5f}")
        
        prev_energy = energy
        engine.step_forward(dt)
        current_time += dt

if __name__ == "__main__":
    analyze_frames()
