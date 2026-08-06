import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import cefe_py as ce
import time

def create_cefe_gif():
    print("Initializing Causal Entropic Field Engine (CEFE) for GIF generation...")
    
    # 1. Setup the Engine Geometry
    # We use a relatively small grid so it plots quickly
    radius = 3.0
    spacing = 0.5
    ads_metric = ce.AntiDeSitterMetric(L=10.0)
    grid = ce.CausalDiamondGrid(radius=radius, spacing=spacing, metric=ads_metric)
    
    engine = ce.EntropicFieldEngine(grid)
    engine.set_mass(1.0)
    engine.set_interaction_coupling(0.5)
    
    # Initialize a heated thermal state
    engine.initialize_field_state(target_t=0.0, temperature=5.0)

    # 2. Extract Grid Geometry for Visualization
    points = grid.get_points()
    x_coords = [p.x for p in points]
    y_coords = [p.y for p in points]
    is_boundary = [p.is_boundary for p in points]
    
    # We will simulate and record the energy evolution
    frames = 60
    dt = 0.05
    energy_history = []
    time_history = []
    
    print("Simulating engine frames...")
    current_time = 0.0
    for i in range(frames):
        energy_history.append(engine.get_field_energy())
        time_history.append(current_time)
        engine.step_forward(dt)
        current_time += dt

    # 3. Setup Matplotlib Figure for Animation
    fig = plt.figure(figsize=(12, 5))
    fig.canvas.manager.set_window_title("CEFE Thermodynamic Simulation")
    
    # Subplot 1: The Causal Diamond Grid Space
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.set_title("Fractal Causal Diamond (AdS Geometry)")
    ax1.set_xlim(-radius - 1, radius + 1)
    ax1.set_ylim(-radius - 1, radius + 1)
    ax1.set_aspect('equal')
    ax1.set_xlabel("X Space")
    ax1.set_ylabel("Y Space")
    
    # Separate boundary and interior points
    x_int = [x for x, b in zip(x_coords, is_boundary) if not b]
    y_int = [y for y, b in zip(y_coords, is_boundary) if not b]
    x_bound = [x for x, b in zip(x_coords, is_boundary) if b]
    y_bound = [y for y, b in zip(y_coords, is_boundary) if b]
    
    ax1.scatter(x_int, y_int, c='blue', alpha=0.3, label="Interior Bulk")
    ax1.scatter(x_bound, y_bound, c='red', alpha=0.8, label="Entangling Boundary")
    
    # An expanding ring to simulate the "peeling layers" or thermodynamic wave
    wave_ring, = ax1.plot([], [], 'g-', lw=2, alpha=0.7)
    
    # Subplot 2: Energy Backreaction Live Chart
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.set_title("Thermodynamic Field Energy Backreaction")
    ax2.set_xlim(0, frames * dt)
    ax2.set_ylim(min(energy_history) * 0.9, max(energy_history) * 1.1)
    ax2.set_xlabel("Emergent Entropic Time ($t$)")
    ax2.set_ylabel("Total Field Energy")
    
    line_energy, = ax2.plot([], [], 'r-', lw=2)
    time_text = ax2.text(0.05, 0.9, '', transform=ax2.transAxes)

    def init():
        wave_ring.set_data([], [])
        line_energy.set_data([], [])
        time_text.set_text('')
        return wave_ring, line_energy, time_text

    def update(frame):
        # 1. Update the wave ring (expanding through the lattice)
        # Modulo radius to make it repeat
        current_r = (frame * 0.1) % radius
        theta = np.linspace(0, 2*np.pi, 100)
        wave_ring.set_data(current_r * np.cos(theta), current_r * np.sin(theta))
        
        # 2. Update the live energy chart
        x_data = time_history[:frame+1]
        y_data = energy_history[:frame+1]
        line_energy.set_data(x_data, y_data)
        
        time_text.set_text(f"Layer step: {frame}\nEnergy: {energy_history[frame]:.2f}")
        return wave_ring, line_energy, time_text

    print(f"Generating GIF animation ({frames} frames)...")
    ani = animation.FuncAnimation(fig, update, frames=frames, init_func=init, blit=True, interval=100)
    
    # Save the animation using the pillow writer (standard, no ffmpeg required)
    output_filename = "cefe_thermodynamic_wind.gif"
    ani.save(output_filename, writer='pillow', fps=10)
    
    print(f"Successfully saved GIF to {output_filename}")

if __name__ == "__main__":
    create_cefe_gif()
