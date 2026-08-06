import argparse
import sys
import cefe_py as ce
from PyQt6.QtWidgets import QApplication
from render_engine import CefeRenderEngine

def run_cli():
    parser = argparse.ArgumentParser(description="CEFE Quantum Simulator CLI")
    parser.add_argument("--mode", type=str, choices=['qm', 'holographic', 'qft'], default='qm', help="Engine simulation mode")
    parser.add_argument("--potential", type=str, choices=['free', 'double_slit', 'harmonic'], default='free', help="Potential energy surface")
    parser.add_argument("--render", type=str, choices=['2d', '3d', 'none'], default='none', help="Rendering mode")
    parser.add_argument("--steps", type=int, default=100, help="Number of time steps")
    parser.add_argument("--dt", type=float, default=0.01, help="Time step size")
    
    args = parser.parse_args()
    print(f"Starting CEFE in {args.mode.upper()} mode...")
    
    if args.mode == 'qm':
        grid = ce.CausalDiamondGrid(radius=2.0, spacing=0.1)
        engine = ce.QuantumWaveEngine(grid)
        engine.set_mass(1.0)
        engine.set_hbar(1.0)
        engine.initialize_state(target_t=0.0)
        engine.set_gaussian_packet(x0=-1.0, y0=0.0, sigma=0.2, px=5.0, py=0.0)
        
        if args.potential == 'double_slit':
            engine.set_double_slit_potential(slit_width=0.4, slit_separation=1.0, barrier_thickness=0.2, barrier_x=0.0)
            
        engine.build_hamiltonian(target_t=0.0)
        engine.prepare_crank_nicolson(args.dt)
        
        if args.render != 'none':
            app = QApplication(sys.argv)
            renderer = CefeRenderEngine(engine, mode='qm')
            if args.render == '3d':
                renderer.render_3d_wave_animation(steps=args.steps, dt=args.dt)
            elif args.render == '2d':
                renderer.render_2d_heatmap()
            sys.exit(app.exec())
        else:
            print(f"Running headless for {args.steps} steps...")
            for i in range(args.steps):
                engine.step_forward()
            print("Simulation complete. Final total probability:", engine.get_total_probability())
            
    elif args.mode == 'holographic':
        print("Initializing Holographic Emergent Gravity Engine...")
        radius = 3.0
        spacing = 0.5
        ads_metric = ce.AntiDeSitterMetric(L=10.0)
        grid = ce.CausalDiamondGrid(radius=radius, spacing=spacing, metric=ads_metric)
        
        engine = ce.EntropicFieldEngine(grid)
        engine.set_mass(1.0)
        engine.set_interaction_coupling(0.5)
        engine.initialize_field_state(target_t=0.0, temperature=5.0)
        
        if args.render != 'none':
            app = QApplication(sys.argv)
            renderer = CefeRenderEngine(engine, mode='heg')
            if args.render == '3d':
                renderer.render_3d_wave_animation(steps=args.steps, dt=args.dt)
            elif args.render == '2d':
                renderer.render_2d_heatmap()
            sys.exit(app.exec())
        else:
            print(f"Running headless for {args.steps} steps...")
            for i in range(args.steps):
                engine.step_forward(args.dt)
            print("Simulation complete. Final total energy:", engine.get_field_energy())
            
    elif args.mode == 'qft':
        print("Initializing QFT Interacting C++ Engine...")
        grid = ce.CausalDiamondGrid(radius=2.0, spacing=0.1)
        engine = ce.LatticeQFTEngine(grid)
        engine.set_mass(1.0)
        engine.set_coupling(10.0)
        engine.initialize_state()
        engine.set_opposing_packets(1.0, 0.2, 5.0, 5.0)
        engine.build_laplacian()
        
        if args.render != 'none':
            app = QApplication(sys.argv)
            renderer = CefeRenderEngine(engine, mode='qft')
            if args.render == '3d':
                renderer.render_3d_wave_animation(steps=args.steps, dt=args.dt)
            elif args.render == '2d':
                renderer.render_2d_heatmap()
            sys.exit(app.exec())
        else:
            print(f"Running headless for {args.steps} steps...")
            for i in range(args.steps):
                engine.step_forward(args.dt)
            print("Simulation complete.")

if __name__ == "__main__":
    run_cli()
