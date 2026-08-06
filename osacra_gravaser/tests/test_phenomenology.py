import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
from osacra_gravaser import HolographicEmergentGravity, OsacraGravaserHoTT, BackreactionEngine

def run_tests():
    print("=========================================================")
    print(" Holographic Emergent Gravity (HEG) Framework Tests")
    print("=========================================================\n")
    
    # 1. Thermodynamics & Emergent Gravity Tests
    print("--- 1. Thermodynamics & Spacetime Interface (Dynamically Derived from CEFE) ---")
    
    import cefe_py as ce
    
    # Initialize an underlying CEFE grid and engine state to compute topological bounds
    grid = ce.CausalDiamondGrid(radius=4.0, spacing=0.5)
    engine = ce.EntropicFieldEngine(grid)
    engine.set_mass(1.0)
    engine.set_interaction_coupling(0.5)
    engine.initialize_field_state(target_t=0.0, temperature=5.0) # Thermal interacting state
    
    heg = HolographicEmergentGravity(engine=engine)
    
    delta_omega = heg.muon_anomalous_precession_shift()
    print(f"Muon Anomalous Precession Shift (delta omega_HEG): {delta_omega:.5e} [Dynamically Computed]")
    
    t_CMB_to_now = 13.8e9 * 365 * 24 * 3600
    delta_chi = heg.cosmic_birefringence_angle(t_CMB_to_now)
    print(f"Cosmic Birefringence Angle over 13.8B years: {delta_chi:.5e} radians [Dynamically Computed]")
    
    k_val = 1e6
    w_plus, w_minus, shift = heg.photonic_dispersion(k_val)
    print(f"Photonic Dispersion Exact Analytical Splitting: {2*shift:.5e} [Dynamically Computed]\n")

    # 2. Osacra-Gravaser HoTT Formalization Tests
    print("--- 2. Osacra-Gravaser HoTT Formalization ---")
    hott = OsacraGravaserHoTT()
    
    vec_S = hott.entropic_time_accumulation(1000)
    print(f"Accumulated Vectorized Entropy for N=1000 layers: {vec_S:.5f}")
    
    winding_numbers = [1, -2, 1, 3] # Mock windings
    J_max, allowed_spins = hott.topological_matter_spin(winding_numbers)
    print(f"Topological Matter Spin (J_max) for given winding numbers: {J_max:.5e} J s")
    print(f"Allowed SU(2) Spin Channels: {allowed_spins}\n")

    # 3. Non-Linear Backreaction & Fractal Horizons Tests
    print("--- 3. Backreaction, Fractal Horizons & MERA Networks ---")
    bre = BackreactionEngine(UV_cutoff_eV=1e9) 
    
    print(f"RG Flow: Ran IR Cosmological Mass ({bre.rg_flow.M_IR} eV) -> UV Hadronic Base Mass: {bre.m0_eV:.5e} eV")
    M_N = bre.entropic_rest_mass(1000)
    print(f"Entropic Rest Mass Accumulation (M_N) after 1000 layers: {M_N:.5e} kg")
    
    # Knot Tensor Mock
    R2_mock = np.array([[1.0, 0.1, 0.0, 0.0],
                        [0.1, -1.0, 0.2, 0.0],
                        [0.0, 0.2, 1.0, 0.0],
                        [0.0, 0.0, 0.0, -1.0]])
    g_mock = np.eye(4)
    u_mock = np.array([1.0, 0.0, 0.0, 0.0])
    
    K_tensor = bre.knot_tensor_fractional_diff(R2_mock, g_mock, u_mock, delta_x_frac=1e-35)
    print(f"Calculated Topological Knot Tensor K^(N)_mu_nu (alpha=0.75):\n{K_tensor}\n")
    
    # CPTP Density Matrix Test
    rho_mock = np.eye(4) / 4.0
    evolved_rho, winding = bre.mera_tensor_network_simulate_layer(10, rho_mock)
    print(f"CPTP Trace Conservation after 10 layers: Tr(rho) = {np.trace(evolved_rho):.5f}")

if __name__ == "__main__":
    run_tests()
