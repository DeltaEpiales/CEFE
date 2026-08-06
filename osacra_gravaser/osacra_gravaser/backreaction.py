import numpy as np
import scipy.linalg
from .constants import M_SCALE_EV, eV_to_Joules

class RenormalizationGroup:
    def __init__(self, M_IR=M_SCALE_EV):
        """
        Handles the Piecewise Renormalization Group (RG) running coupling.
        Integrates across the electron, muon, and QCD confinement thresholds.
        """
        self.M_IR = M_IR
        
        # Define the mass thresholds (in eV) where new degrees of freedom "wake up"
        # and the corresponding beta function coefficients for those energy bands.
        self.thresholds = [
            (0, 511e3, 0.012),       # Below Electron mass (Dark sector + photons)
            (511e3, 105e6, 0.015),   # Electron threshold to Muon threshold
            (105e6, 200e6, 0.018),   # Muon threshold to QCD confinement
            (200e6, float('inf'), 0.022) # Post-QCD (Hadronic vacuum polarization active)
        ]

    def run_mass_to_scale(self, mu_target_eV):
        """
        Piecewise integration of the RG running mass.
        """
        if mu_target_eV <= self.M_IR:
            return self.M_IR
            
        current_mass = self.M_IR
        current_mu = self.M_IR
        
        for lower, upper, beta in self.thresholds:
            if current_mu >= upper:
                continue
                
            integration_upper = min(mu_target_eV, upper)
            
            # Integrate this piece: m(mu_up) = m(mu_low) * (1 + beta * ln(mu_up / mu_low))
            current_mass = current_mass * (1 + beta * np.log(integration_upper / current_mu))
            current_mu = integration_upper
            
            if current_mu >= mu_target_eV:
                break
                
        return current_mass


class BackreactionEngine:
    def __init__(self, UV_cutoff_eV=1e9):
        """
        Initialize the Non-linear Backreaction Framework with Piecewise RG flow.
        """
        self.rg_flow = RenormalizationGroup()
        self.m0_eV = self.rg_flow.run_mass_to_scale(UV_cutoff_eV)
        self.m0_kg = (self.m0_eV * eV_to_Joules) / (299792458.0**2) 

    def entropic_rest_mass(self, N_layers):
        return self.m0_kg * N_layers * (2 * np.log(2))

    def knot_tensor_fractional_diff(self, R2_tensor, g_tensor, u_vec, delta_x_frac):
        D_H = 1.5 
        alpha = D_H / 2.0 # Yields 0.75 via Riemann-Liouville
        
        trace_R2 = np.trace(np.dot(g_tensor, R2_tensor))
        trace_free_pol = R2_tensor - 0.25 * g_tensor * trace_R2
        
        fractional_op_scalar = np.sum(u_vec) / (delta_x_frac ** alpha) 
        K_tensor = fractional_op_scalar * trace_free_pol
        return K_tensor

    def mera_tensor_network_simulate_layer(self, layer_index, initial_density_matrix):
        """
        Simulates the peeling layers using a true Completely Positive Trace-Preserving 
        (CPTP) Kraus operator formulation for generalized phase damping decoherence.
        """
        # Determine the decay probability p based on layer_index
        p = min(0.01 * layer_index, 0.99)
        
        dim = initial_density_matrix.shape[0]
        I = np.eye(dim, dtype=complex)
        
        # We will use a generalized phase damping channel using standard Pauli Z logic 
        # extended to higher dimensions via a diagonal alternating matrix.
        Z = np.diag([1 if i % 2 == 0 else -1 for i in range(dim)])
        
        # Construct Kraus Operators
        K0 = np.sqrt(1 - p) * I
        K1 = np.sqrt(p) * Z
        
        # Verify Completeness: K0^dagger K0 + K1^dagger K1 = I (Mathematically guaranteed)
        
        # Apply the CPTP channel: rho' = K0 rho K0^dagger + K1 rho K1^dagger
        rho_prime = K0 @ initial_density_matrix @ K0.conj().T + K1 @ initial_density_matrix @ K1.conj().T
        
        # The trace is strictly preserved by the algebra of Kraus operators
        mock_winding_number = int(np.round(np.real(np.trace(rho_prime)) * 10)) % 3
        return rho_prime, mock_winding_number
