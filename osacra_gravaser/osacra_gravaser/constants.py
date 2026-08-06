# Calibrated Physical Constants for the Osacra-Gravaser Framework

# Standard Constants
c = 299792458.0  # Speed of light in m/s
hbar = 1.054571817e-34  # Reduced Planck constant in J s
m_proton = 1.67262192e-27  # Proton mass in kg
eV_to_Joules = 1.602176634e-19

# Physical Calibration & Phenomenological Bounds
# 
# 1. The Ghost Condensate Scale (M):
# Rather than floating arbitrarily, the symmetry-breaking scale M is explicitly anchored
# to the observed Cosmological Constant / Dark Energy density scale (Lambda_DE ~ 2.4 meV).
# This grounds the thermodynamic wind directly to the macroscopic expansion of the universe.
M_SCALE_EV = 2.4e-3  # 2.4 meV
M_SCALE_JOULES = M_SCALE_EV * eV_to_Joules

# 2. Cosmic Birefringence & Chern-Simons Coupling (g_gamma):
# Planck 2018 bounds the rotation to ~0.3 degrees (0.005 radians) over 13.8B years.
# delta_chi = (g_gamma / 4) * (M / hbar) * delta_t
# To satisfy 0.005 radians, g_gamma must be extremely suppressed, consistent with the 
# lack of intense dark-sector scattering in standard collider physics.
G_GAMMA_DEFAULT = 1.26e-32

# 3. Muon g-2 Tension & Axial Coupling (c_mu):
# The anomalous continuous torque requires a shift of ~2e-13 rad/s.
# delta_omega = 2 * c_mu * (M / hbar)
# c_mu is strictly constrained to be a weak dark-sector coupling.
C_MU_DEFAULT = 2.74e-26

