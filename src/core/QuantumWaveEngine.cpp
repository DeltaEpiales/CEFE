#include "cefe/core/QuantumWaveEngine.hpp"
#include <iostream>
#include <cmath>

namespace cefe {
namespace core {

QuantumWaveEngine::QuantumWaveEngine(std::shared_ptr<geometry::CausalDiamondGrid> grid_ptr)
    : grid(grid_ptr), mass(1.0), hbar(1.0), current_dt(0.0) {
}

void QuantumWaveEngine::initialize_state(double target_t) {
    auto slice = grid->build_spatial_laplacian(target_t);
    slice_indices = slice.original_indices;
    
    std::size_t n = slice_indices.size();
    psi = Eigen::VectorXcd::Zero(n);
    potential = Eigen::VectorXd::Zero(n);
}

void QuantumWaveEngine::set_gaussian_packet(double x0, double y0, double sigma, double px, double py) {
    std::size_t n = slice_indices.size();
    const auto& points = grid->get_points();
    
    double norm_sum = 0.0;
    
    for (std::size_t i = 0; i < n; ++i) {
        const auto& p = points[slice_indices[i]];
        double dx = p.x - x0;
        double dy = p.y - y0;
        
        double exponent = -(dx * dx + dy * dy) / (4.0 * sigma * sigma);
        double phase = (px * dx + py * dy) / hbar;
        
        std::complex<double> val = std::exp(exponent) * std::polar(1.0, phase);
        psi(i) = val;
        
        norm_sum += std::norm(val);
    }
    
    // Normalize
    if (norm_sum > 0.0) {
        psi /= std::sqrt(norm_sum);
    }
}

void QuantumWaveEngine::set_double_slit_potential(double slit_width, double slit_separation, double barrier_thickness, double barrier_x) {
    std::size_t n = slice_indices.size();
    const auto& points = grid->get_points();
    
    double v_max = 1e6; // Infinite barrier approximation
    
    for (std::size_t i = 0; i < n; ++i) {
        const auto& p = points[slice_indices[i]];
        
        // Check if inside the barrier along x
        if (std::abs(p.x - barrier_x) < barrier_thickness / 2.0) {
            
            // Check for slits
            double y_upper = slit_separation / 2.0;
            double y_lower = -slit_separation / 2.0;
            
            bool in_slit = false;
            if (std::abs(p.y - y_upper) < slit_width / 2.0) in_slit = true;
            if (std::abs(p.y - y_lower) < slit_width / 2.0) in_slit = true;
            
            if (!in_slit) {
                potential(i) = v_max;
            }
        }
    }
}

void QuantumWaveEngine::set_atomic_potential(double Z, double softening) {
    std::size_t n = slice_indices.size();
    const auto& points = grid->get_points();
    
    for (std::size_t i = 0; i < n; ++i) {
        const auto& p = points[slice_indices[i]];
        double r2 = p.x * p.x + p.y * p.y + p.z * p.z;
        double r = std::sqrt(r2 + softening * softening);
        potential(i) = -Z / r;
    }
}

void QuantumWaveEngine::build_hamiltonian(double target_t) {
    auto slice = grid->build_spatial_laplacian(target_t);
    const auto& L = slice.laplacian;
    
    std::size_t n = slice_indices.size();
    H_sparse.resize(n, n);
    
    std::vector<Eigen::Triplet<std::complex<double>>> triplets;
    
    double coeff = -(hbar * hbar) / (2.0 * mass);
    
    for (int k = 0; k < L.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(L, k); it; ++it) {
            triplets.emplace_back(it.row(), it.col(), coeff * it.value());
        }
    }
    
    for (std::size_t i = 0; i < n; ++i) {
        triplets.emplace_back(i, i, potential(i));
    }
    
    H_sparse.setFromTriplets(triplets.begin(), triplets.end());
}

void QuantumWaveEngine::prepare_crank_nicolson(double dt) {
    current_dt = dt;
    std::size_t n = H_sparse.rows();
    
    Eigen::SparseMatrix<std::complex<double>> I(n, n);
    I.setIdentity();
    
    std::complex<double> factor(0.0, dt / (2.0 * hbar));
    
    CN_lhs = I + factor * H_sparse;
    CN_rhs = I - factor * H_sparse;
    
    // Explicitly enforce Dirichlet boundary conditions for barriers
    for (int k = 0; k < CN_lhs.outerSize(); ++k) {
        for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(CN_lhs, k); it; ++it) {
            int r = it.row();
            int c = it.col();
            if (potential(r) >= 1e5 || potential(c) >= 1e5) {
                if (r == c) {
                    it.valueRef() = 1.0;
                } else {
                    it.valueRef() = 0.0;
                }
            }
        }
    }
    
    for (int k = 0; k < CN_rhs.outerSize(); ++k) {
        for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(CN_rhs, k); it; ++it) {
            int r = it.row();
            int c = it.col();
            if (potential(r) >= 1e5 || potential(c) >= 1e5) {
                it.valueRef() = 0.0; // RHS is 0 so psi_new = 0 for boundary points
            }
        }
    }
    
    solver.compute(CN_lhs);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Decomposition failed for Crank-Nicolson!" << std::endl;
    }
}

void QuantumWaveEngine::step_forward() {
    if (current_dt == 0.0) {
        std::cerr << "Must call prepare_crank_nicolson before stepping!" << std::endl;
        return;
    }
    
    Eigen::VectorXcd rhs_vec = CN_rhs * psi;
    psi = solver.solve(rhs_vec);
}

double QuantumWaveEngine::get_total_probability() const {
    return psi.squaredNorm();
}

Eigen::VectorXd QuantumWaveEngine::get_probability_density() const {
    return psi.cwiseAbs2();
}

} // namespace core
} // namespace cefe
