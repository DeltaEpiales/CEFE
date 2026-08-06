#include "cefe/core/LatticeQFTEngine.hpp"
#include <cmath>
#include <iostream>

namespace cefe {
namespace core {

LatticeQFTEngine::LatticeQFTEngine(std::shared_ptr<geometry::CausalDiamondGrid> grid)
    : grid_(grid), mass_(1.0), lambda_(0.0) {
}

void LatticeQFTEngine::set_mass(double m) {
    mass_ = m;
}

void LatticeQFTEngine::set_coupling(double lambda) {
    lambda_ = lambda;
}

void LatticeQFTEngine::initialize_state() {
    std::cout << "[LatticeQFTEngine] Using NEW Leapfrog Engine!" << std::endl;
    auto slice = grid_->build_spatial_laplacian(0.0);
    slice_indices = slice.original_indices;
    int N = slice_indices.size();
    
    phi_current_ = Eigen::VectorXd::Zero(N);
    pi_current_ = Eigen::VectorXd::Zero(N);
}

void LatticeQFTEngine::set_gaussian_packet(double x0, double y0, double sigma, double amplitude, double phase) {
    int N = slice_indices.size();
    const auto& points = grid_->get_points();
    
    for (int i = 0; i < N; ++i) {
        const auto& p = points[slice_indices[i]];
        double dx = p.x - x0;
        double dy = p.y - y0;
        double r2 = dx*dx + dy*dy;
        
        double envelope = amplitude * std::exp(-r2 / (2.0 * sigma * sigma));
        
        phi_current_[i] += envelope;
    }
}

void LatticeQFTEngine::set_opposing_packets(double x_offset, double sigma, double amplitude, double momentum) {
    int N = slice_indices.size();
    const auto& points = grid_->get_points();
    
    double v = momentum / mass_;
    
    for (int i = 0; i < N; ++i) {
        const auto& p = points[slice_indices[i]];
        
        // Packet 1 (Right moving)
        double dx1 = p.x - (-x_offset);
        double dy1 = p.y - 0.0;
        double r1_2 = dx1*dx1 + dy1*dy1;
        double env1 = amplitude * std::exp(-r1_2 / (2.0 * sigma * sigma));
        
        // Packet 2 (Left moving)
        double dx2 = p.x - (x_offset);
        double dy2 = p.y - 0.0;
        double r2_2 = dx2*dx2 + dy2*dy2;
        double env2 = amplitude * std::exp(-r2_2 / (2.0 * sigma * sigma));
        
        phi_current_[i] += env1 + env2;
        
        // Analytical time derivative pi = d(phi)/dt
        double d_env1_dt = env1 * v * dx1 / (sigma * sigma);
        double d_env2_dt = -env2 * v * dx2 / (sigma * sigma);
        
        pi_current_[i] += d_env1_dt + d_env2_dt;
    }
}

void LatticeQFTEngine::build_laplacian() {
    auto slice = grid_->build_spatial_laplacian(0.0);
    // The spatial laplacian built by CausalDiamondGrid already approximates +nabla^2
    // (negative diagonals, positive off-diagonals).
    laplacian_ = slice.laplacian; // Now laplacian_ = +nabla^2
}

void LatticeQFTEngine::step_forward(double dt) {
    int N = phi_current_.size();
    const auto& points = grid_->get_points();
    
    double h = grid_->get_spacing();
    double phi_max = phi_current_.cwiseAbs().maxCoeff();
    double omega_max_sq = 12.0 / (h * h) + (mass_ * mass_) + 0.5 * std::abs(lambda_) * (phi_max * phi_max);
    double dt_stable = 2.0 / std::sqrt(omega_max_sq);
    double dt_actual = std::min(0.0005, dt_stable * 0.8); // 0.8 safety factor
    
    int num_substeps = std::ceil(dt / dt_actual);
    double sub_dt = dt / num_substeps;
    
    for (int step = 0; step < num_substeps; ++step) {
        Eigen::VectorXd laplacian_phi = laplacian_ * phi_current_;
        Eigen::VectorXd force = Eigen::VectorXd::Zero(N);
        
        for (int i = 0; i < N; ++i) {
            const auto& p = points[slice_indices[i]];
            if (!p.is_boundary) {
                double phi = phi_current_[i];
                // EOM: ddot{phi} = nabla^2 phi - m^2 phi - (lambda/6) phi^3
                // laplacian_ is now +nabla^2, so laplacian_phi[i] = (nabla^2 phi)_i
                force[i] = laplacian_phi[i] - (mass_ * mass_) * phi - (lambda_ / 6.0) * (phi * phi * phi);
                // Velocity Verlet half-step p_{n+1/2} = p_n + dt/2 F_n
                pi_current_[i] += (sub_dt / 2.0) * force[i];
                // Full step phi_{n+1} = phi_n + dt p_{n+1/2}
                phi_current_[i] += sub_dt * pi_current_[i];
            }
        }
        
        // New force F_{n+1} for second half-step
        laplacian_phi = laplacian_ * phi_current_;
        for (int i = 0; i < N; ++i) {
            const auto& p = points[slice_indices[i]];
            if (!p.is_boundary) {
                double phi = phi_current_[i];
                double force_next = laplacian_phi[i] - (mass_ * mass_) * phi - (lambda_ / 6.0) * (phi * phi * phi);
                // Velocity Verlet second half-step p_{n+1} = p_{n+1/2} + dt/2 F_{n+1}
                pi_current_[i] += (sub_dt / 2.0) * force_next;
            }
        }
    }
}

std::vector<double> LatticeQFTEngine::get_field_amplitude() const {
    return std::vector<double>(phi_current_.data(), phi_current_.data() + phi_current_.size());
}

double LatticeQFTEngine::get_total_energy() const {
    double energy = 0.0;
    
    // laplacian_ is now +nabla^2. Gradient energy = (1/2)(nabla phi)^2 = -(1/2) phi nabla^2 phi
    Eigen::VectorXd lap_phi = laplacian_ * phi_current_;
    const auto& points = grid_->get_points();
    
    for (int i = 0; i < phi_current_.size(); ++i) {
        if (!points[slice_indices[i]].is_boundary) {
            double phi = phi_current_[i];
            double kinetic = 0.5 * pi_current_[i] * pi_current_[i];
            
            // Gradient energy: (1/2)(nabla phi)^2 = -(1/2) phi * (nabla^2 phi)
            double gradient = -0.5 * phi * lap_phi[i];
            double potential = 0.5 * (mass_ * mass_) * (phi * phi) + (lambda_ / 24.0) * (phi * phi * phi * phi);
            
            energy += kinetic + gradient + potential;
        }
    }
    return energy;
}

} // namespace core
} // namespace cefe
