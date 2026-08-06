#include "DynamicTensorEngine.hpp"
#include <iostream>
#include <cmath>
#include <limits>

namespace cefe {
namespace dynamic_tensors {

DynamicTensorEngine::DynamicTensorEngine(std::shared_ptr<geometry::CausalDiamondGrid> grid)
    : grid_ptr_(grid) {
}

Eigen::Matrix4d DynamicTensorEngine::get_local_metric(std::size_t node_index) const {
    // Simplified getter for the background metric at the node.
    Eigen::Matrix4d g = Eigen::Matrix4d::Zero();
    g(0,0) = -1.0;
    g(1,1) = 1.0;
    g(2,2) = 1.0;
    g(3,3) = 1.0;
    return g;
}

void DynamicTensorEngine::compute_stress_energy_tensor(const std::vector<double>& field_phi, double mass) {
    std::size_t num_nodes = grid_ptr_->get_points().size();
    T_mu_nu_field_.resize(num_nodes, Eigen::Matrix4d::Zero());

    const auto& nodes = grid_ptr_->get_points();

    for (std::size_t i = 0; i < num_nodes; ++i) {
        // 1. Calculate discrete gradients del_mu phi
        Eigen::Vector4d del_phi = Eigen::Vector4d::Zero();
        
        // Find nearest neighbor to approximate spatial gradient
        double min_dist = std::numeric_limits<double>::max();
        std::size_t nearest = i;
        
        for (std::size_t j = 0; j < num_nodes; ++j) {
            if (i == j) continue;
            double dx = nodes[j].x - nodes[i].x;
            double dy = nodes[j].y - nodes[i].y;
            double dz = nodes[j].z - nodes[i].z;
            double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (dist > 1e-9 && dist < min_dist) {
                min_dist = dist;
                nearest = j;
            }
        }
        
        if (nearest != i) {
            double d_phi = field_phi[nearest] - field_phi[i];
            // Approximate spatial radial gradient
            double ri = std::sqrt(nodes[i].x*nodes[i].x + nodes[i].y*nodes[i].y + nodes[i].z*nodes[i].z);
            double rn = std::sqrt(nodes[nearest].x*nodes[nearest].x + nodes[nearest].y*nodes[nearest].y + nodes[nearest].z*nodes[nearest].z);
            double dr = rn - ri;
            if (std::abs(dr) > 1e-9) {
                del_phi(1) = d_phi / dr; // Radial gradient assigned to spatial index 1
            }
        }

        Eigen::Matrix4d g = get_local_metric(i);
        
        // 2. Kinetic term: del^lambda phi del_lambda phi
        double kinetic = 0.0;
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                double g_inv = (g(mu, nu) != 0.0) ? 1.0 / g(mu, nu) : 0.0;
                if (mu == nu) {
                    kinetic += g_inv * del_phi(mu) * del_phi(nu);
                }
            }
        }

        // 3. Potential term: m^2 phi^2
        double potential = mass * mass * field_phi[i] * field_phi[i];
        
        // 4. Construct T_mu_nu
        for (int mu = 0; mu < 4; ++mu) {
            for (int nu = 0; nu < 4; ++nu) {
                T_mu_nu_field_[i](mu, nu) = del_phi(mu) * del_phi(nu) - 0.5 * g(mu, nu) * (kinetic + potential);
            }
        }
    }
}

void DynamicTensorEngine::linearized_einstein_update(double time_step, double G_constant) {
    if (T_mu_nu_field_.empty()) {
        return;
    }

    auto& nodes = grid_ptr_->get_mutable_points();
    
    // Use 3.14159265358979323846 directly instead of M_PI for MSVC
    double coupling = 8.0 * 3.14159265358979323846 * G_constant * time_step;

    for (std::size_t i = 0; i < nodes.size(); ++i) {
        if (i >= T_mu_nu_field_.size()) {
            break;
        }
        double energy_density = T_mu_nu_field_[i](0, 0); // T_00
        
        double warp_factor = coupling * energy_density;
        if (warp_factor > 0.1) warp_factor = 0.1;
        
        // Scale x, y, z proportionally inwards
        nodes[i].x -= warp_factor * nodes[i].x;
        nodes[i].y -= warp_factor * nodes[i].y;
        nodes[i].z -= warp_factor * nodes[i].z;
    }
    
    std::cout << "[DynamicTensorEngine] Linearized Einstein Update applied. Lattice Warped." << std::endl;
}

} // namespace dynamic_tensors
} // namespace cefe
