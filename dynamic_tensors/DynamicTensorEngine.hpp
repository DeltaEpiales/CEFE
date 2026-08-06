#ifndef CEFE_DYNAMIC_TENSORS_DYNAMIC_TENSOR_ENGINE_HPP
#define CEFE_DYNAMIC_TENSORS_DYNAMIC_TENSOR_ENGINE_HPP

#include <vector>
#include <Eigen/Dense>
#include "cefe/geometry/CausalDiamond.hpp"
#include <Eigen/StdVector>

namespace cefe {
namespace dynamic_tensors {

class DynamicTensorEngine {
public:
    DynamicTensorEngine(std::shared_ptr<geometry::CausalDiamondGrid> grid);

    // Calculates T_mu_nu = del_mu phi del_nu phi - 1/2 g_mu_nu (del^lambda phi del_lambda phi + m^2 phi^2)
    void compute_stress_energy_tensor(const std::vector<double>& field_phi, double mass);

    // Warps the underlying CausalDiamondGrid by solving the Linearized Einstein Equations:
    // delta G_mu_nu = 8 * pi * G * T_mu_nu
    void linearized_einstein_update(double time_step, double G_constant = 1.0);

private:
    std::shared_ptr<geometry::CausalDiamondGrid> grid_ptr_;
    
    // Stores the discrete stress energy tensor components for each spatial node
    std::vector<Eigen::Matrix4d, Eigen::aligned_allocator<Eigen::Matrix4d>> T_mu_nu_field_;
    
    // Helper to extract local metric tensor
    Eigen::Matrix4d get_local_metric(std::size_t node_index) const;
};

} // namespace dynamic_tensors
} // namespace cefe

#endif // CEFE_DYNAMIC_TENSORS_DYNAMIC_TENSOR_ENGINE_HPP
