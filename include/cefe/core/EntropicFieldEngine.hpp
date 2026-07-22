#pragma once

#include "cefe/geometry/CausalDiamond.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <memory>

namespace cefe {
namespace core {

enum class BoundaryType {
    DIRICHLET,
    NEUMANN
};

struct FieldState {
    Eigen::VectorXd phi;
    Eigen::VectorXd pi;
};

class EntropicFieldEngine {
private:
    std::shared_ptr<geometry::CausalDiamondGrid> grid;
    Eigen::MatrixXd covariance_C;
    Eigen::MatrixXd covariance_P;
    std::vector<std::size_t> slice_indices;
    double mass;
    BoundaryType boundary_condition;
    
    // Dynamic field state
    FieldState current_state;
    Eigen::SparseMatrix<double> W_sparse;

public:
    EntropicFieldEngine(std::shared_ptr<geometry::CausalDiamondGrid> grid_ptr);
    
    void set_mass(double m);
    void set_boundary_condition(BoundaryType bc);
    
    void evolve_to_slice(double target_t);
    double compute_entanglement_entropy(double subregion_radius);
    
    void initialize_field_state(double target_t);
    void step_forward(double dt);
    
    int get_state_size() const { return covariance_C.rows(); }
    double get_field_energy() const;
};

} // namespace core
} // namespace cefe
