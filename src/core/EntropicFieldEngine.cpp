#include "cefe/core/EntropicFieldEngine.hpp"
#include <iostream>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace cefe {
namespace core {

EntropicFieldEngine::EntropicFieldEngine(std::shared_ptr<geometry::CausalDiamondGrid> grid_ptr)
    : grid(grid_ptr), mass(1.0), boundary_condition(BoundaryType::DIRICHLET) {
}

void EntropicFieldEngine::set_mass(double m) {
    mass = m;
}

void EntropicFieldEngine::set_boundary_condition(BoundaryType bc) {
    boundary_condition = bc;
}

void EntropicFieldEngine::evolve_to_slice(double target_t) {
    auto slice = grid->build_spatial_laplacian(target_t);
    slice_indices = slice.original_indices;
    
    std::size_t N = slice.laplacian.rows();
    if (N == 0) return;
    
    Eigen::MatrixXd W = Eigen::MatrixXd(slice.laplacian);
    
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        W(i, i) += mass * mass;
    }
    
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(W);
    Eigen::VectorXd sqrt_evals = es.eigenvalues().cwiseSqrt();
    Eigen::MatrixXd M = es.eigenvectors() * sqrt_evals.asDiagonal() * es.eigenvectors().transpose();
    
    Eigen::VectorXd inv_sqrt_evals = es.eigenvalues().cwiseInverse().cwiseSqrt();
    Eigen::MatrixXd M_inv = es.eigenvectors() * inv_sqrt_evals.asDiagonal() * es.eigenvectors().transpose();
    
    covariance_C = 0.5 * M_inv;
    covariance_P = 0.5 * M;
}

void EntropicFieldEngine::initialize_field_state(double target_t) {
    auto slice = grid->build_spatial_laplacian(target_t);
    std::size_t N = slice.laplacian.rows();
    
    W_sparse = slice.laplacian;
    for (int k=0; k<W_sparse.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(W_sparse, k); it; ++it) {
            if (it.row() == it.col()) {
                it.valueRef() += mass * mass;
            }
        }
    }
    
    current_state.phi = Eigen::VectorXd::Zero(N);
    current_state.pi = Eigen::VectorXd::Zero(N);
    if (N > 0) {
        current_state.phi(N / 2) = 1.0; // Initial perturbation
    }
}

void EntropicFieldEngine::step_forward(double dt) {
    if (W_sparse.rows() == 0) return;
    
    // Leapfrog step
    current_state.pi -= dt * (W_sparse * current_state.phi);
    current_state.phi += dt * current_state.pi;
}

double EntropicFieldEngine::get_field_energy() const {
    if (W_sparse.rows() == 0) return 0.0;
    double kinetic = 0.5 * current_state.pi.dot(current_state.pi);
    double potential = 0.5 * current_state.phi.dot(W_sparse * current_state.phi);
    return kinetic + potential;
}

double EntropicFieldEngine::compute_entanglement_entropy(double subregion_radius) {
    const auto& points = grid->get_points();
    std::vector<std::size_t> sub_A;
    
    for (std::size_t i = 0; i < slice_indices.size(); ++i) {
        const auto& p = points[slice_indices[i]];
        double r = std::sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
        if (r <= subregion_radius) {
            sub_A.push_back(i);
        }
    }
    
    std::size_t NA = sub_A.size();
    if (NA == 0) return 0.0;
    
    Eigen::MatrixXd C_A(NA, NA);
    Eigen::MatrixXd P_A(NA, NA);
    
    #pragma omp parallel for
    for (int i = 0; i < NA; ++i) {
        for (int j = 0; j < NA; ++j) {
            C_A(i, j) = covariance_C(sub_A[i], sub_A[j]);
            P_A(i, j) = covariance_P(sub_A[i], sub_A[j]);
        }
    }
    
    Eigen::MatrixXd CP = C_A * P_A;
    Eigen::EigenSolver<Eigen::MatrixXd> solver(CP);
    double entropy = 0.0;
    
    for (int k = 0; k < CP.rows(); ++k) {
        double nu2 = solver.eigenvalues()[k].real();
        if (nu2 < 0.25) nu2 = 0.25; 
        double nu = std::sqrt(nu2);
        
        double x1 = nu + 0.5;
        double x2 = nu - 0.5;
        if (x2 < 1e-12) x2 = 1e-12; 
        
        entropy += x1 * std::log(x1) - x2 * std::log(x2);
    }
    
    return entropy;
}

} // namespace core
} // namespace cefe
