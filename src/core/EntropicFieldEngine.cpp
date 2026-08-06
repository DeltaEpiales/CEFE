#include "cefe/core/EntropicFieldEngine.hpp"
#include <iostream>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>

namespace cefe {
namespace core {

EntropicFieldEngine::EntropicFieldEngine(std::shared_ptr<geometry::CausalDiamondGrid> grid_ptr)
    : grid(grid_ptr), mass(1.0), boundary_condition(BoundaryType::DIRICHLET), interaction_coupling(0.0) {
}

void EntropicFieldEngine::set_mass(double m) {
    mass = m;
}

void EntropicFieldEngine::set_boundary_condition(BoundaryType bc) {
    boundary_condition = bc;
}

void EntropicFieldEngine::evolve_to_slice(double target_t, double temperature) {
    auto slice = grid->build_spatial_laplacian(target_t);
    slice_indices = slice.original_indices;
    
    std::size_t N = slice.laplacian.rows();
    if (N == 0) return;
    
    Eigen::MatrixXd W = -Eigen::MatrixXd(slice.laplacian);
    
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        W(i, i) += mass * mass;
    }
    
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(W);
    Eigen::VectorXd evals = es.eigenvalues();
    
    Eigen::VectorXd c_evals(N);
    Eigen::VectorXd p_evals(N);
    
    for(int i=0; i<N; ++i) {
        double w = std::sqrt(evals(i));
        double factor = 1.0;
        if (temperature > 0.0) {
            double beta = 1.0 / temperature;
            factor = 1.0 / std::tanh(beta * w / 2.0); // coth(beta w / 2)
        }
        c_evals(i) = 0.5 * (1.0 / w) * factor;
        p_evals(i) = 0.5 * w * factor;
    }
    
    covariance_C = es.eigenvectors() * c_evals.asDiagonal() * es.eigenvectors().transpose();
    covariance_P = es.eigenvectors() * p_evals.asDiagonal() * es.eigenvectors().transpose();
}

void EntropicFieldEngine::evolve_to_slice_sparse(double target_t, int num_modes, double temperature) {
    auto slice = grid->build_spatial_laplacian(target_t);
    slice_indices = slice.original_indices;
    
    std::size_t N = slice.laplacian.rows();
    if (N == 0) return;
    
    Eigen::SparseMatrix<double> W = -slice.laplacian;
    for (int k=0; k<W.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(W, k); it; ++it) {
            if (it.row() == it.col()) {
                it.valueRef() += mass * mass;
            }
        }
    }
    
    int ncv = std::min((int)N, 2 * num_modes + 1);
    num_modes = std::min((int)N - 1, num_modes);
    if (num_modes <= 0) return;
    
    Spectra::SparseSymShiftSolve<double> op(W);
    Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<double>> eigs(op, num_modes, ncv, 0.0);
    
    eigs.init();
    eigs.compute(Spectra::SortRule::LargestMagn);
    
    Eigen::VectorXd evals = eigs.eigenvalues();
    Eigen::MatrixXd evecs = eigs.eigenvectors();
    
    Eigen::VectorXd c_evals(num_modes);
    Eigen::VectorXd p_evals(num_modes);
    
    for(int i=0; i<num_modes; ++i) {
        double w = std::sqrt(evals(i));
        double factor = 1.0;
        if (temperature > 0.0) {
            double beta = 1.0 / temperature;
            factor = 1.0 / std::tanh(beta * w / 2.0);
        }
        c_evals(i) = 0.5 * (1.0 / w) * factor;
        p_evals(i) = 0.5 * w * factor;
    }
    
    covariance_C = evecs * c_evals.asDiagonal() * evecs.transpose();
    covariance_P = evecs * p_evals.asDiagonal() * evecs.transpose();
}

void EntropicFieldEngine::initialize_field_state(double target_t, double temperature) {
    auto slice = grid->build_spatial_laplacian(target_t);
    std::size_t N = slice.laplacian.rows();
    
    W_sparse = -slice.laplacian;
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
        int center_idx = 0;
        double min_dist = 1e9;
        const auto& points = grid->get_points();
        for (std::size_t i = 0; i < N; ++i) {
            const auto& p = points[slice.original_indices[i]];
            double dist = std::sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
            if (dist < min_dist) {
                min_dist = dist;
                center_idx = i;
            }
        }
        current_state.phi(center_idx) = 1.0; // Initial perturbation
    }
}

void EntropicFieldEngine::step_forward(double dt) {
    if (W_sparse.rows() == 0) return;
    
    // Leapfrog step with lambda phi^4 interaction
    Eigen::VectorXd force = - (W_sparse * current_state.phi);
    if (interaction_coupling > 0.0) {
        for(int i=0; i<force.size(); ++i) {
            double phi_val = current_state.phi(i);
            force(i) -= interaction_coupling * phi_val * phi_val * phi_val;
        }
    }
    
    current_state.pi += dt * force;
    current_state.phi += dt * current_state.pi;
}

double EntropicFieldEngine::get_field_energy() const {
    if (W_sparse.rows() == 0) return 0.0;
    double kinetic = 0.5 * current_state.pi.dot(current_state.pi);
    double potential = 0.5 * current_state.phi.dot(W_sparse * current_state.phi);
    
    double interaction_energy = 0.0;
    if (interaction_coupling > 0.0) {
        for(int i=0; i<current_state.phi.size(); ++i) {
            double phi_val = current_state.phi(i);
            interaction_energy += 0.25 * interaction_coupling * phi_val * phi_val * phi_val * phi_val;
        }
    }
    
    return kinetic + potential + interaction_energy;
}

double EntropicFieldEngine::compute_entanglement_entropy_indices(const std::vector<int>& subregion_indices) {
    std::size_t NA = subregion_indices.size();
    if (NA == 0) return 0.0;
    
    Eigen::MatrixXd C_A(NA, NA);
    Eigen::MatrixXd P_A(NA, NA);
    
    #pragma omp parallel for
    for (int i = 0; i < NA; ++i) {
        for (int j = 0; j < NA; ++j) {
            C_A(i, j) = covariance_C(subregion_indices[i], subregion_indices[j]);
            P_A(i, j) = covariance_P(subregion_indices[i], subregion_indices[j]);
        }
    }
    
    Eigen::MatrixXd CP = C_A * P_A;
    Eigen::EigenSolver<Eigen::MatrixXd> solver(CP);
    double entropy = 0.0;
    
    for (int k = 0; k < CP.rows(); ++k) {
        double nu2 = solver.eigenvalues()[k].real();
        double nu = std::sqrt(std::max(0.25, nu2)); // assure physical minimum
        
        double x1 = nu + 0.5;
        double x2 = std::max(0.0, nu - 0.5); // x2 can be perfectly zero
        
        if (x1 > 0.0) entropy += x1 * std::log(x1);
        if (x2 > 0.0) entropy -= x2 * std::log(x2);
    }
    
    return entropy;
}

double EntropicFieldEngine::compute_entanglement_entropy(double subregion_radius) {
    const auto& points = grid->get_points();
    std::vector<int> sub_A;
    
    for (std::size_t i = 0; i < slice_indices.size(); ++i) {
        const auto& p = points[slice_indices[i]];
        double r = std::sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
        if (r <= subregion_radius) {
            sub_A.push_back((int)i);
        }
    }
    return compute_entanglement_entropy_indices(sub_A);
}

} // namespace core
} // namespace cefe
