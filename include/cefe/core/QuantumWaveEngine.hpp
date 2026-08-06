#pragma once

#include "cefe/geometry/CausalDiamond.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <memory>
#include <complex>

namespace cefe {
namespace core {

class QuantumWaveEngine {
private:
    std::shared_ptr<geometry::CausalDiamondGrid> grid;
    Eigen::VectorXcd psi; 
    Eigen::VectorXd potential; 
    double mass;
    double hbar;
    
    Eigen::SparseMatrix<std::complex<double>> H_sparse;
    std::vector<std::size_t> slice_indices;
    
    Eigen::SparseMatrix<std::complex<double>> CN_lhs;
    Eigen::SparseMatrix<std::complex<double>> CN_rhs;
    Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> solver;
    double current_dt;

public:
    QuantumWaveEngine(std::shared_ptr<geometry::CausalDiamondGrid> grid_ptr);
    
    void set_mass(double m) { mass = m; }
    void set_hbar(double h) { hbar = h; }
    
    void initialize_state(double target_t);
    void set_gaussian_packet(double x0, double y0, double sigma, double px, double py);
    void set_double_slit_potential(double slit_width, double slit_separation, double barrier_thickness, double barrier_x);
    void set_atomic_potential(double Z, double softening);
    
    void build_hamiltonian(double target_t);
    void prepare_crank_nicolson(double dt);
    void step_forward();
    
    double get_total_probability() const;
    Eigen::VectorXd get_probability_density() const;
    
    const std::vector<std::size_t>& get_slice_indices() const { return slice_indices; }
};

} // namespace core
} // namespace cefe
