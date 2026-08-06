#pragma once

#include "cefe/geometry/CausalDiamond.hpp"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>

#include <memory>

namespace cefe {
namespace core {

class LatticeQFTEngine {
public:
    LatticeQFTEngine(std::shared_ptr<geometry::CausalDiamondGrid> grid);
    
    // Physics Parameters
    void set_mass(double m);
    void set_coupling(double lambda);
    
    // Initialization
    void initialize_state();
    void set_gaussian_packet(double x0, double y0, double sigma, double amplitude, double phase);
    void set_opposing_packets(double x_offset, double sigma, double amplitude, double momentum);
    
    // Precomputation
    void build_laplacian();
    
    // Time Evolution
    void step_forward(double dt);
    
    // State Retrieval
    std::vector<double> get_field_amplitude() const;
    double get_total_energy() const;

private:
    std::shared_ptr<geometry::CausalDiamondGrid> grid_;
    std::vector<std::size_t> slice_indices;
    
    double mass_;
    double lambda_;
    
    Eigen::VectorXd phi_current_;
    Eigen::VectorXd pi_current_;
    Eigen::SparseMatrix<double> laplacian_;
};

} // namespace core
} // namespace cefe
