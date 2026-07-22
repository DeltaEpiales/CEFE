#define _USE_MATH_DEFINES
#include "cefe/geometry/CausalDiamond.hpp"
#include <cmath>

namespace cefe {
namespace geometry {

CausalDiamondGrid::CausalDiamondGrid(double radius, double spacing, std::shared_ptr<Metric> m) 
    : R(radius), ds(spacing), metric(m), active_dof(0) {
    generate_grid();
}

void CausalDiamondGrid::generate_grid() {
    int steps = std::ceil(R / ds);
    
    for (int t_step = -steps; t_step <= steps; ++t_step) {
        double t = t_step * ds;
        double r_max = R - std::abs(t);
        
        if (r_max < 0) continue;
        
        int spatial_steps = std::ceil(r_max / ds);
        
        for (int x_step = -spatial_steps; x_step <= spatial_steps; ++x_step) {
            double x = x_step * ds;
            for (int y_step = -spatial_steps; y_step <= spatial_steps; ++y_step) {
                double y = y_step * ds;
                for (int z_step = -spatial_steps; z_step <= spatial_steps; ++z_step) {
                    double z = z_step * ds;
                    
                    if (metric->is_inside_causal_diamond(t, x, y, z, R)) {
                        // Very rough boundary estimate
                        bool is_boundary = (std::abs(std::sqrt(x*x + y*y + z*z) - r_max) <= ds / 2.0);
                        points.push_back({t, x, y, z, is_boundary, active_dof});
                        active_dof++;
                    }
                }
            }
        }
    }
}

double CausalDiamondGrid::get_area_bound(double subregion_radius) const {
    // Bekenstein Bound: S <= A / 4 (in Planck units G=1)
    double area = 4.0 * M_PI * subregion_radius * subregion_radius;
    return area / 4.0; 
}

CausalDiamondGrid::SpatialSlice CausalDiamondGrid::build_spatial_laplacian(double target_t) const {
    SpatialSlice slice;
    for (std::size_t i = 0; i < points.size(); ++i) {
        if (std::abs(points[i].t - target_t) < ds / 2.0) {
            slice.original_indices.push_back(i);
        }
    }
    
    std::size_t N = slice.original_indices.size();
    slice.laplacian = Eigen::SparseMatrix<double>(N, N);
    
    std::vector<Eigen::Triplet<double>> triplets;
    std::vector<double> diag(N, 0.0);
    
    for (std::size_t i = 0; i < N; ++i) {
        const auto& p1 = points[slice.original_indices[i]];
        
        for (std::size_t j = i + 1; j < N; ++j) {
            const auto& p2 = points[slice.original_indices[j]];
            double dist2 = (p1.x - p2.x)*(p1.x - p2.x) + 
                           (p1.y - p2.y)*(p1.y - p2.y) + 
                           (p1.z - p2.z)*(p1.z - p2.z);
                           
            // If nearest neighbor on grid
            if (std::abs(dist2 - ds * ds) < 1e-5) {
                double proper_ds = metric->get_spatial_distance(p1.x, p1.y, p1.z, p2.x, p2.y, p2.z);
                double inv_ds2 = 1.0 / (proper_ds * proper_ds);
                
                triplets.push_back(Eigen::Triplet<double>(i, j, -inv_ds2));
                triplets.push_back(Eigen::Triplet<double>(j, i, -inv_ds2));
                diag[i] += inv_ds2;
                diag[j] += inv_ds2;
            }
        }
    }
    
    for (std::size_t i = 0; i < N; ++i) {
        triplets.push_back(Eigen::Triplet<double>(i, i, diag[i]));
    }
    
    slice.laplacian.setFromTriplets(triplets.begin(), triplets.end());
    return slice;
}

} // namespace geometry
} // namespace cefe
