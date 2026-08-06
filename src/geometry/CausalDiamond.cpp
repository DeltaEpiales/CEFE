#define _USE_MATH_DEFINES
#include "cefe/geometry/CausalDiamond.hpp"
#include <cmath>
#include <unordered_map>
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
    
    std::unordered_map<long long, std::size_t> spatial_hash;
    auto get_hash = [](int x, int y, int z) -> long long {
        return (static_cast<long long>(x + 100000) << 40) | 
               (static_cast<long long>(y + 100000) << 20) | 
               (static_cast<long long>(z + 100000));
    };

    for (std::size_t i = 0; i < N; ++i) {
        const auto& p = points[slice.original_indices[i]];
        int xi = std::round(p.x / ds);
        int yi = std::round(p.y / ds);
        int zi = std::round(p.z / ds);
        spatial_hash[get_hash(xi, yi, zi)] = i;
    }
    
    for (std::size_t i = 0; i < N; ++i) {
        const auto& p1 = points[slice.original_indices[i]];
        int xi = std::round(p1.x / ds);
        int yi = std::round(p1.y / ds);
        int zi = std::round(p1.z / ds);
        
        // Only check positive offsets to avoid double-counting edges
        int offsets[3][3] = {{1,0,0}, {0,1,0}, {0,0,1}};
        for (int d = 0; d < 3; ++d) {
            long long neighbor_key = get_hash(xi + offsets[d][0], yi + offsets[d][1], zi + offsets[d][2]);
            auto it = spatial_hash.find(neighbor_key);
            if (it != spatial_hash.end()) {
                std::size_t j = it->second;
                const auto& p2 = points[slice.original_indices[j]];
                
                double proper_ds = metric->get_spatial_distance(p1.t, p1.x, p1.y, p1.z, p2.x, p2.y, p2.z);
                double inv_ds2 = 1.0 / (proper_ds * proper_ds);
                
                triplets.push_back(Eigen::Triplet<double>(i, j, inv_ds2));
                triplets.push_back(Eigen::Triplet<double>(j, i, inv_ds2));
                diag[i] += inv_ds2;
                diag[j] += inv_ds2;
            }
        }
    }
    
    for (std::size_t i = 0; i < N; ++i) {
        triplets.push_back(Eigen::Triplet<double>(i, i, -diag[i]));
    }
    
    slice.laplacian.setFromTriplets(triplets.begin(), triplets.end());
    return slice;
}

} // namespace geometry
} // namespace cefe
