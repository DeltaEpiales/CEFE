#pragma once

#include <vector>
#include <cstddef>
#include <memory>
#include <Eigen/Sparse>
#include "cefe/geometry/Metric.hpp"

namespace cefe {
namespace geometry {

struct GridPoint {
    double t;
    double x;
    double y;
    double z;
    bool is_boundary;
    std::size_t index;
};

class CausalDiamondGrid {
private:
    double R;
    double ds;
    std::shared_ptr<Metric> metric;
    std::vector<GridPoint> points;
    std::size_t active_dof;

    void generate_grid();

public:
    CausalDiamondGrid(double radius, double spacing, std::shared_ptr<Metric> m = std::make_shared<MinkowskiMetric>());
    
    double get_radius() const { return R; }
    double get_spacing() const { return ds; }
    std::size_t get_dof() const { return active_dof; }
    const std::vector<GridPoint>& get_points() const { return points; }
    double get_area_bound(double subregion_radius) const;
    
    struct SpatialSlice {
        std::vector<std::size_t> original_indices;
        Eigen::SparseMatrix<double> laplacian;
    };

    SpatialSlice build_spatial_laplacian(double target_t) const;
};

} // namespace geometry
} // namespace cefe
