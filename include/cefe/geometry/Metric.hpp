#pragma once

#include <cmath>

namespace cefe {
namespace geometry {

class Metric {
public:
    virtual ~Metric() = default;
    
    // Determines if a point is within the causal diamond of size R
    virtual bool is_inside_causal_diamond(double t, double x, double y, double z, double R) const = 0;
    
    // Computes the proper spatial distance between two nearby points
    virtual double get_spatial_distance(double x1, double y1, double z1, double x2, double y2, double z2) const = 0;
};

class MinkowskiMetric : public Metric {
public:
    bool is_inside_causal_diamond(double t, double x, double y, double z, double R) const override {
        return std::sqrt(x*x + y*y + z*z) <= R - std::abs(t);
    }
    
    double get_spatial_distance(double x1, double y1, double z1, double x2, double y2, double z2) const override {
        return std::sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
    }
};

class SchwarzschildMetric : public Metric {
private:
    double Rs; // Schwarzschild radius
public:
    SchwarzschildMetric(double Rs) : Rs(Rs) {}
    
    bool is_inside_causal_diamond(double t, double x, double y, double z, double R) const override {
        double r = std::sqrt(x*x + y*y + z*z);
        if (r <= Rs) return false; // Exclude interior of horizon
        
        // Approximate tortoise coordinate lightcone deformation
        double local_c = 1.0 - Rs / r; 
        return r <= R - std::abs(t) * local_c;
    }
    
    double get_spatial_distance(double x1, double y1, double z1, double x2, double y2, double z2) const override {
        double r1 = std::sqrt(x1*x1 + y1*y1 + z1*z1);
        double dr2 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
        
        if (r1 <= Rs) return std::sqrt(dr2);
        // Spatial distance is dilated near the horizon
        return std::sqrt(dr2) / std::sqrt(1.0 - Rs / r1);
    }
};

} // namespace geometry
} // namespace cefe
