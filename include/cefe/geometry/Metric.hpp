#pragma once

#include <cmath>

namespace cefe {
namespace geometry {

class Metric {
public:
    virtual ~Metric() = default;
    
    // Determines if a point is within the causal diamond of size R
    virtual bool is_inside_causal_diamond(double t, double x, double y, double z, double R) const = 0;
    
    // Computes the proper spatial distance between two nearby points at time t
    virtual double get_spatial_distance(double t, double x1, double y1, double z1, double x2, double y2, double z2) const = 0;
};

class MinkowskiMetric : public Metric {
public:
    bool is_inside_causal_diamond(double t, double x, double y, double z, double R) const override {
        return std::sqrt(x*x + y*y + z*z) <= R - std::abs(t);
    }
    
    double get_spatial_distance(double t, double x1, double y1, double z1, double x2, double y2, double z2) const override {
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
    
    double get_spatial_distance(double t, double x1, double y1, double z1, double x2, double y2, double z2) const override {
        double r1 = std::sqrt(x1*x1 + y1*y1 + z1*z1);
        double dr2 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
        
        if (r1 <= Rs) return std::sqrt(dr2);
        // Spatial distance is dilated near the horizon
        return std::sqrt(dr2) / std::sqrt(1.0 - Rs / r1);
    }
};

class AntiDeSitterMetric : public Metric {
private:
    double L; // AdS curvature radius
public:
    AntiDeSitterMetric(double L) : L(L) {}
    
    bool is_inside_causal_diamond(double t, double x, double y, double z, double R) const override {
        double r = std::sqrt(x*x + y*y + z*z);
        // Light travels faster in coordinate r at large r in AdS
        double max_r = R * std::cosh(std::abs(t)/L); // simplistic causal bound
        return r <= max_r;
    }
    
    double get_spatial_distance(double t, double x1, double y1, double z1, double x2, double y2, double z2) const override {
        double r1 = std::sqrt(x1*x1 + y1*y1 + z1*z1);
        double dr2 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
        return std::sqrt(dr2) / std::sqrt(1.0 + (r1*r1)/(L*L));
    }
};

class FLRWMetric : public Metric {
public:
    FLRWMetric() {}
    
    bool is_inside_causal_diamond(double t, double x, double y, double z, double R) const override {
        double r = std::sqrt(x*x + y*y + z*z);
        // Matter-dominated expansion a(t) ~ t^(2/3)
        // Assume t_0 = 1 for current conformal time scale
        double a_t = std::pow(std::abs(t) + 1.0, 2.0/3.0);
        return r <= (R - std::abs(t)) / a_t;
    }
    
    double get_spatial_distance(double t, double x1, double y1, double z1, double x2, double y2, double z2) const override {
        double dr2 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
        double a_t = std::pow(std::abs(t) + 1.0, 2.0/3.0);
        return std::sqrt(dr2) * a_t;
    }
};

class KerrMetric : public Metric {
private:
    double M; // Mass
    double a; // Spin parameter (J/M)
public:
    KerrMetric(double M, double a) : M(M), a(a) {}
    
    bool is_inside_causal_diamond(double t, double x, double y, double z, double R) const override {
        double r = std::sqrt(x*x + y*y + z*z);
        double horizon = M + std::sqrt(M*M - a*a);
        if (r <= horizon) return false;
        
        // Approximation of lightcone dragging
        double local_c = 1.0 - (2.0 * M * r) / (r*r + a*a); 
        if (local_c < 0) local_c = 0.1; // Ergosphere simplistic bound
        return r <= R - std::abs(t) * local_c;
    }
    
    double get_spatial_distance(double t, double x1, double y1, double z1, double x2, double y2, double z2) const override {
        double r1 = std::sqrt(x1*x1 + y1*y1 + z1*z1);
        double dr2 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
        
        double horizon = M + std::sqrt(M*M - a*a);
        if (r1 <= horizon) return std::sqrt(dr2);
        
        double delta = r1*r1 - 2.0*M*r1 + a*a;
        double rho2 = r1*r1; // Equatorial approx z=0
        return std::sqrt(dr2) * std::sqrt(rho2 / delta);
    }
};

} // namespace geometry
} // namespace cefe
