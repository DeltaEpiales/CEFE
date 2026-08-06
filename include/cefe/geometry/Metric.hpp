#pragma once

#include <cmath>
#include <algorithm>

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
        if (r <= Rs || R <= Rs) return false; // Exclude interior of horizon
        
        // Exact tortoise coordinate lightcone deformation
        double r_star = r + Rs * std::log(std::abs(r/Rs - 1.0));
        double R_star = R + Rs * std::log(std::abs(R/Rs - 1.0));
        
        return (R_star - r_star) >= std::abs(t);
    }
    
    double get_spatial_distance(double t, double x1, double y1, double z1, double x2, double y2, double z2) const override {
        // Simpson's rule integration along the path for proper distance
        double dist = 0.0;
        int steps = 10;
        double dx = (x2 - x1) / steps;
        double dy = (y2 - y1) / steps;
        double dz = (z2 - z1) / steps;
        double dl2 = dx*dx + dy*dy + dz*dz;
        
        for (int i = 0; i < steps; ++i) {
            double cx = x1 + (i + 0.5) * dx;
            double cy = y1 + (i + 0.5) * dy;
            double cz = z1 + (i + 0.5) * dz;
            double cr = std::sqrt(cx*cx + cy*cy + cz*cz);
            if (cr <= Rs) return std::sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
            
            double r_start = std::sqrt(std::pow(x1 + i*dx, 2) + std::pow(y1 + i*dy, 2) + std::pow(z1 + i*dz, 2));
            double r_end = std::sqrt(std::pow(x1 + (i+1)*dx, 2) + std::pow(y1 + (i+1)*dy, 2) + std::pow(z1 + (i+1)*dz, 2));
            double dr = r_end - r_start;
            
            double step_dist = std::sqrt(dl2 + (Rs / (cr - Rs)) * dr * dr);
            dist += step_dist;
        }
        return dist;
    }
};

class AntiDeSitterMetric : public Metric {
private:
    double L; // AdS curvature radius
public:
    AntiDeSitterMetric(double L) : L(L) {}
    
    bool is_inside_causal_diamond(double t, double x, double y, double z, double R) const override {
        double r = std::sqrt(x*x + y*y + z*z);
        // Exact causal bounds using AdS time
        double r_star = L * std::atan(r / L);
        double R_star = L * std::atan(R / L);
        return (R_star - r_star) >= std::abs(t);
    }
    
    double get_spatial_distance(double t, double x1, double y1, double z1, double x2, double y2, double z2) const override {
        double dist = 0.0;
        int steps = 10;
        double dx = (x2 - x1) / steps;
        double dy = (y2 - y1) / steps;
        double dz = (z2 - z1) / steps;
        double dl2 = dx*dx + dy*dy + dz*dz;
        
        for (int i = 0; i < steps; ++i) {
            double cx = x1 + (i + 0.5) * dx;
            double cy = y1 + (i + 0.5) * dy;
            double cz = z1 + (i + 0.5) * dz;
            double cr = std::sqrt(cx*cx + cy*cy + cz*cz);
            
            double r_start = std::sqrt(std::pow(x1 + i*dx, 2) + std::pow(y1 + i*dy, 2) + std::pow(z1 + i*dz, 2));
            double r_end = std::sqrt(std::pow(x1 + (i+1)*dx, 2) + std::pow(y1 + (i+1)*dy, 2) + std::pow(z1 + (i+1)*dz, 2));
            double dr = r_end - r_start;
            
            double step_dist = std::sqrt(std::max(0.0, dl2 - (cr*cr / (L*L + cr*cr)) * dr * dr));
            dist += step_dist;
        }
        return dist;
    }
};

class FLRWMetric : public Metric {
public:
    FLRWMetric() {}
    
    bool is_inside_causal_diamond(double t, double x, double y, double z, double R) const override {
        double r = std::sqrt(x*x + y*y + z*z);
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
        double horizon = M + std::sqrt(std::max(0.0, M*M - a*a));
        if (r <= horizon || R <= horizon) return false;
        
        // Approximate tortoise coordinate for Kerr
        double r_star = r + (2.0*M*horizon)/(horizon - (M - std::sqrt(std::max(0.0, M*M - a*a)))) * std::log(std::abs(r - horizon));
        double R_star = R + (2.0*M*horizon)/(horizon - (M - std::sqrt(std::max(0.0, M*M - a*a)))) * std::log(std::abs(R - horizon));
        
        return (R_star - r_star) >= std::abs(t);
    }
    
    double get_spatial_distance(double t, double x1, double y1, double z1, double x2, double y2, double z2) const override {
        double dist = 0.0;
        int steps = 10;
        double dx = (x2 - x1) / steps;
        double dy = (y2 - y1) / steps;
        double dz = (z2 - z1) / steps;
        
        double horizon = M + std::sqrt(std::max(0.0, M*M - a*a));
        
        for (int i = 0; i < steps; ++i) {
            double cx = x1 + (i + 0.5) * dx;
            double cy = y1 + (i + 0.5) * dy;
            double cz = z1 + (i + 0.5) * dz;
            double cr = std::sqrt(cx*cx + cy*cy + cz*cz);
            if (cr <= horizon) return std::sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
            
            double rho2 = cr*cr + a*a*cz*cz/(cr*cr); 
            double delta = cr*cr - 2.0*M*cr + a*a;
            if (delta < 1e-6) delta = 1e-6; // Avoid division by zero very close to horizon
            
            double r_start = std::sqrt(std::pow(x1 + i*dx, 2) + std::pow(y1 + i*dy, 2) + std::pow(z1 + i*dz, 2));
            double r_end = std::sqrt(std::pow(x1 + (i+1)*dx, 2) + std::pow(y1 + (i+1)*dy, 2) + std::pow(z1 + (i+1)*dz, 2));
            double dr = r_end - r_start;
            
            double dl2 = dx*dx + dy*dy + dz*dz;
            double step_dist = std::sqrt(std::max(0.0, dl2 + (rho2/delta - 1.0)*dr*dr));
            dist += step_dist;
        }
        return dist;
    }
};

} // namespace geometry
} // namespace cefe
