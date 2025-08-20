#ifndef KERNEL_HPP
#define KERNEL_HPP

#include <cmath>

class Kernel{
    double c = 1;
    double ro = 1;
public: // K(t,s)
    Kernel(double c, double ro) : c(c), ro(ro) {}

    double fun(double t, double s) const{
        if (s < 0 || t < 0) {
            return 0.0; // Return 0 for negative inputs
        }
        return (t>s)? c * exp(-ro * (t-s)) : 0.0;
    }
    // integral of second variable
    double s_integral(double t_fixed, double start, double end) const{
        if (start < 0 || end < 0 || start >= end) {
            return 0.0; // Return 0 for invalid intervals
        }
        return c * (exp(-ro * (t_fixed-end)) - exp(-ro * (t_fixed - start))) / ro;
    } 
    // integral of first variable
    double t_integral(double s_fixed, double start, double end) const{
        if (start < 0 || end < 0 || start >= end) {
            return 0.0; // Return 0 for invalid intervals
        }
        return c * (exp(-ro * (start-s_fixed)) - exp(-ro* (end-s_fixed))) / ro;
    }
};

#endif // KERNEL_HPP