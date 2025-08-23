#ifndef KERNEL_HPP
#define KERNEL_HPP

#include <cmath>

class Kernel{
    double c = 1;
    double ro = 1;
public: // K(t,s)
    Kernel(double c, double ro) : c(c), ro(ro) {}

    double fun(double t, double s) const{
        if (s < 0 || t < 0 || s >= t) {
            return 0.0; // Return 0 for negative inputs
        }
        return c * exp(-ro * (t-s));
    }
    // integral of second variable
    double s_integral(double t_fixed, double start, double end) const{
        if (start < 0 || end < 0 || start >= end || t_fixed <= start) {
            return 0.0; // Return 0 for invalid intervals
        }
        // Adjust end to t_fixed if it exceeds t_fixed
        if (t_fixed < end) end = t_fixed;

        return c * (exp(-ro * (t_fixed-end)) - exp(-ro * (t_fixed - start))) / ro;
    } 
    // integral of first variable
    double t_integral(double s_fixed, double start, double end) const{
        if (start < 0 || end < 0 || start >= end || s_fixed >= end) {
            return 0.0; // Return 0 for invalid intervals
        }
        // Adjust start to s_fixed if it is less than s_fixed
        if (s_fixed > start) start = s_fixed;
        return c * (exp(-ro * (start-s_fixed)) - exp(-ro* (end-s_fixed))) / ro;
    }
};

#endif // KERNEL_HPP