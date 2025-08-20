#ifndef KERNEL_HPP
#define KERNEL_HPP

#include <cmath>

class Kernel{
    double c = 1;
    double ro = 1;
public: // K(t,s)
    Kernel(double c, double ro) : c(c), ro(ro) {}
    double fun(double t, double s) const;
    double s_integral(double t_fixed, double start, double end) const; // integral of second variable
    double t_integral(double s_fixed, double start, double end) const; // integral of first variable
};

#endif // KERNEL_HPP