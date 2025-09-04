#include "Kernel.hpp"

// ExpKernel implementation
ExpKernel::ExpKernel(double c, double ro) : c(c), ro(ro) {}

double ExpKernel::fun(double t, double s) const {
    if (s < 0 || t < 0 || s >= t) {
        return 0.0; // Return 0 for negative inputs
    }
    return c * exp(-ro * (t-s));
}

double ExpKernel::s_integral(double t_fixed, double start, double end) const {
    if (!check_integral_params<IntegralType::S_INTEGRAL>(t_fixed, start, end)) {
        return 0.0; // Return 0 for invalid intervals
    }
    // Adjust end to t_fixed if it exceeds t_fixed
    if (t_fixed < end) end = t_fixed;

    return c * (exp(-ro * (t_fixed-end)) - exp(-ro * (t_fixed - start))) / ro;
}

double ExpKernel::t_integral(double s_fixed, double start, double end) const {
    if (!check_integral_params<IntegralType::T_INTEGRAL>(s_fixed, start, end)) {
        return 0.0; // Return 0 for invalid intervals
    }
    // Adjust start to s_fixed if it is less than s_fixed
    if (s_fixed > start) start = s_fixed;
    return c * (exp(-ro * (start-s_fixed)) - exp(-ro* (end-s_fixed))) / ro;
}

// FracKernel implementation
FracKernel::FracKernel(double c, double alpha) : c(c), alpha(alpha) {}

double FracKernel::fun(double t, double s) const {
    if (s < 0 || t < 0 || s >= t) {
        return 0.0; // Return 0 for negative inputs
    }
    return c * std::pow(t - s, alpha - 1);
}

double FracKernel::s_integral(double t_fixed, double start, double end) const {
    if (!check_integral_params<IntegralType::S_INTEGRAL>(t_fixed, start, end)) {
        return 0.0; // Return 0 for invalid intervals
    }
    // Adjust end to t_fixed if it exceeds t_fixed
    if (t_fixed < end) end = t_fixed;
    return c * (std::pow(t_fixed - start, alpha) - std::pow(t_fixed - end, alpha)) / alpha;
}

double FracKernel::t_integral(double s_fixed, double start, double end) const {
    if (!check_integral_params<IntegralType::T_INTEGRAL>(s_fixed, start, end)) {
        return 0.0; // Return 0 for invalid intervals
    }
    // Adjust start to s_fixed if it is less than s_fixed
    if (s_fixed > start) start = s_fixed;
    return c * (std::pow(start - s_fixed, alpha) - std::pow(end - s_fixed, alpha)) / alpha;
}
