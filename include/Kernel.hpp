#ifndef KERNEL_HPP
#define KERNEL_HPP

#include <cmath>

enum class IntegralType {
    S_INTEGRAL,  // Integral over second variable (s), with t fixed
    T_INTEGRAL   // Integral over first variable (t), with s fixed
};

class Kernel{
public:
    virtual double fun(double t, double s) const = 0;
    virtual double s_integral(double t_fixed, double start, double end) const = 0; // integral of second variable
    virtual double t_integral(double s_fixed, double start, double end) const = 0; // integral of first variable
    virtual ~Kernel() = default;

protected:
    /**
     * @brief Template method to check input parameters for integral functions
     * 
     * @tparam type The type of integral (S_INTEGRAL or T_INTEGRAL)
     * @param fixed The fixed variable (t_fixed for s_integral, s_fixed for t_integral)
     * @param start Start of integration interval
     * @param end End of integration interval
     * @return true if parameters are valid, false otherwise
     */
    template<IntegralType type>
    bool check_integral_params(double fixed, double start, double end) const {
        // Common checks for both types
        if (start < 0 || end < 0 || start >= end) {
            return false;
        }
        // Type-specific checks
        if constexpr (type == IntegralType::S_INTEGRAL) {
            // For s_integral: t_fixed should be greater than start
            // because we integrate s from start to end with s < t_fixed
            return fixed > start;
        } else if constexpr (type == IntegralType::T_INTEGRAL) {
            // For t_integral: s_fixed should be less than end
            // because we integrate t from start to end with t > s_fixed
            return fixed < end;
        }
        return false;
    }
    
};

class ExpKernel : public Kernel{
    double c = 1;
    double ro = 1;
public: // K(t,s)
    ExpKernel(double c, double ro) : c(c), ro(ro) {}

    double fun(double t, double s) const{
        if (s < 0 || t < 0 || s >= t) {
            return 0.0; // Return 0 for negative inputs
        }
        return c * exp(-ro * (t-s));
    }
    // integral of second variable
    double s_integral(double t_fixed, double start, double end) const{
        if (!check_integral_params<IntegralType::S_INTEGRAL>(t_fixed, start, end)) {
            return 0.0; // Return 0 for invalid intervals
        }
        // Adjust end to t_fixed if it exceeds t_fixed
        if (t_fixed < end) end = t_fixed;

        return c * (exp(-ro * (t_fixed-end)) - exp(-ro * (t_fixed - start))) / ro;
    } 
    // integral of first variable
    double t_integral(double s_fixed, double start, double end) const{
        if (!check_integral_params<IntegralType::T_INTEGRAL>(s_fixed, start, end)) {
            return 0.0; // Return 0 for invalid intervals
        }
        // Adjust start to s_fixed if it is less than s_fixed
        if (s_fixed > start) start = s_fixed;
        return c * (exp(-ro * (start-s_fixed)) - exp(-ro* (end-s_fixed))) / ro;
    }
};

class FracKernel : public Kernel{
    double c = 1;
    double alpha = 0.5;
public: // K(t,s)
    FracKernel(double c, double alpha) : c(c), alpha(alpha) {}
    double fun(double t, double s) const{
        if (s < 0 || t < 0 || s >= t) {
            return 0.0; // Return 0 for negative inputs
        }
        return c * std::pow(t - s, alpha - 1);
    }
    // integral of second variable
    double s_integral(double t_fixed, double start, double end) const{
        if (!check_integral_params<IntegralType::S_INTEGRAL>(t_fixed, start, end)) {
            return 0.0; // Return 0 for invalid intervals
        }
        // Adjust end to t_fixed if it exceeds t_fixed
        if (t_fixed < end) end = t_fixed;
        return c * (std::pow(t_fixed - start, alpha) - std::pow(t_fixed - end, alpha)) / alpha;
    }
    // integral of first variable
    double t_integral(double s_fixed, double start, double end) const{
        if (!check_integral_params<IntegralType::T_INTEGRAL>(s_fixed, start, end)) {
            return 0.0; // Return 0 for invalid intervals
        }
        // Adjust start to s_fixed if it is less than s_fixed
        if (s_fixed > start) start = s_fixed;
        return c * (std::pow(start - s_fixed, alpha) - std::pow(end - s_fixed, alpha)) / alpha;
    }
};

#endif // KERNEL_HPP