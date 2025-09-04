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

class ExpKernel : public Kernel {
private:
    double c;
    double ro;

public:
    ExpKernel(double c, double ro);
    double fun(double t, double s) const override;
    double s_integral(double t_fixed, double start, double end) const override;
    double t_integral(double s_fixed, double start, double end) const override;
};

class FracKernel : public Kernel {
private:
    double c;
    double alpha;

public:
    FracKernel(double c, double alpha);
    double fun(double t, double s) const override;
    double s_integral(double t_fixed, double start, double end) const override;
    double t_integral(double s_fixed, double start, double end) const override;
};

#endif // KERNEL_HPP