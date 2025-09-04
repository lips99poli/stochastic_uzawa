#ifndef STOCHASTIC_UZAWA_SOLVER_HPP
#define STOCHASTIC_UZAWA_SOLVER_HPP

#include <memory>
#include <string>
#include <stdexcept>
#include "CommonTypes.hpp"

// Include order: fundamental types first, then dependent classes
#include "Parameters.hpp"
#include "Kernel.hpp"
#include "Constraints.hpp"
#include "OUSimulator.hpp"
#include "LSMCR.hpp"
#include "NystromScheme.hpp"

class StochasticUzawaSolver {
private:
    // Parameters and related objects
    Parameters params;
    const NumericSchemeParams& numeric_params; // Reference to the numeric scheme parameters
    const OUParams& ou_params; // OU process parameters
    const ConstraintsParams& constraints_params; // Constraints parameters
    const KernelParams& kernel_params; // Kernel parameters

    // Time parameters
    const double time_delta; // Time step size
    const std::size_t N; // Number of time steps
    const Vector time_grid; // Time grid

    // Simulation parameter
    const std::size_t M; // Number of Monte Carlo paths

    // Control problem design
    std::unique_ptr<Kernel> kernel; // Kernel object
    Constraints constraints_obj; // 4 Constraints matrix
    MatVec constraints;
    
    // Signal related objects
    MatVec Gamma; // Expectation of phi
    Matrix alpha; // Conditional expectation of integral of finite variation part of the price
    MatVec R; // Conditional expectations of alpha
    Matrix price; // Price matrix

    // State variables
    std::unique_ptr<Matrix> u; // Control vector
    std::unique_ptr<Matrix> X_u; // State vector
    Matrix Z_u; // Transient price impact
    std::vector<std::unique_ptr<Matrix>> lambda; // 4 Lagrange multipliers
    
    // Cycle tools
    std::size_t n;
    std::vector<double> slackness; // Slackness variables for convergence check

    // Private helper methods
    void gradient_update(std::size_t n);
    double compute_slackness();

public:
    // Constructor that takes a Parameters object directly
    explicit StochasticUzawaSolver(const Parameters& p);

    // Simulate signal
    Matrix simulate_signal();

    // Implement the Uzawa algorithm to solve the control problem
    void solve();

    // Getter methods for the Interface class
    const Matrix& get_price() const;
    const Matrix& get_u() const;
    const Matrix& get_X() const;
    const Matrix& get_lambda1() const;
    const Matrix& get_lambda2() const;
    const Matrix& get_lambda3() const;
    const Matrix& get_lambda4() const;
    const Vector& get_time_grid() const;
    std::size_t iterations() const;
};

#endif // STOCHASTIC_UZAWA_SOLVER_HPP