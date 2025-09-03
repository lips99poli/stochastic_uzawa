#ifndef STOCHASTIC_UZAWA_SOLVER_HPP
#define STOCHASTIC_UZAWA_SOLVER_HPP

#include <memory>
#include <iostream>
#include <fstream>
#include <limits>
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
    // Note that outputs u, X_u and lambda are stored in the heap
    
    // Cycle tools
    std::size_t n;
    std::vector<double> slackness; // Slackness variables for convergence check

    void gradient_update(std::size_t n) {
        double adaptive_learning = numeric_params.delta / std::pow(n + 1, numeric_params.beta);
        *lambda[0] += adaptive_learning * (constraints[0] - *u);
        *lambda[1] += adaptive_learning * (*u - constraints[1]);
        *lambda[2] += adaptive_learning * (constraints[2] - *X_u);
        *lambda[3] += adaptive_learning * (*X_u - constraints[3]);
        
        // Project lambda to non-negative values (KKT conditions)
        *lambda[0] = lambda[0]->cwiseMax(0.0);
        *lambda[1] = lambda[1]->cwiseMax(0.0);
        *lambda[2] = lambda[2]->cwiseMax(0.0);
        *lambda[3] = lambda[3]->cwiseMax(0.0);
    }

    auto compute_slackness() {
        MatVec slackness_m(4,Matrix::Zero(M, N));
        slackness_m[0] = (constraints[0] - *u).cwiseProduct(*lambda[0]);
        slackness_m[1] = (*u - constraints[1]).cwiseProduct(*lambda[1]);
        slackness_m[2] = (constraints[2] - *X_u).cwiseProduct(*lambda[2]);
        slackness_m[3] = (*X_u - constraints[3]).cwiseProduct(*lambda[3]);

        std::vector<Vector> slackness_avg_paths(4, Vector::Zero(N));
        // Each slackness variable is a vector containing the sum of all the rows of the corresponding slackness_m matrix
        for (std::size_t i = 0; i < slackness_m.size(); ++i) {
            slackness_avg_paths[i] = slackness_m[i].colwise().sum() / static_cast<double>(M);
        }
        
        // Integrate the slackness average paths
        for (std::size_t i = 0; i < slackness_avg_paths.size(); ++i) {
            //slackness[i] = std::abs(slackness_avg_paths[i].sum() * time_delta); // Left rectangle rule
            slackness[i] = slackness_avg_paths[i].sum() * time_delta; // Left rectangle rule
        }
        return *std::max_element(slackness.begin(), slackness.end()); // Return the maximum slackness iterator
    }


public:
    // Note that variables are of time length N, whereas multipliers are of time length N+1
    StochasticUzawaSolver(const std::string& filename = "data/Parameters.pot") :
        params(filename),

        // Split the params object into its components for easier access
        numeric_params(params.get_numeric_params()),
        ou_params(params.get_ou_params()),
        constraints_params(params.get_constraints_params()),
        kernel_params(params.get_kernel_params()),

        // Constraint class
        constraints_obj(constraints_params.X0, 
            constraints_params.u_min, constraints_params.u_max, 
            constraints_params.X_u_min, constraints_params.X_u_max, 
            constraints_params.terminal_liquidation, 
            constraints_params.stop_trad_price_lb, constraints_params.price_lb),

        // Time parameters
        time_delta(numeric_params.T / (numeric_params.N)),
        N(numeric_params.N), // Number of time steps
        time_grid(Vector::LinSpaced(N+1, 0, numeric_params.T)), //see (3.1) paper

        // Number of Monte Carlo paths
        M(numeric_params.M),

        // Control Variables
        u( std::make_unique<Matrix>(Matrix::Zero(M, N)) ),
        X_u( std::make_unique<Matrix>(Matrix::Constant(M, N, constraints_params.X0)) ),
        Z_u( Matrix::Zero(M, N) ),

        // Lagrange Multipliers
        lambda(4), // Initialize lambda vector with 4 elements
        slackness(std::vector<double>(4, std::numeric_limits<double>::infinity())) // Initialize slackness variables to infinity

    {
        // Create the appropriate kernel based on the type string
        if (kernel_params.kernel_type == "exp") {
            if (kernel_params.kernel_parameters.size() != 2) {
                throw std::invalid_argument("Exponential kernel requires exactly 2 parameters: {c, ro}");
            }
            kernel = std::make_unique<ExpKernel>(kernel_params.kernel_parameters[0], kernel_params.kernel_parameters[1]);
        } else if (kernel_params.kernel_type == "frac") {
            if (kernel_params.kernel_parameters.size() != 2) {
                throw std::invalid_argument("Fractional kernel requires exactly 2 parameters: {c, alpha}");
            }
            kernel = std::make_unique<FracKernel>(kernel_params.kernel_parameters[0], kernel_params.kernel_parameters[1]);
        } else {
            throw std::invalid_argument("Invalid kernel type: " + kernel_params.kernel_type + 
                                      ". Supported types are: 'exp' and 'frac'");
        }
        
        // Initialize lambda vector elements
        for(std::size_t i = 0; i < 4; ++i) {
            lambda[i] = std::make_unique<Matrix>(Matrix::Zero(M, N));
        }
        
        // Initialize gamma to vector of size M of matrices (N,N)
        Gamma.resize(M, Matrix::Zero(N, N));
    }

    // Simulate signal
    Matrix simulate_signal(){
        // Launch the simulation of the OU process
        OUSimulator ou_simulator(ou_params, time_grid, M);
        // Fill the alpha matrix based on the OU parameters
        alpha = ou_simulator.getAlpha();
        // Fill the R matrix based on the OU parameters
        R = ou_simulator.getR();
        price = ou_simulator.getPrice();
        constraints = constraints_obj.compute_constraints(price);
        return price;
    }


    // Implement the Uzawa algorithm to solve the control problem
    // This is a placeholder for the actual Uzawa algorithm implementation
    // The algorithm will iteratively update u, X_u, and lambda based on the R, U, and L matrices

    // Cicle steps are the following:
    // 1. Gradient update lambda
    // 2. Estimate conditional expectations through Least Squares Monte Carlo Regression (LSMCR)
    // 3. Update control variable u and state variables X_u, Z_u through Nystrom Scheme and left rectangle rule
    // 4. Update Slackness Variables to check if convergence is reached
    void solve(){
        // Initialize LSMCR structure
        LSMCR lsmcr(numeric_params.d1, numeric_params.d2, N, time_delta, M, alpha, Z_u, X_u, lambda);

        // Initialize Nystrom Scheme structure
        NystromScheme nystrom_scheme(M, N, time_grid, constraints_params.X0, *kernel, u, Z_u, X_u, R, Gamma);

        // Initialize slackness check variable and check it is correctly initialized
        auto max_slackness = *std::max_element(slackness.begin(), slackness.end());

        // Main loop for the Uzawa algorithm
        n = 0; // Iteration counter
        while (max_slackness > numeric_params.epsilon && n < numeric_params.D){
            // 1. Gradient update lambda
            gradient_update(n);

            // 2. Estimate conditional expectations through Least Squares Monte Carlo Regression (LSMCR)
            lsmcr.update_regressor();
            Gamma = lsmcr.estimate_gamma();
            
            // 3. Update control variable u and state variables X_u, Z_u through Nystrom Scheme and left rectangle rule
            nystrom_scheme.nystrom_update();

            // 4. Update Slackness Variables to check if convergence is reached
            max_slackness = compute_slackness();
           
            // Increment the iteration counter
            ++n;
        }
    }
};


#endif // STOCHASTIC_UZAWA_SOLVER_HPP