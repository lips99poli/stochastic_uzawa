#ifndef STOCHASTIC_UZAWA_SOLVER_HPP
#define STOCHASTIC_UZAWA_SOLVER_HPP

#include <vector>
#include <memory>
#include <iostream>
#include <fstream>
#include <limits>
#include <Eigen/Dense>
//create aliases for Eigen types
using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using Ref = Eigen::Ref<Matrix>; // Specify the template argument for Eigen::Ref

using Mat_Vec = std::vector<Matrix>; // Map to store conditional expectations for each time step couple (i,j)


#include "NumericSchemeParams.hpp"
#include "Kernel.hpp"
#include "OUParams.hpp"
#include "OUSimulator.hpp"
#include "LSMCR.hpp"
#include "NystromScheme.hpp"


class StochasticUzawaSolver {
    private:
    std::ostream* output_stream; // Output stream for writing results

    public:

    double X0;

    const NumericSchemeParams params; // Reference to the numeric scheme parameters
    const double time_delta; // Time step size
    const std::size_t N; // Number of time steps
    const Vector time_grid; // Time grid

    const Kernel& kernel; // Kernel object

    const OUParams ou_params; // OU process parameters
    
    Matrix alpha;
    Mat_Vec R; // Signal matrix
    
    Mat_Vec Gamma; // Expectation of phi

    const Mat_Vec constraints; // 4 Constraints matrix

    // The output u, X_u and lamda are stored in the heap
    std::unique_ptr<Matrix> u; // Control vector
    std::unique_ptr<Matrix> X_u; // State vector
    Matrix Z_u; // Transient price impact
    std::vector<std::unique_ptr<Matrix>> lambda; // 4 Lagrange multipliers
    
    std::vector<double> slackness; // Slackness variables for convergence check

    std::string mode;

    Mat_Vec compute_constraints(const int M) { //questo èda estendere: i constraint possono essere di vario tipo: vorrei implementare i 4 del paper e poi lasciare lo user libero di dare una csustom function che calcola ai_t come funzione del prezzo nei tempi da 0 a t
        // Compute the constraints matrix: each of the M rows corresponds to the constraints path for the same row in the price simulation matrix
        Mat_Vec constraints(4, Matrix::Zero(M, N)); // 4 constraints, each of size MxN
        for (int m = 0; m < M; ++m) {
            for (int i = 0; i < N; ++i) {
                constraints[0](m, i) = -5.0;
                constraints[1](m, i) = 5.0;
                if(i<N-1){
                    constraints[2](m, i) = -100.0;
                    constraints[3](m, i) = 100.0;
                }
            }
        }
        return constraints;
    }

    void gradient_update(std::size_t n) {
        double adaptive_learning = params.delta / std::pow(n + 1, params.beta);
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
        Mat_Vec slackness_m(4,Matrix::Zero(params.M, N));
        slackness_m[0] = (constraints[0] - *u).cwiseProduct(*lambda[0]);
        slackness_m[1] = (*u - constraints[1]).cwiseProduct(*lambda[1]);
        slackness_m[2] = (constraints[2] - *X_u).cwiseProduct(*lambda[2]);
        slackness_m[3] = (*X_u - constraints[3]).cwiseProduct(*lambda[3]);

        std::vector<Vector> slackness_avg_paths(4, Vector::Zero(N));
        // Each slackness variable is a vector containing the sum of all the rows of the corresponding slackness_m matrix
        for (std::size_t i = 0; i < slackness_m.size(); ++i) {
            slackness_avg_paths[i] = slackness_m[i].colwise().sum() / static_cast<double>(params.M);
        }
        
        // Debug print slackness_avg_paths
        if(mode=="verbose"){
            *output_stream << "Debug - Slackness average paths:" << std::endl;
            for (std::size_t i = 0; i < slackness_avg_paths.size(); ++i) {
                *output_stream << "slackness_avg_paths[" << i << "]: " << slackness_avg_paths[i].transpose() << std::endl;
            }
            *output_stream << std::endl;
        }
        
        // Integrate the slackness average paths
        for (std::size_t i = 0; i < slackness_avg_paths.size(); ++i) {
            //slackness[i] = std::abs(slackness_avg_paths[i].sum() * time_delta); // Left rectangle rule
            slackness[i] = slackness_avg_paths[i].sum() * time_delta; // Left rectangle rule
        }

        if(mode=="verbose"){
            // print the slackness values
            *output_stream << "Current slackness values: ";
            for (const auto& s : slackness) {
                *output_stream << s << " ";
            }
            *output_stream << std::endl;
        }

        return *std::max_element(slackness.begin(), slackness.end()); // Return the maximum slackness iterator
    }

    void print_variables(std::size_t iteration) {
        *output_stream << "=== ITERATION " << iteration << " VARIABLES ===" << std::endl;
        
        *output_stream << "u matrix:" << std::endl;
        *output_stream << *(u) << std::endl << std::endl;
        
        *output_stream << "X_u matrix:" << std::endl;
        *output_stream << *(X_u) << std::endl << std::endl;
        
        *output_stream << "Lambda matrices:" << std::endl;
        for (size_t i = 0; i < lambda.size(); ++i) {
            *output_stream << "Lambda[" << i << "]:" << std::endl;
            *output_stream << *(lambda[i]) << std::endl << std::endl;
        }
        
        *output_stream << "Gamma matrices:" << std::endl;
        for (size_t i = 0; i < Gamma.size(); ++i) {
            *output_stream << "Gamma[" << i << "]:" << std::endl;
            *output_stream << Gamma[i] << std::endl << std::endl;
        }
        
        *output_stream << "========================" << std::endl << std::endl;
    }


    public:
    // Note that variables are of time length N, whereas multipliers are of time length N+1
    StochasticUzawaSolver(const NumericSchemeParams& params, const Kernel& k, const OUParams& ou_params, std::ostream* output = &std::cout, const double X0 = -5, std::string mode = "verbose") :
        // Output stream
        output_stream(output)
        ,mode(mode)
        // Numeric scheme parameters
        ,params(params)
        ,time_delta(params.T / (params.N))
        ,N(params.N) // Number of time steps
        ,time_grid(Vector::LinSpaced(params.N+1, 0, params.T)) //see (3.1) paper
        ,kernel(k)
        ,ou_params(ou_params)
        // Control problem parameters
        ,constraints(compute_constraints(params.M))
        // Control Variables and Lagrange Multipliers
        ,u( std::make_unique<Matrix>(Matrix::Zero(params.M, params.N)) )
        ,X0(X0)
        ,X_u( std::make_unique<Matrix>(Matrix::Constant(params.M, params.N,X0)) )
        ,Z_u( Matrix::Zero(params.M, params.N) )
        ,lambda(4) // Initialize lambda vector with 4 elements
        ,slackness(std::vector<double>(4, std::numeric_limits<double>::infinity())) // Initialize slackness variables to infinity
        // Cycle tools

    {
        // Initialize lambda vector elements
        for(std::size_t i = 0; i < 4; ++i) {
            lambda[i] = std::make_unique<Matrix>(Matrix::Zero(params.M, params.N));
        }
        
        //launch the simulation of the OU process
        OUSimulator ou_simulator(ou_params, time_grid, params.M);
        // Fill the alpha matrix based on the OU parameters
        alpha = ou_simulator.getAlpha();
        // Fill the R matrix based on the OU parameters
        R = ou_simulator.getR();
        // Initialize gamma to vector of size M of matrices (N,N)
        Gamma.resize(params.M, Matrix::Zero(params.N, params.N));
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
        const std::size_t d1 = 1; // Upper bound for the sum of the degrees of the Laguerre polynomials used in the first regression
        const std::size_t d2 = 1; // Upper bound for the sum of the degrees of the Laguerre polynomials used in the second regression
        LSMCR lsmcr(d1, d2, N, time_delta, params.M, alpha, Z_u, X_u, lambda);

        // Initialize Nystrom Scheme structure
        NystromScheme nystrom_scheme(params.M, params.N, time_grid, X0, kernel, u, Z_u, X_u, R, Gamma);
        if(mode=="verbose"){
            nystrom_scheme.print_operators(*output_stream);
        }

        // Initialize slackness check variable and check it is correctly initialized
        auto max_slackness = *std::max_element(slackness.begin(), slackness.end());

        // Main loop for the Uzawa algorithm
        std::size_t n = 0; // Iteration counter
        while (max_slackness > params.epsilon && n < params.D){
            // 1. Gradient update lambda
            gradient_update(n);

            // 2. Estimate conditional expectations through Least Squares Monte Carlo Regression (LSMCR)
            lsmcr.update_regressor();
            Gamma = lsmcr.estimate_gamma();
            
            // 3. Update control variable u and state variables X_u, Z_u through Nystrom Scheme and left rectangle rule
            nystrom_scheme.nystrom_update();

            // 4. Update Slackness Variables to check if convergence is reached
            max_slackness = compute_slackness();

            if(mode=="big"){
                if(n%10==0) {
                    // Debug print slackness
                    *output_stream << "Iteration " << n << ", slackness values: ";
                    for (size_t i = 0; i < slackness.size(); ++i) {
                        *output_stream << slackness[i];
                        if (i < slackness.size() - 1) *output_stream << ", ";
                    }
                    *output_stream << std::endl;
                    
                    // Print last column of first 10 rows of X_u
                    *output_stream << "Last column of first 10 rows of X_u: ";
                    for (int m = 0; m < std::min(10, (int)params.M); ++m) {
                        *output_stream << (*X_u)(m, params.N-1);
                        if (m < std::min(10, (int)params.M) - 1) *output_stream << " ";
                    }
                    *output_stream << std::endl;
                }
            }
            if(mode=="verbose"){
            // Print all variables for verbose output
                print_variables(n);
            }
           
            // Increment the iteration counter
            ++n;
        }
        // print iterations
        *output_stream << "Performed " << n << " iterations." << std::endl;
        
    }
};


#endif // STOCHASTIC_UZAWA_SOLVER_HPP