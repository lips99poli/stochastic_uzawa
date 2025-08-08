#ifndef NUMERIC_SCHEME_PARAMS_HPP

#define NUMERIC_SCHEME_PARAMS_HPP

#include <string>
struct NumericSchemeParams{
    double T=0; // Total time duration for the simulation
    std::size_t N = 0; // T/N is the time step size
    std::size_t D = 0; // Max iterations
    std::size_t M = 0; // Number of sample paths
    double epsilon=0; // Tolerance for convergence
    double delta=0; // Adaptive learning
    double beta=0; // Adaptive learning rate
};

NumericSchemeParams load_params(const std::string& file_name);

#endif