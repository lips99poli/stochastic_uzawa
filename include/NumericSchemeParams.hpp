#ifndef NUMERIC_SCHEME_PARAMS_HPP

#define NUMERIC_SCHEME_PARAMS_HPP

#include <string>
#include "GetPot"
#include <stdexcept>

struct NumericSchemeParams{
    double T=0; // Total time duration for the simulation
    std::size_t N = 0; // T/N is the time step size
    std::size_t D = 0; // Max iterations
    std::size_t M = 0; // Number of sample paths
    double epsilon=0; // Tolerance for convergence
    double delta=0; // Adaptive learning
    double beta=0; // Adaptive learning rate
};

NumericSchemeParams load_params(const std::string& file_name){
    GetPot datafile(file_name.c_str());
    NumericSchemeParams params;
    params.T = datafile("T", 1);
    params.N = datafile("N", 100);
    params.D = datafile("D", 10000);
    params.M = datafile("M", 1000);
    params.epsilon = datafile("epsilon", 1e-6);
    params.delta = datafile("delta", 0.1);
    params.beta = datafile("beta", 0.9);

    if (params.N == 0 || params.D == 0 || params.M == 0) {
        throw std::runtime_error("Invalid parameters: N, D, and M must be greater than zero.");
    }
    return params;
};

#endif