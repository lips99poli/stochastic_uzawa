#include "NumericSchemeParams.hpp"
#include "GetPot"
#include <stdexcept>

NumericSchemeParams load_params(const std::string& path){
    GetPot datafile(path.c_str());
    NumericSchemeParams params;
    params.T = datafile("T", 1);
    params.N = datafile("N", 101);
    params.D = datafile("D", 10000);
    params.M = datafile("M", 1000);
    params.epsilon = datafile("epsilon", 1e-6);
    params.delta = datafile("delta", 0.1);
    params.beta = datafile("beta", 0.9);

    if (params.N == 0 || params.D == 0 || params.M == 0) {
        throw std::runtime_error("Invalid parameters: N, D, and M must be greater than zero.");
    }
    return params;
}