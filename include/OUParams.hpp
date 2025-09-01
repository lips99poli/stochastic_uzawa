#ifndef PRICE_PARAMS_HPP
#define PRICE_PARAMS_HPP
#include <cmath> // for M_PI
#include "GetPot"

struct OUParams {
    // Price process
    double sigma = 0.1; // Volatility of the price process
    double S0 = 100; // Initial price level
    // Ornstein-Uhlenbeck process parameters: dIt = (A(t)-k*It)dt + psi*dW_t where A(t) = theta * sin(omega*t + phi)
    double I0 = -2;         // Initial value of the OU process
    double theta = -20;       // Mean reversion level
    double omega = 0;       // Frequency of the oscillation
    double phi = M_PI/2;    // Phase shift of the oscillation
    double k = 1;           // Mean reversion speed
    double psi = 4;       // Volatility of the process
    double T;               // Total time duration for the simulation
    std::size_t N;               // T/N is the time step size
    OUParams(double T, double N) 
        : T(T), N(N) {} // Constructor with default values
};

OUParams OU_load_params(double T, std::size_t N, const std::string& file_name){
    GetPot datafile(file_name.c_str());
    OUParams params(T,N);
    params.sigma = datafile("sigma", 1);
    params.S0 = datafile("S0", 100);
    params.I0 = datafile("I0", 10000);
    params.theta = datafile("theta", 1000);
    params.omega = datafile("omega", 1e-6);
    params.phi = datafile("phi", 0.1);
    params.k = datafile("k", 0.9);
    params.psi = datafile("psi", 0.9);
    return params;
};

#endif // OU_PARAMS_HPP