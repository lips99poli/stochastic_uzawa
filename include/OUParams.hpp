#ifndef OU_PARAMS_HPP
#define OU_PARAMS_HPP
#include <cmath> // for M_PI

struct OUParams {
    // Ornstein-Uhlenbeck process parameters: dIt = (A(t)-k*It)dt + psi*dW_t where A(t) = theta * sin(omega*t + phi)
    double theta = 1;       // Mean reversion level
    double omega = 0;       // Frequency of the oscillation
    double phi = M_PI/2;    // Phase shift of the oscillation
    double k = 1;           // Mean reversion speed
    double psi = 0.2;       // Volatility of the process
    double T;               // Total time duration for the simulation
    double N;               // Time steps
    OUParams(double T = 1.0, double N = 1000) 
        : T(T), N(N) {} // Constructor with default values
};

#endif // OU_PARAMS_HPP