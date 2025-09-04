#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <string>
#include <iostream>
#include "CommonTypes.hpp"

struct NumericSchemeParams{
    double T; // Total time duration for the simulation
    std::size_t N; // T/N is the time step size
    std::size_t D; // Max iterations
    std::size_t M; // Number of sample paths
    double epsilon; // Tolerance for convergence
    double delta; // Adaptive learning
    double beta; // Adaptive learning rate
    std::size_t d1; // Degree for LSMCR first regression
    std::size_t d2; // Degree for LSMCR second regression
};

struct OUParams {
    // Price process
    double sigma; // Volatility of the price process
    double S0; // Initial price level
    // Ornstein-Uhlenbeck process parameters: dIt = (A(t)-k*It)dt + psi*dW_t where A(t) = theta * sin(omega*t + phi)
    double I0;         // Initial value of the OU process
    double theta;       // Mean reversion level
    double omega;       // Frequency of the oscillation
    double phi;    // Phase shift of the oscillation
    double k;           // Mean reversion speed
    double psi;       // Volatility of the process
    double T;               // Total time duration for the simulation
    std::size_t N;               // T/N is the time step size
    int seed1; // Seed for OU process generation
    int seed2; // Seed for price process generation
};

struct ConstraintsParams{
    double X0; // Initial inventory
    double u_min; // Minimum trading rate
    double u_max;  // Maximum trading rate
    double X_u_min; // Minimum inventory
    double X_u_max;  // Maximum inventory
    bool terminal_liquidation; // Whether to enforce terminal liquidation constraint
    bool stop_trad_price_lb; // Whether to stop trading if price hits lower bound
    double price_lb; // Price lower bound for stop trading
};

struct KernelParams{
    std::string kernel_type;
    std::vector<double> kernel_parameters;
};

/**
 * @brief Convert degrees to radians
 * @param degrees Angle in degrees (0, 45, 90, 180, 270, 360, etc.)
 * @return Angle in radians
 */
inline double degrees_to_radians(double degrees);

class Parameters{
private:
    std::string config_file;
    NumericSchemeParams numeric_params;
    OUParams ou_params;
    ConstraintsParams constraints_params;
    KernelParams kernel_params;

public:
    /**
     * @brief Constructor that reads parameters from a file
     * @param filename Path to the parameter file (default: "data/Parameters.pot")
     */
    Parameters(const std::string& filename = "data/Parameters.pot");

    /**
     * @brief Read all parameters from the configuration file
     */
    void read_parameters();

    /**
     * @brief Print all parameters for debugging
     */
    void print_parameters() const;

    /**
     * @brief Get the numeric scheme parameters
     * @return Reference to NumericSchemeParams
     */
    const NumericSchemeParams& get_numeric_params() const;

    /**
     * @brief Get the OU parameters
     * @return Reference to OUParams
     */
    const OUParams& get_ou_params() const;

    /**
     * @brief Get the constraints parameters
     * @return Reference to ConstraintsParams
     */
    const ConstraintsParams& get_constraints_params() const;

    /**
     * @brief Get the kernel parameters
     * @return Reference to KernelParams
     */
    const KernelParams& get_kernel_params() const;
};

#endif // PARAMETERS_HPP