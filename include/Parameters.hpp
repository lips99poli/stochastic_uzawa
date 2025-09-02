#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "GetPot"

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using Ref = Eigen::Ref<Matrix>; // Specify the template argument for Eigen::Ref

using Mat_Vec = std::vector<Matrix>; // Map to store conditional expectations for each time step couple (i,j)

struct NumericSchemeParams{
    double T; // Total time duration for the simulation
    std::size_t N; // T/N is the time step size
    std::size_t D; // Max iterations
    std::size_t M; // Number of sample paths
    double epsilon; // Tolerance for convergence
    double delta; // Adaptive learning
    double beta; // Adaptive learning rate
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
};

struct ConstraintsParams{
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
    std::vector<double> kernel_params;
};

/**
 * @brief Convert degrees to radians
 * @param degrees Angle in degrees (0, 45, 90, 180, 270, 360, etc.)
 * @return Angle in radians
 */
double degrees_to_radians(double degrees) {
    return degrees * M_PI / 180.0;
}

class Parameters{
private:
    std::string config_file;

public:
    NumericSchemeParams numeric_params;
    OUParams ou_params;
    ConstraintsParams constraints_params;

    /**
     * @brief Constructor that reads parameters from a file
     * @param filename Path to the parameter file (default: "data/Parameters.pot")
     */
    Parameters(const std::string& filename = "data/Parameters.pot") : config_file(filename) {
        read_parameters();
    }

    /**
     * @brief Read all parameters from the configuration file
     */
    void read_parameters() {
        GetPot input_file(config_file.c_str());

        // Read Numeric Scheme Parameters
        numeric_params.T = input_file("T", 1.0);
        numeric_params.N = input_file("N", 100);
        numeric_params.D = input_file("D", 50);
        numeric_params.M = input_file("M", 10);
        numeric_params.epsilon = input_file("epsilon", 1e-3);
        numeric_params.delta = input_file("delta", 1);
        numeric_params.beta = input_file("beta", 0.8);
        
        // Read Ornstein-Uhlenbeck and Price Parameters
        ou_params.sigma = input_file("sigma", 0.0);
        ou_params.S0 = input_file("S0", 100.0);
        ou_params.I0 = input_file("I0", 1.0);
        ou_params.theta = input_file("theta", 0.0);
        ou_params.omega = input_file("omega", 0.0);
        
        // Read phi in degrees and convert to radians
        double phi_degrees = input_file("phi", 90.0);  // Default to 90 degrees (M_PI/2)
        ou_params.phi = degrees_to_radians(phi_degrees);
        
        ou_params.k = input_file("k", 0.0);
        ou_params.psi = input_file("psi", 0.0);
        
        // Set T and N for OU params (these might be computed or set elsewhere)
        ou_params.T = numeric_params.T;
        ou_params.N = numeric_params.N;
        
        // Read Constraints Parameters
        constraints_params.X0 = input_file("X0", -100.0);
        constraints_params.u_min = input_file("u_min", -100.0);
        constraints_params.u_max = input_file("u_max", 100.0);
        constraints_params.X_u_min = input_file("X_u_min", -1e6);
        constraints_params.X_u_max = input_file("X_u_max", 1e6);
        constraints_params.terminal_liquidation = static_cast<bool>(input_file("terminal_liquidation", 0));
        constraints_params.stop_trad_price_lb = static_cast<bool>(input_file("stop_trad_price_lb", 0));
        constraints_params.price_lb = input_file("price_lb", 0.0);
    }

    /**
     * @brief Set the total time T for the simulation
     * @param T Total time duration
     */
    void set_total_time(double T) {
        numeric_params.T = T;
        ou_params.T = T;
    }

    /**
     * @brief Print all parameters for debugging
     */
    void print_parameters() const {
        std::cout << "=== Numeric Scheme Parameters ===" << std::endl;
        std::cout << "T: " << numeric_params.T << std::endl;
        std::cout << "N: " << numeric_params.N << std::endl;
        std::cout << "D: " << numeric_params.D << std::endl;
        std::cout << "M: " << numeric_params.M << std::endl;
        std::cout << "epsilon: " << numeric_params.epsilon << std::endl;
        std::cout << "delta: " << numeric_params.delta << std::endl;
        std::cout << "beta: " << numeric_params.beta << std::endl;

        std::cout << "\n=== Ornstein-Uhlenbeck and Price Parameters ===" << std::endl;
        std::cout << "sigma: " << ou_params.sigma << std::endl;
        std::cout << "S0: " << ou_params.S0 << std::endl;
        std::cout << "I0: " << ou_params.I0 << std::endl;
        std::cout << "theta: " << ou_params.theta << std::endl;
        std::cout << "omega: " << ou_params.omega << std::endl;
        std::cout << "phi: " << ou_params.phi << " radians (" << ou_params.phi * 180.0 / M_PI << " degrees)" << std::endl;
        std::cout << "k: " << ou_params.k << std::endl;
        std::cout << "psi: " << ou_params.psi << std::endl;

        std::cout << "\n=== Constraints Parameters ===" << std::endl;
        std::cout << "u_min: " << constraints_params.u_min << std::endl;
        std::cout << "u_max: " << constraints_params.u_max << std::endl;
        std::cout << "X_u_min: " << constraints_params.X_u_min << std::endl;
        std::cout << "X_u_max: " << constraints_params.X_u_max << std::endl;
        std::cout << "terminal_liquidation: " << constraints_params.terminal_liquidation << std::endl;
        std::cout << "stop_trad_price_lb: " << constraints_params.stop_trad_price_lb << std::endl;
        std::cout << "price_lb: " << constraints_params.price_lb << std::endl;
    }

    /**
     * @brief Get the numeric scheme parameters
     * @return Reference to NumericSchemeParams
     */
    const NumericSchemeParams& get_numeric_params() const {
        return numeric_params;
    }

    /**
     * @brief Get the OU parameters
     * @return Reference to OUParams
     */
    const OUParams& get_ou_params() const {
        return ou_params;
    }

    /**
     * @brief Get the constraints parameters
     * @return Reference to ConstraintsParams
     */
    const ConstraintsParams& get_constraints_params() const {
        return constraints_params;
    }
};

#endif // PARAMETERS_HPP