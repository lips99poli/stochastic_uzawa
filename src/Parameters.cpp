#include "Parameters.hpp"
#include "GetPot"

inline double degrees_to_radians(double degrees) {
    return degrees * M_PI / 180.0;
}

Parameters::Parameters(const std::string& filename) : config_file(filename) {
    read_parameters();
}

void Parameters::read_parameters() {
    GetPot input_file(config_file.c_str());

    // Read Numeric Scheme Parameters
    numeric_params.T = input_file("T", 1.0);
    numeric_params.N = input_file("N", 100);
    numeric_params.D = input_file("D", 50);
    numeric_params.M = input_file("M", 10);
    numeric_params.epsilon = input_file("epsilon", 1e-3);
    numeric_params.delta = input_file("delta", 1.0);
    numeric_params.beta = input_file("beta", 0.8);
    numeric_params.d1 = input_file("d1", 1); // Degree for LSMCR first regression
    numeric_params.d2 = input_file("d2", 2); // Degree for LSMCR second regression
    
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
    ou_params.seed1 = input_file("seed1", -1); // Default to -1 for random seed
    ou_params.seed2 = input_file("seed2", -1); // Default to -1 for random seed
    
    // Read Constraints Parameters
    constraints_params.X0 = input_file("X0", -100.0);
    constraints_params.u_min = input_file("u_min", -100.0);
    constraints_params.u_max = input_file("u_max", 100.0);
    constraints_params.X_u_min = input_file("X_u_min", -1e6);
    constraints_params.X_u_max = input_file("X_u_max", 1e6);
    constraints_params.terminal_liquidation = static_cast<bool>(input_file("terminal_liquidation", 0));
    constraints_params.stop_trad_price_lb = static_cast<bool>(input_file("stop_trad_price_lb", 0));
    constraints_params.price_lb = input_file("price_lb", 0.0);

    // Read Kernel Parameters
    kernel_params.kernel_type = input_file("kernel_type", "exp");
    int num_kernel_params = 2;
    kernel_params.kernel_parameters.clear();
    for (int i = 0; i < num_kernel_params; ++i) {
        std::string param_name = "par" + std::to_string(i+1);
        double param_value = input_file(param_name.c_str(), 1.0); // Default to 1.0
        kernel_params.kernel_parameters.push_back(param_value);
    }
}

void Parameters::print_parameters() const {
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

void Parameters::update_price_parameters(const Parameters& other_params) {
    // Update price-related parameters from the other Parameters object
    const auto& other_numeric = other_params.get_numeric_params();
    const auto& other_ou = other_params.get_ou_params();
    
    // Update numeric parameters that affect price simulation
    numeric_params.T = other_numeric.T;
    numeric_params.N = other_numeric.N;
    numeric_params.M = other_numeric.M;
    
    // Update all OU/price parameters using assignment operator
    ou_params = other_ou;
}

const NumericSchemeParams& Parameters::get_numeric_params() const {
    return numeric_params;
}

const OUParams& Parameters::get_ou_params() const {
    return ou_params;
}

const ConstraintsParams& Parameters::get_constraints_params() const {
    return constraints_params;
}

const KernelParams& Parameters::get_kernel_params() const {
    return kernel_params;
}
