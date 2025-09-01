#include <iostream>
#include <fstream>
#include <memory>
#include "../../include/StochasticUzawaSolver.hpp"
#include "../../include/NumericSchemeParams.hpp"
#include "../../include/OUParams.hpp"
#include "../../include/Kernel.hpp"

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <numeric_params_file> <ou_params_file> <experiment_number>" << std::endl;
        std::cerr << "Example: " << argv[0] << " data/data_big.pot data/OU_data.pot 1" << std::endl;
        return 1;
    }
    
    std::string numeric_params_file = argv[1];
    std::string ou_params_file = argv[2];
    std::string experiment_number = argv[3];
    
    try {
        // Load numeric scheme parameters from command line argument
        NumericSchemeParams params = load_params(numeric_params_file);
        
        // Load OU parameters from command line argument
        OUParams ou_params = OU_load_params(params.T, params.N, ou_params_file);
        
        // Create output file stream with experiment number
        std::string output_filename = "verbose_output_" + experiment_number + ".txt";
        std::ofstream output_file(output_filename);
        if (!output_file) {
            std::cerr << "Error: Could not create output file: " << output_filename << std::endl;
            return 1;
        }
        
        // Print initial parameters to output file
        output_file << "=== VERBOSE TEST - NUMERIC SCHEME PARAMETERS ===" << std::endl;
        output_file << "T = " << params.T << std::endl;
        output_file << "N = " << params.N << std::endl;
        output_file << "D = " << params.D << std::endl;
        output_file << "M = " << params.M << std::endl;
        output_file << "epsilon = " << params.epsilon << std::endl;
        output_file << "delta = " << params.delta << std::endl;
        output_file << "beta = " << params.beta << std::endl;
        output_file << std::endl;
        
        output_file << "=== OU PROCESS PARAMETERS ===" << std::endl;
        output_file << "sigma = " << ou_params.sigma << std::endl;
        output_file << "S0 = " << ou_params.S0 << std::endl;
        output_file << "I0 = " << ou_params.I0 << std::endl;
        output_file << "theta = " << ou_params.theta << std::endl;
        output_file << "omega = " << ou_params.omega << std::endl;
        output_file << "phi = " << ou_params.phi << std::endl;
        output_file << "k = " << ou_params.k << std::endl;
        output_file << "psi = " << ou_params.psi << std::endl;
        output_file << std::endl;
        
        // Initial state value
        double X0 = -5.0;
        output_file << "X0 = " << X0 << std::endl;
        output_file << std::endl;
        
        // Create a kernel instance with default parameters
        Kernel kernel(1.0, 1.0); // c=1.0, ro=1.0
        
        output_file << "=== VERBOSE TEST FOR CONSTRAINT OPTIMIZATION ===" << std::endl;
        output_file << std::endl;
        
        // Create and initialize the Stochastic Uzawa Solver
        StochasticUzawaSolver solver(params, kernel, ou_params, &output_file, X0, "verbose");
        
        // Print initialized variables
        output_file << "=== INITIALIZED VARIABLES ===" << std::endl;
        
        output_file << "Initial u matrix:" << std::endl;
        output_file << *(solver.u) << std::endl << std::endl;
        
        output_file << "Initial X_u matrix:" << std::endl;
        output_file << *(solver.X_u) << std::endl << std::endl;
        
        output_file << "Constraints matrices:" << std::endl;
        for (size_t i = 0; i < solver.constraints.size(); ++i) {
            output_file << "Constraint " << i << ":" << std::endl;
            output_file << solver.constraints[i] << std::endl << std::endl;
        }
        
        output_file << "Alpha matrix:" << std::endl;
        output_file << solver.alpha << std::endl << std::endl;
        
        output_file << "R matrices:" << std::endl;
        for (size_t i = 0; i < solver.R.size(); ++i) {
            output_file << "R[" << i << "]:" << std::endl;
            output_file << solver.R[i] << std::endl << std::endl;
        }
        
        output_file << "Starting optimization with verbose output..." << std::endl;
        solver.solve();
        
        // Print final results
        output_file << std::endl << "=== FINAL RESULTS ===" << std::endl;
        
        output_file << "Final u matrix:" << std::endl;
        output_file << *(solver.u) << std::endl << std::endl;
        
        output_file << "Final X_u matrix:" << std::endl;
        output_file << *(solver.X_u) << std::endl << std::endl;
        
        output_file << "Final Z_u matrix:" << std::endl;
        output_file << solver.Z_u << std::endl;
        
        output_file << "Verbose test completed successfully!" << std::endl;
        output_file.close();
        
        std::cout << "Verbose test completed. Results written to " << output_filename << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
