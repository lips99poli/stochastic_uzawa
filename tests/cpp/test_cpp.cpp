#include "Interface.hpp"
#include "chrono.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>
#include <cstring>
#include <iomanip>
#include <cmath>
#include <algorithm>
// Note: OpenMP removed - using Eigen optimization instead

void print_errors(const std::vector<ParamError>& errors) {
    std::cerr << "Parameter validation errors:\n";
    for (const auto& error : errors) {
        std::cerr << "  " << error.path << ": " << error.message << "\n";
    }
}

void write_parameters_to_file(const std::string& output_dir, const std::string& input_file) {
    std::filesystem::create_directories(output_dir);
    std::string output_filename = output_dir + "/Parameters.txt";
    
    // Copy the input parameter file to the output directory
    try {
        std::filesystem::copy_file(input_file, output_filename, std::filesystem::copy_options::overwrite_existing);
        std::cout << "Parameters copied from " << input_file << " to: " << output_filename << std::endl;
    } catch (const std::filesystem::filesystem_error& e) {
        throw std::runtime_error("Failed to copy parameter file: " + std::string(e.what()));
    }
}

void write_matrices_to_file(const std::string& output_dir, const Interface& interface) {
    std::filesystem::create_directories(output_dir);
    std::string filename = output_dir + "/variables.txt";
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    try {
        // Write iterations first
        file << "ITERATIONS=" << interface.get_iterations() << std::endl;
        
        // Write time grid
        file << "TIME_GRID=" << std::endl;
        const auto& time_grid = interface.get_time_grid();
        file << time_grid.transpose() << std::endl;
        
        // Write price matrix
        file << "PRICE=" << std::endl;
        const auto& price = interface.get_price();
        file << price << std::endl;
        
        // Write u matrix
        file << "U=" << std::endl;
        const auto& u = interface.get_u();
        file << u << std::endl;
        
        // Write X matrix
        file << "X=" << std::endl;
        const auto& X = interface.get_X();
        file << X << std::endl;
        
        // Write lambda matrices
        file << "LAMBDA1=" << std::endl;
        const auto& lambda1 = interface.get_lambda1();
        file << lambda1 << std::endl;
        
        file << "LAMBDA2=" << std::endl;
        const auto& lambda2 = interface.get_lambda2();
        file << lambda2 << std::endl;
        
        file << "LAMBDA3=" << std::endl;
        const auto& lambda3 = interface.get_lambda3();
        file << lambda3 << std::endl;
        
        file << "LAMBDA4=" << std::endl;
        const auto& lambda4 = interface.get_lambda4();
        file << lambda4 << std::endl;
        
    } catch (const std::exception& e) {
        file << "Error writing matrices: " << e.what() << std::endl;
        throw;
    }
    
    file.close();
    std::cout << "Variables written to: " << filename << std::endl;
}

int main(int argc, char* argv[]) {
    std::string output_folder_name;
    std::string input_file;
    
    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--input") == 0) {
            if (i + 1 < argc) {
                input_file = argv[++i];
            } else {
                std::cerr << "Error: " << argv[i] << " requires a parameter file argument\n";
                return 1;
            }
        } else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) {
            if (i + 1 < argc) {
                output_folder_name = argv[++i];
            } else {
                std::cerr << "Error: " << argv[i] << " requires an output folder name\n";
                return 1;
            }
        } else {
            std::cerr << "Error: Unknown option " << argv[i] << "\n";
            return 1;
        }
    }
    
    // Validate required arguments
    if (output_folder_name.empty() || input_file.empty()) {
        std::cerr << "Error: Both input parameter file and output folder name are required\n";
        std::cerr << "Usage: " << argv[0] << " -i <parameter_file> -o <output_folder_name>\n";
        return 1;
    }
    
    // Get the project root directory (3 levels up from build directory)
    std::filesystem::path current_path = std::filesystem::current_path();
    std::filesystem::path project_root = current_path.parent_path().parent_path().parent_path();
    
    // Create full output path
    std::string output_dir = project_root / "outputs" / "cpp" / output_folder_name;
    
    std::cout << "Project root: " << project_root << std::endl;
    std::cout << "Using parameter file: " << input_file << std::endl;
    std::cout << "Output directory: " << output_dir << std::endl;
    
    // Instantiate the interface
    Interface interface;
    
    // Read and validate parameters
    std::cout << "Reading parameters from: " << input_file << std::endl;
    std::vector<ParamError> errors = interface.read_par(input_file);
    
    // Check for validation errors
    if (!errors.empty()) {
        print_errors(errors);
        std::cerr << "Parameter validation failed. Exiting.\n";
        return 1;
    }
    
    std::cout << "Parameters validated successfully." << std::endl;
    
    // Write parameters to output directory (done once)
    write_parameters_to_file(output_dir, input_file);
    
    // Simulate signal (price) - initial run for verification
    std::cout << "Simulating signal..." << std::endl;
    Timings::Chrono signal_timer;
    signal_timer.start();
    Eigen::MatrixXd price_matrix = interface.simulate_price();
    signal_timer.stop();
    std::cout << "Signal simulation completed in: " << signal_timer.wallTime() << " microseconds. Generated price matrix of size: " 
                << price_matrix.rows() << " x " << price_matrix.cols() << std::endl;
    
    // Eigen optimization performance test - Multiple runs for consistent timing
    const int num_runs = 5; // Number of performance test runs
    std::vector<double> solver_times;
    std::vector<double> signal_times; 
    
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "EIGEN OPTIMIZATION PERFORMANCE TESTING" << std::endl;
    std::cout << "Running " << num_runs << " iterations for consistent timing" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    Interface* final_interface = nullptr;
    
    for (int run = 1; run <= num_runs; ++run) {
        std::cout << "\n--- Performance Run " << run << "/" << num_runs << " ---" << std::endl;
        
        // Create a fresh interface for this test
        Interface* test_interface = new Interface();
        test_interface->read_par(input_file);
        
        // Signal simulation timing
        std::cout << "Signal simulation (run " << run << ")..." << std::endl;
        Timings::Chrono run_signal_timer;
        run_signal_timer.start();
        test_interface->simulate_price();
        run_signal_timer.stop();
        signal_times.push_back(run_signal_timer.wallTime());
        std::cout << "Signal simulation: " << run_signal_timer.wallTime() << " μs" << std::endl;
        
        // Solver timing
        std::cout << "Solver execution (run " << run << ")..." << std::endl;
        Timings::Chrono run_solver_timer;
        run_solver_timer.start();
        test_interface->solve();
        run_solver_timer.stop();
        solver_times.push_back(run_solver_timer.wallTime());
        std::cout << "Solver execution: " << run_solver_timer.wallTime() << " μs" << std::endl;
        
        // Keep the last interface for output
        if (run == num_runs) {
            final_interface = test_interface;
        } else {
            delete test_interface;
        }
    }
    
    // Calculate statistics
    auto calc_stats = [](const std::vector<double>& times) {
        double sum = 0.0, min_val = times[0], max_val = times[0];
        for (double t : times) {
            sum += t;
            min_val = std::min(min_val, t);
            max_val = std::max(max_val, t);
        }
        double avg = sum / times.size();
        double variance = 0.0;
        for (double t : times) variance += (t - avg) * (t - avg);
        double stddev = std::sqrt(variance / times.size());
        return std::make_tuple(avg, min_val, max_val, stddev);
    };
    
    auto [sig_avg, sig_min, sig_max, sig_std] = calc_stats(signal_times);
    auto [sol_avg, sol_min, sol_max, sol_std] = calc_stats(solver_times);
    
    // Performance summary
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "EIGEN OPTIMIZATION PERFORMANCE SUMMARY" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    std::cout << std::left << std::setw(15) << "Metric" 
              << std::setw(18) << "Signal Time (μs)" 
              << std::setw(18) << "Solver Time (μs)" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    std::cout << std::left << std::setw(15) << "Average"
              << std::setw(18) << std::fixed << std::setprecision(2) << sig_avg
              << std::setw(18) << std::fixed << std::setprecision(2) << sol_avg << std::endl;
    std::cout << std::left << std::setw(15) << "Minimum"
              << std::setw(18) << std::fixed << std::setprecision(2) << sig_min
              << std::setw(18) << std::fixed << std::setprecision(2) << sol_min << std::endl;
    std::cout << std::left << std::setw(15) << "Maximum"
              << std::setw(18) << std::fixed << std::setprecision(2) << sig_max
              << std::setw(18) << std::fixed << std::setprecision(2) << sol_max << std::endl;
    std::cout << std::left << std::setw(15) << "Std Deviation"
              << std::setw(18) << std::fixed << std::setprecision(2) << sig_std
              << std::setw(18) << std::fixed << std::setprecision(2) << sol_std << std::endl;
    
    std::cout << "\nIndividual Run Results:" << std::endl;
    std::cout << std::left << std::setw(8) << "Run" 
              << std::setw(18) << "Signal Time (μs)" 
              << std::setw(18) << "Solver Time (μs)" << std::endl;
    std::cout << std::string(50, '-') << std::endl;
    for (int i = 0; i < num_runs; ++i) {
        std::cout << std::left << std::setw(8) << (i + 1)
                  << std::setw(18) << std::fixed << std::setprecision(2) << signal_times[i]
                  << std::setw(18) << std::fixed << std::setprecision(2) << solver_times[i] << std::endl;
    }
    
    // Write performance results to file
    std::string perf_file = output_dir + "/eigen_performance_results.txt";
    std::ofstream perf_out(perf_file);
    perf_out << "Eigen Optimization Performance Test Results\n";
    perf_out << "==========================================\n\n";
    perf_out << "Test Configuration:\n";
    perf_out << "- Optimization: Eigen automatic parallelization\n";
    perf_out << "- OpenMP: Disabled\n";
    perf_out << "- BLAS/LAPACK: " << (true ? "Enabled" : "Disabled") << "\n"; // Will be filled by CMake
    perf_out << "- Number of runs: " << num_runs << "\n\n";
    
    perf_out << "Summary Statistics:\n";
    perf_out << std::left << std::setw(15) << "Metric" 
             << std::setw(18) << "Signal_Time_μs" 
             << std::setw(18) << "Solver_Time_μs" << "\n";
    perf_out << std::left << std::setw(15) << "Average"
             << std::setw(18) << sig_avg
             << std::setw(18) << sol_avg << "\n";
    perf_out << std::left << std::setw(15) << "Minimum"
             << std::setw(18) << sig_min
             << std::setw(18) << sol_min << "\n";
    perf_out << std::left << std::setw(15) << "Maximum"
             << std::setw(18) << sig_max
             << std::setw(18) << sol_max << "\n";
    perf_out << std::left << std::setw(15) << "Std_Deviation"
             << std::setw(18) << sig_std
             << std::setw(18) << sol_std << "\n\n";
    
    perf_out << "Individual Results:\n";
    perf_out << std::left << std::setw(8) << "Run" 
             << std::setw(18) << "Signal_Time_μs" 
             << std::setw(18) << "Solver_Time_μs" << "\n";
    for (int i = 0; i < num_runs; ++i) {
        perf_out << std::left << std::setw(8) << (i + 1)
                 << std::setw(18) << signal_times[i]
                 << std::setw(18) << solver_times[i] << "\n";
    }
    perf_out.close();
    
    // Write all matrices to variables file (done once with final result)
    std::cout << "\nWriting variables to output..." << std::endl;
    write_matrices_to_file(output_dir, *final_interface);
    
    // Clean up
    delete final_interface;
    
    std::cout << "\nEigen optimization performance test completed successfully!" << std::endl;
    std::cout << "Results saved to: " << perf_file << std::endl;

    return 0;
}
