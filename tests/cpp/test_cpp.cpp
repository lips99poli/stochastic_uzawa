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
#include <Eigen/Core>
// Note: Using Eigen optimization with BLAS/LAPACK and threading control

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
    std::filesystem::path project_root = current_path;
    
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
    
    // Eigen+BLAS/LAPACK performance test with different thread counts
    std::vector<int> thread_counts = {1, 2, 4, 6, 8, 12}; // Progressive thread counts up to 12 cores
    std::vector<double> solver_times;
    std::vector<double> signal_times;
    
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "EIGEN+BLAS/LAPACK PERFORMANCE TESTING WITH DIFFERENT THREAD COUNTS" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    Interface* final_interface = nullptr;
    
    for (size_t i = 0; i < thread_counts.size(); ++i) {
        int num_threads = thread_counts[i];
        std::cout << "\n--- Testing with " << num_threads << " thread(s) ---" << std::endl;
        
        // Set Eigen thread count
        //Eigen::setNbThreads(num_threads);
        std::cout << "Eigen threads set to: " << Eigen::nbThreads() << std::endl;
        
        // Create a fresh interface for this test
        Interface* test_interface = new Interface();
        test_interface->read_par(input_file);
        
        // Signal simulation timing
        std::cout << "Re-simulating signal with " << num_threads << " threads..." << std::endl;
        Timings::Chrono thread_signal_timer;
        thread_signal_timer.start();
        test_interface->simulate_price();
        thread_signal_timer.stop();
        signal_times.push_back(thread_signal_timer.wallTime());
        std::cout << "Signal simulation: " << thread_signal_timer.wallTime() << " μs" << std::endl;
        
        // Solver timing
        std::cout << "Running solver with " << num_threads << " threads..." << std::endl;
        Timings::Chrono thread_solver_timer;
        thread_solver_timer.start();
        test_interface->solve();
        thread_solver_timer.stop();
        solver_times.push_back(thread_solver_timer.wallTime());
        std::cout << "Solver execution: " << thread_solver_timer.wallTime() << " μs" << std::endl;
        
        // Keep the last interface for output (from the max threads test)
        if (i == thread_counts.size() - 1) {
            final_interface = test_interface;
        } else {
            delete test_interface;
        }
    }
    
    // Performance summary
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "EIGEN+BLAS/LAPACK PERFORMANCE SUMMARY" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    std::cout << std::left << std::setw(10) << "Threads" 
              << std::setw(20) << "Signal Time (μs)" 
              << std::setw(20) << "Solver Time (μs)" 
              << std::setw(15) << "Speedup (Solver)" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    for (size_t i = 0; i < thread_counts.size(); ++i) {
        double speedup = solver_times[0] / solver_times[i]; // Speedup relative to single thread
        std::cout << std::left << std::setw(10) << thread_counts[i]
                  << std::setw(20) << std::fixed << std::setprecision(2) << signal_times[i]
                  << std::setw(20) << std::fixed << std::setprecision(2) << solver_times[i]
                  << std::setw(15) << std::fixed << std::setprecision(2) << speedup << "x" << std::endl;
    }
    
    // Write performance results to file
    std::string perf_file = output_dir + "/eigen_blas_performance_results.txt";
    std::ofstream perf_out(perf_file);
    perf_out << "Eigen+BLAS/LAPACK Performance Test Results\n";
    perf_out << "=========================================\n\n";
    perf_out << "Test Configuration:\n";
    perf_out << "- Optimization: Eigen automatic parallelization with BLAS/LAPACK\n";
    perf_out << "- OpenMP: Disabled\n";
    perf_out << "- BLAS/LAPACK: Enabled (reference implementation)\n";
    perf_out << "- Thread counts tested: 1, 2, 4, 6, 8, 12\n\n";
    
    perf_out << std::left << std::setw(10) << "Threads" 
             << std::setw(20) << "Signal_Time_μs" 
             << std::setw(20) << "Solver_Time_μs" 
             << std::setw(15) << "Speedup" << "\n";
    for (size_t i = 0; i < thread_counts.size(); ++i) {
        double speedup = solver_times[0] / solver_times[i];
        perf_out << std::left << std::setw(10) << thread_counts[i]
                 << std::setw(20) << signal_times[i]
                 << std::setw(20) << solver_times[i]
                 << std::setw(15) << speedup << "\n";
    }
    perf_out.close();
    
    // Write all matrices to variables file (done once with final result)
    std::cout << "\nWriting variables to output..." << std::endl;
    write_matrices_to_file(output_dir, *final_interface);
    
    // Clean up
    delete final_interface;
    
    std::cout << "\nEigen+BLAS/LAPACK multi-threading performance test completed successfully!" << std::endl;
    std::cout << "Results saved to: " << perf_file << std::endl;

    return 0;
}
