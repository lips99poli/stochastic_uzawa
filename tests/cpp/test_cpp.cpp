#include "Interface.hpp"
#include "chrono.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>
#include <cstring>
#include <omp.h>
#include <iomanip>
#include <cmath>

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
    
    // Simulate signal (price) - done once since it's the same for all thread tests
    std::cout << "Simulating signal..." << std::endl;
    Timings::Chrono signal_timer;
    signal_timer.start();
    Eigen::MatrixXd price_matrix = interface.simulate_price();
    signal_timer.stop();
    std::cout << "Signal simulation completed in: " << signal_timer.wallTime() << " microseconds. Generated price matrix of size: " 
                << price_matrix.rows() << " x " << price_matrix.cols() << std::endl;
    
    // Thread count tests - Test with different number of threads (3 runs each for consistency)
    std::vector<int> thread_counts = {1, 2, 4, 6, 8, 12}; // Progressive thread counts up to 12 cores
    const int num_runs_per_config = 3; // Number of runs per thread configuration for averaging
    std::vector<std::vector<double>> solver_times_all(thread_counts.size());
    std::vector<std::vector<double>> signal_times_all(thread_counts.size());
    std::vector<double> solver_times_avg(thread_counts.size());
    std::vector<double> signal_times_avg(thread_counts.size());
    Interface* final_interface = nullptr; // Pointer to the final interface for output
    
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "OPENMP PERFORMANCE TESTING WITH DIFFERENT THREAD COUNTS" << std::endl;
    std::cout << "Running " << num_runs_per_config << " iterations per thread configuration for consistency" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    
    for (size_t i = 0; i < thread_counts.size(); ++i) {
        int num_threads = thread_counts[i];
        std::cout << "\n--- Testing with " << num_threads << " thread(s) (" << num_runs_per_config << " runs) ---" << std::endl;
        
        // Set the number of OpenMP threads
        omp_set_num_threads(num_threads);
        std::cout << "OpenMP threads set to: " << omp_get_max_threads() << std::endl;
        
        // Initialize vectors for this thread configuration
        solver_times_all[i].resize(num_runs_per_config);
        signal_times_all[i].resize(num_runs_per_config);
        
        // Run multiple iterations for this thread count
        for (int run = 0; run < num_runs_per_config; ++run) {
            std::cout << "  Run " << (run + 1) << "/" << num_runs_per_config << ":" << std::endl;
            
            // Create a fresh interface for this test (to reset internal state)
            Interface* thread_interface = new Interface();
            thread_interface->read_par(input_file); // Reload parameters
            
            // Re-simulate signal with this thread count
            std::cout << "    Simulating signal..." << std::endl;
            Timings::Chrono thread_signal_timer;
            thread_signal_timer.start();
            thread_interface->simulate_price();
            thread_signal_timer.stop();
            signal_times_all[i][run] = thread_signal_timer.wallTime();
            std::cout << "    Signal simulation: " << thread_signal_timer.wallTime() << " μs" << std::endl;
            
            // Solve the optimization problem with current thread count
            std::cout << "    Running solver..." << std::endl;
            Timings::Chrono thread_solver_timer;
            thread_solver_timer.start();
            thread_interface->solve();
            thread_solver_timer.stop();
            solver_times_all[i][run] = thread_solver_timer.wallTime();
            std::cout << "    Solver execution: " << thread_solver_timer.wallTime() << " μs" << std::endl;
            
            // Keep the last interface for final output (from the last run of max threads)
            if (i == thread_counts.size() - 1 && run == num_runs_per_config - 1) {
                final_interface = thread_interface;
            } else {
                delete thread_interface; // Clean up intermediate interfaces
            }
        }
        
        // Calculate averages for this thread configuration
        double signal_sum = 0, solver_sum = 0;
        for (int run = 0; run < num_runs_per_config; ++run) {
            signal_sum += signal_times_all[i][run];
            solver_sum += solver_times_all[i][run];
        }
        signal_times_avg[i] = signal_sum / num_runs_per_config;
        solver_times_avg[i] = solver_sum / num_runs_per_config;
        
        std::cout << "  Average for " << num_threads << " threads: Signal=" 
                  << std::fixed << std::setprecision(0) << signal_times_avg[i] 
                  << " μs, Solver=" << solver_times_avg[i] << " μs" << std::endl;
    }
    
    // Performance summary
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "OPENMP PERFORMANCE SUMMARY (AVERAGES OF " << num_runs_per_config << " RUNS)" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    std::cout << std::left << std::setw(10) << "Threads" 
              << std::setw(18) << "Avg Signal (μs)" 
              << std::setw(18) << "Avg Solver (μs)" 
              << std::setw(15) << "Speedup" 
              << std::setw(10) << "Std Dev" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    for (size_t i = 0; i < thread_counts.size(); ++i) {
        double speedup = solver_times_avg[0] / solver_times_avg[i]; // Speedup relative to single thread
        
        // Calculate standard deviation for solver times
        double solver_variance = 0;
        for (int run = 0; run < num_runs_per_config; ++run) {
            double diff = solver_times_all[i][run] - solver_times_avg[i];
            solver_variance += diff * diff;
        }
        double solver_std_dev_pct = sqrt(solver_variance / num_runs_per_config) / solver_times_avg[i] * 100;
        
        std::cout << std::left << std::setw(10) << thread_counts[i]
                  << std::setw(18) << std::fixed << std::setprecision(0) << signal_times_avg[i]
                  << std::setw(18) << std::fixed << std::setprecision(0) << solver_times_avg[i]
                  << std::setw(15) << std::fixed << std::setprecision(2) << speedup << "x"
                  << std::setw(10) << std::fixed << std::setprecision(1) << solver_std_dev_pct << "%" << std::endl;
    }
    
    // Write performance results to file
    std::string perf_file = output_dir + "/openmp_performance_results.txt";
    std::ofstream perf_out(perf_file);
    perf_out << "OpenMP Performance Test Results\n";
    perf_out << "===============================\n\n";
    perf_out << "Test Configuration:\n";
    perf_out << "- Optimization: Manual OpenMP parallelization\n";
    perf_out << "- Eigen: No automatic parallelization (EIGEN_DONT_PARALLELIZE)\n";
    perf_out << "- Thread counts tested: 1, 2, 4, 6, 8, 12\n";
    perf_out << "- Runs per configuration: " << num_runs_per_config << " (for statistical consistency)\n";
    perf_out << "- Problem size: Price matrix " << price_matrix.rows() << "x" << price_matrix.cols() << "\n\n";
    
    perf_out << "Average Performance Results:\n";
    perf_out << std::left << std::setw(10) << "Threads" 
             << std::setw(18) << "Avg_Signal_μs" 
             << std::setw(18) << "Avg_Solver_μs" 
             << std::setw(15) << "Speedup" 
             << std::setw(12) << "StdDev_%" << "\n";
    for (size_t i = 0; i < thread_counts.size(); ++i) {
        double speedup = solver_times_avg[0] / solver_times_avg[i];
        
        // Calculate standard deviation for solver times
        double solver_variance = 0;
        for (int run = 0; run < num_runs_per_config; ++run) {
            double diff = solver_times_all[i][run] - solver_times_avg[i];
            solver_variance += diff * diff;
        }
        double solver_std_dev_pct = sqrt(solver_variance / num_runs_per_config) / solver_times_avg[i] * 100;
        
        perf_out << std::left << std::setw(10) << thread_counts[i]
                 << std::setw(18) << signal_times_avg[i]
                 << std::setw(18) << solver_times_avg[i]
                 << std::setw(15) << speedup 
                 << std::setw(12) << solver_std_dev_pct << "\n";
    }
    
    perf_out << "\nDetailed Results (all runs):\n";
    for (size_t i = 0; i < thread_counts.size(); ++i) {
        perf_out << "\n" << thread_counts[i] << " threads:\n";
        for (int run = 0; run < num_runs_per_config; ++run) {
            perf_out << "  Run " << (run+1) << ": Signal=" << signal_times_all[i][run] 
                     << " μs, Solver=" << solver_times_all[i][run] << " μs\n";
        }
    }
    perf_out.close();
    
    // Write all matrices to variables file (done once with final result)
    std::cout << "\nWriting variables to output..." << std::endl;
    write_matrices_to_file(output_dir, *final_interface);
    
    // Clean up
    delete final_interface;
    
    std::cout << "\nOpenMP multi-threading performance test completed successfully!" << std::endl;

    return 0;
}
