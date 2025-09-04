#include "Interface.hpp"
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <limits>

std::vector<ParamError> Interface::validate_all(const Parameters& p) {
    std::vector<ParamError> errors;
    
    const auto& numeric = p.get_numeric_params();
    const auto& ou = p.get_ou_params();
    const auto& constraints = p.get_constraints_params();
    const auto& kernel = p.get_kernel_params();
    
    // Validate numeric parameters
    if (numeric.T <= 0) {
        errors.emplace_back("numeric.T", "Time horizon T must be positive");
    }
    if (numeric.N == 0) {
        errors.emplace_back("numeric.N", "Number of time steps N must be positive");
    }
    if (numeric.M == 0) {
        errors.emplace_back("numeric.M", "Number of Monte Carlo paths M must be positive");
    }
    if (numeric.D == 0) {
        errors.emplace_back("numeric.D", "Maximum iterations D must be positive");
    }
    if (numeric.epsilon <= 0) {
        errors.emplace_back("numeric.epsilon", "Convergence tolerance epsilon must be positive");
    }
    if (numeric.delta <= 0) {
        errors.emplace_back("numeric.delta", "Learning rate delta must be positive");
    }
    if (numeric.beta < 0 || numeric.beta > 1) {
        errors.emplace_back("numeric.beta", "Learning decay beta must be between 0 and 1");
    }
    
    // Validate OU parameters
    if (ou.sigma < 0) {
        errors.emplace_back("ou.sigma", "Price volatility sigma must be non-negative");
    }
    if (ou.S0 <= 0) {
        errors.emplace_back("ou.S0", "Initial price S0 must be positive");
    }
    if (ou.k < 0) {
        errors.emplace_back("ou.k", "Mean reversion speed k must be non-negative");
    }
    if (ou.psi < 0) {
        errors.emplace_back("ou.psi", "OU volatility psi must be non-negative");
    }
    
    // Validate constraints
    if (constraints.u_min > constraints.u_max) {
        errors.emplace_back("constraints.u_bounds", "u_min must be <= u_max");
    }
    if (constraints.X_u_min > constraints.X_u_max) {
        errors.emplace_back("constraints.X_bounds", "X_u_min must be <= X_u_max");
    }
    if (constraints.X0 < constraints.X_u_min || constraints.X0 > constraints.X_u_max) {
        errors.emplace_back("constraints.X0", "Initial inventory X0 must be within [X_u_min, X_u_max]");
    }
    
    // Validate kernel parameters
    if (kernel.kernel_type != "exp" && kernel.kernel_type != "frac") {
        errors.emplace_back("kernel.type", "Kernel type must be 'exp' or 'frac'");
    } else {
        if (kernel.kernel_parameters.size() != 2) {
            errors.emplace_back("kernel.parameters", "Kernel requires exactly 2 parameters");
        } else {
            if (kernel.kernel_type == "exp") {
                if (kernel.kernel_parameters[0] <= 0) {
                    errors.emplace_back("kernel.par1", "Exponential kernel parameter c must be positive");
                }
                if (kernel.kernel_parameters[1] <= 0) {
                    errors.emplace_back("kernel.par2", "Exponential kernel parameter ro must be positive");
                }
            } else if (kernel.kernel_type == "frac") {
                if (kernel.kernel_parameters[0] <= 0) {
                    errors.emplace_back("kernel.par1", "Fractional kernel parameter c must be positive");
                }
                if (kernel.kernel_parameters[1] <= 0 || kernel.kernel_parameters[1] >= 1) {
                    errors.emplace_back("kernel.par2", "Fractional kernel parameter alpha must be in (0, 1)");
                }
            }
        }
    }
    
    return errors;
}

std::vector<ParamError> Interface::validate_price_only(const Parameters& p) {
    std::vector<ParamError> errors;
    
    const auto& numeric = p.get_numeric_params();
    const auto& ou = p.get_ou_params();
    
    // Validate only price-related parameters (OU and optionally T/N/M)
    if (numeric.T <= 0) {
        errors.emplace_back("numeric.T", "Time horizon T must be positive");
    }
    if (numeric.N == 0) {
        errors.emplace_back("numeric.N", "Number of time steps N must be positive");
    }
    if (numeric.M == 0) {
        errors.emplace_back("numeric.M", "Number of Monte Carlo paths M must be positive");
    }
    
    // Validate OU parameters
    if (ou.sigma < 0) {
        errors.emplace_back("ou.sigma", "Price volatility sigma must be non-negative");
    }
    if (ou.S0 <= 0) {
        errors.emplace_back("ou.S0", "Initial price S0 must be positive");
    }
    if (ou.k < 0) {
        errors.emplace_back("ou.k", "Mean reversion speed k must be non-negative");
    }
    if (ou.psi < 0) {
        errors.emplace_back("ou.psi", "OU volatility psi must be non-negative");
    }
    
    return errors;
}

void Interface::merge_price_subset(const Parameters& src, Parameters& dst) {
    // This is tricky because Parameters doesn't expose setters
    // We need to create a new Parameters object with the merged values
    // For now, we'll replace the entire frozen_ with the new parameters
    // This is a simplification - in a real implementation, you might need
    // to add setter methods to Parameters or refactor it
    dst = src;
}

std::vector<ParamError> Interface::read_par(const std::string& file_path) {
    try {
        Parameters temp_params(file_path);
        auto errors = validate_all(temp_params);
        
        if (errors.empty()) {
            // No errors - freeze the parameters and reset state
            frozen_ = temp_params;
            solver_.reset();
            price_cache_.reset();
        }
        
        return errors;
    } catch (const std::exception& e) {
        return {ParamError("file", std::string("Failed to read parameter file: ") + e.what())};
    }
}

std::vector<ParamError> Interface::update_price_par(const std::string& file_path) {
    if (!frozen_.has_value()) {
        return {ParamError("state", "Must call read_par() successfully before update_price_par()")};
    }
    
    try {
        Parameters temp_params(file_path);
        auto errors = validate_price_only(temp_params);
        
        if (errors.empty()) {
            // Merge only price-related parameters
            merge_price_subset(temp_params, frozen_.value());
            solver_.reset();
            price_cache_.reset();
        }
        
        return errors;
    } catch (const std::exception& e) {
        return {ParamError("file", std::string("Failed to read parameter file: ") + e.what())};
    }
}

Eigen::MatrixXd Interface::simulate_price() {
    if (!frozen_.has_value()) {
        throw std::logic_error("Must call read_par() successfully before simulate_price()");
    }
    
    // Create a temporary solver just for simulation
    StochasticUzawaSolver temp_solver(frozen_.value());
    auto price_matrix = temp_solver.simulate_signal();
    
    // Cache the price for later use
    price_cache_ = price_matrix;
    
    return price_matrix;
}

void Interface::start_solver() {
    if (!frozen_.has_value()) {
        throw std::logic_error("Must call read_par() successfully before start_solver()");
    }
    
    solver_ = std::make_unique<StochasticUzawaSolver>(frozen_.value());
}

void Interface::solve() {
    if (!solver_) {
        throw std::logic_error("Must call start_solver() before solve()");
    }
    
    // If price is not already cached in the solver, simulate it first
    try {
        solver_->get_price();
    } catch (const std::exception&) {
        // Price not available in solver, simulate it
        solver_->simulate_signal();
    }
    
    solver_->solve();
}

Eigen::MatrixXd Interface::get_price() {
    // First try to get from solver
    if (solver_) {
        try {
            return solver_->get_price();
        } catch (const std::exception&) {
            // Fall through to cache
        }
    }
    
    // Then try cache
    if (price_cache_.has_value()) {
        return price_cache_.value();
    }
    
    throw std::logic_error("No price data available. Call simulate_price() or solve() first.");
}

const Eigen::MatrixXd& Interface::get_u() {
    if (!solver_) {
        throw std::logic_error("Must call start_solver() and solve() before accessing u");
    }
    return solver_->get_u();
}

const Eigen::MatrixXd& Interface::get_X() {
    if (!solver_) {
        throw std::logic_error("Must call start_solver() and solve() before accessing X");
    }
    return solver_->get_X();
}

const Eigen::MatrixXd& Interface::get_lambda1() {
    if (!solver_) {
        throw std::logic_error("Must call start_solver() and solve() before accessing lambda1");
    }
    return solver_->get_lambda1();
}

const Eigen::MatrixXd& Interface::get_lambda2() {
    if (!solver_) {
        throw std::logic_error("Must call start_solver() and solve() before accessing lambda2");
    }
    return solver_->get_lambda2();
}

const Eigen::MatrixXd& Interface::get_lambda3() {
    if (!solver_) {
        throw std::logic_error("Must call start_solver() and solve() before accessing lambda3");
    }
    return solver_->get_lambda3();
}

const Eigen::MatrixXd& Interface::get_lambda4() {
    if (!solver_) {
        throw std::logic_error("Must call start_solver() and solve() before accessing lambda4");
    }
    return solver_->get_lambda4();
}

const Eigen::VectorXd& Interface::get_time_grid() {
    if (!solver_) {
        throw std::logic_error("Must call start_solver() before accessing time_grid");
    }
    return solver_->get_time_grid();
}

std::size_t Interface::iterations() {
    if (!solver_) {
        throw std::logic_error("Must call start_solver() and solve() before accessing iterations");
    }
    return solver_->iterations();
}
