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
    if (numeric.beta <= 0) {
        errors.emplace_back("numeric.beta", "Learning decay beta must be positive");
    }
    
    // Validate OU parameters
    if (ou.sigma < 0) {
        errors.emplace_back("ou.sigma", "Price volatility sigma must be non-negative");
    }
    if (ou.S0 < 0) {
        errors.emplace_back("ou.S0", "Initial price S0 must be non-negative");
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
    // Use the new update_price_parameters method from Parameters class
    dst.update_price_parameters(src);
}

std::vector<ParamError> Interface::read_par(const std::string& file_path) {
    try {
        Parameters temp_params(file_path);
        auto errors = validate_all(temp_params);
        
        if (errors.empty()) {
            // No errors - freeze the parameters and create the solver
            user_params = std::move(temp_params);
            solver = std::make_unique<StochasticUzawaSolver>(user_params.value());
        }
        
        return errors;
    } catch (const std::exception& e) {
        return {ParamError("file", std::string("Failed to read parameter file: ") + e.what())};
    }
}

std::vector<ParamError> Interface::update_price_par(const std::string& file_path) {
    if (!user_params.has_value() || !solver) {
        return {ParamError("state", "Must call read_par() successfully before update_price_par()")};
    }
    
    try {
        Parameters temp_params(file_path);
        auto errors = validate_price_only(temp_params);
        
        if (errors.empty()) {
            // Merge only price-related parameters
            user_params->update_price_parameters(temp_params);
            // Update the solver with new parameters instead of recreating it
            solver->update_parameters(user_params.value());
        }
        
        return errors;
    } catch (const std::exception& e) {
        return {ParamError("file", std::string("Failed to read parameter file: ") + e.what())};
    }
}

Eigen::MatrixXd Interface::simulate_price() {
    if (!user_params.has_value() || !solver) {
        throw std::logic_error("Must call read_par() successfully before simulate_price()");
    }
    
    return solver->simulate_signal();
}

void Interface::solve() {
    if (!user_params.has_value() || !solver) {
        throw std::logic_error("Must call read_par() successfully before start_solver()");
    }
    
    solver->solve();
}

Eigen::MatrixXd Interface::get_price() {
    if (!solver) {
        throw std::logic_error("Must call start_solver() and simulate_signal() before accessing price");
    }
    return solver->get_price();}

const Eigen::MatrixXd& Interface::get_u() {
    if (!solver) {
        throw std::logic_error("Must call start_solver() and solve() before accessing u");
    }
    return solver->get_u();
}

const Eigen::MatrixXd& Interface::get_X() {
    if (!solver) {
        throw std::logic_error("Must call start_solver() and solve() before accessing X");
    }
    return solver->get_X();
}

const Eigen::MatrixXd& Interface::get_lambda1() {
    if (!solver) {
        throw std::logic_error("Must call start_solver() and solve() before accessing lambda1");
    }
    return solver->get_lambda1();
}

const Eigen::MatrixXd& Interface::get_lambda2() {
    if (!solver) {
        throw std::logic_error("Must call start_solver() and solve() before accessing lambda2");
    }
    return solver->get_lambda2();
}

const Eigen::MatrixXd& Interface::get_lambda3() {
    if (!solver) {
        throw std::logic_error("Must call start_solver() and solve() before accessing lambda3");
    }
    return solver->get_lambda3();
}

const Eigen::MatrixXd& Interface::get_lambda4() {
    if (!solver) {
        throw std::logic_error("Must call start_solver() and solve() before accessing lambda4");
    }
    return solver->get_lambda4();
}

const Eigen::VectorXd& Interface::get_time_grid() {
    if (!solver) {
        throw std::logic_error("Must call start_solver() before accessing time_grid");
    }
    return solver->get_time_grid();
}

std::size_t Interface::iterations() {
    if (!solver) {
        throw std::logic_error("Must call start_solver() and solve() before accessing iterations");
    }
    return solver->iterations();
}
