#ifndef INTERFACE_HPP
#define INTERFACE_HPP

#include "Parameters.hpp"
#include "StochasticUzawaSolver.hpp"
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <memory>
#include <optional>

struct ParamError {
    std::string path;
    std::string message;
    
    ParamError(const std::string& p, const std::string& m) : path(p), message(m) {}
};

class Interface {
private:
    std::optional<Parameters> user_params;
    std::unique_ptr<StochasticUzawaSolver> solver;
    
    // Internal validation helpers
    std::vector<ParamError> validate_all(const Parameters& p);
    std::vector<ParamError> validate_price_only(const Parameters& p);
    void merge_price_subset(const Parameters& src, Parameters& dst);

public:
    Interface() = default;
    ~Interface() = default;
    
    // Parameter management
    std::vector<ParamError> read_par(const std::string& file_path);
    std::vector<ParamError> update_price_par(const std::string& file_path);
    
    // Simulation and solving
    Eigen::MatrixXd simulate_price();
    void solve();
    
    // Getters
    const Eigen::MatrixXd& get_price() const;
    const Eigen::MatrixXd& get_u() const;
    const Eigen::MatrixXd& get_X() const;
    const Eigen::MatrixXd& get_lambda1() const;
    const Eigen::MatrixXd& get_lambda2() const;
    const Eigen::MatrixXd& get_lambda3() const;
    const Eigen::MatrixXd& get_lambda4() const;
    const Eigen::VectorXd& get_time_grid() const;
    const std::size_t iterations() const;
};

#endif // INTERFACE_HPP
