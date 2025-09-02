#ifndef COMMON_TYPES_HPP
#define COMMON_TYPES_HPP

#include <vector>
#include <Eigen/Dense>
#include <cmath>

// Common type aliases used throughout the project
using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using Ref = Eigen::Ref<Matrix>;
using MatVec = std::vector<Matrix>;

#endif // COMMON_TYPES_HPP
