#ifndef COMMON_TYPES_HPP
#define COMMON_TYPES_HPP

#include <vector>
#include <Eigen/Dense>
#include <cmath>

#include <omp.h>

// Common type aliases used throughout the project
using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using Ref = Eigen::Ref<Matrix>;
using MatVec = std::vector<Matrix>;

using par_for_type = std::ptrdiff_t; // Type for parallel loop indices, compatible with OpenMP

#endif // COMMON_TYPES_HPP
