#ifndef OU_SIMULATOR_HPP
#define OU_SIMULATOR_HPP

#include "Parameters.hpp"
#include "CommonTypes.hpp"

#include "chrono.hpp" 
#include <random>
#include <cstdint>

class OUSimulator {
    // This class will handle the simulation of the Ornstein-Uhlenbeck process
    // and provide methods to generate the signal matrix R based on OU parameters.

private:
    const OUParams& ou_params; // Parameters for the OU process
    const Vector& time_grid; // Time grid for the simulation
    const Vector time_delta; // Time delta grid
    Matrix time_integral;
    
    // Store the actual seeds used for reproducibility (must be before matrices that use them)
    std::uint64_t actual_seed1;
    std::uint64_t actual_seed2;
    
    const Matrix OU; // Ornstein-Uhlenbeck process paths
    const Matrix Pt;
    const Matrix price;
    Matrix alpha;
    MatVec R; // Signal matrices for each path

    // Private helper methods - implementations in OUSimulator.cpp
    static std::uint64_t make_seed(int user_seed);
    Matrix construct_time_integral();
    Matrix generateBrownianMotion(const size_t M, int seed) const;
    Matrix generateOrnsteinUhlenbeck(const std::size_t M, int seed) const;
    Matrix compute_Pt();
    Matrix compute_price(int seed) const;
    void compute_alpha_R(const std::size_t M);

public:
    /**
     * @brief Constructor for OUSimulator - implementation in OUSimulator.cpp
     * @param ou_params Parameters for the Ornstein-Uhlenbeck process
     * @param time_grid Time grid for the simulation
     * @param M Number of Monte Carlo paths
     * 
     * Note: When seed1 or seed2 is -1, a random seed based on current time is generated.
     * This ensures different results for each run unless explicit seeds are provided.
     */
    OUSimulator(const OUParams& ou_params, const Vector& time_grid, const std::size_t M);

    // Getter methods - implementations in OUSimulator.cpp
    Matrix getPrice() const;
    Matrix getAlpha() const;
    MatVec getR() const;
};

#endif
