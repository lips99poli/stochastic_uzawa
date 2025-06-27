#ifndef OU_SIMULATOR_HPP
#define OU_SIMULATOR_HPP

#include <Eigen/Dense>
//create aliases for Eigen types
using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using Ref = Eigen::Ref<Matrix>;

#include "OUParams.hpp"
#include <random>
#include <cmath> // for std::sin and std::sqrt

class OUSimulator {
    // This class will handle the simulation of the Ornstein-Uhlenbeck process
    // and provide methods to generate the signal matrix R based on OU parameters.

    private:
        OUParams ou_params; // Parameters for the OU process
        const Vector& time_grid; // Time grid for the simulation
        Matrix BM; // Brownian motion paths
        Matrix OU; // Ornstein-Uhlenbeck process paths

        Matrix generateBrownianMotion(const size_t M) const{
            // Generate M paths of Brownian motion over the time grid
            Matrix BM(M, time_grid.size());
            int seed = 42; // Fixed seed for reproducibility PASSALO COME INPUT SUCCESSIVAMENTE, leggerei dai parameteri
            std::default_random_engine generator(seed); // Fixed seed for reproducibility
            std::normal_distribution<double> distribution(0.0, 1.0);

            // Precompute sqrt of time increments
            Vector dt = time_grid.tail(time_grid.size() - 1) - time_grid.head(time_grid.size() - 1);
            Vector sqrt_dt = dt.array().sqrt();
            
            for (size_t i = 0; i < M; ++i) {
                for (size_t j = 1; j < time_grid.size(); ++j) {
                    BM(i, j) = BM(i, j - 1) + distribution(generator) * sqrt_dt(j - 1);
                }
            }
            return BM;
        }
        Matrix generateOrnsteinUhlenbeck() const {
            // Generate the Ornstein-Uhlenbeck process paths based on the Brownian motion
            Matrix OU(BM.rows(), BM.cols());

            // Initialize the first column with the initial condition
            for (size_t i = 0; i < BM.rows(); ++i) {
                OU(i, 0) = ou_params.theta * std::sin(ou_params.omega * time_grid(0) + ou_params.phi); // Initial condition
            }
            // Version using SDE
            for (size_t i = 0; i < BM.rows(); ++i) { //fix sample_path
                for(size_t j = 1; j<BM.cols(); ++j){ //fix time instant
                    double A_t = ou_params.theta * std::sin(ou_params.omega * time_grid(j) + ou_params.phi);
                    double dt = time_grid(j) - time_grid(j - 1);
                    OU(i, j) = OU(i, j - 1) + (A_t - ou_params.k * OU(i, j - 1)) * dt + ou_params.psi * (BM(i, j) - BM(i, j - 1));
                }
            }
            return OU;
        }


    public:
        OUSimulator(const OUParams& ou_params, const Vector& time_grid, const double M):
            ou_params(ou_params)
            ,time_grid(time_grid)
            ,BM(generateBrownianMotion(M)) // Generate Brownian motion paths
            ,OU(generateOrnsteinUhlenbeck()) // Generate Ornstein-Uhlenbeck process paths
        {}
    
        Matrix getSignalMatrix() const {
            // Return the OU process paths as the signal matrix R
            return OU;
        }

        void print() const {
            // Print BM matrix for debugging purposes
            std::cout << "Brownian motion paths:\n" << BM << std::endl;
            // Print the OU matrix for debugging purposes
            std::cout << "Ornstein-Uhlenbeck process paths:\n" << OU << std::endl;
        }

};

#endif
