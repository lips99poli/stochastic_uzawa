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
        const Vector time_delta; // Time delta grid
        
        const Matrix BM; // Brownian motion paths
        const Matrix OU; // Ornstein-Uhlenbeck process paths
        const Matrix price;
        const Matrix alpha;

        Matrix generateBrownianMotion(const size_t M, int seed = 42) const{
            // Generate M paths of Brownian motion over the time grid
            Matrix BM(M, time_delta.size()); //time_delta.size()=time_grid.size() - 1
            std::default_random_engine generator(seed); // Fixed seed for reproducibility
            std::normal_distribution<double> distribution(0.0, 1.0);

            // Precompute sqrt of time increments
            Vector sqrt_dt = time_delta.array().sqrt();
            
            for (size_t i = 0; i < M; ++i) {
                for (size_t j = 1; j < time_delta.size(); ++j) {
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
                    OU(i, j) = OU(i, j - 1) + (A_t - ou_params.k * OU(i, j - 1)) * time_delta(j-1) + ou_params.psi * (BM(i, j) - BM(i, j - 1));
                }
            }
            return OU;
        }

        Matrix compute_price() const{
            // Compute the price based on the OU process paths
            Matrix price(OU.rows(), OU.cols());
            Matrix BM_price = generateBrownianMotion(OU.rows()); // Generate Brownian motion for pricing
            
            for (size_t m = 0; m < OU.rows(); ++m) {
                for (size_t j = 0; j < OU.cols(); ++j) {
                    price(m, j) = ou_params.S0 + OU.row(m).head(j).dot(time_delta.head(j)) + ou_params.sigma*BM_price(j);
                }
            }
            return price;
        }

        Matrix compute_alpha() const {
            // Compute the alpha matrix based on the OU process paths
            Vector P_T  = OU * time_delta; // Compute the integral of the OU process
            Matrix alpha(OU.rows(), OU.cols());
            for (std::size_t i = 0; i<OU.cols(); ++i){
                alpha.col(i) = P_T - OU.rightCols(ou_params.N-1-i) * time_delta.tail(ou_params.N-1-i);
            }
            return alpha;
        }


    public:
        OUSimulator(const OUParams& ou_params, const Vector& time_grid, const double M):
            ou_params(ou_params)
            ,time_grid(time_grid)
            ,time_delta(time_grid.tail(time_grid.size() - 1) - time_grid.head(time_grid.size() - 1)) // Compute time step size
            ,BM(generateBrownianMotion(M)) // Generate Brownian motion paths
            ,OU(generateOrnsteinUhlenbeck()) // Generate Ornstein-Uhlenbeck process paths
            ,price(compute_price())
            ,alpha(compute_alpha())
        {}
    
        Matrix getOU() const {
            // Return the OU process paths as the signal matrix R
            return OU;
        }
        Matrix getPrice() const {
            // Return the price paths
            return price;
        }
        Matrix getAlpha() const {
            // Return the alpha matrix
            return alpha;
        }

        void print() const {
            // Print BM matrix for debugging purposes
            std::cout << "Brownian motion paths:\n" << BM << std::endl;
            // Print the OU matrix for debugging purposes
            std::cout << "Ornstein-Uhlenbeck process paths:\n" << OU << std::endl;
        }

};

#endif
