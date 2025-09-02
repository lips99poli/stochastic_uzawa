#ifndef OU_SIMULATOR_HPP
#define OU_SIMULATOR_HPP

#include "Parameters.hpp"

#include "chrono.hpp" 
#include <random>

class OUSimulator {
    // This class will handle the simulation of the Ornstein-Uhlenbeck process
    // and provide methods to generate the signal matrix R based on OU parameters.

    public:
        const OUParams& ou_params; // Parameters for the OU process
        const Vector& time_grid; // Time grid for the simulation
        const Vector time_delta; // Time delta grid
        Matrix time_integral;
        
        // Store the actual seeds used for reproducibility (must be before matrices that use them)
        int actual_seed1;
        int actual_seed2;
        
        const Matrix OU; // Ornstein-Uhlenbeck process paths
        const Matrix Pt;
        const Matrix price;
        Matrix alpha;
        MatVec R; // Signal matrices for each path

        /**
         * @brief Generate a truly random seed using std::random_device
         * @return Random seed value
         */
        static int generate_random_seed() {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<int> dis(1, 1000000);
            return dis(gen);
        }

        Matrix construct_time_integral() {
            // Each row of time_integral contains the weights to write: integral_0^t us ds = sum_{j=0}^{N-1} us * (t_{j+1} - t_j)
            Matrix time_int = Matrix::Zero(time_delta.size(), time_delta.size());
            for (std::size_t i = 1; i < time_delta.size(); ++i) {
                time_int.row(i).head(i) = time_delta.head(i);
            }
            return time_int;
        }

        Matrix generateBrownianMotion(const size_t M, int seed) const{
            // Generate M paths of Brownian motion over the time grid
            Matrix BM(M, time_delta.size()); //time_delta.size()=time_grid.size() - 1
            std::default_random_engine generator(seed); // Fixed seed for reproducibility
            std::normal_distribution<double> distribution(0.0, 1.0);

            // Precompute sqrt of time increments
            Vector sqrt_dt = time_delta.array().sqrt();
            
            for (size_t i = 0; i < M; ++i) {
                for (size_t j = 1; j < time_delta.size(); ++j) {
                    double random_value = distribution(generator);
                    BM(i, j) = BM(i, j - 1) + random_value * sqrt_dt(j - 1);
                }
            }
            return BM;
        }
        
        Matrix generateOrnsteinUhlenbeck(const std::size_t M, int seed) const {

            const Matrix BM = generateBrownianMotion(M, seed); // Generate Brownian motion paths
            // Generate the Ornstein-Uhlenbeck process paths based on the Brownian motion
            Matrix OU(BM.rows(), BM.cols());

            // Initialize the first column with the initial condition
            OU.col(0).setConstant(ou_params.I0);

            // Version using SDE
            for (size_t i = 0; i<BM.rows(); ++i) { //fix sample_path
                for(size_t j = 1; j<BM.cols(); ++j){ //fix time instant
                    double A_t = ou_params.theta * std::sin(ou_params.omega * time_grid(j) + ou_params.phi);
                    OU(i, j) = OU(i, j - 1) + (A_t - ou_params.k * OU(i, j - 1)) * time_delta(j-1) + ou_params.psi * (BM(i, j) - BM(i, j - 1));
                }
            }
            return OU;
        }

        Matrix compute_Pt(){
            return OU * time_integral.transpose();
        }

        Matrix compute_price(int seed) const{
            // Compute the price based on the OU process paths
            Matrix initial_price = Matrix::Constant(OU.rows(), OU.cols(), ou_params.S0);
            Matrix BM_price = ou_params.sigma * generateBrownianMotion(OU.rows(),seed); // Generate Brownian motion for pricing
            
            return initial_price + Pt + BM_price;
        }

        // Matrix compute_alpha() const {
        //     // Compute the alpha matrix based on the OU process paths
        //     Vector P_T  = OU * time_delta; // Compute the final integral of the OU process
        //     Matrix alpha(OU.rows(), OU.cols());
        //     for (std::size_t i = 0; i<OU.cols(); ++i){
        //         alpha.col(i) = P_T - Pt.col(i);
        //     }
        //     return alpha;
        // }

        void compute_alpha_R(const std::size_t M) {
            // Compute the signal matrix R based on the OU process
            R.reserve(M); // Reserve space for M matrices
            for (std::size_t m = 0; m < M; ++m) {
                Matrix R_m = Matrix::Zero(ou_params.N, ou_params.N);
                for (std::size_t i = 0; i < ou_params.N; ++i) {
                    for (std::size_t j = i; j < ou_params.N; ++j) {
                        double t_i = time_grid(i);
                        double t_j = time_grid(j);
                        double omega_ti_phi = ou_params.omega * t_i + ou_params.phi;
                        double omega_tj_phi = ou_params.omega * t_j + ou_params.phi;
                        double omega_T_phi = ou_params.omega * ou_params.T + ou_params.phi;
                        double theta_over_den = ou_params.theta / (ou_params.k * ou_params.k + ou_params.omega * ou_params.omega);
                        R_m(i,j) = ( OU(m,i) - theta_over_den*(ou_params.k*std::sin(omega_ti_phi) - ou_params.omega*std::cos(omega_ti_phi)) ) 
                                * ( exp(-ou_params.k*(t_j-t_i)) - exp(-ou_params.k*(ou_params.T-t_i)) )/ou_params.k;
                        if (ou_params.omega != 0){
                            R_m(i,j) -= theta_over_den * ( 
                                ou_params.k/ou_params.omega * (std::cos(omega_T_phi) - std::cos(omega_tj_phi)) + 
                                std::sin(omega_T_phi) - std::sin(omega_tj_phi) );
                        }
                    }
                }
                alpha.row(m) = R_m.diagonal();
                R.push_back(R_m);
            }
        }

    public:
        /**
         * @brief Constructor for OUSimulator
         * @param ou_params Parameters for the Ornstein-Uhlenbeck process
         * @param time_grid Time grid for the simulation
         * @param M Number of Monte Carlo paths
         * @param seed1 Seed for OU process generation (default: -1 for random)
         * @param seed2 Seed for price process generation (default: -1 for random)
         * 
         * Note: When seed1 or seed2 is -1, a random seed based on current time is generated.
         * This ensures different results for each run unless explicit seeds are provided.
         */
        OUSimulator(const OUParams& ou_params, const Vector& time_grid, const std::size_t M):
            ou_params(ou_params)
            ,time_grid(time_grid)
            ,time_delta(time_grid.tail(time_grid.size() - 1) - time_grid.head(time_grid.size() - 1)) // Compute time step size
            ,time_integral(construct_time_integral())
            ,actual_seed1(ou_params.seed1 == -1 ? generate_random_seed() : ou_params.seed1)
            ,actual_seed2(ou_params.seed2 == -1 ? generate_random_seed() : ou_params.seed2)
            ,OU(generateOrnsteinUhlenbeck(M, actual_seed1)) // Generate Ornstein-Uhlenbeck process paths
            ,Pt(compute_Pt())
            ,price(compute_price(actual_seed2))
            ,alpha(Matrix::Zero(M, ou_params.N)) // Initialize alpha to zero matrix
            ,R()
        {
            compute_alpha_R(M);
        }
    

        Matrix getPrice() const {
            // Return the price paths
            return price;
        }
        Matrix getAlpha() const {
            // Return the alpha matrix
            return alpha;
        }
        MatVec getR() const {
            // Return the signal matrices R
            return R;
        }
};

#endif
