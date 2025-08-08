#ifndef STOCHASTIC_UZAWA_SOLVER_HPP
#define STOCHASTIC_UZAWA_SOLVER_HPP

#include <vector>
#include <memory>
#include <Eigen/Dense>
//create aliases for Eigen types
using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using Ref = Eigen::Ref<Matrix>; // Specify the template argument for Eigen::Ref

#include "NumericSchemeParams.hpp"
#include "Kernel.hpp"
#include "OUParams.hpp"
#include "OUSimulator.hpp"
#include "LSMCR.hpp"


class StochasticUzawaSolver {
    private:

    const NumericSchemeParams& params; // Reference to the numeric scheme parameters
    const double time_delta; // Time step size
    const std::size_t N; // Number of time steps
    const Vector time_grid; // Time grid
    
    std::vector<Matrix> R; // Signal matrix
    Matrix alpha;

    const std::vector<Matrix> constraints; // 4 Constraints matrix
    const std::vector<Matrix> U,L; // Kernel matrixes

    // The output u, X_u and lamda are stored in the heap
    std::unique_ptr<Matrix> u; // Control vector
    std::unique_ptr<Matrix> X_u; // State vector
    std::vector<std::unique_ptr<Matrix>> lambda; // 4 Lagrange multipliers
    std::vector<double> slackness; // Slackness variables for convergence check

    // Cycle tools
    std::function<double(double,double)> gradient_update = [](double d1, double d2) { return d1 - d2; }; // Gradient update function

    std::vector<Matrix> compute_R(const Matrix& OU, const int M, const OUParams& ou_params) {
        // Compute the signal matrix R based on the OU process
        std::vector<Matrix> R;
        R.reserve(M); // Reserve space for M matrices
        for (std::size_t m = 0; m < M; ++m) {
            Matrix R_m = Matrix::Zero(N, N);
            for (std::size_t i = 0; i < N; ++i) {
                for (std::size_t j = i; j < N; ++j) {
                    double t_i = time_grid(i);
                    double t_j = time_grid(j);
                    double omega_ti_phi = ou_params.omega * t_i + ou_params.phi;
                    double omega_tj_phi = ou_params.omega * t_j + ou_params.phi;
                    double omega_T_phi = ou_params.omega * ou_params.T + ou_params.phi;
                    double theta_over_den = ou_params.theta / (ou_params.k * ou_params.k + ou_params.omega * ou_params.omega);
                    double termine_dubbio =(ou_params.omega == 0) ? 0 : ( ou_params.k/ou_params.omega * (std::cos(omega_T_phi) - std::cos(omega_tj_phi)) + std::sin(omega_T_phi) - std::sin(omega_tj_phi) );
                    R_m(i,j) = ( OU(m,i) - theta_over_den*(ou_params.k*std::sin(omega_ti_phi) - ou_params.omega*std::cos(omega_ti_phi)) ) 
                            * ( exp(-ou_params.k*(t_j-t_i)) - exp(-ou_params.k*(ou_params.T-t_i)) )/ou_params.k 
                            - theta_over_den * termine_dubbio;
                }
            }
            R.push_back(R_m);
        }
        return R;
    }
    std::vector<Matrix> compute_constraints(const int M) { //questo èda estendere: i constraint possono essere di vario tipo: vorrei implementare i 4 del paper e poi lasciare lo user libero di dare una csustom function che calcola ai_t come funzione del prezzo nei tempi da 0 a t
        // Compute the constraints matrix: each of the M rows corresponds to the constraints path for the same row in the price simulation matrix
        std::vector<Matrix> constraints(4, Matrix::Zero(M, M)); // 4 constraints, each of size MxM
        for (int m = 0; m < M; ++m) {
            for (int j = 0; j < M; ++j) {
                constraints[0](m, j) = -1; // First constraint: u >= 0
                constraints[1](m, j) = 1; // Second constraint: u <= 0
                constraints[2](m, j) = -100; // Third constraint: X_u >= 0
                constraints[3](m, j) = 100; // Fourth constraint: X_u <= 0
            }
        }
        return constraints;
    }

    void do_gradient_update(std::size_t n) {
        const double adaptive_learning = params.delta / std::pow(n, params.beta);
        // Perform the gradient update for the Lagrange multipliers
        *lambda[0] += adaptive_learning * constraints[0].binaryExpr(*u,gradient_update);
        *lambda[1] += adaptive_learning * (*u).binaryExpr(constraints[1],gradient_update);
        *lambda[2] += adaptive_learning * constraints[2].binaryExpr(*X_u,gradient_update);
        *lambda[3] += adaptive_learning * (*X_u).binaryExpr(constraints[3],gradient_update);
    }


    public:
    // Note that variables are of time length N, whereas multipliers are of time length N+1
    StochasticUzawaSolver(const NumericSchemeParams& params, const Kernel& k, const OUParams& ou_params) :
        // Numeric scheme parameters
        params(params)
        ,time_delta(params.T / (params.N+1))
        ,N(params.N) // Number of time steps
        ,time_grid(Vector::LinSpaced(params.N+1, 0, params.T)) //see (3.1) paper
        // Control problem parameters
        ,constraints(compute_constraints(params.M))
        ,U(std::vector<Matrix>(1,Matrix::Identity(params.N, params.N))) //risolvere dubbio su B eq Fredholm
        ,L(std::vector<Matrix>(1,Matrix::Identity(params.N, params.N))) //risolvere dubbio su B eq Fredholm
        // Control Variables and Lagrange Multipliers
        ,u( std::make_unique<Matrix>(Matrix::Zero(params.M, params.N)) )
        ,X_u( std::make_unique<Matrix>(Matrix::Zero(params.M, params.N)) )
        ,lambda(std::vector<std::unique_ptr<Matrix>>(4, std::make_unique<Matrix>(Matrix::Zero(params.M, params.N+1)))) // Initialize lambda vector with nullptrs
        ,slackness(std::vector<double>(4, std::numeric_limits<double>::infinity())) // Initialize slackness variables to infinity
        // Cycle tools
    {
        //launch the simulation of the OU process
        OUSimulator ou_simulator(ou_params, time_grid, params.M);
        ou_simulator.print();
        // Fill the R matrix based on the OU parameters
        R = compute_R(ou_simulator.getOU(), params.M, ou_params);
        // Fill the alpha matrix based on the OU parameters
        alpha = ou_simulator.getAlpha();
    }

    void solve(){
        // Implement the Uzawa algorithm to solve the control problem
        // This is a placeholder for the actual Uzawa algorithm implementation
        // The algorithm will iteratively update u, X_u, and lambda based on the R, U, and L matrices

        // Cicle steps are the following:
        // 1. Gradient update lambda
        // 2. Estimate conditional expectations through Least Squares Monte Carlo Regression (LSMCR)
        // 3. Update control variable u through Nystrom Scheme
        // 4. Update state variable X_u integrating the control variable u with left rectangle rule
        // 5. Update Slackness Variables to check if convergence is reached

        // Initilize LSMCR structure
        const std::size_t d1 = 3; // Upper bound for the sum of the degrees of the Laguerre polynomials used in the first regression
        const std::size_t d2 = 3; // Upper bound for the sum of the degrees of the Laguerre polynomials used in the second regression
        LSMCR lsmcr(d1,d2);

        auto max_slackness_it = std::max_element(slackness.begin(), slackness.end());
        //check it is not equal to slackness.end()
        if (max_slackness_it == slackness.end()) {
            throw std::runtime_error("Slackness variables are not initialized correctly.");
        }

        // Main loop for the Uzawa algorithm
        std::size_t n = 0; // Iteration counter
        while (*max_slackness_it > params.epsilon && n < params.D){
            // 1. Gradient update lambda
            do_gradient_update(n);

            // 2. Estimate conditional expectations through Least Squares Monte Carlo Regression (LSMCR)
            
            

            max_slackness_it = std::max_element(slackness.begin(), slackness.end());
            
        }




    }




    void verify(){
        // print time grid
        std::cout << "Time grid: " << time_grid.transpose() << std::endl;
        // Verify the dimensions of the matrices
        if (R.empty()) {
            throw std::runtime_error("Signal matrix R is empty.");
        }
        if (u->rows() != lambda[0]->rows() || u->cols() != N) {
            throw std::runtime_error("Control vector u has incorrect dimensions.");
        }
        if (X_u->rows() != lambda[0]->rows() || X_u->cols() != N) {
            throw std::runtime_error("State vector X_u has incorrect dimensions.");
        }

        //verify R matrix is full not zero
        for (const auto& R_m : R) {
            if (R_m.isZero()) {
                throw std::runtime_error("Signal matrix R contains a zero matrix.");
            }
        }
        // Print matrixes R
        std::cout << "Signal matrix R (size): " << R.size() << " matrices of size " << R[0].rows() << "x" << R[0].cols() << std::endl;
        std::cout << "Signal matrix R (first matrix):\n" << R[0] << std::endl;
        std::cout << "Signal matrix R (first matrix):\n" << R[1] << std::endl;
    }


};


#endif // STOCHASTIC_UZAWA_SOLVER_HPP