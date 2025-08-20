#ifndef STOCHASTIC_UZAWA_SOLVER_HPP
#define STOCHASTIC_UZAWA_SOLVER_HPP

#include <vector>
#include <memory>
#include <Eigen/Dense>
//create aliases for Eigen types
using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using Ref = Eigen::Ref<Matrix>; // Specify the template argument for Eigen::Ref

using CondExp_ij = std::vector<Matrix>; // Map to store conditional expectations for each time step couple (i,j)


#include "NumericSchemeParams.hpp"
#include "Kernel.hpp"
#include "OUParams.hpp"
#include "OUSimulator.hpp"
#include "LSMCR.hpp"
#include "NystromScheme.hpp"


class StochasticUzawaSolver {
    private:

    double X0 = 1;

    const NumericSchemeParams params; // Reference to the numeric scheme parameters
    const double time_delta; // Time step size
    const std::size_t N; // Number of time steps
    const Vector time_grid; // Time grid

    const Kernel& kernel; // Kernel object

    const OUParams ou_params; // OU process parameters
    
    std::vector<Matrix> R; // Signal matrix
    Matrix alpha;

    std::vector<Matrix> Gamma; // Expectation of phi

    const std::vector<Matrix> constraints; // 4 Constraints matrix
    const std::vector<Matrix> U,L; // Kernel matrixes

    // The output u, X_u and lamda are stored in the heap
    std::unique_ptr<Matrix> u; // Control vector
    std::unique_ptr<Matrix> X_u; // State vector
    Matrix Z_u; // Transient price impact
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
        std::vector<Matrix> constraints(4, Matrix::Zero(M, N)); // 4 constraints, each of size MxN
        for (int m = 0; m < M; ++m) {
            for (int i = 0; i < N; ++i) {
                constraints[0](m, i) = -1; // First constraint: u >= 0
                constraints[1](m, i) = 1; // Second constraint: u <= 0
                constraints[2](m, i) = -100; // Third constraint: X_u >= 0
                constraints[3](m, i) = 100; // Fourth constraint: X_u <= 0
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

    auto compute_slackness() const {
        std::vector<Matrix> slackness_m(4,Matrix::Zero(params.M, N));
        slackness_m[0] = (constraints[0] - *u).cwiseProduct(lambda[0]->leftCols(N));
        slackness_m[1] = (*u - constraints[1]).cwiseProduct(lambda[1]->leftCols(N));
        slackness_m[2] = (constraints[2] - *X_u).cwiseProduct(lambda[2]->leftCols(N));
        slackness_m[3] = (*X_u - constraints[3]).cwiseProduct(lambda[3]->leftCols(N));

        std::vector<Vector> slackness_avg_paths(4, Vector::Zero(N));
        // Each slackness variable is a vector containing the sum of all the rows of the corresponding slackness_m matrix
        for (std::size_t i = 0; i < slackness_m.size(); ++i) {
            slackness_avg_paths[i] = slackness_m[i].colwise().sum() / static_cast<double>(params.M);
        }
        // Integrate the slackness avarage paths
        std::vector<double> slackness (4, 0.0);
        for (std::size_t i = 0; i < slackness_avg_paths.size(); ++i) {
            slackness[i] = slackness_avg_paths[i].sum() * time_delta; // Left rectangle rule
        }
        return std::max_element(slackness.begin(), slackness.end()); // Return the maximum slackness iterator
    }


    public:
    // Note that variables are of time length N, whereas multipliers are of time length N+1
    StochasticUzawaSolver(const NumericSchemeParams& params, const Kernel& k, const OUParams& ou_params) :
        // Numeric scheme parameters
        params(params)
        ,time_delta(params.T / (params.N+1))
        ,N(params.N) // Number of time steps
        ,time_grid(Vector::LinSpaced(params.N+1, 0, params.T)) //see (3.1) paper
        ,kernel(k)
        ,ou_params(ou_params)
        // Control problem parameters
        ,constraints(compute_constraints(params.M))
        ,U(std::vector<Matrix>(1,Matrix::Identity(params.N, params.N))) //risolvere dubbio su B eq Fredholm
        ,L(std::vector<Matrix>(1,Matrix::Identity(params.N, params.N))) //risolvere dubbio su B eq Fredholm
        // Control Variables and Lagrange Multipliers
        ,u( std::make_unique<Matrix>(Matrix::Zero(params.M, params.N)) )
        ,X_u( std::make_unique<Matrix>(Matrix::Zero(params.M, params.N)) )
        ,Z_u( Matrix::Zero(params.M, params.N) )
        ,lambda(std::vector<std::unique_ptr<Matrix>>(4, std::make_unique<Matrix>(Matrix::Zero(params.M, params.N+1)))) // Initialize lambda vector with nullptrs
        ,slackness(std::vector<double>(4, std::numeric_limits<double>::infinity())) // Initialize slackness variables to infinity
        // Cycle tools

    {
        //launch the simulation of the OU process
        OUSimulator ou_simulator(ou_params, time_grid, params.M);
        // Fill the R matrix based on the OU parameters
        R = compute_R(ou_simulator.getOU(), params.M, ou_params);
        // Fill the alpha matrix based on the OU parameters
        alpha = ou_simulator.getAlpha();
        // Initialize gamma to vector of size M of matrices (N,N)
        Gamma.resize(params.M, Matrix::Zero(params.N, params.N));
    }

    // Implement the Uzawa algorithm to solve the control problem
    // This is a placeholder for the actual Uzawa algorithm implementation
    // The algorithm will iteratively update u, X_u, and lambda based on the R, U, and L matrices

    // Cicle steps are the following:
    // 1. Gradient update lambda
    // 2. Estimate conditional expectations through Least Squares Monte Carlo Regression (LSMCR)
    // 3. Update control variable u and state variables X_u, Z_u through Nystrom Scheme and left rectangle rule
    // 4. Update Slackness Variables to check if convergence is reached
    void solve(){
        // Initialize LSMCR structure
        const std::size_t d1 = 3; // Upper bound for the sum of the degrees of the Laguerre polynomials used in the first regression
        const std::size_t d2 = 3; // Upper bound for the sum of the degrees of the Laguerre polynomials used in the second regression
        LSMCR lsmcr(d1, d2, N, params.M, alpha, u, X_u, lambda);

        // Initialize Nystrom Scheme structure
        NystromScheme nystrom_scheme(params.M, params.N, time_grid, X0, kernel, u, Z_u, X_u, R, Gamma);

        // Initialize slackness check variable and check it is correctly initialized
        auto max_slackness_it = std::max_element(slackness.begin(), slackness.end());
        if (max_slackness_it == slackness.end()) {
            throw std::runtime_error("Slackness variables are not initialized correctly.");
        }

        // Main loop for the Uzawa algorithm
        std::size_t n = 0; // Iteration counter
        while (*max_slackness_it > params.epsilon && n < params.D){
            // 1. Gradient update lambda
            do_gradient_update(n);

            // 2. Estimate conditional expectations through Least Squares Monte Carlo Regression (LSMCR)
            lsmcr.update_regressor();
            Gamma = lsmcr.estimate_gamma();
            
            // 3. Update control variable u and state variables X_u, Z_u through Nystrom Scheme and left rectangle rule
            nystrom_scheme.nystrom_update();

            // 4. Update Slackness Variables to check if convergence is reached
            max_slackness_it = compute_slackness();
           
            // Increment the iteration counter
            ++n;
        }
    }
};


#endif // STOCHASTIC_UZAWA_SOLVER_HPP