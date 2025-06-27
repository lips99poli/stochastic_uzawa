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



class StochasticUzawaSolver {
    private:

    const double time_delta; //Time step size
    const size_t N; //Number of time steps
    const Vector time_grid; //Time grid
    
    std::vector<Matrix> R; //Signal matrix
    //std::vector<Matrix> U,L; //Kernel matrixes
    

    // The output u, X_u and lamda are stored in the heap
    std::unique_ptr<Matrix> u; // Control vector
    std::unique_ptr<Matrix> X_u; // State vector
    std::vector<std::unique_ptr<Matrix>> lambda; // Lagrange multipliers

    std::vector<Matrix> fillSignalMatrixes(const OUParams& ou_params, const NumericSchemeParams& params) {
        // Fill the R matrix based on the OU parameters
        OUSimulator ou_simulator(ou_params, time_grid, params.M);
        ou_simulator.print();
        return compute_R(ou_simulator.getSignalMatrix(), params.M, ou_params);
    }

    std::vector<Matrix> compute_R(const Matrix& OU, const int M, const OUParams& ou_params) {
        // Compute the signal matrix R based on the OU process
        std::vector<Matrix> R;
        R.reserve(M); // Reserve space for M matrices
        for (size_t m = 0; m < M; ++m) {
            Matrix R_m = Matrix::Zero(N, N);
            for (size_t i = 0; i < N; ++i) {
                for (size_t j = i; j < N; ++j) {
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


    public:
    StochasticUzawaSolver(const NumericSchemeParams& params, const Kernel& k, const OUParams& ou_params) : 
        time_delta(params.T / (params.N+1))
        ,N(params.N) // Number of time steps
        ,time_grid(Vector::LinSpaced(params.N+1, 0, params.T)) //see (3.1) paper
        ,R(fillSignalMatrixes(ou_params, params))
        // ,U(Matrix::Zero(params.N, params.N))
        // ,L(Matrix::Zero(params.N, params.N))
        
        ,u( std::make_unique<Matrix>(Matrix::Zero(params.M, params.N)) )
        ,X_u( std::make_unique<Matrix>(Matrix::Zero(params.M, params.N)) )
        ,lambda(4)
    {
        // Initialize the lambda vector with unique_ptrs
        for (auto& l : lambda) {
            l = std::make_unique<Matrix>(Matrix::Zero(params.M, params.N));
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