#ifndef NYSTROM_SCHEME_HPP
#define NYSTROM_SCHEME_HPP

#include <Eigen/Dense>
#include <vector>
#include <memory>
#include "Kernel.hpp"

//create aliases for Eigen types
using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using Ref = Eigen::Ref<Matrix>; // Specify the template argument for Eigen::Ref


class NystromScheme {
//private:
public:
    // Parameters
    const std::size_t N; // Number of time steps
    const Matrix I_N;

    const std::size_t M; // Numbeer of Monte Carlo paths
    const Vector& time_grid; // Time grid of length N+1
    const Matrix X0;
    const Kernel& kernel; // Kernel object
    // State Variables
    const std::unique_ptr<Matrix>& u; // u matrix, of size (M, N)
    Matrix& Z_u; // Z_u matrix, of size (M, N)
    const std::unique_ptr<Matrix>& X_u; // X_u matrix, of size (M, N)

    const std::vector<Matrix>& R; // Signal matrix, vector of size M of matrices (N,N)
    const std::vector<Matrix>& Gamma; // Expectation of phi, vector of size M of matrices (N,N)
    // Signal Matrix
    std::vector<Matrix> R_sum_Gamma; // vector of size M ofmatrices (N,N)

    // Operators
    std::vector<Matrix> inv_Dt; // inverse of (I +K_t +K*_t), vector of size N of matrices (N,N)
    Matrix time_integral;
    Matrix L,U; // Kernel matrixes
    Matrix B;
    Matrix I_B_inv; // (I_N - B)^-1

    Matrix a_operator; // Matrix (N,N) where each row i is the operatore to apply at time t_i


public: // constructor
    NystromScheme(const std::size_t M, const std::size_t N, const Vector& time_grid, const double X0, const Kernel& k, const std::unique_ptr<Matrix>& u, Matrix& Z_u,
                  const std::unique_ptr<Matrix>& X_u, const std::vector<Matrix>& R, const std::vector<Matrix>& Gamma)
        : M(M), N(N),I_N(Matrix::Identity(N, N)), time_grid(time_grid), X0(Matrix::Constant(M, N, X0)), kernel(k), u(u), Z_u(Z_u), X_u(X_u), R(R), Gamma(Gamma)
    {
        // Construct Matrixes for operators that do not depend on the control variable u
        construct_time_integral();
        construct_L_and_U(); // Construct the L and U matrices based on the kernel and time grid
        construct_inv_Dt(); // Construct the inverse of (I +K_t +K*_t) for each time step
        construct_B_and_a_operator();
        // Reserve space for the R_sum_Gamma vector of dimension M
        R_sum_Gamma.reserve(M);
    }

//private: // Construct operators
    void construct_time_integral() {
        // Each row of time_integral contains the weights to write: integral_0^t us ds = sum_{j=0}^{N-1} us * (t_{j+1} - t_j)
        Vector time_delta = time_grid.tail(N) - time_grid.head(N);
        time_integral = Matrix::Zero(N, N);
        for (std::size_t i = 0; i < N; ++i) {
            time_integral.row(i).head(i+1) = time_delta.head(i+1);
        }
    }
    // Construct the L and U matrices based on the kernel and time grid
    // L(i,j) = ∫_{t_j}^{t_{j+1}} K(t_i, s) ds
    // U(i,j) = ∫_{t_j}^{t_{j+1}} K(s, t_i) ds
    void construct_L_and_U() {
        L = Matrix::Zero(N, N);
        U = Matrix::Zero(N, N);
        // Construct the L and U matrices based on the kernel and time grid
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = 0; j < i; ++j) { // int_tj^tj+1 of K(t_i, s)
                L(i, j) = kernel.s_integral(time_grid(i), time_grid(j), time_grid(j + 1));
            }
            for (std::size_t j = i; j < N; ++j) { // int_tj^tj+1 of K(s, t_i)
                U(i, j) = kernel.t_integral(time_grid(i), time_grid(j), time_grid(j + 1));
            }
        }
    }
    // Construct the inverse of (I +K_t +K*_t) for each time step
    void construct_inv_Dt() {
        inv_Dt.reserve(N);
        inv_Dt.emplace_back((I_N + L + U).inverse());
        for(std::size_t i = 1; i < N; ++i) {
            // Matrix Kt_i is matrix L where the first i columns and i rows are zero
            Matrix Kt_i = Matrix::Zero(N, N);
            Kt_i.bottomRightCorner(N-i, N-i) = L.bottomRightCorner(N-i, N-i);

            // Matrix K*_t_i is matrix U where the first i columns and i rows are zero
            Matrix K_star_t_i = Matrix::Zero(N, N);
            K_star_t_i.bottomRightCorner(N-i, N-i) = U.bottomRightCorner(N-i, N-i);

            // Compute the inverse of (I + Kt_i + K*_t_i)
            inv_Dt.emplace_back((I_N + Kt_i + K_star_t_i).inverse());
        }
    }
    // Construct operators B and a_operator
    void construct_B_and_a_operator() {
        B = Matrix::Zero(N, N);
        a_operator = Matrix::Zero(N, N);
        for (std::size_t i = 0; i < N; ++i) {
            a_operator.row(i) = U.row(i) * inv_Dt[i];
            for(std::size_t j = 0; j < i; ++j) {
                B(i,j) = a_operator.row(i) * L.col(j);
            }
        }
        // Compute the control variable u based on the Nystrom scheme
        I_B_inv = (I_N - B).inverse();
    }

    // Update the state variables Z_u and X_u based on the control variable u
    void update_state_variables() {
        // Update the transient price impact Z_u based on the current control variable u: Z_u(t) = integral_0^t u(s)*K(t,s)ds
        Z_u.noalias() = (*u) * L.transpose();
        // Update the inventory X_u based on the current control variable u: X_u(t) = X0 + integral_0^t u(s)ds
        X_u->noalias() = X0 + (*u) * time_integral.transpose();
    }

public:
    // Update control variable u and state variables Z_u and X_u
    void nystrom_update(){
        // Update the R_sum_Gamma vector with R+Gamma at each m<M
        R_sum_Gamma.clear();
        for (std::size_t m = 0; m < M; ++m) {
            R_sum_Gamma.emplace_back(R[m] + Gamma[m]);
        }
        // Update the control variable u based on the Nystrom scheme
        for (std::size_t m = 0; m < M; ++m) {
            Vector a(N);
            for (std::size_t i = 0; i < N; ++i) {
                a[i] = R_sum_Gamma[m](i,i) - a_operator.row(i).dot(R_sum_Gamma[m].row(i));
            }
            u->row(m) = I_B_inv * a;
        }
        update_state_variables();
    }
    
    void print_operators(std::ostream& output_stream) {
        output_stream << "=== NYSTROM SCHEME OPERATORS ===" << std::endl;
        
        output_stream << "time_integral matrix:" << std::endl;
        output_stream << time_integral << std::endl << std::endl;
        
        output_stream << "L matrix:" << std::endl;
        output_stream << L << std::endl << std::endl;
        
        output_stream << "U matrix:" << std::endl;
        output_stream << U << std::endl << std::endl;
        
        output_stream << "inv_Dt matrices:" << std::endl;
        for (size_t i = 0; i < inv_Dt.size(); ++i) {
            output_stream << "inv_Dt[" << i << "]:" << std::endl;
            output_stream << inv_Dt[i] << std::endl << std::endl;
        }
        
        output_stream << "B matrix:" << std::endl;
        output_stream << B << std::endl << std::endl;
        
        output_stream << "I_B_inv matrix:" << std::endl;
        output_stream << I_B_inv << std::endl << std::endl;
        
        output_stream << "a_operator matrix:" << std::endl;
        output_stream << a_operator << std::endl << std::endl;
        
        output_stream << "========================" << std::endl << std::endl;
    }
};


#endif // NYSTROM_SCHEME_HPP