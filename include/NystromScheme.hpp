#ifndef NYSTROM_SCHEME_HPP
#define NYSTROM_SCHEME_HPP

#include <memory>
#include "CommonTypes.hpp"
#include "Kernel.hpp"

class NystromScheme {
private:
    // Parameters
    const std::size_t N; // Number of time steps
    const Matrix I_N;

    const std::size_t M; // Number of Monte Carlo paths
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
    std::vector<Matrix> R_sum_Gamma; // vector of size M of matrices (N,N)

    // Operators
    std::vector<Matrix> inv_Dt; // inverse of (I +K_t +K*_t), vector of size N of matrices (N,N)
    Matrix time_integral;
    Matrix L,U; // Kernel matrices
    Matrix B;
    Matrix I_B_inv; // (I_N - B)^-1
    Matrix a_operator; // Matrix (N,N) where each row i is the operator to apply at time t_i

    // Private helper methods
    void construct_time_integral();
    void construct_L_and_U();
    void construct_inv_Dt();
    void construct_B_and_a_operator();
    inline void update_state_variables();

public:
    // Constructor
    NystromScheme(const std::size_t M, const std::size_t N, const Vector& time_grid, 
        const double X0, 
        const Kernel& k, 
        const std::unique_ptr<Matrix>& u, Matrix& Z_u, const std::unique_ptr<Matrix>& X_u, 
        const std::vector<Matrix>& R, 
        const std::vector<Matrix>& Gamma);

    // Update control variable u and state variables Z_u and X_u
    void nystrom_update();
};

#endif // NYSTROM_SCHEME_HPP