#include "NystromScheme.hpp"

NystromScheme::NystromScheme(const std::size_t M, const std::size_t N, const Vector& time_grid, 
    const double X0, 
    const Kernel& k, 
    const std::unique_ptr<Matrix>& u, Matrix& Z_u, const std::unique_ptr<Matrix>& X_u, 
    const std::vector<Matrix>& R, 
    const std::vector<Matrix>& Gamma)
    : M(M), N(N), I_N(Matrix::Identity(N, N)), time_grid(time_grid), X0(Matrix::Constant(M, N, X0)), 
      kernel(k), u(u), Z_u(Z_u), X_u(X_u), R(R), Gamma(Gamma)
{
    // Construct Matrices for operators that do not depend on the control variable u
    construct_time_integral();
    construct_L_and_U(); // Construct the L and U matrices based on the kernel and time grid
    construct_inv_Dt(); // Construct the inverse of (I +K_t +K*_t) for each time step
    construct_B_and_a_operator();
    // Reserve space for the R_sum_Gamma vector of dimension M
    R_sum_Gamma.reserve(M);
}

void NystromScheme::construct_time_integral() {
    // Each row of time_integral contains the weights to write: integral_0^t us ds = sum_{j=0}^{N-1} us * (t_{j+1} - t_j)
    Vector time_delta = time_grid.tail(N) - time_grid.head(N);
    time_integral = Matrix::Zero(N, N);
    for (std::size_t i = 0; i < N; ++i) {
        time_integral.row(i).head(i+1) = time_delta.head(i+1);
    }
}

void NystromScheme::construct_L_and_U() {
    // L(i,j) = ∫_{t_j}^{t_{j+1}} K(t_i, s) ds
    // U(i,j) = ∫_{t_j}^{t_{j+1}} K(s, t_i) ds
    L = Matrix::Zero(N, N);
    U = Matrix::Zero(N, N);
    
    // Parallelize the construction of L and U matrices
    // #pragma omp parallel for schedule(static)
    for (par_for_type i = 0; i < static_cast<par_for_type>(N); ++i) {
        for (std::size_t j = 0; j < i; ++j) { // int_tj^tj+1 of K(t_i, s)
            L(i, j) = kernel.s_integral(time_grid(i), time_grid(j), time_grid(j + 1));
        }
        for (std::size_t j = i; j < N; ++j) { // int_tj^tj+1 of K(s, t_i)
            U(i, j) = kernel.t_integral(time_grid(i), time_grid(j), time_grid(j + 1));
        }
    }
}

void NystromScheme::construct_inv_Dt() {
    // Construct the inverse of (I +K_t +K*_t) for each time step
    inv_Dt.resize(N);
    inv_Dt[0] = (I_N + L + U).inverse();
    
    // #pragma omp parallel for schedule(static)
    for(par_for_type i = 1; i < static_cast<par_for_type>(N); ++i) {
        // Matrix Kt_i is matrix L where the first i columns and i rows are zero
        Matrix Kt_i = Matrix::Zero(N, N);
        Kt_i.bottomRightCorner(N-i, N-i) = L.bottomRightCorner(N-i, N-i);

        // Matrix K*_t_i is matrix U where the first i columns and i rows are zero
        Matrix K_star_t_i = Matrix::Zero(N, N);
        K_star_t_i.bottomRightCorner(N-i, N-i) = U.bottomRightCorner(N-i, N-i);

        // Compute the inverse of (I + Kt_i + K*_t_i)
        inv_Dt[i] = (I_N + Kt_i + K_star_t_i).inverse();
    }
}

void NystromScheme::construct_B_and_a_operator() {
    // Construct operators B and a_operator
    B = Matrix::Zero(N, N);
    a_operator = Matrix::Zero(N, N);
    
    // Parallelize the construction of a_operator and B
    // #pragma omp parallel for schedule(static)
    for (par_for_type i = 0; i < static_cast<par_for_type>(N); ++i) {
        a_operator.row(i) = U.row(i) * inv_Dt[i];
        for(std::size_t j = 0; j < i; ++j) {
            B(i,j) = a_operator.row(i) * L.col(j);
        }
    }
    // Compute the control variable u based on the Nystrom scheme
    I_B_inv = (I_N - B).inverse();
}

void NystromScheme::update_state_variables() {
    // Update the state variables Z_u and X_u based on the control variable u
    // Update the transient price impact Z_u based on the current control variable u: Z_u(t) = integral_0^t u(s)*K(t,s)ds
    Z_u.noalias() = (*u) * L.transpose();
    // Update the inventory X_u based on the current control variable u: X_u(t) = X0 + integral_0^t u(s)ds
    X_u->noalias() = X0 + (*u) * time_integral.transpose();
}

void NystromScheme::nystrom_update() {
    // Update control variable u and state variables Z_u and X_u
    // Update the R_sum_Gamma vector with R+Gamma at each m<M
    R_sum_Gamma.clear();
    R_sum_Gamma.resize(M);
    
    // Parallelize the R_sum_Gamma computation
    // #pragma omp parallel for schedule(static)
    for (par_for_type m = 0; m < static_cast<par_for_type>(M); ++m) {
        R_sum_Gamma[m] = R[m] + Gamma[m];
    }
    
    // Update the control variable u based on the Nystrom scheme
    // Parallelize over Monte Carlo paths
    // #pragma omp parallel for schedule(static)
    for (par_for_type m = 0; m < static_cast<par_for_type>(M); ++m) {
        Vector a(N);
        for (std::size_t i = 0; i < N; ++i) {
            a[i] = R_sum_Gamma[m](i,i) - a_operator.row(i).dot(R_sum_Gamma[m].row(i));
        }
        u->row(m) = I_B_inv * a;
    }
    update_state_variables();
}
