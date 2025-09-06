#include "LSMCR.hpp"

// Constructor implementation
LSMCR::LSMCR(const std::size_t d1, const std::size_t d2, 
    const std::size_t N, const double time_delta, const std::size_t M, 
    const Matrix& alpha, const Matrix& Z_u, const std::unique_ptr<Matrix>& X_u, const std::vector<std::unique_ptr<Matrix>>& lambda)
    : d1(d1), d2(d2), comb_upto_d1((d1+3)*(d1+2)*(d1+1)/6), comb_upto_d2((d2+3)*(d2+2)*(d2+1)/6), N(N), time_delta(time_delta), M(M), 
      alpha(alpha), Z_u(Z_u), X_u(X_u), lambda(lambda)
    {
    laguerre_alpha.reserve(N);
    laguerre_Z_u.reserve(N);
    laguerre_X_u.reserve(N);
    regressors_i.reserve(N);
    regressors_ij.reserve(N);
    // Allocate container for coefficients
    coeff_i = Matrix::Zero(comb_upto_d1, N); // (d1+1 choose 2) coefficients for each time step i
    coeff_ij.reserve(N-1);
}

// This function will precompute the target vector for the first regression
void LSMCR::precompute_target_i(){
    target_i = Matrix::Zero(M, N); // Initialize the target vector for the first regression
    
    // Compute cumulative sum from right to left for efficiency
    // Parallelize over paths (rows)
    // #pragma omp parallel for schedule(static)
    for (par_for_type m = 0; m < static_cast<par_for_type>(M); ++m) {
        double cumsum = 0.0;
        for (int j = static_cast<int>(N)-1; j >= 0; --j) {
            cumsum += diff_lambda_3_4(m, j); // Add current element to cumulative sum
            target_i(m, j) = cumsum * time_delta; // Store for time step j
        }
    }
}

// This function will precompute the target vector for the second regression
void LSMCR::precompute_target_ij(){
    target_ij = Matrix::Zero(M, N-1);
    
    // Parallelize over time steps
    // #pragma omp parallel for schedule(static)
    for (par_for_type j = 0; j < static_cast<par_for_type>(N-1); ++j){
        // Recall that target_i are already shifted by 1 in time
        target_ij.col(j) = target_i.col(j+1) + diff_lambda_1_2.col(j+1);
    }
}

// This function will precompute the regressors for the first regression
void LSMCR::precompute_regressors_i(){
    regressors_i.clear();
    regressors_i.resize(N);
    
    // Parallelize over time steps
    // #pragma omp parallel for schedule(static)
    for(par_for_type i = 0; i < static_cast<par_for_type>(N); ++i){
        // The matrix PHI_i will be of dimension Mx(d1+3 choose 3), where M is the number of sample paths.
        // Each column of PHI_i will contain all possible combinations of multiplication of Laguerre polynomials of the three variables such that the sum of the degrees is less than or equal to d1.
        regressors_i[i] = Matrix::Zero(M, comb_upto_d1); // Initialize matrix
        Matrix& PHI_i = regressors_i[i]; // Get reference to the matrix
        
        std::size_t col_idx = 0; // Column index for PHI_i
        for (std::size_t x = 0; x <= d1; ++x) {
            for (std::size_t y = 0; y <= d1-x; ++y) {
                for (std::size_t z = 0; z <= d1-x-y; ++z) {
                    // Vectorized multiplication for all paths at once
                    PHI_i.col(col_idx) = laguerre_alpha[i].col(x).array() *
                                        laguerre_Z_u[i].col(y).array() *
                                        laguerre_X_u[i].col(z).array();
                    col_idx++;
                }
            }
        }
    }
}

// This function will precompute the regressors for the second regression
void LSMCR::precompute_regressors_ij(){
    regressors_ij.clear();
    regressors_ij.resize(N-1);
    
    // Parallelize over time steps
    // #pragma omp parallel for schedule(static)
    for(par_for_type i = 0; i < static_cast<par_for_type>(N-1); ++i){
        // The matrix PHI_i will be of dimension Mx(d2+3 choose 3), where M is the number of sample paths.
        // Each column of PHI_i will contain all possible combinations of multiplication of Laguerre polynomials of the three variables such that the sum of the degrees is less than or equal to d2.
        regressors_ij[i] = Matrix::Zero(M, comb_upto_d2); // Initialize matrix
        Matrix& PHI_i = regressors_ij[i]; // Get reference to the matrix
        
        std::size_t col_idx = 0; // Column index for PHI_i
        for (std::size_t x = 0; x <= d2; ++x) {
            for (std::size_t y = 0; y <= d2-x; ++y) {
                for (std::size_t z = 0; z <= d2-x-y; ++z) {
                    // Vectorized multiplication for all paths at once
                    PHI_i.col(col_idx) = laguerre_alpha[i].col(x).array() *
                                        laguerre_Z_u[i].col(y).array() *
                                        laguerre_X_u[i].col(z).array();
                    col_idx++;
                }
            }
        }
    }
}

// Function to compute the coefficients of the first LSMCR regression
void LSMCR::regression_i(){
    // The Linear Regression will solve for each time step i the following equation: PHI_i * coeffs_i = Y_i giving as output coeffs_i of size (d1+3  3) x 1
    // The final output will be the matrix coeffs_i of size ((d1+3  3), N)
    
    // Parallelize over time steps - each QR decomposition is independent
    // #pragma omp parallel for schedule(static)
    for (par_for_type i = 0; i < static_cast<par_for_type>(N); ++i) {
        // Construct matrix PHI_i
        const Matrix& PHI_i = regressors_i[i];

        // Solve the linear regression problem using least squares QR decomposition
        coeff_i.col(i) = PHI_i.colPivHouseholderQr().solve(target_i.col(i));
    }
}

// Function to compute the coefficients of the second LSMCR regression
void LSMCR::regression_ij(){
    // The Linear Regression will solve for each time step couple (i,j) s.t. 0<=i<j<N the following equation: PHI_i * coeffs_ij = Y_j giving as output coeffs_ij of size (d2+3  3) x 1
    // The final output will be the vector of Matrix coeffs_ij of size (N-1, (d2+3  3), N-1-i)
    coeff_ij.clear(); // Clear the previous coefficients
    coeff_ij.resize(N-1);
    
    // Initialize matrices for each i
    for (std::size_t i = 0; i < N-1; ++i) {
        coeff_ij[i] = Matrix::Zero(comb_upto_d2, N-1-i); // Create a new Matrix for each i
    }
    
    // Parallelize over all (i,j) pairs
    // #pragma omp parallel for collapse(2)
    for (par_for_type i = 0; i < static_cast<par_for_type>(N-1); ++i) {
        for (par_for_type j_idx = 0; j_idx < static_cast<par_for_type>(N-1-i); ++j_idx) {
            std::size_t j = i + j_idx;
            coeff_ij[i].col(j_idx) = regressors_ij[i].colPivHouseholderQr().solve(target_ij.col(j));
        }
    }
}

// Function to estimate the conditional expectation E_i[ sum_l_i+1_N {lamda3_l - lambda4_l} ] given the coefficients of the first regression
void LSMCR::estimate_conditional_expectation_i() {
    cond_exp_i = Matrix::Zero(M, N);
    
    // Parallelize over time steps
    // #pragma omp parallel for schedule(static)
    for (par_for_type i = 0; i < static_cast<par_for_type>(N); ++i){
        cond_exp_i.col(i) = regressors_i[i] * coeff_i.col(i);
    }
}

// Function to estimate the conditional expectation E_i[ lamda1_j - lamda2_j + sum_l_j+1_N {lambda3_l - lambda2_4} ] given the coefficients of the second regression
void LSMCR::estimate_conditional_expectation_ij() {
    cond_exp_ij.clear();
    cond_exp_ij.resize(N-1);
    
    // Initialize matrices
    for (std::size_t i=0; i<N-1; ++i){
        cond_exp_ij[i] = Matrix::Zero(M, N-1-i); // Create a new Matrix for each i
    }
    
    // Parallelize over all (i,j) pairs
    // #pragma omp parallel for collapse(2)
    for (par_for_type i = 0; i < static_cast<par_for_type>(N-1); ++i){
        for (par_for_type j_idx = 0; j_idx < static_cast<par_for_type>(N-1-i); ++j_idx){
            cond_exp_ij[i].col(j_idx) = regressors_ij[i] * coeff_ij[i].col(j_idx);
        }
    }
}

// Main update regressor method
void LSMCR::update_regressor(){
    // LAGUERRE POLYNOMIALS
    // Compute Laguerre polynomials of the 3 variables up to degrees max(d1,d2)
    laguerre_alpha.clear();
    laguerre_Z_u.clear();
    laguerre_X_u.clear();
    laguerre_alpha.resize(N);
    laguerre_Z_u.resize(N);
    laguerre_X_u.resize(N);
    
    // Compute the Laguerre polynomials for each variable in parallel
    int degree = std::max(d1, d2);
    // #pragma omp parallel for schedule(static)
    for (par_for_type i = 0; i < static_cast<par_for_type>(N); ++i) {
        laguerre_alpha[i] = LaguerrePolynomial(alpha.col(i), degree);
        laguerre_Z_u[i] = LaguerrePolynomial(Z_u.col(i), degree);
        laguerre_X_u[i] = LaguerrePolynomial((*X_u).col(i), degree);
    }
    
    // LAMBDA1 - LAMBDA2 and LAMBDA3 - LAMBDA4
    // Parallelize matrix operations
    // #pragma omp parallel sections
    {
        // #pragma omp section
        {
            diff_lambda_1_2 = (*lambda[0]) - (*lambda[1]);
        }
        // #pragma omp section
        {
            diff_lambda_3_4 = (*lambda[2]) - (*lambda[3]);
        }
    }

    // TARGETS
    // Precompute the targets
    precompute_target_i();
    precompute_target_ij();

    // REGRESSORS
    precompute_regressors_i();
    precompute_regressors_ij();

    // COEFFICIENTS
    regression_i();
    regression_ij();

    // CONDITIONAL EXPECTATIONS
    estimate_conditional_expectation_i();
    estimate_conditional_expectation_ij();
}

// Estimate gamma method
MatVec LSMCR::estimate_gamma(){
    MatVec Gamma(M, Matrix::Zero(N, N));
    
    // Parallelize over simulation paths
    // #pragma omp parallel for schedule(static)
    for (par_for_type m = 0; m < static_cast<par_for_type>(M); ++m){
        // Fill diagonal
        Gamma[m].diagonal() = cond_exp_i.row(m) + diff_lambda_1_2.row(m);
        // Fill upper triangle
        for (std::size_t i = 0; i < N-1; ++i){
            // Recall cond_exp_ij: fix i as the time of the expectation than we have a matrix of size (M, N-1-i) with the conditional expectation for each j>i up to N-1
            // The rows in the upper triangle (diagonal excluded) are filled for each m with the rows cond_exp_ij[i](m,:)
            for (std::size_t j = i+1; j < N; ++j){
                Gamma[m](i,j) = cond_exp_ij[i](m,j-1-i);
            }
        }
    }
    return Gamma;
}

// Function to compute Laguerre polinomial on a given vector of a given degree.
// For the implementation it used a recursive implementation: L_n(x) = (2n-1-x)L_{n-1}(x) - (n-1)L_{n-2}(x) / n
Matrix LaguerrePolynomial(const Vector& x, const int degree){
    Matrix result(x.size(), degree+1);
    result.col(0) = Vector::Ones(x.size()); // L_0(x) = 1
    if (degree >= 1) {
        result.col(1) = result.col(0) - x; // L_1(x) = 1-x
    }
    
    for (int n = 2; n <= degree; ++n) {
        result.col(n) = ((2 * n - 1 - x.array()) * result.col(n - 1).array() - (n - 1) * result.col(n - 2).array()) / n;
    }
    return result;
}