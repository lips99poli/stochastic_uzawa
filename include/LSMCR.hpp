#ifndef LSMCR_HPP
#define LSMCR_HPP

#include <stdexcept>
#include <memory>
#include <vector>
#include <numeric>
#include <cstddef>

#include "CommonTypes.hpp"

// Forward declaration of LaguerrePolynomial function
Matrix LaguerrePolynomial(const Vector& x, const int degree);

// This class implements 2 LSMCRegressions:
// 1. LSMCR: Least Squares Monte Carlo Regression for the conditional expectation E_i[ sum_l_i+1_N {lamda3_l - lambda4_l} ] with 0<i<=N-1
//  The output of this regression are (d1+1  2) coefficients for each time step i
//  so it is a matrix of size ((d1+1  2), N-1)
// 2. LSMCR2: Least Squares Monte Carlo Regression for the conditional expectation E_i[ lamda1_j - lamda2_j + sum_l_j+1_N {lambda3_l - lambda2_4} ] with 0<i<j<=N-1
// The output of this regression are (d2+1  2) coefficients for each time step couples (i,j)
//  so it is a matrix of size ((d2+1  2), N-2, N-1)

// la time grid è lunga (N+1) e contiene i tempi da 0 a T, ma T a noi non interessa, infatti le variabili di stato sono lunghe N cioè valorizzate sui tempi 0,...,T*(N-1)/N
// la regressione ha senso fino all'ultimo valore delle variabili di stato e con le info fino al precedente -> output indici 1 a N-1 e input da 0 a N-2 (entrambi lunghi N-1)

class LSMCR{
    private:

    public:
    const std::size_t d1; // Upper bound for the sum of the degrees of the Laguerre polynomials used in the first regression
    const std::size_t d2; // Upper bound for the sum of the degrees of the Laguerre polynomials used in the second regression
    const std::size_t comb_upto_d1, comb_upto_d2;

    const std::size_t N; // Number of time steps
    const double time_delta;
    const std::size_t M; // Number of sample paths

    const Matrix& alpha; // Alpha matrix, of size (M, N)
    const Matrix& Z_u; // Z_u matrix, of size (M, N)
    const std::unique_ptr<Matrix>& X_u; // X_u matrix, of size (M, N)
    const std::vector<std::unique_ptr<Matrix>>& lambda; // Lambda vector, of size (4, M, N)

    Matrix diff_lambda_1_2; // Difference lambda1 - lambda2, of size (M, N)
    Matrix diff_lambda_3_4; // Difference lambda3 - lambda4, of size (M, N)

    // Laguerre polynomials for the three variables of state, vector length N containinig matrixes of size (M, max(d1,d2)+1), (M, max(d1,d2)+1), (M, max(d1,d2)+1)
    MatVec laguerre_alpha; // Laguerre polynomials for alpha, of size (M, max(d1,d2)+1)
    MatVec laguerre_Z_u; // Laguerre polynomials for Z_u, of size (M, max(d1,d2)+1)
    MatVec laguerre_X_u; // Laguerre polynomials for X_u, of size (M, max(d1,d2)+1)

    // Targets for the two regressions
    Matrix target_i; // Target matrix for the first regression, of size (N,M): for every timestamp i I have a target vector of size M corresponding to time i+1
    Matrix target_ij; // Target matrix for the second regression, of size (N,M): for every timestamp couple j>i (i indicates the time of regressors) I have a target vector of size M

    // Regressors for the 2 regressions
    MatVec regressors_i; // Vector of length N where each element is the matrix PHI_i of size((d+3  3), M) 
    MatVec regressors_ij; // Vector of length N where each element is the matrix PHI_i of size((d+3  3), M)

    // Coefficients for the two regressions
    Matrix coeff_i; // Coefficients matrix for the first regression, of size ((d1+3  3), N-1)
    MatVec coeff_ij; // Coefficients map for the second regression, of size (N-1,(d2+3  3),N-1-i): for each i going from 0 to N-2 i do regression for target at all j in [i+1,N-1] and i get (d2+3  3) coefficientes 

    // Estimates for the two conditional expectations
    Matrix cond_exp_i;
    MatVec cond_exp_ij;


    // Dubbio: non ho capito da quando inizia la regressione: i=0 o i=1?

    // Note: lambda1 and lambda2 are of length N referring to times 0,...,T*(N-1)/N
    //       lambda3 and lambda4 are of length N referring to times 1,...,T

    // This function will precompute the target vector for the first regression
    void precompute_target_i(){
        target_i = Matrix::Zero(M, N); // Initialize the target vector for the first regression
        // Compute cumulative sum from right to left for efficiency
        Vector cumsum = Vector::Zero(M);
        for (int j = static_cast<int>(N)-1; j >= 0; --j) {
            cumsum += diff_lambda_3_4.col(j); // Add current column to cumulative sum
            target_i.col(j) = cumsum * time_delta; // Store for time step j
        }
    }

    // This function will precompute the target vector for the second regression
    void precompute_target_ij(){
        target_ij = Matrix::Zero(M, N-1);
        for (std::size_t j=0; j<N-1; ++j){
            // Recall that targets are already shifted by 1 in time
            target_ij.col(j) = target_i.col(j+1) + diff_lambda_1_2.col(j+1);
        }
    }
    
    // This function will precompute the regressors for the first regression
    void precompute_regressors_i(){
        regressors_i.clear();
        for(std::size_t i=0; i<N; ++i){
            // The matrix PHI_i will be of dimension Mx(d1+3 choose 3), where M is the number of sample paths.
            // Each column of PHI_i will contain all possible combinations of multiplication of Laguerre polynomials of the three variables such that the sum of the degrees is less than or equal to d1.
            regressors_i.emplace_back(M, comb_upto_d1); // Construct matrix directly in vector
            Matrix& PHI_i = regressors_i.back(); // Get reference to the matrix
            PHI_i.setZero(); // Initialize to zero
            std::size_t col_idx = 0; // Column index for PHI_i
            for (std::size_t x = 0; x <= d1; ++x) {
                for (std::size_t y = 0; y <= d1-x; ++y) {
                    for (std::size_t z = 0; z <= d1-x-y; ++z) {
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
    void precompute_regressors_ij(){
        regressors_ij.clear();
        for(std::size_t i=0; i<N-1; ++i){
            // The matrix PHI_i will be of dimension Mx(d2+3 choose 3), where M is the number of sample paths.
            // Each column of PHI_i will contain all possible combinations of multiplication of Laguerre polynomials of the three variables such that the sum of the degrees is less than or equal to d2.
            regressors_ij.emplace_back(M, comb_upto_d2); // Construct matrix directly in vector
            Matrix& PHI_i = regressors_ij.back(); // Get reference to the matrix
            PHI_i.setZero(); // Initialize to zero
            std::size_t col_idx = 0; // Column index for PHI_i
            for (std::size_t x = 0; x <= d2; ++x) {
                for (std::size_t y = 0; y <= d2-x; ++y) {
                    for (std::size_t z = 0; z <= d2-x-y; ++z) {
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
    void regression_i(){
        // The Linear Regression will solve for each time step i the following equation: PHI_i * coeffs_i = Y_i giving as output coeffs_i of size (d1+3  3) x 1
        // The final output will be the matrix coeffs_i of size ((d1+3  3), N)
        
        // Cicle over the time steps i from 0 to N-1
        for (std::size_t i = 0; i < N; ++i) {
            // Construct matrix PHI_i
            const Matrix& PHI_i = regressors_i[i];

            // Solve the linear regression problem using least squares QR decomposition
            coeff_i.col(i) = PHI_i.colPivHouseholderQr().solve(target_i.col(i));
        }
    }

    // Function to compute the coefficients of the second LSMCR regression
    void regression_ij(){
        // The Linear Regression will solve for each time step couple (i,j) s.t. 0<=i<j<N the following equation: PHI_i * coeffs_ij = Y_j giving as output coeffs_ij of size (d2+3  3) x 1
        // The final output will be the vector of Matrix coeffs_ij of size (N-1, (d2+3  3), N-1-i)
        coeff_ij.clear(); // Clear the previous coefficients
        for (std::size_t i = 0; i < N-1; ++i) {
            coeff_ij.emplace_back(comb_upto_d2, N-1-i); // Create a new Matrix for each i
            // Loop over j from i+1 to N-1 to gill the columns of coeff_ij[i]
            for (std::size_t j=i; j<N-1; ++j) {
                coeff_ij[i].col(j-i) = regressors_ij[i].colPivHouseholderQr().solve(target_ij.col(j));
            }
        }
    }

    // Function to estimate the conditional expectation E_i[ sum_l_i+1_N {lamda3_l - lambda4_l} ] given the coefficients of the first regression
    void estimate_conditional_expectation_i() {
        cond_exp_i = Matrix::Zero(M, N);
        for (std::size_t i=0; i<N; ++i){
            cond_exp_i.col(i) = regressors_i[i] * coeff_i.col(i);
        }
    }
    // Function to estimate the conditional expectation E_i[ lamda1_j - lamda2_j + sum_l_j+1_N {lambda3_l - lambda2_4} ] given the coefficients of the second regression
    void estimate_conditional_expectation_ij() {
        cond_exp_ij.clear(); // Reserve space for N-1 matrices
        for (std::size_t i=0; i<N-1; ++i){
            cond_exp_ij.emplace_back(M, N-1-i); // Create a new Matrix for each i
            for (std::size_t j=i; j<N-1; ++j){
                cond_exp_ij[i].col(j-i) = regressors_ij[i] * coeff_ij[i].col(j-i);
            }
        }
    }

    public:
    LSMCR(const std::size_t d1, const std::size_t d2, const std::size_t N, const double time_delta, const std::size_t M, const Matrix& alpha, const Matrix& Z_u, const std::unique_ptr<Matrix>& X_u, const std::vector<std::unique_ptr<Matrix>>& lambda)
        : d1(d1), d2(d2), comb_upto_d1((d1+3)*(d1+2)*(d1+1)/6), comb_upto_d2((d2+3)*(d2+2)*(d2+1)/6), N(N), time_delta(time_delta), M(M), 
          alpha(alpha), Z_u(Z_u), X_u(X_u), lambda(lambda)
        {
        if (alpha.size() == 0|| Z_u.size() == 0 || !X_u) {
            throw std::invalid_argument("Invalid input for state variables.");
        }
        laguerre_alpha.reserve(N);
        laguerre_Z_u.reserve(N);
        laguerre_X_u.reserve(N);
        regressors_i.reserve(N);
        regressors_ij.reserve(N);
        // Allocate container for coefficients
        coeff_i = Matrix::Zero(comb_upto_d1, N); // (d1+1 choose 2) coefficients for each time step i
        coeff_ij.reserve(N-1);
    }

    void update_regressor(){
        // LAGUERRE POLYNOMIALS
        // Compute Laguerre polynomials of the 3 variables up to degrees max(d1,d2)
        laguerre_alpha.clear();
        laguerre_Z_u.clear();
        laguerre_X_u.clear();
        // Compute the Laguerre polynomials for each variable
        int degree = std::max(d1, d2);
        for (std::size_t i = 0; i < N; ++i) {
            laguerre_alpha.push_back(LaguerrePolynomial(alpha.col(i), degree));
            laguerre_Z_u.push_back(LaguerrePolynomial(Z_u.col(i), degree));
            laguerre_X_u.push_back(LaguerrePolynomial((*X_u).col(i), degree));
        }
        
        // LAMBDA1 - LAMBDA2 and LAMBDA3 - LAMBDA4
        diff_lambda_1_2 = (*lambda[0]) - (*lambda[1]);
        diff_lambda_3_4 = (*lambda[2]) - (*lambda[3]);

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


    // Getter functions to access coefficients for testing/debugging
    const Matrix& get_coeff_i() const { return coeff_i; }
    const MatVec& get_coeff_ij() const { return coeff_ij; }

    MatVec estimate_gamma(){
        MatVec Gamma(M, Matrix::Zero(N, N));
        for (std::size_t m=0; m<M; ++m){
            // Fill diagonal
            Gamma[m].diagonal() = cond_exp_i.row(m) + diff_lambda_1_2.row(m);
            // Fill upper triangle
            for (std::size_t i=0; i<N-1; ++i){
                // Recall cond_exp_ij: fix i as the time of the expectation than we have a matrix of size (M, N-1-i) with the conditional expectation for each j>i up to N-1
                // The rows in the upper triangle (diagonal excluded) are filled for each m with the rows cond_exp_ij[i](m,:)
                for (std::size_t j=i+1; j<N; ++j){
                    Gamma[m](i,j) = cond_exp_ij[i](m,j-1-i);
                }
            }
        }
        return Gamma;
    }

};


// Function to compute Laguerre polinomial on a given vector of a given degree.
// For the implementation it used a recursive implementation: L_n(x) = (2n-1-x)L_{n-1}(x) - (n-1)L_{n-2}(x) / n
Matrix LaguerrePolynomial(const Vector& x, const int degree){
    if (degree < 0) {
        throw std::invalid_argument("Degree must be non-negative");
    }
    
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

#endif // LSMCR_HPP