#ifndef LSMCR_HPP
#define LSMCR_HPP

#include <memory>
#include <numeric>
#include <cstddef>

#include "CommonTypes.hpp"

// Forward declaration of LaguerrePolynomial function
inline Matrix LaguerrePolynomial(const Vector& x, const int degree);

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

    // Private member functions - implementations in LSMCR.cpp
    void precompute_target_i();
    void precompute_target_ij();
    void precompute_regressors_i();
    void precompute_regressors_ij();
    void regression_i();
    void regression_ij();
    void estimate_conditional_expectation_i();
    void estimate_conditional_expectation_ij();

public:
    // Constructor - implementation in LSMCR.cpp
    LSMCR(const std::size_t d1, const std::size_t d2, 
        const std::size_t N, const double time_delta, const std::size_t M, 
        const Matrix& alpha, const Matrix& Z_u, const std::unique_ptr<Matrix>& X_u, const std::vector<std::unique_ptr<Matrix>>& lambda);

    // Public methods - implementations in LSMCR.cpp
    void update_regressor();
    
    // Getter functions to access coefficients for testing/debugging
    const Matrix& get_coeff_i() const { return coeff_i; }
    const MatVec& get_coeff_ij() const { return coeff_ij; }

    MatVec estimate_gamma();

};

#endif // LSMCR_HPP