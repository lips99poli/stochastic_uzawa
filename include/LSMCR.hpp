#ifndef LSMCR_HPP
#define LSMCR_HPP

#include <stdexcept>
#include <memory>
#include <vector>
#include <map>
#include <numeric>
#include <cstddef>
#include <Eigen/Dense>



//create aliases for Eigen types
using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using Ref = Eigen::Ref<Matrix>;

using CondExp_ij = std::vector<Matrix>; // Map to store conditional expectations for each time step couple (i,j)

// Forward declaration of LaguerrePolynomial function
Matrix LaguerrePolynomial(const Vector& x, const int degree);

// This class implements 2 LSMCRegressions:
// 1. LSMCR: Least Squares Monte Carlo Regression for the conditional expectation E_i[ sum_l_i+1_N {lamda3_l - lambda4_l} ] with 0<i<=N-1
//  The output of this regression are (d1+1  2) coefficients for each time step i
//  so it is a matrix of size ((d1+1  2), N-1)
// 2. LSMCR2: Least Squares Monte Carlo Regression for the conditional expectation E_i[ lamda1_j - lamda2_j + sum_l_j+1_N {lambda3_l - lambda2_4} ] with 0<i<j<=N-1
// The output of this regression are (d2+1  2) coefficients for each time step couples (i,j)
//  so it is a matrix of size ((d2+1  2), N-2, N-1)



// NOTE: SALVARE GLI SVILUPPI DI LAGUERRE DELLE 3 VARIABILI, DATO D E D' LI SVILUPPO TUTTI E TRE FINO A MAX(D/3,D'/3,1), QUANDO FARò L'OTTIMIZZAZIONE NON ANDRò MAI OLTRE QUESTO VALORE
// li salvo come membri e poi pesco da li sia per la regressione che per la stima
// li salvo in dei vec di vec/array di dim: con m,i accedo all'espansione fino al grado dim 
// meglio sostituire anche la mappa con questi vec di vec/array con un vec di vec di vec/array coì faccio il tensore. Questo perchè ho accesso più rapido.
// poi ci sono da fare le due funzioni di stima


// la time grid è lunga (N+1) e contiene i tempi da 0 a T, ma T a noi non interessa, infatti le variabili di stato sono lunghe N cioè valorizzate sui tempi 0,...,T*(N-1)/N
// la regressione ha senso fino all'ultimo valore delle variabili di stato e con le info fino al precedente -> output indici 1 a N-1 e input da 0 a N-2 (entrambi lunghi N-1)

class LSMCR{
    private:
    const std::size_t d1; // Upper bound for the sum of the degrees of the Laguerre polynomials used in the first regression
    const std::size_t d2; // Upper bound for the sum of the degrees of the Laguerre polynomials used in the second regression
    const std::size_t comb_upto_d1, comb_upto_d2;

    const std::size_t N; // Number of time steps
    const std::size_t M; // Number of sample paths

    const std::unique_ptr<Matrix>& alpha; // Alpha matrix, of size (M, N)
    const std::unique_ptr<Matrix>& Z_u; // Z_u matrix, of size (M, N)
    const std::unique_ptr<Matrix>& X_u; // X_u matrix, of size (M, N)
    const std::vector<std::unique_ptr<Matrix>>& lambda; // Lambda vector, of size (4, M, N+1)

    Matrix diff_lambda_3_4; // Difference lambda3 - lambda4, of size (M, N+1)

    // Laguerre polynomials for the three variables of state, vector length N containinig matrixes of size (M, max(d1,d2)+1), (M, max(d1,d2)+1), (M, max(d1,d2)+1)
    std::vector<Matrix> laguerre_alpha; // Laguerre polynomials for alpha, of size (M, max(d1,d2)+1)
    std::vector<Matrix> laguerre_Z_u; // Laguerre polynomials for Z_u, of size (M, max(d1,d2)+1)
    std::vector<Matrix> laguerre_X_u; // Laguerre polynomials for X_u, of size (M, max(d1,d2)+1)

    // Targets for the two regressions
    Matrix target_i; // Target matrix for the first regression, of size (N,M): for every timestamp i I have a target vector of size M corresponding to time i+1
    Matrix target_ij; // Target matrix for the second regression, of size (N,M): for every timestamp couple j>i (i indicates the time of regressors) I have a target vector of size M

    // Regressors for the 2 regressions
    std::vector<Matrix> regressors_i; // Vector of length N where each element is the matrix PHI_i of size((d+3  3), M) 
    std::vector<Matrix> regressors_ij; // Vector of length N where each element is the matrix PHI_i of size((d+3  3), M)

    // Coefficients for the two regressions
    Matrix coeff_i; // Coefficients matrix for the first regression, of size ((d1+3  3), N-1)
    CondExp_ij coeff_ij; // Coefficients map for the second regression, of size (N-1,(d2+3  3),N-1-i): for each i going from 0 to N-2 i do regression for target at all j in [i+1,N-1] and i get (d2+3  3) coefficientes 

    // Dubbio: non ho capito da quando inizia la regressione: i=0 o i=1?

    // This function will precompute the target vector for the first regression
    void precompute_target_i(){
        target_i = Matrix::Zero(M, N); // Initialize the target vector for the first regression
        // Compute cumulative sum from right to left for efficiency
        Vector cumsum = Vector::Zero(M);
        for (std::size_t j=N; j>=1; --j) {
            cumsum += diff_lambda_3_4.col(j); // Add current column to cumulative sum
            target_i.col(j-1) = cumsum; // Store for time step i-1
        }
    }

    // This function will precompute the target vector for the second regression
    void precompute_target_ij(){
        target_ij = Matrix::Zero(M, N-1);
        for (std::size_t j=0; j<N-1; ++j){
            // Recall that targets are already shifted by 1 in time
            target_ij.col(j) = target_i.col(j+1) + lambda[0]->col(j+1) - lambda[1]->col(j+1);
        }
    }
    
    // This function will precompute the regressors for the first regression
    void precompute_regressors_i(){
        regressors_i.clear();
        regressors_i.reserve(N); // Reserve space to avoid reallocations
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
        regressors_ij.reserve(N); // Reserve space to avoid reallocations
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
        
        // Initialize the output matrix for coefficients
        coeff_i = Matrix::Zero(comb_upto_d1, N); // (d1+1 choose 2) coefficients for each time step i

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

    public:
    LSMCR(const std::size_t d1, const std::size_t d2, const std::size_t N, const std::size_t M, const std::unique_ptr<Matrix>& alpha, const std::unique_ptr<Matrix>& Z_u, const std::unique_ptr<Matrix>& X_u, const std::vector<std::unique_ptr<Matrix>>& lambda)
        : d1(d1), d2(d2), comb_upto_d1((d1+3)*(d1+2)*(d1+1)/6), comb_upto_d2((d2+3)*(d2+2)*(d2+1)/6), N(N), M(M), 
          alpha(alpha), Z_u(Z_u), X_u(X_u), lambda(lambda)
        {
        if (!alpha || !Z_u || !X_u) {
            throw std::invalid_argument("Invalid input for state variables.");
        }
        regressors_i.reserve(N);
        coeff_ij.reserve(N-1);
    }

    void update_regressor(){
        // LAGUERRE POLYNOMIALS: Compute Laguerre polynomials of the 3 variables up to degrees max(d1,d2)
        laguerre_alpha.clear();
        laguerre_Z_u.clear();
        laguerre_X_u.clear();
        // Reserve space for the vectors
        laguerre_alpha.reserve(N);
        laguerre_Z_u.reserve(N);
        laguerre_X_u.reserve(N);
        // Compute the Laguerre polynomials for each variable
        int degree = std::max(d1, d2);
        for (std::size_t i = 0; i < N; ++i) {
            laguerre_alpha.push_back(LaguerrePolynomial((*alpha).col(i), degree));
            laguerre_Z_u.push_back(LaguerrePolynomial((*Z_u).col(i), degree));
            laguerre_X_u.push_back(LaguerrePolynomial((*X_u).col(i), degree));
        }
        
        // LAMBDA3 - LAMBDA4
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
    }

    // Function to estimate the conditional expectation E_i[ sum_l_i+1_N {lamda3_l - lambda4_l} ] given the coefficients of the first regression
    Matrix estimate_conditional_expectation_i() const {
        Matrix cond_exp_i = Matrix::Zero(M, N);
        for (std::size_t i=0; i<N; ++i){
            cond_exp_i.col(i) = regressors_i[i] * coeff_i.col(i);
        }
        return cond_exp_i;
    }
    // Function to estimate the conditional expectation E_i[ lamda1_j - lamda2_j + sum_l_j+1_N {lambda3_l - lambda2_4} ] given the coefficients of the second regression
    CondExp_ij estimate_conditional_expectation_ij() const {
        CondExp_ij cond_exp_ij;
        cond_exp_ij.reserve(N-1); // Reserve space for N-1 matrices
        for (std::size_t i=0; i<N-1; ++i){
            cond_exp_ij.emplace_back(M, N-1-i); // Create a new Matrix for each i
            for (std::size_t j=i; j<N-1; ++j){
                cond_exp_ij[i].col(j-i) = regressors_ij[i] * coeff_ij[i].col(j-i);
            }
        }
        return cond_exp_ij;
    }

    // Getter functions to access coefficients for testing/debugging
    const Matrix& get_coeff_i() const { return coeff_i; }
    const CondExp_ij& get_coeff_ij() const { return coeff_ij; }

    // Getters fpr debugging
    const std::vector<Matrix>& get_laguerre_alpha() const { return laguerre_alpha; }
    const std::vector<Matrix>& get_laguerre_Z_u() const { return laguerre_Z_u; }
    const std::vector<Matrix>& get_laguerre_X_u() const { return laguerre_X_u; }
    const Matrix& get_target_i() const { return target_i; }
    const Matrix& get_target_ij() const { return target_ij; }
    const std::vector<Matrix>& get_regressors_i() const { return regressors_i; }
    const std::vector<Matrix>& get_regressors_ij() const { return regressors_ij; }
    const Matrix& get_diff_lambda_3_4() const { return diff_lambda_3_4; }

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