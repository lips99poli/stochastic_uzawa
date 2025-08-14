#include <iostream>
#include <memory>
#include <vector>
#include <iomanip>
#include <Eigen/Dense>
#include "include/LSMCR.hpp"

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

static void print_matrix(const std::string& title, const Matrix& A){
    std::cout << "\n-- " << title << " (" << A.rows() << "x" << A.cols() << ") --\n" << A << "\n";
}

int main(){
    std::cout << "=== LSMCR Verbose Second Regression (d1=1, d2=2, N=3, M=4) ===\n";
    const std::size_t d1 = 1;
    const std::size_t d2 = 2;
    const std::size_t N  = 3;
    const std::size_t M  = 4;
    try{
        auto alpha = std::make_unique<Matrix>(M,N);
        auto Z_u   = std::make_unique<Matrix>(M,N);
        auto X_u   = std::make_unique<Matrix>(M,N);

        *alpha << 0.1, 0.2, 0.3,
                  0.4, 0.5, 0.6,
                  0.7, 0.8, 0.9,
                  1.0, 1.1, 1.2;
        *Z_u <<   0.2, 0.3, 0.4,
                   0.5, 0.6, 0.7,
                   0.8, 0.9, 1.0,
                   1.1, 1.2, 1.3;
        *X_u <<   0.3, 0.4, 0.5,
                   0.6, 0.7, 0.8,
                   0.9, 1.0, 1.1,
                   1.2, 1.3, 1.4;

        std::vector<std::unique_ptr<Matrix>> lambda(4);
        for(int k=0;k<4;++k) lambda[k] = std::make_unique<Matrix>(M,N+1);
        *lambda[0] << 1.0, 2.0, 4.0, 4.0,
                      1.1, 2.1, 4.1, 4.1,
                      1.2, 2.2, 4.2, 4.2,
                      1.3, 2.3, 4.3, 4.3;
        *lambda[1] << 0.5, 1.0, 1.5, 2.0,
                      0.6, 1.1, 1.6, 2.1,
                      0.7, 1.2, 1.7, 2.2,
                      0.8, 1.3, 1.8, 2.3;
        *lambda[2] << 2.0, 3.0, 4.0, 5.0,
                      2.1, 3.1, 4.1, 5.1,
                      2.2, 3.2, 4.2, 5.2,
                      2.3, 3.3, 4.3, 5.3;
        *lambda[3] << 1.5, 2.5, 3.5, 4.5,
                      1.6, 2.6, 3.6, 4.6,
                      1.7, 2.7, 3.7, 4.7,
                      1.8, 2.8, 3.8, 4.8;

        LSMCR model(d1,d2,N,M, alpha, Z_u, X_u, lambda);
        model.update_regressor();

        // Print Laguerre expansions
        const auto& L_alpha = model.get_laguerre_alpha();
        const auto& L_Zu    = model.get_laguerre_Z_u();
        const auto& L_Xu    = model.get_laguerre_X_u();
        for(std::size_t i=0;i<N;++i){
            print_matrix("Laguerre alpha col "+std::to_string(i), L_alpha[i]);
            print_matrix("Laguerre Z_u   col "+std::to_string(i), L_Zu[i]);
            print_matrix("Laguerre X_u   col "+std::to_string(i), L_Xu[i]);
        }

        // Print regressors
        const auto& Phi_i  = model.get_regressors_i();
        const auto& Phi_ij = model.get_regressors_ij();
        for(std::size_t i=0;i<N;++i){
            print_matrix("Regressor PHI_i for i="+std::to_string(i), Phi_i[i]);
            print_matrix("Regressor PHI_ij for i="+std::to_string(i), Phi_ij[i]);
        }

        // Print targets
        print_matrix("diff_lambda_3_4", model.get_diff_lambda_3_4());
        print_matrix("target_i", model.get_target_i());
        print_matrix("target_ij", model.get_target_ij());

        // Print coefficients
        print_matrix("coeff_i", model.get_coeff_i());
        const auto& cij = model.get_coeff_ij();
        for(std::size_t i=0;i<cij.size();++i){
            print_matrix("coeff_ij["+std::to_string(i)+"]", cij[i]);
        }

        // Print estimations
        auto cond_i  = model.estimate_conditional_expectation_i();
        auto cond_ij = model.estimate_conditional_expectation_ij();
        print_matrix("cond_i", cond_i);
        for(std::size_t i=0;i<cond_ij.size();++i){
            print_matrix("cond_ij["+std::to_string(i)+"]", cond_ij[i]);
        }
        std::cout << "\n=== Verbose run complete ===\n";
    }catch(const std::exception& e){
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
