#ifndef CONSTRAINTS_HPP
#define CONSTRAINTS_HPP

#include <vector>
#include <Eigen/Dense>

using Matrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using Ref = Eigen::Ref<Matrix>; // Specify the template argument for Eigen::Ref

using Mat_Vec = std::vector<Matrix>; // Map to store conditional expectations for each time step couple (i,j)


class Constraints{
private:
    const double X0, u_min, u_max, X_u_min, X_u_max; 
    const bool terminal_liquidation, stop_trad_price_lb; 
    const double price_lb;
    
public:
    Constraints (const double X0, const double u_min, const double u_max, const double X_u_min, const double X_u_max, 
        const bool terminal_liquidation, 
        const bool stop_trad_price_lb, const double price_lb) 
        : X0(X0), u_min(u_min), u_max(u_max), X_u_min(X_u_min), X_u_max(X_u_max),
          terminal_liquidation(terminal_liquidation), 
          stop_trad_price_lb(stop_trad_price_lb), price_lb(price_lb) {}

    Mat_Vec compute_constraints(const Matrix& price){
        const std::size_t M = price.rows(); // Number of Monte Carlo paths
        const std::size_t N = price.cols(); // Number of time steps
        // Compute the constraints matrix: each of the M rows corresponds to the constraints path for the same row in the price simulation matrix
        Mat_Vec c;
        c.reserve(4);
        c.emplace_back(Matrix::Constant(M, N, u_min)); // u_min constraint
        c.emplace_back(Matrix::Constant(M, N, u_max)); // u_max constraint
        c.emplace_back(Matrix::Constant(M, N, X_u_min)); // X_u_min constraint
        c.emplace_back(Matrix::Constant(M, N, X_u_max)); // X_u_max constraint
        if(terminal_liquidation){
            c[2].col(N-1) = Vector::Zero(M); // terminal liquidation constraint for X_u_min
            c[3].col(N-1) = Vector::Zero(M); // terminal liquidation constraint for X_u_max
        }
        if(stop_trad_price_lb){
            for(std::size_t m=0; m<M; ++m){
                for(std::size_t i=0; i<N; ++i){
                    if(price(m,i) <= price_lb){
                        c[0].row(m).tail(N-i).setZero(); // u_min = 0 if price <= price_lb
                        c[1].row(m).tail(N-i).setZero(); // u_max = 0 if price <= price_lb
                        break; 
                    }
                }
            }
        }
        return c; 
    }

};

#endif // CONSTRAINTS_HPP