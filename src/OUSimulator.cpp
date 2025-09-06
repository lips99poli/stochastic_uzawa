#include "OUSimulator.hpp"

// Static helper method implementation
std::uint64_t OUSimulator::make_seed(int user_seed) {
    if (user_seed != 0) {
        return static_cast<std::uint64_t>(user_seed);
    } else {
        return static_cast<std::uint64_t>(std::random_device{}());
    }
}

// Constructor implementation
OUSimulator::OUSimulator(const OUParams& ou_params, const Vector& time_grid, const std::size_t M):
    ou_params(ou_params)
    ,time_grid(time_grid)
    ,time_delta(time_grid.tail(time_grid.size() - 1) - time_grid.head(time_grid.size() - 1)) // Compute time step size
    ,time_integral(construct_time_integral())
    ,actual_seed1(make_seed(ou_params.seed1))
    ,actual_seed2(make_seed(ou_params.seed2))
    ,OU(generateOrnsteinUhlenbeck(M, actual_seed1)) // Generate Ornstein-Uhlenbeck process paths
    ,Pt(compute_Pt())
    ,price(compute_price(actual_seed2))
    ,alpha(Matrix::Zero(M, ou_params.N)) // Initialize alpha to zero matrix
    ,R()
{
    compute_alpha_R(M);
}

// Private method implementations
Matrix OUSimulator::construct_time_integral() {
    // Each row of time_integral contains the weights to write: integral_0^t us ds = sum_{j=0}^{N-1} us * (t_{j+1} - t_j)
    Matrix time_int = Matrix::Zero(ou_params.N, ou_params.N);
    for (std::size_t i = 1; i < ou_params.N; ++i) {
        time_int.row(i).head(i) = time_delta.head(i);
    }
    return time_int;
}

Matrix OUSimulator::generateBrownianMotion(const size_t M, int seed) const{
    // Generate M paths of Brownian motion over the time grid
    Matrix BM(M, ou_params.N); //time_delta.size()=time_grid.size() - 1
    // Precompute sqrt of time increments
    Vector sqrt_dt = time_delta.array().sqrt();
    
    // Parallelize over paths (rows) - each thread handles complete paths independently
    // #pragma omp parallel for schedule(static)
    for (par_for_type m = 0; m < static_cast<par_for_type>(M); ++m) {
        std::mt19937_64 gen(seed + m); // Different seed for each path
        std::normal_distribution<double> Z(0.0, 1.0);

        for (size_t i = 1; i < ou_params.N; ++i) {
            BM(m, i) = BM(m, i - 1) + Z(gen) * sqrt_dt(i - 1);
        }
    }
    return BM;
}

Matrix OUSimulator::generateOrnsteinUhlenbeck(const std::size_t M, int seed) const {
    const Matrix BM = generateBrownianMotion(M, seed); // Generate Brownian motion paths
    // Generate the Ornstein-Uhlenbeck process paths based on the Brownian motion
    Matrix OU(BM.rows(), BM.cols());

    // Precompute A(t)
    std::vector<double> A_t(ou_params.N);
    // #pragma omp parallel for schedule(static)
    for (par_for_type i = 0; i < static_cast<par_for_type>(ou_params.N); ++i) {
        A_t[i] = ou_params.theta * std::sin(ou_params.omega * time_grid(i) + ou_params.phi);
    }

    // Initialize the first column with the initial condition
    OU.col(0).setConstant(ou_params.I0);

    // Version using SDE - each path computed independently to avoid race conditions
    // #pragma omp parallel for schedule(static)
    for (par_for_type m = 0; m < static_cast<par_for_type>(BM.rows()); ++m) { //fix sample_path
        for(size_t i = 1; i < static_cast<size_t>(BM.cols()); ++i){ //fix time instant
            OU(m, i) = OU(m, i - 1) + (A_t[i] - ou_params.k * OU(m, i - 1)) * time_delta(i-1) + ou_params.psi * (BM(m, i) - BM(m, i - 1));
        }
    }
    return OU;
}

Matrix OUSimulator::compute_Pt(){
    return OU * time_integral.transpose();
}

Matrix OUSimulator::compute_price(int seed) const{
    // Compute the price based on the OU process paths
    Matrix initial_price = Matrix::Constant(OU.rows(), OU.cols(), ou_params.S0);
    Matrix BM_price = ou_params.sigma * generateBrownianMotion(OU.rows(),seed); // Generate Brownian motion for pricing
    
    return initial_price + Pt + BM_price;
}

void OUSimulator::compute_alpha_R(const std::size_t M) {
    // Compute the signal matrixes R and alpha based on the OU process

    // Precompute deterministic terms
    const double omega_T_phi = ou_params.omega * ou_params.T + ou_params.phi;
    const double sin_omega_T_phi = std::sin(omega_T_phi);
    const double cos_omega_T_phi = std::cos(omega_T_phi);
    const double theta_over_den = ou_params.theta / (ou_params.k * ou_params.k + ou_params.omega * ou_params.omega);

    std::vector<double> sin_omega_ti_phi(ou_params.N), cos_omega_ti_phi(ou_params.N), exp_neg_k_T_ti(ou_params.N);
    // #pragma omp parallel for schedule(static)
    for (par_for_type i = 0; i < static_cast<par_for_type>(ou_params.N); ++i) {
        const double t_i = time_grid(i);
        sin_omega_ti_phi[i] = std::sin(ou_params.omega * t_i + ou_params.phi);
        cos_omega_ti_phi[i] = std::cos(ou_params.omega * t_i + ou_params.phi);
        exp_neg_k_T_ti[i] = std::exp(-ou_params.k * (ou_params.T - t_i));
    }

    // Precompute coeff first row term
    Matrix coeff_first_term_ij = Matrix::Zero(ou_params.N, ou_params.N);
    // #pragma omp parallel for schedule(static)
    for (par_for_type i = 0; i < static_cast<par_for_type>(ou_params.N); ++i) {
        for (par_for_type j = i; j < static_cast<par_for_type>(ou_params.N); ++j) {
            coeff_first_term_ij(i,j) = (std::exp(ou_params.k * (time_grid(i)-time_grid(j))) - exp_neg_k_T_ti[j]) / ou_params.k;
        }
    }

    // Precompute second row term
    std::vector<double> second_row_term(ou_params.N, 0.0);
    if(ou_params.omega != 0){
        // #pragma omp parallel for schedule(static)
        for (par_for_type i = 0; i < static_cast<par_for_type>(ou_params.N); ++i) {
            second_row_term[i] = theta_over_den * ( 
                ou_params.k/ou_params.omega * (std::cos(omega_T_phi) - cos_omega_ti_phi[i]) + 
                sin_omega_T_phi - sin_omega_ti_phi[i] );
        }
    }

    // Initialize R
    R.clear();
    R.resize(M);
    // Compute R and alpha
    // #pragma omp parallel for schedule(static)
    for (par_for_type m = 0; m < static_cast<par_for_type>(M); ++m) { //fix path
        Matrix R_m = Matrix::Zero(ou_params.N, ou_params.N);

        for (std::size_t i = 0; i < ou_params.N; ++i) { //fix smaller time instant
            const double first_term = OU(m,i) - theta_over_den*(ou_params.k * sin_omega_ti_phi[i] - ou_params.omega * cos_omega_ti_phi[i]);

            for (std::size_t j = i; j < ou_params.N; ++j) { //fix larger time instant
                R_m(i,j) = first_term * coeff_first_term_ij(i,j);
                if (ou_params.omega != 0){
                    R_m(i,j) -= second_row_term[j];
                }
            }
        }
        alpha.row(m) = R_m.diagonal();
        R[m] = std::move(R_m);
    }
}

// Public getter method implementations
Matrix OUSimulator::getPrice() const {
    // Return the price paths
    return price;
}

Matrix OUSimulator::getAlpha() const {
    // Return the alpha matrix
    return alpha;
}

MatVec OUSimulator::getR() const {
    // Return the signal matrices R
    return R;
}
