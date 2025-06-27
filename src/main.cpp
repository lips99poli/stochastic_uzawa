#include <iostream>
#include "NumericSchemeParams.hpp"
#include "Kernel.hpp"
#include "StochasticUzawaSolver.hpp"
#include "OUParams.hpp"

#include <Eigen/Dense>

int main() {
    // Load parameters from file
    NumericSchemeParams params;
    try {
        std::string file_name = "../data/NumericSchemeParams.pot";
        params = load_params(file_name);
        std::cout << "T = " << params.T << "\n"
                << "N = " << params.N << "\n"
                << "D = " << params.D << "\n"
                << "M = " << params.M << "\n"
                << "epsilon = " << params.epsilon << "\n"
                << "delta = " << params.delta << "\n"
                << "beta = " << params.beta << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Failed to load parameters: " << e.what() << "\n";
        return 1;
    }

    //hardcoded kernel: se riesco faccio parsing anche di questo
    Kernel kernel(1.0, 1.0);

    //signal_generator for price
    //qua o carico il modello da qualche parte oppure leggo dei parametri per l'OU
    OUParams ou_params(params.T, params.N);

    // Ora istanzio l'UzawaSolver: per il momento accetta solo un set di parametri per simulare l'OU e sfruttare la forma esplicita di R sul paper


    // test
    StochasticUzawaSolver uzawa_solver(params, kernel, ou_params);
    uzawa_solver.verify();

    return 0;
}
