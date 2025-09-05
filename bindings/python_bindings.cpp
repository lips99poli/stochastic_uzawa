#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include "Interface.hpp"

namespace py = pybind11;

PYBIND11_MODULE(stochastic_uzawa, m) {
    m.doc() = "Stochastic Uzawa Solver - Python bindings for optimal trading with constraints";
    
    // Bind ParamError struct
    py::class_<ParamError>(m, "ParamError")
        .def(py::init<const std::string&, const std::string&>(),
             "Create a parameter error",
             py::arg("path"), py::arg("message"))
        .def_readwrite("path", &ParamError::path, "Parameter path that caused the error")
        .def_readwrite("message", &ParamError::message, "Error message description")
        .def("__str__", [](const ParamError& e) {
            return "ParamError(path='" + e.path + "', message='" + e.message + "')";
        })
        .def("__repr__", [](const ParamError& e) {
            return "ParamError(path='" + e.path + "', message='" + e.message + "')";
        });
    
    // Bind Interface class
    py::class_<Interface>(m, "Interface")
        .def(py::init<>(), "Create a new Interface instance")
        
        // Parameter management methods
        .def("read_par", &Interface::read_par,
             "Read and validate parameters from a file",
             py::arg("file_path"),
             R"pbdoc(
                Read parameters from a configuration file and validate them.
                
                Args:
                    file_path (str): Path to the parameter file
                    
                Returns:
                    list[ParamError]: List of parameter errors (empty if successful)
                    
                Example:
                    interface = Interface()
                    errors = interface.read_par("data/Parameters.pot")
                    if not errors:
                        print("Parameters loaded successfully!")
                    else:
                        for error in errors:
                            print(f"Error in {error.path}: {error.message}")
             )pbdoc")
        
        .def("update_price_par", &Interface::update_price_par,
             "Update only price-related parameters from a file",
             py::arg("file_path"),
             R"pbdoc(
                Update price-related parameters (T, N, M, OU parameters) from a file.
                Must call read_par() successfully first.
                
                Args:
                    file_path (str): Path to the parameter file
                    
                Returns:
                    list[ParamError]: List of parameter errors (empty if successful)
                    
                Example:
                    errors = interface.update_price_par("new_price_params.pot")
                    if not errors:
                        print("Price parameters updated successfully!")
             )pbdoc")
        
        // Simulation and solving methods
        .def("simulate_price", &Interface::simulate_price,
             "Simulate price paths using current parameters",
             R"pbdoc(
                Generate new price realizations based on current OU parameters.
                Can be called multiple times to try different price scenarios.
                
                Returns:
                    numpy.ndarray: Price matrix of shape (M, N) where M is number of paths and N is time steps
                    
                Example:
                    price1 = interface.simulate_price()
                    price2 = interface.simulate_price()  # Different realization
                    print(f"Price matrix shape: {price1.shape}")
             )pbdoc")
        
        .def("solve", &Interface::solve,
             "Solve the stochastic optimal control problem",
             R"pbdoc(
                Solve the constrained optimization problem using the Uzawa algorithm.
                Uses the most recently simulated price paths.
                
                Note:
                    Must call simulate_price() first to generate price paths.
                    This method resets all variables and solves from initial conditions.
                    
                Example:
                    interface.simulate_price()  # Generate price scenario
                    interface.solve()           # Solve optimization
                    print(f"Converged in {interface.iterations()} iterations")
             )pbdoc")
        
        // Getter methods for results
        .def("get_price", &Interface::get_price,
             "Get the price matrix used in the last solve",
             R"pbdoc(
                Returns:
                    numpy.ndarray: Price matrix of shape (M, N)
             )pbdoc")
        
        .def("get_u", &Interface::get_u,
             "Get the optimal control (trading rate)",
             py::return_value_policy::reference_internal,
             R"pbdoc(
                Returns:
                    numpy.ndarray: Control matrix u of shape (M, N) representing optimal trading rates
             )pbdoc")
        
        .def("get_X", &Interface::get_X,
             "Get the state variable (inventory)",
             py::return_value_policy::reference_internal,
             R"pbdoc(
                Returns:
                    numpy.ndarray: State matrix X of shape (M, N) representing inventory levels
             )pbdoc")
        
        .def("get_lambda1", &Interface::get_lambda1,
             "Get Lagrange multiplier lambda1 (lower bound on u)",
             py::return_value_policy::reference_internal,
             R"pbdoc(
                Returns:
                    numpy.ndarray: Lagrange multiplier matrix lambda1 of shape (M, N)
             )pbdoc")
        
        .def("get_lambda2", &Interface::get_lambda2,
             "Get Lagrange multiplier lambda2 (upper bound on u)",
             py::return_value_policy::reference_internal,
             R"pbdoc(
                Returns:
                    numpy.ndarray: Lagrange multiplier matrix lambda2 of shape (M, N)
             )pbdoc")
        
        .def("get_lambda3", &Interface::get_lambda3,
             "Get Lagrange multiplier lambda3 (lower bound on X)",
             py::return_value_policy::reference_internal,
             R"pbdoc(
                Returns:
                    numpy.ndarray: Lagrange multiplier matrix lambda3 of shape (M, N)
             )pbdoc")
        
        .def("get_lambda4", &Interface::get_lambda4,
             "Get Lagrange multiplier lambda4 (upper bound on X)",
             py::return_value_policy::reference_internal,
             R"pbdoc(
                Returns:
                    numpy.ndarray: Lagrange multiplier matrix lambda4 of shape (M, N)
             )pbdoc")
        
        .def("get_time_grid", &Interface::get_time_grid,
             "Get the time grid used in the simulation",
             py::return_value_policy::reference_internal,
             R"pbdoc(
                Returns:
                    numpy.ndarray: Time grid vector of length N+1 from 0 to T
             )pbdoc")
        
        .def("iterations", &Interface::iterations,
             "Get the number of iterations used in the last solve",
             R"pbdoc(
                Returns:
                    int: Number of Uzawa algorithm iterations
             )pbdoc")
        
        // Python-friendly string representation
        .def("__repr__", [](const Interface&) {
            return "<stochastic_uzawa.Interface>";
        });
    
    // Module-level information
    m.attr("__version__") = "1.0.0";
    m.attr("__author__") = "Filippo Lipari";
}
