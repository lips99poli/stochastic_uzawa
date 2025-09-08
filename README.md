# Stochastic Uzawa — Overview and Quick Start

This project implements a **stochastic Uzawa** method to simulate and solve constrained stochastic optimal control problems—motivated by **optimal execution** and **battery storage**. The core engine is **C++20 (Eigen3)**; plotting and the interactive workflow are in **Python** via a **pybind11** module.

**Background (paper).** The method follows the formulation in *Trading with propagators and constraints: applications to optimal execution and battery storage* (Abi Jaber, De Carvalho, Pham): a stochastic convex program with pathwise inequality constraints, **KKT** optimality conditions leading to a **stochastic Fredholm equation of the second kind**, and a numerical resolution via **stochastic Uzawa** (projected ascent on multipliers) with **LSMC** to approximate conditional expectations—exactly the workflow reproduced here.

**Paper reference:** arXiv:2409.12098v1 [math.OC], 18 Sep 2024. See `doc/` for the PDF.

## Branch Structure

This repository implements two optimization strategies for performance comparison:

- **`main`**: Stable branch with OpenMP parallelization (production-ready)
- **`openmp-parallel`**: Development branch for OpenMP optimization features
- **`eigen-optimal`**: Development branch for Eigen+BLAS/LAPACK optimization

The dual-branch architecture enables systematic comparison of parallelization strategies while maintaining clean, optimized implementations of each approach.

## Prerequisites

- **C++**: C++20 compiler, CMake 3.16+, Eigen3, OpenMP (for `main`/`openmp-parallel`) or BLAS/LAPACK (for `eigen-optimal`)
- **Python**: Python 3.8+ (for bindings and interactive features)

## Quick Start

```bash
# Interactive experiment (recommended)
./scripts/run_app.sh soft data/Parameters.pot

# C++ test (fast native performance)
./tests/cpp/test_cpp.sh data/Parameters.pot

# Python bindings test
./tests/python/test_bindings.sh soft data/Parameters.pot
```

---

## Project Structure

### Core Implementation
- **`/src/`, `/include/`**: C++20 implementation with Eigen3 matrix operations and OpenMP parallelization
- **`/bindings/`**: pybind11 Python bindings for seamless Python integration
- **`CMakeLists.txt`**: Branch-specific CMake configuration (OpenMP vs Eigen+BLAS/LAPACK optimization)
- **`pyproject.toml`**: Modern Python build; `pip install -e .` invokes CMake as needed

### Applications & Testing
- **`/apps/`**: Interactive applications (`interactive_experiment.py`)
- **`/tests/cpp/`**: C++ test executable - validates core algorithm, exports matrices
- **`/tests/python/`**: Python binding tests - same validation with automatic plotting

### Utilities & Data
- **`/tools/`**: Utilities (`plotter.py` - dual-purpose plotting tool)
- **`/data/`**: Configuration files with template/working copy structure
  - `Parameters.pot`: Working parameter file (read-write)
  - `Parameters_template.pot`: Reference configuration (read-only)
- **`/doc/`**: Paper PDF (arXiv reference) and project reports

### Automation
- **`/scripts/`**: Complete automation scripts for all workflows
  - `run_app.sh`: Launch interactive experiment
  - `test_cpp.sh`: Execute C++ tests with optional plotting
  - `test_bindings.sh`: Execute Python binding tests
  - `python_setup.sh`: Environment setup and package installation
  - `create_venv.sh`: Virtual environment management
  - `cleaner.sh`: Project cleanup utilities

---

## Detailed Usage

All scripts support three setup modes and accept parameters for custom output folders. For complete details, see each folder's README.

### Setup Modes
- **`soft`** — reuse venv/ if present; build/install only if needed
- **`hard`** — recreate/refresh environment and rebuild the module  
- **`lib`** — rebuild the extension without touching the rest

### 1) Interactive Application (Recommended)

Terminal-driven workflow: parameter validation → 10-path price preview → optional solve → save plots.

```bash
./scripts/run_app.sh <soft|hard|lib> <parameter_file> [output_folder_name]

# Example with custom output folder:
./scripts/run_app.sh soft data/Parameters.pot exp_$(date +%Y%m%d_%H%M%S)
```

### 2) C++ Native Test (Fast Performance)

Builds C++ executable, runs algorithm, exports matrices, generates plots via Python plotter.

```bash
./tests/cpp/test_cpp.sh <parameter_file> [output_folder_name]

# Example:
./tests/cpp/test_cpp.sh data/Parameters.pot my_test_results
```

### 3) Python Bindings Test

Validates Python module imports and interface functionality with automatic plotting.

```bash
./tests/python/test_bindings.sh <soft|hard|lib> <parameter_file> [output_folder_name]

# Example:
./tests/python/test_bindings.sh soft data/Parameters.pot binding_test
```

## Quick usage (via `.sh` scripts)

The scripts handle flags, environment, and paths. For full details, see each folder’s README.

### 1) Interactive application (Python + bindings + plots)

Terminal-driven workflow: parameter validation → 10-path price preview (allows price-only edits) → optional solve → save plots and parameter snapshot.

```bash
./scripts/run_app.sh <soft|hard|lib> <parameter_file> [output_folder_name]
# soft — reuse venv/ if present; build/install only if needed
# hard — recreate/refresh environment and rebuild the module
# lib  — rebuild the extension without touching the rest

# Example:
./scripts/run_app.sh soft data/Parameters.pot exp_$(date +%Y%m%d_%H%M%S)

# Results:
# outputs/app/<output_folder>/
```

### 2) C++-only end-to-end test (fast native loop) 
Builds and runs the native C++ test, writes `.txt` matrices, then uses the Python plotter to generate figures. 
```bash 
./tests/cpp/test_cpp.sh <parameter_file> [output_folder_name]
```
What it does: 
* Configures root CMake with bindings disabled (`-DBUILD_BINDINGS=OFF`), 
* Builds the C++ test executable with parallel compilation (`make -j`),
* Runs: `./test_cpp -i <parameter_file> -o <output_folder>`, 
* Activates `venv/` and calls `python tools/plotter.py "outputs/cpp/<folder>/variables.txt"`, 
* Results in: `outputs/cpp/<output_folder>/`.


### 3) Python bindings test 
Ensures the module imports and validates algorithm functionality with automatic plotting.
```bash
./tests/python/test_bindings.sh <soft|hard|lib> <parameter_file> [output_folder_name]
```
Internally ensures `venv/` exists and `stochastic_uzawa` module is importable; rebuilds only if needed. 
Saves results under `outputs/python/<output_folder>/`.


---

## Output Structure

All results are organized under the `outputs/` directory:
- **`outputs/app/`**: Interactive experiment results (full workflow outputs)
- **`outputs/cpp/`**: C++ test results (matrices in text format)
- **`outputs/python/`**: Python test results (plots + parameters)

## Key Features

- **Dual optimization strategies** with systematic performance comparison (OpenMP vs Eigen+BLAS/LAPACK)
- **Dual C++/Python workflows** with consistent interfaces
- **Interactive parameter adjustment** with real-time feedback  
- **Automated environment management** (virtual environments, dependencies)
- **Safe parameter experimentation** (template/working file structure)
- **Comprehensive test coverage** with automatic plotting and performance analysis
- **Modern packaging** (CMake + scikit-build-core with parallel compilation)
- **Branch-specific optimizations** for research and production use

## Getting Help

For detailed information on each component, see the README files in individual directories:
- `/scripts/README.md`: Complete script documentation
- `/apps/README.md`: Interactive application details
- `/tests/cpp/README.md` & `/tests/python/README.md`: Test suite documentation
- `/data/README.md`: Parameter file structure and usage
- `/tools/README.md`: Plotting utilities