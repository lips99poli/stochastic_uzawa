# C++ Test

This directory contains the C++ test executable that validates the core Stochastic Uzawa algorithm implementation.

## What the Test Does

The `test_cpp.cpp` program:

1. **Validates parameters** from a given input file
2. **Simulates price signals** using Ornstein-Uhlenbeck process
3. **Solves the optimization problem** using the Uzawa algorithm
4. **Exports all results** to text files for analysis

### Test Output

Results are saved in `outputs/cpp/{folder_name}/`:
- `Parameters.txt` - Copy of input parameters used
- `variables.txt` - All computed matrices (time grid, price, u, X, λ1-λ4) and iteration count
- Plot files (when using shell script with plotting)

## Running the Test

### Method 1: Automated Workflow (Recommended)

Use the provided shell script for complete workflow including plotting:

```bash
# Required: parameter file, optional: output folder name
./test_cpp.sh data/Parameters.pot

# With custom output folder
./test_cpp.sh data/Parameters.pot my_experiment

# Using different parameter file
./test_cpp.sh path/to/custom_params.pot validation_test
```

**What the script does:**
1. Configures CMake with `BUILD_BINDINGS=OFF`
2. Builds the executable in `build/` directory
3. Runs the test with specified parameters
4. Generates plots using Python plotter

**Arguments:**
- `parameter_file` (required): Path to parameter file
- `output_folder_name` (optional): Output folder name, defaults to timestamped folder

### Method 2: Manual CMake Build

For direct control over the build process:

```bash
# From project root
mkdir -p tests/cpp/build
cd tests/cpp/build

# Configure (without Python bindings)
cmake ../../.. -DBUILD_BINDINGS=OFF

# Build
make

# Run executable (requires both input and output arguments)
./test_cpp -i ../../../data/Parameters.pot -o my_experiment_folder
```

**Required arguments:**
- `-i, --input`: Path to parameter file
- `-o, --output`: Output folder name (created under `outputs/cpp/`)

### Plotting Results

After running the test, generate plots:

```bash
# From project root (requires active venv)
source venv/bin/activate
python tools/plotter.py outputs/cpp/{my_experiment_folder}/variables.txt
```

## Examples

```bash
# Test with default parameter file
./test_cpp.sh data/Parameters.pot

# Test with custom output folder
./test_cpp.sh data/Parameters.pot validation_test

# Test with custom parameter file
./test_cpp.sh path/to/custom_params.pot experiment

# Manual build and test
cd build
cmake ../../.. -DBUILD_BINDINGS=OFF
make
./test_cpp -i ../../../data/Parameters.pot -o manual_test
```

## Notes

- The executable expects to be run from the `build/` directory for proper path resolution
- Both input parameter file and output folder name are mandatory
- Results include full matrix data suitable for analysis and plotting
- The shell script handles environment setup and plotting automatically