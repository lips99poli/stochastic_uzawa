# Python Bindings Test Suite

This directory contains the Python test script that validates the stochastic_uzawa Python bindings and core algorithm functionality.

## What the Test Does

The `test_bindings.py` script:

1. **Validates parameters** from a given input file using Python bindings
2. **Simulates price signals** using the Ornstein-Uhlenbeck process
3. **Solves the optimization problem** using the Uzawa algorithm
4. **Generates plots** directly from computed matrices
5. **Exports all results** for analysis

### Test Output

Results are saved in `outputs/python/{folder_name}/`:
- `Parameters.txt` - Copy of input parameters used
- Plot files for each variable (price, u, X, λ1-λ4) - PNG format
- All plots generated automatically during test execution

## Running the Test

### Method 1: Automated Workflow (Recommended)

Use the provided shell script for complete workflow including environment setup:

```bash
# Required: setup option and parameter file, optional: output folder name
./test_bindings.sh soft data/Parameters.pot

# With custom output folder
./test_bindings.sh soft data/Parameters.pot my_experiment

# Force clean environment setup
./test_bindings.sh hard data/Parameters.pot validation_test
```

**What the script does:**
1. Sets up Python environment using specified option
2. Activates virtual environment
3. Runs the Python test with specified parameters
4. Handles all dependencies and plotting automatically

**Arguments:**
- `setup_option` (required): Environment setup mode
  - `soft`: Use existing environment if available
  - `hard`: Force recreate environment and reinstall everything
  - `lib`: Keep environment, reinstall library only
- `parameter_file` (required): Path to parameter file
- `output_folder_name` (optional): Output folder name, defaults to timestamped folder

### Method 2: Direct Python Execution

For direct control when environment is already set up:

```bash
# From project root (requires active venv with stochastic_uzawa installed)
source venv/bin/activate
python tests/python/test_bindings.py data/Parameters.pot my_experiment
```

**Required arguments:**
- `parameter_path`: Path to parameter file
- `output_folder_name`: Output folder name (created under `outputs/python/`)

## Examples

```bash
# Test with default parameter file and timestamped output
./test_bindings.sh soft data/Parameters.pot

# Test with custom output folder
./test_bindings.sh soft data/Parameters.pot validation_test

# Test with custom parameter file
./test_bindings.sh soft path/to/custom_params.pot experiment

# Force clean setup (useful for troubleshooting)
./test_bindings.sh hard data/Parameters.pot clean_test

# Direct execution (environment must be ready)
./scripts/python_setup.sh 
python tests/python/test_bindings.py data/Parameters.pot direct_test
```

## Notes

- The script automatically copies the parameter file to the output directory
- All plots are generated automatically using the integrated plotter
- The shell script handles all environment setup and dependency management
- Results include both numerical data and visual plots for analysis
- Environment setup options provide flexibility for different usage scenarios