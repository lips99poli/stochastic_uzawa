# Apps Directory

This directory contains interactive applications for the Stochastic Uzawa project.

## Interactive Experiment (`interactive_experiment.py`)

A user-guided application that allows interactive experimentation with the Stochastic Uzawa algorithm. The app provides step-by-step control over the simulation process with real-time feedback and parameter adjustment capabilities.

### What it Does

The interactive experiment provides a guided workflow with three main phases:

#### 1. Parameter Validation Phase
- **Loads parameters** from specified file
- **Validates all parameters** using the algorithm's validation system
- **Shows detailed error messages** for any invalid parameters
- **Allows parameter correction**: Pauses for user to fix issues and retry
- **Continues only when all parameters are valid**

#### 2. Price Simulation Phase (Interactive Loop)
- **Simulates price paths** using Ornstein-Uhlenbeck process
- **Generates price plot** for immediate visual feedback
- **Asks for user satisfaction** with the generated price scenario
- **If not satisfied**: Allows modification of price parameters (T, N, M, OU parameters)
- **Re-simulates** with updated parameters until user approves
- **Preserves approved price scenario** for final optimization

#### 3. Final Optimization Phase
- **Solves the complete optimization problem** using approved price paths
- **Generates all result plots** (PRICE, U, X, LAMBDA1-4)
- **Saves all results** to specified output directory

### Key Features

- **Smart error reporting**: Shows exact parameter issues with paths and descriptions
- **Visual feedback**: Immediate plot generation for price assessment
- **Iterative refinement**: Modify parameters until satisfied with simulation
- **User-friendly prompts**: Clear y/n validation with error handling
- **Complete result export**: Both numerical data and visual plots

### Output Structure

Results are saved in `outputs/app/{folder_name}/`:
- Individual PNG plots for each variable (price, u, x, lambda1-4)
- All plots generated with consistent formatting and quality

## Running the Interactive Experiment

### Method 1: Automated Workflow (Recommended)

Use the provided shell script for complete automation:

```bash
# From scripts directory
./run_app.sh <setup_option> <parameter_file> [output_folder]

# Examples:
./run_app.sh soft data/Parameters.pot                    # Timestamped output folder
./run_app.sh soft data/Parameters.pot my_experiment      # Custom output folder
./run_app.sh hard data/Parameters.pot clean_test         # Force clean environment
```

**What `run_app.sh` does:**
1. **Sets up Python environment** using specified option (soft/hard/lib)
2. **Activates virtual environment** automatically
3. **Launches interactive experiment** with specified parameters
4. **Handles all dependencies** and error reporting

**Arguments:**
- `setup_option` (required): Environment setup mode
  - `soft`: Use existing environment if available
  - `hard`: Force recreate environment and reinstall everything  
  - `lib`: Keep environment, reinstall library only
- `parameter_file` (required): Path to parameter file
- `output_folder` (optional): Output folder name, defaults to timestamped

### Method 2: Manual Execution

For direct control when environment is already set up:

```bash
# From project root (requires active venv with stochastic_uzawa installed)
./scripts/python_setup.sh
python apps/interactive_experiment.py <parameter_file> <output_folder>

# Example:
python apps/interactive_experiment.py data/Parameters.pot my_experiment
```

**Required arguments:**
- `parameter_file`: Path to parameter file (must exist)
- `output_folder`: Output folder name (created under `outputs/app/`)

## Interactive Usage Flow

### Typical Session:
1. **Start application** with parameter file and output folder
2. **Fix parameter errors** (if any) - edit file when prompted
3. **Review price simulation** - examine generated price plot
4. **Adjust price parameters** (if needed) - modify T, N, M, or OU parameters
5. **Approve final price scenario** when satisfied
6. **Wait for optimization** to complete
7. **Review all results** in output directory

### User Prompts:
- **Parameter errors**: `"Please modify {file} and press Enter when done..."`
- **Price satisfaction**: `"Are you satisfied with the price simulation? (y/yes/n/no):"`
- **Parameter modification**: `"Please modify the price parameters or N,M,T in the parameter file."`

## Examples

```bash
# Quick start with default parameters
./run_app.sh soft data/Parameters.pot

# Experiment with custom output folder
./run_app.sh soft data/Parameters.pot price_sensitivity_test

# Clean environment setup for troubleshooting
./run_app.sh hard data/Parameters.pot debug_session

# Manual execution (environment ready)
./scripts/python_setup.sh 
python apps/interactive_experiment.py data/custom_params.pot manual_test
```

## Use Cases

- **Parameter exploration**: Test different price scenarios interactively
- **Algorithm validation**: Verify results with known parameter sets
- **Research experiments**: Generate results for specific market conditions
- **Educational purposes**: Understand algorithm behavior through guided interaction
- **Debugging**: Isolate issues with specific parameter combinations

## Requirements

- **Python environment** with stochastic_uzawa library installed
- **Interactive terminal** for user prompts and confirmations
- **Parameter file** in valid format (.pot files)
- **Write permissions** for output directory creation