# Scripts Directory

This directory contains automation scripts for the Stochastic Uzawa project. These scripts handle environment setup, application execution, and project maintenance.

## Scripts Overview

### 🚀 **run_app.sh** - Interactive Experiment Launcher
Launches the interactive experiment application with automatic environment setup.

**Usage:** `./run_app.sh [setup_option] [parameter_file] [output_folder]`

**Parameters:**
- `setup_option` (default: `soft`): Environment setup mode
  - `soft`: Use existing environment if available
  - `hard`: Force recreate environment and reinstall everything
  - `lib`: Keep environment, reinstall library only
- `parameter_file` (default: `data/Parameters.pot`): Path to parameter file
- `output_folder` (default: `exp_YYYYMMDD_HHMMSS`): Output folder name (created under `outputs/app/`)

**Examples:**
```bash
./run_app.sh                                    # Default: soft setup, default params, timestamped folder
./run_app.sh soft data/Parameters.pot my_exp   # Custom experiment name
./run_app.sh hard                               # Force clean setup
```

---

### 🐍 **python_setup.sh** - Python Environment Manager
Sets up the complete Python environment with virtual environment and stochastic_uzawa library.

**Usage:** `./python_setup.sh [option]`

**Options:**
- `soft` (default): Skip if venv and module already exist and work
- `hard`: Recreate venv and reinstall module from scratch
- `lib`: Keep existing venv, reinstall module only

**What it does:**
1. Creates/manages virtual environment via `create_venv.sh`
2. Activates the environment
3. Installs/reinstalls the stochastic_uzawa Python module

---

### 📦 **create_venv.sh** - Virtual Environment Creator
Creates a Python virtual environment with required dependencies.

**Usage:** `./create_venv.sh [option]`

**Options:**
- `soft` (default): Skip if venv already exists
- `hard`: Remove existing venv and create new one

**Dependencies installed:**
- `numpy==1.24.3`
- `matplotlib`
- `pybind11`
- `scikit-build-core`

---

### 🧹 **cleaner.sh** - Project Maintenance Tool
Comprehensive cleaning tool for managing project artifacts and builds.

**Usage:** `./cleaner.sh [option|path]`

**Cleaning Options:**
- `cpp`: Remove `outputs/cpp` directory
- `python`: Remove `outputs/python` directory  
- `build_cpp`: Remove `tests/cpp/build` directory
- `build_py`: Remove all `__pycache__` directories (excluding venv)
- `params`: Remove all files in `data/` except `Parameters.pot` and `README.md`
- `venv`: Remove virtual environment
- `su_lib`: Uninstall stochastic_uzawa library from venv
- `outputs`: Remove entire `outputs` directory
- `all`: Clean everything except venv and library
- `reset`: Clean everything including venv and library
- `[path]`: Remove specific folder or file path
- `help`: Show help message

**Examples:**
```bash
./cleaner.sh help           # Show all options
./cleaner.sh cpp            # Clean C++ outputs only
./cleaner.sh all            # Clean everything but keep environment
./cleaner.sh reset          # Full reset - clean everything
./cleaner.sh outputs/cpp/exp_20250905_123456  # Remove specific experiment
```

---

## Workflow Examples

### Quick Start (First Time Setup)
```bash
# Setup environment and run interactive experiment
./run_app.sh hard

# Or step by step:
./python_setup.sh hard
python ../apps/interactive_experiment.py data/Parameters.pot my_experiment
```

### Development Workflow
```bash
# Clean previous results and run new experiment
./cleaner.sh outputs
./run_app.sh soft data/Parameters.pot new_experiment

# After code changes, reinstall library and test
./python_setup.sh lib
./run_app.sh soft
```

### Maintenance
```bash
# Clean specific outputs
./cleaner.sh cpp              # Remove C++ test results
./cleaner.sh python           # Remove Python test results

# Clean builds (after compilation issues)
./cleaner.sh build_cpp
./cleaner.sh build_py

# Full cleanup and restart
./cleaner.sh reset
./run_app.sh hard
```

---

## Script Dependencies

- **run_app.sh** → **python_setup.sh** → **create_venv.sh**
- All scripts are designed to be run from the scripts directory
- Scripts automatically detect project root directory
- Virtual environment activation is handled automatically where needed

---

## Notes

- All scripts include error handling and status reporting
- Default parameters are designed for typical usage scenarios
- Scripts create necessary directories automatically
- The `soft` option is generally recommended for repeated usage
- Use `hard` options when troubleshooting environment issues