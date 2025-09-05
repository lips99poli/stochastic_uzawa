# Tools Directory

This directory contains utility tools for the Stochastic Uzawa project.

## Plotter (`plotter.py`)

A versatile plotting utility for visualizing Stochastic Uzawa algorithm results. The plotter can be used both as a standalone Python script and as an importable module, making it suitable for both C++ and Python workflows.

### What it Does

The plotter generates time-series plots for all algorithm variables:
- **PRICE** - Price signal paths
- **U** - Control variable (trading strategy)
- **X** - State variable (position/inventory)
- **LAMBDA1-LAMBDA4** - Lagrange multipliers for constraints

### Key Features

- **Automatic path limiting**: Plots up to 10 paths maximum for readability
- **Smart time grid handling**: Excludes last time point for proper visualization
- **Consistent formatting**: All plots use same style with grid, legend, and markers
- **High-quality output**: 150 DPI PNG files with tight layout
- **Iteration tracking**: Shows algorithm iteration count in plot titles

### Dual Usage

#### 1. As Standalone Script (C++ Workflow)

Used by C++ tests that export data to `variables.txt` files:

```bash
# From project root (requires active venv)
python tools/plotter.py path/to/variables.txt

# Example from test_cpp.sh
python tools/plotter.py "outputs/cpp/experiment_folder/variables.txt"
```

**Input:** `variables.txt` file with format:
```
ITERATIONS=150
TIME_GRID=
0.0 0.1 0.2 ... 1.0
PRICE=
path1_values...
path2_values...
U=
control_values...
...
```

**Output:** Individual PNG files for each variable in the same directory as input file.

#### 2. As Importable Module (Python Workflow)

Used by Python scripts that have direct access to matrices:

```python
from tools.plotter import Plotter

# Initialize with output directory
plotter = Plotter("outputs/python/my_experiment")

# Plot from matrices (Python bindings)
plotter.plot_from_matrices(
    time_grid=time_grid,
    price=price_matrix,
    u=u_matrix,
    x=x_matrix,
    lambda1=lambda1_matrix,
    lambda2=lambda2_matrix,
    lambda3=lambda3_matrix,
    lambda4=lambda4_matrix,
    iterations=iteration_count
)

# Or plot individual variables
plotter.plot_variable("PRICE", time_grid, price_matrix, iterations)
```

### Manual Usage Examples

```bash
# Plot results from C++ test
python tools/plotter.py outputs/cpp/exp_20250905_123456/variables.txt

# Plot results from specific experiment
python tools/plotter.py outputs/cpp/my_experiment/variables.txt

# From different directory (use full path)
python /path/to/project/tools/plotter.py /path/to/variables.txt
```

### Integration in Workflows

- **C++ Tests**: `test_cpp.sh` automatically calls plotter on generated `variables.txt`
- **Python Tests**: `test_bindings.py` uses plotter as module with direct matrix input
- **Interactive App**: `interactive_experiment.py` uses plotter for both price preview and final results

### Output Structure

All plots are saved as PNG files in the same directory as the input/specified output:
- `price_plot.png`
- `u_plot.png`
- `x_plot.png`
- `lambda1_plot.png`
- `lambda2_plot.png`
- `lambda3_plot.png`
- `lambda4_plot.png`

### Requirements

- Python environment with `numpy`, `matplotlib`, and `pathlib`
- For standalone usage: Virtual environment activation
- For module usage: Import path setup (handled automatically in project structure)