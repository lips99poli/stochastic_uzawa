# Data Directory

This directory contains data files and configurations for the Stochastic Uzawa project.

## File Organization

### Template vs Working Files

The data directory follows a template/working copy structure for safety:

- **Template Files** (`*_template.*`): Read-only reference configurations
  - `Parameters_template.pot`: Original parameter configuration (read-only)
  - Use as reference for parameter meanings and default values
  - Should NOT be modified directly

- **Working Files**: User-modifiable copies for experimentation
  - `Parameters.pot`: Working copy for experiments and modifications
  - This is what scripts will use by default
  - Safe to modify without losing reference configuration

### Usage Guidelines

1. **To experiment with parameters**: Modify `Parameters.pot` directly
2. **To reset to defaults**: Copy from template: `cp Parameters_template.pot Parameters.pot`  
3. **To understand parameters**: Check both the template file and the documentation below

## Parameter File Format

The `.pot` files use GetPot format with key-value pairs:

```
# Comments start with #
parameter_name = value
section/nested_parameter = value
```

## Key Parameters

### Simulation Parameters
- `T`: Final time for simulation
- `n_timesteps`: Number of time steps
- `n_realizations`: Number of Monte Carlo realizations

### Spatial Discretization
- `n_dofs`: Number of degrees of freedom
- Grid and mesh parameters for spatial domain

### Algorithmic Parameters
- `lambda`: Uzawa penalty parameter
- `rho`: Stochastic step size parameter
- Convergence tolerances

### Output Control
- File output settings
- Verbosity levels
- Plotting options

## Safety Features

- Template file is read-only (permissions: 444)
- Working file can be safely modified
- Easy restoration from template
- Clear separation between reference and experimental configurations

## Integration

All scripts in `/scripts/` directory use `Parameters.pot` (working copy) by default:
- `run_app.sh`: Interactive experiment application
- `test_cpp.sh`: C++ test executable  
- `test_bindings.sh`: Python binding tests

This ensures consistent behavior while protecting the reference configuration.
