#!/usr/bin/env python3

import sys
import os
from pathlib import Path
import numpy as np
import stochastic_uzawa
sys.path.append(str(Path(__file__).parents[2] / "tools"))
from plotter import Plotter

def test_bindings(parameter_path, output_folder_name):
    """Test stochastic_uzawa Python bindings."""
    
    # Setup paths
    project_root = Path(__file__).parents[2]
    output_dir = project_root / "outputs" / "python" / output_folder_name
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Parameter file: {parameter_path}")
    print(f"Output directory: {output_dir}")
    
    # Create interface and read parameters
    interface = stochastic_uzawa.Interface()
    try:
        interface.read_par(str(parameter_path))
        print("Parameters read successfully")
    except Exception as e:
        print(f"Error reading parameters: {e}")
        return False
    
    # Simulate and solve
    print("Simulating price...")
    interface.simulate_price()
    
    print("Solving...")
    interface.solve()
    
    # Get all matrices
    print("Retrieving results...")
    time_grid = interface.get_time_grid()
    price = interface.get_price()
    u = interface.get_u()
    x = interface.get_X()
    lambda1 = interface.get_lambda1()
    lambda2 = interface.get_lambda2()
    lambda3 = interface.get_lambda3()
    lambda4 = interface.get_lambda4()
    it = interface.get_iterations()
    
    # Copy parameter file to output directory
    import shutil
    shutil.copy2(parameter_path, output_dir / "Parameters.txt")
    print(f"Parameters copied to: {output_dir / 'Parameters.txt'}")
    
    # Generate plots directly from matrices
    print("Generating plots...")
    plotter = Plotter(output_dir)
    plotter.plot_from_matrices(
        time_grid=time_grid,
        price=price,
        u=u,
        x=x,
        lambda1=lambda1,
        lambda2=lambda2,
        lambda3=lambda3,
        lambda4=lambda4,
        iterations=it
    )
    
    print(f"Plots saved in: {output_dir}")
    print("Test completed successfully!")
    return True

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python test_bindings.py <parameter_path> <output_folder_name>")
        sys.exit(1)
    
    parameter_path = sys.argv[1]
    output_folder_name = sys.argv[2]
    
    success = test_bindings(parameter_path, output_folder_name)
    sys.exit(0 if success else 1)