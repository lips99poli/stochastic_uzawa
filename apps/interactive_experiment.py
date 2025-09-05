#!/usr/bin/env python3
"""
Interactive Experiment Script
Usage: python interactive_experiment.py <parameter_file> <output_folder>
"""

import sys
import os
from pathlib import Path
import stochastic_uzawa

# Add project root to path for imports
sys.path.append(str(Path(__file__).parent.parent))
from tools.plotter import Plotter


def get_user_confirmation(prompt="Are you satisfied? (y/yes/n/no): "):
    """Get user confirmation with proper validation"""
    while True:
        response = input(prompt).strip().lower()
        if response in ['y', 'yes']:
            return True
        elif response in ['n', 'no']:
            return False
        else:
            print("Please answer with 'y', 'yes', 'n', or 'no'")


def wait_for_user_modification(param_file):
    """Wait for user to modify parameter file"""
    input(f"Please modify {param_file} and press Enter when done...")


def main():
    # Check command line arguments
    if len(sys.argv) != 3:
        print("Usage: python interactive_experiment.py <parameter_file> <output_folder>")
        sys.exit(1)
    
    param_file = sys.argv[1]
    output_folder = sys.argv[2]
    
    # Validate parameter file exists
    if not os.path.exists(param_file):
        print(f"Error: Parameter file '{param_file}' not found")
        sys.exit(1)
    
    # Create output directory
    output_dir = Path("outputs/app") / output_folder
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Interactive Experiment")
    print(f"Parameter file: {param_file}")
    print(f"Output directory: {output_dir}")
    print("-" * 50)
    
    # Initialize Interface
    interface = stochastic_uzawa.Interface()
    
    # Parameter reading loop - continue until no errors
    while True:
        print(f"Reading parameters from {param_file}...")
        errors = interface.read_par(param_file)
        
        if not errors:  # Empty list means success
            print("Parameters loaded successfully!")
            break
        else:
            print("Errors found in parameters:")
            for error in errors:
                print(f"  - {error.path}: {error.message}")
            wait_for_user_modification(param_file)
    
    # Initialize plotter
    plotter = Plotter(str(output_dir))
    
    # Price simulation loop - continue until user is satisfied
    while True:
        print("Simulating price signal...")
        price_matrix = interface.simulate_price()
        
        print("Plotting price simulation...")
        time_grid = interface.get_time_grid()
        plotter.plot_variable("PRICE", time_grid, price_matrix)
        print(f"Price plot saved to {output_dir}/price_plot.png")
        
        if get_user_confirmation("Are you satisfied with the price simulation? (y/yes/n/no): "):
            break
        else:
            print("Please modify the price parameters or N,M,T in the parameter file.")
            wait_for_user_modification(param_file)
            
            # Re-read parameters after modification
            print(f"Re-reading parameters from {param_file}...")
            errors = interface.read_par(param_file)
            if errors:
                print("Warning: Errors found in modified parameters:")
                for error in errors:
                    print(f"  - {error.path}: {error.message}")
                print("Continuing with previous parameter values...")
    
    # Final solve and complete plotting
    print("Running final solve...")
    interface.solve()
    
    print("Generating all plots...")
    time_grid = interface.get_time_grid()
    price = interface.get_price()
    u = interface.get_u()
    x = interface.get_X()
    lambda1 = interface.get_lambda1()
    lambda2 = interface.get_lambda2()
    lambda3 = interface.get_lambda3()
    lambda4 = interface.get_lambda4()
    iterations = interface.get_iterations()
    
    plotter.plot_from_matrices(
        time_grid=time_grid,
        price=price,
        u=u,
        x=x,
        lambda1=lambda1,
        lambda2=lambda2,
        lambda3=lambda3,
        lambda4=lambda4,
        iterations=iterations
    )
    print(f"Plots saved in: {output_dir}")
    print("Test completed successfully!")


if __name__ == "__main__":
    main()
