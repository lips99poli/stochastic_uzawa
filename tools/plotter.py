import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path


class Plotter:
    def __init__(self, variables_file_path):
        """Initialize plotter with path to variables.txt file."""
        self.variables_file = Path(variables_file_path)
        self.output_dir = self.variables_file.parent
        self.data = {}
        self.iterations = 0
        self._parse_file()
    
    def _parse_file(self):
        """Parse the variables.txt file."""
        with open(self.variables_file, 'r') as f:
            lines = f.readlines()
        
        current_matrix = None
        matrix_data = []
        
        for line in lines:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith("ITERATIONS="):
                self.iterations = int(line.split("=")[1])
            elif line.endswith("=") and line[:-1] in ["TIME_GRID", "PRICE", "U", "X", "LAMBDA1", "LAMBDA2", "LAMBDA3", "LAMBDA4"]:
                if current_matrix and matrix_data:
                    self.data[current_matrix] = np.array(matrix_data)
                current_matrix = line[:-1]  # Remove the '=' at the end
                matrix_data = []
            else:
                # Parse numeric data
                try:
                    row = [float(x) for x in line.split()]
                    if row:
                        matrix_data.append(row)
                except ValueError:
                    continue
        
        # Add the last matrix
        if current_matrix and matrix_data:
            self.data[current_matrix] = np.array(matrix_data)
    
    def plot_all_variables(self):
        """Generate plots for all variables and save them."""
        time_grid = self.data.get("TIME_GRID", np.array([[]]))[0]  # First row
        
        print("Plotting variables...")
        
        # Plot each variable
        for var_name in ["PRICE", "U", "X", "LAMBDA1", "LAMBDA2", "LAMBDA3", "LAMBDA4"]:
            if var_name in self.data:
                matrix = self.data[var_name]
                # Adjust time grid to match matrix columns
                plot_time_grid = time_grid[:matrix.shape[1]] if len(time_grid) > matrix.shape[1] else time_grid
                self._plot_variable(var_name, matrix, plot_time_grid)
        
        print(f"Plots saved in {self.output_dir}")
    
    def _plot_variable(self, var_name, matrix, time_grid):
        """Plot a single variable matrix."""
        plt.figure(figsize=(10, 6))
        
        # Plot only the first 10 paths (rows) in the matrix
        num_paths = min(10, matrix.shape[0])
        for i in range(num_paths):
            plt.plot(time_grid, matrix[i], label=f'Path {i+1}', marker='o')
        
        plt.xlabel('Time')
        plt.ylabel(var_name)
        plt.title(f'{var_name} - Iterations: {self.iterations}')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Save plot
        output_file = self.output_dir / f'{var_name.lower()}_plot.png'
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python plotter.py <variables.txt_path>")
        sys.exit(1)
    
    plotter = Plotter(sys.argv[1])
    plotter.plot_all_variables()
