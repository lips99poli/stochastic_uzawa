import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


class Plotter:
    def __init__(self, output_dir):
        """Initialize plotter with required output directory."""
        if not output_dir:
            raise ValueError("Output directory is required")
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def plot_variable(self, var_name, time_grid, matrix, iterations=0):
        """Plot a variable and save to output directory."""
        time_values = time_grid[:-1]  # Exclude last value
        matrix = np.array(matrix)
        
        plt.figure(figsize=(10, 6))
        num_paths = min(10, matrix.shape[0])
        for i in range(num_paths):
            plt.plot(time_values, matrix[i, :len(time_values)], label=f'Path {i+1}', marker='o')
        
        plt.xlabel('Time')
        plt.ylabel(var_name)
        plt.title(f'{var_name} - Iterations: {iterations}')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        output_file = self.output_dir / f'{var_name.lower()}_plot.png'
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
    
    def plot_from_file(self, variables_file):
        """Parse variables.txt and plot all variables."""
        data = self._parse_file(variables_file)
        time_grid = data['TIME_GRID'][0]  # First row
        iterations = data.get('iterations', 100)
        
        for var_name in ['PRICE', 'U', 'X', 'LAMBDA1', 'LAMBDA2', 'LAMBDA3', 'LAMBDA4']:
            if var_name in data:
                self.plot_variable(var_name, time_grid, data[var_name], iterations)
    
    def plot_from_matrices(self, time_grid, price, u, x, lambda1, lambda2, lambda3, lambda4, iterations=100):
        """Plot all variables from given matrices."""
        variables = {
            'PRICE': price,
            'U': u,
            'X': x,
            'LAMBDA1': lambda1,
            'LAMBDA2': lambda2,
            'LAMBDA3': lambda3,
            'LAMBDA4': lambda4
        }
        
        for var_name, matrix in variables.items():
            self.plot_variable(var_name, time_grid, matrix, iterations)
    
    def _parse_file(self, variables_file):
        """Parse variables.txt file and return data dict."""
        data = {}
        iterations = 100
        
        with open(variables_file, 'r') as f:
            lines = f.readlines()
        
        current_matrix = None
        matrix_data = []
        
        for line in lines:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith("ITERATIONS="):
                iterations = int(line.split("=")[1])
            elif line.endswith("=") and line[:-1] in ["TIME_GRID", "PRICE", "U", "X", "LAMBDA1", "LAMBDA2", "LAMBDA3", "LAMBDA4"]:
                if current_matrix and matrix_data:
                    data[current_matrix] = np.array(matrix_data)
                current_matrix = line[:-1]
                matrix_data = []
            else:
                try:
                    row = [float(x) for x in line.split()]
                    if row:
                        matrix_data.append(row)
                except ValueError:
                    continue
        
        if current_matrix and matrix_data:
            data[current_matrix] = np.array(matrix_data)
        
        data['iterations'] = iterations
        return data


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python plotter.py <variables_file_path>")
        sys.exit(1)
    
    variables_file = sys.argv[1]
    output_dir = Path(variables_file).parent
    plotter = Plotter(output_dir)
    plotter.plot_from_file(variables_file)
    print(f"Plots saved in {output_dir}")
