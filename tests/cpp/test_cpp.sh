#!/bin/bash

# test_cpp.sh - Complete workflow for C++ test
# Usage: ./test_cpp.sh [path_to_par] [output_folder_name]

set -e  # Exit on any error

# Default values
DEFAULT_PAR_FILE="data/Parameters.pot"
DEFAULT_OUTPUT_FOLDER="exp_$(date +%Y%m%d_%H%M%S)"

# Parse arguments
PAR_FILE="${1:-$DEFAULT_PAR_FILE}"
OUTPUT_FOLDER="${2:-$DEFAULT_OUTPUT_FOLDER}"

echo "=== Stochastic Uzawa C++ Test Workflow ==="
echo "Parameter file: $PAR_FILE"
echo "Output folder: $OUTPUT_FOLDER"

# Get project root (script is in tests/cpp/)
PROJECT_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
BUILD_DIR="$PROJECT_ROOT/tests/cpp/build"

# Make paths absolute
ABS_PAR_FILE="$PROJECT_ROOT/$PAR_FILE"
ABS_OUTPUT_FOLDER="$OUTPUT_FOLDER"  # This will be relative to project root

echo "Project root: $PROJECT_ROOT"

# Step 1: Configure CMake without bindings
echo "=== Step 1: Configuring CMake ==="
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"
cmake "$PROJECT_ROOT" -DBUILD_BINDINGS=OFF

# Step 2: Build executable
echo "=== Step 2: Building executable ==="
make

# Step 3: Run test
echo "=== Step 3: Running test ==="
./test_cpp -i "$ABS_PAR_FILE" -o "$ABS_OUTPUT_FOLDER"

# Step 4: Generate plots
echo "=== Step 4: Generating plots ==="
cd "$PROJECT_ROOT"
source venv/bin/activate
python tools/plotter.py "outputs/$ABS_OUTPUT_FOLDER/variables.txt"

echo "=== Workflow completed successfully! ==="
echo "Results in: $PROJECT_ROOT/outputs/$ABS_OUTPUT_FOLDER"
