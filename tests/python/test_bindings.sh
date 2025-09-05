#!/bin/bash

# Python bindings test workflow script
# Usage: ./test_bindings.sh [soft|hard|lib] [parameter_file] [output_folder_name]

OPTION=${1:-soft}
PARAMETER_FILE=${2:-data/Parameters.pot}
OUTPUT_FOLDER=${3:-exp_$(date +%Y%m%d_%H%M%S)}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

echo "=== Python Bindings Test Workflow ==="
echo "Setup option: $OPTION"
echo "Parameter file: $PARAMETER_FILE"
echo "Output folder: $OUTPUT_FOLDER"
echo "Project root: $PROJECT_ROOT"

# Setup Python environment
echo "=== Step 1: Setting up Python environment ==="
"$PROJECT_ROOT/scripts/python_setup.sh" "$OPTION"

# Activate virtual environment and run test
echo "=== Step 2: Running Python bindings test ==="
cd "$PROJECT_ROOT"
source venv/bin/activate
python tests/python/test_bindings.py "$PARAMETER_FILE" "$OUTPUT_FOLDER"

if [ $? -eq 0 ]; then
    echo "=== Workflow completed successfully! ==="
    echo "Results in: $PROJECT_ROOT/outputs/python/$OUTPUT_FOLDER"
else
    echo "=== Workflow failed! ==="
    exit 1
fi