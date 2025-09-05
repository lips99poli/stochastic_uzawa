#!/bin/bash

# Interactive Experiment App workflow script
# Usage: ./run_app.sh <soft|hard|lib> <parameter_file> [output_folder_name]

# Check for required arguments
if [[ $# -lt 2 ]]; then
    echo "Error: Missing required arguments"
    echo "Usage: ./run_app.sh <soft|hard|lib> <parameter_file> [output_folder_name]"
    echo "  soft|hard|lib: Environment setup option (required)"
    echo "  parameter_file: Path to parameter file (required)"
    echo "  output_folder_name: Output folder name (optional, default: timestamped)"
    exit 1
fi

OPTION="$1"
PARAMETER_FILE="$2"
OUTPUT_FOLDER=${3:-exp_$(date +%Y%m%d_%H%M%S)}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

echo "=== Interactive Experiment App Workflow ==="
echo "Setup option: $OPTION"
echo "Parameter file: $PARAMETER_FILE"
echo "Output folder: $OUTPUT_FOLDER"
echo "Project root: $PROJECT_ROOT"

# Setup Python environment
echo "=== Step 1: Setting up Python environment ==="
"$PROJECT_ROOT/scripts/python_setup.sh" "$OPTION"

# Activate virtual environment and run interactive experiment
echo "=== Step 2: Running Interactive Experiment ==="
cd "$PROJECT_ROOT"
source venv/bin/activate
python apps/interactive_experiment.py "$PARAMETER_FILE" "$OUTPUT_FOLDER"

if [ $? -eq 0 ]; then
    echo "=== Workflow completed successfully! ==="
    echo "Results in: $PROJECT_ROOT/outputs/app/$OUTPUT_FOLDER"
else
    echo "=== Workflow failed! ==="
    exit 1
fi
