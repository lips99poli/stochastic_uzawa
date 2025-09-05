#!/bin/bash

# Setup complete Python environment for Stochastic Uzawa project
# Usage: ./python_setup.sh [soft|hard|lib]
# soft: skip if venv and module exist
# hard: recreate venv and reinstall module  
# lib: keep venv, reinstall module only

OPTION=${1:-soft}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Handle venv based on option
if [ "$OPTION" = "lib" ]; then
    # Keep existing venv, just check it exists
    if [ ! -d "$PROJECT_ROOT/venv" ]; then
        echo "No venv found. Creating one..."
        "$SCRIPT_DIR/create_venv.sh" soft
    else
        echo "Using existing venv."
    fi
else
    # Create/setup virtual environment for soft/hard
    "$SCRIPT_DIR/create_venv.sh" "$OPTION"
fi

# Activate virtual environment
source "$PROJECT_ROOT/venv/bin/activate"

# Handle module installation based on option
if [ "$OPTION" = "soft" ]; then
    # Check if module works, skip if yes
    if python -c "import stochastic_uzawa" 2>/dev/null; then
        echo "Module stochastic_uzawa already installed and working."
        echo "Setup complete."
        exit 0
    else
        echo "Installing module..."
        pip install -e .
    fi
else
    # For hard or lib: always reinstall module
    echo "Reinstalling module..."
    pip uninstall -y stochastic_uzawa 2>/dev/null || true
    pip install -e .
fi

echo "Python setup complete. Module ready to use."