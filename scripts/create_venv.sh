#!/bin/bash

# Create Python virtual environment for Stochastic Uzawa project
# Usage: ./create_venv.sh [soft|hard]

OPTION=${1:-soft}
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VENV_PATH="$PROJECT_ROOT/venv"

if [ -d "$VENV_PATH" ]; then
    if [ "$OPTION" = "soft" ]; then
        echo "Virtual environment already exists. Use 'hard' to recreate."
        exit 0
    elif [ "$OPTION" = "hard" ]; then
        rm -rf "$VENV_PATH"
    else
        echo "Invalid option. Use 'soft' or 'hard'."
        exit 1
    fi
fi

cd "$PROJECT_ROOT"
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install numpy==1.24.3 matplotlib pybind11 scikit-build-core

echo "Virtual environment created. Activate with: source venv/bin/activate"