#!/bin/bash

# Repository cleaner script
# Usage: ./cleaner.sh [option|path]

show_help() {
    echo "Usage: ./cleaner.sh [option|path]"
    echo ""
    echo "Cleaning options:"
    echo "  cpp        - Remove outputs/cpp"
    echo "  python     - Remove outputs/python"
    echo "  build_cpp  - Remove tests/cpp/build"
    echo "  build_py   - Remove all __pycache__ directories"
    echo "  venv       - Remove venv"
    echo "  su_lib     - Uninstall stochastic_uzawa library"
    echo "  outputs    - Remove entire outputs directory"
    echo "  all        - Clean everything except venv and library"
    echo "  reset      - Clean everything including venv and library"
    echo "  [path]     - Remove specific experiment folder"
    echo "  help       - Show this help"
}

OPTION=${1:-help}
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

cd "$PROJECT_ROOT"

case "$OPTION" in
    "cpp")
        echo "Removing outputs/cpp..."
        rm -rf outputs/cpp
        echo "Done."
        ;;
    "python")
        echo "Removing outputs/python..."
        rm -rf outputs/python
        echo "Done."
        ;;
    "build_cpp")
        echo "Removing tests/cpp/build..."
        rm -rf tests/cpp/build
        echo "Done."
        ;;
    "build_py")
        echo "Removing __pycache__ directories in tools/..."
        rm -rf tools/__pycache__/
        echo "Done."
        ;;
    "venv")
        echo "Removing venv..."
        rm -rf venv
        echo "Done."
        ;;
    "su_lib")
        echo "Uninstalling stochastic_uzawa library..."
        if [ -d "venv" ]; then
            source venv/bin/activate
            pip uninstall -y stochastic_uzawa 2>/dev/null || echo "Library not installed"
        else
            echo "No venv found"
        fi
        echo "Done."
        ;;
    "outputs")
        echo "Removing entire outputs directory..."
        rm -rf outputs
        echo "Done."
        ;;
    "all")
        echo "Cleaning everything except venv and library..."
        echo "Removing outputs directory..."
        rm -rf outputs
        echo "Removing tests/cpp/build..."
        rm -rf tests/cpp/build
        echo "Removing __pycache__ directories in tools/..."
        rm -rf tools/__pycache__/
        echo "All cleaning completed."
        ;;
    "reset")
        echo "Full reset - cleaning everything including venv and library..."
        echo "Removing outputs directory..."
        rm -rf outputs
        echo "Removing tests/cpp/build..."
        rm -rf tests/cpp/build
        echo "Removing __pycache__ directories in tools/..."
        rm -rf tools/__pycache__/
        echo "Uninstalling stochastic_uzawa library..."
        if [ -d "venv" ]; then
            source venv/bin/activate
            pip uninstall -y stochastic_uzawa 2>/dev/null || echo "Library not installed"
        else
            echo "No venv found"
        fi
        echo "Removing venv..."
        rm -rf venv
        echo "Full reset completed."
        ;;
    "help")
        show_help
        ;;
    *)
        if [ -d "$OPTION" ]; then
            echo "Removing specific path: $OPTION"
            rm -rf "$OPTION"
            echo "Done."
        else
            echo "Error: Unknown option or path not found: $OPTION"
            echo ""
            show_help
            exit 1
        fi
        ;;
esac
