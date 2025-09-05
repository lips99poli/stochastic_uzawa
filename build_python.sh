#!/bin/bash
# Build script for Stochastic Uzawa Solver Python bindings

set -e  # Exit on any error

echo "🔨 Building Stochastic Uzawa Solver Python Bindings"
echo "=================================================="

# Check if we're in the right directory
if [ ! -f "setup.py" ]; then
    echo "❌ Error: setup.py not found. Please run from the project root directory."
    exit 1
fi

# Check if pybind11 is available
echo "🔍 Checking for pybind11..."
python3 -c "import pybind11" 2>/dev/null || {
    echo "❌ pybind11 not found. Installing..."
    pip3 install pybind11 numpy
}

# Method 1: Try building with setup.py (recommended)
echo -e "\n🚀 Building Python module with setup.py..."
python3 setup.py build_ext --inplace

if [ $? -eq 0 ]; then
    echo "✅ Python module built successfully!"
    
    # Test if the module can be imported
    echo -e "\n🧪 Testing module import..."
    python3 -c "import stochastic_uzawa; print('✅ Module import successful!')" || {
        echo "⚠️  Module built but import failed. Checking for potential issues..."
        
        # Try to find the built module
        echo "🔍 Looking for built module files..."
        find . -name "stochastic_uzawa*.so" -o -name "stochastic_uzawa*.pyd" | head -5
    }
    
    # Run the test script if it exists and module imports successfully
    if [ -f "test_python_bindings.py" ]; then
        echo -e "\n🧪 Running comprehensive tests..."
        python3 test_python_bindings.py || {
            echo "⚠️  Tests failed, but module was built. Check test output above."
        }
    fi
    
else
    echo "❌ Setup.py build failed. Trying alternative approaches..."
    
    # Method 2: Try manual compilation
    echo -e "\n🔧 Attempting manual compilation..."
    
    # Get pybind11 includes
    PYBIND11_INCLUDES=$(python3 -m pybind11 --includes)
    PYTHON_SUFFIX=$(python3-config --extension-suffix)
    
    echo "Using pybind11 includes: $PYBIND11_INCLUDES"
    echo "Python extension suffix: $PYTHON_SUFFIX"
    
    # Compile manually
    g++ -O2 -Wall -shared -std=c++20 -fPIC \
        $PYBIND11_INCLUDES \
        -I./include \
        -I/u/sw/toolchains/gcc-glibc/11.2.0/pkgs/eigen/3.3.9/include/eigen3 \
        -fopenmp \
        src/python_bindings.cpp \
        src/Interface.cpp \
        src/Parameters.cpp \
        src/StochasticUzawaSolver.cpp \
        src/Kernel.cpp \
        src/NystromScheme.cpp \
        src/LSMCR.cpp \
        src/OUSimulator.cpp \
        -o stochastic_uzawa$PYTHON_SUFFIX \
        -fopenmp
    
    if [ $? -eq 0 ]; then
        echo "✅ Manual compilation successful!"
        
        echo -e "\n🧪 Testing manually compiled module..."
        python3 -c "import stochastic_uzawa; print('✅ Manual compilation - Module import successful!')"
    else
        echo "❌ Manual compilation also failed."
        echo "💡 Suggestions:"
        echo "   - Check that all dependencies are installed"
        echo "   - Verify Eigen path: /u/sw/toolchains/gcc-glibc/11.2.0/pkgs/eigen/3.3.9/include/eigen3"
        echo "   - Try: pip3 install pybind11[global] numpy"
        echo "   - Check compiler version supports C++20"
        exit 1
    fi
fi

echo -e "\n🎉 Build process completed!"
echo "📖 Usage:"
echo "   python3 -c 'import stochastic_uzawa; interface = stochastic_uzawa.Interface()'"
echo "   python3 test_python_bindings.py  # Run full test suite"

echo -e "\n📋 Next steps:"
echo "   1. Run: python3 test_python_bindings.py"
echo "   2. Check the test output for functionality verification"
echo "   3. Use the module in your Python scripts!"
