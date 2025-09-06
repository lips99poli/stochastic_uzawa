cmake_minimum_required(VERSION 3.16)
project(StochasticUzawa CXX)

# Option to build Python bindings
option(BUILD_BINDINGS "Build Python bindings" OFF)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compiler flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Default to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Set binary directory to tests/cpp/build
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/tests/cpp/build)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# Find required packages
find_package(Eigen3 QUIET)
if(NOT Eigen3_FOUND)
    # Try to find Eigen3 in common locations
    find_path(EIGEN3_INCLUDE_DIR 
        NAMES Eigen/Dense
        PATHS 
            /u/sw/toolchains/gcc-glibc/11.2.0/pkgs/eigen/3.3.9/include/eigen3
            /usr/include/eigen3
            /usr/local/include/eigen3
            /opt/homebrew/include/eigen3
            /usr/include
            /usr/local/include
    )
    
    if(EIGEN3_INCLUDE_DIR)
        message(STATUS "Found Eigen3 at: ${EIGEN3_INCLUDE_DIR}")
        add_library(Eigen3::Eigen INTERFACE IMPORTED)
        target_include_directories(Eigen3::Eigen INTERFACE ${EIGEN3_INCLUDE_DIR})
    else()
        message(FATAL_ERROR "Eigen3 not found. Please install libeigen3-dev package.")
    endif()
endif()

find_package(OpenMP QUIET)  # Optional now
if(OpenMP_FOUND)
    message(STATUS "OpenMP found but NOT using it - Eigen optimization enabled")
else()
    message(STATUS "OpenMP not found - using Eigen optimization")
endif()

# Enable Eigen's automatic parallelization and optimizations
add_compile_definitions(EIGEN_USE_BLAS)
add_compile_definitions(EIGEN_USE_LAPACKE)
# Remove EIGEN_DONT_PARALLELIZE to enable Eigen threading
# Remove EIGEN_NO_AUTOMATIC_RESIZING to allow Eigen optimizations

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# All source files (library sources + test source)
file(GLOB_RECURSE ALL_SOURCES "${CMAKE_SOURCE_DIR}/src/*.cpp")

if(BUILD_BINDINGS)
    # Build Python bindings
    message(STATUS "Building with Python bindings")
    
    # Find Python and pybind11
    find_package(Python COMPONENTS Interpreter Development REQUIRED)
    
    # Try to find pybind11, fallback to Python module import
    find_package(pybind11 QUIET)
    if(NOT pybind11_FOUND)
        message(STATUS "pybind11 not found via CMake, trying Python import")
        execute_process(
            COMMAND ${Python_EXECUTABLE} -m pybind11 --cmakedir
            OUTPUT_VARIABLE pybind11_DIR
            OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE pybind11_RESULT
        )
        if(pybind11_RESULT EQUAL 0)
            message(STATUS "Found pybind11 via Python at: ${pybind11_DIR}")
            find_package(pybind11 REQUIRED PATHS ${pybind11_DIR})
        else()
            message(FATAL_ERROR "pybind11 not found. Install with: pip install pybind11")
        endif()
    endif()
    
    # Create Python module
    pybind11_add_module(stochastic_uzawa 
        ${ALL_SOURCES}
        ${CMAKE_SOURCE_DIR}/bindings/python_bindings.cpp
    )
    
    # Link libraries to the Python module
    target_link_libraries(stochastic_uzawa PRIVATE 
        Eigen3::Eigen
    )
    
    # Add BLAS/LAPACK if available for maximum performance
    find_package(BLAS QUIET)
    find_package(LAPACK QUIET)
    if(BLAS_FOUND AND LAPACK_FOUND)
        target_link_libraries(stochastic_uzawa PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
        message(STATUS "Linked with BLAS/LAPACK for optimal Eigen performance")
    endif()
    
    # Install the Python module
    install(TARGETS stochastic_uzawa DESTINATION .)
    
    # Set properties for the Python module
    target_compile_definitions(stochastic_uzawa PRIVATE VERSION_INFO=${PROJECT_VERSION})
    
else()
    # Build test executable
    message(STATUS "Building test executable only")
    list(APPEND ALL_SOURCES "${CMAKE_SOURCE_DIR}/tests/cpp/test_cpp.cpp")
    
    add_executable(test_cpp ${ALL_SOURCES})
    target_link_libraries(test_cpp Eigen3::Eigen)
    
    # Add BLAS/LAPACK if available for maximum performance
    find_package(BLAS QUIET)
    find_package(LAPACK QUIET)
    if(BLAS_FOUND AND LAPACK_FOUND)
        target_link_libraries(test_cpp ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
        message(STATUS "Linked with BLAS/LAPACK for optimal Eigen performance")
    endif()
    
    # Link filesystem library for std::filesystem (needed in C++17)
    target_link_libraries(test_cpp stdc++fs)
endif()
