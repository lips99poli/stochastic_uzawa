cmake_minimum_required(VERSION 3.16)
project(StochasticUzawa CXX)

# Option to build Python bindings
option(BUILD_BINDINGS "Build Python bindings" OFF)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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

find_package(OpenMP REQUIRED)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# All source files (library sources + test source)
file(GLOB_RECURSE ALL_SOURCES "${CMAKE_SOURCE_DIR}/src/*.cpp")

if(BUILD_BINDINGS)
    # Build Python bindings
    message(STATUS "Building with Python bindings")
    # TODO: Add pybind11 and bindings compilation here
else()
    # Build test executable
    message(STATUS "Building test executable only")
    list(APPEND ALL_SOURCES "${CMAKE_SOURCE_DIR}/tests/cpp/test_cpp.cpp")
    
    add_executable(test_cpp ${ALL_SOURCES})
    target_link_libraries(test_cpp Eigen3::Eigen OpenMP::OpenMP_CXX)
    
    # Link filesystem library for std::filesystem (needed in C++17)
    target_link_libraries(test_cpp stdc++fs)
endif()
