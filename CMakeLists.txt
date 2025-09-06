cmake_minimum_required(VERSION 3.16)
project(StochasticUzawa LANGUAGES CXX)

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

# --- Eigen ---
find_package(Eigen3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
  find_path(EIGEN3_INCLUDE_DIR
    NAMES Eigen/Dense
    PATHS
      /u/sw/toolchains/gcc-glibc/11.2.0/pkgs/eigen/3.3.9/include/eigen3
      /usr/include/eigen3 /usr/local/include/eigen3
      /opt/homebrew/include/eigen3
      /usr/include /usr/local/include
  )
  if(NOT EIGEN3_INCLUDE_DIR)
    message(FATAL_ERROR "Eigen3 not found")
  endif()
  message(STATUS "Found Eigen3 at: ${EIGEN3_INCLUDE_DIR}")
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  target_include_directories(Eigen3::Eigen INTERFACE ${EIGEN3_INCLUDE_DIR})
endif()

# --- Your custom BLAS/LAPACK (toolchain-built) ---
set(OPENBLAS_LIB "$ENV{HOME}/local/openblas-tc/lib/libopenblas.so")
if(NOT EXISTS "${OPENBLAS_LIB}")
  message(WARNING "OpenBLAS not found at ${OPENBLAS_LIB}, trying system BLAS")
  # Fallback to system BLAS if custom build not available
  find_package(BLAS)
  if(NOT BLAS_FOUND)
    message(FATAL_ERROR "No BLAS library found")
  endif()
endif()

set(LAPACKE_INC "$ENV{HOME}/local/lapack-tc/include")
set(LAPACKE_LIB "$ENV{HOME}/local/lapack-tc/lib/liblapacke.a")
if(NOT EXISTS "${LAPACKE_LIB}")
  message(WARNING "Custom LAPACKE not found at ${LAPACKE_LIB}, trying system LAPACK")
  find_package(LAPACK)
  if(NOT LAPACK_FOUND)
    message(FATAL_ERROR "No LAPACK library found")
  endif()
endif()

# --- Project includes & sources ---
include_directories(${CMAKE_SOURCE_DIR}/include)
file(GLOB_RECURSE LIB_SOURCES "${CMAKE_SOURCE_DIR}/src/*.cpp")

# --- Threading ---
find_package(Threads REQUIRED)

if(BUILD_BINDINGS)
    # Build Python bindings
    message(STATUS "Building with Python bindings and Eigen+BLAS/LAPACK optimization")
    
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
        ${LIB_SOURCES}
        ${CMAKE_SOURCE_DIR}/bindings/python_bindings.cpp
    )
    
    # Enable Eigen BLAS/LAPACK hooks for Python module
    target_compile_definitions(stochastic_uzawa PRIVATE
        EIGEN_USE_BLAS
        EIGEN_USE_LAPACKE
        EIGEN_NO_DEBUG
        NDEBUG
        VERSION_INFO=${PROJECT_VERSION}
    )
    
    # Optimization flags for Python module
    target_compile_options(stochastic_uzawa PRIVATE -O3 -march=native -mtune=native)
    
    # Include custom LAPACKE headers if available
    if(EXISTS "${LAPACKE_INC}")
        target_include_directories(stochastic_uzawa PRIVATE ${LAPACKE_INC})
    endif()
    
    # Link libraries to the Python module
    target_link_libraries(stochastic_uzawa PRIVATE Eigen3::Eigen)
    
    # Link custom BLAS/LAPACK if available, otherwise system libraries
    if(EXISTS "${LAPACKE_LIB}" AND EXISTS "${OPENBLAS_LIB}")
        message(STATUS "Using custom BLAS/LAPACK libraries for Python module")
        target_link_libraries(stochastic_uzawa PRIVATE
            ${LAPACKE_LIB}
            ${OPENBLAS_LIB}
            Threads::Threads
            m dl
        )
        # Set RPATH for custom OpenBLAS
        set_target_properties(stochastic_uzawa PROPERTIES
            BUILD_RPATH "$ENV{HOME}/local/openblas-tc/lib"
            INSTALL_RPATH "$ENV{HOME}/local/openblas-tc/lib"
        )
    else()
        message(STATUS "Using system BLAS/LAPACK libraries for Python module")
        if(BLAS_FOUND AND LAPACK_FOUND)
            target_link_libraries(stochastic_uzawa PRIVATE 
                ${LAPACK_LIBRARIES} 
                ${BLAS_LIBRARIES}
                Threads::Threads
            )
        endif()
    endif()
    
    # Install the Python module
    install(TARGETS stochastic_uzawa DESTINATION .)
    
else()
    # Build test executable
    message(STATUS "Building C++ test executable with Eigen+BLAS/LAPACK optimization")
    
    set(TEST_MAIN "${CMAKE_SOURCE_DIR}/tests/cpp/test_cpp.cpp")
    add_executable(test_cpp ${LIB_SOURCES} ${TEST_MAIN})
    
    target_link_libraries(test_cpp PRIVATE Eigen3::Eigen)
    
    # Enable Eigen BLAS/LAPACK hooks for test executable
    target_compile_definitions(test_cpp PRIVATE
        EIGEN_USE_BLAS
        EIGEN_USE_LAPACKE
        EIGEN_NO_DEBUG
        NDEBUG
    )
    
    # Optimization flags for test executable
    target_compile_options(test_cpp PRIVATE -O3 -march=native -mtune=native)
    
    # Include custom LAPACKE headers if available
    if(EXISTS "${LAPACKE_INC}")
        target_include_directories(test_cpp PRIVATE ${LAPACKE_INC})
    endif()
    
    # Link custom BLAS/LAPACK if available, otherwise system libraries
    if(EXISTS "${LAPACKE_LIB}" AND EXISTS "${OPENBLAS_LIB}")
        message(STATUS "Using custom BLAS/LAPACK libraries for test executable")
        target_link_libraries(test_cpp PRIVATE
            ${LAPACKE_LIB}
            ${OPENBLAS_LIB}
            Threads::Threads
            m dl
        )
        # Set RPATH for custom OpenBLAS
        set_target_properties(test_cpp PROPERTIES
            BUILD_RPATH "$ENV{HOME}/local/openblas-tc/lib"
            INSTALL_RPATH "$ENV{HOME}/local/openblas-tc/lib"
        )
    else()
        message(STATUS "Using system BLAS/LAPACK libraries for test executable")
        if(BLAS_FOUND AND LAPACK_FOUND)
            target_link_libraries(test_cpp PRIVATE 
                ${LAPACK_LIBRARIES} 
                ${BLAS_LIBRARIES}
                Threads::Threads
            )
        endif()
        # Link filesystem library for std::filesystem (needed in C++17)
        target_link_libraries(test_cpp stdc++fs)
    endif()
    
endif()