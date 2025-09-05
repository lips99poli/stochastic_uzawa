cmake_minimum_required(VERSION 3.16)  # 3.16+ is recommended for C++20 features

project(StochasticUzawa VERSION 0.1 LANGUAGES CXX)

# Enable C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)  # portable, avoids compiler-specific extensions

# Find required packages
find_package(OpenMP REQUIRED)

# Include header files
include_directories(${CMAKE_SOURCE_DIR}/include)

# Define source files (excluding main.cpp for library)
set(LIB_SOURCES
    src/Interface.cpp
    src/Parameters.cpp
    src/StochasticUzawaSolver.cpp
    src/Kernel.cpp
    src/NystromScheme.cpp
    src/LSMCR.cpp
    src/OUSimulator.cpp
)

# Define the main executable
add_executable(uzawa_exec src/main.cpp ${LIB_SOURCES})

# Include header files for executable
target_include_directories(uzawa_exec PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    /u/sw/toolchains/gcc-glibc/11.2.0/pkgs/eigen/3.3.9/include/eigen3
)

# Link OpenMP for executable
target_link_libraries(uzawa_exec PRIVATE OpenMP::OpenMP_CXX)

# Optional: Add pybind11 target if pybind11 is available
find_package(pybind11 QUIET)
if(pybind11_FOUND)
    message(STATUS "pybind11 found - building Python bindings")
    
    # Create the Python module
    pybind11_add_module(stochastic_uzawa 
        src/python_bindings.cpp
        ${LIB_SOURCES}
    )
    
    # Include directories for Python module
    target_include_directories(stochastic_uzawa PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        /u/sw/toolchains/gcc-glibc/11.2.0/pkgs/eigen/3.3.9/include/eigen3
    )
    
    # Compiler options for Python module
    target_compile_definitions(stochastic_uzawa PRIVATE VERSION_INFO="1.0.0")
    target_compile_options(stochastic_uzawa PRIVATE -O2 -w)
    
    # Link OpenMP for Python module
    target_link_libraries(stochastic_uzawa PRIVATE OpenMP::OpenMP_CXX)
    
    # Set properties
    set_target_properties(stochastic_uzawa PROPERTIES
        CXX_VISIBILITY_PRESET "hidden"
        INTERPROCEDURAL_OPTIMIZATION TRUE
    )
else()
    message(STATUS "pybind11 not found - skipping Python bindings")
endif()

# Compiler warning suppression
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")  # Suppress all warnings
# Alternative: suppress specific warnings instead of all warnings
# add_compile_options(
#     -Wno-deprecated-enum-enum-conversion
#     -Wno-unused-variable
#     -Wno-unused-parameter
#     -Wno-sign-compare
#     -Wno-deprecated-declarations
#     -Wno-maybe-uninitialized
#     -Wno-unused-but-set-variable
# )
