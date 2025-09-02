cmake_minimum_required(VERSION 3.16)  # 3.16+ is recommended for C++20 features

project(StochasticUzawa VERSION 0.1 LANGUAGES CXX)

# Enable C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)  # portable, avoids compiler-specific extensions

# Include header files
include_directories(${CMAKE_SOURCE_DIR}/include)

# Gather source files
file(GLOB SOURCES "src/*.cpp")
set(MAIN_SRC "${CMAKE_SOURCE_DIR}/src/main.cpp")

# Define the executable
add_executable(uzawa_exec ${SOURCES} ${MAIN_SRC})

# Include header files
target_include_directories(uzawa_exec PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    /u/sw/toolchains/gcc-glibc/11.2.0/pkgs/eigen/3.3.9/include/eigen3
)

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
