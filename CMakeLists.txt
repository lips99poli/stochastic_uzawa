cmake_minimum_required(VERSION 3.10)

project(StochasticUzawa VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(include)

file(GLOB SOURCES "src/*.cpp")

add_executable(uzawa_exec ${SOURCES})