cmake_minimum_required(VERSION 3.0)
project(MyProject)

# Set C++11 standard
set(CMAKE_CXX_STANDARD 11)

include_directories(include)
include_directories(include/ac_fixed)
# Add source files
set(SOURCES
    src/maxPool.cpp
    src/reshape.cpp
    src/percetron.cpp
)

# Add test files
set(TESTS
    tests/test_maxPool.cpp
    tests/test_reshape.cpp
    tests/test_perceptron.cpp
    )

# Create executables
add_executable(test_maxPool src/maxPool.cpp tests/test_maxPool.cpp)
# add_executable(test_perceptron src/perceptron.cpp tests/test_perceptron.cpp)
# Set output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build)

