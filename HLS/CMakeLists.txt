cmake_minimum_required(VERSION 3.0)
project(MyProject)

# Set C++11 standard
set(CMAKE_CXX_STANDARD 11)

# Add source files
set(SOURCES
    source/maxpool.cpp
    source/reshape.cpp
    source/convolution.cpp
)

# Add test files
set(TESTS
    tests/main.cpp
)

# Create executables
add_executable(MyProject ${SOURCES} ${TESTS})

# Set output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build)

