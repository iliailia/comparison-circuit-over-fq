cmake_minimum_required(VERSION 3.11 FATAL_ERROR)

## Use -std=c++17 as default.
set(CMAKE_CXX_STANDARD 17)
## Disable C++ extensions
set(CMAKE_CXX_EXTENSIONS OFF)
## Require full C++ standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(comparison_circuit
       	LANGUAGES CXX)

find_package(helib REQUIRED)

include_directories(${PROJECT_SOURCE_DIR})

add_executable(comparison_circuit comparison_circuit.cpp comparator.cpp tools.cpp)
add_executable(sorting_circuit sorting_circuit.cpp comparator.cpp tools.cpp)
add_executable(min_max_circuit min_max_circuit.cpp comparator.cpp tools.cpp)

target_link_libraries(comparison_circuit helib)
target_link_libraries(sorting_circuit helib)
target_link_libraries(min_max_circuit helib)

