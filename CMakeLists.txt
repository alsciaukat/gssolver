cmake_minimum_required(VERSION 4.0)
project(gssolver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS=ON)
set(CMAKE_BUILD_TYPE=Debug)

include_directories(${PROJECT_SOURCE_DIR}/include)
add_library(base OBJECT ./src/base.cpp)
add_executable(solver ./src/solver.cpp)
target_link_libraries(solver PRIVATE base)

add_custom_target(debug COMMAND gdb solver DEPENDS solver WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
