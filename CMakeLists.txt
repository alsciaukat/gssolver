cmake_minimum_required(VERSION 4.1)
project(gssolver LANGUAGES CXX)

include_directories(${PROJECT_SOURCE_DIR}/include)
add_library(base OBJECT ./src/base.cpp)
add_executable(solver ./src/solver.cpp)
target_link_libraries(solver PRIVATE base)

set(CMAKE_EXPORT_COMPILE_COMMANDS=ON)
set(CMAKE_BUILD_TYPE=Debug)

add_custom_target(debug COMMAND gdb solver DEPENDS solver WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
