cmake_minimum_required(VERSION 3.12)
project(gssolver LANGUAGES CXX)

# GCC: `--std=c++20`
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# generate `compile_commands.json`
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# GCC: `-g`
set(CMAKE_BUILD_TYPE Debug)

# GCC: `-I<PROJECT_ROOT>/include`
include_directories(${PROJECT_SOURCE_DIR}/include)

add_library(base STATIC ./src/base.cpp)
add_library(dataio STATIC ./src/dataio.cpp)

find_package(PkgConfig)
if(PKG_CONFIG_FOUND)
  pkg_check_modules(NETCDF_CXX netcdf-cxx4)
endif()

if(NETCDF_CXX_FOUND)
  target_include_directories(dataio PRIVATE ${NETCDF_CXX_INCLUDE_DIRS})
  target_link_libraries(dataio PRIVATE ${NETCDF_CXX_LIBRARIES})
  target_compile_definitions(dataio PUBLIC USE_NETCDF)
  message(STATUS "NetCDF supported")
else()
  message(WARNING "NetCDF not supported")
endif()

add_executable(solver ./src/solver.cpp)
target_link_libraries(solver PRIVATE base dataio)

add_custom_target(debug COMMAND gdb solver DEPENDS solver WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
