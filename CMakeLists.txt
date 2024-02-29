cmake_minimum_required(VERSION 3.10)
project(monolayer_film_tools)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# Include GROMACS libraries
message(STATUS "Locating GROMACS package...")
if (WIN32)
else ()
    find_package(GROMACS 2021 REQUIRED)
    gromacs_check_compiler(CXX)
    message(STATUS "GROMACS package at:${GROMACS_INCLUDE_DIRS}")
endif ()

# To build, we need gromacs. For some trajectory analysis functions, we require headers only available in the submodule
include_directories(${GROMACS_INCLUDE_DIRS} submodules/gromacs/src/)

#Include OpenMP capabilities
find_package(OpenMP)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../planar_bin)

set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/../planar_bin")

add_subdirectory(src)
