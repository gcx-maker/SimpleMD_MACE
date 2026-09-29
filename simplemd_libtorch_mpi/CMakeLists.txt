cmake_minimum_required(VERSION 3.18)

project(simplemd_libtorch_mpi)


# ==========================
# C++ Standard
# ==========================

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# ==========================
# LibTorch
# ==========================

# 用户需要提供:
# -DTorch_DIR=/path/to/libtorch/share/cmake/Torch

find_package(Torch REQUIRED)



# ==========================
# MPI
# ==========================

find_package(MPI REQUIRED)



# ==========================
# Source files
# ==========================

add_executable(main

    main.cpp

    system.cpp
    Vector.cpp
    Random.cpp
    input.cpp
    xyzreader.cpp

    neighbor.cpp

    integrator.cpp
    thermostat.cpp

    output.cpp
    setup.cpp


    ljforcefield.cpp

    maceforcefield.cpp
    macedatabuilder.cpp

)



# ==========================
# Include directories
# ==========================

target_include_directories(main PRIVATE

    ${MPI_CXX_INCLUDE_DIRS}

)



# ==========================
# Link libraries
# ==========================

target_link_libraries(main PRIVATE

    ${TORCH_LIBRARIES}

    MPI::MPI_CXX

)



# ==========================
# Torch settings
# ==========================

set_property(TARGET main PROPERTY CXX_STANDARD 17)



# ==========================
# Output directory
# ==========================

set_target_properties(main PROPERTIES

    RUNTIME_OUTPUT_DIRECTORY
    ${CMAKE_BINARY_DIR}/bin

)