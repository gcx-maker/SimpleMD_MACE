cmake_minimum_required(VERSION 3.18)

project(simplemd_libtorch_mpi)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


set(Torch_DIR
"D:/Software/libtorch/libtorch210/libtorch/share/cmake/Torch"
)

find_package(Torch REQUIRED)


set(MPI_HOME "D:/Software/MPI")


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


target_include_directories(main PRIVATE

    ${MPI_HOME}/Include

)


target_link_directories(main PRIVATE

    ${MPI_HOME}/Lib/x64

)


target_link_libraries(main PRIVATE

    "${TORCH_LIBRARIES}"

    msmpi.lib

)


set_target_properties(main PROPERTIES

    RUNTIME_OUTPUT_DIRECTORY
    ${CMAKE_SOURCE_DIR}

)