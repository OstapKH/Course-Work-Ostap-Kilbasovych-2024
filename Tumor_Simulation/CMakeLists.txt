# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

cmake_minimum_required(VERSION 3.1...3.23)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpic++")
#set(CMAKE_C_COMPILER "gcc")
#set(CMAKE_CXX_COMPILER "g++")

project(Tumor_Simulation)

# MPI and its components
find_package(MPI REQUIRED)
if (NOT MPI_FOUND)
    message(STATUS "This project requires the OpenMPI or MPICH libraries, and will not be compiled.")
    return()
endif ()

# CGAL and its components
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)
if (NOT CGAL_FOUND)
    message(STATUS "This project requires the CGAL library, and will not be compiled.")
    return()
endif ()


# Boost and its components
find_package(Boost REQUIRED)
if (NOT Boost_FOUND)
    message(STATUS "This project requires the Boost library, and will not be compiled.")
    return()
endif ()

# include for local directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# include for local package

# Creating entries for all C++ files with "main" routine
# ##########################################################
#create_single_source_cgal_program( "automata_2d.cpp" )
#create_single_source_cgal_program( "main.cpp" )

set(SOURCES
        #  main.cpp
        #  automata_2d.cpp
        Automata.hpp
        Automata.cpp
        Cell.hpp
        Cell.cpp
        DiffusionReaction2D.hpp
        DiffusionReaction2D.cpp
        DiffusionReaction3D.hpp
        DiffusionReaction3D.cpp
        AutomataMPI.cpp
        AutomataMPI.hpp
        AutomataMPI3D.cpp
        AutomataMPI3D.hpp
        )
add_executable(main main.cpp ${SOURCES})
target_link_libraries(main PUBLIC CGAL::CGAL Boost::boost ${MPI_CXX_LIBRARIES})

add_executable(automata_2d automata_2d.cpp ${SOURCES})
target_link_libraries(automata_2d PUBLIC CGAL::CGAL Boost::boost ${MPI_CXX_LIBRARIES})

add_executable(automataMPI main_automataMPI.cpp ${SOURCES})
target_link_libraries(automataMPI PUBLIC CGAL::CGAL Boost::boost ${MPI_CXX_LIBRARIES})

add_executable(automataMPI3D main_automataMPI3D.cpp ${SOURCES})
target_link_libraries(automataMPI3D PUBLIC CGAL::CGAL Boost::boost ${MPI_CXX_LIBRARIES})

add_executable(mainDiffusionReaction mainDiffusionReaction.cpp ${SOURCES})
target_link_libraries(mainDiffusionReaction PUBLIC CGAL::CGAL Boost::boost ${MPI_CXX_LIBRARIES})
