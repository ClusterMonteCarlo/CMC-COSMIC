# CMC 

# This minimum version is mostly set in order to get a newer version
# of the FindMPI check.  Note that you can easily install a newer cmake version
# using conda or pip.
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

foreach(policy
    CMP0048
    CMP0074
    CMP0077
  )
  if(POLICY ${policy})
    cmake_policy(SET ${policy} NEW)
  endif()
endforeach()

project(CMC VERSION 1.0.0 LANGUAGES C Fortran)

# The version number.
set (CMC_VERSION_MAJOR 1)
set (CMC_VERSION_MINOR 0)
set (CMCPRETTYNAME "\"ClusterMonteCarlo\"")
 
# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/CMCConfig.h.in"
  "${PROJECT_BINARY_DIR}/CMCConfig.h"
  )

install (FILES "${PROJECT_BINARY_DIR}/CMCConfig.h"        
         DESTINATION include)

install(DIRECTORY include/bse_wrap DESTINATION include)
install(DIRECTORY include/cmc DESTINATION include)
install(DIRECTORY include/common DESTINATION include)
install(DIRECTORY include/fewbody-0.24 DESTINATION include)

# add the binary tree to the search path for include files
# so that we will find CMCConfig.h
include_directories("${PROJECT_BINARY_DIR}")

# Auxiliary files
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Use GNUInstallDirs to install libraries into correct locations on all
# platforms.
include(GNUInstallDirs)

# Build defaults
#include(BuildType)

# We are building libraries that will eventually be linked into shared
# modules.  All code should be built with PIC.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# External packages

find_package(MPI REQUIRED)

find_package(GSL REQUIRED)

find_package(CFITSIO REQUIRED)

find_package(ZLIB REQUIRED)

#find_package(HDF5 COMPONENTS C HL REQUIRED)

# enable USE_MPI
option(WITH_MPI "Option description" ON)
add_definitions(-DUSE_MPI)

# Internal products

add_subdirectory(src)
