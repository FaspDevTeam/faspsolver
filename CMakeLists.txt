# CMakeLists for FASP 
# 
# Author            Date           Action
# -----------------------------------------------------------------------------------
# Ludmil Zikatanov  Aug/08/2015    Create CMakeLists
# Chensong Zhang    Oct/20/2017    Update to new FASP code
# Ludmil Zikatanov  Aug/08/2019    Fix installation problem for cmake 3.14+
# Chensong Zhang    Mar/13/2021    Fix OpenMP problems
#
# Some sample usages:
#   cmake .                                 // build in default configuration
#   cmake -DCMAKE_BUILD_TYPE=Release .      // build in Release configuration
#   cmake -DCMAKE_BUILD_TYPE=Debug .        // build in Debug configuration
#   cmake -DCMAKE_VERBOSE_MAKEFILE=ON .     // build with verbose on
#   cmake -DUSE_UMFPACK=ON .                // build with UMFPACK package support

#####################################################################################
## General environment setting
#####################################################################################

# Minimum cmake version needed
cmake_minimum_required (VERSION 2.8.12)

# Helper modules
include(CheckFunctionExists)
include(CheckIncludeFile)

set(GDB         1 CACHE BOOL "debugging or not")
set(OPENMP      0 CACHE BOOL "Openmp use")
set(USE_MUMPS   0 CACHE BOOL "MUMPS use")
set(USE_UMFPACK 0 CACHE BOOL "UMFPACK use")
set(USE_SUPERLU 0 CACHE BOOL "SUPERLU use")
set(USE_PARDISO 0 CACHE BOOL "PARDISO use")
set(USE_DOXYGEN 0 CACHE BOOL "Doxygen use")

# Search for C compilers in the specified order
find_program(THE_C NAMES $ENV{CC} gcc icc clang)
set(CMAKE_C_COMPILER ${THE_C} CACHE INTERNAL "the C compiler" FORCE)

# Name the project at this moment and this will set the compiler id:
project(FASP C)

set(CMAKE_MODULE_PATH ${FASP_SOURCE_DIR}/modules)

# Default configuration types and flags
set(CMAKE_CONFIGURATION_TYPES "Release;Debug" CACHE STRING "Configs" FORCE)
if (MSVC)
    # Compiler flags for Windows MSVC, skip Fortran
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /O3")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /W4 /WX")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O3")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /W4 /WX")
else()
    # Compiler flags for Linux or Mac OS X
    include(get_cxx_fc_compilers) # find matching Fortran compiler
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -g")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g")
endif()

# Set default build type to Release
set(build_type "Release")
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "${build_type}" CACHE STRING "the build type" FORCE)
endif()
message(STATUS "Set build type to ${CMAKE_BUILD_TYPE} configuration")

########### check the prefix for installation #########################
include(install_prefix)
 
#######################################################################
# OpenMP : defined on command line in the top Makefile
# For cmake since 3.9, we can use find_package(OpenMP REQUIRED)
# We keep this for older cmake versions!!!
#######################################################################

if(USE_OPENMP)
    if(CMAKE_VERSION VERSION_GREATER 3)
        cmake_policy(SET "CMP0054" NEW)
    endif()

    find_package(OpenMP)
    if(OPENMP_FOUND)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")	 
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
    else(OPENMP_FOUND)
        message("   WARNING: OpenMP was requested but disabled!")
    endif(OPENMP_FOUND)
endif(USE_OPENMP)

#####################################################################################
## Project specific parameters
#####################################################################################

if(ADD_CFLAGS)
    set(CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE} 
        "${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}} ${ADD_CFLAGS}")	 
endif(ADD_CFLAGS)

if(ADD_CXXFLAGS)
    set(CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE} 
        "${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}} ${ADD_CXXFLAGS}")
endif(ADD_CXXFLAGS)

if(ADD_FFLAGS)
    set(CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE} 
        "${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}} ${ADD_FFLAGS}")	 
endif(ADD_FFLAGS)

# FOR MAC OS X to find shared libs in install location
set(CMAKE_INSTALL_NAME_DIR 
    ${CMAKE_INSTALL_PREFIX}/lib CACHE PATH "path for Mac OS X")

# FOR LINUX to find shared libs in install location
set(CMAKE_INSTALL_RPATH
    ${CMAKE_INSTALL_NAME_DIR} CACHE PATH "path for Linux")

# Set source files required for FASP
file(GLOB FASP_C_SOURCES RELATIVE ${FASP_SOURCE_DIR} 
	  		                      ${FASP_SOURCE_DIR}/base/src/*.c
	  		                      ${FASP_SOURCE_DIR}/base/src/*.inl )

# Auto generating headers if needed
add_custom_target(headers ${FASP_SOURCE_DIR}/util/fhead.sh "${FASP_SOURCE_DIR}/base"
         	      DEPENDS ${FASP_C_SOURCES}
                  WORKING_DIRECTORY "${FASP_SOURCE_DIR}/util"
                  COMMENT "Auto-generating header file for C functions..."
                  VERBATIM)

#####################################################################################
## External libraries (Optional)
#####################################################################################
include(external)

#####################################################################################
## Installtion targets for lib and exe files
#####################################################################################
set(CMAKE_INSTALL_PREFIX "${FASP_SOURCE_DIR}" )
set(FASPLIB_BASE_PATH "base" CACHE PATH "base path to FASPLIB")
set(FASP_INSTALL TRUE)

if(SHARED) 
    set(FASP_LIBRARY_TYPE SHARED)
else(SHARED)
    set(FASP_LIBRARY_TYPE STATIC)
endif(SHARED)   

set(FASP_LIBRARY_NAME fasp)

# Uninstall target
if (NOT TARGET uninstall)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/modules/cmake_uninstall.cmake.in"
                   "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
                   IMMEDIATE @ONLY)
    add_custom_target(uninstall
        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()

# Recursively look for CMakeLists.txt in the library subdir
add_subdirectory("base")
add_subdirectory("test")      # Requires Fortran compiler
add_subdirectory("tutorial")  # Requires Fortran compiler
add_subdirectory("benchmark") # Only needed for benchmark tests

include(write_config)
