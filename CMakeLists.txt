############################################################
# CMakeLists for FASP 
# Modified   2015-08-08   --Ludmil Zikatanov
# Modified   2017-10-20   --Chensong Zhang
# Modified   2019-08-08   --Ludmil Zikatanov & Chensong Zhang
############################################################

cmake_minimum_required (VERSION 2.8)

# Helper modules
include(CheckFunctionExists)
include(CheckIncludeFile)

set(CMAKE_VERBOSE_MAKEFILE 1) 
set(GDB 1 CACHE BOOL "debugging or not")
set(OPENMP 0 CACHE BOOL "Openmp use")
set(USE_MUMPS 0 CACHE BOOL "MUMPS use")

#######################################################################
# Search for C compilers in the specified order.
#######################################################################
if(DEFINED ENV{CC}) 
	find_program(THE_C NAMES $ENV{CC} gcc-mp-5 gcc-mp-4.9 gcc-mp-4.8
	             gcc-mp-4.6 gcc46 gcc45 gcc44 gcc icc clang)
else(DEFINED ENV{CC}) 
	find_program(THE_C NAMES gcc-mp-5 gcc-mp-4.9 gcc-mp-4.8 gcc-mp-4.6
	             gcc46 gcc45 gcc44 gcc icc clang)
endif(DEFINED ENV{CC}) 

if( THE_C )
    set(CMAKE_C_COMPILER ${THE_C} CACHE INTERNAL "the C compiler" FORCE)
endif( THE_C )

#######################################################################
# name the project at this moment and this will set the compiler id:
project(FASP C)
#######################################################################
set(CMAKE_MODULE_PATH ${FASP_SOURCE_DIR}/modules)

########### try to find matching c++ and fortran compilers ############
# Fortran compiler is needed for Fortran examples and demos
#######################################################################
include(get_cxx_fc_compilers)

########### check the prefix for installation #########################
include(install_prefix)
 
#######################################################################
# OpenMP : defined on command line in the top Makefile
#######################################################################
if(USE_OPENMP)
    if(CMAKE_VERSION VERSION_GREATER 3)
        cmake_policy(SET "CMP0054" NEW)
    endif()
    find_package(OpenMP)

    if(OPENMP_FOUND)
        set (CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE} 
             "${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}} ${OpenMP_C_FLAGS}")	 
        set (CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE} 
	     "${CMAKE_CXX_FLAGS_${CMAKE_CBUILD_TYPE}} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE} 
	     "${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}} ${OpenMP_C_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    else(OPENMP_FOUND)
        message(" WARNING: OpenMP was requested but not supported!")
    endif(OPENMP_FOUND)
endif(USE_OPENMP)

########## Additional compiler flags (not defined by the build
#
if(ADD_CFLAGS)
    set (CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE} 
	 "${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}} ${ADD_CFLAGS}")	 
endif(ADD_CFLAGS)
#
if(ADD_CXXFLAGS)
    set (CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE} 
	 "${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}} ${ADD_CXXFLAGS}")
endif(ADD_CXXFLAGS)
#
if(ADD_FFLAGS)
    set (CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE} 
	 "${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}} ${ADD_FFLAGS}")	 
endif(ADD_FFLAGS)
#
    set(CMAKE_INSTALL_PREFIX "${FASP_SOURCE_DIR}" )
    set(FASPLIB_BASE_PATH "base" CACHE PATH "base path to FASPLIB")
    set(FASP_INSTALL TRUE)
#
    if(SHARED) 
	set(FASP_LIBRARY_TYPE SHARED)
    else(SHARED)
	set(FASP_LIBRARY_TYPE STATIC)
    endif(SHARED)   
#

# FOR MAC OS X to find shared libs in install location
    set(CMAKE_INSTALL_NAME_DIR 
	${CMAKE_INSTALL_PREFIX}/lib CACHE PATH "path for Mac OS X")
# FOR LINUX to find shared libs in install location
    set(CMAKE_INSTALL_RPATH
        ${CMAKE_INSTALL_NAME_DIR} CACHE PATH "path for Linux")

file(GLOB FASP_C_SOURCES RELATIVE ${FASP_SOURCE_DIR} 
	  		          ${FASP_SOURCE_DIR}/base/src/*.c 
	  		          ${FASP_SOURCE_DIR}/base/src/*.inl )

# keep here for auto generating headers if needed
add_custom_target(headers ${FASP_SOURCE_DIR}/util/fhead.sh "${FASP_SOURCE_DIR}/base"
         	  DEPENDS ${FASP_C_SOURCES}
                  WORKING_DIRECTORY "${FASP_SOURCE_DIR}/util"
                  COMMENT 
		  "FASP: Autogenerating header file with C functions..."
                  VERBATIM)

##################################################################
# For UMFPACK
##################################################################
if(USE_UMFPACK)

   # set some path to the UMFPACK pacakge
   # metis is not part of suitesparse, so theremay be also some other metis dir. 
   set(METIS_DIR "${SUITESPARSE_DIR}")

   find_package(UMFPACK)     
   if (UMFPACK_FOUND)
      add_definitions("-DWITH_UMFPACK=1")
      include_directories(${UMFPACK_INCLUDE_DIRS})
   else(UMFPACK_FOUND)
      message("   WARNING: UMFPACK was requested but not supported!")
   endif(UMFPACK_FOUND)

endif(USE_UMFPACK)

##################################################################
# For SuperLU
##################################################################
if(USE_SUPERLU)

   # set the path to find specific modules
   set(CMAKE_MODULE_PATH "${FASP_SOURCE_DIR}/modules")

   set(SUPERLU_DIR "${SUPERLU_DIR}")

   # try to find SuperLU
   find_package(SUPERLU)

   if (SUPERLU_FOUND)
      add_definitions("-DWITH_SuperLU=1")
      include_directories(${SUPERLU_INCLUDE_DIRS})
   else(SUPERLU_FOUND)
      message("   WARNING: SuperLU was requested but not supported!")
   endif(SUPERLU_FOUND)

endif(USE_SUPERLU)

##################################################################
# For MUMPS 
##################################################################
if(USE_MUMPS)

   # set the path to find specific modules
   set(CMAKE_MODULE_PATH "${FASP_SOURCE_DIR}/modules")

   set(MUMPS_DIR "${MUMPS_DIR}")

   # try to find MUMPS
   find_package(MUMPS)

   if (MUMPS_FOUND)
      add_definitions("-DWITH_MUMPS=1")
      include_directories(${MUMPS_INCLUDE_DIRS})
   else(MUMPS_FOUND)
      message("   WARNING: MUMPS was requested but not supported!")
   endif(MUMPS_FOUND)

endif(USE_MUMPS)

##################################################################
# For Intel MKL PARDISO
##################################################################
if(USE_PARDISO)

   # set the path to find specific modules
   set(CMAKE_MODULE_PATH "${FASP_SOURCE_DIR}/modules")

   set(MKL_DIR "${MKL_DIR}")

   # try to find MKL
   find_package(MKL)

   if (MKL_FOUND)
       add_definitions("-DWITH_PARDISO=1")
       include_directories(${MKL_INCLUDE_DIRS})
   else(MKL_FOUND)
       message("   WARNING: Intel MKL was requested but not supported!")
   endif(MKL_FOUND)

endif(USE_PARDISO)

##################################################################
# For Doxygen
##################################################################
# target generating the FASP documentation with doxygen
find_package(Doxygen)

if(DOXYGEN_FOUND)
    if(EXISTS ${FASP_SOURCE_DIR}/doc/fasp.Doxygen.cnf.in)
    configure_file(
	${FASP_SOURCE_DIR}/doc/fasp.Doxygen.cnf.in
        ${CMAKE_CURRENT_BINARY_DIR}/fasp.Doxygen.cnf @ONLY)
    set(DOXY_EXEC "${DOXYGEN_EXECUTABLE}")
    if(DOXYWIZARD)
        find_program(WIZARD doxywizard)
        if(APPLE AND (NOT WIZARD) )
           find_program(WIZARD 
	       /Applications/Doxygen.app/Contents/MacOS/Doxywizard)
        endif()
        if(WIZARD)
	       set(DOXY_EXEC "${WIZARD}")
	    endif()
    endif(DOXYWIZARD)
    add_custom_target(docs ${DOXY_EXEC}
 	       ${CMAKE_CURRENT_BINARY_DIR}/fasp.Doxygen.cnf
	       WORKING_DIRECTORY 
	       "${CMAKE_CURRENT_BINARY_DIR}"
	       COMMENT 
	       "Generating FASP documentation (Doxygen)" 
	       VERBATIM)
    else(EXISTS ${FASP_SOURCE_DIR}/doc/fasp.Doxygen.cnf.in)
        message("   WARNING: Doxygen configuration file cannot be found!")
    endif(EXISTS ${FASP_SOURCE_DIR}/doc/fasp.Doxygen.cnf.in)
endif(DOXYGEN_FOUND)
###############end documentation with doxygen

#MISMATCH CHECK:
if(NOT (CMAKE_Fortran_COMPILER_ID STREQUAL CMAKE_C_COMPILER_ID))
  message("
** WARNING: Mismatching CC and FC mismatch. CC=\"${CMAKE_C_COMPILER_ID}${CMAKE_C_COMPILER_VERSION}\"; FC=\"${CMAKE_Fortran_COMPILER_ID}${CMAKE_Fortran_COMPILER_VERSION}\"!
"
) 
endif(NOT (CMAKE_Fortran_COMPILER_ID STREQUAL CMAKE_C_COMPILER_ID))

##################################################################
# set FASP library name for the library
    set(FASP_LIBRARY_NAME fasp)
##################################################################
# Recursively look for CMakeLists.txt in the library subdir.
    add_subdirectory("base")
    add_subdirectory("test")
    add_subdirectory("tutorial")
# write the configuration
    include(write_config)
    
##################################################################
# Only the library will be installed with "make". The rest of the
# targets must be installed separately, so we do EXCLUDE_FROM_ALL
##################################################################
# We no longer treat "test" and "tutorial" separately due to
# cmake 3.14+ changes its behavoir on EXCLUDE_FROM_ALL!!! It
# makes the compilation results unpredictable.
##################################################################
#    add_subdirectory("test" EXCLUDE_FROM_ALL)
#    add_subdirectory("tutorial" EXCLUDE_FROM_ALL)
