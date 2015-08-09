############################################################
# CMakeLists for FASP 
# Modified   2015-08-08   --ltz
############################################################
cmake_minimum_required (VERSION 2.8)
# Helper modules.
include(CheckFunctionExists)
include(CheckIncludeFile)

set(CMAKE_VERBOSE_MAKEFILE 1) 
set(GDB 1 CACHE BOOL "debugging or not")
set(OPENMP 0 CACHE BOOL "Openmp use")
set(USE_MUMPS 0 CACHE BOOL "MUMPS use")

# Search for C compilers in the specified order. That will determine the rest.
if(DEFINED ENV{CC}) 
	find_program(THE_C NAMES $ENV{CC} gcc-mp-4.9 gcc-mp-4.8 gcc-mp-4.6 gcc46 gcc45 gcc44 gcc icc clang)
else(DEFINED ENV{CC}) 
	find_program(THE_C NAMES gcc-mp-4.9 gcc-mp-4.8 gcc-mp-4.6 gcc46 gcc45 gcc44 gcc icc clang)
endif(DEFINED ENV{CC}) 

if( THE_C )
    set(CMAKE_C_COMPILER ${THE_C} CACHE INTERNAL   "the C compiler" FORCE)
endif( THE_C )

# name the project at this moment and this will set the compiler id:
project(FASP C )
##message("C compiler ID and Version: ${CMAKE_C_COMPILER_ID} and ${CMAKE_C_COMPILER_VERSION}" )
#message("The system is ${CMAKE_HOST_SYSTEM}; apple=${APPLE}; unix=${UNIX}")
set(REAL_C "${CMAKE_C_COMPILER_ID}${CMAKE_C_COMPILER_VERSION}")
############################
if(${REAL_C} MATCHES "GNU.*" AND  ${THE_C} MATCHES "gcc.*") 
     		string(REPLACE "gcc" "g++" C_XX ${THE_C} )
		string(REPLACE "gcc" "gfortran" F_C ${THE_C} )
		find_program(THE_CXX NAMES ${C_XX})
		find_program(THE_FC NAMES ${F_C})
#		message("GNU: ${THE_C} ${THE_CXX} ${THE_FC}" )
       elseif( ${REAL_C} MATCHES "Intel.*" AND ${THE_C} MATCHES "icc.*" )
		find_program(THE_CXX NAMES icpc)
		find_program(THE_FC NAMES ifort)
#		message("Intel: ${THE_CXX} ${THE_FC}" )
       elseif( ${REAL_C} MATCHES "Clang.*")
	 find_program(THE_C NAMES clang)
	find_program(THE_CXX NAMES clang++ clang)
		# clang is seems to work with gnu compilers if
		# dragonegg plugin is installed; it is not clear how
		# to include this dragonegg plugin on Linux in the
		# CMAKEs, so we will issue a warning at the end if the 
		#C compiler and fortran compilers are different. 
		#As compiler gfortran works with clang in many
		# cases, so we assume this standard behavior
		# here as well
		if(APPLE)
#		    message(SATATUS "MAC OS X and Clang. Looking for dragonegg install")
#		  find_program(THE_FC NAMES dragonegg-3.4-gfortran-mp-4.6  gfortran gfortran-mp-4.9 gfortran-mp-4.8  gfortran-mp-4.6 g95 g77)
		  find_program(THE_FC NAMES gfortran gfortran-mp-4.9 gfortran-mp-4.8  gfortran-mp-4.6 g95 g77)
		else(APPLE)
		  find_program(THE_FC NAMES gfortran gfortran-mp-4.9 gfortran-mp-4.8  gfortran-mp-4.6 g95 g77)
		endif(APPLE)
#		message("Clang: ${THE_CXX} ${THE_FC}" )
	else()       
	   message("WARNING: ${THE_C} did not match any of the preset C compilers" )
	   message("Continuing with the default compiler: ${CMAKE_C_COMPILER}" )
   	   set(THE_C "0")
   	   set(THE_CXX "0")	      
   	   set(THE_FC "0")
    endif(${REAL_C} MATCHES "GNU.*" AND  ${THE_C} MATCHES "gcc.*") 
#
	if( THE_C AND THE_CXX AND THE_FC )
	     	 set(CMAKE_C_COMPILER ${THE_C} CACHE INTERNAL   "the C   compiler" FORCE) 
	    	 set(CMAKE_CXX_COMPILER ${THE_CXX} CACHE INTERNAL   "the C++ compiler" FORCE)
                 set(CMAKE_Fortran_COMPILER ${THE_FC} CACHE INTERNAL    "the F compiler" FORCE)
	endif( THE_C AND THE_CXX AND THE_FC )
# END COMPILERS SET UP................ 
##################################
	enable_language(CXX) 
	enable_language(Fortran) 

#
# OpenMP : defined on command line in the top Makefile
#
    if(USE_OPENMP)
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

	if(ADD_CFLAGS)
          set (CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE} 
	  	  "${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}} ${ADD_CFLAGS}")	 
	endif(ADD_CFLAGS)
	if(ADD_CXXFLAGS)
          set (CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE} 
	  	  "${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}} ${ADD_CXXFLAGS}")	 
	endif(ADD_CXXFLAGS)
	if(ADD_FFLAGS)
          set (CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE} 
	  	  "${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}} ${ADD_FFLAGS}")	 
	endif(ADD_FFLAGS)

#       message(" CHECKC:  ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}") 
#       message(" CHECKCXX:  ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}") 
#       message(" CHECKF:  ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}") 
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
# FOR MAC OS X to find shared libs in install location
        set(CMAKE_INSTALL_NAME_DIR 
		${CMAKE_INSTALL_PREFIX}/lib CACHE PATH "path for apple")
# FOR LINUX to find shared libs in install location
        set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_NAME_DIR} CACHE PATH "path for Linux")

file(GLOB FASP_C_SOURCES RELATIVE ${FASP_SOURCE_DIR} 
	  		    ${FASP_SOURCE_DIR}/base/src/*.c 
	  		    ${FASP_SOURCE_DIR}/base/src/*.inl )

## keep here for auto generating headers if needed
 add_custom_target(headers ./fhead.sh "${FASP_SOURCE_DIR}/base"
		   DEPENDS ${FASP_C_SOURCES}
                   WORKING_DIRECTORY "${FASP_SOURCE_DIR}/util"
                   COMMENT 
		   "FASP: Autogenerating header file with C functions..."
                   VERBATIM)

##################################################################
# For UMFPACK
##################################################################
if (USE_UMFPACK)
  
   # set the path to find specific modules	  
   set(CMAKE_MODULE_PATH "${FASP_SOURCE_DIR}/modules")

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
if (USE_SUPERLU)

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
if (USE_MUMPS)

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
# For Doxygen
##################################################################
# target generating the FASP documentation with doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
    if(EXISTS ${FASP_SOURCE_DIR}/doc/fasp.Doxygen.cnf.in)
    configure_file(
		${FASP_SOURCE_DIR}/doc/fasp.Doxygen.cnf.in
	    ${CMAKE_CURRENT_BINARY_DIR}/fasp.Doxygen.cnf @ONLY
		)
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

# Add include directories.
    include_directories(${FASPLIB_BASE_PATH}/include)
# 
# Recursively look for CMakeLists.txt in subdirs.
    add_subdirectory("base")
    add_subdirectory("test")
    add_subdirectory("tutorial")


##############################below are debugs and tries out
#enable_testing()
#add_test(headers ${FASP_SOURCE_DIR}/util/fhead.sh ${CMAKE_CURRENT_SOURCE_DIR})
## add_custom_command(OUTPUT include/fasp_functs.h
##                   COMMAND ./fhead.sh ARGS ${CMAKE_CURRENT_SOURCE_DIR}
##		   DEPENDS ${FASP_BASE_SOURCES}
##                   WORKING_DIRECTORY "${FASP_SOURCE_DIR}/util"
##                   COMMENT "Autogenerating header files first..."
##                   VERBATIM)


#MISMATCH CHECK:
if(NOT (CMAKE_Fortran_COMPILER_ID STREQUAL CMAKE_C_COMPILER_ID))
  message("
** WARNING: CC and FC mismatch. C compiler=\"${CMAKE_C_COMPILER_ID}${CMAKE_C_COMPILER_VERSION}\" ; Fortran compiler=\"${CMAKE_Fortran_COMPILER_ID}${CMAKE_Fortran_COMPILER_VERSION}\"
"
) 
endif(NOT (CMAKE_Fortran_COMPILER_ID STREQUAL CMAKE_C_COMPILER_ID))

