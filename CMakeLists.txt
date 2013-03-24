	cmake_minimum_required (VERSION 2.6)


# INIT_CONFIG load...
	include ( FASP_init.cmake OPTIONAL )


#
	###################################################

	project(FASP C CXX Fortran )

##       message(" CHECK:  ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}") 
##       message(" CHECK1:  ${CMAKE_C_FLAGS}") 

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
          message(WARNING "OpenMP was requested but not supported!")
       endif(OPENMP_FOUND)
    endif(USE_OPENMP)
############## Additional compiler flags (not defined by the build
################type

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


	set(FASP_MAJOR_VERSION 1)
	set(FASP_MINOR_VERSION 2)
	set(FASP_PATCH_VERSION 1)
	set(FASP_VERSION
		${FASP_MAJOR_VERSION}.${FASP_MINOR_VERSION}.${FASP_PATCH_VERSION})
#
        set(CMAKE_INSTALL_PREFIX "${FASP_SOURCE_DIR}" )
	set(FASPLIB_BASE_PATH "base" CACHE PATH "base path to FASPLIB")
        set(FASP_INSTALL TRUE)
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
   set(SUITESPARSE_CONFIG_DIR "${SUITESPARSE_DIR}")
   set(AMD_DIR "${SUITESPARSE_DIR}")
   set(CHOLMOD_DIR "${SUITESPARSE_DIR}")
   set(COLAMD_DIR "${SUITESPARSE_DIR}")
   set(CAMD_DIR "${SUITESPARSE_DIR}")
   set(METIS_DIR "${SUITESPARSE_DIR}")
   set(CCOLAMD_DIR "${SUITESPARSE_DIR}")
   set(UMFPACK_DIR "${SUITESPARSE_DIR}")   

   find_package(UMFPACK)     
   if (UMFPACK_FOUND)
      add_definitions("-DWITH_UMFPACK=1")
      include_directories(${UMFPACK_INCLUDE_DIRS})
   else(UMFPACK_FOUND)
      message(WARNING "UMFPACK was requested but not supported!")
   endif(UMFPACK_FOUND)
endif(USE_UMFPACK)

############## documentation with doxygen
# target generating the FASP documentation with doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
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
##        message( "DOXY exec is ${DOXY_EXEC}")
    endif(DOXYWIZARD)
    add_custom_target(docs ${DOXY_EXEC}
 	       ${CMAKE_CURRENT_BINARY_DIR}/fasp.Doxygen.cnf
	       WORKING_DIRECTORY 
	       "${CMAKE_CURRENT_BINARY_DIR}"
	       COMMENT 
	       "Generating FASP documentation (Doxygen)" 
	       VERBATIM)
else(DOXYGEN_FOUND)
message(" WARNING: Doxygen was requested but is not found ")
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
