	cmake_minimum_required (VERSION 2.8)


# INIT_CONFIG load...
	include(fasp_init_config.cmake OPTIONAL )

#
	###################################################
#
	project(FASP C CXX Fortran )
#
# OpenMP : defined on command line in the top Makefile
#  if(OPENMP)
#	include(FindOpenMP)
#	if(OPENMP_FOUND)
#		message("OPENMP FOUND")
#		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} 
#			   ${OpenMP_EXE_LINKER_FLAGS}")
#	else()
#	       message(WARNING "OpenMP was requested but support was not found")
#        endif(OPENMP_FOUND)
#    endif(OPENMP)
#

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

############## documentation with doxygen
# target generating the FASP documentation with doxygen
find_package(Doxygen)
##message(" WAS doxygen found: ${DOXYGEN_FOUND} ")
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
        message( "DOXY exec is ${DOXY_EXEC}")
    endif(DOXYWIZARD)
    add_custom_target(docs ${DOXY_EXEC}
 	       ${CMAKE_CURRENT_BINARY_DIR}/fasp.Doxygen.cnf
	       WORKING_DIRECTORY 
	       "${CMAKE_CURRENT_BINARY_DIR}"
	       COMMENT 
	       "Generating FASP documentation (Doxygen)" 
	       VERBATIM)
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
