## build test library

set(FASP_TLIB_DIR ${PROJECT_SOURCE_DIR}/test)
set(FASP_TINC_DIR ${PROJECT_SOURCE_DIR}/test/include)
set(RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/test)

## requires math library
if(UNIX OR LINUX)
	set(ADD_STDLIBS m)
endif()

include_directories(${FASP_TINC_DIR}
		    ${PROJECT_SOURCE_DIR}/base/include)

file(GLOB FASP_TLIB_SRC src/*.c src/*.f)

# Build libfasptest as STATIC and do not install it. If we want
# SHARED, then we need to figure out how apple works with multiple
# rpaths put the shared lib in the same dir as the executables.

add_library(fasptest STATIC ${FASP_TLIB_SRC}) 

if(MSVC)
    file( GLOB FASP_TSRC RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/main" 
          main/*.c )
else()
    file( GLOB FASP_TSRC RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/main" 
          main/*.c main/*.cpp main/*.f main/*.f90 )
endif()

foreach( srcfile ${FASP_TSRC} )
    string( REGEX REPLACE "\\.[^\\.]*$" ".ex" execfile ${srcfile})
    add_executable( ${execfile} main/${srcfile})
# 
    target_link_libraries(${execfile} fasptest 
                          ${FASP_LIBRARY_NAME}
                          ${UMFPACK_LIBRARIES} 
                          ${SUPERLU_LIBRARIES} 
                          ${MUMPS_LIBRARIES}
                          ${MKL_LIBRARIES}
                          ${ADD_STDLIBS})
#
    install(TARGETS ${execfile} 
            DESTINATION ${PROJECT_SOURCE_DIR}/test)
endforeach( srcfile ${FASP_TSRC} )