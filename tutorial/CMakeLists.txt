## build tutorial programs

file( GLOB FASP_TUTSRC RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/main" main/*.c main/*.f main/*.f90 )
#message("FASP_TUTSRC=${FASP_TUTSRC} ")

foreach( srcfile ${FASP_TUTSRC} )
 	 string( REGEX MATCH "\\.c$|\\.f90$|\\.cpp$" X ${srcfile} )
#        message( "X=${X} ${srcfile} ")
    if("${X}" STREQUAL ".c")
        string( REGEX REPLACE "\\.c$" "-c.ex" execfile ${srcfile} )
    elseif("${X}" STREQUAL ".f90")
        string( REGEX REPLACE "\\.f90$" "-f.ex" execfile ${srcfile} )
    elseif("${X}" STREQUAL ".cpp")
        string( REGEX REPLACE "\\.cpp$" "-cpp.ex" execfile ${srcfile} )
    endif()
    if(execfile)
      add_executable( ${execfile} main/${srcfile} )
      target_link_libraries( ${execfile} fasp ${UMFPACK_LIBRARIES}) 
      install(TARGETS ${execfile} DESTINATION ${FASP_SOURCE_DIR}/tutorial)
#      message( "${execfile} ${srcfile} ")
    endif(execfile)
endforeach( srcfile ${FASP_TUTSRC} )

