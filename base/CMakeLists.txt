#
set(FASP_BASE_SRC_DIR "src")
#
set(FASP_EXTRA_SRC_DIR "extra")
#
file(GLOB FASP_EXTRA_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}  
	  		     	  ${FASP_EXTRA_SRC_DIR}/sparsekit/ilu.f 
				      ${FASP_EXTRA_SRC_DIR}/interface/interface_umfpack.c
				      ${FASP_EXTRA_SRC_DIR}/interface/interface_superlu.c
				      ${FASP_EXTRA_SRC_DIR}/interface/interface_mumps.c)
#
file(GLOB FASP_BASE_HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
	  		    	     ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h)

file(GLOB FASP_BASE_INL     RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
	  		    	     ${CMAKE_CURRENT_SOURCE_DIR}/src/*.inl)
#
file(GLOB FASP_BASE_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
	  		    	     ${FASP_BASE_SRC_DIR}/*.c 
	  		    	     ${FASP_BASE_SRC_DIR}/*.f 
	  		    	     ${FASP_EXTRA_SOURCES} )
#
##
# Build libfasp
# depends also on the source, header and inl files 
add_library(fasp ${FASP_LIBRARY_TYPE} ${FASP_BASE_SOURCES} ${FASP_BASE_INL} 
	  		    	     ${FASP_BASE_HEADERS}) 

# install libfasp
    install(TARGETS fasp
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION lib
    ARCHIVE DESTINATION lib)

