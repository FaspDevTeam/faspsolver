set(REAL_C "${CMAKE_C_COMPILER_ID}${CMAKE_C_COMPILER_VERSION}")
###################################################XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
if(${REAL_C} MATCHES "GNU.*" AND ${THE_C} MATCHES "gcc.*") 
    string(REPLACE "gcc" "g++" C_XX ${THE_C} )
	string(REPLACE "gcc" "gfortran" F_C ${THE_C} )
	find_program(THE_CXX NAMES ${C_XX})
	find_program(THE_FC NAMES ${F_C})
elseif( ${REAL_C} MATCHES "Intel.*" AND ${THE_C} MATCHES "icc.*" )
    message("${THE_C} and ${REAL_C} is 3")
    find_program(THE_CXX NAMES icpc)
    find_program(THE_FC NAMES ifort)
elseif( ${REAL_C} MATCHES "Clang.*")
	find_program(THE_C NAMES clang)
  	find_program(THE_CXX NAMES clang++ clang)
		# clang seems to work with gnu compilers if
		# dragonegg plugin is installed; it is not clear how
		# to include this dragonegg plugin on Linux in the
		# CMAKEs, so we will issue a warning at the end if the 
		# C compiler and fortran compilers are different. 
		# As compiler gfortran works with clang in many
		# cases, but it might fail on Apple M1 chips
    find_program(THE_FC NAMES $ENV{FC} gfortran g95 g77)
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
###################################################XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
enable_language(CXX) 
enable_language(Fortran) 
###################################################XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
