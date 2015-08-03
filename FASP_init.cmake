# Helper modules.
include(CheckFunctionExists)
include(CheckIncludeFile)

set(CMAKE_VERBOSE_MAKEFILE 1) 
set(GDB 1 CACHE BOOL "debugging or not")
set(OPENMP 0 CACHE BOOL "Openmp use")
set(USE_MUMPS 0 CACHE BOOL "MUMPS use")

# For which compilers we shall search (if none found then the 
# default cmake compiler detection will be invoked. 
	set(F_C "dragonegg-3.4-gfortran-mp-4.6")
        if( THE_C AND C_XX AND F_C )
	     find_program(THE_CXX ${C_XX} )
	     find_program(THE_F ${F_C} )
	     if( THE_F AND THE_CXX ) 
	    	 set(CMAKE_CXX_COMPILER ${THE_CXX} CACHE INTERNAL   "the C++ compiler" FORCE)
             set(CMAKE_Fortran_COMPILER ${THE_F} CACHE INTERNAL    "the F compiler" FORCE)
	     endif( THE_F AND THE_CXX )
	endif( THE_C AND C_XX AND F_C )


# END COMPILERS SET UP................ 
