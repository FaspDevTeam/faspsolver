        if( THE_C AND C_XX AND F_C )
	     find_program(THE_CXX ${C_XX} )
	     find_program(THE_F ${F_C} )
	     if( THE_F AND THE_CXX ) 
	    	 set(CMAKE_CXX_COMPILER ${THE_CXX} CACHE INTERNAL   "the C++ compiler" FORCE)
             set(CMAKE_Fortran_COMPILER ${THE_F} CACHE INTERNAL    "the F compiler" FORCE)
	     endif( THE_F AND THE_CXX )
	endif( THE_C AND C_XX AND F_C )


# END COMPILERS SET UP................ 
