#  
#  Writing configuration so that it can be used by regular makefiles 
#
#  Modified   20151017   --ltz
###################################################################
set(CONFIGMK ${FASP_SOURCE_DIR}/Config.mk)
message("-- Writing FASP configuration to ${CONFIGMK}")
file(WRITE  ${CONFIGMK} 
"
\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\# Automatically generated \#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#
\# Fast Auxiliary Space Preconditioners (FASP) 
\#
\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#
\# This file is rewritten when \"make config\" is run.
\# It is (s)included by test/Makefile and tutorial/Makefile
\#
\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#
fasp_prefix=${FASP_INSTALL_PREFIX}
CC=${CMAKE_C_COMPILER}
CXX=${CMAKE_CXX_COMPILER}
FC=${CMAKE_Fortran_COMPILER}
\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#
")
