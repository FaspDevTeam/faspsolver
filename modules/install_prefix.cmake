## - Try to see if the prefix is OK. 
#  
#  FASP_INSTALL_PREFIX is set here. 
#
#  Modified   20151016   --ltz
if(DEFINED FASP_INSTALL_PREFIX)
  string(LENGTH "${FASP_INSTALL_PREFIX}" Z)
  if(Z)
    get_filename_component(FASP_INSTALL_PREFIX ${FASP_INSTALL_PREFIX} REALPATH)
  else(Z)
   unset(FASP_INSTALL_PREFIX)
  endif(Z)
endif(DEFINED FASP_INSTALL_PREFIX)
####
if(DEFINED FASP_INSTALL_PREFIX)
if((EXISTS ${FASP_INSTALL_PREFIX}) AND (IS_DIRECTORY ${FASP_INSTALL_PREFIX}))
  else((EXISTS ${FASP_INSTALL_PREFIX}) AND (IS_DIRECTORY ${FASP_INSTALL_PREFIX}))
  message("
** WARNING: Installation dir=\"${FASP_INSTALL_PREFIX}\" does not exist or is not a directory ; 
   	    Configured install in the default location \"${FASP_SOURCE_DIR}\"
"
) 
 set(FASP_INSTALL_PREFIX ${FASP_SOURCE_DIR})
 endif((EXISTS ${FASP_INSTALL_PREFIX})  AND  (IS_DIRECTORY ${FASP_INSTALL_PREFIX}))
else(DEFINED FASP_INSTALL_PREFIX)
#  message("
#** WARNING: Installation prefix \"${FASP_INSTALL_PREFIX}\"
#   does not seem to be a valid name; 
#**          Configured install in the default location \"${FASP_SOURCE_DIR}\"
#"
#) 
   set(FASP_INSTALL_PREFIX ${FASP_SOURCE_DIR})
endif(DEFINED FASP_INSTALL_PREFIX)
   message( "Info: FASP library/headers installation dir =\"${FASP_INSTALL_PREFIX}\"")
