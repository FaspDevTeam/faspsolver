# - Try to find CAMD
#  
#  OUTPUT:
#  CAMD_FOUND        - system has CAMD
#  CAMD_INCLUDE_DIRS - include directories for CAMD
#  CAMD_LIBRARIES    - libraries for CAMD
#
#  Xiaozhe Hu 
#  02/27/2013

message(STATUS "Checking for package 'CAMD'")

# Check for header file
find_path(CAMD_INCLUDE_DIRS camd.h
 HINTS ${CAMD_DIR}/include ${CAMD_DIR}/CAMD/include $ENV{CAMD_DIR}/include $ENV{CAMD_DIR}/CAMD/include
 PATH_SUFFIXES suitesparse ufsparse
 DOC "Directory where the CAMD header is located"
 )
mark_as_advanced(CAMD_INCLUDE_DIRS)

# Check for CAMD library
find_library(CAMD_LIBRARIES camd
  HINTS ${CAMD_DIR}/lib ${CAMD_DIR}/CAMD/lib $ENV{CAMD_DIR}/lib $ENV{CAMD_DIR}/CAMD/lib
  DOC "The CAMD library"
  )
mark_as_advanced(CAMD_LIBRARIES)

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CAMD
  "CAMD could not be found. Be sure to set CAMD_DIR."
  CAMD_LIBRARIES CAMD_INCLUDE_DIRS)
