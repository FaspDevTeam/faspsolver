# - Try to find COLAMD
# 
#  OUTPUT:
#  COLAMD_FOUND        - system has COLAMD
#  COLAMD_INCLUDE_DIRS - include directories for COLAMD
#  COLAMD_LIBRARIES    - libraries for COLAMD
#
#  Xiaozhe Hu
#  02/27/2013

message(STATUS "Checking for package 'COLAMD'")

# Check for header file
find_path(COLAMD_INCLUDE_DIRS colamd.h
 HINTS ${COLAMD_DIR}/include ${COLAMD_DIR}/COLAMD/include $ENV{COLAMD_DIR}/include $ENV{COLAMD_DIR}/COLAMD/include
 PATH_SUFFIXES suitesparse ufsparse
 DOC "Directory where the COLAMD header is located"
 )
mark_as_advanced(COLAMD_INCLUDE_DIRS)

# Check for COLAMD library
find_library(COLAMD_LIBRARIES colamd
  HINTS ${COLAMD_DIR}/lib ${COLAMD_DIR}/COLAMD/lib $ENV{COLAMD_DIR}/lib $ENV{COLAMD_DIR}/COLAMD/lib
  DOC "The COLAMD library"
  )
mark_as_advanced(COLAMD_LIBRARIES)

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(COLAMD
  "COLAMD could not be found. Be sure to set COLAMD_DIR."
  COLAMD_LIBRARIES COLAMD_INCLUDE_DIRS)
