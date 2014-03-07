# - Try to find METIS
# 
#  OUTPUT:
#  METIS_FOUND        - system has METIS
#  METIS_INCLUDE_DIRS - include directories for METIS
#  METIS_LIBRARIES    - libraries for METIS
#
#  Xiaozhe Hu
#  02/27/2013

message(STATUS "Checking for package 'METIS'")

# Check for header file
find_path(METIS_INCLUDE_DIRS metis.h
 HINTS ${METIS_DIR}/include ${METIS_DIR}/lib ${METIS_DIR}/METIS/include ${METIS_DIR}/METIS/lib $ENV{METIS_DIR}/include $ENV{METIS_DIR}/lib $ENV{METIS_DIR}/METIS/include $ENV{METIS_DIR}/METIS/lib
 PATH_SUFFIXES suitesparse ufsparse
 DOC "Directory where the METIS header is located"
 )
mark_as_advanced(METIS_INCLUDE_DIRS)

# Check for METIS library
find_library(METIS_LIBRARIES metis
  HINTS ${METIS_DIR}/include ${METIS_DIR}/lib ${METIS_DIR}/METIS/include ${METIS_DIR}/METIS/lib $ENV{METIS_DIR}/include $ENV{METIS_DIR}/lib $ENV{METIS_DIR}/METIS/include $ENV{METIS_DIR}/METIS/lib
  DOC "The METIS library"
  )
mark_as_advanced(METIS_LIBRARIES)

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(METIS
  "METIS could not be found. Be sure to set METIS_DIR."
  METIS_LIBRARIES METIS_INCLUDE_DIRS)
