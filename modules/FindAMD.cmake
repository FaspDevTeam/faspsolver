# - Try to find AMD
# 
#  OUTPUT:
#  AMD_FOUND        - system has AMD
#  AMD_INCLUDE_DIRS - include directories for AMD
#  AMD_LIBRARIES    - libraries for AMD
#
#  Xiaozhe Hu
#  02/27/2013

message(STATUS "Checking for package 'AMD'")

# Check for header file
find_path(AMD_INCLUDE_DIRS amd.h
 HINTS ${AMD_DIR}/include ${AMD_DIR}/AMD/include $ENV{AMD_DIR}/include $ENV{AMD_DIR}/AMD/include
 PATH_SUFFIXES suitesparse ufsparse
 DOC "Directory where the AMD header is located"
 )
mark_as_advanced(AMD_INCLUDE_DIRS)

# Check for AMD library
find_library(AMD_LIBRARIES amd
  HINTS ${AMD_DIR}/lib ${AMD_DIR}/AMD/lib $ENV{AMD_DIR}/lib $ENV{AMD_DIR}/AMD/lib
  DOC "The AMD library"
  )
mark_as_advanced(AMD_LIBRARIES)

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(AMD
  "AMD could not be found. Be sure to set AMD_DIR."
  AMD_LIBRARIES AMD_INCLUDE_DIRS)
