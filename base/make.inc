# This is the included file for the main Makefile -- put your localization definitions here.
# But do not check in your changes to the repository. :D

# User compilers
# CC=icc
# CPP=icpc
# FC=ifort

# User compiling options
BOPT=-g -O2 -funroll-loops -ftree-vectorize

# Uncomment the following line to turn on OpenMP version
   BOPT+=-fopenmp
#OMP=1

##iifdef OMP
#   BOPT+=-fopenmp
#   CDEFS=-D FASP_USE_OPENMP=1
#   FDEFS=-D FASP_USE_OPENMP=1
##ielse
#   CDEFS=-D FASP_USE_OPENMP=0
#   FDEFS=-D FASP_USE_OPENMP=0
#endif

# User preprocessing definitions
# CDEFS+=xx
# FDEFS+=xx
