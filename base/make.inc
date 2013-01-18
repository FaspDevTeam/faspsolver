# This is the included file by the main Makefile. Please put your localization 
# definitions here. But do not check in your changes to the repository. :D

#==============================================================================#
# User compilers                                                               #
# FASP has been tested with many different compilers; see README for details.  #
# Note: This will overwrite the seting in Makefile!!!                          #
#==============================================================================#
# CC=icc
# CPP=icpc
# FC=ifort

#==============================================================================#
# User compiling options                                                       #
# For reasonable speed, uncomment the following line; you may try -O3 instead  #
# of -O2. Also, "-ffast-math", "-msse2" or "-maltivec" can be applied to allow #
# vectorization.                                                               #
# Note: This will overwrite the seting in Makefile!!!                          #
#==============================================================================#
BOPT=-O2 -funroll-loops -ftree-vectorize

#==============================================================================#
# OpenMP Support                                                               #
# Uncomment the following line to turn on OpenMP version                       #
#==============================================================================#
# BOPT+=-fopenmp

#==============================================================================#
# User preprocessing definitions                                               #
#==============================================================================#
# CDEFS+=
# FDEFS+=

# Rules for compiling the source files
.SUFFIXES: .c .cc .cpp .C .for .f .f77 .f90 .f95
#
FSRC := $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*.for))
FSRC += $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*.f))
FSRC += $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*.f77))
FSRC += $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*.f90))
FSRC += $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*.f95))
CSRC := $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*.c))
CSRC += $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*.cc))
CSRC += $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*.cpp))
CSRC += $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*.C))
FSRC += $(foreach DIR,$(EXTDIR),$(wildcard $(DIR)/*.f))
CSRC += $(foreach DIR,$(EXTDIR),$(wildcard $(DIR)/*.c))
#
OBJSF := $(patsubst %.for,%.o,$(FSRC))
OBJSF += $(patsubst %.f,%.o,$(FSRC))
OBJSF += $(patsubst %.f77,%.o,$(FSRC))
OBJSF += $(patsubst %.f90,%.o,$(FSRC))
OBJSF += $(patsubst %.f95,%.o,$(FSRC))
OBJSC := $(patsubst %.c,%.o,$(CSRC))
OBJSC += $(patsubst %.cc,%.o,$(CSRC))
OBJSC += $(patsubst %.cpp,%.o,$(CSRC))
OBJSC += $(patsubst %.C,%.o,$(CSRC))
