########################################################################
# Fast Auxiliary Space Preconditioners (FASP)
# 
# This is the Makefile for FASP test suite! 
# Do NOT put your personal definitions in this general file. 
# Use ../core/make.inc to define user settings instead.
#
########################################################################

########################################################################
# Compilers and dependences
########################################################################
AR=ar ruc
CC=gcc
CPP=g++
FC=gfortran

CSRCDIR=./src
FSRCDIR=./src
INCLUDE=../core/include
FASPLIB=../core/lib/libfasp.a

########################################################################      
# Compiling options                                                             
########################################################################        
BOPT=-g -pg -O0 -Wall

COPTS=$(BOPT)
CDEFS=
CINCLUDES=-I$(INCLUDE)
CFLAGS=$(CDEFS) $(COPTS) $(CINCLUDES)

FOPTS=$(BOPT)
FDEFS=$(CDEFS)
FINCLUDES=$(CINCLUDES)
FFLAGS=$(FDEFS) $(FOPTS) $(FINCLUDES)

########################################################################
# Set libraries
########################################################################
LIBS=$(TESTLIB) $(FASPLIB)

########################################################################
# Load user-defined parameters
########################################################################
include ../core/make.inc

########################################################################
# Link options
########################################################################
LINKOPTS=$(BOPT)
CLFLAGS=-lstdc++ $(LINKOPTS) $(LIBS)
FLFLAGS=-lm $(LINKOPTS) $(LIBS)

########################################################################
# Rules for compiling the source files
########################################################################
.SUFFIXES: .c .cc .cpp .for .f .f77 .f90 .f95
#
FSRC := $(foreach dir,$(FSRCDIR),$(wildcard $(FSRCDIR)/*.for))
FSRC += $(foreach dir,$(FSRCDIR),$(wildcard $(FSRCDIR)/*.f))
FSRC += $(foreach dir,$(FSRCDIR),$(wildcard $(FSRCDIR)/*.f77))
FSRC += $(foreach dir,$(FSRCDIR),$(wildcard $(FSRCDIR)/*.f90))
FSRC += $(foreach dir,$(FSRCDIR),$(wildcard $(FSRCDIR)/*.f95))
CSRC := $(foreach dir,$(CSRCDIR),$(wildcard $(CSRCDIR)/*.c))
CSRC += $(foreach dir,$(EXTRDIR),$(wildcard $(EXTRDIR)/*.c))
#
OBJSF := $(patsubst %.for,%.o,$(FSRC))
OBJSF += $(patsubst %.f,%.o,$(FSRC))
OBJSF += $(patsubst %.f77,%.o,$(FSRC))
OBJSF += $(patsubst %.f90,%.o,$(FSRC))
OBJSF += $(patsubst %.f95,%.o,$(FSRC))
OBJSC := $(patsubst %.c,%.o,$(CSRC))
#
.for.o:
	$(FC) -c $< -o $@ $(FFLAGS) 
	$(AR) $(TESTLIB) $@
#
.f.o:
	$(FC) -c $< -o $@ $(FFLAGS) 
	$(AR) $(TESTLIB) $@
#
.f77.o:
	$(FC) -c $< -o $@ $(FFLAGS) 
	$(AR) $(TESTLIB) $@
#
.f90.o:
	$(FC) -c $< -o $@ $(FFLAGS) 
	$(AR) $(TESTLIB) $@
#
.f95.o:
	$(FC) -c $< -o $@ $(FFLAGS) 
	$(AR) $(TESTLIB) $@
#
.c.o:
	$(CC) -c $< -o $@ $(CFLAGS) 
	$(AR) $(TESTLIB) $@
#
.cc.o:
	$(CPP) -c $< -o $@ $(CFLAGS) 
	$(AR) $(TESTLIB) $@
#
.cpp.o:
	$(CPP) -c $< -o $@ $(CFLAGS) 
	$(AR) $(TESTLIB) $@
#
########################################################################
# List of all programs to be compiled
########################################################################

# Everything
ALLPROG=$(TESTLIB) poisson-amg

# Test for solvers
TESTPROG=$(TESTLIB) main/poisson-amg.o 

########################################################################
# Link
########################################################################

all: $(ALLPROG)

Default: 
	poisson-amg

#$(TESTLIB): $(OBJSC) $(OBJSF)
#	ranlib $(TESTLIB)

########################################################################
# Some test problems
########################################################################

poisson-amg: 
	$(CC) $(CFLAGS) -c main/poisson-amg.c -o main/poisson-amg.o
	$(FC) -o poisson-amg.ex main/poisson-amg.o $(FLFLAGS)

########################################################################
# Clean up
########################################################################

.PHONY : clean allclean help

clean:
	rm -f $(CSRCDIR)/*.o
	rm -f $(FSRCDIR)/*.o
	rm -f main/*.o

allclean:
	make clean
	rm -f lib/*.a
	rm -f *~ *.ex *.out
	rm -f $(CSRCDIR)/*~
	rm -f $(FSRCDIR)/*~
	rm -f main/*~

help:
	@echo "======================================================"
	@echo " Fast Auxiliary Space Preconditioners (FASP)"
	@echo "======================================================"
	@echo " "
	@echo " make             : build all examples "
	@echo " make poisson-amg : build the AMG solver for Poisson "
	@echo " make clean       : clean all obj files "
	@echo " make allclean    : clean all obj, exe, bak, out files "
	@echo " make help        : show this screen "
	@echo " "
