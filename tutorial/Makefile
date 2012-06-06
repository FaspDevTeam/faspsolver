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

FASPDIR=../../faspsolver
SRCDIR=./src
INCLUDE=-I${FASPDIR}/core/include
FASPLIB=${FASPDIR}/core/lib/libfasp.a

########################################################################      
# Compiling options                                                             
########################################################################        
BOPT=-g -O0 -Wall

COPTS=$(BOPT)
CDEFS=
CINCLUDES=$(INCLUDE)
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
include ${FASPDIR}/core/make.inc

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
FSRC := $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*.for))
FSRC += $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*.f))
FSRC += $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*.f77))
FSRC += $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*.f90))
FSRC += $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*.f95))
CSRC := $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*.c))
CSRC += $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*.cc))
CSRC += $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*.cpp))
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
ALLPROG=$(TESTLIB)
ALLPROG+=poisson-amg.c poisson-amg.f90
ALLPROG+=poisson-its.c
ALLPROG+=poisson-pcg.c poisson-pcg.f90
ALLPROG+=spe01-its.c

# Test for solvers
TESTPROG=$(TESTLIB) main/poisson-amg.o 

########################################################################
# Link
########################################################################

all: $(ALLPROG)

Default: 
	poisson-amg.c

########################################################################
# Some test problems
########################################################################

poisson-amg.c:
	$(CC) -c main/poisson-amg.c -o main/poisson-amg.o $(CFLAGS)
	$(FC) -o poisson-amg-c.ex main/poisson-amg.o $(FLFLAGS)

poisson-amg.f90:
	$(FC) -c main/poisson-amg.f90 -o main/poisson-amg.o $(FFLAGS)
	$(FC) -o poisson-amg-f.ex main/poisson-amg.o $(FLFLAGS)

poisson-its.c:
	$(CC) -c main/poisson-its.c -o main/poisson-its.o $(CFLAGS)
	$(FC) -o poisson-its-c.ex main/poisson-its.o $(FLFLAGS)

poisson-pcg.c:
	$(CC) -c main/poisson-pcg.c -o main/poisson-pcg.o $(CFLAGS)
	$(FC) -o poisson-pcg-c.ex main/poisson-pcg.o $(FLFLAGS)

poisson-pcg.f90:
	$(FC) -c main/poisson-pcg.f90 -o main/poisson-pcg.o $(FFLAGS)
	$(FC) -o poisson-pcg-f.ex main/poisson-pcg.o $(FLFLAGS)

spe01-its.c:
	$(CC) -c main/spe01-its.c -o main/spe01-its.o $(CFLAGS)
	$(FC) -o spe01-its-c.ex main/spe01-its.o $(FLFLAGS)

########################################################################
# Clean up
########################################################################

.PHONY : clean allclean help

clean:
	rm -f $(SRCDIR)/*.o
	rm -f main/*.o

allclean:
	make clean
	rm -f lib/*.a
	rm -f *~ *.ex *.out
	rm -f $(SRCDIR)/*~
	rm -f main/*~

help:
	@echo "========================================================="
	@echo "||     Fast Auxiliary Space Preconditioners (FASP)     ||"
	@echo "========================================================="
	@echo "                                                         "
	@echo " make                 : build all examples               "
	@echo " make poisson-amg.c   : build the C version AMG example  "
	@echo " make poisson-amg.f90 : build the F version AMG example  "
	@echo " make poisson-its.c   : build the C version ITS example  "
	@echo " make poisson-pcg.c   : build the C version PCG example  "
	@echo " make poisson-pcg.f90 : build the F version PCG example  "
	@echo " make clean           : clean all obj files              "
	@echo " make allclean        : clean all obj, exe, out files    "
	@echo " make help            : show this screen                 "
	@echo "                                                         "
