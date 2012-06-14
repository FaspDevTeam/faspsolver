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
INCLUDE=-I$(FASPDIR)/core/include -I./include
FASPLIB=$(FASPDIR)/core/lib/libfasp.a
TESTLIB=lib/libfasptest.a

########################################################################      
# Compiling options                                                             
########################################################################        
BOPT=-g -pg -O0 -Wall

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
ALLPROG=$(TESTLIB) test

# Regression Tests
REGRESSION=$(TESTLIB) main/regression.o 

# Test for solvers
TESTPROG=$(TESTLIB) main/test.o 

# Test for matrix properties
TESTMATPROG=$(TESTLIB) main/testmat.o

########################################################################
# Link
########################################################################

all: $(ALLPROG)

Default: 
	test 

$(TESTLIB): $(OBJSC) $(OBJSF)
	ranlib $(TESTLIB)

lib: $(OBJSC) $(OBJSF)
	ranlib $(TESTLIB)

########################################################################
# Some test problems
########################################################################

test:
	$(CC) -c main/test.c -o main/test.o $(CFLAGS)
	$(FC) -o test.ex main/test.o $(FLFLAGS)

regression: 
	$(CC) -c main/regression.c -o main/regression.o $(CFLAGS)
	$(FC) -o regression.ex main/regression.o $(FLFLAGS)

testmm:
	$(CC) -c main/testmm.c -o main/testmm.o $(CFLAGS)
	$(FC) -o testmm.ex main/testmm.o $(FLFLAGS)

testf: 
	$(FC) -c main/testf.f90 -o main/testf.o $(CFLAGS)
	$(FC) -o testf.ex main/testf.o $(FLFLAGS)

testbsr:
	$(CC) -c main/testbsr.c -o main/testbsr.o $(CFLAGS)
	$(FC) -o testbsr.ex main/testbsr.o $(FLFLAGS)

testmat: 
	$(CC) -c main/testmat.c -o main/testmat.o $(CFLAGS)
	$(FC) -o testmat.ex main/testmat.o $(FLFLAGS)

testrap: 
	$(CC) -c main/testrap.c -o main/testrap.o $(CFLAGS)
	$(FC) -o testrap.ex main/testrap.o $(FLFLAGS)

testfem:
	$(CC) -c main/testfem.c -o main/testfem.o $(CFLAGS)
	$(FC) -o testfem.ex main/testfem.o $(FLFLAGS)

testheat:
	$(CC) -c main/testheat.c -o main/testheat.o $(CFLAGS)
	$(FC) -o testheat.ex main/testheat.o $(FLFLAGS)

testfdm2d :
	${CPP} -c main/testfdm2d.cpp -o main/testfdm2d.o $(CFLAGS)
	$(CPP) -o testfdm2d.ex main/testfdm2d.o $(FLFLAGS)

testfdm3d :
	${CPP} -c main/testfdm3d.cpp -o main/testfdm3d.o $(CFLAGS)
	$(CPP) -o testfdm3d.ex main/testfdm3d.o $(FLFLAGS)

testomp: 
	$(CC) -c main/testomp.c -o main/testomp.o $(CFLAGS)
	$(FC) -o testomp.ex main/testomp.o $(FLFLAGS)

test_rsamg_omp: 
	$(CC) -c main/test_rsamg_omp.c -o main/test_rsamg_omp.o $(CFLAGS)
	$(FC) -o test_rsamg_omp.ex main/test_rsamg_omp.o $(FLFLAGS)

fasp: 
	$(CC) -c main/fasp.c -o main/fasp.o $(CFLAGS)
	$(FC) -o fasp.ex main/fasp.o $(FLFLAGS)

testall:
	make 
	make regression
	make testmm
	make testf
	make testbsr
	make testmat
	make testrap
	make testfem
	make testfdm2d
	make testfdm3d
	make testheat
	make testomp
	make test_rsamg_omp
	nake fasp

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
	rm -f out/mat_*
	rm -f out/rhs_*
	rm -f out/sol_*
	rm -f out/sp_*
	rm -f out/matrix.*

help:
	@echo "======================================================"
	@echo "||   Fast Auxiliary Space Preconditioners (FASP)    ||"
	@echo "======================================================"
	@echo "                                                      "
	@echo " make            : build the lib and simple test      "
	@echo " make test       : build the simple test              "
	@echo " make regression : build the regression test          "
	@echo " make testf      : build the F90 interface test       "
	@echo " make testfem    : build the FEM assembling test      "
	@echo " make testfdm2d  : build the FDM 2D assembling test   "
	@echo " make testfdm3d  : build the FDM 3D assembling test   "
	@echo " make testheat   : build the FEM matrix for Heat eqn  "
	@echo " make testmat    : build the matrix checking test     "
	@echo " make testrap    : build the comparison of RAP        "
	@echo " make testall    : build all test cases               "
	@echo " make clean      : clean all obj files                "
	@echo " make allclean   : clean all obj, bak, out, exe files "
	@echo " make help       : show this help screen              "
	@echo "                                                      "
