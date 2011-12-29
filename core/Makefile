########################################################################
# Fast Auxiliary Space Preconditioners (FASP)
# 
# This is the Makefile for FASP core! 
# Do NOT put your personal definitions in this general file. 
# Use make.inc to define user settings instead.
#
########################################################################

########################################################################
# Compilers and dependences
########################################################################
AR=ar ruc
CC=gcc 
CPP=g++ 
FC=gfortran 

CC=icc -openmp
CPP=icpc -openmp
FC=ifort -openmp
CSRCDIR=./src
FSRCDIR=./src
EXTRDIR= ./extra
INCLUDE=./include
AMGLIB=lib/libfasp.a

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
LIBS=$(AMGLIB)

########################################################################
# Load user-defined parameters
########################################################################
include make.inc

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
FSRC += $(foreach dir,$(EXTRDIR),$(wildcard $(EXTRDIR)/*.f))
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
	$(AR) $(AMGLIB) $@
#
.f.o:
	$(FC) -c $< -o $@ $(FFLAGS) 
	$(AR) $(AMGLIB) $@
#
.f77.o:
	$(FC) -c $< -o $@ $(FFLAGS) 
	$(AR) $(AMGLIB) $@
#
.f90.o:
	$(FC) -c $< -o $@ $(FFLAGS) 
	$(AR) $(AMGLIB) $@
#
.f95.o:
	$(FC) -c $< -o $@ $(FFLAGS) 
	$(AR) $(AMGLIB) $@
#
.c.o:
	$(CC) -c $< -o $@ $(CFLAGS) 
	$(AR) $(AMGLIB) $@
#
.cc.o:
	$(CPP) -c $< -o $@ $(CFLAGS) 
	$(AR) $(AMGLIB) $@
#
.cpp.o:
	$(CPP) -c $< -o $@ $(CFLAGS) 
	$(AR) $(AMGLIB) $@
#
########################################################################
# List of all programs to be compiled
########################################################################

# Everything
ALLPROG=$(AMGLIB)

########################################################################
# Libraries
########################################################################

all: $(ALLPROG)

Default: 
	lib

headers: 
	/bin/cat $(CSRCDIR)/*.c \
	| awk -v name="fasp_functs.h" -f ../bin/mkheaders.awk > $(INCLUDE)/fasp_functs.h

$(AMGLIB): $(OBJSC) $(OBJSF)
	ranlib $(AMGLIB)

lib: $(OBJSC) $(OBJSF)
	ranlib $(AMGLIB)

########################################################################
# Clean up
########################################################################

.PHONY : clean allclean help

clean:
	rm -f $(CSRCDIR)/*.o
	rm -f $(FSRCDIR)/*.o
	rm -f $(EXTRDIR)/*.o

libclean:
	rm -f include/fasp_functs.h
	rm -f lib/libfasp.a

allclean:
	make clean
	rm -f lib/libfasp.a
	rm -f *~ *.ex *.out
	rm -f $(CSRCDIR)/*~
	rm -f $(FSRCDIR)/*~
	rm -f $(EXTRDIR)/*~

help:
	@echo "======================================================"
	@echo " Fast Auxiliary Space Preconditioners (FASP)"
	@echo "======================================================"
	@echo " "
	@echo " make            : build all exe files "
	@echo " make headers    : build the header file automatically"
	@echo " make lib        : build the library "
	@echo " make clean      : clean all obj files "
	@echo " make libclean   : clean all lib files "
	@echo " make allclean   : clean all obj, exe, bak, out files "
	@echo " make help       : show this screen "
	@echo " "
