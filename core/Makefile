########################################################################
# Fast Auxiliary Space Preconditioners (FASP)
# 
# This is the Makefile for FASP core! 
#
########################################################################

# !!! IMPORTANT !!!
# Do NOT modify this Makefile and put your personal definitions here.
# Use make.inc to define user settings instead. If you really need to
# revise this file, use "chmod u+w Makefiel" first.  

########################################################################
# Compilers and dependences
########################################################################
AR=ar ruc

CC=gcc 
CPP=g++
FC=gfortran
#CC=icc
#CPP=icpc
#FC=ifort

SRCDIR=./src
EXTDIR=./extra ./extra/interface ./extra/sparsekit ./extra/dlmalloc ./extra/nedmalloc
INCLUDE=-I./include
AMGLIB=lib/libfasp.a

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
LIBS=$(AMGLIB)

########################################################################
# Load user-defined parameters
########################################################################
include make.inc

########################################################################
# Rules for compiling the source files
########################################################################
.SUFFIXES: .for .f .f77 .f90 .f95 .c .cc .cpp
#
FSRC := $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*.for))
FSRC += $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*.f))
FSRC += $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*.f77))
FSRC += $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*.f90))
FSRC += $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*.f95))
FSRC += $(foreach DIR,$(EXTDIR),$(wildcard $(DIR)/*.f))
CSRC := $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*.c))
CSRC += $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*.cc))
CSRC += $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*.cpp))
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
	/bin/cat $(SRCDIR)/*.c \
	| awk -v name="fasp_functs.h" -f ../bin/mkheaders.awk > include/fasp_functs.h

$(AMGLIB): $(OBJSC) $(OBJSF)
	ranlib $(AMGLIB)

lib: $(OBJSC) $(OBJSF)
	ranlib $(AMGLIB)

########################################################################
# Clean up
########################################################################

.PHONY: clean allclean help

clean:
	@-rm -f $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*.o))
	@-rm -f $(foreach DIR,$(EXTDIR),$(wildcard $(DIR)/*.o))

libclean:
	rm -f include/fasp_functs.h
	rm -f lib/libfasp.a

allclean:
	make clean
	@-rm -f $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*~))
	@-rm -f $(foreach DIR,$(EXTDIR),$(wildcard $(DIR)/*~))
	rm -f lib/libfasp.a
	rm -f *~ *.ex *.out

help:
	@echo "======================================================"
	@echo "||   Fast Auxiliary Space Preconditioners (FASP)    ||"
	@echo "======================================================"
	@echo "                                                      "
	@echo " make            : build all exe files                "
	@echo " make lib        : build the library                  "
	@echo " make headers    : build the header file automatically"
	@echo " make clean      : clean all obj files                "
	@echo " make libclean   : clean all lib files                "
	@echo " make allclean   : clean all obj, exe, bak, out files "
	@echo " make help       : show this help screen              "
	@echo "                                                      "
