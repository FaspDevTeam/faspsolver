########################################################################
# Fast Auxiliary Space Preconditioners (FASP)
# 
# This is the Makefile for FASP base! 
#
########################################################################

# !!! IMPORTANT !!!
# Do NOT modify this Makefile and put your personal definitions here.
# Use make.inc to define user settings instead. If you really need to
# revise this file, use "chmod u+w Makefiel" first.  

########################################################################
# Location of sources and libraries 
########################################################################
SRCDIR=./src
EXTDIR=./extra ./extra/interface 
EXTDIR+=./extra/sparsekit ./extra/dlmalloc ./extra/nedmalloc
INCLUDE=-I./include
AMGLIB=lib/libfasp.a

########################################################################
# Compilers and dependences
########################################################################
AR=ar ruc
CC=gcc
CPP=g++
FC=gfortran

########################################################################      
# Compiling options                                                             
########################################################################        
BOPT=-O0 -Wall -g -pg

COPTS=$(BOPT)
#CDEFS=-D FASP_USE_OPENMP=0
CINCLUDES=$(INCLUDE)
CFLAGS=$(COPTS) $(CINCLUDES) $(CDEFS)

FOPTS=$(BOPT)
#FDEFS=-D FASP_USE_OPENMP=0
FINCLUDES=$(CINCLUDES)
FFLAGS=$(FOPTS) $(FINCLUDES) $(FDEFS)

########################################################################
# Set libraries
########################################################################

LIBS=$(AMGLIB)

########################################################################
# Load user-defined parameters
########################################################################

include make.inc

########################################################################
# Options for building libraries
########################################################################

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
.C.o:
	$(CPP) -c $< -o $@ $(CFLAGS) 
	$(AR) $(AMGLIB) $@

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
	| awk -v name="fasp_functs.h" -f ../util/mkheaders.awk > include/fasp_functs.h

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
	@-rm -f $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*.o))
	@-rm -f $(foreach DIR,$(SRCDIR),$(wildcard $(DIR)/*~))
	@-rm -f $(foreach DIR,$(EXTDIR),$(wildcard $(DIR)/*.o))
	@-rm -f $(foreach DIR,$(EXTDIR),$(wildcard $(DIR)/*~))
	rm -f lib/libfasp.a *~

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
