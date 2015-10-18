#######################################################################
# Fast Auxiliary Space Preconditioners (FASP)
#
########################################################################
# Makefile for building the test suite. If the file ../Config.mk
# exists it will take the values of configuration variables from
# there. Otherwise, they have to be user supplied below.
########################################################################

fasp_prefix = not-defined-yet
fasp_library = not-defined-yet
CC = not-defined-yet
FC = not-defined-yet
CXX = not-defined-yet

# include the configuration written by CMake at config time if found
sinclude ../Config.mk

ifeq ($(fasp_prefix),not-defined-yet)
    fasp_prefix = ..
endif

ifeq ($(fasp_library),not-defined-yet)
    fasp_library = libfasp.a
endif

ifeq ($(CC),not-defined-yet) 
   CC=gcc 
endif

ifeq ($(FC),not-defined-yet)
   FC=gfortran
endif

ifeq ($(CXX),not-defined-yet)
   CXX=g++
endif

CFLAGS = -I$(fasp_prefix)/include/ -I./include 
CFLAGS+=-O3
LINKER = $(FC)
LFLAGS = $(LINK_OPT) -L$(fasp_prefix)/lib -lfasp 

fasp_lib_file=$(fasp_prefix)/lib/$(fasp_library) 

DSRC=$(shell find ./src -name '*.c')
DOBJ=$(DSRC:.c=.o)

.PHONY: all clean distclean

all:	regression.ex regression_mf.ex regression_mm.ex \
	test.ex testbamg.ex testbcsr.ex testbsr.ex testgmg.ex \
	testlab.ex testmat.ex testsmat.ex \
	testrap.ex testfem.ex testheat.ex

testrap.ex: main/testrap.o $(DOBJ) $(fasp_lib_file)
	@$(LINKER) main/testrap.o $(DOBJ) -o $@ $(LFLAGS) 

testfem.ex: main/testfem.o $(DOBJ) $(fasp_lib_file)
	@$(LINKER) main/testfem.o $(DOBJ) -o $@ $(LFLAGS) 

testheat.ex: main/testheat.o $(DOBJ) $(fasp_lib_file)
	@$(LINKER) main/testheat.o $(DOBJ) -o $@ $(LFLAGS) 

%.ex: main/%.c $(fasp_lib_file)
	@$(CC) -c $(CFLAGS) -o main/$@.o $< 
	@$(LINKER) main/$@.o -o $@ $(LFLAGS) 
	@echo 'Building executable file $@' ... Done

%.o: %.c $(fasp_lib_file)
	@$(CC) $(CFLAGS) -c -o $@ $<

$(fasp_lib_file):	
	$(error The FASP library $@ is not found)

clean:
	@+find . -name '*.o' -exec rm {} \;

distclean: clean
	@+rm -f *.ex *~
