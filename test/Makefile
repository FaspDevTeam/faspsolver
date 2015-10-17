#######################################################################
# Fast Auxiliary Space Preconditioners (FASP)
#
########################################################################
# Makefile for building the test suite. If the file ../Config.mk
# exists it will take the values of configuration variables from
# there. Otherwise, they have to be user supplied below.

fasp_prefix = not-defined-yet
CC = not-defined-yet
FC = not-defined-yet
CXX = not-defined-yet

sinclude ../Config.mk

ifeq ($(fasp_prefix),not-defined-yet)
    fasp_prefix = ..
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

#no fflags here
#FFLAGS= -I$(fasp_prefix)/include/ #-I./include
#FFLAGS+=-O3

#LCFLAGS=$(LINK_OPT) -L$(fasp_prefix)/lib -lgfortran -lfasp 
LINKER = $(FC)
LFLAGS = $(LINK_OPT) -L$(fasp_prefix)/lib -lfasp 

DSRC=$(shell find ./src -name '*.c')
DOBJ=$(DSRC:.c=.o)

.PHONY: all clean distclean

all:	regression.ex regression_mf.ex regression_mm.ex \
	test.ex testbamg.ex testbcsr.ex testbsr.ex testgmg.ex \
	testlab.ex testmat.ex testsmat.ex \
	testrap.ex testfem.ex testheat.ex

testrap.ex: main/testrap.o $(DOBJ)
	@$(LINKER) main/testrap.o $(DOBJ) -o $@ $(LFLAGS) 
	@echo 'Building executable file $@' ... Done

testfem.ex: main/testfem.o $(DOBJ)
	@$(LINKER) main/testfem.o $(DOBJ) -o $@ $(LFLAGS) 
	@echo 'Building executable file $@' ... Done

testheat.ex: main/testheat.o $(DOBJ)
	@$(LINKER) main/testheat.o $(DOBJ) -o $@ $(LFLAGS) 
	@echo 'Building executable file $@' ... Done

%.o: %.c
	@$(CC) $(CFLAGS) -c -o $@ $<
#	@echo 'Compiling source file $<' ... Done

%.ex: main/%.c
	@$(CC) -c $(CFLAGS) -o $@.o $< 
#	@echo 'Compiling source file $<' ... Done
	@$(LINKER) $@.o -o $@ $(LFLAGS) 
	@echo 'Building executable file $@' ... Done

clean:
		+find . -name '*.o' -exec rm {} \;

distclean:	clean
		+rm -f *.ex *~
