#######################################################################
# Fast Auxiliary Space Preconditioners (FASP)
# 
########################################################################
# Set base variables for building the tutorial executables. If the
# file ../Config.mk exists it will take the values of these variables
# from there.

fasp_prefix = not-set
CC = not-set
FC = not-set
CXX = not-set

sinclude ../Config.mk

ifeq ($(fasp_prefix),not-set)
    fasp_prefix = ..
endif

ifeq ($(CC),not-set) 
   CC=gcc 
endif

ifeq ($(FC),not-set)
   FC=gfortran
endif

ifeq ($(CXX),not-set)
   CXX=g++
endif

CFLAGS = -I$(fasp_prefix)/include/ #-I./include 
CFLAGS+=-O3

FFLAGS= -I$(fasp_prefix)/include/ #-I./include
FFLAGS+=-O3

#LCFLAGS=$(LINK_OPT) -L$(fasp_prefix)/lib -lgfortran -lfasp 
LINKER = $(FC)
LFLAGS = $(LINK_OPT) -L$(fasp_prefix)/lib -lfasp 


examples = poisson-amg-c.ex poisson-its-c.ex poisson-pcg-c.ex \
           poisson-gmg-c.ex spe01-its-c.ex \
           poisson-amg-f.ex poisson-pcg-f.ex

examples_f = $(filter %-f.ex,$(examples))
examples_c = $(filter-out %-f.ex,$(examples))

.PHONY:	all clean

all:	$(examples_c) $(examples_f)

%-c.ex:	main/%.c
	@$(CC) -c $(CFLAGS) -o $@.o $<
#	@echo 'Compiling source file $<' ... Done
	@$(LINKER) -o $@ $@.o $(LFLAGS)
	@echo 'Building executable file $@' ... Done

%-f.ex:	main/%.f90	
	@$(FC) $(FFLAGS) $(LFLAGS) -o $@ $< 


clean:
		+rm -f *.o main/*.o  *~

distclean:	clean
		+rm -f poisson-amg-c.ex poisson-its-c.ex poisson-pcg-c.ex \
		poisson-gmg-c.ex spe01-its-c.ex \
		poisson-amg-f.ex spe01-pcg-f.ex 

