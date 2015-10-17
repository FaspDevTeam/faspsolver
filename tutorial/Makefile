#######################################################################
# Fast Auxiliary Space Preconditioners (FASP)
########################################################################

CSRC=$(shell find ./main -name '*.c')
COBJ=$(CSRC:.c=.o)

FASP_DIR=..

CC=gcc
FC=gfortran
LINKER=gfortran

CFLAGS=-I./include -I$(FASP_DIR)/include/
CFLAGS+=-O3

FFLAGS=-I./include -I$(FASP_DIR)/include/
FFLAGS+=-O3

LFLAGS=$(LINK_OPT) -L$(FASP_DIR)/lib -lfasp
LFLAGS+=-O3

.PHONY: all objs clean distclean help

all: objs poisson-amg-c poisson-its-c poisson-pcg-c \
     poisson-gmg-c spe01-its-c

objs: $(COBJ)

poisson-amg-c: main/poisson-amg.o
	@$(LINKER) main/poisson-amg.o -o $@.ex $(LFLAGS)
	@echo 'Building executable file $@.ex' ... Done

poisson-its-c: main/poisson-its.o
	@$(LINKER) main/poisson-its.o -o $@.ex $(LFLAGS)
	@echo 'Building executable file $@.ex' ... Done

poisson-pcg-c: main/poisson-pcg.o
	@$(LINKER) main/poisson-pcg.o -o $@.ex $(LFLAGS)
	@echo 'Building executable file $@.ex' ... Done

poisson-gmg-c: main/poisson-gmg.o
	@$(LINKER) main/poisson-gmg.o -o $@.ex $(LFLAGS)
	@echo 'Building executable file $@.ex' ... Done

spe01-its-c: main/spe01-its.o
	@$(LINKER) main/spe01-its.o -o $@.ex $(LFLAGS)
	@echo 'Building executable file $@.ex' ... Done

%.o: %.c
	@$(CC) $(CFLAGS) -c -o $@ $^ 
	@echo 'Building C object file $@' ... Done

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $^ 
	@echo 'Building F object file $@' ... Done

clean:
	@find . -name '*.o' -exec rm {} \;

distclean:
	@find . -name '*.o' -exec rm {} \;
	@rm -f *.ex *~

.PHONY: clean distclean
