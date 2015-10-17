#######################################################################
# Fast Auxiliary Space Preconditioners (FASP)
########################################################################

DSRC=$(shell find ./src -name '*.c')
DOBJ=$(DSRC:.c=.o)

PSRC=$(shell find ./main -name '*.c')
POBJ=$(PSRC:.c=.o)

FASP_DIR=..

LINKER=gfortran

CFLAGS=-I./include -I$(FASP_DIR)/include/
CFLAGS+=-O3

LFLAGS=$(LINK_OPT) -L$(FASP_DIR)/lib -lfasp
LFLAGS+=-O3

.PHONY: all objs clean distclean help

all: objs test testbamg testbcsr testbsr testgmg \
     testrap testlab testmat testrap testsmat    \
     regression regression_mf regression_mm      \
     testfem testheat

regression: main/regression.o
	@$(LINKER) main/regression.o -o $@.ex $(LFLAGS) 
	@echo 'Building executable file $@.ex' ... Done

regression_mf: main/regression_mf.o
	@$(LINKER) main/regression_mf.o -o $@.ex $(LFLAGS) 
	@echo 'Building executable file $@.ex' ... Done

regression_mm: main/regression_mm.o
	@$(LINKER) main/regression_mm.o -o $@.ex $(LFLAGS) 
	@echo 'Building executable file $@.ex' ... Done

test: main/test.o
	@$(LINKER) main/test.o -o $@.ex $(LFLAGS) 
	@echo 'Building executable file $@.ex' ... Done

testbamg: main/testbamg.o
	@$(LINKER) main/testbamg.o -o $@.ex $(LFLAGS) 
	@echo 'Building executable file $@.ex' ... Done

testbcsr: main/testbcsr.o
	@$(LINKER) main/testbcsr.o -o $@.ex $(LFLAGS) 
	@echo 'Building executable file $@.ex' ... Done

testbsr: main/testbsr.o
	@$(LINKER) main/testbsr.o -o $@.ex $(LFLAGS) 
	@echo 'Building executable file $@.ex' ... Done

testgmg: main/testgmg.o
	@$(LINKER) main/testgmg.o -o $@.ex $(LFLAGS) 
	@echo 'Building executable file $@.ex' ... Done

testlab: main/testlab.o
	@$(LINKER) main/testlab.o -o $@.ex $(LFLAGS) 
	@echo 'Building executable file $@.ex' ... Done

testmat: main/testmat.o
	@$(LINKER) main/testmat.o -o $@.ex $(LFLAGS) 
	@echo 'Building executable file $@.ex' ... Done

testsmat: main/testsmat.o
	@$(LINKER) main/testsmat.o -o $@.ex $(LFLAGS) 
	@echo 'Building executable file $@.ex' ... Done

testrap: main/testrap.o $(DOBJ)
	@$(LINKER) main/testrap.o $(DOBJ) -o $@.ex $(LFLAGS) 
	@echo 'Building executable file $@.ex' ... Done

testfem: main/testfem.o $(DOBJ)
	@$(LINKER) main/testfem.o $(DOBJ) -o $@.ex $(LFLAGS) 
	@echo 'Building executable file $@.ex' ... Done

testheat: main/testheat.o $(DOBJ)
	@$(LINKER) main/testheat.o $(DOBJ) -o $@.ex $(LFLAGS) 
	@echo 'Building executable file $@.ex' ... Done

objs: $(DOBJ) $(POBJ)

%.o: %.c
	@$(CC) $(CFLAGS) -c -o $@ $^ 
	@echo 'Building C object file $@' ... Done

clean:
	@find . -name '*.o' -exec rm {} \;

distclean:
	@find . -name '*.o' -exec rm {} \;
	@rm -f *.ex *~

.PHONY: clean distclean
