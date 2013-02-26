########################################################################
# Fast Auxiliary Space Preconditioners (FASP) 
#
# Top level Makefile: Calls cmake to configure and build the library
# and the test suite. modeled on METIS's Makefile.
#
########################################################################

openmp     = no	    # openmp support
shared     = no     # shared or static libs
doxywizard = no     # use the GUI for doxygen 
		    # (if there is one in a standard location)
verbose = no	    # verbose compiling


############## COMPLER FLAGS 
cflags="-O3 -funroll-all-loops"
cxxflags="-O3 -funroll-all-loops"
fflags="-O3 -funroll-all-loops"

########################CHANGE UP TO HERE ##############

# Let cmake do the configuration. Set up a build dir
cpu0=$(shell uname -m | sed -e 's/[[:space:]][[:space:]]*/_/g')
sys0=$(shell uname -s)
build_dir = BUILD_$(cpu0)-$(sys0)

CONFIG_FLAGS = -DCMAKE_VERBOSE_MAKEFILE=0
ifneq ($(verbose), no)
    CONFIG_FLAGS = -DCMAKE_VERBOSE_MAKEFILE=1
endif
ifneq ($(openmp), no)
    CONFIG_FLAGS += -DOPENMP=$(OPENMP)
endif
ifneq ($(shared), no)
    CONFIG_FLAGS += -DSHARED=$(shared)
endif
ifneq ($(doxywizard), no)
    CONFIG_FLAGS += -DDOXYWIZARD=$(doxywizard)
endif
    CONFIG_FLAGS += -DCMAKE_C_FLAGS=$(cflags) -DCMAKE_CXX_FLAGS=$(cflags) -DCMAKE_Fortran_FLAGS=$(fflags)


VER0=1.2.1
PKG0=fasp-$(VER0)

all clean install docs headers:
	@if [ ! -f $(build_dir)/Makefile ]; then \
		echo "Configuration not found! Please perform configuration first."; \
		echo "See the following help screen for usage ..."; \
		echo " "; \
		cat INSTALL; \
	else \
	  	make -C $(build_dir) $@ $(MAKEFLAGS); \
	fi

config: distclean
	mkdir -p $(build_dir)
	cd $(build_dir) && cmake $(CURDIR) $(CONFIG_FLAGS)


uninstall:
	@if [ ! -f $(build_dir)/install_manifest.txt ]; then \
		echo "Installation manifest not found! Nothing to uninstall."; \
		echo "See the following help screen for usages ..."; \
		echo " "; \
		cat INSTALL; \
	else \
		xargs rm < $(build_dir)/install_manifest.txt; \
		rm -rf $(build_dir)/install_manifest.txt \
			doc/htdocs \
			doc/doc.zip; \
	fi

distclean:
	@-rm -rf $(build_dir)
	@-find . -name '*~' -exec rm {} \;

help:
	@cat INSTALL

backup:
	@-rm -f faspsolver.zip
	@-zip -r faspsolver.zip README INSTALL License Makefile \
		                    *.txt *.cmake doc/userguide.pdf \
		                    base data test tutorial vs08

.PHONY: config distclean all clean install docs headers uninstall backup help
