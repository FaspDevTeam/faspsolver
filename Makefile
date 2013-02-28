########################################################################
# Fast Auxiliary Space Preconditioners (FASP) 
#
# Top level Makefile: Calls cmake to configure and build the library
# and the test suite.
#
########################################################################

####################   User Defined Options   ##########################
#
# The default setting for vebosity level for FASP is verbose=no. If you
# want to increase verbosity level, uncomment the next line:
#
# verbose=yes
#
# By default, FASP generates static libraries. If you need to generate 
# shared libs instead of static libs, uncomment the next line:
#
# shared=yes
#
# If you want to compile with OpenMP support, uncomment the next line:
#
# openmp=yes
#
# By default, FASP uses the command-line Doxygen to generate a reference
# manual. If you want to use the GUI of Doxgen instead of command-line
# (if there is one installed on your system), uncomment the next line:
#
# doxywizard=yes
#
# If you want to use the UMFPACK, uncomment the next line: 
# 
# umfpack=yes
#
# If you want to specify the path to SparseSuite, uncomment the next 
# line and give the correct path to it:
#
# suitesparse_dir="/Users/XiaozheHu/Research/Package/SuiteSparse"
#
####################  User Defined Compiler Flags  #####################
#
cflags="-O3 -funroll-all-loops"
cxxflags="-O3 -funroll-all-loops"
fflags="-O3 -funroll-all-loops"
#
# If you need to generate debug information during building, you should 
# uncomment the next few lines: 
#
# cflags="-Wall -g -pg"
# cxxflags="-Wall -g -pg"
# fflags="-Wall -g -pg"
#
####################  User Changes UP TO HERE   ########################

VER0=1.2.1
PKG0=fasp-$(VER0)

# Let cmake do the configuration. Set up a build dir
cpu0=$(shell uname -m | sed -e 's/[[:space:]][[:space:]]*/_/g')
sys0=$(shell uname -s)
build_dir=BUILD_$(cpu0)-$(sys0)

CONFIG_FLAGS=-DCMAKE_VERBOSE_MAKEFILE=OFF -DCMAKE_RULE_MESSAGES=ON
ifeq ($(verbose),yes)
    CONFIG_FLAGS=-DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_RULE_MESSAGES=ON
endif
ifeq ($(shared),yes)
    CONFIG_FLAGS+=-DSHARED=$(shared)
endif
ifeq ($(openmp),yes)
    CONFIG_FLAGS+=-DUSE_OPENMP=$(openmp)
endif
ifeq ($(doxywizard),yes)
    CONFIG_FLAGS+=-DDOXYWIZARD=$(doxywizard)
endif
ifeq ($(umfpack), yes)
    CONFIG_FLAGS+=-DUSE_UMFPACK=$(umfpack)
    CONFIG_FLAGS+=-DSUITESPARSE_DIR=$(suitesparse_dir)
endif
CONFIG_FLAGS+=-DCMAKE_C_FLAGS=$(cflags)
CONFIG_FLAGS+=-DCMAKE_CXX_FLAGS=$(cxxflags)
CONFIG_FLAGS+=-DCMAKE_Fortran_FLAGS=$(fflags)

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
		       doc/htdocs; \
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
