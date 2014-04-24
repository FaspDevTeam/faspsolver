########################################################################
# Fast Auxiliary Space Preconditioners (FASP) 
#
# Top level Makefile: Calls cmake to configure and build the library
# and the test suite.
#
########################################################################

####################   User Defined Options   ##########################
#
# The default setting for build type for FASP is RELEASE. The compiler 
# options then include by default "-Wall -g" as well as "-DDEBUG_MODE". 
# The RELEASE build type by default has the "-O3". If you want to work 
# with build type DEBUG, then uncomment the next line:
#
# debug=yes
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
# These user options can also be applied as make command line options.
# For example, to enforce the debug compiling options:
#
# make config debug=yes
#
#-------------------------------------------------------------------------
#
# By default, FASP uses the command-line Doxygen to generate a reference
# manual. If you want to use the GUI of Doxgen instead of command-line
# (if there is one installed on your system), uncomment the next line:
#
# doxywizard=yes
#
#-------------------------------------------------------------------------
# If you want to use UMFPACK (part of SparseSuite), uncomment the next 
# line:
# 
# umfpack=yes
#
# If you want to specify the path to SparseSuite, uncomment the next 
# line and give the correct path to SparseSuite here. For example:
#
# suitesparse_dir="/dir/to/SuiteSparse"
#-------------------------------------------------------------------------
# If you want to use SuperLU, uncomment the next line:
#
# superlu=yes
#
# If you want to specify the path to SuperLU, uncomment the next line
# and give the correct path to SuperLU here. For example:
#
# superlu_dir="/dir/to/SuperLU"
#-------------------------------------------------------------------------
# If you want to use MUMPS, uncomment the next line:
#
# mumps=yes
#
# If you want to specify the path to MUMPS, uncomment the next line
# and give the correct path to MUMPS here. For example:
#
# mumps_dir="/dir/to/MUMPS"
#
####################  User Defined Compiler Flags  #####################
ifeq ($(debug),yes)
	cflags="-Wall -g -DDEBUG_MODE"
	cxxflags="-Wall -g -DDEBUG_MODE"
	fflags="-Wall -g -DDEBUG_MODE"
else
	cflags="-O3 -funroll-loops"
	cxxflags="-O3 -funroll-loops"
	fflags="-O3 -funroll-loops"
endif
####################  User Changes UP TO HERE   ########################

# Let cmake do the configuration. Set up a build dir
#cpu0=$(shell uname -m | sed -e 's/[[:space:]][[:space:]]*/_/g')
#sys0=$(shell uname -s)
#build_dir=BUILD_$(cpu0)-$(sys0)
build_dir=BUILD_FASP

CONFIG_FLAGS=-DCMAKE_RULE_MESSAGES=ON

ifeq ($(verbose),yes)
    CONFIG_FLAGS+=-DCMAKE_VERBOSE_MAKEFILE=ON
else
    CONFIG_FLAGS+=-DCMAKE_VERBOSE_MAKEFILE=OFF
endif

ifeq ($(debug),yes)
    CONFIG_FLAGS+=-DCMAKE_BUILD_TYPE=DEBUG
else
    CONFIG_FLAGS+=-DCMAKE_BUILD_TYPE=RELEASE
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

ifeq ($(superlu), yes)
    CONFIG_FLAGS+=-DUSE_SUPERLU=$(superlu)
    CONFIG_FLAGS+=-DSUPERLU_DIR=$(superlu_dir)
endif

ifeq ($(mumps), yes)
    CONFIG_FLAGS+=-DUSE_MUMPS=$(mumps)
    CONFIG_FLAGS+=-DMUMPS_DIR=$(mumps_dir)
endif

CONFIG_FLAGS+=-DADD_CFLAGS=$(cflags) 
CONFIG_FLAGS+=-DADD_CXXFLAGS=$(cxxflags) 
CONFIG_FLAGS+=-DADD_FFLAGS=$(fflags) 

all clean install docs headers:
	@if [ ! -f $(build_dir)/Makefile ] ; then \
		echo "Configuration not found! Please perform configuration first."; \
		echo "See the following help screen for usage ..."; \
		echo " "; \
		cat INSTALL; \
	else \
	  	make -C $(build_dir) $@ ; \
	fi

config: distclean
	mkdir -p $(build_dir)
	cd $(build_dir) && cmake $(CURDIR) $(CONFIG_FLAGS)

uninstall:
	@if [ ! -f $(build_dir)/install_manifest.txt ]; then \
		echo "Installation manifest not found! Nothing to uninstall."; \
		echo "See the following help screen for usage ..."; \
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
	@-zip -r faspsolver.zip README INSTALL License Makefile VERSION     \
	                        base data test tutorial *.txt *.cmake *.tcl \
                            doc/userguide.pdf doc/refman.pdf vs08 vs10

version:
	@-hg -q id > VERSION
	@-hg log -r "." --template 'FASP {latesttag}.{latesttagdistance}:'
	@-cat VERSION

.PHONY: all backup config clean distclean install uninstall docs headers help version
