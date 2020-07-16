#######################################################################
# Fast Auxiliary Space Preconditioners (FASP) 
#
########################################################################
#
# TOP LEVEL FASP Makefile: Calls cmake to configure and build
# the library and the test suite. 
# 
#   Hopefully you will *NOT NEED TO CHANGE* this top level Makefile.
#
#   FOR USER DEFINED OPTIONS, GO IN "FASP.mk" which is included by this
#   Makefile. Copy "FASP.mk.example file to "FASP.mk" and edit it to
#   adjust the settings to your liking, and then type "make help" to
#   see how to configure/build FASP.
# 
#  Modified   2015-10-18   --ltz
#  Modified   2017-01-10   --zcs
#  Modified   2019-08-08   --zcs
#  Modified   2020-05-06   --zcs
########################################################################
prefix = no-prefix

ifneq ($(wildcard FASP.mk),)
	include ./FASP.mk
endif

ifeq ($(debug),yes)
	cflags="-Wall -g"
	cxxflags="-Wall -g"
	fflags="-Wall -g"
endif
#
ifeq ($(debug),some)
	cflags="-Wall -g -DDEBUG_MODE=1"
	cxxflags="-Wall -g -DDEBUG_MODE=1"
	fflags="-Wall -g -DDEBUG_MODE=1"
endif
#
ifeq ($(debug),all)
	cflags="-Wall -g -DDEBUG_MODE=3 -DCHMEM_MODE=1"
	cxxflags="-Wall -g -DDEBUG_MODE=3 -DCHMEM_MODE=1"
	fflags="-Wall -g -DDEBUG_MODE=3 -DCHMEM_MODE=1"
endif
####################  User Changes UP TO HERE   ########################

# Let cmake do the configuration. Set up a build dir
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

ifneq ($(prefix),no-prefix)
    CONFIG_FLAGS+=-DFASP_INSTALL_PREFIX=$(prefix)
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
    CONFIG_FLAGS+=-DMETIS_DIR=$(metis_dir)
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

ifeq ($(pardiso), yes)
    CONFIG_FLAGS+=-DUSE_PARDISO=$(pardiso)
    CONFIG_FLAGS+=-DMKL_DIR=$(mkl_dir)
endif

CONFIG_FLAGS+=-DADD_CFLAGS=$(cflags)
CONFIG_FLAGS+=-DADD_CXXFLAGS=$(cxxflags)
CONFIG_FLAGS+=-DADD_FFLAGS=$(fflags)

fasp headers docs clean: 
	@if [ ! -f $(build_dir)/Makefile ] ; then \
		echo "*=======================================================================*"; \
		echo "* WARNING: Configuration not found! Please perform 'make config' first. *"; \
		echo "* See the following help screen for usage ...                           *"; \
		echo "*=======================================================================*"; \
		echo " "; \
		cat INSTALL; \
	else \
		make -C $(build_dir) $@ ; \
	fi

install: 
	@if [ ! -f $(build_dir)/Makefile ] ; then \
		echo "*=======================================================================*"; \
		echo "* WARNING: Configuration not found! Please perform 'make config' first. *"; \
		echo "* See the following help screen for usage ...                           *"; \
		echo "*=======================================================================*"; \
		echo " "; \
		cat INSTALL; \
	else \
		umaskz=`umask`;  umask 0022; \
		make -C $(build_dir) $@ ; umask $$umaskz ; \
	fi

test tutorial:
	@if [ ! -f $(build_dir)/$@/Makefile ] ; then \
		echo "*=======================================================================*"; \
		echo "* WARNING: Configuration not found! Please perform 'make config' first. *"; \
		echo "* See the following help screen for usage ...                           *"; \
		echo "*=======================================================================*"; \
		echo " "; \
		cat INSTALL; \
	else \
		make -C $(build_dir)/$@ install ; \
	fi

config:	distclean
	@if [ ! -f ./FASP.mk ] ; then \
		echo "*=======================================================================*"; \
		echo "* WARNING: FASP.mk is missing from the current directory!               *"; \
		echo "* Using the DEFAULT configuration instead ...                           *"; \
		echo "*=======================================================================*"; \
	fi
	@mkdir -p $(build_dir) ; 
	@cd $(build_dir) && cmake $(CURDIR) $(CONFIG_FLAGS) 

uninstall:
	@if [ ! -f $(build_dir)/install_manifest.txt ]; then \
		echo "*=======================================================================*"; \
		echo "* WARNING: Installation manifest not found! Nothing to be uninstalled.  *"; \
		echo "* Type \"make help\" for help on usage ...                                *"; \
		echo "*=======================================================================*"; \
		echo " "; \
	else \
		xargs rm < $(build_dir)/install_manifest.txt; \
		rm -rf $(build_dir)/install_manifest.txt \
		       doc/htdocs; \
		echo "FASP library and header files have been successfully uninstalled."; \
	fi

distclean:
	@-rm -rf Config.mk
	@-rm -rf $(build_dir)   
	@-find . -name '*~' -exec rm {} \;

help:
	@-clear
	@-cat INSTALL

backup:
	@-rm -f faspsolver.zip
	@-zip -r faspsolver.zip README INSTALL License Makefile VERSION FASP* \
	                        base data test tutorial modules log util vs10 \
                                vs15 vs19 *.txt *.tcl doc/*.pdf doc/*.in doc/FAQ 

version:
	@-git describe --tags > VERSION
	@-cat VERSION

.PHONY:	backup config fasp install tutorial test clean distclean uninstall docs headers help version 
