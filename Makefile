########################################################################
# Fast Auxiliary Space Preconditioners (FASP) 
#
# Top level Makefile: Calls cmake to configure and build the library
# and the test suite. modeled on METIS's Makefile.
#
########################################################################

openmp = no

# Let cmake do the configuration. Set up a build dir

cpu0=$(shell uname -m | sed -e 's/[[:space:]][[:space:]]*/_/g')
sys0=$(shell uname -s)
build_dir = build/$(cpu0)-$(sys0)

CONFIG_FLAGS = -DCMAKE_VERBOSE_MAKEFILE=1
 
ifneq ($(openmp), no)
    CONFIG_FLAGS += -DOPENMP=$(openmp)
endif

VER0=1.2.1
PKG0=fasp-$(VER0)

define do-cmake-build
mkdir -p $(build_dir)
cd $(build_dir) && cmake $(CURDIR) $(CONFIG_FLAGS)
endef

all clean install docs headers:
	@if [ ! -f $(build_dir)/Makefile ]; then \
		echo "Configuration not found! Please erform configuration first."; \
		echo "See the following help screen for usages ..."; \
		echo " "; \
		cat INSTALL; \
	else \
	  	make -C $(build_dir) $@ $(MAKEFLAGS); \
	fi

config: distclean
	$(do-cmake-build)

uninstall:
	@if [ ! -f $(build_dir)/install_manifest.txt ]; then \
		echo "Installation manifest not found! Nothing to uninstall."; \
		echo "See the following help screen for usages ..."; \
		echo " "; \
		cat INSTALL; \
	else \
		xargs rm < $(build_dir)/install_manifest.txt; \
	fi

distclean:
	rm -rf $(build_dir)

help:
	@cat INSTALL

backup:
	@-rm -f faspsolver.zip
	@-zip -r faspsolver.zip README INSTALL License Makefile \
		                    *.txt *.cmake doc/userguide.pdf \
		                    base data test tutorial vs08

.PHONY: config distclean all clean install docs headers uninstall backup help
