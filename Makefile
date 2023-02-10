.POSIX:
.SUFFIXES:

# CC = cc
# HDF5_CFLAGS = -fPIC
# HDF5_MAJOR = 1
# HDF5_MINOR = 12
# HDF5_PATCH = 1
# HDF5_VERSION = $(HDF5_MAJOR).$(HDF5_MINOR).$(HDF5_PATCH)
# HDF5_PREFIX = $(PWD)/third_party/hdf5

# UNAME = $(shell uname)
# ifeq ($(UNAME), Darwin)
#   # NPROC = $(shell sysctl -n hw.ncpu)
#   SHARED_EXT = dylib
# else
#   # NPROC = $(shell nproc --all)
#   SHARED_EXT = so
# endif

.PHONY: all
all: haskell

# PREFIX = $(PWD)
# PACKAGE = lattice-symmetries-haskell
# GIT_COMMIT = $(shell git rev-parse --short HEAD)
# DIST = $(PACKAGE)-$(GIT_COMMIT)

HALIDE_PATH ?= $(PWD)/third_party/Halide
BIN_DIR = $(PWD)/kernels/build

TRUE_HALIDE_PATH := $(realpath $(HALIDE_PATH))
$(info $(TRUE_HALIDE_PATH))

UNAME = $(shell uname)
ifeq ($(UNAME), Darwin)
  SHARED_EXT = dylib
else
  SHARED_EXT = so
endif

.PHONY: haskell
haskell: cabal.project.local kernels
	cabal build

.PHONY: kernels
kernels: $(BIN_DIR)/libkernels.a

$(BIN_DIR)/libkernels.a:
	$(MAKE) -C kernels BIN_DIR=$(BIN_DIR) HALIDE_PATH=$(TRUE_HALIDE_PATH)

cabal.project.local:
	$(MAKE) -C kernels BIN_DIR=$(BIN_DIR) HALIDE_PATH=$(TRUE_HALIDE_PATH) ../cabal.project.local

# Determine the GHC version with which the library was built
GHC_VERSION := $(shell ghc --version | sed -e 's/[^0-9]*//')
# Find the shared library
LIBRARY_NAME = liblattice_symmetries_haskell.$(SHARED_EXT)
HASKELL_LIBRARY := $(shell find dist-newstyle/ -type f -name $(LIBRARY_NAME) | grep $(GHC_VERSION))

.PHONY: bindist
bindist: haskell
	export HASKELL_LIBRARY=$$(find dist-newstyle/ -type f -name $(LIBRARY_NAME) | grep $(GHC_VERSION)) && \
	mkdir -p bundle/include && \
	install -m644 cbits/lattice_symmetries_haskell.h bundle/include/ && \
	mkdir -p bundle/lib/haskell && \
	ldd $$HASKELL_LIBRARY | \
		grep $(SHARED_EXT) | \
		sed -e '/^[\^t]/d' | \
		sed -e 's/\t//' | \
		sed -e 's/ (0.*)//' | \
		grep -E '(libHS|libffi)' | \
		sed -e 's/.* => //' | \
		xargs -I '{}' install -m 644 '{}' bundle/lib/haskell/ && \
	install -m644 $$HASKELL_LIBRARY bundle/lib/ && \
	patchelf --set-rpath '$$ORIGIN/haskell' bundle/lib/$(LIBRARY_NAME) && \
	find bundle/lib/haskell -type f -exec patchelf --set-rpath '$$ORIGIN' {} \;

# .PHONY: release
# release: haskell
# 	mkdir -p $(DIST)/include
# 	mkdir -p $(DIST)/lib
# 	install -m644 -C cbits/lattice_symmetries_haskell.h $(DIST)/include/
# 	# install -m644 -C kernels/build/liblattice_symmetries_core.$(SHARED_EXT) $(DIST)/lib/
# 	find dist-newstyle -name "liblattice_symmetries_haskell.$(SHARED_EXT)" \
# 	  -exec install -m644 -C {} $(DIST)/lib/ \;
# ifeq ($(UNAME), Linux)
# 	LIBFFI=`ldd $(DIST)/lib/liblattice_symmetries_haskell.$(SHARED_EXT) | grep libffi | sed -E 's:.*=>\s+(.*/libffi.$(SHARED_EXT).[6-8]).*:\1:'`; cp -d $$LIBFFI* $(DIST)/lib/
# 	patchelf --set-rpath '$$ORIGIN' $(DIST)/lib/liblattice_symmetries_haskell.$(SHARED_EXT)
# endif
# 	tar -cf $(DIST).tar $(DIST)
# 	bzip2 $(DIST).tar
# ifneq ($(realpath $(PREFIX)), $(PWD))
# 	install -m644 -C $(DIST).tar.bz2 $(PREFIX)
# endif
# 	rm -r $(DIST)


# .PHONY: haskell
# haskell: kernels
# 	cabal v2-build

# .PHONY: kernels
# kernels: cabal.project.local
# kernels/build/liblattice_symmetries_core.$(SHARED_EXT)

# LS_HS_PATH = $(dir $(shell find dist-newstyle -name liblattice_symmetries_haskell.$(SHARED_EXT)))

# test_01: test/test_01.c cbits/lattice_symmetries_haskell.h
# 	$(CC) -I cbits -o $@ $< -L $(LS_HS_PATH) -Wl,-rpath=$(LS_HS_PATH) -llattice_symmetries_haskell

# kernels/build/liblattice_symmetries_core.$(SHARED_EXT): kernels/generator.cpp kernels/*.c cbits/lattice_symmetries_haskell.h
# 	cd kernels && $(MAKE) Halide
# 	cd kernels && $(MAKE)

# .PHONY: hdf5
# hdf5: $(HDF5_PREFIX) $(HDF5_PREFIX)/lib/pkgconfig/hdf5.pc

# third_party/hdf5-$(HDF5_VERSION).tar.bz2:
# 	@mkdir -p $(@D)
# 	cd third_party && \
# 	wget -q https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-$(HDF5_MAJOR).$(HDF5_MINOR)/hdf5-$(HDF5_VERSION)/src/hdf5-$(HDF5_VERSION).tar.bz2

# third_party/hdf5-$(HDF5_VERSION)-src: third_party/hdf5-$(HDF5_VERSION).tar.bz2
# 	cd third_party && \
# 	tar -xf hdf5-$(HDF5_VERSION).tar.bz2 && \
# 	mv hdf5-$(HDF5_VERSION) hdf5-$(HDF5_VERSION)-src

# $(HDF5_PREFIX): third_party/hdf5-1.12.1-src
# 	cd third_party/hdf5-1.12.1-src && \
# 	CC=$(CC) CFLAGS="$(HDF5_CFLAGS)" ./configure --prefix="$(HDF5_PREFIX)" \
# 	    --disable-java --disable-fortran --disable-cxx \
# 	    --disable-tests --disable-tools \
# 	    --disable-shared && \
# 	make -j $(NPROC) install

# $(HDF5_PREFIX)/lib/pkgconfig/hdf5.pc:
# 	@mkdir -p $(@D)
# 	cd third_party/hdf5/lib/pkgconfig && \
# 	echo "# Package Information for pkg-config" >>hdf5.pc && \
# 	echo "" >>hdf5.pc && \
# 	echo "prefix=$(HDF5_PREFIX)" >>hdf5.pc && \
# 	echo "libdir=\$${prefix}/lib" >>hdf5.pc && \
# 	echo "includedir=\$${prefix}/include" >>hdf5.pc && \
# 	echo "" >>hdf5.pc && \
# 	echo "Name: HDF5" >>hdf5.pc && \
# 	echo "Description: HDF5" >>hdf5.pc && \
# 	echo "Version: $(HDF5_MAJOR).$(HDF5_MINOR).$(HDF5_PATCH)" >>hdf5.pc && \
# 	echo "Libs: -L\$${libdir} -lhdf5_hl -lhdf5 -lz" >>hdf5.pc && \
# 	echo "Cflags: -I\$${includedir}" >>hdf5.pc

.PHONY: clean
clean:
	@$(MAKE) -C kernels BIN_DIR=$(BIN_DIR) HALIDE_PATH=$(TRUE_HALIDE_PATH) clean
