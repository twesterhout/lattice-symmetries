.POSIX:
.SUFFIXES:

# NOTE: You probably want to override the following
export HALIDE_PATH = third_party/Halide
export BIN_DIR = kernels-build

# Command to use to run docker
SUDO = sudo

.PHONY: all
all: haskell

# If HALIDE_PATH is not an absolute path, make it so, because otherwise our
# sub-makefile will fail
ifneq ($(HALIDE_PATH), $(realpath $(HALIDE_PATH)))
  TRUE_HALIDE_PATH := $(realpath $(HALIDE_PATH))
  export HALIDE_PATH = $(TRUE_HALIDE_PATH)
  $(info "Setting HALIDE_PATH to $(HALIDE_PATH) ...")
endif

ifneq ($(BIN_DIR), $(abspath $(BIN_DIR)))
  TRUE_BIN_DIR := $(abspath $(BIN_DIR))
  export BIN_DIR = $(TRUE_BIN_DIR)
  $(info "Setting BIN_DIR to $(HALIDE_PATH) ...")
endif

UNAME = $(shell uname)
ifeq ($(UNAME), Darwin)
  SHARED_EXT = dylib
else
  SHARED_EXT = so
endif

.PHONY: haskell
haskell: cabal.project.local kernels
	cabal build

.PHONY: check
check: haskell
	cabal test --test-show-details=direct

.PHONY: kernels
kernels: $(BIN_DIR)/libkernels.a

$(BIN_DIR)/libkernels.a:
	$(MAKE) -C kernels
	# $(MAKE) -C kernels BIN_DIR=$(BIN_DIR) HALIDE_PATH=$(TRUE_HALIDE_PATH)

cabal.project.local:
	$(MAKE) -C kernels ../cabal.project.local
	# $(MAKE) -C kernels BIN_DIR=$(BIN_DIR) HALIDE_PATH=$(TRUE_HALIDE_PATH) ../cabal.project.local

# Determine the GHC version with which the library was built
GHC_VERSION := $(shell ghc --version | sed -e 's/[^0-9]*//')
# Find the shared library
LIBRARY_NAME = liblattice_symmetries_haskell.$(SHARED_EXT)
HASKELL_LIBRARY := $(shell find dist-newstyle/ -type f -name $(LIBRARY_NAME) | grep $(GHC_VERSION))

.PHONY: bundle
bundle: haskell
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

.PHONY: bundle-docker
bundle-docker: 
	mkdir -p $@
	mkdir -p kernel-build-docker
	mkdir -p dist-newstyle-docker
	WORKDIR=/work/lattice-symmetries-haskell && \
	$(SUDO) docker run --rm \
	  -v $$PWD/src:$$WORKDIR/src:ro \
	  -v $$PWD/lib:$$WORKDIR/lib:ro \
	  -v $$PWD/test:$$WORKDIR/test:ro \
	  -v $$PWD/cbits:$$WORKDIR/cbits:ro \
	  -v $$PWD/kernels:$$WORKDIR/kernels:ro \
	  -v $$PWD/LICENSE:$$WORKDIR/LICENSE:ro \
	  -v $$PWD/README.md:$$WORKDIR/README.md:ro \
	  -v $$PWD/cabal.project:$$WORKDIR/cabal.project:ro \
	  -v $$PWD/lattice-symmetries-haskell.cabal:$$WORKDIR/lattice-symmetries-haskell.cabal:ro \
	  -v $$PWD/Makefile:$$WORKDIR/Makefile:ro \
	  -v $$PWD/$@:$$WORKDIR/bundle:z \
	  -v $$PWD/kernel-build-docker:$$WORKDIR/kernel-build:z \
	  -v $$PWD/dist-newstyle-docker:$$WORKDIR/dist-newstyle:z \
	  twesterhout/lattice-symmetries-haskell \
	  bash -c 'make HALIDE_PATH=/opt/Halide bundle'
	$(SUDO) chown -R $$USER:$$USER $@
	$(SUDO) chown -R $$USER:$$USER kernel-build-docker
	$(SUDO) chown -R $$USER:$$USER dist-newstyle-docker

.PHONY: conda-package
conda-package:
	@rm -rf python/build python/*.egg-info python/lattice_symmetries/__pycache__ python/lattice_symmetries/*.so
	conda build python/conda

.PHONY: clean
clean:
	@$(MAKE) -C kernels clean
	cabal clean
	rm -rf dist-newstyle-docker kernel-build-docker bundle-docker

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

