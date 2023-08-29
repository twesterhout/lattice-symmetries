.RULES:
.POSIX:

OPTIMIZATION ?= --debug
COMPILER ?= llvm
PREFIX = /usr/local
# LS_HS_PATH = $(PWD)/lattice-symmetries-haskell
# PRIMME_PATH = /home/tom/src/primme

SUDO = sudo
UNAME = $(shell uname)
CFLAGS = $(OPTIMIZATION) --target-compiler=$(COMPILER)
# TODO: enable later to get rid of __primitive calls: -s cPtrToLogicalValue=true
LDFLAGS = -llattice_symmetries_haskell
# -llattice_symmetries_core
# ifeq ($(UNAME), Linux)
#   LDFLAGS += -lnuma
# endif

# LS_HS_CFLAGS = -I$(LS_HS_PATH)/include
# LS_HS_LDFLAGS = -L$(LS_HS_PATH)/lib -llattice_symmetries_haskell

PRIMME_CFLAGS = -I$(PRIMME_PATH)/include
PRIMME_LDFLAGS = -L$(PRIMME_PATH)/lib -lprimme -llapacke -lopenblas -lm -lgomp -lpthread

ifeq ($(UNAME), Linux)
  HDF5_CFLAGS ?= $(shell pkg-config --cflags hdf5)
  HDF5_LDFLAGS ?= -lhdf5_hl $(shell pkg-config --libs hdf5)
  # CONDA_CC ?= $(shell conda run -n ci_devel bash -c "which \$${CC}")
  # CONDA_PREFIX ?= $(shell conda run -n ci_devel bash -c "echo \$${CONDA_PREFIX}")
  SHARED_EXT = so
  # SHARED_FLAG = -shared -rdynamic
else
  HDF5_CFLAGS ?=
  HDF5_LDFLAGS ?= -lhdf5_hl -lhdf5
  # CONDA_CC = $(CC)
  # CONDA_PREFIX = 
  SHARED_EXT = dylib
  # SHARED_FLAG = -dynamiclib
endif

# ifeq ($(UNAME), Linux)
#   CHPL_LIBS = LD_LIBRARY_PATH=$(PWD)/third_party/lib:$$LD_LIBRARY_PATH
# endif
# ifeq ($(UNAME), Darwin)
#   CHPL_LIBS = DYLD_LIBRARY_PATH=$(PWD)/third_party/lib:$$DYLD_LIBRARY_PATH
# endif
# CHPL_ARGS =


PREFIX = $(PWD)
PACKAGE = lattice-symmetries-chapel
GIT_COMMIT = $(shell git rev-parse --short HEAD)
DIST = $(PACKAGE)-$(GIT_COMMIT)

# MODULES = src/ApplyOperator.chpl src/StatesEnumeration.chpl src/helper.c
LIB_MODULES = src/LatticeSymmetries.chpl \
              src/FFI.chpl \
              src/ForeignTypes.chpl \
              src/Vector.chpl \
              src/ConcurrentAccessor.chpl \
              src/StatesEnumeration.chpl \
              src/BatchedOperator.chpl \
              src/DistributedMatrixVector.chpl \
              src/CommonParameters.chpl \
              src/Utils.chpl

APP_MODULES = $(LIB_MODULES) \
              src/HashedToBlock.chpl \
              src/BlockToHashed.chpl \
              src/MyHDF5.chpl

.PHONY: all
all: lib
# examples

.PHONY: examples
examples: bin/Example02 bin/Example05

# .PHONY: test
# test: bin/TestStatesEnumeration bin/TestMatrixVectorProduct

.PHONY: check
check: check-states-enumeration check-matrix-vector-product

.PHONY: benchmark-states-enumeration
benchmark-states-enumeration: bin/TestStatesEnumeration
	# $(CHPL_LIBS) $< $(CHPL_ARGS) --kBasis data/heisenberg_pyrochlore_2x2x2.yaml --kRepresentatives data/large-scale/construction/heisenberg_pyrochlore_2x2x2.h5
	# $(CHPL_LIBS) $< $(CHPL_ARGS) --kBasis data/heisenberg_kagome_36.yaml --kRepresentatives data/large-scale/construction/heisenberg_kagome_36.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_square_6x6.yaml --kRepresentatives data/large-scale/matvec/heisenberg_square_6x6.h5

.PHONY: check-states-enumeration
check-states-enumeration: bin/TestStatesEnumeration
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_4.yaml --kRepresentatives data/matvec/heisenberg_chain_4.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_6.yaml --kRepresentatives data/matvec/heisenberg_chain_6.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_8.yaml --kRepresentatives data/matvec/heisenberg_chain_8.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_10.yaml --kRepresentatives data/matvec/heisenberg_chain_10.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_12.yaml --kRepresentatives data/matvec/heisenberg_chain_12.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_16.yaml --kRepresentatives data/matvec/heisenberg_chain_16.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_20.yaml --kRepresentatives data/matvec/heisenberg_chain_20.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_24.yaml --kRepresentatives data/matvec/heisenberg_chain_24.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_24_symm.yaml --kRepresentatives data/matvec/heisenberg_chain_24_symm.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_kagome_12.yaml --kRepresentatives data/matvec/heisenberg_kagome_12.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_kagome_12_symm.yaml --kRepresentatives data/matvec/heisenberg_kagome_12_symm.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_kagome_16.yaml --kRepresentatives data/matvec/heisenberg_kagome_16.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_square_4x4.yaml --kRepresentatives data/matvec/heisenberg_square_4x4.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_square_5x5.yaml --kRepresentatives data/matvec/heisenberg_square_5x5.h5

.PHONY: benchmark-matrix-vector-product
benchmark-matrix-vector-product: bin/TestMatrixVectorProduct
	# $(CHPL_LIBS) $< $(CHPL_ARGS) --kBasis data/heisenberg_pyrochlore_2x2x2.yaml --kRepresentatives data/heisenberg_pyrochlore_2x2x2.h5
	# $(CHPL_LIBS) $< $(CHPL_ARGS) --kBasis data/heisenberg_kagome_36.yaml --kRepresentatives data/heisenberg_kagome_36.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_square_6x6.yaml --kVectors data/large-scale/matvec/heisenberg_square_6x6.h5

.PHONY: check-matrix-vector-product
check-matrix-vector-product: bin/TestMatrixVectorProduct
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/issue_01.yaml --kVectors data/matvec/issue_01.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_4.yaml --kVectors data/matvec/heisenberg_chain_4.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_6.yaml --kVectors data/matvec/heisenberg_chain_6.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_8.yaml --kVectors data/matvec/heisenberg_chain_8.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_10.yaml --kVectors data/matvec/heisenberg_chain_10.h5
	# $(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_12.yaml --kVectors data/matvec/heisenberg_chain_12.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_16.yaml --kVectors data/matvec/heisenberg_chain_16.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_20.yaml --kVectors data/matvec/heisenberg_chain_20.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_24.yaml --kVectors data/matvec/heisenberg_chain_24.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_chain_24_symm.yaml --kVectors data/matvec/heisenberg_chain_24_symm.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_kagome_12.yaml --kVectors data/matvec/heisenberg_kagome_12.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_kagome_12_symm.yaml --kVectors data/matvec/heisenberg_kagome_12_symm.h5
	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_kagome_16.yaml --kVectors data/matvec/heisenberg_kagome_16.h5


TEST_DATA_URL = https://surfdrive.surf.nl/files/index.php/s/OK5527Awfgl1hT2/download

.PHONY: data/construction
data/construction:
	mkdir -p data && cd data && \
	wget -q -O tmp.zip $(TEST_DATA_URL)?path=%2Fdata%2Fconstruction && \
	unzip tmp.zip && rm tmp.zip

.PHONY: data/matvec
data/matvec:
	mkdir -p data && cd data && \
	wget -q -O tmp.zip $(TEST_DATA_URL)?path=%2Fdata%2Fmatvec && \
	unzip tmp.zip && rm tmp.zip

.PHONY: data/large-scale
data/large-scale:
	mkdir -p data && cd data && \
	wget -q -O tmp.zip $(TEST_DATA_URL)?path=%2Fdata%2Flarge-scale && \
	unzip tmp.zip && rm tmp.zip

lib: lib/liblattice_symmetries_chapel.$(SHARED_EXT)

src/library.o: src/library.c
	$(CC) -c $^ -fPIC -o $@

lib/liblattice_symmetries_chapel.$(SHARED_EXT): $(LIB_MODULES) src/library.o
	@mkdir -p $(@D)
	@rm -rf $(@D)/*
	CHPL_LIB_PIC=pic chpl --print-commands --devel $(CFLAGS) $(CHPL_CFLAGS) --library --dynamic -o lattice_symmetries_chapel $^ \
		$(CHPL_LDFLAGS) -llattice_symmetries_haskell

install: lib
	@mkdir -p $(PREFIX)/lib
	# @mkdir -p $(PREFIX)/include
	install -m 755 lib/liblattice_symmetries_chapel.$(SHARED_EXT) $(PREFIX)/lib/
	# install -m 644 lib/lattice_symmetries_chapel.h $(PREFIX)/include/

.PHONY: bundle
bundle: lib
	mkdir -p $@
	cp -v -r $(LS_HS_PATH)/* $@/
	install -m644 -C lib/liblattice_symmetries_chapel.$(SHARED_EXT) $@/lib/
ifeq ($(UNAME), Linux)
	patchelf --set-rpath '$$ORIGIN' $@/lib/liblattice_symmetries_chapel.$(SHARED_EXT)
endif

.PHONY: bundle-docker
bundle-docker: 
	mkdir -p $@
	WORKDIR=/work/distributed-matvec && \
	docker run --rm \
	  -v $$PWD/src:$$WORKDIR/src:ro \
	  -v $$PWD/Makefile:$$WORKDIR/Makefile:ro \
	  -v $(LS_HS_PATH):$$WORKDIR/lattice-symmetries-haskell:ro \
	  -v $$PWD/$@:$$WORKDIR/bundle:z \
	  twesterhout/haskell-chapel-halide \
	  bash -c 'cd /work/distributed-matvec && make PATH=/opt/chapel/bin:$$PATH OPTIMIZATION=--fast bundle'
	$(SUDO) chown -R $$USER:$$USER $@


# ifeq ($(UNAME), Darwin)
# 	# install_name_tool -id lib/liblattice_symmetries_chapel.$(SHARED_EXT) lib/liblattice_symmetries_core.$(SHARED_EXT)
# else
# 	# chpl $(CFLAGS) --library --static -o lattice_symmetries_chapel $(LIB_MODULES) $(LDFLAGS)
# 	# $(CONDA_CC) $(SHARED_FLAG) -o lib/liblattice_symmetries_chapel.$(SHARED_EXT) src/library.c lib/liblattice_symmetries_chapel.a `$$CHPL_HOME/util/config/compileline --libraries` $(LDFLAGS)
# 	# rm lib/liblattice_symmetries_chapel.a
# endif

.PHONY: release
release: lib
	mkdir -p $(DIST)/include
	mkdir -p $(DIST)/lib
	install -m644 -C third_party/include/*.h $(DIST)/include/
	# NOTE: Only copy liblattice_symmetries_haskell;
	# it is assumed that libffi will be installed via Conda
	install -m644 -C third_party/lib/liblattice_symmetries_*.$(SHARED_EXT) $(DIST)/lib/
	install -m644 -C lib/liblattice_symmetries_chapel.* $(DIST)/lib/
	rm -f $(DIST)/lib/liblattice_symmetries_chapel.a
ifeq ($(UNAME), Linux)
	find $(DIST)/lib/ -name "*.$(SHARED_EXT)" -exec patchelf --set-rpath '$$ORIGIN' {} \;
endif
	tar -cf $(DIST).tar $(DIST)
	rm -f $(DIST).tar.bz2
	bzip2 $(DIST).tar
ifneq ($(realpath $(PREFIX)), $(PWD))
	install -m644 -C $(DIST).tar.bz2 $(PREFIX)
endif
	rm -r $(DIST)

bin/TestStatesEnumeration: test/TestStatesEnumeration.chpl $(APP_MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) $(CHPL_CFLAGS) $(HDF5_CFLAGS) -o $@ --main-module $(@F) $^ $(HDF5_LDFLAGS) $(CHPL_LDFLAGS) $(LDFLAGS)

bin/TestMatrixVectorProduct: test/TestMatrixVectorProduct.chpl $(APP_MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) $(CHPL_CFLAGS) $(HDF5_CFLAGS) -o $@ --main-module $(@F) $^ $(HDF5_LDFLAGS) $(CHPL_LDFLAGS) $(LDFLAGS)

bin/BenchmarkStatesEnumeration: benchmark/BenchmarkStatesEnumeration.chpl $(LIB_MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) $(CHPL_CFLAGS) -o $@ --main-module $(@F) $^ $(CHPL_LDFLAGS) $(LDFLAGS)

bin/BenchmarkMatrixVectorProduct: benchmark/BenchmarkMatrixVectorProduct.chpl $(LIB_MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) $(CHPL_CFLAGS) -o $@ --main-module $(@F) $^ $(CHPL_LDFLAGS) $(LDFLAGS)

bin/BenchmarkBlockHashed: benchmark/BenchmarkBlockHashed.chpl $(LIB_MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) $(CHPL_CFLAGS) -o $@ --main-module $(@F) $^ $(CHPL_LDFLAGS) $(LDFLAGS)

bin/Example01: example/Example01.chpl $(APP_MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) $(HDF5_CFLAGS) -o $@ --main-module $(@F) $^ $(HDF5_LDFLAGS) $(LDFLAGS)

bin/Example02: example/Example02.chpl $(APP_MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) $(HDF5_CFLAGS) -o $@ --main-module $(@F) $^ $(HDF5_LDFLAGS) $(LDFLAGS)

bin/Example03: example/Example03.chpl $(APP_MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) $(HDF5_CFLAGS) -o $@ --main-module $(@F) $^ $(HDF5_LDFLAGS) $(LDFLAGS)

bin/Example04: example/Example04.chpl $(APP_MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) $(HDF5_CFLAGS) -o $@ --main-module $(@F) $^ $(HDF5_LDFLAGS) $(LDFLAGS)

bin/Example05: example/Example05.chpl $(APP_MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) $(HDF5_CFLAGS) -o $@ --main-module $(@F) $^ $(HDF5_LDFLAGS) $(LDFLAGS)

bin/TestRoseTree: test/TestRoseTree.chpl $(APP_MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) $(CHPL_CFLAGS) $(HDF5_CFLAGS) -o $@ --main-module $(@F) $^ $(HDF5_LDFLAGS) $(CHPL_LDFLAGS) $(LDFLAGS)

bin/dummy: src/dummy.chpl
	@mkdir -p $(@D)
	chpl $(CFLAGS) -o $@ $<

bin/atomic_queue: src/ConcurrentQueue.chpl
	@mkdir -p $(@D)
	chpl $(CFLAGS) -o $@ $<

bin/primme: src/PRIMME.chpl
	@mkdir -p $(@D)
	chpl $(CFLAGS) $(PRIMME_CFLAGS) -o $@ $< $(PRIMME_LDFLAGS)

bin/Diagonalize: src/Diagonalize.chpl src/PRIMME.chpl $(APP_MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) $(HDF5_CFLAGS) $(PRIMME_CFLAGS) -o $@ --main-module $(@F) $^ \
		$(PRIMME_LDFLAGS) \
		$(HDF5_LDFLAGS) \
		$(LDFLAGS)

# Dummy file we use to reproduce internal compiler errors in Chapel for
# submitting issues.
.PHONY: error
error: src/error.chpl
	chpl -o $@ $^

.PHONY: clean
clean:
	rm -rf bin/
