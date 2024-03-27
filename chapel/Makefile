.RULES:
.POSIX:

PREFIX ?= /usr/local
TEST_DATA ?= $(PWD)
CHPL_ARGS ?=
CHPL_CFLAGS ?= --fast
CHPL_LDFLAGS = -llattice_symmetries_haskell

PRIMME_CFLAGS = -I$(PRIMME_PATH)/include
PRIMME_LDFLAGS = -L$(PRIMME_PATH)/lib -lprimme -llapacke -lopenblas -lm -lgomp -lpthread

HALIDE_CFLAGS = -I$(HALIDE_PATH)/include

ifeq ($(shell uname -s), Linux)
  HDF5_CFLAGS ?= $(shell pkg-config --cflags hdf5)
  HDF5_LDFLAGS ?= -lhdf5_hl $(shell pkg-config --libs hdf5)
  SHARED_EXT = so
else
  HDF5_CFLAGS ?=
  HDF5_LDFLAGS ?= -lhdf5_hl -lhdf5
  SHARED_EXT = dylib
endif

LIB_MODULES = src/library.c \
              src/CommonParameters.chpl \
              src/FFI.chpl \
              src/Utils.chpl \
              src/ForeignTypes.chpl \
              src/Timing.c src/Timing.chpl \
              src/Vector.chpl \
              src/BatchedOperator.chpl \
              src/ConcurrentQueue.c src/ConcurrentQueue.chpl \
              src/StatesEnumeration.chpl \
              src/CSR.chpl \
              src/MatrixVectorProduct.chpl \
              src/LatticeSymmetries.chpl

.PHONY: all
all: lib

lib: lib/liblattice_symmetries_chapel.$(SHARED_EXT)

lib/liblattice_symmetries_chapel.$(SHARED_EXT): $(LIB_MODULES)
	@mkdir -p $(@D)
	@rm -rf $(@D)/*
	CHPL_LIB_PIC=pic chpl \
		$(CHPL_CFLAGS) \
		$(HALIDE_CFLAGS) \
		--library --dynamic -o lattice_symmetries_chapel \
		--ccflags -DLS_CHPL_LIBRARY=1 \
		$^ \
		$(CHPL_LDFLAGS) -llattice_symmetries_haskell

install: lib
	@mkdir -p $(PREFIX)/lib
	@mkdir -p $(PREFIX)/include
	install -m 755 lib/liblattice_symmetries_chapel.$(SHARED_EXT) $(PREFIX)/lib/
	install -m 644 src/lattice_symmetries_chapel.h $(PREFIX)/include/

# .PHONY: benchmark-states-enumeration
# benchmark-states-enumeration: bin/TestStatesEnumeration
# 	# $(CHPL_LIBS) $< $(CHPL_ARGS) --kBasis data/heisenberg_pyrochlore_2x2x2.yaml --kRepresentatives data/large-scale/construction/heisenberg_pyrochlore_2x2x2.h5
# 	# $(CHPL_LIBS) $< $(CHPL_ARGS) --kBasis data/heisenberg_kagome_36.yaml --kRepresentatives data/large-scale/construction/heisenberg_kagome_36.h5
# 	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_square_6x6.yaml --kRepresentatives data/large-scale/matvec/heisenberg_square_6x6.h5
# 
# .PHONY: benchmark-matrix-vector-product
# benchmark-matrix-vector-product: bin/TestMatrixVectorProduct
# 	# $(CHPL_LIBS) $< $(CHPL_ARGS) --kBasis data/heisenberg_pyrochlore_2x2x2.yaml --kRepresentatives data/heisenberg_pyrochlore_2x2x2.h5
# 	# $(CHPL_LIBS) $< $(CHPL_ARGS) --kBasis data/heisenberg_kagome_36.yaml --kRepresentatives data/heisenberg_kagome_36.h5
# 	$(CHPL_LIBS) $< $(CHPL_ARGS) --kHamiltonian data/heisenberg_square_6x6.yaml --kVectors data/large-scale/matvec/heisenberg_square_6x6.h5

bin/TestConcurrentQueue: test/TestConcurrentQueue.chpl src/ConcurrentQueue.chpl src/ConcurrentQueue.c
	@mkdir -p $(@D)
	chpl $(CHPL_CFLAGS) -o $@ --main-module $(@F) $^ $(CHPL_LDFLAGS) $(LDFLAGS)

bin/TestUtils: test/TestUtils.chpl src/CommonParameters.chpl src/FFI.chpl src/Utils.chpl src/Timing.c src/Timing.chpl src/ForeignTypes.chpl src/Vector.chpl
	@mkdir -p $(@D)
	chpl $(CHPL_CFLAGS) $(HALIDE_CFLAGS) -o $@ --main-module $(@F) $^ $(CHPL_LDFLAGS) $(LDFLAGS)

bin/TestBatchedOperator: test/TestBatchedOperator.chpl \
                         src/FFI.chpl \
                         src/Utils.chpl \
                         src/ForeignTypes.c src/ForeignTypes.chpl \
                         src/Timing.c src/Timing.chpl \
                         src/CommonParameters.chpl \
                         src/BatchedOperator.chpl \
                         src/ConcurrentQueue.c src/ConcurrentQueue.chpl \
                         src/MatrixVectorProduct.chpl
	@mkdir -p $(@D)
	chpl $(CHPL_CFLAGS) $(HALIDE_CFLAGS) -I src/ -o $@ --main-module $(@F) $^ $(CHPL_LDFLAGS) $(LDFLAGS)

bin/TestSingleLocaleMatrixVector: test/TestSingleLocaleMatrixVector.chpl $(LIB_MODULES)
	@mkdir -p $(@D)
	chpl $(CHPL_CFLAGS) $(HALIDE_CFLAGS) $(HDF5_CFLAGS) -I src/ -o $@ --main-module $(@F) $^ $(HDF5_CFLAGS) $(CHPL_LDFLAGS) $(LDFLAGS)


bin/TestToCsr: test/TestToCsr.chpl $(LIB_MODULES)
	@mkdir -p $(@D)
	chpl $(CHPL_CFLAGS) $(HALIDE_CFLAGS) $(HDF5_CFLAGS) -I src/ -o $@ --main-module $(@F) $^ $(HDF5_CFLAGS) $(CHPL_LDFLAGS) $(LDFLAGS)

bin/TestStatesEnumeration: test/TestStatesEnumeration.chpl $(LIB_MODULES)
	@mkdir -p $(@D)
	chpl $(CHPL_CFLAGS) $(HALIDE_CFLAGS) -I src/ -o $@ --main-module $(@F) $^ $(CHPL_LDFLAGS) $(LDFLAGS)

bin/TestMatrixVectorProduct: test/TestMatrixVectorProduct.chpl $(APP_MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) $(CHPL_CFLAGS) $(HDF5_CFLAGS) -o $@ --main-module $(@F) $^ $(HDF5_LDFLAGS) $(CHPL_LDFLAGS) $(LDFLAGS)

bin/BenchmarkStatesEnumeration: benchmark/BenchmarkStatesEnumeration.chpl $(LIB_MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) $(CHPL_CFLAGS) -o $@ --main-module $(@F) $^ $(CHPL_LDFLAGS) $(LDFLAGS)

bin/BenchmarkMatrixVectorProduct: benchmark/BenchmarkMatrixVectorProduct.chpl $(LIB_MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) $(CHPL_CFLAGS) -o $@ --main-module $(@F) $^ $(CHPL_LDFLAGS) $(LDFLAGS)

bin/BenchmarkSingleLocaleMatrixVector: benchmark/BenchmarkSingleLocaleMatrixVector.chpl $(LIB_MODULES)
	@mkdir -p $(@D)
	chpl $(CHPL_CFLAGS) $(HALIDE_CFLAGS) -I src/ -o $@ --main-module $(@F) $^ $(CHPL_LDFLAGS) $(LDFLAGS)

bin/BenchmarkBlockHashed: benchmark/BenchmarkBlockHashed.chpl $(LIB_MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) $(CHPL_CFLAGS) -o $@ --main-module $(@F) $^ $(CHPL_LDFLAGS) $(LDFLAGS)

.PHONY: check-states-enumeration
check-states-enumeration: bin/TestStatesEnumeration $(TEST_DATA)/data/matvec
	for yaml_file in $$(cd data/ && ls *.yaml); do \
		if test -f "$(TEST_DATA)/data/matvec/$${yaml_file%.yaml}.h5"; then \
			echo "Testing $$yaml_file ..."; \
			$(CHPL_LIBS) $< $(CHPL_ARGS) \
				--kHamiltonian "data/$$yaml_file" \
				--kRepresentatives "$(TEST_DATA)/data/matvec/$${yaml_file%.yaml}.h5"; \
		fi; \
	done

.PHONY: check-matrix-vector-product
check-matrix-vector-product: bin/TestSingleLocaleMatrixVector
	find ../test/random_quspin -name "*_expr.json" | xargs -P 4 -I {} ./bin/TestSingleLocaleMatrixVector --kHamiltonian {}


TEST_DATA_URL = https://surfdrive.surf.nl/files/index.php/s/OK5527Awfgl1hT2/download

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


# .PHONY: bundle
# bundle: lib
# 	mkdir -p $@
# 	cp -v -r $(LS_HS_PATH)/* $@/
# 	install -m644 -C lib/liblattice_symmetries_chapel.$(SHARED_EXT) $@/lib/
# ifeq ($(UNAME), Linux)
# 	patchelf --set-rpath '$$ORIGIN' $@/lib/liblattice_symmetries_chapel.$(SHARED_EXT)
# endif

# .PHONY: bundle-docker
# bundle-docker: 
# 	mkdir -p $@
# 	WORKDIR=/work/distributed-matvec && \
# 	docker run --rm \
# 	  -v $$PWD/src:$$WORKDIR/src:ro \
# 	  -v $$PWD/Makefile:$$WORKDIR/Makefile:ro \
# 	  -v $(LS_HS_PATH):$$WORKDIR/lattice-symmetries-haskell:ro \
# 	  -v $$PWD/$@:$$WORKDIR/bundle:z \
# 	  twesterhout/haskell-chapel-halide \
# 	  bash -c 'cd /work/distributed-matvec && make PATH=/opt/chapel/bin:$$PATH OPTIMIZATION=--fast bundle'
# 	$(SUDO) chown -R $$USER:$$USER $@


# ifeq ($(UNAME), Darwin)
# 	# install_name_tool -id lib/liblattice_symmetries_chapel.$(SHARED_EXT) lib/liblattice_symmetries_core.$(SHARED_EXT)
# else
# 	# chpl $(CFLAGS) --library --static -o lattice_symmetries_chapel $(LIB_MODULES) $(LDFLAGS)
# 	# $(CONDA_CC) $(SHARED_FLAG) -o lib/liblattice_symmetries_chapel.$(SHARED_EXT) src/library.c lib/liblattice_symmetries_chapel.a `$$CHPL_HOME/util/config/compileline --libraries` $(LDFLAGS)
# 	# rm lib/liblattice_symmetries_chapel.a
# endif

# .PHONY: release
# release: lib
# 	mkdir -p $(DIST)/include
# 	mkdir -p $(DIST)/lib
# 	install -m644 -C third_party/include/*.h $(DIST)/include/
# 	# NOTE: Only copy liblattice_symmetries_haskell;
# 	# it is assumed that libffi will be installed via Conda
# 	install -m644 -C third_party/lib/liblattice_symmetries_*.$(SHARED_EXT) $(DIST)/lib/
# 	install -m644 -C lib/liblattice_symmetries_chapel.* $(DIST)/lib/
# 	rm -f $(DIST)/lib/liblattice_symmetries_chapel.a
# ifeq ($(UNAME), Linux)
# 	find $(DIST)/lib/ -name "*.$(SHARED_EXT)" -exec patchelf --set-rpath '$$ORIGIN' {} \;
# endif
# 	tar -cf $(DIST).tar $(DIST)
# 	rm -f $(DIST).tar.bz2
# 	bzip2 $(DIST).tar
# ifneq ($(realpath $(PREFIX)), $(PWD))
# 	install -m644 -C $(DIST).tar.bz2 $(PREFIX)
# endif
# 	rm -r $(DIST)


# bin/Example01: example/Example01.chpl $(APP_MODULES)
# 	@mkdir -p $(@D)
# 	chpl $(CFLAGS) $(HDF5_CFLAGS) -o $@ --main-module $(@F) $^ $(HDF5_LDFLAGS) $(LDFLAGS)

# bin/Example02: example/Example02.chpl $(APP_MODULES)
# 	@mkdir -p $(@D)
# 	chpl $(CFLAGS) $(HDF5_CFLAGS) -o $@ --main-module $(@F) $^ $(HDF5_LDFLAGS) $(LDFLAGS)

# bin/Example03: example/Example03.chpl $(APP_MODULES)
# 	@mkdir -p $(@D)
# 	chpl $(CFLAGS) $(HDF5_CFLAGS) -o $@ --main-module $(@F) $^ $(HDF5_LDFLAGS) $(LDFLAGS)

# bin/Example04: example/Example04.chpl $(APP_MODULES)
# 	@mkdir -p $(@D)
# 	chpl $(CFLAGS) $(HDF5_CFLAGS) -o $@ --main-module $(@F) $^ $(HDF5_LDFLAGS) $(LDFLAGS)

# bin/Example05: example/Example05.chpl $(APP_MODULES)
# 	@mkdir -p $(@D)
# 	chpl $(CFLAGS) $(HDF5_CFLAGS) -o $@ --main-module $(@F) $^ $(HDF5_LDFLAGS) $(LDFLAGS)

# bin/TestRoseTree: test/TestRoseTree.chpl $(APP_MODULES)
# 	@mkdir -p $(@D)
# 	chpl $(CFLAGS) $(CHPL_CFLAGS) $(HDF5_CFLAGS) -o $@ --main-module $(@F) $^ $(HDF5_LDFLAGS) $(CHPL_LDFLAGS) $(LDFLAGS)

# bin/dummy: src/dummy.chpl
# 	@mkdir -p $(@D)
# 	chpl $(CFLAGS) -o $@ $<

# bin/atomic_queue: src/ConcurrentQueue.chpl
# 	@mkdir -p $(@D)
# 	chpl $(CFLAGS) -o $@ $<

# bin/primme: src/PRIMME.chpl
# 	@mkdir -p $(@D)
# 	chpl $(CFLAGS) $(PRIMME_CFLAGS) -o $@ $< $(PRIMME_LDFLAGS)

bin/Diagonalize: src/Diagonalize.chpl src/PRIMME.chpl $(APP_MODULES)
	@mkdir -p $(@D)
	chpl $(CFLAGS) $(HDF5_CFLAGS) $(PRIMME_CFLAGS) -o $@ --main-module $(@F) $^ \
		$(PRIMME_LDFLAGS) \
		$(HDF5_LDFLAGS) \
		$(LDFLAGS)

# Dummy file we use to reproduce internal compiler errors in Chapel for
# submitting issues.
# .PHONY: error
# error: src/error.chpl
# 	chpl -o $@ $^

.PHONY: clean
clean:
	rm -rf bin/
