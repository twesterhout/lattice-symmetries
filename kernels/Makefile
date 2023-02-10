MAKEFLAGS += --no-builtin-rules
.SUFFIXES:

UNAME = $(shell uname)
UNAME_M = $(shell uname -m)

ifeq ($(UNAME), Darwin)
  HL_OS = osx
  SHARED_EXT = dylib
else
  HL_OS = linux
  SHARED_EXT = so
endif

# HALIDE_PATH =
# $(PWD)/../third_party/Halide
# HL_VERSION = 14.0.0
# HL_HASH = 6b9ed2afd1d6d0badf04986602c943e287d44e46
# HALIDE_DISTRIB_PATH ?= $(PWD)/../third_party/Halide-$(HL_VERSION)-$(HL_ARCH)-$(HL_BITS)-$(HL_OS)
# HL_TARGET_SUFFIX ?= -no_asserts-no_bounds_query-no_runtime
# HL_TARGET_SUFFIX ?= -no_runtime

HL_TARGET_FEATURES ?= -no_asserts-no_bounds_query-no_runtime
# HL_TARGET_FEATURES ?= -no_runtime

ifeq ($(UNAME_M),x86_64)
  _HL_ALL_TARGETS = $(HL_OS)-x86-64-sse41$(HL_TARGET_FEATURES) \
                   $(HL_OS)-x86-64-avx$(HL_TARGET_FEATURES) \
                   $(HL_OS)-x86-64-avx2-fma$(HL_TARGET_FEATURES) \
                   $(HL_OS)-x86-64$(HL_TARGET_FEATURES)
else
  _HL_ALL_TARGETS = host$(HL_TARGET_FEATURES)
endif
space := $() $()
comma := ,
HL_TARGETS ?= $(subst $(space),$(comma),$(_HL_ALL_TARGETS))

# ifeq ($(UNAME), Linux)
# CONDA_CC ?= $(shell conda run -n ci_devel bash -c "which \$${CC}")
# CONDA_PREFIX ?= $(shell conda run -n ci_devel bash -c "echo \$${CONDA_PREFIX}")
# $(info CONDA_CC=$(CONDA_CC))
# $(info CONDA_PREFIX=$(CONDA_PREFIX))
# else
# CONDA_CC = $(CC)
# CONDA_PREFIX = 
# endif

ifeq ($(UNAME_M),x86_64)
  CFLAGS ?= -O3 -DNDEBUG -march=nehalem -mtune=haswell
else
  CFLAGS ?= -O3 -DNDEBUG
endif

#
# SHELL = bash
# BIN_DIR ?= build
# GENERATOR_OUTPUTS ?= static_library,c_header,stmt_html,assembly

# CFLAGS += -I ../cbits -isystem $(HALIDE_DISTRIB_PATH)/include/ -isystem $(HALIDE_DISTRIB_PATH)/share/Halide/tools/ -Wall -Wextra -Wpedantic -Wconversion -Wdouble-promotion -Wformat=2 -Wformat-overflow=2 -Winit-self -Wimplicit-fallthrough -Wunused -Wuninitialized -Wstrict-overflow=2 -Wduplicated-branches -Wfloat-equal -Wshadow -Wundef -Wunused-macros -Wmissing-field-initializers

HL_CXXFLAGS ?= -std=c++17 -I $(HALIDE_PATH)/include/ -I $(HALIDE_PATH)/share/Halide/tools/ \
               -Wall -Wno-unused-function -Wcast-qual -Wignored-qualifiers -Wno-comment \
               -Wsign-compare -Wno-unknown-warning-option -Wno-psabi
ifeq ($(UNAME), Darwin)
  HL_CXXFLAGS += -fvisibility=hidden
endif

HL_LDFLAGS ?= -Wl,-rpath,$(HALIDE_PATH)/lib -L $(HALIDE_PATH)/lib \
              -Wl,-rpath,$(HALIDE_PATH)/lib64 -L $(HALIDE_PATH)/lib64 -lHalide -lz
ifeq ($(UNAME), Linux)
  HL_LDFLAGS += -ldl -lpthread
endif

all: $(BIN_DIR)/libkernels.a

../cabal.project.local:
	@echo "Creating cabal.project.local ..."
	@echo "-- DO NOT MODIFY (Generated automatically by kernels/Makefile)" >$@
	@echo "package lattice-symmetries-haskell" >>$@
	@echo "  extra-lib-dirs: $(realpath $(BIN_DIR))" >>$@

# ifeq ($(UNAME), Darwin)
#   SHARED_FLAG = -dynamiclib
# else
#   SHARED_FLAG = -shared
# endif

# LIB_HALIDE_STATIC = $(HALIDE_DISTRIB_PATH)/lib/libHalide.a
# LIB_HALIDE = $(HALIDE_PATH)/lib/libHalide.$(SHARED_EXT)

# GENERATOR_DEPS ?= $(LIB_AUTOSCHEDULER) \
# 		  $(LIB_HALIDE) \
# 		  $(HALIDE_DISTRIB_PATH)/include/Halide.h \
# 		  $(HALIDE_DISTRIB_PATH)/share/Halide/tools/GenGen.cpp
# GENERATOR_DEPS_STATIC ?= $(LIB_HALIDE_STATIC) \
# 			 $(HALIDE_DISTRIB_PATH)/include/Halide.h \
# 			 $(HALIDE_DISTRIB_PATH)/share/Halide/tools/GenGen.cpp

# Generators which use autoscheduler plugin need to specify the linker where to find libHalide.so required by the plugin.
# LIBHALIDE_LDFLAGS ?= -Wl,-rpath,$(dir $(LIB_HALIDE)) -L $(dir $(LIB_HALIDE)) -lHalide -lz $(LDFLAGS)
# LIBHALIDE_LDFLAGS_STATIC ?= $(LIB_HALIDE_STATIC) -lz $(LDFLAGS)

# Autoschedulers.
# AUTOSCHEDULER ?= adams2019
# ifneq ($(AUTOSCHEDULER),)
#   LIB_AUTOSCHEDULER ?= $(HALIDE_DISTRIB_PATH)/lib/libautoschedule_$(AUTOSCHEDULER).$(SHARED_EXT)
#   ifeq ($(UNAME), Darwin)
#     LIBHALIDE_LDFLAGS += -Wl,-force_load $(HALIDE_DISTRIB_PATH)/lib/libautoschedule_$(AUTOSCHEDULER).$(SHARED_EXT)
#   else
#     LIBHALIDE_LDFLAGS += -Wl,--no-as-needed -lautoschedule_$(AUTOSCHEDULER) -Wl,--as-needed
#   endif
# endif


# OPS = $(addsuffix /state_info_symmetric.a,$(addprefix $(BIN_DIR)/,$(HL_ALL_TARGETS))) \
#       $(addsuffix /state_info_antisymmetric.a,$(addprefix $(BIN_DIR)/,$(HL_ALL_TARGETS))) \
#       $(addsuffix /state_info_general.a,$(addprefix $(BIN_DIR)/,$(HL_ALL_TARGETS))) \
#       $(addsuffix /is_representative_symmetric.a,$(addprefix $(BIN_DIR)/,$(HL_ALL_TARGETS))) \
#       $(addsuffix /is_representative_antisymmetric.a,$(addprefix $(BIN_DIR)/,$(HL_ALL_TARGETS))) \
#       $(addsuffix /is_representative_general.a,$(addprefix $(BIN_DIR)/,$(HL_ALL_TARGETS))) \
#       $(addsuffix /apply_off_diagonal.a,$(addprefix $(BIN_DIR)/,$(HL_ALL_TARGETS)))

# all: $(BIN_DIR)/kernels.o
# $(OPS) $(BIN_DIR)/runtime.a $(BIN_DIR)/kernels.h

# all: ../cabal.project.local
# $(BIN_DIR)/liblattice_symmetries_core.$(SHARED_EXT) 

# HS_OBJECT_FILES = $(OPS) $(BIN_DIR)/runtime.a
# HS_LINK_FLAGS = $(foreach obj,$(HS_OBJECT_FILES),$(realpath $(obj)))
# HS_INCLUDES = $(realpath $(BIN_DIR)) \
#               $(HALIDE_DISTRIB_PATH)/include \
#               $(HALIDE_DISTRIB_PATH)/share/Halide/tools

# $(BIN_DIR)/liblattice_symmetries_core.$(SHARED_EXT): $(BIN_DIR)/reference.o $(BIN_DIR)/indexing.o $(BIN_DIR)/kernels.o $(OPS) $(BIN_DIR)/runtime.a 
# 	@mkdir -p $(@D)
# 	$(CONDA_CC) $(SHARED_FLAG) -o $(BIN_DIR)/liblattice_symmetries_core.$(SHARED_EXT) $^ $(LDFLAGS)
# ifeq ($(UNAME), Darwin)
# 	install_name_tool -id $(CURDIR)/$(BIN_DIR)/liblattice_symmetries_core.$(SHARED_EXT) $(BIN_DIR)/liblattice_symmetries_core.$(SHARED_EXT)
# endif

# $(BIN_DIR)/reference.o: reference.c ../cbits/lattice_symmetries_haskell.h
# 	@mkdir -p $(@D)
# 	$(CC) -fPIC $(CFLAGS) -I$(PWD)/../cbits -o $@ -c $<

# $(BIN_DIR)/indexing.o: indexing.c ../cbits/lattice_symmetries_haskell.h
# 	@mkdir -p $(@D)
# 	$(CC) -fPIC $(CFLAGS) -I$(PWD)/../cbits -o $@ -c $<

# $(BIN_DIR)/binary_search.o: binary_search.c ../cbits/lattice_symmetries_haskell.h
# 	@mkdir -p $(@D)
# 	$(CC) -fPIC $(CFLAGS) -I$(PWD)/../cbits -o $@ -c $<

# $(BIN_DIR)/kernels.o: kernels.c $(BIN_DIR)/kernels.h ../cbits/lattice_symmetries_haskell.h
# 	@mkdir -p $(@D)
# 	$(CC) -fPIC $(CFLAGS) -I$(PWD)/../cbits -I$(BIN_DIR) -o $@ -c $< 

# $(BIN_DIR)/kernels.h: $(OPS)
# 	@mkdir -p $(@D)
# 	echo "#pragma once" > $(BIN_DIR)/kernels.h
# 	echo "" >> $(BIN_DIR)/kernels.h
# 	find $(BIN_DIR) -type f -name "*.h" | grep --invert-match "$(BIN_DIR)/kernels.h" | sed -E 's:$(BIN_DIR)/(.*):#include "\1":' >> $(BIN_DIR)/kernels.h

# clean_kernel_name = $(shell echo $(1) | sed -E 's/(no_runtime|no_asserts|no_bounds_query|$(HL_OS)|$(HL_ARCH)-$(HL_BITS))//g;s/[-]+/_/g;s/_$$//g')

# GENERATOR_ENVIRONMENT = LD_LIBRARY_PATH=$(HALIDE_PATH)/lib DYLD_LIBRARY_PATH=$(HALIDE_PATH)/lib


# Build the Halide generator for our kernels
$(BIN_DIR)/generator: generator.cpp 
	@echo Building the generator ...
	@mkdir -p $(@D)
	$(CXX) $(HL_CXXFLAGS) $< $(HALIDE_PATH)/share/tools/GenGen.cpp -o $@ $(HL_LDFLAGS)

# Build common runtime
$(BIN_DIR)/runtime.o: $(BIN_DIR)/generator
	@echo Building the runtime ...
	@mkdir -p $(@D)
	$^ -r runtime -o $(@D) -e object target=$(HL_TARGETS)

$(BIN_DIR)/libkernels.a: $(BIN_DIR)/runtime.o \
					     $(BIN_DIR)/kernels.o \
					     $(BIN_DIR)/state_info_symmetric.o \
					     $(BIN_DIR)/state_info_antisymmetric.o \
					     $(BIN_DIR)/state_info_general.o \
					     $(BIN_DIR)/is_representative_symmetric.o \
					     $(BIN_DIR)/is_representative_antisymmetric.o \
					     $(BIN_DIR)/is_representative_general.o
	$(AR) rcs $@ $^

$(BIN_DIR)/kernels.o: kernels.c kernels.h ../cbits/lattice_symmetries_haskell.h \
					  $(BIN_DIR)/state_info_symmetric.o \
					  $(BIN_DIR)/state_info_antisymmetric.o \
					  $(BIN_DIR)/state_info_general.o \
					  $(BIN_DIR)/is_representative_symmetric.o \
					  $(BIN_DIR)/is_representative_antisymmetric.o \
					  $(BIN_DIR)/is_representative_general.o
	$(CC) -fPIC $(CFLAGS) -I../cbits -I$(BIN_DIR) -I$(HALIDE_PATH)/include -o $@ -c $<
	
$(BIN_DIR)/%_symmetric.o: $(BIN_DIR)/generator
	$^ -g $* -o $(@D) -e object,c_header \
		-n $*_symmetric \
		-f ls_internal_$*_symmetric_kernel \
		spin_inversion=1 auto_schedule=false \
		target=$(HL_TARGETS)
	$(LD) -relocatable $(BIN_DIR)/$*_symmetric*.o -o $@

$(BIN_DIR)/%_antisymmetric.o: $(BIN_DIR)/generator
	$^ -g $* -o $(@D) -e object,c_header \
		-n $*_antisymmetric \
		-f ls_internal_$*_antisymmetric_kernel \
		spin_inversion=-1 auto_schedule=false \
		target=$(HL_TARGETS)
	$(LD) -relocatable $(BIN_DIR)/$*_antisymmetric*.o -o $@

$(BIN_DIR)/%_general.o: $(BIN_DIR)/generator
	$^ -g $* -o $(@D) -e object,c_header \
		-n $*_general \
		-f ls_internal_$*_general_kernel \
		spin_inversion=0 auto_schedule=false \
		target=$(HL_TARGETS)
	$(LD) -relocatable $(BIN_DIR)/$*_general*.o -o $@

# $(BIN_DIR)/%/apply_off_diagonal.a: $(BIN_DIR)/generator
# 	@mkdir -p $(@D)
# 	$(GENERATOR_ENVIRONMENT) \
# 	$^ -g apply_off_diagonal -o $(@D) -e $(GENERATOR_OUTPUTS) \
# 		-n apply_off_diagonal \
# 		-f ls_internal_apply_off_diag_kernel$(call clean_kernel_name,$*) \
# 		number_words=1 auto_schedule=false \
# 		target=$*

# $(BIN_DIR)/%/state_info_symmetric.a: $(BIN_DIR)/generator
# 	@mkdir -p $(@D)
# 	$(GENERATOR_ENVIRONMENT) \
# 	$^ -g state_info -o $(@D) -e $(GENERATOR_OUTPUTS) \
# 		-n state_info_symmetric \
# 		-f ls_internal_state_info_symmetric_kernel$(call clean_kernel_name,$*) \
# 		spin_inversion=1 auto_schedule=false \
# 		target=$*

# $(BIN_DIR)/%/state_info_antisymmetric.a: $(BIN_DIR)/generator
# 	@mkdir -p $(@D)
# 	$(GENERATOR_ENVIRONMENT) \
# 	$^ -g state_info -o $(@D) -e $(GENERATOR_OUTPUTS) \
# 		-n state_info_antisymmetric \
# 		-f ls_internal_state_info_antisymmetric_kernel$(call clean_kernel_name,$*) \
# 		spin_inversion=-1 auto_schedule=false \
# 		target=$*

# $(BIN_DIR)/%/state_info_general.a: $(BIN_DIR)/generator
# 	@mkdir -p $(@D)
# 	$(GENERATOR_ENVIRONMENT) \
# 	$^ -g state_info -o $(@D) -e $(GENERATOR_OUTPUTS) \
# 		-n state_info_general \
# 		-f ls_internal_state_info_general_kernel$(call clean_kernel_name,$*) \
# 		spin_inversion=0 auto_schedule=false \
# 		target=$*

# $(BIN_DIR)/%/is_representative_symmetric.a: $(BIN_DIR)/generator
# 	@mkdir -p $(@D)
# 	$(GENERATOR_ENVIRONMENT) \
# 	$^ -g is_representative -o $(@D) -e $(GENERATOR_OUTPUTS) \
# 		-n is_representative_symmetric \
# 		-f ls_internal_is_representative_symmetric_kernel$(call clean_kernel_name,$*) \
# 		spin_inversion=1 auto_schedule=false \
# 		target=$*

# $(BIN_DIR)/%/is_representative_antisymmetric.a: $(BIN_DIR)/generator
# 	@mkdir -p $(@D)
# 	$(GENERATOR_ENVIRONMENT) \
# 	$^ -g is_representative -o $(@D) -e $(GENERATOR_OUTPUTS) \
# 		-n is_representative_antisymmetric \
# 		-f ls_internal_is_representative_antisymmetric_kernel$(call clean_kernel_name,$*) \
# 		spin_inversion=-1 auto_schedule=false \
# 		target=$*

# $(BIN_DIR)/%/is_representative_general.a: $(BIN_DIR)/generator
# 	@mkdir -p $(@D)
# 	$(GENERATOR_ENVIRONMENT) \
# 	$^ -g is_representative -o $(@D) -e $(GENERATOR_OUTPUTS) \
# 		-n is_representative_general \
# 		-f ls_internal_is_representative_general_kernel$(call clean_kernel_name,$*) \
# 		spin_inversion=0 auto_schedule=false \
# 		target=$*


# Halide: $(PWD)/../third_party/Halide-$(HL_VERSION)-$(HL_ARCH)-$(HL_BITS)-$(HL_OS)

# $(PWD)/../third_party/Halide-$(HL_VERSION)-$(HL_ARCH)-$(HL_BITS)-$(HL_OS)-$(HL_HASH).$(HL_ARCHIVE):
# 	@echo Downloading Halide ...
# 	@mkdir -p $(@D)
# 	cd ../third_party && wget -q https://github.com/halide/Halide/releases/download/v$(HL_VERSION)/Halide-$(HL_VERSION)-$(HL_ARCH)-$(HL_BITS)-$(HL_OS)-$(HL_HASH).$(HL_ARCHIVE)
 
# $(PWD)/../third_party/Halide-$(HL_VERSION)-$(HL_ARCH)-$(HL_BITS)-$(HL_OS): $(PWD)/../third_party/Halide-$(HL_VERSION)-$(HL_ARCH)-$(HL_BITS)-$(HL_OS)-$(HL_HASH).$(HL_ARCHIVE)
# 	cd ../third_party && tar -xf Halide-$(HL_VERSION)-$(HL_ARCH)-$(HL_BITS)-$(HL_OS)-$(HL_HASH).$(HL_ARCHIVE)

.PHONY: clean
clean:
	rm -rf $(BIN_DIR) ../cabal.project.local
