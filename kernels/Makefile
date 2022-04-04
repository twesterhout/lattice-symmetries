MAKEFLAGS += --no-builtin-rules
.SUFFIXES:

HL_ARCH ?= x86
HL_BITS ?= 64

UNAME = $(shell uname)
ifeq ($(UNAME), Linux)
HL_OS = linux
else
ifeq ($(UNAME), Darwin)
HL_OS = osx
else # This will probably never happen :P
HL_OS = windows
endif
endif

HALIDE_DISTRIB_PATH ?= $(realpath ../third_party/Halide-13.0.4-$(HL_ARCH)-$(HL_BITS)-$(HL_OS))
ifeq (,$(HALIDE_DISTRIB_PATH))
$(error Could not find a Halide distribution. Before building, set HALIDE_DISTRIB_PATH to the root of an existing Halide distribution)
endif

# HL_TARGET_SUFFIX ?= -no_asserts-no_bounds_query-no_runtime
HL_TARGET_SUFFIX ?= -no_runtime
HL_BASE_TARGET ?= $(HL_OS)-$(HL_ARCH)-$(HL_BITS)$(HL_TARGET_SUFFIX)

UNAME_M = $(shell uname -m)
ifeq ($(UNAME_M),x86_64) #-----------
  HL_ALL_TARGETS ?= $(HL_OS)-x86-64-sse41$(HL_TARGET_SUFFIX) \
                    $(HL_OS)-x86-64-avx$(HL_TARGET_SUFFIX) \
                    $(HL_OS)-x86-64-avx2-fma$(HL_TARGET_SUFFIX)
else #-------------------------------
  HL_ALL_TARGETS ?= host$(HL_TARGET_SUFFIX)
endif #------------------------------
HL_ALL_TARGETS += $(HL_BASE_TARGET)

CXX ?= c++
LDFLAGS ?=
SHELL = bash
BIN_DIR ?= build
SANITIZER_FLAGS ?=
GENERATOR_OUTPUTS ?= static_library,c_header,stmt_html,assembly

LDFLAGS += -ldl -lpthread -lz
CFLAGS += -I $(HALIDE_DISTRIB_PATH)/include/ -I $(HALIDE_DISTRIB_PATH)/share/Halide/tools/
CXXFLAGS += -std=c++17 -I $(HALIDE_DISTRIB_PATH)/include/ -I $(HALIDE_DISTRIB_PATH)/share/Halide/tools/ $(SANITIZER_FLAGS) -Wall -Werror -Wno-unused-function -Wcast-qual -Wignored-qualifiers -Wno-comment -Wsign-compare -Wno-unknown-warning-option -Wno-psabi

ifeq ($(UNAME), Darwin)
  CXXFLAGS += -fvisibility=hidden
endif

ifeq ($(UNAME), Darwin)
  SHARED_EXT = dylib
else
  SHARED_EXT = so
endif

ifeq ($(UNAME), Darwin)
  SHARED_FLAG = -dynamiclib
else
  SHARED_FLAG = -shared
endif

# Autoschedulers. Mullapudi2016 is currently the default, because it's fast.
AUTOSCHEDULER ?= mullapudi2016
ifneq ($(AUTOSCHEDULER),)
  LIB_AUTOSCHEDULER ?= $(HALIDE_DISTRIB_PATH)/lib/libautoschedule_$(AUTOSCHEDULER).$(SHARED_EXT)
ifeq ($(UNAME), Darwin)
  LIBHALIDE_LDFLAGS += -Wl,-force_load $(HALIDE_DISTRIB_PATH)/lib/libautoschedule_$(AUTOSCHEDULER).$(SHARED_EXT)
else
  LIBHALIDE_LDFLAGS += -Wl,--no-as-needed -lautoschedule_$(AUTOSCHEDULER) -Wl,--as-needed
endif
endif

LIB_HALIDE_STATIC = $(HALIDE_DISTRIB_PATH)/lib/libHalide.a
LIB_HALIDE = $(HALIDE_DISTRIB_PATH)/lib/libHalide.$(SHARED_EXT)

GENERATOR_DEPS ?= $(LIB_AUTOSCHEDULER) $(LIB_HALIDE) $(HALIDE_DISTRIB_PATH)/include/Halide.h $(HALIDE_DISTRIB_PATH)/share/Halide/tools/GenGen.cpp
GENERATOR_DEPS_STATIC ?= $(LIB_HALIDE_STATIC) $(HALIDE_DISTRIB_PATH)/include/Halide.h $(HALIDE_DISTRIB_PATH)/share/Halide/tools/GenGen.cpp

# Generators which use autoscheduler plugin need to specify the linker where to find libHalide.so required by the plugin.
LIBHALIDE_LDFLAGS ?= -Wl,-rpath,$(dir $(LIB_HALIDE)) -L $(dir $(LIB_HALIDE)) -lHalide $(LDFLAGS)
LIBHALIDE_LDFLAGS_STATIC ?= $(LIB_HALIDE_STATIC) $(LDFLAGS)


OPS = $(addsuffix /state_info_symmetric.a,$(addprefix $(BIN_DIR)/,$(HL_ALL_TARGETS))) \
      $(addsuffix /state_info_antisymmetric.a,$(addprefix $(BIN_DIR)/,$(HL_ALL_TARGETS))) \
      $(addsuffix /state_info_general.a,$(addprefix $(BIN_DIR)/,$(HL_ALL_TARGETS))) \
      $(addsuffix /is_representative_symmetric.a,$(addprefix $(BIN_DIR)/,$(HL_ALL_TARGETS))) \
      $(addsuffix /is_representative_antisymmetric.a,$(addprefix $(BIN_DIR)/,$(HL_ALL_TARGETS))) \
      $(addsuffix /is_representative_general.a,$(addprefix $(BIN_DIR)/,$(HL_ALL_TARGETS)))

# all: $(BIN_DIR)/kernels.o
# $(OPS) $(BIN_DIR)/runtime.a $(BIN_DIR)/kernels.h

all: $(BIN_DIR)/liblattice_symmetries_core.$(SHARED_EXT)

$(BIN_DIR)/liblattice_symmetries_core.$(SHARED_EXT): $(BIN_DIR)/reference.o $(BIN_DIR)/indexing.o $(BIN_DIR)/kernels.o $(OPS)
	@mkdir -p $(@D)
	$(CC) $(SHARED_FLAG) -o $(BIN_DIR)/liblattice_symmetries_core.$(SHARED_EXT) $^
ifeq ($(UNAME), Darwin)
	install_name_tool -id $(CURDIR)/$(BIN_DIR)/liblattice_symmetries_core.$(SHARED_EXT) $(BIN_DIR)/liblattice_symmetries_core.$(SHARED_EXT)
endif

$(BIN_DIR)/reference.o: ../cbits/reference.c ../cbits/lattice_symmetries_haskell.h
	@mkdir -p $(@D)
	$(CC) -fPIC $(CFLAGS) -I$(PWD)/../cbits `pkg-config --cflags lattice_symmetries` -o $@ -c $<

$(BIN_DIR)/indexing.o: ../cbits/indexing.c ../cbits/lattice_symmetries_haskell.h
	@mkdir -p $(@D)
	$(CC) -fPIC $(CFLAGS) -I$(PWD)/../cbits `pkg-config --cflags lattice_symmetries` -o $@ -c $<

$(BIN_DIR)/kernels.o: kernels.c $(BIN_DIR)/kernels.h
	@mkdir -p $(@D)
	$(CC) -fPIC $(CFLAGS) -I$(PWD)/../cbits `pkg-config --cflags lattice_symmetries` -I$(BIN_DIR) -o $@ -c $< 

$(BIN_DIR)/kernels.h: $(OPS)
	@mkdir -p $(@D)
	echo "#pragma once" > $(BIN_DIR)/kernels.h
	echo "" >> $(BIN_DIR)/kernels.h
	find $(BIN_DIR) -type f -name "*.h" | grep --invert-match "$(BIN_DIR)/kernels.h" | sed -E 's:$(BIN_DIR)/(.*):#include "\1":' >> $(BIN_DIR)/kernels.h

clean_kernel_name = $(shell echo $(1) | sed -E 's/(no_runtime|no_asserts|no_bounds_query|$(HL_OS)|$(HL_ARCH)-$(HL_BITS))//g;s/[-]+/_/g;s/_$$//g')

$(BIN_DIR)/%/state_info_symmetric.a: $(BIN_DIR)/generator
	@mkdir -p $(@D)
	$^ -g state_info -o $(@D) -e $(GENERATOR_OUTPUTS) \
		-n state_info_symmetric \
		-f ls_internal_state_info_symmetric_kernel$(call clean_kernel_name,$*) \
		spin_inversion=1 auto_schedule=false \
		target=$*

$(BIN_DIR)/%/state_info_antisymmetric.a: $(BIN_DIR)/generator
	@mkdir -p $(@D)
	$^ -g state_info -o $(@D) -e $(GENERATOR_OUTPUTS) \
		-n state_info_antisymmetric \
		-f ls_internal_state_info_antisymmetric_kernel$(call clean_kernel_name,$*) \
		spin_inversion=-1 auto_schedule=false \
		target=$*

$(BIN_DIR)/%/state_info_general.a: $(BIN_DIR)/generator
	@mkdir -p $(@D)
	$^ -g state_info -o $(@D) -e $(GENERATOR_OUTPUTS) \
		-n state_info_general \
		-f ls_internal_state_info_general_kernel$(call clean_kernel_name,$*) \
		spin_inversion=0 auto_schedule=false \
		target=$*

$(BIN_DIR)/%/is_representative_symmetric.a: $(BIN_DIR)/generator
	@mkdir -p $(@D)
	$^ -g is_representative -o $(@D) -e $(GENERATOR_OUTPUTS) \
		-n is_representative_symmetric \
		-f ls_internal_is_representative_symmetric_kernel$(call clean_kernel_name,$*) \
		spin_inversion=1 auto_schedule=false \
		target=$*

$(BIN_DIR)/%/is_representative_antisymmetric.a: $(BIN_DIR)/generator
	@mkdir -p $(@D)
	$^ -g is_representative -o $(@D) -e $(GENERATOR_OUTPUTS) \
		-n is_representative_antisymmetric \
		-f ls_internal_is_representative_antisymmetric_kernel$(call clean_kernel_name,$*) \
		spin_inversion=-1 auto_schedule=false \
		target=$*

$(BIN_DIR)/%/is_representative_general.a: $(BIN_DIR)/generator
	@mkdir -p $(@D)
	$^ -g is_representative -o $(@D) -e $(GENERATOR_OUTPUTS) \
		-n is_representative_general \
		-f ls_internal_is_representative_general_kernel$(call clean_kernel_name,$*) \
		spin_inversion=0 auto_schedule=false \
		target=$*

# Build common runtime
$(BIN_DIR)/runtime.a: $(BIN_DIR)/generator
	@echo Building runtime ...
	@mkdir -p $(@D)
	$^ -r runtime -o $(@D) target=$(HL_BASE_TARGET)

# Build the Halide generator for our kernels
$(BIN_DIR)/generator: generator.cpp $(GENERATOR_DEPS_STATIC)
	@echo Building Generator ...
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(filter-out %.h,$^) -o $@ $(LIBHALIDE_LDFLAGS_STATIC)
