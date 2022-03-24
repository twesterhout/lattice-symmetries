MAKEFLAGS += --no-builtin-rules
.POSIX:
.SUFFIXES:

CC     = cc
CFLAGS = -fPIC -O3 -march=native -W -Wall -Wextra -pedantic
UNAME = $(shell uname)
BIN_DIR = build


ifeq ($(UNAME), Darwin)
SHARED_FLAG = -dynamiclib
else
SHARED_FLAG = -shared
endif

ifeq ($(OS), Windows_NT)
    $(error Windows is not supported)
else
ifeq ($(UNAME), Darwin)
    SHARED_EXT = dylib
else
    SHARED_EXT = so
endif
endif

all: $(BIN_DIR)/liblattice_symmetries_core.$(SHARED_EXT)

$(BIN_DIR)/reference.o: reference.c ../cbits/lattice_symmetries_haskell.h
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -I../cbits `pkg-config --cflags lattice_symmetries` -o $@ -c $<

$(BIN_DIR)/liblattice_symmetries_core.$(SHARED_EXT): $(BIN_DIR)/reference.o
	@mkdir -p $(@D)
	$(CC) $(SHARED_FLAG) -o $(BIN_DIR)/liblattice_symmetries_core.$(SHARED_EXT) $^
ifeq ($(UNAME), Darwin)
	install_name_tool -id $(CURDIR)/$(BIN_DIR)/liblattice_symmetries_core.$(SHARED_EXT) $(BIN_DIR)/liblattice_symmetries_core.$(SHARED_EXT)
endif

clean:
	rm -fr $(BIN_DIR)
