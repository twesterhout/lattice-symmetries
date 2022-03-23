MAKEFLAGS += --no-builtin-rules
.POSIX:
.SUFFIXES:

CC     = cc
CFLAGS = -fPIC -O3 -march=native -W -Wall -Wextra -pedantic
UNAME = $(shell uname)
BIN_DIR = build

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
	$(CXX) -shared -o $(BIN_DIR)/liblattice_symmetries_core.$(SHARED_EXT) $^

clean:
	rm -fr $(BIN_DIR)
