.POSIX:
.SUFFIXES:
CC     = cc
CFLAGS = -fPIC -O3 -march=native -W -Wall -Wextra -pedantic
LDLIBS =
# -lpthread -lm

$(shell mkdir -p build)

all: build/liblattice_symmetries_core.so

build/reference.o: reference.c ../cbits/lattice_symmetries_haskell.h
	$(CC) $(CFLAGS) -I../cbits `pkg-config --cflags lattice_symmetries` -o $@ -c $<

build/liblattice_symmetries_core.so: build/reference.o
	$(CC) -shared -o $@ $^

clean:
	rm -fr build/*
