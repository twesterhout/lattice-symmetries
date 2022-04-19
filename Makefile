.POSIX:
.SUFFIXES:

CC = cc
CFLAGS = -fPIC
HDF5_MAJOR = 1
HDF5_MINOR = 12
HDF5_PATCH = 1
HDF5_VERSION = $(HDF5_MAJOR).$(HDF5_MINOR).$(HDF5_PATCH)
PREFIX = $(PWD)/third_party/hdf5

UNAME = $(shell uname)
ifeq ($(UNAME), Darwin)
  NPROC = $(shell sysctl -n hw.ncpu)
else
  NPROC = $(shell nproc --all)
endif

all:
hdf5: $(PREFIX)

third_party/hdf5-$(HDF5_VERSION).tar.bz2:
	@mkdir -p $(@D)
	cd third_party && \
	wget -q https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-$(HDF5_MAJOR).$(HDF5_MINOR)/hdf5-$(HDF5_VERSION)/src/hdf5-$(HDF5_VERSION).tar.bz2

third_party/hdf5-$(HDF5_VERSION)-src: third_party/hdf5-$(HDF5_VERSION).tar.bz2
	cd third_party && \
	tar -xf hdf5-$(HDF5_VERSION).tar.bz2 && \
	mv hdf5-$(HDF5_VERSION) hdf5-$(HDF5_VERSION)-src

$(PREFIX): third_party/hdf5-1.12.1-src
	cd third_party/hdf5-1.12.1-src && \
	CC=$(CC) CFLAGS="$(CFLAGS)" ./configure --prefix="$(PREFIX)" \
	    --disable-java --disable-fortran --disable-cxx \
	    --disable-tests --disable-tools \
	    --disable-shared && \
	make -j $(NPROC) install

.PHONY: clean
clean:
	rm -rf $(PREFIX) \
	       third_party/hdf5-$(HDF5_VERSION)-src
