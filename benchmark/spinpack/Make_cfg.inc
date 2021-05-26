# Generated automatically by configure 2.57c on Thu Dec 20 11:53:33 CET 2018.
# found CC=g++
# tested successfully: MPICC=mpicxx -DMPICH_IGNORE_CXX_SEEK -DCONFIG_MPI=1
# use option -pg or -gprof for gprof-analyzes
# best performance: CFLAGS+= -march=native -mtune=native
CC=gcc
#MPICC=mpicc -DMPICH_IGNORE_CXX_SEEK -DCONFIG_MPI=16
#MPICC=g++
MPICC=gcc
MPIEXEC?=mpiexec -np 2
CFLAGS= -g -Wall -Wno-unused-function -Wno-unknown-pragmas
CFLAGS+= -O2 -funroll-loops -ftree-vectorize
CFLAGS+= -march=native -mtune=native
#CFLAGS+= -march=atom -mtune=atom
#CFLAGS+= -march=i586 -mtune=i586 -m32
DFLAGS= -D_REENTRANT -DCONFIG_PTHREAD=64
DFLAGS += -DCFG_CPUSET=2 -D_GNU_SOURCE
DFLAGS += -DNZXMAX=75
DFLAGS += -DNUM_AHEAD_LINES=2048
#DFLAGS+= -DNOSZSYM=1
#DFLAGS += -DHRMAX=0
#DFLAGS += -DCFG_DM=1
#DFLAGS+= -DCFG_FTLM=1
DFLAGS += -DVecType=1
#DFLAGS+= -DDBG_NHM=1
LIB= -g  -lm
LIB += -lpthread
LIB += -static
# tested by JueSch1905: openblas failed
#LIB += ${MKL_LIB}
