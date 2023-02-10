#!/bin/bash

mkdir -p bundle-docker
mkdir -p dist-newstyle-docker
make clean
sudo docker run -it \
  --rm \
  -v $PWD/src:/work/lattice-symmetries-haskell/src:ro \
  -v $PWD/lib:/work/lattice-symmetries-haskell/lib:ro \
  -v $PWD/test:/work/lattice-symmetries-haskell/test:ro \
  -v $PWD/cbits:/work/lattice-symmetries-haskell/cbits:ro \
  -v $PWD/kernels:/work/lattice-symmetries-haskell/kernels:ro \
  -v $PWD/LICENSE:/work/lattice-symmetries-haskell/LICENSE:ro \
  -v $PWD/README.md:/work/lattice-symmetries-haskell/README.md:ro \
  -v $PWD/cabal.project:/work/lattice-symmetries-haskell/cabal.project:ro \
  -v $PWD/lattice-symmetries-haskell.cabal:/work/lattice-symmetries-haskell/lattice-symmetries-haskell.cabal:ro \
  -v $PWD/Makefile:/work/lattice-symmetries-haskell/Makefile:ro \
  -v $PWD/bundle-docker:/work/lattice-symmetries-haskell/bundle \
  -v $PWD/dist-newstyle-docker:/work/lattice-symmetries-haskell/dist-newstyle \
  twesterhout/lattice-symmetries-haskell \
  bash -c 'mkdir build && make HALIDE_PATH=/opt/Halide BIN_DIR=$PWD/build bindist'
sudo chown -R $USER:$USER bundle-docker
sudo chown -R $USER:$USER dist-newstyle-docker
# tar -czf lattice-symmetries-haskell.tar.gz lattice-symmetries-haskell
# rm -r lattice-symmetries-haskell
