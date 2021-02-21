#!/bin/bash

set -ex

mkdir -p build.shared
pushd build.shared
rm -rf -- *
cmake -GNinja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DBUILD_SHARED_LIBS=ON \
  -DLatticeSymmetries_ENABLE_UNIT_TESTING=OFF \
  -DLatticeSymmetries_ENABLE_CLANG_TIDY=OFF \
  -DLatticeSymmetries_ENABLE_CPPCHECK=OFF \
  -DLatticeSymmetries_ENABLE_CODE_COVERAGE=OFF \
  -DLatticeSymmetries_LINK_STDLIB_STATICALLY=OFF \
  ..
cmake --build .
find . -name "*.so*" -maxdepth 1 -type f | while read sofile; do
  echo "Setting rpath of $sofile to" '$ORIGIN'
  patchelf --set-rpath '$ORIGIN' --force-rpath $sofile
  patchelf --print-rpath $sofile
done
cmake --build . --target install
popd

pushd python
$PYTHON setup.py install
