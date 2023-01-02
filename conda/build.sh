#!/bin/bash

set -ex

mkdir -p build.shared
pushd build.shared
rm -rf -- *
cmake -GNinja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DBUILD_SHARED_LIBS=ON \
  -DLatticeSymmetries_ENABLE_UNIT_TESTING=OFF \
  -DLatticeSymmetries_ENABLE_CLANG_TIDY=OFF \
  -DLatticeSymmetries_ENABLE_CPPCHECK=OFF \
  -DLatticeSymmetries_ENABLE_CODE_COVERAGE=OFF \
  -DLatticeSymmetries_LINK_STDLIB_STATICALLY=OFF \
  -DCMAKE_OSX_DEPLOYMENT_TARGET=11.3 \
  ..
cmake --build .
if [[ $(uname) == "Linux" ]]; then
  find . -name "*.so*" -maxdepth 1 -type f | while read -r sofile; do
    # shellcheck disable=SC2016
    echo "Setting rpath of $sofile to" '$ORIGIN'
    # shellcheck disable=SC2016
    patchelf --set-rpath '$ORIGIN' --force-rpath "$sofile"
    patchelf --print-rpath "$sofile"
  done
fi
cmake --build . --target install
popd

pushd python
$PYTHON setup.py install
