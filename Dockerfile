FROM quay.io/centos/centos:stream8 as base

# Set an encoding to make things work smoothly.
ENV LANG=en_US.UTF-8

# gcc gcc-c++
# clang-tools-extra-14.0.6
# python3-devel 
RUN dnf -y install \
    bash readline \
    m4 perl python3 make gawk git cmake \
    ncurses-compat-libs gmp gmp-devel \
    llvm-devel-14.0.6 \
    clang-14.0.6 clang-devel-14.0.6 clang-libs-14.0.6 \
    lld-14.0.6 lld-devel-14.0.6
RUN dnf clean all && rm -rf /var/cache/dnf/* /var/cache/yum/*

RUN mkdir /work

## Install Halide
################################################################################
RUN cd /work && \
    git clone --depth=1 --single-branch --branch v14.0.0 https://github.com/halide/Halide.git

RUN cd /work/Halide && \
    cmake -DCMAKE_BUILD_TYPE=Release \
    -DHalide_SHARED_LLVM=ON \
    -DWITH_PYTHON_BINDINGS=OFF \
    -DWITH_DOCS=OFF \
    -DWITH_TESTS=OFF \
    -DWITH_TUTORIALS=OFF \
    -DTARGET_HEXAGON=OFF \
    -DTARGET_POWERPC=OFF \
    -DTARGET_RISCV=OFF \
    -DTARGET_WEBASSEMBLY=OFF \
    -DTARGET_METAL=OFF \
    -DTARGET_D3D12COMPUTE=OFF \
    -S . -B build && \
    cmake --build build -- -j4 && \
    cmake --install build --prefix /opt/Halide && \
    cd / && \
    rm -rf /work/Halide

## Install Chapel
################################################################################
RUN cd /opt && \
    git clone --depth=1 --single-branch --branch 1.29.0 https://github.com/chapel-lang/chapel.git

ENV CHPL_HOME=/opt/chapel
ENV GASNET_BACKTRACE=1
ENV CHPL_COMM=none
ENV CHPL_LLVM=system
ENV CHPL_TARGET_CPU=none
ENV CHPL_TASKS=qthreads
ENV CHPL_MEM=jemalloc
ENV CHPL_ATOMICS=cstdlib
ENV CHPL_UNWIND=bundled
ENV CHPL_HWLOC=bundled
ENV CHPL_RE2=none
ENV CHPL_GMP=none
ENV CHPL_LIB_PIC=pic

RUN cd ${CHPL_HOME} && \
    make -j4

## Install Haskell
################################################################################
ENV GHCUP_INSTALL_BASE_PREFIX=/opt
ENV CABAL_DIR=/opt/cabal

RUN curl -s -L https://downloads.haskell.org/~ghcup/x86_64-linux-ghcup > /usr/bin/ghcup && \
    chmod +x /usr/bin/ghcup

ENV PATH /opt/.ghcup/bin:$PATH

RUN dnf install -y readline ncurses-compat-libs gmp gmp-devel

RUN ghcup install ghc 9.2.5 && \
    ghcup set ghc 9.2.5 && \
    ghcup install cabal 3.6.2.0 && \
    ghcup set cabal 3.6.2.0 && \
    cabal update

## Install patchelf
################################################################################
RUN cd /usr && \
    curl -s -L -O https://github.com/NixOS/patchelf/releases/download/0.17.2/patchelf-0.17.2-x86_64.tar.gz && \
    tar -x ./bin/patchelf -f patchelf-0.17.2-x86_64.tar.gz -C /usr && \
    rm patchelf-0.17.2-x86_64.tar.gz && \
    chmod +x /usr/bin/patchelf

FROM base as dependencies

RUN mkdir -p /work/lattice-symmetries-haskell
ADD cabal.project                     /work/lattice-symmetries-haskell/
ADD lattice-symmetries-haskell.cabal  /work/lattice-symmetries-haskell/

WORKDIR /work/lattice-symmetries-haskell

RUN cabal build --dependencies-only && \
    rm cabal.project lattice-symmetries-haskell.cabal

RUN dnf install -y zlib-devel
# Install basic requirements.
# RUN yum update -y dnf
# 
# RUN yum install -y bash centos-release-scl
# 
# RUN yum install -y devtoolset-11-gcc*
# 
# RUN yum install -y devtoolset-11-gcc*
#     scl-utils \
#     scl-utils-build \
#     llvm-toolset-11.0 \
#     epel-release

# RUN ln -sf /bin/bash /bin/sh
# 
# RUN scl --list
# 
# RUN find /etc/scl/conf && \
#     scl enable llvm-toolset-11.0 bash
# 
# RUN echo 'source scl_source enable devtoolset-11' >> ~/.bashrc && \
#     yum install -y \
#     bash \
#     cmake3 \
#     gawk \
#     git \
#     gcc \
#     gcc-c++ \
#     gmp \
#     gmp-devel \
#     make \
#     m4 \
#     ncurses \
#     ncurses-compat-libs \
#     perl \
#     python3 \
#     tar \
#     xz \
#     which && \
#     echo 'export CMAKE=cmake3' >> ~/.bashrc && \
#     yum clean all && \
#     rm -rf /var/cache/yum/*
# 
# RUN mkdir /work
# 
# WORKDIR /work
# 
# # Install Chapel
# RUN git clone --depth=1 --single-branch --branch 1.29.0 https://github.com/chapel-lang/chapel.git
# 
# WORKDIR /work/chapel
# 
# RUN export CHPL_HOME=/work/chapel && \
#     export CHPL_LLVM=system && \
#     source util/printchplenv --all
# 
# 
# ENV GHCUP_INSTALL_BASE_PREFIX /opt
# ENV CABAL_DIR /opt/cabal
# 
# RUN curl -s -L https://downloads.haskell.org/~ghcup/x86_64-linux-ghcup > /usr/bin/ghcup && \
#     chmod +x /usr/bin/ghcup
# 
# ENV PATH /opt/.ghcup/bin:$PATH
# 
# RUN ghcup install ghc 9.2.5 && \
#     ghcup set ghc 9.2.5 && \
#     ghcup install cabal 3.6.2.0 && \
#     ghcup set cabal 3.6.2.0
# 
# RUN curl -s -L -O https://github.com/NixOS/patchelf/releases/download/0.17.2/patchelf-0.17.2-x86_64.tar.gz && \
#     tar -x ./bin/patchelf -f patchelf-0.17.2-x86_64.tar.gz -C /usr && \
#     rm patchelf-0.17.2-x86_64.tar.gz && \
#     chmod +x /usr/bin/patchelf
# 
# RUN mkdir /work
# 
# WORKDIR /work
# 
# RUN git clone --depth=1 --single-branch --branch llvmorg-15.0.7 https://github.com/llvm/llvm-project.git
# 
# RUN cmake -G Ninja -S llvm-project -B llvm-build \
#     -DCMAKE_BUILD_TYPE=Release \
#     "-DCMAKE_OSX_ARCHITECTURES=arm64;x86_64" \
#     "-DLLVM_TARGETS_TO_BUILD=X86;ARM;AArch64" \
#     "-DLLVM_ENABLE_PROJECTS=clang;lld" \
#     -DLLVM_ENABLE_ASSERTIONS=ON \
#     -DLLVM_ENABLE_RTTI=ON \
#     -DLLVM_ENABLE_EH=ON \
#     -DLLVM_ENABLE_LIBXML2=OFF \
#     -DLLVM_ENABLE_TERMINFO=OFF \
#     -DLLVM_ENABLE_ZSTD=OFF \
#     -DLLVM_ENABLE_ZLIB=OFF \
#     -DLLVM_ENABLE_OCAMLDOC=OFF \
#     -DLLVM_ENABLE_BINDINGS=OFF \
#     -DLLVM_ENABLE_IDE=OFF
# 
# RUN cabal update
# 
# # FROM base as dependencies
# # 
# # RUN mkdir -p /work/lattice-symmetries-haskell
# # ADD cabal.project                     /work/lattice-symmetries-haskell/
# # ADD lattice-symmetries-haskell.cabal  /work/lattice-symmetries-haskell/
# # 
# # WORKDIR /work/lattice-symmetries-haskell
# # 
# # RUN cabal build --dependencies-only && \
# #     rm cabal.project lattice-symmetries-haskell.cabal
# 
CMD [ "/bin/bash", "-i" ]
