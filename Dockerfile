FROM twesterhout/haskell-chapel-halide:latest as dependencies

RUN mkdir -p /work/lattice-symmetries-haskell
ENV WORKDIR /work/lattice-symmetries-haskell
WORKDIR ${WORKDIR}

ADD cabal.project                     ${WORKDIR}/
ADD lattice-symmetries-haskell.cabal  ${WORKDIR}/

RUN cabal build --dependencies-only --disable-tests && \
    rm cabal.project lattice-symmetries-haskell.cabal

CMD [ "/bin/bash", "-i" ]
