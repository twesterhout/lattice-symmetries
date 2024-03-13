{ lib
, stdenv
, fetchurl
, gfortran
, python3
, blas
, lapack
, mpi
, openssh
, petsc
, gnumake
}:

stdenv.mkDerivation rec {
  pname = "slepc";
  version = "3.20.1";

  src = fetchurl {
    url = "https://slepc.upv.es/download/distrib/slepc-${version}.tar.gz";
    sha256 = "sha256-Wja2ZIlYgdOFjQZE9Wv3u5Ir2rcNcy+hHL9kQv7BGAY=";
  };

  patches = [ ./filter_mpi_warnings.patch ];

  strictDeps = true;
  nativeBuildInputs = [ python3 gnumake gfortran ]
    ++ lib.optional petsc.mpiSupport mpi
    ++ lib.optional (petsc.mpiSupport && mpi.pname == "openmpi") openssh;
  buildInputs = [ blas lapack petsc ];

  preConfigure = ''
    export PETSC_DIR=${petsc} PETSC_ARCH=""
  '';
  configureScript = "python ./configure";

  preBuild = ''
    export SLEPC_DIR=$PWD
  '';

  checkPhase = "make check";

  enableParallelBuilding = true;
  doCheck = stdenv.hostPlatform == stdenv.buildPlatform;

  meta = with lib; {
    description = "Scalable Library for Eigenvalue Problem Computations";
    homepage = "https://slepc.upv.es";
    license = licenses.bsd2;
    maintainers = with maintainers; [ twesterhout ];
  };
}
