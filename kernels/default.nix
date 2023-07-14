{ stdenv
, halide
, version
}:

stdenv.mkDerivation {
  pname = "lattice-symmetries-kernels";
  inherit version;
  src = ./.;
  dontConfigure = true;
  makeFlags = [ "HALIDE_PATH=${halide}" "PREFIX=$(out)" ];
}
