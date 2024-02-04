{ stdenv
, gcc
, cmake
, halide
, libffcall
, version
}:

stdenv.mkDerivation {
  pname = "lattice-symmetries-kernels2";
  inherit version;
  src = ./.;
  buildInputs = [
    halide
    libffcall
  ];
  nativeBuildInputs = [
    cmake
    gcc
  ];
}
