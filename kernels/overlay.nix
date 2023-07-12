{ version
}:

final: prev: {
  lattice-symmetries = (prev.lattice-symmetries or { }) // {
    kernels = final.callPackage ./. { inherit version; };
  };
}
