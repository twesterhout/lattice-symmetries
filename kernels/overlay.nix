{ version
}:

final: prev: {
  lattice-symmetries = (prev.lattice-symmetries or { }) // {
    # kernels = final.callPackage ./. { inherit version; };
    kernels_v2 = final.callPackage ./v2 { inherit version; };
  };
}
