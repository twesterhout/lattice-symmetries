{ lib
, doInstallForeignLibs
, doEnableRelocatedStaticLibs
}:

lib.composeManyExtensions [
  (doEnableRelocatedStaticLibs "ghc962")
  (doEnableRelocatedStaticLibs "ghc963")
  (doEnableRelocatedStaticLibs "ghc964")

  (self: super: rec {
    haskell = super.haskell // {
      packageOverrides = lib.composeExtensions super.haskell.packageOverrides
        (hself: hsuper: {
          lattice-symmetries-haskell =
            doInstallForeignLibs
              { headers = [ "lattice_symmetries_functions.h" ]; }
              (hself.callCabal2nix "lattice-symmetries-haskell" ./. { inherit (super.lattice-symmetries) kernels_v2; });
        });
    };
    lattice-symmetries = (super.lattice-symmetries or { }) // {
      haskell = haskell.packages.ghc96.lattice-symmetries-haskell.lib;
    };
  })
]
