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
      haskell =
        let original = self.haskell.packages.ghc96.lattice-symmetries-haskell;
        in original.lib;
      # self.stdenv.mkDerivation {
      #   pname = original.pname;
      #   version = original.version;
      #   dontUnpack = true;
      #   dontConfigure = true;
      #   dontBuild = true;
      #   dontCheck = true;
      #   installPhase = ''
      #     mkdir -p $out/lib
      #     for f in $(find ${original}/lib/ghc-*/lib -maxdepth 1 -type f -regex '.*\.\(so\|dylib\)'); do
      #       install -v -Dm 755 "$f" $out/lib/
      #     done

      #     mkdir -p $out/include
      #     install -v -Dm 644 "${original}/include/lattice_symmetries_functions.h" $out/include/
      #   '';
      # };
    };
  })
]
