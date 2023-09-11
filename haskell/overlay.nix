{ lib
, stdenv
, withPic
}:

self: super: rec {
  haskell = super.haskell // {
    packageOverrides = lib.composeExtensions super.haskell.packageOverrides
      (hself: hsuper: {
        lattice-symmetries-haskell =
          (hself.callCabal2nix "lattice-symmetries-haskell" ./. {
            inherit (super.lattice-symmetries) kernels;
          }).overrideAttrs
            (attrs: {
              outputs = (attrs.outputs or [ ]) ++ [ "lib" ];
              postInstall = ''
                ${attrs.postInstall or ""}

                echo "Installing foreign library to $lib ..."
                mkdir -p $lib/lib
                install -v -Dm 755 \
                  $out/lib/ghc-*/lib/liblattice_symmetries_haskell.so \
                  $lib/lib/

                mkdir -p $out/include
                install -v -Dm 644 \
                  lattice_symmetries_functions.h \
                  $out/include/
              '';
            });
      }
      # We only need to rebuild all packages with fPIC on Linux, Darwin does it by default (right?)
      // lib.optionalAttrs (stdenv.isLinux && withPic && (lib.versionAtLeast hsuper.ghc.version "9.6.1")) {
        # Ensure that all Haskell packages are built with -fPIC
        mkDerivation = args: (hsuper.mkDerivation args).overrideAttrs (attrs: {
          configureFlags = (attrs.configureFlags or [ ]) ++ [
            "--ghc-option=-fPIC"
            "--ghc-option=-fexternal-dynamic-refs"
          ];
        });
      }
      // lib.optionalAttrs (lib.versionAtLeast hsuper.ghc.version "9.6.1") {
        # Loosen constraints to make them build with 9.6.1
        tagged = super.haskell.lib.doJailbreak hsuper.tagged;
        zigzag = super.haskell.lib.doJailbreak hsuper.zigzag;
        typerep-map = super.haskell.lib.doJailbreak hsuper.typerep-map;
        relude = super.haskell.lib.dontCheck (super.haskell.lib.doJailbreak hsuper.relude);
        bytebuild = super.haskell.lib.doJailbreak hsuper.bytebuild;
        chronos = super.haskell.lib.doJailbreak hsuper.chronos;
      });
  };
  lattice-symmetries = (super.lattice-symmetries or { }) // {
    haskell = haskell.packages.ghc96.lattice-symmetries-haskell;
  };
}
