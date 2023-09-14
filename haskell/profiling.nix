{ lib
}:

self: super: {
  haskell = super.haskell // {
    packageOverrides = lib.composeExtensions super.haskell.packageOverrides
      (hself: hsuper: {
        mkDerivation = builtins.trace "enabling profiling" (args: hsuper.mkDerivation (args // { enableLibraryProfiling = true; }));
      });
  };
}
