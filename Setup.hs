{-# OPTIONS_GHC -Wall #-}

module Main (main) where

import Data.Graph (topSort)
import Distribution.InstalledPackageInfo
  ( extraLibraries,
    hsLibraries,
    libraryDirs,
  )
import Distribution.Package
  ( PackageInstalled,
    getHSLibraryName,
  )
import Distribution.PackageDescription (PackageDescription ())
import Distribution.Simple
  ( UserHooks (..),
    defaultMainWithHooks,
    simpleUserHooks,
  )
import Distribution.Simple.BuildPaths (autogenPackageModulesDir)
import Distribution.Simple.LocalBuildInfo
import Distribution.Simple.PackageIndex
import Distribution.Simple.Setup
  ( ConfigFlags (configProfLib, configVerbosity),
    fromFlag,
    fromFlagOrDefault,
  )
import Distribution.Simple.Utils
  ( createDirectoryIfMissingVerbose,
    rewriteFileEx,
  )
import Distribution.Verbosity (Verbosity)
import System.FilePath

main :: IO ()
main =
  defaultMainWithHooks
    simpleUserHooks
      { postConf = \args flags pkg lbi -> do
          generateBuildModule (fromFlag (configVerbosity flags)) pkg lbi
          postConf simpleUserHooks args flags pkg lbi
      }

-- | Generate a part of a Makefile which contains all libraries and
-- include locations used by the Cabal library.
generateBuildModule ::
  Verbosity -> PackageDescription -> LocalBuildInfo -> IO ()
generateBuildModule verbosity pkgDesc lbi = do
  let autodir = autogenPackageModulesDir lbi
  createDirectoryIfMissingVerbose verbosity True autodir
  let installDirs = absoluteInstallDirs pkgDesc lbi NoCopyDest

  withLibLBI pkgDesc lbi $ \_ libLBI -> do
    let thisLib = getHSLibraryName (componentUnitId libLBI)
    let pkgs = orderedPackagesList (installedPkgs lbi)
        libdirs = libdir installDirs : concatMap libraryDirs pkgs
        libNames = thisLib : map threadedVersion (concatMap hsLibraries pkgs)
        mkLibName x
          | fromFlagOrDefault
              False
              (configProfLib (configFlags lbi)) =
            x ++ "_p"
          | otherwise = x

    rewriteFileEx verbosity (autodir </> "HS_LIBRARIES_LIST") $
      unlines $
        map mkLibName libNames

    rewriteFileEx verbosity (autodir </> "HS_LIBRARY_PATHS_LIST") $
      unlines libdirs

    rewriteFileEx verbosity (autodir </> "EXTRA_LIBRARIES_LIST") $
      unlines $
        extraLibraries =<< pkgs

orderedPackagesList :: PackageInstalled a => PackageIndex a -> [a]
orderedPackagesList pkgs = lookupVertex <$> topSort g
  where
    (g, lookupVertex, _findVertex) = dependencyGraph pkgs

-- We needed the threaded run-time so that SIGINT can be handled
-- cleanly when C code has called into Haskell
threadedVersion :: String -> String
threadedVersion lib =
  case lib of
    "Cffi" -> "Cffi_thr"
    "HSrts" -> "HSrts_thr"
    _ -> lib
