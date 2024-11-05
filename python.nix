{ version
}:

final: prev: {
  pythonPackagesExtensions = prev.pythonPackagesExtensions ++ [
    (python-final: python-prev: {

      # grip = python-prev.grip.overrideAttrs (attrs: {
      #   src = final.fetchFromGitHub {
      #     owner = "Antonio-R1";
      #     repo = "grip";
      #     rev = "d2efd3c6a896c01cfd7624b6504107e7b3b4b20f";
      #     hash = "sha256-0wgIM7Ll5WELvAOiu1TLyoNSrhJ22Y1SRbWqa3BDF3k=";
      #   };
      #   checkPhase = "true";
      #   installCheckPhase = "true";
      # });

      parallel-sparse-tools = python-final.buildPythonPackage {
        pname = "parallel-sparse-tools";
        version = "0.2.3";
        src = final.fetchFromGitHub {
          owner = "QuSpin";
          repo = "parallel-sparse-tools";
          rev = "ef72b076a0f50bab56afdec709f284ed63aa4dbf";
          hash = "sha256-/Lf6vVrXB80tLOKFdMxGMdsNj3w3HScW4tmgGLIBoQk=";
        };
        pyproject = true;
        propagatedBuildInputs = with python-final; [ numpy scipy ];
        pythonRelaxDeps = [ "numpy" ];
        buildInputs = [ ];
        nativeBuildInputs = with python-final; [ setuptools cython ];
        nativeCheckInputs = with python-final; [ pip pytestCheckHook ];
      };

      quspin-extensions = python-final.buildPythonPackage {
        pname = "quspin-extensions";
        version = "0.1.6";
        src = final.fetchFromGitHub {
          owner = "QuSpin";
          repo = "quspin-extensions";
          rev = "c27603d7c6c3ffe39c8cf31e1a053fe22d35a313";
          hash = "sha256-cqFWz/4kB+X/QS/lrEg8MHWpCJBCwzY8YX2rwqXgI0A=";
        };
        pyproject = true;

        postPatch = ''
          export BOOST_ROOT=${final.boost}
        '';
        propagatedBuildInputs = with python-final; [ numpy scipy gmpy2 ];
        pythonRelaxDeps = [ "gmpy2" ];
        buildInputs = [ final.boost ];
        nativeBuildInputs = with python-final; [ setuptools cython ];
      };

      quspin = python-final.buildPythonPackage rec {
        pname = "quspin";
        version = "1.0.0";
        src = final.fetchFromGitHub {
          owner = "QuSpin";
          repo = "QuSpin";
          rev = "v${version}";
          hash = "sha256-Etu45rhkeLMsWObasJoFL0AFLQ+o7fC4Lg1msSntLKA=";
        };
        pyproject = true;

        propagatedBuildInputs = with python-final; [
          quspin-extensions
          parallel-sparse-tools
          numpy
          dill
          scipy
          matplotlib
          numexpr
          numba
          six
          joblib
        ];
        pythonRelaxDeps = [ "numpy" "numexpr" ];
        nativeBuildInputs = with python-final; [ pdm-backend ];
        nativeCheckInputs = with python-final; [ pytestCheckHook ];
      };

      petsc4py =
        assert final.petsc.version == "3.21.3";
        python-final.buildPythonPackage rec {
          pname = "petsc4py";
          version = final.petsc.version;
          src = final.fetchurl {
            url = "https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc4py-${version}.tar.gz";
            hash = "sha256-HDZk1bUnNUFxB3yJxLH+899KQb5xltErynSydZx+Jkg=";
          };
          preConfigure = ''
            export PETSC_DIR=${final.petsc} PETSC_ARCH=""
            rm conf/epydoc*
          '';
          strictDeps = true;
          propagatedBuildInputs = with python-final; [ numpy ];
          buildInputs = with final; [ petsc ];
          nativeBuildInputs = with final; with python-final;
            [ cython ]
              ++ lib.optional petsc.mpiSupport mpi
              ++ lib.optional (petsc.mpiSupport && mpi.pname == "openmpi") openssh;
          nativeCheckInputs = with python-final; [ pytestCheckHook ];
          pytestFlagsArray = [
            "test/"
          ];
          disabledTestPaths = [
            "test/test_stdout.py"
          ];
        };

      slepc4py =
        assert final.slepc.version == "3.21.2";
        python-final.buildPythonPackage rec {
          pname = "slepc4py";
          version = final.slepc.version;
          src = python-final.fetchPypi {
            inherit pname version;
            hash = "sha256-9hH/dOR0nyFEWyNp29Dt9ATN9jnuyv1UGH0KKGXVIaA=";
          };
          preConfigure = ''
            export PETSC_DIR=${final.petsc} PETSC_ARCH="" SLEPC_DIR=${final.slepc}
            rm conf/epydoc*
          '';
          strictDeps = true;
          propagatedBuildInputs = with python-final; [ numpy petsc4py ];
          buildInputs = with final; [ petsc slepc ];
          nativeBuildInputs = with final; with python-final;
            [ cython ]
              ++ lib.optional petsc.mpiSupport mpi
              ++ lib.optional (petsc.mpiSupport && mpi.pname == "openmpi") openssh;
          nativeCheckInputs = with python-final; [ pytestCheckHook ];
          pytestFlagsArray = [
            "test/"
          ];
        };

      dynamite = python-final.buildPythonPackage rec {
        pname = "dynamite";
        version = "0.4.0";
        src = final.fetchFromGitHub {
          owner = "GregDMeyer";
          repo = "dynamite";
          rev = "v${version}";
          hash = "sha256-bFH0H/Asc7yokA/jqqDWL/UYgHU4HKSFOnIkv60X7qE=";
        };
        pyproject = true;
        postPatch = ''
          substituteInPlace setup.py \
            --replace-fail "['git', 'describe', '--always']" "['echo', '${version}']" \
            --replace-fail "['git', 'rev-parse', '--abbrev-ref', 'HEAD']" "['echo', 'master']"
        '';
        preConfigure = ''
          export PETSC_DIR=${final.petsc} PETSC_ARCH="" SLEPC_DIR=${final.slepc}
        '';
        pythonRelaxDeps = [ "numpy" "petsc4py" "slepc4py" ];
        propagatedBuildInputs = with final; with python-final;
          [ numpy petsc4py scipy slepc4py threadpoolctl ]
            ++ lib.optional petsc.mpiSupport mpi4py;
        nativeBuildInputs = with final; with python-final;
          [ cython setuptools ]
            ++ lib.optional petsc.mpiSupport mpi
            ++ lib.optional (petsc.mpiSupport && mpi.pname == "openmpi") openssh;
        nativeCheckInputs = with python-final; [ pytestCheckHook ];
        pytestFlagsArray = [
          "tests/unit"
        ];
      };

      lattice-symmetries = python-final.buildPythonPackage rec {
        pname = "lattice-symmetries";
        inherit version;
        src = ./.;
        pyproject = true;

        # buildInputs = with final; [
        #   lattice-symmetries.kernels_v2
        #   lattice-symmetries.haskell
        #   lattice-symmetries.chapel
        # ];
        propagatedBuildInputs = with python-final; [
          cffi
          loguru
          numpy
          scipy
          sympy
          igraph
          halide
          more-itertools
          lark
        ] ++ lib.optionals (lib.versionAtLeast python-final.python.version "3.12") [
          quspin
          dynamite
        ];

        buildInputs = [
          final.lattice-symmetries-chapel
        ];

        nativeBuildInputs = with python-final; [
          setuptools
          final.tree
          final.ocl-icd
        ];

        nativeCheckInputs = with python-final; [
          pip
          pytestCheckHook
          pythonOutputDistHook
          hypothesis
          scalene
          # igraph
        ];

        postPatch = ''
          cp -v ${final.lattice-symmetries-chapel}/lib/liblattice_symmetries_chapel.* lattice_symmetries/
        '';

        preInstall = ''
          pushd dist/
          WHEEL_FILE=$(ls *.whl)
          wheel unpack $WHEEL_FILE
          rm -v $WHEEL_FILE

          pushd lattice_symmetries-${version}
          tree
          patchelf --debug --remove-rpath lattice_symmetries/liblattice_symmetries_chapel.*
          patchelf --debug --set-rpath '$ORIGIN' lattice_symmetries/_ls.*
          popd

          wheel pack lattice_symmetries-${version}
          rm -r lattice_symmetries-${version}
          popd
        '';

        # NEW_RPATH=$(patchelf --print-rpath $out/${python-final.python.sitePackages}/lattice_symmetries/_ls.* | sed -E 's;${final.lattice-symmetries-chapel}/lib;$ORIGIN;')
        # patchelf --debug --set-rpath "$NEW_RPATH" $out/${python-final.python.sitePackages}/lattice_symmetries/_ls.*
        #   # - Let CFFI see chunks between `python-cffi: START` and `python-cffi: STOP`
        #   # - Hide `LS_HS_ATOMIC` from CFFI
        #   awk '/python-cffi: START/{flag=1;next}/python-cffi: STOP/{flag=0}flag' \
        #     ${final.lattice-symmetries.kernels_v2}/include/lattice_symmetries_types.h \
        #     | sed -E 's/LS_HS_ATOMIC\(([^)]+)\)/\1/' \
        #     >lattice_symmetries/extracted_declarations.h
        #   awk '/python-cffi: START/{flag=1;next}/python-cffi: STOP/{flag=0}flag' \
        #     ${final.lattice-symmetries.haskell}/include/lattice_symmetries_functions.h \
        #     >>lattice_symmetries/extracted_declarations.h
        #   awk '/python-cffi: START/{flag=1;next}/python-cffi: STOP/{flag=0}flag' \
        #     ${final.lattice-symmetries.chapel}/include/lattice_symmetries_chapel.h \
        #     >>lattice_symmetries/extracted_declarations.h
        # '';

        preCheck = "rm -rf lattice_symmetries";

        checkPhase = ''
          runHook preCheck
          python3 -m pytest --color=yes --capture=no test/test_api.py | tee output.txt
          grep -q -E '(FAILURES|failed)' output.txt && exit 1
          runHook postCheck
        '';

        doCheck = false;

        shellHook = ''
          if test -e setup.py; then
            tmp_path="$PWD/.pip-install"
            rm -rf $tmp_path build/ lattice_symmetries/*.so
            ${postPatch}

            mkdir -p "$tmp_path"
            export PYTHONPATH="$tmp_path/${python-final.python.sitePackages}:$PYTHONPATH"
            python -m pip install -e . --prefix $tmp_path --no-deps # --no-build-isolation --config-settings editable_mode=compat
            export NIX_PYTHONPATH="$tmp_path/${python-final.python.sitePackages}:$${NIX_PYTHONPATH-}"
          fi
          export OCL_ICD_PATH=${final.ocl-icd}
        '';
      };

    } // (
      let
        disableTests = drv: drv.overridePythonAttrs (attrs: { doCheck = false; });
      in
      final.lib.optionalAttrs (python-prev.python.pythonOlder "3.11") {
        scipy = disableTests python-prev.scipy;
        django = disableTests python-prev.django;
        tifffile = disableTests python-prev.tifffile;
        pytest-django = disableTests python-prev.pytest-django;
        zarr = disableTests python-prev.zarr;
        geoip2 = disableTests python-prev.geoip2;
        aiohttp = disableTests python-prev.aiohttp;
        pillow-heif = disableTests python-prev.pillow-heif;
        astropy = disableTests python-prev.astropy;
        fsspec = disableTests python-prev.fsspec;
        dask = disableTests python-prev.dask;
        scikit-image = disableTests python-prev.scikit-image;
        patsy = disableTests python-prev.patsy;
        statsmodels = disableTests python-prev.statsmodels;
        imageio = python-prev.imageio.overridePythonAttrs (attrs: {
          optional-dependencies = {
            bsdf = [ ];
            dicom = [ ];
            feisem = [ ];
            ffmpeg = [ ];
            fits = [ ];
            freeimage = [ ];
            lytro = [ ];
            numpy = [ ];
            pillow = [ ];
            simpleitk = [ ];
            spe = [ ];
            swf = [ ];
            tifffile = [ ];
            pyav = [ ];
            heif = [ ];
          };
          doCheck = false;
          nativeCheckInputs = [
            python-prev.psutil
            python-prev.pytestCheckHook
          ];
        });
      }
    ))
  ];
  # lattice-symmetries = (prev.lattice-symmetries or { }) // {
  #   python = final.python3Packages.lattice-symmetries;
  #   apptainer-python-minimal = final.singularity-tools.buildImage {
  #     name = "lattice-symmetries-apptainer";
  #     contents = with final; [
  #       (python3.withPackages (ps: with ps; [
  #         lattice-symmetries
  #         loguru
  #         scipy
  #         numpy
  #         sympy
  #       ]))
  #       coreutils
  #       less # for displaying docs in the Python interpreter
  #     ];
  #     diskSize = 10240;
  #     memSize = 5120;
  #   };
  # };
}
