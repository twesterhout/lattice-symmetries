{ version
}:

final: prev: {
  pythonPackagesExtensions = prev.pythonPackagesExtensions ++ [
    (python-final: python-prev: {
      lattice-symmetries = python-final.buildPythonPackage rec {
        pname = "lattice-symmetries";
        inherit version;
        src = ./.;
        pyproject = true;
        dependencies = with python-final; [ loguru numpy scipy sympy igraph halide more-itertools lark quspin final.simde cffi ];
        nativeBuildInputs = with python-final; [ final.tree final.ocl-icd setuptools ipython ];
        nativeCheckInputs = with python-final; [ pip pytestCheckHook pythonOutputDistHook hypothesis jax jaxlib jax-cuda12-plugin ];
        # preInstall = ''
        #   pushd dist/
        #   WHEEL_FILE=$(ls *.whl)
        #   wheel unpack $WHEEL_FILE
        #   rm -v $WHEEL_FILE

        #   pushd lattice_symmetries-${version}
        #   tree
        #   patchelf --debug --remove-rpath lattice_symmetries/liblattice_symmetries_chapel.*
        #   patchelf --debug --set-rpath '$ORIGIN' lattice_symmetries/_ls.*
        #   popd

        #   wheel pack lattice_symmetries-${version}
        #   rm -r lattice_symmetries-${version}
        #   popd
        # '';
        preCheck = "rm -rf lattice_symmetries";
        checkPhase = ''
          runHook preCheck
          python3 -m pytest --color=yes --capture=no test/test_api.py | tee output.txt
          grep -q -E '(FAILURES|failed)' output.txt && exit 1
          runHook postCheck
        '';
        doCheck = true;

        shellHook = ''
          if test -e setup.py; then
            tmp_path="$PWD/.pip-install"
            rm -rf $tmp_path build/ lattice_symmetries/*.so

            mkdir -p "$tmp_path"
            export PYTHONPATH="$tmp_path/${python-final.python.sitePackages}:$PYTHONPATH"
            python -m pip install -e . --prefix $tmp_path --no-deps # --no-build-isolation --config-settings editable_mode=compat
            export NIX_PYTHONPATH="$tmp_path/${python-final.python.sitePackages}:$${NIX_PYTHONPATH-}"
          fi
          export OCL_ICD_PATH=${final.ocl-icd}
        '';
      };
    })
  ];
}
