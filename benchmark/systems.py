import sys
import os
import time
from loguru import logger
import numpy as np

sys.path.insert(0, "/home/tom/src/lattice-symmetries/python")


def get_processor_name():
    import subprocess
    import json

    result = subprocess.run(["lscpu", "-J"], check=False, capture_output=True)
    if result.returncode != 0:
        logger.warn(
            "Failed to get processor name: {} returned error code {}: {}",
            result.args,
            result.returncode,
            result.stderr,
        )
        return None
    for obj in json.loads(result.stdout)["lscpu"]:
        if obj["field"].startswith("Model name"):
            return obj["data"]


def square_lattice_symmetries(L_x, L_y, sectors=dict()):
    sites = np.arange(L_y * L_x, dtype=np.int32)
    x = sites % L_x
    y = sites // L_x

    symmetries = []
    if L_x > 1:
        T_x = (x + 1) % L_x + L_x * y  # translation along x-direction
        symmetries.append(("T_x", T_x, sectors.get("T_x", 0)))
        P_x = (L_x - 1 - x) + L_x * y  # reflection over y-axis
        symmetries.append(("P_x", P_x, sectors.get("P_x", 0)))
    if L_y > 1:
        T_y = x + L_x * ((y + 1) % L_y)  # translation along y-direction
        symmetries.append(("T_y", T_y, sectors.get("T_y", 0)))
        P_y = x + L_x * (L_y - 1 - y)  # reflection around x-axis
        symmetries.append(("P_y", P_y, sectors.get("P_y", 0)))
    if L_x == L_y and L_x > 1:  # Rotations are valid only for square samples
        R = np.rot90(sites.reshape(L_y, L_x), k=-1).reshape(-1)
        symmetries.append(("R", R, sectors.get("R", 0)))
    if L_x * L_y % 2 == 0:
        symmetries.append(("I", None, sectors.get("I", 0)))
    return symmetries


def chain_edges(n):
    return [(i, (i + 1) % n) for i in range(n)]


def _quspin_make_basis(symmetries, number_spins, hamming_weight=None, build=True, **kwargs):
    import quspin.basis

    def transform(t):
        x0, x1, x2 = t
        if x1 is None:
            return (x0, (-(np.arange(number_spins, dtype=np.int32) + 1), x2))
        return x0, (x1, x2)

    basis = quspin.basis.spin_basis_general(
        N=number_spins,
        Nup=hamming_weight,
        make_basis=build,
        **kwargs,
        **dict(transform(s) for s in symmetries)
    )
    return basis


def _ls_make_basis(symmetries, number_spins, hamming_weight=None, build=True):
    import lattice_symmetries

    def transform(t):
        x0, x1, x2 = t
        if x1 is None:
            return lattice_symmetries.Symmetry(
                np.arange(number_spins, dtype=np.int32), flip=True, sector=x2
            )
        return lattice_symmetries.Symmetry(x1, sector=x2)

    basis = lattice_symmetries.SpinBasis(
        lattice_symmetries.Group([transform(s) for s in symmetries]),
        number_spins=number_spins,
        hamming_weight=hamming_weight,
    )
    if build:
        basis.build()
    return basis


def make_basis(*args, backend="ls", **kwargs):
    if backend == "ls":
        return _ls_make_basis(*args, **kwargs)
    elif backend == "quspin":
        return _quspin_make_basis(*args, **kwargs)
    else:
        raise ValueError("invalid backend: {}; expected either 'ls' or 'quspin'".format(backend))
