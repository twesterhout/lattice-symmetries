import sys
import os
import time
from loguru import logger
import numpy as np

try:
    import lattice_symmetries as ls
except ImportError:
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "python"))
    import lattice_symmetries as ls


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
    assert L_x > 0 and L_y > 0
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


def square_lattice_edges(L_x, L_y):
    assert L_x > 0 and L_y > 0
    # Example 4x6 square to illustrate the ordering of sites:
    #  [[ 0,  1,  2,  3,  4,  5],
    #   [ 6,  7,  8,  9, 10, 11],
    #   [12, 13, 14, 15, 16, 17],
    #   [18, 19, 20, 21, 22, 23]])
    sites = np.arange(L_y * L_x, dtype=np.int32).reshape(L_y, L_x)

    def generate_nearest_neighbours():
        for y in range(L_y):
            for x in range(L_x):
                if L_x > 1:
                    yield (sites[y, x], sites[y, (x + 1) % L_x])
                if L_y > 1:
                    yield (sites[y, x], sites[(y + 1) % L_y, x])

    def generate_next_nearest_neighbours():
        if L_x == 1 or L_y == 1:
            return
        for y in range(L_y):
            for x in range(L_x):
                yield (sites[y, x], sites[(y + 1) % L_y, (x + L_x - 1) % L_x])
                yield (sites[y, x], sites[(y + 1) % L_y, (x + 1) % L_x])

    return list(generate_nearest_neighbours()), list(generate_next_nearest_neighbours())


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
    spin_inversion = None
    processed_symmetries = []
    for s in symmetries:
        _, x1, x2 = s
        if x1 is None:
            assert x2 == 0 or x2 == 1
            spin_inversion = 1 if x2 == 0 else -1
        else:
            processed_symmetries.append(ls.Symmetry(x1, sector=x2))

    group = ls.Group(processed_symmetries)
    logger.info("Symmetry group contains {} elements", len(group))
    basis = ls.SpinBasis(
        group,
        number_spins=number_spins,
        hamming_weight=hamming_weight,
        spin_inversion=spin_inversion,
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


def _quspin_make_heisenberg(basis, nearest, next_nearest=None, j2=None, dtype=np.float64, matrix=False):
    from quspin.operators import quantum_LinearOperator, hamiltonian

    static = [
        ["+-", [[0.5, i, j] for (i, j) in nearest]],
        ["-+", [[0.5, i, j] for (i, j) in nearest]],
        ["zz", [[1.0, i, j] for (i, j) in nearest]],
    ]
    if next_nearest is not None:
        assert j2 is not None
        static += [
            ["+-", [[0.5 * j2, i, j] for (i, j) in next_nearest]],
            ["-+", [[0.5 * j2, i, j] for (i, j) in next_nearest]],
            ["zz", [[1.0 * j2, i, j] for (i, j) in next_nearest]],
        ]
    if matrix:
        return hamiltonian(static, [], basis=basis, dtype=dtype)
    return quantum_LinearOperator(static, basis=basis, dtype=dtype)


def _ls_make_heisenberg(basis, nearest, next_nearest=None, j2=None, dtype=None):
    matrix = [[1, 0, 0, 0], [0, -1, 2, 0], [0, 2, -1, 0], [0, 0, 0, 1]]
    interactions = [ls.Interaction(matrix, nearest)]
    if next_nearest is not None:
        assert j2 is not None
        interactions.append(ls.Interaction(j2 * matrix, next_nearest))
    return ls.Operator(basis, interactions)


def make_heisenberg(*args, backend="ls", **kwargs):
    if backend == "ls":
        return _ls_make_heisenberg(*args, **kwargs)
    elif backend == "quspin":
        return _quspin_make_heisenberg(*args, **kwargs)
    else:
        raise ValueError("invalid backend: {}; expected either 'ls' or 'quspin'".format(backend))
