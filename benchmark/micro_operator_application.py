import time
import sys
import os
import timeit
import subprocess
from loguru import logger
import numpy as np
import quspin
from quspin.basis import spin_basis_general

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "python"))
import lattice_symmetries as ls

from systems import *


def benchmark_rectangles(backend, block_size=1):
    assert backend in {"ls", "quspin"}
    name = "Lattice Symmetries" if backend == "ls" else "QuSpin"
    # r = dict()
    for L_y, L_x in [(1, 34)]:
        logger.info("Benchmarking {} for {}x{}...", name, L_y, L_x)
        symmetries = square_lattice_symmetries(L_x, L_y)
        nearest, _ = square_lattice_edges(L_x, L_y)
        if backend == "quspin" and (L_y, L_x) == (5, 5):
            # quspin incorrectly computes Ns_block_est
            extra_args = {"Ns_block_est": 26500}
        elif backend == "quspin" and (L_y, L_x) == (6, 6):
            # quspin incorrectly computes Ns_block_est
            extra_args = {"Ns_block_est": 15804956}
        else:
            extra_args = dict()

        basis = make_basis(
            symmetries,
            backend=backend,
            number_spins=L_x * L_y,
            hamming_weight=(L_x * L_y) // 2,
            **extra_args
        )
        logger.info("Constructed the basis.")
        hamiltonian = make_heisenberg(basis, nearest, backend=backend)

        n = basis.number_states if backend == "ls" else basis.Ns
        x = 0.5 - np.random.rand(n, block_size)
        x = np.asfortranarray(x)

        if backend == "ls":
            func = lambda: hamiltonian(x)
        else:
            assert hamiltonian.matmat(x).shape[1] == block_size
            if block_size == 1:
                func = lambda: hamiltonian.matvec(x)
            else:
                func = lambda: hamiltonian.matmat(x)

        ts = timeit.repeat(func, repeat=3, number=1)
        logger.info("  -> {:.2f} ± {:.2f}", np.mean(ts), np.std(ts))
        yield ("${} \\times {}$".format(L_y, L_x), (np.mean(ts), np.std(ts)))
    # return r


if __name__ == "__main__":
    dict(benchmark_rectangles(backend="ls", block_size=1))
