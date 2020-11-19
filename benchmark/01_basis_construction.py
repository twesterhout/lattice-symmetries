import time
import sys
import os
import timeit
from loguru import logger
import numpy as np
import quspin
from quspin.basis import spin_basis_general

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "python"))
import lattice_symmetries as ls

from systems import make_basis, chain_symmetries


def benchmark_chains(backend):
    assert backend in {"ls", "quspin"}
    name = "Lattice Symmetries" if backend == "ls" else "QuSpin"
    r = dict()
    for L in range(20, 26, 2):
        logger.info("Benchmarking {} for L={}...", name, L)
        symmetries = chain_symmetries(L)
        func = lambda: make_basis(
            symmetries, backend=backend, number_spins=L, hamming_weight=L // 2
        )
        ts = timeit.repeat(func, repeat=3, number=1)
        r["$1 \\times {}$".format(L)] = (np.mean(ts), np.std(ts))
    return r

def benchmark_squares(backend):
    assert backend in {"ls", "quspin"}
    pass


def main():
    output_file = "01_basis_construction.dat"
    with open(output_file, "w") as output:
        output.write("system\trelative\terror\n")
    rs = benchmark_chains(backend="ls")
    for (key, (mean, _)) in benchmark_chains(backend="quspin").items():
        with open(output_file, "a") as output:
            output.write("{}\t{}\t{}\n".format(key, mean / rs[key][0], 0.0))


if __name__ == "__main__":
    main()
