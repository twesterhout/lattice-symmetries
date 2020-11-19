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

from systems import get_processor_name, make_basis, square_lattice_symmetries


def benchmark_rectangles(backend):
    assert backend in {"ls", "quspin"}
    name = "Lattice Symmetries" if backend == "ls" else "QuSpin"
    r = dict()
    for L_y, L_x in [(1, 24), (1, 28), (4, 6), (5, 5), (5, 6)]:
        logger.info("Benchmarking {} for {}x{}...", name, L_y, L_x)
        symmetries = square_lattice_symmetries(L_x, L_y)
        if backend == "quspin" and (L_y, L_x) == (5, 5):
            # quspin incorrectly computes Ns_block_est
            extra_args = {"Ns_block_est": 26500}
        else:
            extra_args = dict()

        func = lambda: make_basis(
            symmetries, backend=backend, number_spins=L_x * L_y,
            hamming_weight=(L_x * L_y) // 2,
            **extra_args
        )
        ts = timeit.repeat(func, repeat=3, number=1)
        r["${} \\times {}$".format(L_y, L_x)] = (np.mean(ts), np.std(ts))
    return r


def main():
    output_file = "01_basis_construction.dat"
    with open(output_file, "w") as output:
        output.write("# Date: {}\n".format(time.asctime()))
        cpu = get_processor_name()
        if cpu is not None:
            output.write("#  CPU: {}\n".format(cpu))
        output.write("system\trelative\terror\n")
    rs = benchmark_rectangles(backend="ls")
    for (key, (mean, _)) in benchmark_rectangles(backend="quspin").items():
        with open(output_file, "a") as output:
            output.write("{}\t{}\t{}\n".format(key, mean / rs[key][0], 0.0))


if __name__ == "__main__":
    main()
