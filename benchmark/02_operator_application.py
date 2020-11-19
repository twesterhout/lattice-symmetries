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
    r = dict()
    # for L_y, L_x in [(1, 24), (1, 28), (4, 6), (5, 5), (5, 6)]:
    for L_y, L_x in [(1, 28)]:
        logger.info("Benchmarking {} for {}x{}...", name, L_y, L_x)
        symmetries = square_lattice_symmetries(L_x, L_y)
        nearest, _ = square_lattice_edges(L_x, L_y)
        if backend == "quspin" and (L_y, L_x) == (5, 5):
            # quspin incorrectly computes Ns_block_est
            extra_args = {"Ns_block_est": 26500}
        else:
            extra_args = dict()

        basis = make_basis(
            symmetries,
            backend=backend,
            number_spins=L_x * L_y,
            hamming_weight=(L_x * L_y) // 2,
            **extra_args
        )
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

        # if backend == "ls":
        #     E, V = ls.diagonalize(hamiltonian)
        #     print(E)
        # else:
        #     E, V = hamiltonian.eigsh(k=1, which="SA")
        #     print(E)

        ts = timeit.repeat(func, repeat=1, number=1)
        r["${} \\times {}$".format(L_y, L_x)] = (np.mean(ts), np.std(ts))
    return r


def main():
    block_size = 1
    output_file = "_02_operator_application_{}.dat".format(block_size)
    with open(output_file, "w") as output:
        output.write("# Date: {}\n".format(time.asctime()))
        cpu = get_processor_name()
        if cpu is not None:
            output.write("#  CPU: {}\n".format(cpu))
        output.write("system\tls\tquspin\trelative\terror\n")
    rs = benchmark_rectangles("ls", block_size)
    for (key, (mean, _)) in benchmark_rectangles("quspin", block_size).items():
        with open(output_file, "a") as output:
            output.write("{}\t{}\t{}\t{}\t{}\n".format(key, rs[key][0], mean, mean / rs[key][0], 0.0))


if __name__ == "__main__":
    main()
