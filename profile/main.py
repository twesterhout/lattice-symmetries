import time
import sys
import os
import timeit
import subprocess
from loguru import logger
import numpy as np

import lattice_symmetries as ls


def make_basis(L_x, L_y, sectors=dict()):
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

    hamming_weight = (L_x * L_y) // 2
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
    basis = ls.SpinBasis(
        group,
        number_spins=L_x * L_y,
        hamming_weight=hamming_weight,
        spin_inversion=spin_inversion,
    )
    basis.build()
    return basis


def make_operator(L_x, L_y, basis):
    assert L_x > 0 and L_y > 0
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

    edges = list(generate_nearest_neighbours())
    matrix = np.array(
        [[1, 0, 0, 0], [0, -1, 2, 0], [0, 2, -1, 0], [0, 0, 0, 1]], dtype=np.complex128
    )
    operator = ls.Operator(basis, [ls.Interaction(matrix, edges)])
    return operator


def do_benchmark(L_x, L_y):
    basis = make_basis(L_x, L_y)
    operator = make_operator(L_x, L_y, basis)

    block_size = 1
    n = basis.number_states
    x = 0.5 - np.random.rand(n, block_size)
    x = np.asfortranarray(x)

    tick = time.time()
    for _ in range(10):
        x = operator(x)
    tock = time.time()
    print(tock - tick)


def main():
    do_benchmark(6, 4)


if __name__ == "__main__":
    main()
