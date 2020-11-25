# The following three lines can be commented out
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "..", "python"))

import numpy as np
import scipy.sparse.linalg
import lattice_symmetries as ls


def main():
    number_spins = 10  # System size
    hamming_weight = number_spins // 2  # Hamming weight (i.e. number of spin ups)

    # Constructing symmetries
    symmetries = []
    sites = np.arange(number_spins)
    # Momentum in x direction with eigenvalue π
    T = (sites + 1) % number_spins
    symmetries.append(ls.Symmetry(T, sector=number_spins // 2))
    # Parity with eigenvalue π
    P = sites[::-1]
    symmetries.append(ls.Symmetry(P, sector=1))
    # Global spin inversion with eigenvalue π
    # Z = sites  # no permutation
    # symmetries.append(ls.Symmetry(Z, flip=True, sector=1))

    # Constructing the group
    symmetry_group = ls.Group(symmetries)
    print("Symmetry group contains {} elements".format(len(symmetry_group)))

    # Constructing the basis
    basis = ls.SpinBasis(
        symmetry_group, number_spins=number_spins, hamming_weight=hamming_weight, spin_inversion=-1
    )
    basis.build()  # Build the cache, we need it since we're doing ED
    print("Hilbert space dimension is {}".format(basis.number_states))

    # Heisenberg Hamiltonian
    # fmt: off
    σ_x = np.array([ [0, 1]
                   , [1, 0] ])
    σ_y = np.array([ [0 , -1j]
                   , [1j,   0] ])
    σ_z = np.array([ [1,  0]
                   , [0, -1] ])
    # fmt: on
    σ_p = σ_x + 1j * σ_y
    σ_m = σ_x - 1j * σ_y

    matrix = 0.5 * (np.kron(σ_p, σ_m) + np.kron(σ_m, σ_p)) + np.kron(σ_z, σ_z)
    edges = [
        (i, (i + 1) % number_spins) for i in range(number_spins)
    ]  # periodic boundary conditions
    interaction = ls.Interaction(matrix, edges)
    hamiltonian = ls.Operator(basis, [interaction])

    # Diagonalize the Hamiltonian using ARPACK
    eigenvalues, eigenstates = ls.diagonalize(hamiltonian, k=1)
    print("Ground state energy is {}".format(eigenvalues[0]))
    assert np.isclose(eigenvalues[0], -18.06178542)


if __name__ == "__main__":
    main()
