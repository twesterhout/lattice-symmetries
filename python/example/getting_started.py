from functools import reduce
import lattice_symmetries as ls
import numpy as np
import operator
import scipy.sparse.linalg


def main():
    number_spins = 10  # System size
    hamming_weight = number_spins // 2  # Hamming weight (i.e. number of spin ups)

    # Constructing symmetries
    sites = np.arange(number_spins)
    # Momentum in x direction with phase φ=1/2 (i.e., with the eigenvalue λ=exp(-2πiφ)=-1)
    T = (ls.Permutation((sites + 1) % number_spins), ls.Rational(1, 2))
    # Parity with eigenvalue λ=-1
    P = (ls.Permutation(sites[::-1]), ls.Rational(1, 2))

    # Constructing the basis
    basis = ls.SpinBasis(
        number_spins=number_spins,
        hamming_weight=hamming_weight,
        spin_inversion=-1,
        symmetries=[T, P],
    )
    basis.build()  # Build the list of representatives, we need it since we're doing ED
    print("Hilbert space dimension is {}".format(basis.number_states))

    edges = [(i, (i + 1) % number_spins) for i in range(number_spins)]
    expr = ls.Expr("2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁", sites=edges)
    print("Expression:", expr)

    # Construct the Hamiltonian
    hamiltonian = ls.Operator(expr, basis)

    # Diagonalize the Hamiltonian using ARPACK
    eigenvalues, eigenstates = scipy.sparse.linalg.eigsh(hamiltonian, k=1, which="SA")
    print("Ground state energy is {}".format(eigenvalues[0]))
    assert np.isclose(eigenvalues[0], -18.06178542)

    # We can also explicitly build the Hamiltonian matrix first, before diagonalizing:
    hamiltonian.build_matrix()
    eigenvalues, _ = scipy.sparse.linalg.eigsh(hamiltonian, k=1, which="SA")
    assert np.isclose(eigenvalues[0], -18.06178542)


if __name__ == "__main__":
    main()
