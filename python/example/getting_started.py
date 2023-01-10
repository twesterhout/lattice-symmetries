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
    # Momentum in x direction with eigenvalue π
    T = ls.Symmetry((sites + 1) % number_spins, sector=number_spins // 2)
    # Parity with eigenvalue π
    P = ls.Symmetry(sites[::-1], sector=1)

    # Constructing the group
    symmetries = ls.Symmetries([T, P])
    print("Symmetry group contains {} elements".format(len(symmetries)))

    # Constructing the basis
    basis = ls.SpinBasis(
        number_spins=number_spins,
        # NOTE: we don't actually need to specify hamming_weight when spin_inversion
        # is set. The library will guess that hamming_weight = number_spins / 2.
        spin_inversion=-1,
        symmetries=symmetries,
    )
    basis.build()  # Build the list of representatives, we need it since we're doing ED
    print("Hilbert space dimension is {}".format(basis.number_states))

    edges = [(i, (i + 1) % number_spins) for i in range(number_spins)]
    expr = ls.Expr("2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁", sites=edges)
    print("Expression:", expr)

    # Alternatively, we can create expr in an algebraic way:
    two_site_expr = ls.Expr("2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁")
    many_exprs = [two_site_expr.replace_indices({0: i, 1: j}) for (i, j) in edges]
    expr2 = reduce(operator.add, many_exprs)
    assert expr == expr2

    # Construct the Hamiltonian
    hamiltonian = ls.Operator(basis, expr)

    # Diagonalize the Hamiltonian using ARPACK
    eigenvalues, eigenstates = scipy.sparse.linalg.eigsh(hamiltonian, k=1, which="SA")
    print("Ground state energy is {}".format(eigenvalues[0]))
    assert np.isclose(eigenvalues[0], -18.06178542)


if __name__ == "__main__":
    main()
