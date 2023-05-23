from collections import defaultdict

import lattice_symmetries as ls
import numpy as np
import scipy.sparse.linalg


def test_symmetry():
    a = ls.Symmetry([0, 1, 2], sector=0)
    assert a.sector == 0
    assert len(a) == 3
    assert a.permutation.tolist() == [0, 1, 2]
    b = ls.Symmetry(a._payload, None)
    assert b.sector == 0
    assert len(b) == 3
    assert b.permutation.tolist() == [0, 1, 2]
    del b
    del a


def test_symmetries():
    a = ls.Symmetry([1, 2, 3, 0], sector=0)
    b = ls.Symmetry([3, 2, 1, 0], sector=0)
    c = ls.Symmetries([a, b])


def test_index():
    basis = ls.SpinBasis(4)
    basis.build()
    assert np.array_equal(basis.index(basis.states), basis.states)
    assert np.array_equal(basis.index(basis.states[2]), 2)


def test_kagome_symmetries():
    expr = ls.Expr(
        "1.0 × σᶻ₀ σᶻ₁ + 1.0 × σᶻ₀ σᶻ₃ + 1.0 × σᶻ₀ σᶻ₈ + 1.0 × σᶻ₀ σᶻ₁₀ + 2.0 × σ⁺₀ σ⁻₁ + 2.0 × σ⁺₀ σ⁻₃ + 2.0 × σ⁺₀ σ⁻₈ + 2.0 × σ⁺₀ σ⁻₁₀ + 2.0 × σ⁻₀ σ⁺₁ + 2.0 × σ⁻₀ σ⁺₃ + 2.0 × σ⁻₀ σ⁺₈ + 2.0 × σ⁻₀ σ⁺₁₀ + 1.0 × σᶻ₁ σᶻ₂ + 0.8 × σᶻ₁ σᶻ₃ + 0.8 × σᶻ₁ σᶻ₉ + 2.0 × σ⁺₁ σ⁻₂ + 1.6 × σ⁺₁ σ⁻₃ + 1.6 × σ⁺₁ σ⁻₉ + 2.0 × σ⁻₁ σ⁺₂ + 1.6 × σ⁻₁ σ⁺₃ + 1.6 × σ⁻₁ σ⁺₉ + 1.0 × σᶻ₂ σᶻ₄ + 1.0 × σᶻ₂ σᶻ₉ + 1.0 × σᶻ₂ σᶻ₁₀ + 2.0 × σ⁺₂ σ⁻₄ + 2.0 × σ⁺₂ σ⁻₉ + 2.0 × σ⁺₂ σ⁻₁₀ + 2.0 × σ⁻₂ σ⁺₄ + 2.0 × σ⁻₂ σ⁺₉ + 2.0 × σ⁻₂ σ⁺₁₀ + 1.0 × σᶻ₃ σᶻ₅ + 0.8 × σᶻ₃ σᶻ₁₁ + 2.0 × σ⁺₃ σ⁻₅ + 1.6 × σ⁺₃ σ⁻₁₁ + 2.0 × σ⁻₃ σ⁺₅ + 1.6 × σ⁻₃ σ⁺₁₁ + 0.8 × σᶻ₄ σᶻ₆ + 1.0 × σᶻ₄ σᶻ₇ + 0.8 × σᶻ₄ σᶻ₁₀ + 1.6 × σ⁺₄ σ⁻₆ + 2.0 × σ⁺₄ σ⁻₇ + 1.6 × σ⁺₄ σ⁻₁₀ + 1.6 × σ⁻₄ σ⁺₆ + 2.0 × σ⁻₄ σ⁺₇ + 1.6 × σ⁻₄ σ⁺₁₀ + 1.0 × σᶻ₅ σᶻ₆ + 1.0 × σᶻ₅ σᶻ₈ + 1.0 × σᶻ₅ σᶻ₁₁ + 2.0 × σ⁺₅ σ⁻₆ + 2.0 × σ⁺₅ σ⁻₈ + 2.0 × σ⁺₅ σ⁻₁₁ + 2.0 × σ⁻₅ σ⁺₆ + 2.0 × σ⁻₅ σ⁺₈ + 2.0 × σ⁻₅ σ⁺₁₁ + 1.0 × σᶻ₆ σᶻ₇ + 0.8 × σᶻ₆ σᶻ₈ + 2.0 × σ⁺₆ σ⁻₇ + 1.6 × σ⁺₆ σ⁻₈ + 2.0 × σ⁻₆ σ⁺₇ + 1.6 × σ⁻₆ σ⁺₈ + 1.0 × σᶻ₇ σᶻ₉ + 1.0 × σᶻ₇ σᶻ₁₁ + 2.0 × σ⁺₇ σ⁻₉ + 2.0 × σ⁺₇ σ⁻₁₁ + 2.0 × σ⁻₇ σ⁺₉ + 2.0 × σ⁻₇ σ⁺₁₁ + 0.8 × σᶻ₈ σᶻ₁₀ + 1.6 × σ⁺₈ σ⁻₁₀ + 1.6 × σ⁻₈ σ⁺₁₀ + 0.8 × σᶻ₉ σᶻ₁₁ + 1.6 × σ⁺₉ σ⁻₁₁ + 1.6 × σ⁻₉ σ⁺₁₁"
    )
    # top_shift = ls.Symmetry([5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 11, 10], sector=0)
    right_shift = ls.Symmetry([2, 10, 0, 4, 3, 7, 11, 5, 9, 8, 1, 6], sector=1)
    assert expr == expr.replace_indices(dict(zip(range(12), right_shift.permutation)))
    symmetries = ls.Symmetries([right_shift])
    basis = ls.SpinBasis(
        symmetries=symmetries, number_spins=12, hamming_weight=6, spin_inversion=None
    )
    basis.build()
    print(basis.states)
    print(basis.state_info(basis.states))
    hamiltonian = ls.Operator(basis, expr)
    energy, state = scipy.sparse.linalg.eigsh(hamiltonian, k=1, which="SA")


def test_issue01():
    def create_operator(basis):
        a = ls.Expr("n↑₀ n↓₀", [(0,), (1,)])
        b = ls.Expr("c↑₀ c†↑₁", [(0, 1)])
        c = ls.Expr("c↓₀ c†↓₁", [(0, 1)])
        return ls.Operator(basis, 2 * a + b + b.adjoint() + c + c.adjoint())

    basis = ls.SpinfulFermionBasis(2)
    basis.build()
    op = create_operator(basis)

    for x in basis.states:
        print("Applying to {} ...".format(x))
        op.apply_diag_to_basis_state(x)


def create_hubbard_hamiltonian(basis, t, U):
    operator = -t * ls.Expr("c†↑₀ c↑₁", [(0, 1)])
    operator -= t * ls.Expr("c†↑₁ c↑₀", [(0, 1)])
    operator -= t * ls.Expr("c†↓₁ c↓₀", [(0, 1)])
    operator -= t * ls.Expr("c†↓₀ c↓₁", [(0, 1)])
    operator += U * ls.Expr("n↑₀ n↓₀", [(0,)])
    operator += U * ls.Expr("n↑₁ n↓₁", [(1,)])

    return ls.Operator(basis, operator)


def calculate_hamiltonian_matrix(operator, basis_states):
    hamiltonian_matrix = np.zeros((len(basis_states), len(basis_states)))
    for ket_state in basis_states:
        operator_dict = defaultdict(list)
        coeff = operator.apply_diag_to_basis_state(ket_state)
        if coeff != 0:
            operator_dict[coeff].append(ket_state)
        for coeff, state in operator.apply_off_diag_to_basis_state(ket_state):
            if coeff != 0:
                operator_dict[coeff].append(state)

        for amplitude, states in operator_dict.items():
            for bra_state in basis_states:
                if bra_state in states:
                    c = amplitude * states.count(bra_state)
                    assert c.imag == 0
                    hamiltonian_matrix[basis_states.index(bra_state)][
                        basis_states.index(ket_state)
                    ] = c.real

    return hamiltonian_matrix


def test_constructed_basis_states():
    basis = ls.SpinfulFermionBasis(2)
    basis.build()
    assert np.array_equal(
        basis.states, np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])
    )


def test_constructed_basis_states_with_one_particle():
    basis = ls.SpinfulFermionBasis(2, 1)
    basis.build()
    assert np.array_equal(basis.states, np.array([1, 2, 4, 8]))


def test_constructed_operator():
    basis = ls.SpinfulFermionBasis(2)
    basis.build()
    t = 1
    U = 2
    op = create_hubbard_hamiltonian(basis, t, U)
    assert (
        "2.0 × n↑₀ n↓₀ + -1.0 × c†↑₀ c↑₁ + 1.0 × c↑₀ c†↑₁ + 2.0 × n↑₁ n↓₁ + -1.0 × c†↓₀ c↓₁ + 1.0 × c↓₀ c†↓₁"
        == str(op.expression)
    )


def test_hubbard_hamiltonian_matrix():
    basis = ls.SpinfulFermionBasis(2)
    basis.build()
    t = 1
    U = 2
    op = create_hubbard_hamiltonian(basis, t, U)
    hamiltonian_matrix = calculate_hamiltonian_matrix(op, list(basis.states))
    valid_matrix = np.array(
        [
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 2.0, -1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, -1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, -1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 2.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0],
        ]
    )

    assert np.array_equal(hamiltonian_matrix, valid_matrix)


def test_hubbard_hamiltonian_matrix_with_defined_particle_number():
    basis = ls.SpinfulFermionBasis(2, 2)
    basis.build()
    t = 1
    U = 2
    op = create_hubbard_hamiltonian(basis, t, U)
    hamiltonian_matrix = calculate_hamiltonian_matrix(op, list(basis.states))
    valid_matrix = np.array(
        [
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 2.0, -1.0, -1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0, 0.0, -1.0, 0.0],
            [0.0, -1.0, 0.0, 0.0, -1.0, 0.0],
            [0.0, 0.0, -1.0, -1.0, 2.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ]
    )

    assert np.array_equal(hamiltonian_matrix, valid_matrix)
