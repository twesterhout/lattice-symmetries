import numpy as np
from collections import defaultdict
import lattice_symmetries


def create_hubbard_hamiltonian(basis, t, U):
    operator = -t * lattice_symmetries.Operator(basis, "c†↑₀ c↑₁", [(0, 1)])
    operator -= t * lattice_symmetries.Operator(basis, "c†↑₁ c↑₀", [(0, 1)])
    operator -= t * lattice_symmetries.Operator(basis, "c†↓₁ c↓₀", [(0, 1)])
    operator -= t * lattice_symmetries.Operator(basis, "c†↓₀ c↓₁", [(0, 1)])

    operator += U * lattice_symmetries.Operator(basis, "n↑₀ n↓₀", [(0)])
    operator += U * lattice_symmetries.Operator(basis, "n↑₁ n↓₁", [(1)])

    return operator


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
                    hamiltonian_matrix[basis_states.index(bra_state)][
                        basis_states.index(ket_state)] = amplitude * states.count(bra_state)

    return hamiltonian_matrix


def test_constructed_basis_states():
    basis = lattice_symmetries.SpinfulFermionBasis(2)
    basis.build()
    assert np.array_equal(basis.states, np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]))


def test_constructed_basis_states_with_one_particle():
    basis = lattice_symmetries.SpinfulFermionBasis(2, 1)
    basis.build()
    assert np.array_equal(basis.states, np.array([1, 2, 4, 8]))


def test_constructed_operator():
    basis = lattice_symmetries.SpinfulFermionBasis(2)
    basis.build()
    t = 1
    U = 2
    op = create_hubbard_hamiltonian(basis, t, U)
    assert "2.0 × n↑₀ n↓₀ + -1.0 × c†↑₀ c↑₁ + 1.0 × c↑₀ c†↑₁ + 2.0 × n↑₁ n↓₁ + -1.0 × c†↓₀ c↓₁ + 1.0 × c↓₀ c†↓₁" == \
           str(op)


def test_hubbard_hamiltonian_matrix():
    basis = lattice_symmetries.SpinfulFermionBasis(2)
    basis.build()
    t = 1
    U = 2
    op = create_hubbard_hamiltonian(basis, t, U)
    hamiltonian_matrix = calculate_hamiltonian_matrix(op, list(basis.states))
    valid_matrix = np.array([
        [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
        [0., 0., -1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
        [0., -1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
        [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
        [0., 0., 0., 0., 0., 0., 0., 0., -1., 0., 0., 0., 0., 0., 0., 0.],
        [0., 0., 0., 0., 0., 2., -1., 0., 0., -1., 0., 0., 0., 0., 0., 0.],
        [0., 0., 0., 0., 0., -1., 0., 0., 0., 0., -1., 0., 0., 0., 0., 0.],
        [0., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0., -1., 0., 0., 0., 0.],
        [0., 0., 0., 0., -1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
        [0., 0., 0., 0., 0., -1., 0., 0., 0., 0., -1., 0., 0., 0., 0., 0.],
        [0., 0., 0., 0., 0., 0., -1., 0., 0., -1., 2., 0., 0., 0., 0., 0.],
        [0., 0., 0., 0., 0., 0., 0., -1., 0., 0., 0., 2., 0., 0., 0., 0.],
        [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
        [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2., -1., 0.],
        [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 2., 0.],
        [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 4.]])

    assert np.array_equal(hamiltonian_matrix, valid_matrix)


def test_hubbard_hamiltonian_matrix_with_defined_particle_number():
    basis = lattice_symmetries.SpinfulFermionBasis(2, 2)
    basis.build()
    t = 1
    U = 2
    op = create_hubbard_hamiltonian(basis, t, U)
    hamiltonian_matrix = calculate_hamiltonian_matrix(op, list(basis.states))
    valid_matrix = np.array(
        [[0., 0., 0., 0., 0., 0.],
         [0., 2., -1., -1., 0., 0.],
         [0., -1., 0., 0., -1., 0.],
         [0., -1., 0., 0., -1., 0.],
         [0., 0., -1., -1., 2., 0.],
         [0., 0., 0., 0., 0., 0.]])

    assert np.array_equal(hamiltonian_matrix, valid_matrix)


