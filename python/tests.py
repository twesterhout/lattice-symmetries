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
        for coeff, state in operator.apply_off_diag_to_basis_state(ket_state):
            if coeff != 0:
                operator_dict[coeff].append(state)
        coeff = operator.apply_diag_to_basis_state(ket_state)
        if coeff != 0:
            operator_dict[coeff].append(ket_state)
        for amplitude, states in operator_dict.items():
            for bra_state in basis_states:
                if bra_state in states:
                    hamiltonian_matrix[basis_states.index(bra_state)][
                        basis_states.index(ket_state)] = amplitude * states.count(bra_state)

    return hamiltonian_matrix


if __name__ == '__main__':
    basis = lattice_symmetries.SpinfulFermionBasis(2, 2)
    basis.build()
    t = 1
    U = 2
    op = create_hubbard_hamiltonian(basis, t, U)
    print(op)
    h_matrix = calculate_hamiltonian_matrix(op, list(basis.states))
    print(h_matrix)
