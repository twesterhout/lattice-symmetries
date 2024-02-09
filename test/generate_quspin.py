import quspin
from quspin.operators import quantum_LinearOperator, hamiltonian
from quspin.basis import spin_basis_1d
import numpy as np
import json

rng = np.random.default_rng(seed=42)


def random_M_conserving_term(order: int, number_sites: int):
    get_index = lambda: rng.integers(0, number_sites)

    quspin_str = ""
    quspin_indices = []
    ls_strs = []

    def spsm():
        nonlocal quspin_str, quspin_indices, ls_strs
        i = get_index()
        j = get_index()
        quspin_str += "+-"
        quspin_indices.append(i)
        quspin_indices.append(j)
        ls_strs.append("σ⁺_{} σ⁻_{}".format(i, j))
        return 2

    def sz():
        nonlocal quspin_str, quspin_indices, ls_strs
        i = get_index()
        quspin_str += "z"
        quspin_indices.append(i)
        ls_strs.append("σᶻ_{}".format(i))
        return 1

    # def id():
    #     nonlocal quspin_str, quspin_indices, ls_strs
    #     i = get_index()
    #     quspin_str += "I"
    #     quspin_indices.append(i)
    #     ls_strs.append("I")
    #     return 1

    while order > 0:
        choices = [sz]  # , id]
        if order >= 2:
            choices.append(spsm)
        n = choices[rng.integers(0, len(choices))]()
        order -= n

    coeff = rng.random() - 0.5
    quspin_term = [quspin_str, [[coeff] + quspin_indices]]
    ls_term = str(coeff) + " " + " ".join(ls_strs)
    print(quspin_term)
    print(ls_term)

    return quspin_term, ls_term


def random_M_conserving_operator(
    number_terms: int, max_order: int, number_sites: int, hamming_weight: int
):
    terms = [
        random_M_conserving_term(rng.integers(1, max_order), number_sites)
        for _ in range(number_terms)
    ]

    basis = spin_basis_1d(L=number_sites, Nup=number_sites - hamming_weight, pauli=-1)
    quspin_op = hamiltonian(
        static_list=[t for (t, _) in terms],
        dynamic_list=[],
        basis=basis,
        dtype=np.float64,
        check_herm=False,
    )
    # print(quspin_op.toarray())

    ls_str = ""
    for _, s in terms:
        if len(ls_str) > 0:
            ls_str += " " if s.startswith("-") else " + "
        ls_str += s

    ls_expr = {
        "basis": dict(number_spins=number_sites, hamming_weight=hamming_weight),
        "expression": {"expression": ls_str, "particle": "spin-1/2"},
    }
    return quspin_op, ls_expr


def reverse_bits(x, n_bits):
    x = np.array(x)
    x_reversed = np.zeros_like(x)
    for _ in range(n_bits):
        x_reversed = (x_reversed << 1) | x & 1
        x >>= 1
    return x_reversed


def generate_test_file(name, quspin_op, ls_expr):
    basis = quspin_op.basis
    x = np.random.rand(basis.Ns) + np.random.rand(basis.Ns) * 1j - (0.5 + 0.5j)
    y = quspin_op.dot(x)

    with open("{}_expr.json".format(name), "w") as out:
        json.dump(ls_expr, out)

    # QuSpin orders basis states dirrefently:
    #   - It stores the spin 0 in the most significant bit, and lattice-symmetries in the least significant bit
    basis_states = reverse_bits(basis.states, basis.L)
    #   - It uses 1 to represent ↑, and 0 to represent ↓, but lattice symmetries does the inverse 😭
    basis_states = basis_states ^ ((1 << basis.L) - 1)

    permutation = np.argsort(basis_states)
    x = x[permutation]
    y = y[permutation]

    arr_info = {
        "x_real": [float(c.real) for c in x],
        "x_imag": [float(c.imag) for c in x],
        "y_real": [float(c.real) for c in y],
        "y_imag": [float(c.imag) for c in y],
    }
    with open("{}_arrays.json".format(name), "w") as out:
        json.dump(arr_info, out)


def generate_M_conserving_tests():
    for number_sites in [1, 2, 3, 4, 5, 10, 16]:
        number_terms = 4 * number_sites
        max_order = 10
        for hamming_weight in range(0, number_sites + 1):
            generate_test_file(
                f"test_{number_sites}_{hamming_weight}_{max_order}_{number_terms}",
                *random_M_conserving_operator(
                    number_terms, max_order, number_sites, hamming_weight
                ),
            )


if __name__ == "__main__":
    generate_M_conserving_tests()
