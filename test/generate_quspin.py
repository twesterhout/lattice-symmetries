import quspin
from quspin.operators import quantum_LinearOperator, hamiltonian
from quspin.basis import spin_basis_1d
import numpy as np
from typing import Any, List, Optional
import json

rng = np.random.default_rng(seed=42)


class TermBuilder:
    def __init__(self):
        self.quspin_str = ""
        self.quspin_indices = []
        self.ls_strs = []

    def spsm(self, i, j):
        self.quspin_str += "+-"
        self.quspin_indices.append(i)
        self.quspin_indices.append(j)
        self.ls_strs.append("σ⁺_{} σ⁻_{}".format(i, j))
        return 2

    def smsp(self, i, j):
        self.quspin_str += "-+"
        self.quspin_indices.append(i)
        self.quspin_indices.append(j)
        self.ls_strs.append("σ⁻_{} σ⁺_{}".format(i, j))
        return 2

    def sp(self, i):
        self.quspin_str += "+"
        self.quspin_indices.append(i)
        self.ls_strs.append("σ⁺_{}".format(i))
        return 1

    def sm(self, i):
        self.quspin_str += "-"
        self.quspin_indices.append(i)
        self.ls_strs.append("σ⁻_{}".format(i))
        return 1

    def sz(self, i):
        self.quspin_str += "z"
        self.quspin_indices.append(i)
        self.ls_strs.append("σᶻ_{}".format(i))
        return 1

    def id(self, i):
        self.quspin_str += "I"
        self.quspin_indices.append(i)
        self.ls_strs.append("I")
        return 1

    def build(self, coeff):
        quspin_term = [self.quspin_str, [[coeff] + self.quspin_indices]]
        ls_term = str(coeff) + " " + " ".join(self.ls_strs)
        return [(quspin_term, ls_term)]


class SpinInversionTermBuilder:
    def __init__(self):
        self.up = TermBuilder()
        self.down = TermBuilder()
        self.coeff = 1

    def spsm(self, i, j):
        self.up.spsm(i, j)
        self.down.smsp(i, j)
        return 2

    def smsp(self, i, j):
        self.up.smsp(i, j)
        self.down.spsm(i, j)
        return 2

    def sp(self, i):
        self.up.sp(i)
        self.down.sm(i)
        return 1

    def sm(self, i):
        self.up.sm(i)
        self.down.sp(i)
        return 1

    def sz(self, i):
        self.up.sz(i)
        self.down.sz(i)
        self.coeff *= -1
        return 1

    def id(self, i):
        self.up.id(i)
        self.down.id(i)
        return 1

    def build(self, coeff):
        a = self.up.build(coeff)
        b = self.down.build(self.coeff * coeff)
        return a + b


def random_term(
    order: int,
    number_sites: int,
    hamming_weight: Optional[int],
    spin_inversion: Optional[int],
    complex_coeff=False,
):
    get_index = lambda: rng.integers(0, number_sites)
    builder = TermBuilder() if spin_inversion is None else SpinInversionTermBuilder()

    while order > 0:
        choices: List[Any] = [
            lambda: builder.sz(get_index()),
            lambda: builder.id(get_index()),
        ]
        if hamming_weight is None:
            choices += [
                lambda: builder.sp(get_index()),
                lambda: builder.sm(get_index()),
            ]
        elif order >= 2:
            choices += [lambda: builder.spsm(get_index(), get_index())]

        k = rng.integers(0, len(choices))
        n = choices[k]()
        order -= n

    coeff = rng.random() - 0.5
    if complex_coeff:
        coeff = coeff + (rng.random() - 0.5) * 1j

    return builder.build(coeff)


def random_operator(
    number_terms: int,
    max_order: int,
    number_sites: int,
    hamming_weight: Optional[int],
    spin_inversion: Optional[int],
):
    if spin_inversion is not None:
        number_terms = number_terms // 2

    terms = []
    for _ in range(number_terms):
        order = rng.integers(1, max_order)
        terms += random_term(order, number_sites, hamming_weight, spin_inversion)

    ls_str = ""
    for _, s in terms:
        if len(ls_str) > 0:
            ls_str += " " if s.startswith("-") else " + "
        ls_str += s

    ls_expr = {
        "basis": dict(
            number_spins=number_sites,
            hamming_weight=hamming_weight,
            spin_inversion=spin_inversion,
        ),
        "expression": {"expression": ls_str, "particle": "spin-1/2"},
    }

    Nup = number_sites - hamming_weight if hamming_weight is not None else None
    if spin_inversion is not None:
        blocks = dict(zblock=spin_inversion)
    else:
        blocks = dict()
    basis = spin_basis_1d(L=number_sites, Nup=Nup, pauli=-1, **blocks)
    quspin_op = hamiltonian(
        static_list=[t for (t, _) in terms],
        dynamic_list=[],
        basis=basis,
        dtype=np.float64,
        check_herm=False,
        check_symm=False,
    )

    return quspin_op, ls_expr


def reverse_bits(x, n_bits):
    x = np.array(x)
    x_reversed = np.zeros_like(x)
    for _ in range(n_bits):
        x_reversed = (x_reversed << 1) | x & 1
        x >>= 1
    return x_reversed


def generate_test_file(
    number_terms: int,
    max_order: int,
    number_sites: int,
    hamming_weight: Optional[int],
    spin_inversion: Optional[int],
):
    quspin_op, ls_expr = random_operator(
        number_terms, max_order, number_sites, hamming_weight, spin_inversion
    )

    basis = quspin_op.basis
    x = np.random.rand(basis.Ns) + np.random.rand(basis.Ns) * 1j - (0.5 + 0.5j)
    y = quspin_op.dot(x)

    name = "test_{}_{}_{}".format(number_sites, hamming_weight, spin_inversion)
    with open("{}_expr.json".format(name), "w") as out:
        json.dump(ls_expr, out)

    # QuSpin orders basis states dirrefently:
    #   - It stores the spin 0 in the most significant bit, and lattice-symmetries in the least significant bit
    basis_states = reverse_bits(basis.states, basis.L)
    #   - It uses 1 to represent ↑, and 0 to represent ↓, but lattice symmetries does the inverse 😭
    basis_states = basis_states ^ ((1 << basis.L) - 1)

    if "zblock" in basis.blocks:
        mask = (1 << basis.L) - 1
        need_flipping = (basis_states ^ mask) < basis_states
        x[need_flipping] *= basis.blocks["zblock"]
        y[need_flipping] *= basis.blocks["zblock"]
        basis_states = np.minimum(basis_states, basis_states ^ mask)

    permutation = np.argsort(basis_states)
    basis_states = basis_states[permutation]
    x = x[permutation]
    y = y[permutation]

    arr_info = {
        "states": [int(c) for c in basis_states],
        "x_real": [float(c.real) for c in x],
        "x_imag": [float(c.imag) for c in x],
        "y_real": [float(c.real) for c in y],
        "y_imag": [float(c.imag) for c in y],
    }
    with open("{}_arrays.json".format(name), "w") as out:
        json.dump(arr_info, out)


def generate_tests():
    for number_sites in [1, 2, 3, 4, 5, 10]:
        number_terms = 4 * number_sites
        max_order = 10

        for spin_inversion in [None, 1, -1]:
            generate_test_file(
                number_terms,
                max_order,
                number_sites,
                hamming_weight=None,
                spin_inversion=spin_inversion,
            )

        for hamming_weight in range(0, number_sites + 1):
            generate_test_file(
                number_terms,
                max_order,
                number_sites,
                hamming_weight=hamming_weight,
                spin_inversion=None,
            )
            if 2 * hamming_weight == number_sites:
                for spin_inversion in [1, -1]:
                    generate_test_file(
                        number_terms,
                        max_order,
                        number_sites,
                        hamming_weight=hamming_weight,
                        spin_inversion=spin_inversion,
                    )


if __name__ == "__main__":
    generate_tests()
