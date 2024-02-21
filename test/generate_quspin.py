import quspin
from quspin.operators import quantum_LinearOperator, hamiltonian
from quspin.basis import spin_basis_general
import numpy as np
from typing import Any, List, Optional, Tuple
import json

rng = np.random.default_rng(seed=42)


def quspin_static_to_ls(terms) -> str:
    ls_str = ""
    mapping = {
        "I": lambda _i: "I",
        "+": lambda i: "σ⁺_{}".format(i),
        "-": lambda i: "σ⁻_{}".format(i),
        "z": lambda i: "σᶻ_{}".format(i),
    }
    for quspin_string, _term in terms:
        for coeff, *indices in _term:
            term_str = " ".join(
                [mapping[ch](i) for (ch, i) in zip(quspin_string, indices)]
            )
            if len(ls_str) == 0:
                ls_str += "{} {}".format(coeff, term_str)
            else:
                if coeff > 0:
                    ls_str += " + {} {}".format(coeff, term_str)
                else:
                    ls_str += " {} {}".format(coeff, term_str)
    return ls_str


class TermBuilder:
    def __init__(self, complex_coeff=False):
        self.coeff = rng.random() - 0.5
        if complex_coeff:
            self.coeff = self.coeff + (rng.random() - 0.5) * 1j

        self.quspin_str = ""
        self.quspin_indices = []
        # self.ls_strs = []

    def spsm(self, i, j):
        self.quspin_str += "+-"
        self.quspin_indices.append(i)
        self.quspin_indices.append(j)
        # self.ls_strs.append("σ⁺_{} σ⁻_{}".format(i, j))
        return 2

    def smsp(self, i, j):
        self.quspin_str += "-+"
        self.quspin_indices.append(i)
        self.quspin_indices.append(j)
        # self.ls_strs.append("σ⁻_{} σ⁺_{}".format(i, j))
        return 2

    def sp(self, i):
        self.quspin_str += "+"
        self.quspin_indices.append(i)
        # self.ls_strs.append("σ⁺_{}".format(i))
        return 1

    def sm(self, i):
        self.quspin_str += "-"
        self.quspin_indices.append(i)
        # self.ls_strs.append("σ⁻_{}".format(i))
        return 1

    def sz(self, i):
        self.quspin_str += "z"
        self.quspin_indices.append(i)
        # self.ls_strs.append("σᶻ_{}".format(i))
        return 1

    def id(self, i):
        self.quspin_str += "I"
        self.quspin_indices.append(i)
        # self.ls_strs.append("I")
        return 1

    def build(self):
        return [self.quspin_str, [[self.coeff] + self.quspin_indices]]
        # ls_term = str(coeff) + " " + " ".join(self.ls_strs)
        # return [(quspin_term, ls_term)]
        # return [quspin_term]

    def invert_spin(self):
        mapping = {
            "+": (1, "-"),
            "-": (1, "+"),
            "z": (-1, "z"),
            "I": (1, "I"),
        }
        coeff = 1
        chars = []
        for z, c in (mapping[c] for c in self.quspin_str):
            coeff *= z
            chars.append(c)

        other = TermBuilder()
        other.coeff = coeff * self.coeff
        other.quspin_str = "".join(chars)
        other.quspin_indices = self.quspin_indices
        return other

    def permute(self, permutation):
        other = TermBuilder()
        other.coeff = self.coeff
        other.quspin_str = self.quspin_str
        other.quspin_indices = [
            int(np.where(permutation == i)[0][0]) for i in self.quspin_indices
        ]
        print(other.quspin_indices)
        return other


# class SpinInversionTermBuilder:
#     def __init__(self):
#         self.up = TermBuilder()
#         self.down = TermBuilder()
#         self.coeff = 1
#
#     def spsm(self, i, j):
#         self.up.spsm(i, j)
#         self.down.smsp(i, j)
#         return 2
#
#     def smsp(self, i, j):
#         self.up.smsp(i, j)
#         self.down.spsm(i, j)
#         return 2
#
#     def sp(self, i):
#         self.up.sp(i)
#         self.down.sm(i)
#         return 1
#
#     def sm(self, i):
#         self.up.sm(i)
#         self.down.sp(i)
#         return 1
#
#     def sz(self, i):
#         self.up.sz(i)
#         self.down.sz(i)
#         self.coeff *= -1
#         return 1
#
#     def id(self, i):
#         self.up.id(i)
#         self.down.id(i)
#         return 1
#
#     def build(self, coeff):
#         a = self.up.build(coeff)
#         b = self.down.build(self.coeff * coeff)
#         return a + b
#
#     def permute(self, permutation):
#         other = SpinInversionTermBuilder()


def random_term(
    order: int,
    number_sites: int,
    hamming_weight: Optional[int],
    complex_coeff: bool,
):
    get_index = lambda: rng.integers(0, number_sites)
    builder = TermBuilder(complex_coeff=complex_coeff)
    # if spin_inversion is None else SpinInversionTermBuilder()

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

    return builder


def random_operator(
    number_terms: int,
    max_order: int,
    number_sites: int,
    hamming_weight: Optional[int],
    spin_inversion: Optional[int],
    symmetry,
    complex_coeff,
):
    if spin_inversion is not None:
        number_terms = number_terms // 2

    terms = []
    for _ in range(number_terms):
        order = rng.integers(1, max_order)
        terms.append(random_term(order, number_sites, hamming_weight, complex_coeff))

    additional = []
    if symmetry is not None:
        _name, permutation, sector = symmetry
        for p, _, _ in build_cyclic_group(permutation, sector):
            if not np.array_equal(p, np.arange(p.size)):
                additional += [t.permute(p) for t in terms]
    terms += additional

    if spin_inversion is not None:
        terms += [t.invert_spin() for t in terms]

    terms = [t.build() for t in terms]
    print(terms)
    ls_str = quspin_static_to_ls(terms)
    # print(ls_str)
    # ls_str = ""
    # for _, s in terms:
    #     if len(ls_str) > 0:
    #         ls_str += " " if s.startswith("-") else " + "
    #     ls_str += s

    ls_expr = {
        "basis": dict(
            number_spins=number_sites,
            hamming_weight=hamming_weight,
            spin_inversion=spin_inversion,
            symmetries=[dict(permutation=symmetry[1], sector=symmetry[2])]
            if symmetry is not None
            else [],
        ),
        "expression": {"expression": ls_str, "particle": "spin-1/2"},
    }

    Nup = number_sites - hamming_weight if hamming_weight is not None else None
    if spin_inversion is not None:
        blocks = dict(zblock=([-(i + 1) for i in range(number_sites)], spin_inversion))
    else:
        blocks = dict()
    if symmetry is not None:
        name, permutation, sector = symmetry
        blocks[name] = (permutation, sector)

    basis = spin_basis_general(N=number_sites, Nup=Nup, pauli=-1, **blocks)
    quspin_op = hamiltonian(
        static_list=terms,
        dynamic_list=[],
        basis=basis,
        dtype=np.float64 if not complex_coeff else np.complex128,
        check_herm=False,
        check_symm=True,  # False,
    )

    return quspin_op, ls_expr


def reverse_bits(x, n_bits):
    x = np.array(x)
    x_reversed = np.zeros_like(x)
    for _ in range(n_bits):
        x_reversed = (x_reversed << 1) | x & 1
        x >>= 1
    return x_reversed


def permute_bits(x, permutation):
    x = np.array(x)
    x_permuted = np.zeros_like(x)
    for i in range(permutation.size):
        x_permuted |= ((x >> permutation[i]) & 1) << i
    print(x, permutation, x_permuted)
    return x_permuted


def build_cyclic_group(permutation, sector):
    permutation = np.asarray(permutation)
    identity = np.arange(permutation.size)

    elements = [identity]
    sectors = [0]

    p = permutation
    s = sector
    n = 1
    while not np.array_equal(p, identity):
        elements.append(p)
        sectors.append(s)

        p = p[permutation]
        s += sector
        n += 1

    sectors = [s % n for s in sectors]
    characters = [np.exp(-2j * np.pi * s / n) for s in sectors]
    return zip(elements, sectors, characters)


def generate_test_file(
    number_terms: int,
    max_order: int,
    number_sites: int,
    hamming_weight: Optional[int],
    spin_inversion: Optional[int],
    symmetry=None,
    complex_coeff=False,
):
    quspin_op, ls_expr = random_operator(
        number_terms,
        max_order,
        number_sites,
        hamming_weight,
        spin_inversion,
        symmetry,
        complex_coeff,
    )

    basis = quspin_op.basis
    x = np.random.rand(basis.Ns) + np.random.rand(basis.Ns) * 1j - (0.5 + 0.5j)
    y = quspin_op.dot(x)

    name = "test_{}_{}_{}".format(number_sites, hamming_weight, spin_inversion)
    with open("{}_expr.json".format(name), "w") as out:
        json.dump(ls_expr, out)

    # QuSpin orders basis states dirrefently:
    #   - It stores the spin 0 in the most significant bit, and lattice-symmetries in the least significant bit
    basis_states = reverse_bits(basis.states, basis.N)
    #   - It uses 1 to represent ↑, and 0 to represent ↓, but lattice symmetries does the inverse 😭
    basis_states = basis_states ^ ((1 << basis.N) - 1)

    if symmetry is not None:
        assert spin_inversion is None
        _, permutation, sector = symmetry
        representatives = basis_states
        characters = np.ones(basis_states.size)
        for p, _, c in build_cyclic_group(permutation, sector):
            permuted = permute_bits(basis_states, p)
            representatives = np.minimum(representatives, permuted)
            characters[permuted < representatives] = c

        basis_states = representatives
        x *= characters
        y *= characters

    if "zblock" in basis.blocks:
        (_, spin_inversion) = basis.blocks["zblock"]
        mask = (1 << basis.L) - 1
        need_flipping = (basis_states ^ mask) < basis_states
        x[need_flipping] *= spin_inversion
        y[need_flipping] *= spin_inversion
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
    generate_test_file(
        3,
        4,
        4,
        hamming_weight=None,
        spin_inversion=None,
        symmetry=("Tx", [1, 2, 3, 0], 0),
        complex_coeff=False,
    )
    # for number_sites in [1, 2, 3, 4, 5, 10]:
    #     number_terms = 4 * number_sites
    #     max_order = 10

    #     for spin_inversion in [None, 1, -1]:
    #         generate_test_file(
    #             number_terms,
    #             max_order,
    #             number_sites,
    #             hamming_weight=None,
    #             spin_inversion=spin_inversion,
    #         )

    #     for hamming_weight in range(0, number_sites + 1):
    #         generate_test_file(
    #             number_terms,
    #             max_order,
    #             number_sites,
    #             hamming_weight=hamming_weight,
    #             spin_inversion=None,
    #         )
    #         if 2 * hamming_weight == number_sites:
    #             for spin_inversion in [1, -1]:
    #                 generate_test_file(
    #                     number_terms,
    #                     max_order,
    #                     number_sites,
    #                     hamming_weight=hamming_weight,
    #                     spin_inversion=spin_inversion,
    #                 )


if __name__ == "__main__":
    generate_tests()
