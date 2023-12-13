import glob
import os
import lattice_symmetries as ls
import math
import numpy as np
import scipy.sparse.linalg
from pytest import approx, raises


def sum1(xs):
    s = None
    for x in xs:
        if s is not None:
            s += x
        else:
            s = x
    return s


def test_symmetry():
    a = ls.Symmetry([0, 1, 2], sector=0)
    assert a.sector == 0
    assert len(a) == 3
    assert a.permutation.tolist() == [0, 1, 2]
    del a

    a = ls.Symmetry(list(range(10)), sector=123)
    assert a.phase == 0

    a = ls.Symmetry(list(range(10)), sector=-1)
    assert a.phase == 0


def test_symmetries():
    a = ls.Symmetry([1, 2, 3, 0], sector=0)
    b = ls.Symmetry([3, 2, 1, 0], sector=0)
    c = ls.Symmetries([a, b])
    assert len(c) == 8

    with raises(SystemError):
        ls.Symmetries([ls.Symmetry([1, 2, 0], sector=1), ls.Symmetry([1, 2, 0], sector=2)])


def test_index():
    basis = ls.SpinBasis(4)
    basis.build()
    assert np.array_equal(basis.index(basis.states), basis.states)
    assert np.array_equal(basis.index(basis.states[2]), 2)


def test_kagome_symmetries():
    expr = ls.Expr(
        "1.0 σᶻ₀ σᶻ₁ + 1.0 σᶻ₀ σᶻ₃ + 1.0 σᶻ₀ σᶻ₈ + 1.0 σᶻ₀ σᶻ₁₀ + 2.0 σ⁺₀ σ⁻₁ + 2.0 σ⁺₀ σ⁻₃ + 2.0 σ⁺₀ σ⁻₈ + 2.0 σ⁺₀ σ⁻₁₀ + 2.0 σ⁻₀ σ⁺₁ + 2.0 σ⁻₀ σ⁺₃ + 2.0 σ⁻₀ σ⁺₈ + 2.0 σ⁻₀ σ⁺₁₀ + 1.0 σᶻ₁ σᶻ₂ + 0.8 σᶻ₁ σᶻ₃ + 0.8 σᶻ₁ σᶻ₉ + 2.0 σ⁺₁ σ⁻₂ + 1.6 σ⁺₁ σ⁻₃ + 1.6 σ⁺₁ σ⁻₉ + 2.0 σ⁻₁ σ⁺₂ + 1.6 σ⁻₁ σ⁺₃ + 1.6 σ⁻₁ σ⁺₉ + 1.0 σᶻ₂ σᶻ₄ + 1.0 σᶻ₂ σᶻ₉ + 1.0 σᶻ₂ σᶻ₁₀ + 2.0 σ⁺₂ σ⁻₄ + 2.0 σ⁺₂ σ⁻₉ + 2.0 σ⁺₂ σ⁻₁₀ + 2.0 σ⁻₂ σ⁺₄ + 2.0 σ⁻₂ σ⁺₉ + 2.0 σ⁻₂ σ⁺₁₀ + 1.0 σᶻ₃ σᶻ₅ + 0.8 σᶻ₃ σᶻ₁₁ + 2.0 σ⁺₃ σ⁻₅ + 1.6 σ⁺₃ σ⁻₁₁ + 2.0 σ⁻₃ σ⁺₅ + 1.6 σ⁻₃ σ⁺₁₁ + 0.8 σᶻ₄ σᶻ₆ + 1.0 σᶻ₄ σᶻ₇ + 0.8 σᶻ₄ σᶻ₁₀ + 1.6 σ⁺₄ σ⁻₆ + 2.0 σ⁺₄ σ⁻₇ + 1.6 σ⁺₄ σ⁻₁₀ + 1.6 σ⁻₄ σ⁺₆ + 2.0 σ⁻₄ σ⁺₇ + 1.6 σ⁻₄ σ⁺₁₀ + 1.0 σᶻ₅ σᶻ₆ + 1.0 σᶻ₅ σᶻ₈ + 1.0 σᶻ₅ σᶻ₁₁ + 2.0 σ⁺₅ σ⁻₆ + 2.0 σ⁺₅ σ⁻₈ + 2.0 σ⁺₅ σ⁻₁₁ + 2.0 σ⁻₅ σ⁺₆ + 2.0 σ⁻₅ σ⁺₈ + 2.0 σ⁻₅ σ⁺₁₁ + 1.0 σᶻ₆ σᶻ₇ + 0.8 σᶻ₆ σᶻ₈ + 2.0 σ⁺₆ σ⁻₇ + 1.6 σ⁺₆ σ⁻₈ + 2.0 σ⁻₆ σ⁺₇ + 1.6 σ⁻₆ σ⁺₈ + 1.0 σᶻ₇ σᶻ₉ + 1.0 σᶻ₇ σᶻ₁₁ + 2.0 σ⁺₇ σ⁻₉ + 2.0 σ⁺₇ σ⁻₁₁ + 2.0 σ⁻₇ σ⁺₉ + 2.0 σ⁻₇ σ⁺₁₁ + 0.8 σᶻ₈ σᶻ₁₀ + 1.6 σ⁺₈ σ⁻₁₀ + 1.6 σ⁻₈ σ⁺₁₀ + 0.8 σᶻ₉ σᶻ₁₁ + 1.6 σ⁺₉ σ⁻₁₁ + 1.6 σ⁻₉ σ⁺₁₁"
    )
    # top_shift = ls.Symmetry([5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 11, 10], sector=0)
    right_shift = ls.Symmetry([2, 10, 0, 4, 3, 7, 11, 5, 9, 8, 1, 6], sector=1)
    assert right_shift.periodicity == 2
    assert right_shift.phase == 0.5
    assert expr == expr.replace_indices(dict(zip(range(12), right_shift.permutation)))
    symmetries = ls.Symmetries([right_shift])
    basis = ls.SpinBasis(
        symmetries=symmetries, number_spins=12, hamming_weight=6, spin_inversion=None
    )
    basis.build()
    # print(basis.states)
    # print(basis.state_info(basis.states))
    hamiltonian = ls.Operator(basis, expr)
    assert basis.is_real
    assert expr.is_real
    assert hamiltonian.is_real
    energy, _ = scipy.sparse.linalg.eigsh(hamiltonian, k=1, which="SA")
    assert energy == approx(-19.95338528)


def test_operator_apply():
    basis = ls.SpinBasis(2)
    basis.build()
    expr = ls.Expr("1.0 σᶻ₀ σᶻ₁ + 2.0 σ⁺₀ σ⁻₁ + 2.0 σ⁻₀ σ⁺₁")
    hamiltonian = ls.Operator(basis, expr)
    assert len(hamiltonian.apply_off_diag_to_basis_state(basis.states[1])) == 1


def test_simple_spin_expr():
    sp = np.array([[0, 1], [0, 0]])
    sm = np.array([[0, 0], [1, 0]])
    sz = np.diag([1, -1])
    s0 = np.eye(2)

    def check(number_spins, expression, matrix_ref):
        basis = ls.SpinBasis(number_spins)
        basis.build()
        operator = ls.Operator(basis, ls.Expr(expression))
        cols = []
        for i in range(basis.number_states):
            v = np.zeros(basis.number_states)
            v[i] = 1
            cols.append(operator @ v)
        matrix = np.vstack(cols).T
        assert matrix.tolist() == matrix_ref.tolist()
        matrix = operator.to_csr().todense()
        assert matrix.tolist() == matrix_ref.tolist()

    check(1, "σ⁻₀", sm)
    check(1, "σ⁺₀", sp)
    check(1, "σᶻ₀", sz)
    check(2, "σ⁺₀ σ⁻₁", np.kron(sm, sp))
    check(2, "σ⁺₀ σᶻ₁", np.kron(sz, sp))
    check(2, "σ⁻₁", np.kron(sm, s0))
    check(3, "σ⁺₀ σᶻ₁ σ⁻₂", np.kron(sm, np.kron(sz, sp)))
    check(3, "σ⁺₀ σ⁻₂", np.kron(sm, np.kron(s0, sp)))

    check(1, "σʸ₀", -1j * sp + 1j * sm)
    check(3, "3im σ⁺₀ σᶻ₁ σʸ₂", 3j * np.kron((1j * sm - 1j * sp), np.kron(sz, sp)))


def test_prepare_hphi():
    basis = ls.SpinBasis(2)
    expr = ls.Expr("σᶻ₀ σᶻ₁ + 2 σ⁺₀ σ⁻₁ + 2 σ⁻₀ σ⁺₁")
    hamiltonian = ls.Operator(basis, expr)
    hamiltonian.prepare_inputs_for_hphi("/tmp/lattice-symmetries-python/hphi")


def test_prepare_mvmc():
    # basis = ls.SpinBasis(4)
    # expr = ls.Expr("σᶻ₀ σᶻ₁ + 2 σ⁺₀ σ⁻₁ + 2 σ⁻₀ σ⁺₁")
    # expr = sum((expr.replace_indices({0: i, 1: j}) for (i, j) in [(0, 1), (1, 2), (2, 3), (3, 0)]), ls.Expr(""))
    basis = ls.SpinfulFermionBasis(4, number_particles=4)

    hopping = ls.Expr("- (c†₁↑ c₀↑ + c†₀↑ c₁↑ + c†₁↓ c₀↓ + c†₀↓ c₁↓)")
    coulomb = ls.Expr("4.0 n₀↑ n₀↓")
    expr = hopping + coulomb
    for i, j in [(1, 2), (2, 3), (3, 0)]:
        expr += hopping.replace_indices({0: i, 1: j})
    for i in [1, 2, 3]:
        expr += coulomb.replace_indices({0: i})
    print(expr)
    hamiltonian = ls.Operator(basis, expr)
    hamiltonian.prepare_inputs_for_mvmc("/tmp/lattice-symmetries-python/mvmc")


def test_anisotropic_kagome_9():
    # fmt: off
    nearest = [
        (0, 1), (1, 2), (2, 0),
        (3, 4), (4, 5), (5, 3),
        (6, 7), (7, 8), (8, 6),
    ]
    next_nearest = [
        (0, 5), (0, 7),
        (1, 3), (1, 8),
        (2, 4), (2, 6),
        (3, 8),
        (4, 6),
        (5, 7),
    ]
    # fmt: on

    basis = ls.SpinfulFermionBasis(number_sites=9, number_particles=3)
    basis.build()
    hopping = lambda i, j: ls.Expr("c†₁↑ c₀↑ + c†₀↑ c₁↑ + c†₁↓ c₀↓ + c†₀↓ c₁↓").replace_indices(
        {0: i, 1: j}
    )
    coulomb = lambda i: ls.Expr("n₀↑ n₀↓").replace_indices({0: i})

    t1 = -0.3251
    t2 = 0.0845
    U = 2.8
    expr = (
        t1 * sum1((hopping(i, j) for i, j in nearest))
        + t2 * sum1((hopping(i, j) for i, j in next_nearest))
        + U * sum1((coulomb(i) for i in range(9)))
    )
    print(expr)
    hamiltonian = ls.Operator(basis, expr)
    hamiltonian.prepare_inputs_for_hphi("/tmp/kagome")
    energy, state = scipy.sparse.linalg.eigsh(hamiltonian, k=1, which="SA", tol=1e-6)
    print(energy)


def notest_vs_hphi():
    prefix = "../../test"
    for folder in os.listdir(prefix):
        print(folder)
        config = ls.load_yaml_config(os.path.join(prefix, folder, "hamiltonian.yaml"))
        config.basis.build()
        energy, state = scipy.sparse.linalg.eigsh(config.hamiltonian, k=1, which="SA", tol=1e-6)
        with open(os.path.join(prefix, folder, "HPhi", "output", "zvo_energy.dat")) as f:
            for line in f.readlines():
                if line.startswith("Energy"):
                    ref_energy = float(line.strip().split(" ")[-1])
        print(energy, ref_energy)
        assert ref_energy is not None
        assert energy == approx(ref_energy)


def test_apply_off_diag_projected():
    basis1 = ls.SpinBasis(4)
    basis1.build()
    expr = ls.Expr(
        "σᶻ₀ σᶻ₁ + 2 σ⁺₀ σ⁻₁ + 2 σ⁻₀ σ⁺₁",
        sites=[(0, 1), (1, 2), (2, 3), (3, 0)],
    )
    hamiltonian1 = ls.Operator(basis1, expr)
    assert hamiltonian1.apply_off_diag_to_basis_state(int("0101", base=2)) == [
        ((2 + 0j), 6),
        ((2 + 0j), 3),
        ((2 + 0j), 12),
        ((2 + 0j), 9),
    ]

    group = ls.Symmetries([ls.Symmetry([1, 2, 3, 0], sector=0)])
    basis = ls.SpinBasis(4, symmetries=group)
    basis.build()
    hamiltonian = ls.Operator(basis, expr)
    for c, x in hamiltonian.apply_off_diag_to_basis_state(int("0101", base=2)):
        assert x == 3
        assert c == approx(math.sqrt(2))

    group = ls.Symmetries([ls.Symmetry([1, 2, 3, 0], sector=2)])
    basis = ls.SpinBasis(4, symmetries=group)
    basis.build()
    hamiltonian = ls.Operator(basis, expr)
    assert hamiltonian.apply_off_diag_to_basis_state(int("0101", base=2)) == [
        (approx(-math.sqrt(2)), 3),
        (approx(math.sqrt(2)), 3),
        (approx(math.sqrt(2)), 3),
        (approx(-math.sqrt(2)), 3),
    ]


# def test_matvec():
#     if os.en


def get_csr_hamiltonian(hamiltonian):
    basis = hamiltonian.basis
    states, coeffs, row_idxs = hamiltonian.apply_off_diag_to_basis_state(basis.states)
    columns = basis.index(states)
    off_diagonal_matrix = scipy.sparse.csr_matrix(
        (coeffs, columns, row_idxs),
        shape=(basis.number_states, basis.number_states),
    )
    diagonal = scipy.sparse.diags(hamiltonian.apply_diag_to_basis_state(basis.states))
    return off_diagonal_matrix + diagonal


def test_correct_imag_part():
    expr = ls.Expr(
        # "σᶻ₀ σᶻ₁ +
        "2 σ⁺₀ σ⁻₁"
        # + 2 σ⁻₀ σ⁺₁"
        ,
        sites=[(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)],
    )
    basis = ls.SpinBasis(5, hamming_weight=2)
    basis.build()
    hamiltonian = ls.Operator(basis, expr)

    x = np.zeros(basis.number_states, dtype=complex)
    #                  43210
    x[basis.index(int("00011", base=2))] += 1.0
    x[basis.index(int("00110", base=2))] += np.exp(+2j * np.pi * 1 / 5)
    x[basis.index(int("01100", base=2))] += np.exp(+2j * np.pi * 2 / 5)
    x[basis.index(int("11000", base=2))] += np.exp(+2j * np.pi * 3 / 5)
    x[basis.index(int("10001", base=2))] += np.exp(+2j * np.pi * 4 / 5)
    x /= np.linalg.norm(x)

    y = np.zeros(basis.number_states, dtype=complex)
    y[basis.index(int("00101", base=2))] += 1.0
    y[basis.index(int("01010", base=2))] += np.exp(+2j * np.pi * 1 / 5)
    y[basis.index(int("10100", base=2))] += np.exp(+2j * np.pi * 2 / 5)
    y[basis.index(int("01001", base=2))] += np.exp(+2j * np.pi * 3 / 5)
    y[basis.index(int("10010", base=2))] += np.exp(+2j * np.pi * 4 / 5)
    y /= np.linalg.norm(y)

    matrix = np.array(
        [
            [np.vdot(x, hamiltonian @ x), np.vdot(x, hamiltonian @ y)],
            [np.vdot(y, hamiltonian @ x), np.vdot(y, hamiltonian @ y)],
        ]
    )

    group = ls.Symmetries([ls.Symmetry([4, 0, 1, 2, 3], sector=1)])
    basis = ls.SpinBasis(5, hamming_weight=2, symmetries=group)
    hamiltonian = ls.Operator(basis, expr)
    basis.build()

    cols = []
    for i in range(basis.number_states):
        v = np.zeros(basis.number_states)
        v[i] = 1
        cols.append(hamiltonian @ v)
    matrix2 = np.vstack(cols).T
    matrix3 = hamiltonian.to_csr().todense()
    print(matrix)
    print(matrix2)
    print(matrix3)
    assert np.allclose(matrix, matrix2)
    assert np.allclose(matrix, matrix3)


def test_convert_to_csr():
    def check(basis, expr):
        assert expr == expr.adjoint()
        hamiltonian = ls.Operator(basis, expr)
        basis.build()
        matrix = hamiltonian.to_csr()
        assert matrix.has_canonical_format
        assert (matrix != get_csr_hamiltonian(hamiltonian)).nnz == 0
        np.random.seed(42)
        for i in range(5):
            x = np.random.rand(basis.number_states)  # + 1j * np.random.rand(basis.number_states)
            y1 = hamiltonian @ x
            y2 = matrix @ x
            if not np.allclose(y1, y2):
                print(y1.tolist()[:10])
                print(y2.tolist()[:10])
            assert np.allclose(y1, y2)

    basis = ls.SpinBasis(2)
    expr = ls.Expr("σᶻ₀ σᶻ₁ + 2 σ⁺₀ σ⁻₁ + 2 σ⁻₀ σ⁺₁")
    check(basis, expr)

    group = ls.Symmetries([ls.Symmetry([1, 2, 3, 0], sector=0)])
    basis = ls.SpinBasis(4, symmetries=group)
    expr = ls.Expr(
        "σᶻ₀ σᶻ₁ + 2 σ⁺₀ σ⁻₁ + 2 σ⁻₀ σ⁺₁",
        sites=[(0, 1), (1, 2), (2, 3), (3, 0)],
    )
    check(basis, expr)

    basis = ls.Basis.from_json(
        '{"hamming_weight":null,"number_spins":9,"particle":"spin-1/2","spin_inversion":null,"symmetries":[{"permutation":[0,1,2,3,4,5,6,7,8],"sector":0},{"permutation":[1,2,0,4,5,3,7,8,6],"sector":0},{"permutation":[2,0,1,5,3,4,8,6,7],"sector":0},{"permutation":[3,4,5,6,7,8,0,1,2],"sector":0},{"permutation":[4,5,3,7,8,6,1,2,0],"sector":0},{"permutation":[5,3,4,8,6,7,2,0,1],"sector":0},{"permutation":[6,7,8,0,1,2,3,4,5],"sector":0},{"permutation":[7,8,6,1,2,0,4,5,3],"sector":0},{"permutation":[8,6,7,2,0,1,5,3,4],"sector":0}]}'
    )
    expr = ls.Expr(
        "-σᶻ₀ σᶻ₁ - σᶻ₀ σᶻ₂ - σᶻ₀ σᶻ₃ - σᶻ₀ σᶻ₆ - 2.0 σ⁺₀ σ⁻₁ - 2.0 σ⁺₀ σ⁻₂ - 2.0 σ⁺₀ σ⁻₃ - 2.0 σ⁺₀ σ⁻₆ - 2.0 σ⁻₀ σ⁺₁ - 2.0 σ⁻₀ σ⁺₂ - 2.0 σ⁻₀ σ⁺₃ - 2.0 σ⁻₀ σ⁺₆ - σᶻ₁ σᶻ₂ - σᶻ₁ σᶻ₄ - σᶻ₁ σᶻ₇ - 2.0 σ⁺₁ σ⁻₂ - 2.0 σ⁺₁ σ⁻₄ - 2.0 σ⁺₁ σ⁻₇ - 2.0 σ⁻₁ σ⁺₂ - 2.0 σ⁻₁ σ⁺₄ - 2.0 σ⁻₁ σ⁺₇ - σᶻ₂ σᶻ₅ - σᶻ₂ σᶻ₈ - 2.0 σ⁺₂ σ⁻₅ - 2.0 σ⁺₂ σ⁻₈ - 2.0 σ⁻₂ σ⁺₅ - 2.0 σ⁻₂ σ⁺₈ - σᶻ₃ σᶻ₄ - σᶻ₃ σᶻ₅ - σᶻ₃ σᶻ₆ - 2.0 σ⁺₃ σ⁻₄ - 2.0 σ⁺₃ σ⁻₅ - 2.0 σ⁺₃ σ⁻₆ - 2.0 σ⁻₃ σ⁺₄ - 2.0 σ⁻₃ σ⁺₅ - 2.0 σ⁻₃ σ⁺₆ - σᶻ₄ σᶻ₅ - σᶻ₄ σᶻ₇ - 2.0 σ⁺₄ σ⁻₅ - 2.0 σ⁺₄ σ⁻₇ - 2.0 σ⁻₄ σ⁺₅ - 2.0 σ⁻₄ σ⁺₇ - σᶻ₅ σᶻ₈ - 2.0 σ⁺₅ σ⁻₈ - 2.0 σ⁻₅ σ⁺₈ - σᶻ₆ σᶻ₇ - σᶻ₆ σᶻ₈ - 2.0 σ⁺₆ σ⁻₇ - 2.0 σ⁺₆ σ⁻₈ - 2.0 σ⁻₆ σ⁺₇ - 2.0 σ⁻₆ σ⁺₈ - σᶻ₇ σᶻ₈ - 2.0 σ⁺₇ σ⁻₈ - 2.0 σ⁻₇ σ⁺₈"
    )
    check(basis, expr)

    basis = ls.Basis.from_json(
        '{"hamming_weight":null,"number_spins":14,"particle":"spin-1/2","spin_inversion":null,"symmetries":[{"permutation":[0,1,2,3,4,5,6,7,8,9,10,11,12,13],"sector":0},{"permutation":[1,2,3,4,5,6,7,8,9,10,11,12,13,0],"sector":11},{"permutation":[2,3,4,5,6,7,8,9,10,11,12,13,0,1],"sector":4},{"permutation":[3,4,5,6,7,8,9,10,11,12,13,0,1,2],"sector":5},{"permutation":[4,5,6,7,8,9,10,11,12,13,0,1,2,3],"sector":1},{"permutation":[5,6,7,8,9,10,11,12,13,0,1,2,3,4],"sector":13},{"permutation":[6,7,8,9,10,11,12,13,0,1,2,3,4,5],"sector":5},{"permutation":[7,8,9,10,11,12,13,0,1,2,3,4,5,6],"sector":1},{"permutation":[8,9,10,11,12,13,0,1,2,3,4,5,6,7],"sector":2},{"permutation":[9,10,11,12,13,0,1,2,3,4,5,6,7,8],"sector":1},{"permutation":[10,11,12,13,0,1,2,3,4,5,6,7,8,9],"sector":6},{"permutation":[11,12,13,0,1,2,3,4,5,6,7,8,9,10],"sector":9},{"permutation":[12,13,0,1,2,3,4,5,6,7,8,9,10,11],"sector":3},{"permutation":[13,0,1,2,3,4,5,6,7,8,9,10,11,12],"sector":3}]}'
    )
    expr = ls.Expr(
        "-σᶻ₀ σᶻ₁ - σᶻ₀ σᶻ₁₃ - 2.0 σ⁺₀ σ⁻₁ - 2.0 σ⁺₀ σ⁻₁₃ - 2.0 σ⁻₀ σ⁺₁ - 2.0 σ⁻₀ σ⁺₁₃ - σᶻ₁ σᶻ₂ - 2.0 σ⁺₁ σ⁻₂ - 2.0 σ⁻₁ σ⁺₂ - σᶻ₂ σᶻ₃ - 2.0 σ⁺₂ σ⁻₃ - 2.0 σ⁻₂ σ⁺₃ - σᶻ₃ σᶻ₄ - 2.0 σ⁺₃ σ⁻₄ - 2.0 σ⁻₃ σ⁺₄ - σᶻ₄ σᶻ₅ - 2.0 σ⁺₄ σ⁻₅ - 2.0 σ⁻₄ σ⁺₅ - σᶻ₅ σᶻ₆ - 2.0 σ⁺₅ σ⁻₆ - 2.0 σ⁻₅ σ⁺₆ - σᶻ₆ σᶻ₇ - 2.0 σ⁺₆ σ⁻₇ - 2.0 σ⁻₆ σ⁺₇ - σᶻ₇ σᶻ₈ - 2.0 σ⁺₇ σ⁻₈ - 2.0 σ⁻₇ σ⁺₈ - σᶻ₈ σᶻ₉ - 2.0 σ⁺₈ σ⁻₉ - 2.0 σ⁻₈ σ⁺₉ - σᶻ₉ σᶻ₁₀ - 2.0 σ⁺₉ σ⁻₁₀ - 2.0 σ⁻₉ σ⁺₁₀ - σᶻ₁₀ σᶻ₁₁ - 2.0 σ⁺₁₀ σ⁻₁₁ - 2.0 σ⁻₁₀ σ⁺₁₁ - σᶻ₁₁ σᶻ₁₂ - 2.0 σ⁺₁₁ σ⁻₁₂ - 2.0 σ⁻₁₁ σ⁺₁₂ - σᶻ₁₂ σᶻ₁₃ - 2.0 σ⁺₁₂ σ⁻₁₃ - 2.0 σ⁻₁₂ σ⁺₁₃"
    )
    check(basis, expr)


def test_csr_matvec():
    basis = ls.SpinBasis(2)
    basis.build()
    expr = ls.Expr("σᶻ₀ σᶻ₁ + 2 σ⁺₀ σ⁻₁ + 2 σ⁻₀ σ⁺₁")
    hamiltonian = ls.Operator(basis, expr)
    matrix = hamiltonian.to_csr()
    x = np.random.rand(basis.number_states).astype(np.complex128)
    assert np.allclose(hamiltonian @ x, ls.matrix_vector_product_csr(matrix, x))


def test_abelian_representations():
    basis = ls.SpinBasis(3)
    basis.build()
    expr = ls.Expr("σᶻ₀ σᶻ₁ + 2 σ⁺₀ σ⁻₁ + 2 σ⁻₀ σ⁺₁", sites=[[0, 1], [1, 2], [2, 0]])
    hamiltonian = ls.Operator(basis, expr)
    # print(hamiltonian.to_csr().todense().real)
    for r in hamiltonian.abelian_representations():
        print(r)


def test_ground_state_in_abelian_representations():
    basis = ls.SpinBasis(10)
    basis.build()
    expr = ls.Expr(
        "σᶻ₀ σᶻ₁ + 2 σ⁺₀ σ⁻₁ + 2 σ⁻₀ σ⁺₁",
        sites=[(i, (i + 1) % basis.number_bits) for i in range(basis.number_bits)],
    )
    hamiltonian = ls.Operator(basis, expr)
    real_energy, _ = scipy.sparse.linalg.eigsh(hamiltonian, k=1, which="SA")

    energies = []
    for r in hamiltonian.abelian_representations():
        symm_basis = ls.SpinBasis(basis.number_bits, symmetries=r)
        symm_basis.build()
        symm_hamiltonian = ls.Operator(symm_basis, expr)
        e, _ = scipy.sparse.linalg.eigsh(symm_hamiltonian, k=1, which="SA")
        energies.append(e)

    assert np.any(np.isclose(real_energy, energies))
