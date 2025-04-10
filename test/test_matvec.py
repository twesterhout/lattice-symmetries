import time, numpy as np, igraph as ig, sympy, lattice_symmetries as ls
import pytest, hypothesis, hypothesis.strategies as st, hypothesis.extra.numpy
from dataclasses import dataclass
from sympy import S, simplify, Rational
from sympy.combinatorics import Permutation

phases = (hypothesis.Phase.explicit, hypothesis.Phase.reuse, hypothesis.Phase.generate)
rng = np.random.default_rng(seed=123)

def _off_diag(m): return m - np.diag(np.diag(m))

def test_ring3_no_symm():
    rng = np.random.default_rng(5)
    b = ls.B(ls.BasisInfo(bits=3))
    e = ls.heisenberg(ig.Graph.Ring(3, circular=True), h=0.125)
    o = ls.O(e, b)

    b.build()
    x, n = rng.random(128, dtype=np.float64).view(np.complex128), b.number_states

    ref = _off_diag(e.to_dense()) @ x[:n]
    out = o._prepare_Matvec()._off_diag64(b.states, b.norms, x)
    np.testing.assert_allclose(out[:n], ref)


def test_ring3_symm():
    rng = np.random.default_rng(5)
    b = ls.B(ls.BasisInfo(bits=3, symmetries=[(Permutation([1, 2, 0]), Rational(0))]))
    e = ls.heisenberg(ig.Graph.Ring(3, circular=True), h=0.125)

    b.build()
    assert b.states.tolist() == [0, 1, 3, 7]
    assert b.norms.tolist() == [3, 1, 1, 3]
    x = rng.random(128, dtype=np.float64).view(np.complex128)
    matvec = ls.O(e, b)._prepare_Matvec()

    out = matvec._off_diag64(b.states, b.norms, x.real)
    assert out[0] == pytest.approx(0.0)
    assert out[1] == pytest.approx(4 * x.real[1])
    assert out[2] == pytest.approx(4 * x.real[2])
    assert out[3] == pytest.approx(0.0)

    out = matvec._off_diag64(b.states, b.norms, x)
    assert out[0] == pytest.approx(0.0)
    assert out[1] == pytest.approx(4 * x[1])
    assert out[2] == pytest.approx(4 * x[2])
    assert out[3] == pytest.approx(0.0)

    out = matvec._diag64(b.states, x.real)
    assert out[0] == pytest.approx(x.real[0] * (3 - 3 * 0.125))
    assert out[1] == pytest.approx(x.real[1] * (-1 - 0.125))
    assert out[2] == pytest.approx(x.real[2] * (-1 + 0.125))
    assert out[3] == pytest.approx(x.real[3] * (3 + 3 * 0.125))

    out = matvec._diag64(b.states, x)
    assert out[0] == pytest.approx(x[0] * (3 - 3 * 0.125))
    assert out[1] == pytest.approx(x[1] * (-1 - 0.125))
    assert out[2] == pytest.approx(x[2] * (-1 + 0.125))
    assert out[3] == pytest.approx(x[3] * (3 + 3 * 0.125))


def test_ring4():
    rng = np.random.default_rng(5)
    b = ls.B(ls.BasisInfo(bits=4, symmetries=[(Permutation([1, 2, 3, 0]), Rational(1, 2))]))
    e = ls.heisenberg(ig.Graph.Ring(4, circular=True), h=0.125)

    b.build()
    assert b.states.tolist() == [1, 3, 5, 7]
    assert b.norms.tolist() == [1, 1, 2, 1]

    matvec = ls.O(e, b)._prepare_Matvec()
    x = rng.random(128, dtype=np.float64).view(np.complex128)

    out = matvec._off_diag64(b.states, b.norms, x.real)
    assert out[0] == pytest.approx(-4 * x.real[0])
    assert out[1] == pytest.approx(0.0)
    assert out[2] == pytest.approx(0.0)
    assert out[3] == pytest.approx(-4 * x.real[3])

    out = matvec._diag64(b.states, x.real)
    assert out[0] == pytest.approx(x.real[0] * (0 - 2 * 0.125))
    assert out[1] == pytest.approx(x.real[1] * 0)
    assert out[2] == pytest.approx(x.real[2] * (-4))
    assert out[3] == pytest.approx(x.real[3] * (0 + 2 * 0.125))


def quspin_heisenberg(info, graph, symmetries=None):
    import quspin.basis, quspin.operators
    if symmetries is None: symmetries = info.symmetries
    blocks = {
        chr(ord("a") + i) + "block": (p.array_form, int(k * p.order()))
        for i, (p, k) in enumerate(symmetries)
    }
    # print(blocks)
    basis = quspin.basis.spin_basis_general(info.bits, pauli=-1, **blocks)
    nearest = [e.tuple for e in graph.es]
    static = [
        ["+-", [[2, i, j] for (i, j) in nearest]],
        ["-+", [[2, i, j] for (i, j) in nearest]],
        ["zz", [[1, i, j] for (i, j) in nearest]],
        ["z", [[-0.125, i] for i in range(info.bits)]]
    ]
    # QuSpin doesn't detect that the Hamiltonian is real if some characters are -1...
    hamiltonian = quspin.operators.hamiltonian(static, [], basis=basis, dtype=np.complex128, check_symm=False, check_herm=False)
    return basis, hamiltonian

@pytest.mark.parametrize("number_bits,k", [(3, 0), (3, 1), (3, 2), (4, 0), (4, 1), (4, 2), (5, 0), (6, 0), (6, 1), (6, 2), (6, 3), (10, 5)])
def test_ringX_translation(number_bits, k):
    import quspin, quspin.operators
    p, k = Permutation(np.roll(np.arange(number_bits), shift=-1)), k % number_bits
    symmetries = [(p, Rational(k, p.order()))]
    b = ls.Basis(ls.BasisInfo(bits=number_bits, symmetries=symmetries))
    g = ig.Graph.Ring(number_bits, circular=True)
    h = ls.heisenberg(g, h=0.125)

    def invert(p): inv = np.empty_like(p); inv[p] = np.arange(len(inv), dtype=inv.dtype); return inv
    quspin_basis, quspin_hamiltonian = quspin_heisenberg(b.i, g, symmetries=symmetries)
    order = np.argsort(quspin_basis.states ^ (2**b.i.bits - 1))
    # print((quspin_basis.states ^ (2**i.bits - 1))[order])
    # print(quspin_basis.normalization(quspin_basis.states)[order])

    b.build()
    o = ls.O(h, b)

    rng = np.random.default_rng(5)
    if k == 0 or 2 * k == number_bits: # characters are real
        x = rng.random(b.states.size, dtype=np.float64)
        out = o @ x
        out_ref = quspin_hamiltonian.dot(x[invert(order)])[order].real
        np.testing.assert_allclose(out, out_ref, rtol=1e-8, atol=1e-10)

    x = rng.random(2 * b.states.size, dtype=np.float64).view(np.complex128)
    out = o @ x
    out_ref = quspin_hamiltonian.dot(x[invert(order)])[order]
    np.testing.assert_allclose(out, out_ref, rtol=1e-8, atol=1e-10)

@pytest.mark.parametrize("number_bits,k", [(3, 0), (3, 1), (4, 0), (4, 1), (5, 0), (6, 0), (6, 1), (10, 0), (10, 1)])
def test_ringX_reflection(number_bits, k):
    import quspin, quspin.operators
    p, k = Permutation(np.arange(number_bits)[::-1]), k % number_bits
    symmetries = [(p, Rational(k, p.order()))]
    b = ls.Basis(ls.BasisInfo(bits=number_bits, symmetries=symmetries))
    g = ig.Graph.Ring(number_bits, circular=True)
    h = ls.heisenberg(g, h=0.125)

    def invert(p): inv = np.empty_like(p); inv[p] = np.arange(len(inv), dtype=inv.dtype); return inv
    quspin_basis, quspin_hamiltonian = quspin_heisenberg(b.i, g, symmetries=symmetries)
    order = np.argsort(quspin_basis.states ^ (2**b.i.bits - 1))

    b.build(); o = ls.O(h, b)

    rng = np.random.default_rng(5)
    # real
    x = rng.random(b.states.size, dtype=np.float64)
    out = o @ x
    out_ref = quspin_hamiltonian.dot(x[invert(order)])[order].real
    np.testing.assert_allclose(out, out_ref)
    # complex
    x = rng.random(2 * b.states.size, dtype=np.float64).view(np.complex128)
    out = o @ x
    out_ref = quspin_hamiltonian.dot(x[invert(order)])[order]
    np.testing.assert_allclose(out, out_ref, rtol=1e-8, atol=1e-10)

@pytest.mark.parametrize("number_bits,k1,k2", [(3, 0, 0), (4, 0, 0), (4, 2, 1), (8, 0, 0), (8, 4, 1), (8, 0, 1)])
def test_ringX_both(number_bits, k1, k2):
    p1 = Permutation(np.roll(np.arange(number_bits), shift=-1))
    p2 = Permutation(np.arange(number_bits)[::-1])
    symmetries = [(p1, Rational(k1, p1.order())), (p2, Rational(k2, p2.order()))]
    b = ls.Basis(ls.BasisInfo(bits=number_bits, symmetries=symmetries))
    g = ig.Graph.Ring(number_bits, circular=True)
    h = ls.heisenberg(g, h=0.125)

    def invert(p): inv = np.empty_like(p); inv[p] = np.arange(len(inv), dtype=inv.dtype); return inv
    quspin_basis, quspin_hamiltonian = quspin_heisenberg(b.i, g, symmetries=symmetries)
    order = np.argsort(quspin_basis.states ^ (2**b.i.bits - 1))

    b.build(); o = ls.O(h, b)
    assert (quspin_basis.states ^ (2**b.i.bits - 1))[order].tolist() == b.states.tolist()
    
    rng = np.random.default_rng(5)
    # real
    x = rng.random(b.states.size, dtype=np.float64)
    out = o @ x
    out_ref = quspin_hamiltonian.dot(x[invert(order)])[order].real
    np.testing.assert_allclose(out, out_ref)
    # complex
    x = rng.random(2 * b.states.size, dtype=np.float64).view(np.complex128)
    out = o @ x
    out_ref = quspin_hamiltonian.dot(x[invert(order)])[order]
    np.testing.assert_allclose(out, out_ref, rtol=1e-8, atol=1e-10)


@st.composite
def rnd_pauli_term(draw, bits, max_order=3):
    import operator
    from functools import reduce
    p = st.sampled_from([ls.SigmaX, ls.SigmaY, ls.SigmaZ, ls.SigmaPlus, ls.SigmaMinus])
    i = st.integers(min_value=0, max_value=bits - 1)
    c = draw(st.complex_numbers(min_magnitude=1e-3, max_magnitude=10, allow_subnormal=False))
    t = draw(st.lists(st.builds(lambda a, b: a(b), p, i), min_size=0, max_size=max_order))
    return c * reduce(operator.mul, t, S.One)

@st.composite
def rnd_pauli_expr(draw, max_bits=8, max_terms=10, max_order=5):
    bits = draw(st.integers(min_value=1, max_value=max_bits - 1))
    elements = draw(st.lists(rnd_pauli_term(bits=bits, max_order=max_order), min_size=1,
        max_size=max_terms))
    return sum(elements)

@hypothesis.given(rnd_pauli_expr(max_bits=5))
@hypothesis.example(sympy.Float(0.001))
@hypothesis.example(-9.999 - 2.2250738585072e-309 * sympy.I - 10.0 * sympy.I * ls.SigmaX(1))
@hypothesis.example(ls.SigmaPlus(0))
@hypothesis.example((0.5 + 5.96046447753906e-8 * sympy.I) * ls.SigmaPlus(1) + (-1.40129846432482e-45 - 5.80934524614774 * sympy.I) * ls.SigmaZ(2))
@hypothesis.example((-0.32052911031739 + 9.99486173438328 * sympy.I) * ls.SigmaPlus(0) * ls.SigmaPlus(1) + 10.0 * ls.SigmaX(0))
@hypothesis.settings(max_examples=20, deadline=None, phases=phases)
def test_pauli_small(e):
    rng = np.random.default_rng(123)
    o = ls.O(ls.Expr(e))
    o.b.build()

    x = rng.uniform(0, 1, size=o.b.states.size)
    ref = o.e.to_dense() @ x
    out = o @ x
    np.testing.assert_allclose(out, ref, rtol=1e-8, atol=1e-10)

