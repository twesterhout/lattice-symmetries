import time, numpy as np, igraph as ig, lattice_symmetries as ls
import pytest, hypothesis, hypothesis.strategies as st, hypothesis.extra.numpy
from sympy import S, simplify, Rational
from sympy.combinatorics import Permutation

phases = (hypothesis.Phase.explicit,) # hypothesis.Phase.reuse, hypothesis.Phase.generate)
rng = np.random.default_rng(seed=123)
kernels = ls.compiler.build_kernels()

def reverse_bits(x, n):
    x = np.asarray(x).copy(); r = np.zeros_like(x)
    for _ in range(n): r <<= 1; r |= x & 1; x >>= 1
    return r

@st.composite
def rnd_symmetries(draw, n: int):
    p = Permutation(draw(st.permutations(list(range(n)))))
    k = st.integers(min_value=0, max_value=p.order() - 1)
    return [(p, Rational(k, p.order()))]

@pytest.mark.parametrize("number_bits", list(range(1, 4)))
def test_no_symmetries(number_bits):
    import quspin
    b = ls.Basis(ls.BasisInfo(number_bits))
    b.build()
    quspin_basis = quspin.basis.spin_basis_general(b.i.bits)
    states_ref, norms_ref = quspin_basis.states, quspin_basis.normalization(quspin_basis.states)
    states_ref = reverse_bits(states_ref, b.i.bits)
    order = np.argsort(states_ref); states_ref, norms_ref = states_ref[order], norms_ref[order]
    assert b.states.tolist() == states_ref.tolist()
    assert b.norms.tolist() == norms_ref.tolist()

@pytest.mark.parametrize("number_bits", list(range(1, 10)))
def test_just_id(number_bits):
    b_ref = ls.Basis(ls.BasisInfo(number_bits)); b_ref.build()
    i = ls.BasisInfo(bits=number_bits, symmetries=[(Permutation(np.arange(number_bits)), Rational(0))])
    b = ls.Basis(i); b.build()
    assert b.states.tolist() == b_ref.states.tolist()
    assert b.norms.tolist() == b_ref.norms.tolist()

@pytest.mark.parametrize("number_bits", list(range(2, 10)))
def test_ring_cyclic(number_bits):
    import quspin
    p, k = Permutation(np.roll(np.arange(number_bits), shift=-1)), 0
    i = ls.BasisInfo(bits=number_bits, symmetries=[(p, Rational(k, p.order()))])
    b = ls.Basis(i); b.build()
    group_size = len(i.symmetries)
    quspin_basis = quspin.basis.spin_basis_general(i.bits, pxblock=(p.array_form, k))
    states_ref, norms_ref = quspin_basis.states, quspin_basis.normalization(quspin_basis.states)
    states_ref, norms_ref = states_ref ^ (2**i.bits - 1), norms_ref // group_size
    order = np.argsort(states_ref); states_ref, norms_ref = states_ref[order], norms_ref[order]
    assert b.states.tolist() == states_ref.tolist()
    assert b.norms.tolist() == norms_ref.tolist()

@hypothesis.given(st.integers(min_value=2, max_value=10), st.integers(min_value=0, max_value=10))
@hypothesis.example(3, 1)
@hypothesis.settings(deadline=None, phases=phases)
def test_ring_cyclic_nontrivial(number_bits, k):
    import quspin
    p = Permutation(np.roll(np.arange(number_bits), shift=-1))
    k = k % number_bits
    i = ls.BasisInfo(bits=number_bits, symmetries=[(p, Rational(k, p.order()))])
    b = ls.Basis(i); b.build()
    group_size = len(i.symmetries)
    quspin_basis = quspin.basis.spin_basis_general(i.bits, pxblock=(p.array_form, k))
    states_ref, norms_ref = quspin_basis.states, quspin_basis.normalization(quspin_basis.states)
    states_ref, norms_ref = states_ref ^ (2**i.bits - 1), norms_ref // group_size
    order = np.argsort(states_ref); states_ref, norms_ref = states_ref[order], norms_ref[order]
    assert b.states.tolist() == states_ref.tolist()
    assert b.norms.tolist() == norms_ref.tolist()
