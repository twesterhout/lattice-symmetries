import time, numpy as np, igraph as ig, lattice_symmetries as ls
import pytest, hypothesis, hypothesis.strategies as st, hypothesis.extra.numpy
from sympy import S, simplify, Rational
from sympy.combinatorics import Permutation

phases = (hypothesis.Phase.explicit, hypothesis.Phase.reuse, hypothesis.Phase.generate)
rng = np.random.default_rng(seed=123)

@pytest.mark.parametrize("number_bits", list(range(2, 10)))
def test_just_id(number_bits):
    b = ls.Basis(ls.BasisInfo(bits=number_bits, symmetries=[(Permutation(np.arange(number_bits)), Rational(0))]))
    xs = np.array([0, 1, 2, 3, 48, 100, 143], dtype=np.uint64)
    rep, idx = b.state_info(xs)
    assert rep.tolist() == xs.tolist()
    assert idx.tolist() == [0] * xs.size

def test_ring3():
    p, k = Permutation(np.roll(np.arange(3), shift=-1)), 0
    b = ls.Basis(ls.BasisInfo(bits=3, symmetries=[(p, Rational(k, p.order()))]))
    xs = np.arange(8, dtype=np.uint64)
    rep, idx = b.state_info(xs)
    assert list(zip(rep.tolist(), idx.tolist())) == [
        (0, 0), # 000
        (1, 0), # 001
        (1, 1), # 010
        (3, 0), # 011
        (1, 2), # 100
        (3, 2), # 101
        (3, 1), # 110
        (7, 0), # 111
    ]
