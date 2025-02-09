import lattice_symmetries as ls
import importlib, hypothesis, numpy as np, random, sympy, pytest
from sympy.combinatorics import Permutation, PermutationGroup
from pytest import raises, approx
import hypothesis.strategies as st

our_phases = (hypothesis.Phase.explicit, hypothesis.Phase.reuse, hypothesis.Phase.generate)

def test_permutation_to_benes_network_examples():
    benes = ls.perm2benes(Permutation([0, 1, 2, 3]))
    assert benes(0b0100) == 0b0100
    benes = ls.perm2benes(Permutation([1, 2, 0, 3]))
    assert benes(0b0100) == 0b0010
    assert benes(0b0101) == 0b0110
    benes = ls.perm2benes(Permutation([1, 2, 3, 0, 5, 6, 7, 4, 9, 10, 11, 8]))
    assert benes(0b100111010010) == 0b110011100001
    benes = ls.perm2benes(Permutation([0]))
    assert benes(0b0) == 0b0
    assert benes(0b1) == 0b1
    benes = ls.perm2benes(Permutation([]))
    assert benes(0b0) == 0b0
    assert benes(0b10100110) == 0b10100110

@st.composite
def permutations(draw, min_size=1, max_size=1000):
    size = draw(st.integers(min_value=min_size, max_value=max_size))
    arr = draw(st.permutations(list(range(size))))
    return Permutation(arr)

def reference_permute_bits(p: Permutation, bits: int):
    assert bits < 2**p.size
    s = ("{:0" + str(p.size) + "b}").format(bits)[::-1]
    s = "".join(s[int(i)] for i in p.array_form)[::-1]
    return int(s, base=2)

@hypothesis.given(permutations())
@hypothesis.settings(max_examples=10, deadline=None, phases=our_phases)
def test_permutation_to_benes_network(permutation):
    random.seed(42)
    benes = ls.perm2benes(permutation)
    for i in range(10):
        bits = random.randint(0, 2**permutation.size - 1)
        assert benes(bits) == reference_permute_bits(permutation, bits)
