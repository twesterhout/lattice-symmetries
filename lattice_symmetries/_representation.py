import itertools
from sympy import Rational
from sympy.combinatorics import Permutation, PermutationGroup


def generate_representation(symmetries):
    symmetries = list(symmetries)
    if len(symmetries) == 0:
        return []
    generators = [g for g, _ in symmetries]
    representation = dict()

    def add(p, r):
        if p not in representation:
            representation[p] = r % 1
        elif not (r % 1).equals(representation[p]):
            msg = f"conflicting phase factors for the permutation {p.array_form}: {r} and {representation[p]}"
            raise ValueError(msg)

    for g, r in symmetries:
        add(g, r)

    interior = set()
    boundary = set(generators)
    while len(boundary) > 0:
        interior |= boundary
        next_boundary = set()
        for h, g in itertools.product(boundary, generators):
            p = h * g
            r = (representation[h] + representation[g]) % 1
            add(p, r)
            if p not in interior:
                next_boundary.add(p)
        boundary = next_boundary

    for g in PermutationGroup(generators).elements:
        if g not in representation:
            raise ValueError(f"missing the representation of {g.array_form}")
    return sorted(list(representation.items()), key=lambda x: x[0].array_form)
