import dataclasses
import math
import itertools
from copy import deepcopy
from functools import reduce

import more_itertools
from sympy.combinatorics import Permutation


@dataclasses.dataclass(frozen=True)
class Edge:
    focus: str
    left: tuple[int, int]
    right: tuple[int, int]


@dataclasses.dataclass(frozen=True)
class Swap:
    focus: str
    left: int
    right: int


class InvertiblePermutation:
    permutation: any
    inverse: any

    def __init__(self, permutation):
        self.permutation = permutation
        self.inverse = ~self.permutation

    def index(self, v) -> int:
        return int(self.inverse.apply(v))

    def value(self, i) -> int:
        return int(self.permutation.apply(i))


def _bit_permute_step(x, mask, delta):
    y = ((x >> delta) ^ x) & mask
    return (x ^ y) ^ (y << delta)


@dataclasses.dataclass(frozen=True)
class BenesNetwork:
    masks: list[int]
    shifts: list[int]

    def __call__(self, bits: int):
        assert isinstance(bits, int)
        for m, d in zip(self.masks, self.shifts):
            bits = _bit_permute_step(bits, m, d)
        return bits


def _get_cycle(
    d: int, src: InvertiblePermutation, tgt: InvertiblePermutation, edge: Edge
) -> list[Edge]:
    is_smaller = lambda i: i % (2 * d) < d  # noqa
    get_neighbor = lambda i: i + d if is_smaller(i) else i - d  # noqa

    edge0 = edge
    loc = edge0.focus
    x = edge0.right[1]
    cycle = [edge0]
    while True:
        loc = "target" if loc == "source" else "source"
        p = dict(source=src, target=tgt)[loc]
        i = p.index(x)
        xi = p.value(i)
        j = get_neighbor(i)
        xj = p.value(j)

        if is_smaller(i):
            edge = Edge(focus=loc, left=(i, xi), right=(j, xj))
        else:
            edge = Edge(focus=loc, left=(j, xj), right=(i, xi))

        if edge == edge0:
            cycle.append(edge0)
            return cycle
        cycle.append(edge)
        x = {xi: xj, xj: xi}[x]


def _solve_cycle(edges: list[Edge]) -> list[Swap]:
    edges0 = deepcopy(edges)
    edges = deepcopy(edges)

    def should_swap(edge1, edge2):
        if edge1.left[1] == edge2.left[1] or edge1.right[1] == edge2.right[1]:
            return False
        if edge1.left[1] == edge2.right[1] or edge1.right[1] == edge2.left[1]:
            return True
        assert False, f"invalid cycle {edges0}"

    swaps = []
    for k in range(len(edges) - 1):
        e0 = edges[k]
        e1 = edges[k + 1]
        if should_swap(e0, e1):
            swaps.append(Swap(e1.focus, e1.left[0], e1.right[0]))
            edges[k + 1] = Edge(e1.focus, e1.right, e1.left)
    if len(edges) > 0:
        e0 = edges[-1]
        assert e0.left[0] < e0.right[0], f"unsolvable cycle {edges0}"

    return swaps[::-1]


def _initial_edges(d, src):
    n = src.permutation.size
    mk_edge = lambda a, b: Edge("source", (a, src.value(a)), (b, src.value(b)))  # noqa
    return [mk_edge(i + j, i + j + d) for i in range(0, n - d + 1, 2 * d) for j in range(d)]


def _stage_cycles(d, src, tgt):
    visited = set()
    cycles = []
    edges = _initial_edges(d, src)
    for e in edges:
        if e not in visited:
            cycle = _get_cycle(d, src, tgt, e)
            cycles.append(cycle)
            visited |= set(cycle)
    return cycles


def _apply_swaps(p, swaps):
    return reduce(lambda acc, s: Permutation(s.left, s.right) * acc, swaps, p)


def _get_mask(swaps):
    return reduce(lambda acc, s: acc | (1 << min(s.left, s.right)), swaps, 0)


def _solve_stage(d, src, tgt):
    swaps = itertools.chain(*(_solve_cycle(cycle) for cycle in _stage_cycles(d, src, tgt)))
    tgt_swaps, src_swaps = more_itertools.partition(lambda s: s.focus == "source", swaps)
    src_swaps = list(src_swaps)
    tgt_swaps = list(tgt_swaps)
    src = InvertiblePermutation(_apply_swaps(src.permutation, src_swaps))
    tgt = InvertiblePermutation(_apply_swaps(tgt.permutation, tgt_swaps))
    return _get_mask(src_swaps), _get_mask(tgt_swaps), (src, tgt)


def _solve(src, tgt):
    n = src.permutation.size
    d = 1
    stages = []
    while 2 * d <= n:
        (src_mask, tgt_mask, (src, tgt)) = _solve_stage(d, src, tgt)
        stages.append((src_mask, tgt_mask, d))
        d *= 2

    if len(stages) == 0:
        return [], []
    else:
        src_masks, tgt_masks, shifts = tuple(map(list, zip(*stages)))
        assert src_masks[-1] == 0, f"invalid masks: {src_masks} and {tgt_masks}"
        return src_masks[:-1] + tgt_masks[::-1], shifts[:-1] + shifts[::-1]


def permutation_to_benes_network(p: Permutation) -> BenesNetwork:
    if p.size == 0:
        return BenesNetwork(masks=[], shifts=[])
    n = 2 ** math.ceil(math.log2(p.size))
    p = p.resize(n)
    solution = _solve(InvertiblePermutation(Permutation(list(range(n)))), InvertiblePermutation(p))
    return BenesNetwork(*solution)
