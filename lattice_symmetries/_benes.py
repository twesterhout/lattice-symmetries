import dataclasses, math, numbers, numpy as np
from itertools import islice, accumulate
from more_itertools import iterate
from functools import reduce, lru_cache
from numpy.typing import NDArray
from sympy.combinatorics import Permutation, PermutationGroup
from typing import Callable

def _bit_permute_step(x, mask, delta):
    y = ((x >> delta) ^ x) & mask
    return (x ^ y) ^ (y << delta)

@dataclasses.dataclass(frozen=True)
class BenesNetwork:
    masks: list[int]
    shifts: list[int]
    swaps: list[NDArray] | None = None
    def __call__(self, bits: int | NDArray) -> int | NDArray:
        bits = int(bits) if isinstance(bits, numbers.Integral) else np.asarray(bits)
        return reduce(lambda acc, t: _bit_permute_step(acc, *t), zip(self.masks, self.shifts), bits)

def _stage(d, n, r):
    smaller = np.arange(n) % (2 * d) < d; neighbor = np.arange(n) + np.where(smaller, d, -d)
    # building the cycles
    def f(acc): loc, v = acc; return loc ^ 1, r[loc, 0, neighbor[r[loc, 1, v]]]
    es0 = (np.arange(0, n - d + 1, 2 * d).reshape(-1, 1) + np.arange(d)).reshape(-1, 1)
    i = r[np.arange(2 * n) % 2, 1, np.hstack([v for _, v in islice(iterate(f, (0, r[0, 0, es0])), 2 * n)])]
    si = i[:, ::2]; ni = np.where(~smaller[si], neighbor[si], si)
    # determine unique cycles
    seen = np.zeros(ni.shape, dtype=np.bool_); seen[np.ogrid[:ni.shape[0], :ni.shape[1]][0], ni] = 1
    j = np.unique(seen, axis=0, return_index=True)[1] # TODO: this is probably taking a while...
    i, ni = i[j], ni[j]
    # finding swaps
    m = np.bitwise_xor.accumulate(smaller[i[:, :-1]] == smaller[i[:, 1:]], axis=1)
    swaps, locs = i[:, 1:][m], np.tile([False, True], (m.shape[0], n))[:, 1:][m]
    src_swaps, dst_swaps = np.unique(swaps[~locs]), np.unique(swaps[locs])
    # apply swaps to src
    k, j = src_swaps, neighbor[src_swaps]
    a, b = r[0, 0, k], r[0, 0, j]
    r[0, 0, k], r[0, 0, j] = b, a
    r[0, 1, a], r[0, 1, b] = r[0, 1, b], r[0, 1, a]
    src_swaps = np.minimum(k, j)
    # apply swaps to dst
    k, j = dst_swaps, neighbor[dst_swaps]
    a, b = r[1, 0, k], r[1, 0, j]
    r[1, 0, k], r[1, 0, j] = b, a
    r[1, 1, a], r[1, 1, b] = r[1, 1, b], r[1, 1, a]
    dst_swaps = np.minimum(k, j)
    return src_swaps, dst_swaps

def _solve(p: Permutation) -> tuple[list, NDArray, list]:
    n = p.size; log2n = (n - 1).bit_length()
    r = np.array([[np.arange(n), np.arange(n)], [p.array_form, (~p).array_form]])
    steps = 2 * log2n - 1; swaps, masks = [None] * steps, [None] * steps
    shifts = 1 << np.arange(log2n); shifts = np.hstack([shifts[:-1], shifts[::-1]]).tolist()
    def mask(x): return reduce(lambda acc, s: acc | (1 << int(s)), x, 0)
    for i, d in enumerate(1 << np.arange(log2n)):
        src, dst = _stage(d, n, r)
        swaps[i], masks[i] = src, mask(src)
        swaps[-(i + 1)], masks[-(i + 1)] = dst, mask(dst)
    return BenesNetwork(masks=masks, shifts=shifts, swaps=swaps)

@lru_cache(maxsize=1000)
def perm2benes(p: Permutation) -> BenesNetwork:
    if p.size == 0: return BenesNetwork(masks=[], shifts=[], swaps=[])
    return _solve(p.resize(2**(p.size - 1).bit_length()))

def unpack2cycle(g: PermutationGroup) -> Callable:
    networks = list(map(perm2benes, g.elements))
    def kernel(bits: int | NDArray):
        if isinstance(bits, numbers.Integral): bits = int(bits); return [b(bits) for b in networks]
        bits = np.asarray(bits); return np.vstack([b(bits) for b in networks])
    return kernel
