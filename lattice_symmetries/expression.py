from functools import reduce
# from lattice_symmetries._kernels import PauliNonbranchingTerm
from numpy.typing import NDArray
from sympy import Add, Mul, Number, Atom
from sympy.core.singleton import S
from sympy.physics.quantum import pauli, represent
from sympy.combinatorics import Permutation, PermutationGroup
from sympy.physics.quantum.pauli import SigmaOpBase, qsimplify_pauli, SigmaX, SigmaY, SigmaZ, SigmaPlus, SigmaMinus
import dataclasses
import itertools
from itertools import groupby
import more_itertools
import numbers
import numpy as np
import operator
import sympy
from typing import Mapping
# from lattice_symmetries.basis import SpinBasis
from lattice_symmetries._parser import parse_expr

try:
    import igraph
except ImportError:
    igraph = None


def _visit_pauli_matrices(e: sympy.Expr, x, y, z):
    if isinstance(e, SigmaX): return x(e)
    if isinstance(e, SigmaY): return y(e)
    if isinstance(e, SigmaZ): return z(e)
    raise ValueError(f"expected 'e' to be a either SigmaX, SigmaY or SigmaZ, but got {type(e)}")

def simplify_pauli_expression(e: sympy.Expr):
    if e.is_zero: return e
    e = e.replace(lambda x: isinstance(x, SigmaPlus), lambda x: (SigmaX(x.name) + sympy.I * SigmaY(x.name)) / 2)
    e = e.replace(lambda x: isinstance(x, SigmaMinus), lambda x: (SigmaX(x.name) - sympy.I * SigmaY(x.name)) / 2)
    e = e.expand()
    if isinstance(e, Mul):
        c, nc = e.args_cnc()
        others, paulis = more_itertools.partition(lambda x: isinstance(x, SigmaOpBase), nc)
        paulis = sorted(list(paulis), key=lambda x: x.name)
        paulis = Mul(*(qsimplify_pauli(Mul(*g)) for _, g in groupby(paulis, lambda x: x.name))).expand()
        c = Mul(*c, *others)
        return Add(*[c * t for t in paulis.args]) if isinstance(paulis, sympy.Add) else c * paulis
    if isinstance(e, Add): return e.func(*map(simplify_pauli_expression, e.args))
    return e

def _collect_indices(e):
    if isinstance(e, (Add, Mul)): yield from itertools.chain.from_iterable(map(_collect_indices, e.args))
    if isinstance(e, SigmaOpBase): yield e.name

def _replace_indices(e: sympy.Expr, m: Mapping[int, int]):
    if isinstance(e, SigmaOpBase): return e.func(m[int(e.name) if isinstance(e.name, numbers.Integral) else e.name])
    if isinstance(e, sympy.Expr): return e.func(*map(lambda x: _replace_indices(x, m), e.args)) if len(e.args) > 0 else e
    raise ValueError(f"expected 'e' to be a SymPy expression, but got {type(e)}")

def pauli_expression_to_matrix(e: sympy.Expr, simplify=True, number_bits=None):
    if simplify: e = simplify_pauli_expression(e)
    if number_bits is None: number_bits = 1 + max(_collect_indices(e), default=0)
    f = lambda x: pauli_expression_to_matrix(x, simplify=False, number_bits=number_bits)
    if isinstance(e, sympy.Atom): return e * sympy.eye(2**number_bits)
    if isinstance(e, SigmaOpBase): return represent(e)
    if isinstance(e, Add): return reduce(operator.add, map(f, e.args))
    if isinstance(e, Mul):
        c, nc = e.args_cnc()
        terms = [sympy.eye(2)] * number_bits
        for t in sorted(nc, key=lambda t: t.name): terms[number_bits - 1 - t.name] = f(t)
        return Mul(*c, sympy.kronecker_product(*terms))
    raise ValueError(f"Unsupported expression type {type(e)}")

def _pauli2edges(e: sympy.Expr):
    # {a, b, c} means a, b, and c are interconnected
    # [a, b, c] means that edges from a, b, and c should be merged
    def f(e):
        if isinstance(e, Add): return map(f, e.args)
        if isinstance(e, Mul): return reduce(operator.or_, map(f, e.args_cnc()[1]), set())
        if isinstance(e, SigmaOpBase): return {e.name}
        if isinstance(e, sympy.Atom): return set()
        raise ValueError(f"Unsupported expression type {type(e)}")
    def expand(xs):
        if isinstance(xs, set): return itertools.combinations(xs, 2) if len(xs) >= 2 else []
        else: return itertools.chain.from_iterable(map(expand, xs))
    return set(expand(f(e)))

def _d_nd(e):
    """Split into diagonal and non-diagonal parts"""
    def is_off_diag(e):
        if isinstance(e, (SigmaX, SigmaY)): return True
        if isinstance(e, Mul): return any(map(is_off_diag, e.args))
        if isinstance(e, (SigmaZ, Atom)): return False
        assert False, f"shouldn't have been reached: {e}"
    if isinstance(e, Add): d, nd = more_itertools.partition(is_off_diag, e.args); return Add(*d), Add(*nd)
    return (S.Zero, e) if is_off_diag(e) else (e, S.Zero)


class Expr:
    raw: sympy.Expr

    def __init__(self, expression: sympy.Expr | str, sites: list[list[int]] | None = None, particle: str | None = None, simplify: bool = True):
        if isinstance(expression, str): expression = parse_expr(expression)
        if simplify: expression = simplify_pauli_expression(expression)
        if particle is not None and particle != "spin-1/2": raise ValueError(f"expected 'particle' to be 'spin-1/2', but got {particle}")
        self.raw = expression if sites is None else expression.on(sites).raw
    def __str__(self): return str(self.raw)
    def __repr__(self): return repr(self.raw)
    def __add__(self, other: "Expr") -> "Expr": return Expr(self.raw + other.raw) if isinstance(other, Expr) else Expr(self.raw + other)
    def __sub__(self, other: "Expr") -> "Expr": return Expr(self.raw - other.raw) if isinstance(other, Expr) else Expr(self.raw - other)
    def __mul__(self, other: "Expr") -> "Expr":  return Expr(self.raw * other.raw) if isinstance(other, Expr) else Expr(self.raw * other)
    def __rmul__(self, other: "Expr") -> "Expr": return Expr(other.raw * self.raw) if isinstance(other, Expr) else Expr(other * self.raw)
    def __neg__(self) -> "Expr": return Expr(-self.raw)
    def __eq__(self, other) -> bool: return (self - other).is_zero

    @property
    def number_spins(self) -> int: return 1 + max(_collect_indices(self.raw), default=0)
    @property
    def number_sites(self) -> int: return self.number_spins
    @property
    def particle_type(self) -> str: return "spin-1/2"
    @property
    def is_zero(self): return self.raw.is_zero
    @property
    def is_hermitian(self) -> bool: return self == self.adjoint()
    @property
    def is_identity(self) -> bool: return self.is_one

    def adjoint(self) -> "Expr": return Expr(self.raw.adjoint())
    def to_dense(self, number_bits: int | None = None, evalf: bool = True):
        m = pauli_expression_to_matrix(self.raw, number_bits=number_bits or self.number_spins)
        if evalf: m = sympy.matrix2numpy(m, dtype=float if all(x.is_real for x in m.flat()) else complex)
        return m
    def replace_indices(self, mapping: Mapping[int, int], simplify: bool = True) -> "Expr":
        return Expr(_replace_indices(self.raw, mapping), simplify=simplify)
    def on(self, graph) -> "Expr":
        if self.is_zero: return self
        indices = sorted(frozenset(_collect_indices(self.raw)))
        sites = (edge.tuple for edge in graph.es) if igraph is not None and isinstance(graph, igraph.Graph) else graph
        args = (_replace_indices(self.raw, dict(zip(indices, edge))) for edge in sites)
        return Expr(sympy.Add(*args))
    def is_invariant_under(self, p: Permutation): return self == self.replace_indices(dict(enumerate(p.array_form)), simplify=False)
    def permutation_group(self) -> PermutationGroup:
        if igraph is None: raise NotImplementedError("😭")
        g = igraph.Graph(); g.add_vertices(range(self.number_sites)); g.add_edges(_pauli2edges(self.raw))
        pg = PermutationGroup(list(map(Permutation, g.get_isomorphisms_vf2())))
        return PermutationGroup(list(filter(self.is_invariant_under, pg.strong_gens)))

            

        

    # def nonbranching_terms(self): return pauli_expression_to_nonbranching_terms(self.raw)
    # def abelian_permutation_group(self) -> PermutationGroup: raise NotImplementedError("😭")
    # def to_json(self) -> str: raise NotImplementedError("😭")
    # @property
    # def conserves_number_particles(self) -> bool:
    #     raise NotImplementedError("😭")
    # @property
    # def spin_inversion_invariant(self) -> bool:
    #     raise NotImplementedError("😭")
    # @property
    # def conserves_total_spin(self) -> bool:
    #     raise NotImplementedError("😭")
    # def full_basis(self):
    #     return SpinBasis(number_spins=self.number_spins)
    # def hilbert_space_sectors(self):
    #     raise NotImplementedError("😭")
    # def ground_state_sectors(self):
    #     raise NotImplementedError("😭")





def _pre_g(g):
    g = [e.tuple for e in g.es] if igraph is not None and isinstance(g, igraph.Graph) else list(g)
    n = 1 + max(itertools.chain.from_iterable(g), default=0)
    return g, ((i,) for i in range(n))

def ising(g, J=1, h=0):
    g, i = _pre_g(g)
    if (J == 0 and h == 0) or len(g) == 0: return Expr(S.Zero)
    return J * Expr("σᶻ₀ σᶻ₁").on(g) - h * Expr("σˣ₀").on(i)

def heisenberg(g, J=1, h=0):
    g, i = _pre_g(g)
    if (J == 0 and h == 0) or len(g) == 0: return Expr(S.Zero)
    return Expr("2 (σ⁺₀ σ⁻₁ + σ⁺₁ σ⁻₀) + σᶻ₀ σᶻ₁").on(g) - h * Expr("σᶻ₀").on(i)


@dataclasses.dataclass(frozen=True)
class PauliNonbranchingTerm:
    v: sympy.Expr
    x: int
    s: int

    @property
    def is_diagonal(self) -> bool: return self.x == 0
    def act_on_ket(self, ket: int) -> tuple[sympy.Expr, int]:
        sign = 1 - 2 * ((ket & self.s).bit_count() % 2)
        return sign * self.v, ket ^ self.x
    def act_on_bra(self, bra: int) -> tuple[sympy.Expr, int]:
        coeff, beta = self.act_on_ket(bra)
        return coeff.conjugate(), beta
    def __mul__(self, other):
        if not isinstance(other, PauliNonbranchingTerm):
            return NotImplemented
        return PauliNonbranchingTerm(v=self.v * other.v, x=self.x ^ other.x, s=self.s ^ other.s)

def _pauli2nbts(e: sympy.Expr) -> list:
    if isinstance(e, Add): return list(itertools.chain.from_iterable(map(_pauli2nbts, e.args)))
    if isinstance(e, Mul):
        c, nc = e.args_cnc()
        if len(nc) > 0:
            nbt = reduce(operator.mul, (more_itertools.one(_pauli2nbts(t)) for t in nc))
            return [dataclasses.replace(nbt, v=sympy.Mul(*c, nbt.v))]
        return [PauliNonbranchingTerm(v=sympy.Mul(*c), x=0, s=0)]
    if isinstance(e, SigmaOpBase):
        i = int(e.name)
        return [_visit_pauli_matrices(e,
            lambda t: PauliNonbranchingTerm(v=S.One, x=2**i, s=0),
            lambda t: PauliNonbranchingTerm(v=-sympy.I, x=2**i, s=2**i),
            lambda t: PauliNonbranchingTerm(v=S.One, x=0, s=2**i)
        )]
    return [PauliNonbranchingTerm(v=e, x=0, s=0)]

def pauli2nbts(e: sympy.Expr, simplify=True) -> list:
    """Convert an expression with Pauli matrices into a list of nonbranching terms."""
    if simplify: e = simplify_pauli_expression(e)
    key = lambda t: (t.x, t.s)
    term = lambda arg: PauliNonbranchingTerm(v=Add(*(t.v for t in arg[1])), x=arg[0][0], s=arg[0][1])
    return list(filter(lambda t: not t.v.is_zero, map(term, groupby(sorted(_pauli2nbts(e), key=key), key))))

class PauliLoweredTerms:
    has_diag: bool
    has_off_diag: bool
    is_real: bool
    is_imag: bool
    s_diag: NDArray[np.uint64]
    v_re_diag: NDArray[np.float64]
    v_im_diag: NDArray[np.float64]
    s_2d: NDArray[np.uint64]
    v_re_2d: NDArray[np.float64]
    v_im_2d: NDArray[np.float64]
    mask: NDArray[np.uint64]

    def __init__(self, terms: list[PauliNonbranchingTerm]):
        v = np.asarray([complex(t.v) for t in terms], dtype=np.complex128)
        x = np.asarray([t.x for t in terms], dtype=np.uint64)
        s = np.asarray([t.s for t in terms], dtype=np.uint64)
        n_diag = np.sum(x == 0)
        self.has_diag = n_diag != 0
        self.has_off_diag = n_diag < len(x)
        self.is_real = np.all(v.imag == 0)
        self.is_imag = np.all(v.real == 0) and not self.is_real
        self.s_diag = s[:n_diag] if self.has_diag else None
        self.v_re_diag = np.ascontiguousarray(v[:n_diag].real) if not self.is_imag else None
        self.v_im_diag = np.ascontiguousarray(v[:n_diag].imag) if not self.is_real else None
        if self.has_off_diag:
            unique_xs, counts = np.unique(x[n_diag:], return_counts=True)
            n_terms, n_reduced = counts.size, np.max(counts)
            self.s_2d = np.zeros((n_terms, n_reduced), dtype=np.uint64)
            self.v_re_2d = np.zeros((n_terms, n_reduced), dtype=np.float64) if not self.is_imag else None
            self.v_im_2d = np.zeros((n_terms, n_reduced), dtype=np.float64) if not self.is_real else None
            offsets = n_diag + np.pad(np.cumsum(counts), ((1, 0),))
            for i in range(n_terms):
                self.s_2d[i, :counts[i]] = s[offsets[i] : offsets[i] + counts[i]]
                if self.v_re_2d is not None: self.v_re_2d[i, :counts[i]] = v[offsets[i] : offsets[i] + counts[i]].real
                if self.v_im_2d is not None: self.v_im_2d[i, :counts[i]] = v[offsets[i] : offsets[i] + counts[i]].imag
            self.mask = unique_xs
        else:
            self.s_2d = None
            self.v_re_2d = None
            self.v_im_2d = None
            self.mask = None

    def diag(self, alpha: NDArray[np.uint64]) -> NDArray:
        if not self.has_diag: raise ValueError("don't call diag() when has_diag is False")
        assert alpha.dtype == np.dtype("uint64")
        sign = np.bitwise_count(alpha.reshape(1, -1) & self.s_diag.reshape(-1, 1)).astype(np.uint64) << 63
        re = (self.v_re_diag.reshape(-1, 1).view(np.uint64) ^ sign).view(np.float64).sum(axis=0) if not self.is_imag else None
        im = (self.v_im_diag.reshape(-1, 1).view(np.uint64) ^ sign).view(np.float64).sum(axis=0) if not self.is_real else None
        if self.is_real: return re
        if self.is_imag: return 1j * im
        return re + 1j * im

    def off_diag(self, alpha: NDArray[np.uint64]) -> tuple[NDArray, NDArray]:
        if not self.has_off_diag: raise ValueError("don't call off_diag() when has_off_diag is False")
        assert alpha.dtype == np.dtype("uint64")
        sign = np.bitwise_count(alpha[None, None, :] & self.s_2d[:, :, None]).astype(np.uint64) << 63
        re = (self.v_re_2d[:, :, None].view(np.uint64) ^ sign).view(np.float64).sum(axis=1) if not self.is_imag else None
        im = (self.v_im_2d[:, :, None].view(np.uint64) ^ sign).view(np.float64).sum(axis=1) if not self.is_real else None
        beta = alpha[None, :] ^ self.mask[:, None]
        if self.is_real: return beta, re
        if self.is_imag: return beta, 1j * im
        return beta, re + 1j * im
