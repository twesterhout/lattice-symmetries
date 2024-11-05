from functools import reduce
from lattice_symmetries._kernels import PauliNonbranchingTerm
from numpy.typing import NDArray
from sympy.core.singleton import S
from sympy.physics.quantum import pauli, represent
from sympy.combinatorics import PermutationGroup
import dataclasses
import itertools
import more_itertools
import numbers
import numpy as np
import operator
import sympy
from typing import Mapping
from lattice_symmetries.basis import SpinBasis
from lattice_symmetries._parser import parse_expr

try:
    import igraph
except ImportError:
    igraph = None


def _visit_pauli_matrices(e: sympy.Expr, x, y, z):
    if isinstance(e, pauli.SigmaX):
        return x(e)
    if isinstance(e, pauli.SigmaY):
        return y(e)
    if isinstance(e, pauli.SigmaZ):
        return z(e)
    else:
        msg = f"expected 'e' to be a either SigmaX, SigmaY or SigmaZ, but got {type(e)}"
        raise ValueError(msg)

def simplify_pauli_expression(e: sympy.Expr):
    e = e.replace(
        lambda x: isinstance(x, pauli.SigmaPlus),
        lambda x: (pauli.SigmaX(x.name) + sympy.I * pauli.SigmaY(x.name)) / 2,
    )
    e = e.replace(
        lambda x: isinstance(x, pauli.SigmaMinus),
        lambda x: (pauli.SigmaX(x.name) - sympy.I * pauli.SigmaY(x.name)) / 2,
    )
    e = e.expand()

    # def _simplify_re_im(e):
    #     if isinstance(e, sympy.re) and isinstance(e.args[0], pauli.SigmaOpBase):
    #         return _visit_pauli_matrices(e.args[0], lambda x: x, lambda y: sympy.S.Zero, lambda z: z)
    #     if isinstance(e, sympy.im) and isinstance(e.args[0], pauli.SigmaOpBase):
    #         return _visit_pauli_matrices(e.args[0], lambda x: sympy.S.Zero, lambda y: -sympy.I * y, lambda z: sympy.S.Zero)
    #     return e

    # e = e.replace(lambda x: isinstance(x, (sympy.re, sympy.im)), _simplify_re_im)

    if isinstance(e, sympy.Mul):
        c, nc = e.args_cnc()
        others, paulis = more_itertools.partition(lambda x: isinstance(x, pauli.SigmaOpBase), nc)
        # Sort Pauli matrices
        paulis = sorted(list(paulis), key=lambda x: x.name)
        # Simplify on each site
        paulis = [
            pauli.qsimplify_pauli(sympy.Mul(*g))
            for _, g in itertools.groupby(paulis, lambda x: x.name)
        ]
        paulis = sympy.Mul(*paulis).expand()
        factor = sympy.Mul(*c, *others)
        if isinstance(paulis, sympy.Add):
            return sympy.Add(*[factor * t for t in paulis.args])
        return factor * paulis
    if isinstance(e, sympy.Add):
        return e.func(*map(simplify_pauli_expression, e.args))
    return e


def _pauli_expression_to_nonbranching_terms(e: sympy.Expr) -> list[PauliNonbranchingTerm]:
    f = lambda x: _pauli_expression_to_nonbranching_terms(x)
    if isinstance(e, sympy.Add):
        return list(itertools.chain.from_iterable(map(f, e.args)))
    if isinstance(e, sympy.Mul):
        c, nc = e.args_cnc()
        if len(nc) > 0:
            term = reduce(operator.mul, (more_itertools.one(f(t)) for t in nc))
            return [dataclasses.replace(term, v=sympy.Mul(*c, term.v))]
        else:
            return [PauliNonbranchingTerm(v=sympy.Mul(*c), x=0, s=0)]
    if isinstance(e, pauli.SigmaOpBase):
        if not isinstance(e.name, numbers.Integral):
            raise ValueError(f"expected 'e.name' to be a site index, but got {e.name}")
        i = int(e.name)
        if isinstance(e, pauli.SigmaX):
            return [PauliNonbranchingTerm(v=S.One, x=2**i, s=0)]
        elif isinstance(e, pauli.SigmaY):
            return [PauliNonbranchingTerm(v=-sympy.I, x=2**i, s=2**i)]
        elif isinstance(e, pauli.SigmaZ):
            return [PauliNonbranchingTerm(v=S.One, x=0, s=2**i)]
        else:
            msg = f"expected 'e' to be a either SigmaX, SigmaY or SigmaZ, but got {type(e)}"
            raise ValueError(msg)
    return [PauliNonbranchingTerm(v=e, x=0, s=0)]


def pauli_expression_to_nonbranching_terms(
    e: sympy.Expr, simplify=True
) -> list[PauliNonbranchingTerm]:
    if simplify:
        e = simplify_pauli_expression(e)
    terms = sorted(_pauli_expression_to_nonbranching_terms(e), key=lambda t: (t.x, t.s))
    terms = (
        PauliNonbranchingTerm(v=sympy.Add(*(t.v for t in g)), x=k[0], s=k[1])
        for k, g in itertools.groupby(terms, lambda x: (x.x, x.s))
    )
    return [t for t in terms if t.v != 0]


def _collect_indices(e):
    if isinstance(e, sympy.Add):
        for i in itertools.chain.from_iterable(map(_collect_indices, e.args)):
            yield i
    if isinstance(e, sympy.Mul):
        for i in itertools.chain.from_iterable(map(_collect_indices, e.args)):
            yield i
    if isinstance(e, pauli.SigmaOpBase):
        yield e.name

def _replace_indices(e: sympy.Expr, mapping: Mapping[int, int]):
    if isinstance(e, pauli.SigmaOpBase):
        i = e.name
        if isinstance(i, numbers.Integral):
            i = int(i)
        return e.func(mapping[i])
    if isinstance(e, sympy.Expr):
        if len(e.args) > 0:
            return e.func(*map(lambda x: _replace_indices(x, mapping), e.args))
        else:
            return e
    else:
        raise ValueError(f"expected 'e' to be a SymPy expression, but got {type(e)}")

def pauli_expression_to_matrix(e: sympy.Expr, simplify=True, number_bits=None):
    if simplify:
        e = simplify_pauli_expression(e)
    if number_bits is None:
        number_bits = 1 + max(_collect_indices(e), default=0)

    f = lambda x: pauli_expression_to_matrix(x, simplify=False, number_bits=number_bits)
    if isinstance(e, sympy.Number):
        m = sympy.eye(2**number_bits)
        return e * m
    if isinstance(e, sympy.Add):
        return reduce(operator.add, map(f, e.args))
    if isinstance(e, sympy.Mul):
        c, nc = e.args_cnc()
        terms = [sympy.Matrix([[1, 0], [0, 1]])] * number_bits
        for t in sorted(nc, key=lambda t: t.name):
            terms[number_bits - 1 - t.name] = f(t)
        return sympy.Mul(*c, sympy.kronecker_product(*terms))
    if isinstance(e, pauli.SigmaOpBase):
        return represent(e)
    else:
        raise ValueError(f"expected 'e' to be a either SigmaX, SigmaY or SigmaZ, but got {type(e)}")

class Expr:
    raw: sympy.Expr

    def __init__(self, expression: sympy.Expr | str, sites: list[list[int]] | None = None, particle: str | None = None, simplify: bool = True):
        if isinstance(expression, str):
            expression = parse_expr(expression)
        if simplify:
            expression = simplify_pauli_expression(expression)
        if particle is not None and particle != "spin-1/2":
            raise ValueError(f"expected 'particle' to be 'spin-1/2', but got {particle}")
        self.raw = expression
        if sites is not None:
            self.raw = self.on(sites).raw

    def __str__(self):
        """Get the string representation of the underlying expression."""
        return str(self.raw)

    def __repr__(self):
        return repr(self.raw)

    @property
    def number_spins(self) -> int:
        return 1 + max(_collect_indices(self.raw), default=0)

    @property
    def number_sites(self) -> int:
        return self.number_spins

    def to_dense(self, number_bits: int | None = None, evalf: bool = True):
        if number_bits is None:
            number_bits = self.number_spins
        m = pauli_expression_to_matrix(self.raw, number_bits=number_bits)
        if evalf:
            is_real = all(x.is_real for x in m.flat())
            m = sympy.matrix2numpy(m, dtype=float if is_real else complex)
        return m

    def nonbranching_terms(self):
        return pauli_expression_to_nonbranching_terms(self.raw)

    def permutation_group(self) -> PermutationGroup:
        raise NotImplementedError("😭")

    def abelian_permutation_group(self) -> PermutationGroup:
        raise NotImplementedError("😭")

    def to_json(self) -> str:
        raise NotImplementedError("😭")

    @property
    def is_real(self) -> bool:
        return self == Expr(sympy.re(self.raw))

    @property
    def is_identity(self) -> bool:
        return self == sympy.S.One

    @property
    def is_hermitian(self) -> bool:
        return self == self.adjoint()

    def replace_indices(self, mapping: Mapping[int, int], simplify: bool = True) -> "Expr":
        return Expr(_replace_indices(self.raw, mapping), simplify=simplify)

    def on(self, graph) -> "Expr":
        indices = sorted(frozenset(_collect_indices(self.raw)))
        if igraph is not None and isinstance(graph, igraph.Graph):
            sites = (edge.tuple for edge in graph.es)
        else:
            sites = graph
        args = (_replace_indices(self.raw, dict(zip(indices, edge))) for edge in sites)
        return Expr(sympy.Add(*args))

    def adjoint(self) -> "Expr":
        return Expr(self.raw.adjoint())

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Expr):
            return NotImplemented
        eq = self.raw.equals(other.raw)
        if eq is not None:
            return eq
        eq = simplify_pauli_expression(self.raw - other.raw).equals(sympy.S.Zero)
        if eq is not None:
            return eq
        raise ValueError(f"could not determine equality of {self} and {other}")

    def __add__(self, other: "Expr") -> "Expr":
        if isinstance(other, Expr):
            return Expr(self.raw + other.raw)
        return Expr(self.raw + other)

    def __sub__(self, other: "Expr") -> "Expr":
        if isinstance(other, Expr):
            return Expr(self.raw - other.raw)
        return Expr(self.raw - other)

    def __mul__(self, other: "Expr") -> "Expr":
        if isinstance(other, Expr):
            return Expr(self.raw * other.raw)
        return Expr(self.raw * other)

    def __neg__(self) -> "Expr":
        return Expr(-self.raw)

    def __rmul__(self, other) -> "Expr":
        if isinstance(other, Expr):
            return Expr(other.raw * self.raw)
        return Expr(other * self.raw)

    @property
    def particle_type(self) -> str:
        return "spin-1/2"

    @property
    def conserves_number_particles(self) -> bool:
        raise NotImplementedError("😭")

    @property
    def spin_inversion_invariant(self) -> bool:
        raise NotImplementedError("😭")

    @property
    def conserves_total_spin(self) -> bool:
        raise NotImplementedError("😭")

    def full_basis(self):
        return SpinBasis(number_spins=self.number_spins)

    def hilbert_space_sectors(self):
        raise NotImplementedError("😭")

    def ground_state_sectors(self):
        raise NotImplementedError("😭")