from functools import reduce
from lattice_symmetries._kernels import PauliNonbranchingTerm
from numpy.typing import NDArray
from sympy.core.singleton import S
from sympy.physics.quantum import pauli, represent
import dataclasses
import itertools
import more_itertools
import numbers
import numpy as np
import operator
import sympy


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
