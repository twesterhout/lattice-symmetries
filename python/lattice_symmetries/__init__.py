# Copyright (c) 2022-2024, Tom Westerhout
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import json
import typing
import weakref
from functools import cached_property
from typing import Any, Callable, Dict, Optional, Tuple

import numpy as np
from lattice_symmetries._ls_hs import ffi, lib
from loguru import logger
from numpy.typing import NDArray
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import LinearOperator
from sympy.combinatorics import Permutation, PermutationGroup
from sympy.core import Rational

try:
    import igraph as ig
except:
    ig = None

__version__ = "2.2.0"


class _RuntimeInitializer:
    def __init__(self):
        logger.trace("Initializing Haskell runtime...")
        lib.ls_hs_init()
        logger.trace("Initializing Chapel runtime...")
        lib.ls_chpl_init()

    def __del__(self):
        # NOTE: The order of these should actually be reversed, but ls_chpl_finalize calls exit(0) :/
        logger.trace("Deinitializing Haskell runtime...")
        lib.ls_hs_exit()
        logger.trace("Deinitializing Chapel runtime...")
        # sys.stdout.flush()
        # sys.stderr.flush()
        lib.ls_chpl_finalize()


# NOTE: even though _runtime_init is never accessed, we need the variable to
# force initialization of Chapel and Haskell runtimes
_runtime_init = _RuntimeInitializer()


class LatticeSymmetriesEncoder(json.JSONEncoder):
    """Extension of the json.JSONEncoder class to support encoding types required by the Haskell code."""

    def default(self, o):
        if isinstance(o, Rational):
            name = "Rational"
        else:
            name = type(o).__name__

        try:
            encoder = getattr(self, f"encode_{name}")
        except AttributeError:
            super().default(o)
        else:
            encoded = encoder(o)
            if isinstance(encoded, dict):
                encoded["__type__"] = name
            return encoded

    def encode_Permutation(self, obj: Permutation):
        return obj.list()

    def encode_Rational(self, obj: Rational):
        return f"{obj.numerator}/{obj.denominator}"

    def encode_complex(self, obj: complex):
        return {"real": obj.real, "imag": obj.imag}


class LatticeSymmetriesDecoder(object):
    def object_hook(self, obj):
        if isinstance(obj, dict):
            if "Left" in obj or "Right" in obj:
                return self.decode_Result(obj)
            if "__type__" in obj:
                name = obj["__type__"]
                decoder = getattr(self, f"decode_{name}", lambda x: x)
                return decoder(obj)
        return obj

    def decode_Result(self, obj):
        if "Left" in obj:
            raise ValueError(obj["Left"])
        elif "Right" in obj:
            return self.object_hook(obj["Right"])
        else:
            return obj

    # def decode_Symmetry(self, obj):
    #     return Symmetry(permutation=obj["permutation"], sector=obj["sector"])

    def decode_Complex(self, obj):
        return complex(obj["real"], obj["imag"])


def _to_json(obj) -> bytes:
    return json.dumps(obj, cls=LatticeSymmetriesEncoder).encode("utf-8")


def _from_haskell_string(c_str) -> str:
    s = ffi.string(c_str).decode("utf-8")
    lib.ls_hs_destroy_string(c_str)
    return s


def _from_json(obj) -> Any:
    if isinstance(obj, ffi.CData):
        obj = _from_haskell_string(obj)
    return json.loads(obj, object_hook=LatticeSymmetriesDecoder().object_hook)


def _assert_subtype(variable, required_type):
    if not isinstance(variable, required_type):
        raise TypeError("expected a '{}', but got '{}'".format(required_type, type(variable)))


class ExternalArrayWrapper:
    payload: ffi.CData
    keep_alive: Optional[Any]
    finalizer: Optional[weakref.finalize]

    def __init__(
        self,
        arr: ffi.CData,
        typestr: str,
        keep_alive: Optional[Any] = None,
        finalizer=lib.ls_hs_internal_destroy_external_array,
    ):
        self.payload = arr
        self.keep_alive = keep_alive

        if finalizer is not None:
            self.finalizer = weakref.finalize(self, finalizer, self.payload)
        else:
            self.finalizer = None

        n = arr.num_elts
        p = int(ffi.cast("uintptr_t", arr.elts))
        self.__array_interface__ = {
            "version": 3,
            "typestr": typestr,
            # TODO: NumPy is buggy and doesn't work if p == 0
            "data": (p if p != 0 else -1, True),
            "shape": (n,),
        }


class HsWrapper(object):
    _payload: ffi.CData
    _finalizer: Optional[weakref.finalize]

    def __init__(self, payload: ffi.CData, finalizer: Optional[Callable[[ffi.CData], None]]):
        assert isinstance(payload, ffi.CData)
        assert payload != 0
        self._payload = payload
        self._finalizer = (
            weakref.finalize(self, finalizer, self._payload) if finalizer is not None else None
        )


class Basis(HsWrapper):
    def __init__(
        self,
        json_string: str = "",
        payload: ffi.CData | int = 0,
        finalizer: Callable[[ffi.CData], None] = lib.ls_hs_destroy_basis,
    ):
        if len(json_string) > 0:
            payload = _from_json(lib.ls_hs_basis_from_json(json_string.encode("utf-8")))
        elif payload == 0:
            raise ValueError(
                "incompatible arguments: either 'payload' and 'finalizer' "
                "or 'json_string' may be specified"
            )

        if isinstance(payload, int):
            payload = ffi.cast("ls_hs_basis *", payload)

        super().__init__(payload=payload, finalizer=finalizer)
        # NOTE: Important!! Pre-initialize the basis info
        lib.ls_hs_init_basis_info(self._payload)

    def _get_info(self):
        return self._payload.info

    @property
    def hamming_weight(self) -> Optional[int]:
        h = self._get_info().hamming_weight
        return int(h) if h != -1 else None

    @property
    def spin_inversion(self) -> Optional[int]:
        i = self._get_info().spin_inversion
        return int(i) if i != 0 else None

    @property
    def min_state_estimate(self) -> int:
        return int(self._get_info().min_state_estimate)

    @property
    def max_state_estimate(self) -> int:
        return int(self._get_info().max_state_estimate)

    @property
    def number_sites(self) -> int:
        return int(self._get_info().number_sites)

    @property
    def number_bits(self) -> int:
        return int(self._get_info().number_bits)

    @property
    def number_words(self) -> int:
        return int(self._get_info().number_words)

    @property
    def has_permutation_symmetries(self) -> bool:
        return bool(self._get_info().has_permutation_symmetries)

    @property
    def requires_projection(self) -> bool:
        return bool(self._get_info().requires_projection)

    @property
    def is_state_index_identity(self) -> bool:
        return bool(self._get_info().is_state_index_identity)

    @property
    def is_real(self) -> bool:
        return bool(self._get_info().is_real)

    @property
    def is_built(self) -> bool:
        return bool(self._payload.is_built)

    @cached_property
    def symmetries(self) -> list[tuple[Permutation, Rational]]:
        elements = _from_json(lib.ls_hs_basis_permutation_group(self._payload))
        return [(Permutation(p), Rational(n, m)) for p, n, m in elements]

    @cached_property
    def permutation_group(self) -> PermutationGroup:
        return PermutationGroup([p for p, _ in self.symmetries])

    def check_is_built(self):
        if not self.is_built:
            raise ValueError(
                "basis states have not been built yet; "
                "if you wish to do so, use the basis.build() function"
            )

    def build(self) -> None:
        """Generate a list of representatives.

        These can later be accessed using the `number_states` and `states` attributes.
        """
        if not self.is_built:
            lib.ls_chpl_local_enumerate_states(
                self._payload, self.min_state_estimate, self.max_state_estimate
            )
        assert self.is_built

    @property
    def number_states(self) -> int:
        self.check_is_built()
        return self._payload.local_representatives.num_elts

    @property
    def states(self) -> NDArray[np.uint64]:
        self.check_is_built()
        arr = ExternalArrayWrapper(
            arr=self._payload.local_representatives,
            typestr="u8",
            keep_alive=self,
            finalizer=None,
        )
        return np.array(arr, copy=False)

    @property
    def norms(self) -> NDArray[np.uint16]:
        self.check_is_built()
        arr = ExternalArrayWrapper(
            arr=self._payload.local_norms,
            typestr="u2",
            keep_alive=self,
            finalizer=None,
        )
        return np.array(arr, copy=False)

    def is_representative(self, x: int | NDArray[np.uint64]) -> int | NDArray[np.uint16]:
        is_scalar = False
        x = np.asarray(x, dtype=np.uint64, order="C")
        if x.ndim == 0:
            is_scalar = True
            x = np.expand_dims(x, axis=0)
        if x.ndim != 1:
            raise ValueError(f"'x' has invalid shape: {x.shape}; expected a one-dimensional array")

        count = x.size
        norms = np.empty(x.size, dtype=np.uint16)
        lib.ls_chpl_is_representative(
            self._payload,
            count,
            ffi.from_buffer("uint64_t const*", x, require_writable=False),
            ffi.from_buffer("uint16_t *", norms, require_writable=True),
        )

        if is_scalar:
            return int(norms[0])
        else:
            return norms

    def index(self, x: int | NDArray[np.uint64]) -> int | NDArray[np.int64]:
        """Return the index of a basis state or a batch of basis states."""
        if self.number_bits > 64:
            raise ValueError(
                "it is impractical to compute indices of states with more than 64 bits"
            )

        is_scalar = False
        x = np.asarray(x, dtype=np.uint64, order="C")
        if x.ndim == 0:
            is_scalar = True
            x = np.expand_dims(x, axis=0)
        if x.ndim != 1:
            raise ValueError(f"'x' has invalid shape: {x.shape}; expected a one-dimensional array")

        if self.is_state_index_identity:
            indices = x.astype(np.int64)
        else:
            if self.has_permutation_symmetries:
                self.check_is_built()

            count = x.shape[0]
            indices = np.zeros(count, dtype=np.int64)
            x_ptr = ffi.from_buffer("uint64_t const*", x, require_writable=False)
            indices_ptr = ffi.from_buffer("int64_t *", indices, require_writable=True)
            lib.ls_chpl_basis_state_to_index(self._payload, count, x_ptr, indices_ptr)

        if is_scalar:
            return int(indices[0])
        else:
            return indices


#
#     def build(self) -> None:
#         """Generate a list of representatives.
#
#         These can later be accessed using the `number_states` and `states` attributes.
#         """
#         if not self.is_built:
#             lib.ls_hs_basis_build(self._payload)
#         assert self.is_built
#
#     @property
#     def number_states(self) -> int:
#         self.check_is_built()
#         return int(self._payload.representatives.num_elts)
#
#     @property
#     def states(self) -> NDArray[np.uint64]:
#         n = self.number_states
#         p = int(ffi.cast("uintptr_t", self._payload.representatives.elts))
#         # int(ffi.cast("uintptr_t", lib.ls_hs_basis_states(self._payload)))
#
#         class StatesWrapper(object):
#             def __init__(self, wrapped):
#                 self.wrapped = wrapped
#                 self.__array_interface__ = {
#                     "version": 3,
#                     "typestr": "u8",
#                     "data": (p, True),
#                     "shape": (n,),
#                 }
#
#         return np.array(StatesWrapper(self), copy=False)
#


class SpinBasis(Basis):
    def __init__(
        self,
        number_spins: int,
        hamming_weight: int | None = None,
        spin_inversion: int | None = None,
        symmetries: list[Tuple[Permutation, Rational]] = [],
    ):
        """Create a Hilbert space basis for `number_spins` spin-1/2 particles."""
        super().__init__(
            json_string=json.dumps(
                {
                    "particle": "spin-1/2",
                    "number_sites": number_spins,
                    "hamming_weight": hamming_weight,
                    "spin_inversion": spin_inversion,
                    "symmetries": symmetries,
                },
                cls=LatticeSymmetriesEncoder,
            )
        )


class Expr(HsWrapper):
    def __init__(
        self,
        expression: str = "",
        sites: Any = None,
        particle: Optional[str] = None,
        payload: ffi.CData | int = 0,
        finalizer: Callable[[ffi.CData], None] = lib.ls_hs_destroy_expr,
    ):
        if len(expression) > 0:
            json_object: Dict[str, Any] = {"expression": expression}
            if sites is not None:
                if ig is not None and isinstance(sites, ig.Graph):
                    json_object["sites"] = [list(edge.tuple) for edge in sites.es]
                else:
                    json_object["sites"] = sites
            if particle is not None:
                json_object["particle"] = particle
            payload = _from_json(lib.ls_hs_expr_from_json(json.dumps(json_object).encode("utf-8")))
        else:
            if not (len(expression) == 0 and sites is None and particle is None):
                raise ValueError(
                    "incompatible arguments: either 'payload' and 'finalizer' "
                    "or 'expression', 'sites', and 'particle' may be specified"
                )
        if isinstance(payload, int):
            payload = ffi.cast("ls_hs_expr *", payload)
        super().__init__(payload=payload, finalizer=finalizer)

    def permutation_group(self) -> PermutationGroup:
        obj = _from_json(lib.ls_hs_expr_permutation_group(self._payload))
        return PermutationGroup([Permutation(p) for p in obj["permutations"]])

    def abelian_permutation_group(self) -> PermutationGroup:
        obj = _from_json(lib.ls_hs_expr_abelian_permutation_group(self._payload))
        g = PermutationGroup([Permutation(p) for p in obj["permutations"]])
        assert g.is_abelian
        return g

    def to_json(self) -> str:
        c_str = lib.ls_hs_expr_to_json(self._payload)
        s = ffi.string(c_str).decode("utf-8")
        lib.ls_hs_destroy_string(c_str)
        return s

    def __str__(self) -> str:
        """Get the string representation of the underlying expression."""
        return _from_haskell_string(lib.ls_hs_expr_to_string(self._payload))

    @property
    def is_real(self) -> bool:
        return lib.ls_hs_expr_is_real(self._payload) != 0

    @property
    def is_identity(self) -> bool:
        return lib.ls_hs_expr_is_identity(self._payload) != 0

    @property
    def is_hermitian(self) -> bool:
        return lib.ls_hs_expr_is_hermitian(self._payload) != 0

    def replace_indices(self, mapping: Dict[int, int]) -> "Expr":
        return Expr(
            payload=_from_json(
                lib.ls_hs_replace_indices(
                    self._payload,
                    _to_json([(int(k), int(v)) for k, v in mapping.items()]),
                )
            )
        )

    def on(self, graph) -> "Expr":
        if ig is not None and isinstance(graph, ig.Graph):
            return Expr(expression=str(self), sites=graph)
        return Expr(expression=str(self), sites=[list(edge) for edge in graph])

    def adjoint(self) -> "Expr":
        return Expr(payload=lib.ls_hs_expr_adjoint(self._payload))

    def scale(self, coeff: complex) -> "Expr":
        coeff = complex(coeff)
        return Expr(payload=lib.ls_hs_expr_scale(coeff.real, coeff.imag, self._payload))

    def __eq__(self, other: object) -> bool:
        _assert_subtype(other, Expr)
        other = typing.cast(Expr, other)
        return lib.ls_hs_expr_equal(self._payload, other._payload) != 0

    def __add__(self, other: "Expr") -> "Expr":
        _assert_subtype(other, Expr)
        return Expr(payload=_from_json(lib.ls_hs_expr_plus(self._payload, other._payload)))

    def __sub__(self, other: "Expr") -> "Expr":
        _assert_subtype(other, Expr)
        return Expr(payload=_from_json(lib.ls_hs_expr_minus(self._payload, other._payload)))

    def __mul__(self, other: "Expr") -> "Expr":
        _assert_subtype(other, Expr)
        return Expr(payload=_from_json(lib.ls_hs_expr_times(self._payload, other._payload)))

    def __neg__(self) -> "Expr":
        return Expr(payload=lib.ls_hs_expr_negate(self._payload))

    def __rmul__(self, other: complex) -> "Expr":
        if np.isscalar(other):
            return self.scale(typing.cast(complex, other))
        else:
            return NotImplemented

    @property
    def particle_type(self) -> str:
        return _from_haskell_string(lib.ls_hs_expr_particle_type(self._payload))

    @cached_property
    def conserves_number_particles(self) -> bool:
        return bool(lib.ls_hs_expr_conserves_number_particles(self._payload))

    @cached_property
    def spin_inversion_invariant(self) -> bool:
        return bool(lib.ls_hs_expr_spin_inversion_invariant(self._payload))

    @cached_property
    def conserves_total_spin(self) -> bool:
        if self.particle_type == "spin-1/2":
            s = Expr(
                "σˣ₀ σˣ₁ + σʸ₀ σʸ₁ + σᶻ₀ σᶻ₁",
                sites=[(i, j) for i in range(self.number_sites) for j in range(self.number_sites)],
            )
            return self * s == s * self
        else:
            return NotImplemented

    @property
    def number_sites(self) -> int:
        return int(lib.ls_hs_expr_number_sites(self._payload))

    def full_basis(self) -> Basis:
        return Basis(
            json_string=json.dumps(
                dict(particle=self.particle_type, number_sites=self.number_sites)
            )
        )

    def hilbert_space_sectors(self):
        r = _from_json(lib.ls_hs_expr_hilbert_space_sectors(self._payload))
        return [Basis(payload=p) for p in r]

    def ground_state_sectors(self):
        r = _from_json(lib.ls_hs_expr_ground_state_sectors(self._payload))
        return [Basis(payload=p) for p in r]

    def to_dense(self) -> NDArray:
        m = Operator(self)
        m.basis.build()
        return m.to_dense()


class Operator(HsWrapper, LinearOperator):
    _basis: Basis
    _expression: Expr
    _off_diag_csr: csr_matrix | None
    _diag_arr: NDArray | None

    def __init__(
        self,
        expression: Optional[Expr] = None,
        basis: Optional[Basis] = None,
        payload: ffi.CData | int = 0,
        finalizer: Callable[[ffi.CData], None] = lib.ls_hs_destroy_operator,
    ):
        if expression is not None:
            _assert_subtype(expression, Expr)
            if basis is None:
                basis = expression.full_basis()
            _assert_subtype(basis, Basis)

            payload = _from_json(lib.ls_hs_create_operator(basis._payload, expression._payload))
        else:
            if basis is not None or payload == 0:
                raise ValueError(
                    "incompatible arguments: either 'payload' and 'finalizer' "
                    "or 'expression' and 'basis' may be specified"
                )

        if isinstance(payload, int):
            payload = ffi.cast("ls_hs_operator *", payload)

        # Manually increase refcounts, because self.basis and self.expression
        # Python objects will automatically with decrease them on destruction
        lib.ls_hs_internal_object_inc_ref_count(ffi.addressof(payload.basis, "base"))  # type: ignore
        lib.ls_hs_internal_object_inc_ref_count(ffi.addressof(payload.expr, "base"))  # type: ignore
        self._basis = Basis(payload=payload.basis)  # type: ignore
        self._expression = Expr(payload=payload.expr)  # type: ignore
        self._off_diag_csr = None
        self._diag_arr = None

        super().__init__(payload=payload, finalizer=finalizer)

    @staticmethod
    def from_json(json_string: str) -> "Operator":
        return Operator(
            payload=_from_json(lib.ls_hs_operator_from_json(json_string.encode("utf-8")))
        )

    @property
    def shape(self) -> Tuple[int, int]:
        n = self.basis.number_states
        return (n, n)

    @property
    def dtype(self):
        if self.basis.is_real and self.expression.is_real:
            return np.dtype("float64")
        else:
            return np.dtype("complex128")

    @property
    def basis(self):
        return self._basis

    @property
    def expression(self):
        return self._expression

    def apply_to_state_vector(self, vector: NDArray, out: Optional[NDArray] = None) -> NDArray:
        self.basis.check_is_built()

        if vector.ndim > 2 or vector.shape[0] != self.basis.number_states:
            raise ValueError(
                "'vector' has invalid shape: {}; expected {}"
                "".format(vector.shape, (self.basis.number_states,))
            )

        operator_is_real = self.basis.is_real and self.expression.is_real
        vector_is_real = np.isrealobj(vector)

        if vector_is_real and operator_is_real:
            vector = np.asfortranarray(vector, dtype=np.float64)
            c_type_str = "double[]"
            fn = lib.ls_chpl_matrix_vector_product_f64
        else:
            vector = np.asfortranarray(vector, dtype=np.complex128)
            c_type_str = "ls_hs_scalar[]"
            fn = lib.ls_chpl_matrix_vector_product_c128

        if out is None:
            out = np.zeros_like(vector)
        else:
            if out.shape != vector.shape:
                raise ValueError(
                    "'out' has invalid shape: {}; expected {}" "".format(out.shape, vector.shape)
                )
            if not out.flags.f_contiguous:
                raise ValueError("'out' must be Fortran contiguous")

        if self._diag_arr is not None and self._off_diag_csr is not None:
            _csr_matvec(self._off_diag_csr, vector, out)
            out += self._diag_arr * vector
        else:
            x_ptr = ffi.from_buffer(c_type_str, vector, require_writable=False)
            y_ptr = ffi.from_buffer(c_type_str, out, require_writable=True)
            num_vectors = out.shape[1] if out.ndim > 1 else 1
            fn(self._payload, num_vectors, x_ptr, y_ptr)
        return out

    def _matvec(self, x):
        return self.apply_to_state_vector(x)

    #
    #     # def adjoint(self):
    #     #     return Operator(self.basis, lib.ls_hs_operator_hermitian_conjugate(self._payload))
    #
    #     # @property
    #     # def is_identity(self) -> bool:
    #     #     return lib.ls_hs_operator_is_identity(self._payload)
    #
    #     # @property
    #     # def is_hermitian(self) -> bool:
    #     #     return lib.ls_hs_operator_is_hermitian(self._payload)
    #     @property
    #     def is_real(self) -> bool:
    #         return self.basis.is_real and self.expression.is_real
    #
    def __add__(self, other):
        """Add two operators."""
        if isinstance(other, Operator):
            return Operator(self.expression + other.expression, self.basis)
        else:
            return NotImplemented

    def __sub__(self, other):
        """Subtract `other` Operator from `self`."""
        if isinstance(other, Operator):
            return Operator(self.expression - other.expression, self.basis)
        else:
            return NotImplemented

    def scale(self, coeff: complex) -> "Operator":
        return Operator(self.expression.scale(coeff), self.basis)

    def __mul__(self, other):
        """Multiply `self` by `other` Operator."""
        if isinstance(other, Operator):
            return Operator(self.expression * other.expression, self.basis)
        else:
            return NotImplemented

    def __rmul__(self, other: complex):
        """Scale `self` by a scalar `other`."""
        if np.isscalar(other):
            return self.scale(typing.cast(complex, other))
        else:
            return NotImplemented

    def __matmul__(self, other: NDArray[Any]) -> NDArray[Any]:
        return self.apply_to_state_vector(other)

    def __repr__(self):
        return "<Operator defined on {}>".format(self.basis.__class__.__name__)

    def diag_to_array(self) -> NDArray:
        if self._diag_arr is not None:
            return self._diag_arr

        self.basis.check_is_built()
        diag = ffi.new("chpl_external_array *")
        lib.ls_chpl_extract_diag_c128(self._payload, diag)
        arr = np.array(ExternalArrayWrapper(arr=diag, typestr="c16"), copy=False)
        # Our Chapel code builds everything using complex128 even when the
        # operator is real, so we manually cast to float64 when it's safe
        if self.dtype == np.dtype("float64"):
            arr = np.ascontiguousarray(arr.real)
        assert arr.dtype == self.dtype
        return arr

    def off_diag_to_triple(
        self,
        states: NDArray[np.uint64] | None = None,
        norms: NDArray[np.uint16] | None = None,
        convert_to_index: bool = True,
    ) -> tuple[NDArray, NDArray[np.uint64] | NDArray[np.int64], NDArray[np.int64]]:
        if states is None and norms is None:
            states = self.basis.states
            norms = self.basis.norms
        elif states is not None:
            if states.ndim != 1:
                raise ValueError(f"'states' has wrong shape: {states.shape}")
            states = np.require(states, dtype=np.uint64, requirements=["C"])
        else:
            raise ValueError("'states' must not be None when 'norms' is not None")

        if norms is None:
            norms = self.basis.is_representative(states)
        else:
            if norms.shape != states.shape:
                raise ValueError(f"'norms' has wrong shape: {norms.shape}; expected {states.shape}")
            norms = np.require(norms, dtype=np.uint16, requirements=["C"])

        c_basis_states = ffi.new("chpl_external_array *")
        c_basis_states.num_elts = states.size
        c_basis_states.elts = ffi.from_buffer("uint64_t[]", states)
        c_norms = ffi.new("chpl_external_array *")
        c_norms.num_elts = norms.size
        c_norms.elts = ffi.from_buffer("uint16_t[]", norms)

        c_coeffs = ffi.new("chpl_external_array *")
        c_offsets = ffi.new("chpl_external_array *")
        c_indices = ffi.new("chpl_external_array *")
        lib.ls_chpl_off_diag_to_csr_c128(
            self._payload,
            c_basis_states,
            c_norms,
            c_coeffs,
            c_offsets,
            c_indices,
            convert_to_index,
        )

        coeffs_arr = np.array(ExternalArrayWrapper(c_coeffs, "c16"), copy=False)
        # Our Chapel code builds everything using complex128 even when the
        # operator is real, so we manually cast to float64 when it's safe
        if self.dtype == np.dtype("float64"):
            coeffs_arr = np.ascontiguousarray(coeffs_arr.real)
        assert coeffs_arr.dtype == self.dtype

        offsets_arr = np.array(ExternalArrayWrapper(c_offsets, "i8"), copy=False)
        if convert_to_index:
            indices_arr = np.array(ExternalArrayWrapper(c_indices, "i8"), copy=False)
        else:
            indices_arr = np.array(ExternalArrayWrapper(c_indices, "u8"), copy=False)

        return (coeffs_arr, indices_arr, offsets_arr)

    def off_diag_to_csr(self) -> csr_matrix:
        if self._off_diag_csr is not None:
            return self._off_diag_csr
        n = self.basis.number_states
        return csr_matrix(self.off_diag_to_triple(), shape=(n, n), dtype=self.dtype)

    def to_partial_csr(
        self,
        states: NDArray[np.uint64],
        norms: NDArray[np.uint16] | None = None,
        add_diagonal=True,
    ) -> csr_matrix:
        """
        Constructs a sparse matrix that is a slice of the Hamiltonian matrix.
        Included rows are determined by ``states``, included columns are ``states``
        and their first neighbors.

        Parameters
        ----------
        states : NDArray[np.uint64]
            The states whose neighbors to include in the matrix.

        norms : NDArray[np.uint16]
            The corresponding norms, as returned by `self.basis.is_representative(states)`
            If None, calculated automatically

        add_diagonal: bool (default: True)
            Whether to add a diagonal

        Returns
        -------
        matrix : csr_matrix
            The sparse matrix.

        all_states : npt.NDArray[np.uint64]
            The sorted array of states and their first neighbors.

        For each ``i`` and ``j``, the following holds
        (provided ``add_diagonal=True``):

        ``matrix[i, j] = <states[i] | H | all_states[j]>``

        or, equivalently,

        ``(H |psi>)[states] = matrix @ (psi[all_states])``
        """
        if states.ndim != 1:
            raise ValueError(f"'states' has wrong shape: {states.shape}")
        states = np.require(states, dtype=np.uint64, requirements=["C"])

        if norms is None:
            norms = self.basis.is_representative(states)
        else:
            if norms.shape != states.shape:
                raise ValueError(f"'norms' has wrong shape: {norms.shape}; expected {states.shape}")
            norms = np.require(norms, dtype=np.uint16, requirements=["C"])

        coeffs, other_states, offsets = self.off_diag_to_triple(
            states=states, norms=norms, convert_to_index=False
        )
        all_states = np.unique(np.concatenate([other_states, states]))
        other_indices = np.searchsorted(all_states, other_states)

        matrix_without_diag = csr_matrix(
            (coeffs, other_indices, offsets),
            shape=(states.size, all_states.size),
        )
        if add_diagonal:
            # TODO: fix me
            diag_coeffs = self.diag_to_array()[self.basis.index(states)]
            diag_indices = np.searchsorted(all_states, states)
            diag_matrix = csr_matrix(
                (diag_coeffs, diag_indices, np.arange(states.size + 1)),
                shape=matrix_without_diag.shape,
            )
            matrix = matrix_without_diag + diag_matrix

        matrix.sum_duplicates()
        return matrix, all_states

    def build_matrix(self):
        if self._off_diag_csr is None:
            self._off_diag_csr = self.off_diag_to_csr()
        if self._diag_arr is None:
            self._diag_arr = self.diag_to_array()

    def to_dense(self) -> NDArray[Any]:
        m = self.off_diag_to_csr().todense()
        m += np.diag(self.diag_to_array())
        return m

    def estimate_max_number_off_diag(self) -> int:
        return lib.ls_hs_operator_max_number_off_diag(self._payload)


def _axpy(alpha: complex, x: NDArray, y: NDArray):
    if x.dtype != y.dtype or x.dtype != np.dtype("complex128"):
        raise ValueError("_axpy currently only supports complex128")
    if x.shape != y.shape:
        raise ValueError(f"shape mismatch: {x.shape} != {y.shape}")
    if not (
        x.flags.c_contiguous
        and y.flags.c_contiguous
        or x.flags.f_contiguous
        and y.flags.f_contiguous
    ):
        raise ValueError("x and y must be contiguous arrays")

    size = x.size
    alpha = complex(alpha)
    x_ptr = ffi.from_buffer("ls_hs_scalar[]", x, require_writable=False)
    y_ptr = ffi.from_buffer("ls_hs_scalar[]", y, require_writable=True)
    lib.ls_chpl_experimental_axpy_c128(size, alpha.real, alpha.imag, x_ptr, y_ptr, y_ptr)


def _csr_matvec(
    matrix: csr_matrix,
    x: NDArray,
    out: None | NDArray = None,
    num_tasks: int = -1,
) -> NDArray:
    if np.iscomplexobj(matrix.data) or np.iscomplexobj(x):
        elt_dtype = np.dtype("complex128")
        elt_dtype_str = "ls_hs_scalar[]"
        elt_suffix = "c128"
    else:
        elt_dtype = np.dtype("float64")
        elt_dtype_str = "double[]"
        elt_suffix = "f64"

    if matrix.indices.dtype == np.dtype("int64"):
        idx_dtype = np.dtype("int64")
        idx_dtype_str = "int64_t[]"
        idx_suffix = "i64"
    else:
        idx_dtype = np.dtype("int32")
        idx_dtype_str = "int32_t[]"
        idx_suffix = "i32"

    data = np.require(matrix.data, dtype=elt_dtype, requirements=["C"])
    indptr = np.require(matrix.indptr, dtype=idx_dtype, requirements=["C"])
    indices = np.require(matrix.indices, dtype=idx_dtype, requirements=["C"])
    matrix_elements = ffi.from_buffer(elt_dtype_str, data, require_writable=False)
    row_offsets = ffi.from_buffer(idx_dtype_str, indptr, require_writable=False)
    col_indices = ffi.from_buffer(idx_dtype_str, indices, require_writable=False)

    x = np.asarray(x, order="C", dtype=elt_dtype)
    if x.shape != (matrix.shape[1],):
        raise ValueError(f"'x' has invalid shape: {x.shape}; expected ({matrix.shape[1]},)")
    x_ptr = ffi.from_buffer(elt_dtype_str, x, require_writable=False)

    if out is None:
        out = np.empty(matrix.shape[0], dtype=elt_dtype)
    else:
        if out.shape != (matrix.shape[0],):
            raise ValueError(f"'out' has invalid shape: {out.shape}; expected ({matrix.shape[0]},)")
        if not out.flags.c_contiguous:
            raise ValueError("'out' must be contiguous")
        out = np.require(out, dtype=elt_dtype, requirements=["C"])
    y_ptr = ffi.from_buffer(elt_dtype_str, out, require_writable=True)

    func_name = f"ls_chpl_matrix_vector_product_csr_{idx_suffix}_{elt_suffix}"
    func = getattr(lib, func_name)

    func(
        matrix.shape[0],
        matrix.shape[1],
        matrix.nnz,
        matrix_elements,
        row_offsets,
        col_indices,
        x_ptr,
        y_ptr,
        num_tasks,
    )
    return out


#
#     # def apply_off_diag_to_basis_states(
#     #     op: ls.Operator, alphas: npt.NDArray[np.uint64]
#     # ) -> tuple[npt.NDArray[np.uint64], npt.NDArray[np.complex128], npt.NDArray[np.int64]]:
#     #     alphas = np.asarray(alphas, order="C", dtype=np.uint64)
#     #     alphas_ptr = ls.ffi.from_buffer("uint64_t[]", alphas, require_writable=False)
#     #
#     #     betas = ls.ffi.new("chpl_external_array *")
#     #     coeffs = ls.ffi.new("chpl_external_array *")
#     #     offsets = ls.ffi.new("chpl_external_array *")
#     #     kernels = ls.lib.ls_hs_internal_get_chpl_kernels()
#     #     kernels.operator_apply_off_diag(
#     #         op._payload, alphas.size, alphas_ptr, betas, coeffs, offsets, 0
#     #     )
#     #
#     #     offsets_arr = ls._chpl_external_array_as_ndarray(offsets, np.int64)
#     #     betas_arr = ls._chpl_external_array_as_ndarray(betas, np.uint64)[: offsets_arr[-1]]
#     #     coeffs_arr = ls._chpl_external_array_as_ndarray(coeffs, np.complex128)[: offsets_arr[-1]]
#     #     return betas_arr, coeffs_arr, offsets_arr
#     #
#     #
#     # def apply_diag_to_basis_states(
#     #     op: ls.Operator, alphas: npt.NDArray[np.uint64]
#     # ) -> npt.NDArray[np.float64]:
#     #     alphas = np.asarray(alphas, order="C", dtype=np.uint64)
#     #     alphas_ptr = ls.ffi.from_buffer("uint64_t[]", alphas, require_writable=False)
#     #
#     #     coeffs = ls.ffi.new("chpl_external_array *")
#     #     kernels = ls.lib.ls_hs_internal_get_chpl_kernels()
#     #     kernels.operator_apply_diag(op._payload, alphas.size, alphas_ptr, coeffs, 0)
#     #
#     #     coeffs_arr = ls._chpl_external_array_as_ndarray(coeffs, np.float64)
#     #     return coeffs_arr
#
#     @overload
#     def apply_diag_to_basis_state(self, alphas: int) -> float:
#         ...
#
#     @overload
#     def apply_diag_to_basis_state(self, alphas: NDArray[np.uint64]) -> NDArray[np.float64]:
#         ...
#
#     def apply_diag_to_basis_state(
#         self, alphas: int | NDArray[np.uint64]
#     ) -> float | NDArray[np.float64]:
#         if isinstance(alphas, (int, np.integer)):
#             is_scalar = True
#             count = 1
#             alphas_ptr = _basis_state_to_array(int(alphas), self.basis.number_words)
#         else:
#             is_scalar = False
#             alphas = np.asarray(alphas, order="C", dtype=np.uint64)
#             alphas_ptr = ffi.from_buffer("uint64_t[]", alphas, require_writable=False)
#             count = alphas.size
#
#         coeffs = ffi.new("chpl_external_array *")
#         kernels = lib.ls_hs_internal_get_chpl_kernels()
#         kernels.operator_apply_diag(self._payload, count, alphas_ptr, coeffs, 0)
#
#         coeffs_arr = _chpl_external_array_as_ndarray(coeffs, np.float64)
#         if is_scalar:
#             return float(coeffs_arr[0])
#         else:
#             return coeffs_arr
#
#     @overload
#     def apply_off_diag_to_basis_state(self, alphas: int) -> List[Tuple[complex, int]]:
#         ...
#
#     @overload
#     def apply_off_diag_to_basis_state(
#         self, alphas: NDArray[np.uint64]
#     ) -> Tuple[NDArray[np.uint64], NDArray[np.complex128], NDArray[np.int64]]:
#         ...
#
#     def apply_off_diag_to_basis_state(
#         self, alphas: int | NDArray[np.uint64]
#     ) -> (
#         List[Tuple[complex, int]]
#         | Tuple[NDArray[np.uint64], NDArray[np.complex128], NDArray[np.int64]]
#     ):
#         if isinstance(alphas, (int, np.integer)):
#             is_scalar = True
#             count = 1
#             alphas_ptr = _basis_state_to_array(int(alphas), self.basis.number_words)
#         else:
#             is_scalar = False
#             alphas = np.asarray(alphas, order="C", dtype=np.uint64)
#             alphas_ptr = ffi.from_buffer("uint64_t[]", alphas, require_writable=False)
#             count = alphas.size
#
#         betas = ffi.new("chpl_external_array *")
#         coeffs = ffi.new("chpl_external_array *")
#         offsets = ffi.new("chpl_external_array *")
#         kernels = lib.ls_hs_internal_get_chpl_kernels()
#         kernels.operator_apply_off_diag(self._payload, count, alphas_ptr, betas, coeffs, offsets, 0)
#
#         offsets_arr = _chpl_external_array_as_ndarray(offsets, np.int64)
#         betas_arr = _chpl_external_array_as_ndarray(betas, np.uint64)[: offsets_arr[-1]]
#         coeffs_arr = _chpl_external_array_as_ndarray(coeffs, np.complex128)[: offsets_arr[-1]]
#         if is_scalar:
#             return list(zip(coeffs_arr, betas_arr))
#         else:
#             return betas_arr, coeffs_arr, offsets_arr
#
#     def to_csr(self):
#         row_offsets = ffi.new("chpl_external_array *")
#         col_indices = ffi.new("chpl_external_array *")
#         matrix_elements = ffi.new("chpl_external_array *")
#         kernels = lib.ls_hs_internal_get_chpl_kernels()
#         kernels.operator_to_csr(self._payload, row_offsets, col_indices, matrix_elements, 0)
#
#         offsets_arr = _chpl_external_array_as_ndarray(row_offsets, np.int32)
#         indices_arr = _chpl_external_array_as_ndarray(col_indices, np.int32)
#         coeffs_arr = _chpl_external_array_as_ndarray(matrix_elements, np.complex128)
#         size = offsets_arr[-1]
#         dim = self.basis.number_states
#
#         return scipy.sparse.csr_matrix(
#             (coeffs_arr[:size], indices_arr[:size], offsets_arr), shape=(dim, dim)
#         )
#
#     def _check_basis_is_built(self, attribute):
#         if not self.basis.is_built:
#             raise AttributeError(
#                 "'Operator' object has no attribute '{}' (did you forget to build the basis?)"
#                 "".format(attribute),
#                 name=attribute,
#                 obj=self,
#             )
#
#     @property
#     def dtype(self):
#         self._check_basis_is_built("dtype")
#         return np.dtype("float64") if self.is_real else np.dtype("complex128")
#
#     @property
#     def shape(self):
#         self._check_basis_is_built("shape")
#         n = self.basis.number_states
#         return (n, n)
#
#     def _matvec(self, x):
#         return self.apply_to_state_vector(x)
#
#     def prepare_inputs_for_hphi(self, folder: str) -> None:
#         lib.ls_hs_prepare_hphi(self._payload, folder.encode("utf-8"))
#
#     def prepare_inputs_for_mvmc(self, folder: str) -> None:
#         lib.ls_hs_prepare_mvmc(self._payload, folder.encode("utf-8"))
#
#     def abelian_representations(self):
#         # with process_symmetries_impl_lock:
#         #     global process_symmetries_impl
#         #     representations = []
#         #     process_symmetries_impl = lambda p: representations.append(Symmetries(p))
#         return _from_json(lib.ls_hs_operator_abelian_representations(self._payload))
#
#
# def load_yaml_config(filename: str):
#     config = lib.ls_hs_load_yaml_config(filename.encode("utf-8"))
#     basis = Basis(payload=lib.ls_hs_clone_basis(config.basis))
#     hamiltonian = None
#     if config.hamiltonian != 0:
#         hamiltonian = Operator(lib.ls_hs_clone_operator(config.hamiltonian))
#     if config.number_observables != 0:
#         raise NotImplementedError
#     lib.ls_hs_destroy_yaml_config(config)
#     Config = namedtuple("Config", ["basis", "hamiltonian", "observables"], defaults=[None, None])
#     return Config(basis, hamiltonian)
#
#
# def matrix_vector_product_csr(
#     matrix: scipy.sparse.csr_matrix,
#     x: NDArray[np.complex128],
#     out: None | NDArray[np.complex128] = None,
# ):
#     data = np.require(matrix.data, dtype=np.complex128, requirements=["C"])
#     indptr = np.require(matrix.indptr, dtype=np.int32, requirements=["C"])
#     indices = np.require(matrix.indices, dtype=np.int32, requirements=["C"])
#     matrix_elements = ffi.from_buffer("ls_hs_scalar[]", matrix.data, require_writable=False)
#     row_offsets = ffi.from_buffer("int32_t[]", matrix.indptr, require_writable=False)
#     col_indices = ffi.from_buffer("int32_t[]", matrix.indices, require_writable=False)
#
#     x = np.asarray(x, order="C", dtype=np.complex128)
#     x_ptr = ffi.from_buffer("ls_hs_scalar[]", x, require_writable=False)
#     if out is None:
#         out = np.empty_like(x)
#     else:
#         out = np.asarray(out, order="C", dtype=np.complex128)
#     y_ptr = ffi.from_buffer("ls_hs_scalar[]", out, require_writable=True)
#
#     kernels = lib.ls_hs_internal_get_chpl_kernels()
#     kernels.matrix_vector_product_csr_i32_c128(
#         matrix.shape[0],
#         matrix.shape[1],
#         matrix.nnz,
#         matrix_elements,
#         row_offsets,
#         col_indices,
#         x_ptr,
#         y_ptr,
#         -1,
#     )
#     return out
