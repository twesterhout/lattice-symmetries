import dataclasses
from typing import Optional
import numpy as np
from numpy.typing import NDArray
from scipy.special import comb
from sympy import Rational
from sympy.combinatorics import Permutation, PermutationGroup

from lattice_symmetries._kernels import (
    BasisInfo,
    CompiledKernel,
    is_representative_kernel,
    xored_state_to_index_kernel,
)
from lattice_symmetries import _ls


def trailing_zeros(x):
    return (x & -x).bit_length() - 1


def fixed_hamming_state_to_index(state: int) -> int:
    state = int(state)
    if state < 0:
        raise ValueError(f"invalid 'state': {state}; expected a non-negative number")
    k = 0
    index = 0
    while state > 0:
        n = trailing_zeros(state)
        k += 1
        if k <= n:
            index += comb(n, k, exact=True)
        state &= state - 1
    return index


def fixed_hamming_index_to_state(index: int, *, number_sites: int, hamming_weight: int) -> int:
    state = 0
    k = hamming_weight
    for n in range(number_sites, 0, -1):
        state <<= 1
        current = 0 if k > n - 1 else comb(n - 1, k)
        if index >= current:
            index -= current
            k -= 1
            state |= 1
    return state


@dataclasses.dataclass
class BasisKernels:
    is_representative: Optional[CompiledKernel] = None
    state_info: Optional[CompiledKernel] = None
    state_to_index: Optional[CompiledKernel] = None


def enumerate_basis_states(
    info: BasisInfo, kernels: BasisKernels | None = None
) -> tuple[NDArray[np.uint64], NDArray[np.uint16]]:
    if info.is_state_index_identity:
        if info.spin_inversion is None:
            max_state = 2**info.number_bits
        else:
            # If spin_inversion is not None, leave the most significant bit as 0
            max_state = 2 ** (info.number_bits - 1)
        basis_states = np.arange(max_state, dtype=np.uint64)
        norms = np.ones(basis_states.size, dtype=np.uint16)
    elif not info.has_permutation_symmetries:
        assert info.hamming_weight is not None
        if info.hamming_weight > info.number_bits:
            msg = f"Hamming weight {info.hamming_weight} exceeds the number of bits {info.number_bits}"
            raise ValueError(msg)

        min_state = 2**info.hamming_weight - 1
        if info.spin_inversion is None:
            max_state = min_state << (info.number_bits - info.hamming_weight)
        else:
            max_state = min_state << (info.number_bits - 1 - info.hamming_weight)

        assert fixed_hamming_state_to_index(min_state) == 0
        number_states = fixed_hamming_state_to_index(max_state) + 1

        basis_states = np.empty(number_states, dtype=np.uint64)
        offsets = np.array([0, number_states], dtype=np.int64)
        values = np.array([min_state], dtype=np.uint64)
        _ls.lib.ls_enumerate_states_fixed_hamming(
            1,
            _ls.ffi.from_buffer("const int64_t*", offsets, require_writable=False),
            _ls.ffi.from_buffer("const uint64_t*", values, require_writable=False),
            _ls.ffi.from_buffer("uint64_t*", basis_states, require_writable=True),
        )

        # All norms are 1, because
        # - either no projection takes place (i.e., info.spin_inversion is None)
        # - or we apply spin inversion, but x is always not equal to invert(x)
        norms = np.ones(number_states, dtype=np.uint16)
    else:
        raise NotImplementedError("😭")

    basis_states.flags.writeable = False
    norms.flags.writeable = False
    return basis_states, norms


@dataclasses.dataclass(frozen=True)
class Basis:
    info: BasisInfo
    kernels: BasisKernels = dataclasses.field(default_factory=BasisKernels)

    states: NDArray[np.uint64] | None = None
    norms: NDArray[np.uint8] | None = None

    @property
    def hamming_weight(self) -> Optional[int]:
        return self.info.hamming_weight

    @property
    def spin_inversion(self) -> Optional[int]:
        return self.info.spin_inversion

    @property
    def number_bits(self) -> int:
        return self.info.number_bits

    @property
    def number_words(self) -> int:
        return (self.number_bits + 7) // 8

    @property
    def symmetries(self) -> list[tuple[Permutation, Rational]]:
        return self.info.symmetries

    @property
    def has_permutation_symmetries(self) -> bool:
        return self.info.has_permutation_symmetries

    @property
    def is_built(self) -> bool:
        return self.states is not None

    def check_is_built(self):
        if not self.is_built:
            msg = "basis states have not been built yet; if you wish to do so, use the basis.build() function"
            raise ValueError(msg)

    def build(self) -> None:
        """Generate a list of representatives.

        These can later be accessed using the `number_states` and `states` attributes.
        """
        if not self.is_built:
            states, norms = enumerate_basis_states(self.info, self.kernels)
            object.__setattr__(self, "states", states)
            object.__setattr__(self, "norms", norms)
        assert self.is_built

    @property
    def number_states(self) -> int:
        self.check_is_built()
        return int(self.states.size)

    def is_representative(self, x: int | NDArray[np.uint64]) -> int | NDArray[np.uint16]:
        is_scalar = False
        x = np.asarray(x, dtype=np.uint64, order="C")
        if x.ndim == 0:
            is_scalar = True
            x = np.expand_dims(x, axis=0)
        if x.ndim != 1:
            raise ValueError(f"'x' has invalid shape: {x.shape}; expected a one-dimensional array")

        if self.kernels.is_representative is None:
            self.kernels.is_representative = is_representative_kernel(self.info)

        norms = np.empty(x.size, dtype=np.uint16)
        self.kernels.is_representative.callable(
            # NOTE: temporary hack until Halide start supporting uint64 NumPy arrays
            x.view(np.int64),
            norms,
        )

        if is_scalar:
            return int(norms[0])
        else:
            return norms

    def index(self, x: int | NDArray[np.uint64]) -> int | NDArray[np.int64]:
        """Return the index of a basis state or a batch of basis states."""
        if self.number_bits > 64:
            msg = "it is impractical to compute indices of states with more than 64 bits"
            raise ValueError(msg)

        is_scalar = False
        x = np.asarray(x, dtype=np.uint64, order="C")
        if x.ndim == 0:
            is_scalar = True
            x = np.expand_dims(x, axis=0)
        if x.ndim != 1:
            raise ValueError(f"'x' has invalid shape: {x.shape}; expected a one-dimensional array")

        if self.kernels.state_to_index is None:
            self.kernels.state_to_index = xored_state_to_index_kernel(self.info)

        indices = np.empty((1, x.size), dtype=np.int64)
        masks = np.zeros(1, dtype=np.int64)
        basis_states = self.states.view(np.int64) if self.is_built else np.empty(1, dtype=np.int64)
        self.kernels.state_to_index.callable(
            # NOTE: temporary hack until Halide start supporting uint64 NumPy arrays
            x.view(np.int64),
            masks,
            basis_states,
            indices,
        )

        if is_scalar:
            return int(indices[0, 0])
        else:
            return indices[0]


class SpinBasis(Basis):
    def __init__(
        self,
        number_spins: int,
        hamming_weight: int | None = None,
        spin_inversion: int | None = None,
    ):
        if number_spins < 0:
            raise ValueError(f"invalid number_spins={number_spins}")
        if hamming_weight is not None and (hamming_weight < 0 or hamming_weight > number_spins):
            raise ValueError(f"invalid hamming_weight={hamming_weight}")
        if spin_inversion is not None and spin_inversion != 1 and spin_inversion != -1:
            raise ValueError(f"invalid spin_inversion={spin_inversion}")
        if spin_inversion is not None and hamming_weight is not None:
            if 2 * hamming_weight != number_spins:
                msg = f"incompatible spin_inversion={spin_inversion} and hamming_weight={hamming_weight}"
                raise ValueError(msg)
        info = BasisInfo(
            number_bits=number_spins, hamming_weight=hamming_weight, spin_inversion=spin_inversion
        )
        super().__init__(info)

    @property
    def number_spins(self) -> int:
        return self.info.number_bits

    def state_to_string(self, state: int) -> str:
        return ("|{:0" + str(self.info.number_bits) + "b}⟩").format(state)
