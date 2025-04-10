import dataclasses, numpy as np, time, lattice_symmetries as ls
from numpy.typing import NDArray
from sympy import Rational
from sympy.combinatorics import Permutation

@dataclasses.dataclass(frozen=True)
class B:
    i: any
    states: NDArray[np.uint64] | None = None
    norms: NDArray[np.uint16] | None = None
    bs_ctx: any = None; search_ctx: any = None

    hamming_weight = property(lambda self: self.i.hamming)
    spin_inversion = property(lambda self: self.i.inversion)
    number_bits = property(lambda self: self.i.bits)
    symmetries = property(lambda self: self.i.symmetries)
    has_permutation_symmetries = property(lambda self: self.i.has_ps)
    is_built = property(lambda self: self.states is not None)
    number_states = property(lambda self: self._check_is_built() or int(self.states.size))

    def index(self, x: int | NDArray[np.uint64]):
        """Return the index of a basis state or a batch of basis states."""
        if self.number_bits > 64: raise ValueError("it is impractical to compute indices of states with more than 64 bits")
        is_scalar, x = False, np.asarray(x, dtype=np.uint64, order="C")
        if x.ndim == 0: is_scalar, x = True, np.expand_dims(x, axis=0)
        if x.ndim != 1: raise ValueError(f"'x' has invalid shape: {x.shape}; expected a one-dimensional array")
        if self.i.is_s2i_id:
            indices = x
        elif self.i.has_ps and self.i.hamming is None:
            self._prepare_search_ctx()
            indices = ls.compiler.state_to_index(x, self.search_ctx)
        else:
            raise NotImplementedError()
        return int(indices[0]) if is_scalar else indices

    def state_info(self, x: int | NDArray[np.uint64]):
        is_scalar, x = False, np.asarray(x, dtype=np.uint64, order="C")
        if x.ndim == 0: is_scalar, x = True, np.expand_dims(x, axis=0)
        if x.ndim != 1: raise ValueError(f"'x' has invalid shape: {x.shape}; expected a one-dimensional array")
        if self.i.is_s2i_id:
            rep, idx = x, np.full(x.size, fill_value=-1, dtype=np.int64)
        elif self.i.has_ps and self.i.hamming is None:
            self._prepare_bs_ctx()
            rep, idx = ls.compiler.state_info(x, self.bs_ctx)
        else:
            raise NotImplementedError()
        return (rep[0], idx[0]) if is_scalar else (rep, idx)

    def build(self):
        if not self.is_built:
            bs_ctx = ls.compiler.bs_ctx_t(self.i)
            states, norms = ls.compiler.enumerate_states(self.i, ctx=bs_ctx)
            object.__setattr__(self, "states", states)
            object.__setattr__(self, "norms", norms)
            object.__setattr__(self, "bs_ctx", bs_ctx)

    def _prepare_bs_ctx(self):
        if self.bs_ctx is None: object.__setattr__(self, "bs_ctx", ls.compiler.bs_ctx_t(self.i))
    def _prepare_search_ctx(self):
        if self.search_ctx is None:
            if not self.i.is_s2i_id:
                self._check_is_built()
                ctx = ls.compiler.search_ctx_t(self.i, self.states, self.norms)
            else:
                ctx = ls.compiler.Ctx(ls.compiler.NULL, None)
            object.__setattr__(self, "search_ctx", ctx)
                
    def _check_is_built(self):
        if not self.is_built: raise ValueError("basis states have not been built yet; if you wish to do so, use the basis.build() function")
Basis = B

    
class SpinBasis(Basis):
    def __init__(self, number_spins: int, hamming_weight: int | None = None, spin_inversion: int | None = None,
            symmetries: list[tuple[Permutation, Rational]] = []):
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
        super().__init__(ls.BasisInfo(bits=number_spins, hamming=hamming_weight, inversion=spin_inversion,
            symmetries=symmetries))
    number_spins = property(lambda self: self.i.bits)
    def state_to_string(self, state: int) -> str: return ("|{:0" + str(self.i.bits) + "b}⟩").format(state)
