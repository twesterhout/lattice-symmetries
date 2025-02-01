import dataclasses
import functools
import contextlib
import itertools
import os
import tempfile
import subprocess
import time
from typing import Callable
from lattice_symmetries import _ls
from lattice_symmetries._benes import BenesNetwork, permutation_to_benes_network
from lattice_symmetries._representation import generate_representation

import cffi
import halide as hl
import numpy as np
from numpy.typing import NDArray
from loguru import logger
from scipy.special import comb
import sympy
from sympy.combinatorics import Permutation
from sympy import Rational

_KERNEL_DEFINITIONS = """
  struct halide_device_interface_t;

  typedef enum halide_type_code_t
  {
      halide_type_int = 0,     ///< signed integers
      halide_type_uint = 1,    ///< unsigned integers
      halide_type_float = 2,   ///< IEEE floating point numbers
      halide_type_handle = 3,  ///< opaque pointer type (void *)
      halide_type_bfloat = 4,  ///< floating point numbers in the bfloat format
  } halide_type_code_t;

  typedef struct halide_dimension_t {
    int32_t min, extent, stride;
    uint32_t flags;
  } halide_dimension_t;

  struct halide_type_t {
    uint8_t code;
    uint8_t bits;
    uint16_t lanes;
  };

  typedef struct halide_buffer_t {
    uint64_t device;
    const struct halide_device_interface_t *device_interface;
    uint8_t *host;
    uint64_t flags;
    struct halide_type_t type;
    int32_t dimensions;
    halide_dimension_t *dim;
    void *padding;
  } halide_buffer_t;

  int is_representative_kernel(struct halide_buffer_t *, struct halide_buffer_t *);

  int xored_state_to_index_kernel(struct halide_buffer_t *alphas, struct halide_buffer_t *mask, struct halide_buffer_t *basis_states, struct halide_buffer_t *indices);
  int state_to_index_kernel(struct halide_buffer_t *alphas, struct halide_buffer_t *basis_states, struct halide_buffer_t *indices);

  int diag_matrix_complex_kernel(struct halide_buffer_t *, struct halide_buffer_t *, struct halide_buffer_t *);
  int diag_matrix_real_kernel(struct halide_buffer_t *, struct halide_buffer_t *);

  int off_diag_matrix_kernel(struct halide_buffer_t *, struct halide_buffer_t *, struct halide_buffer_t *, struct halide_buffer_t *, struct halide_buffer_t *);
"""


@contextlib.contextmanager
def measure_time():
    tick = tock = time.perf_counter() 
    yield lambda: tock - tick
    tock = time.perf_counter() 


@dataclasses.dataclass(frozen=True)
class CompiledKernel:
    ffi_fun_ptr: any
    callable: any


@dataclasses.dataclass(frozen=True)
class BasisInfo:
    number_bits: int
    hamming_weight: int | None = None
    spin_inversion: int | None = None
    symmetries: list[tuple[Permutation, Rational]] = dataclasses.field(default_factory=list)

    @property
    def has_permutation_symmetries(self) -> bool:
        return len(self.symmetries) > 0

    @property
    def is_state_index_identity(self) -> bool:
        return self.hamming_weight is None and not self.has_permutation_symmetries

    @property
    def min_and_max_state_estimate(self) -> int:
        if self.hamming_weight is None:
            min_state = 0
            if self.spin_inversion is None:
                max_state = 2**self.number_bits
            else:
                # If spin_inversion is not None, leave the most significant bit as 0
                max_state = 2 ** (self.number_bits - 1)
        else:
            min_state = 2**self.hamming_weight - 1
            if self.spin_inversion is None:
                max_state = min_state << (self.number_bits - self.hamming_weight)
            else:
                max_state = min_state << (self.number_bits - 1 - self.hamming_weight)
        return min_state, max_state


@dataclasses.dataclass(frozen=True)
class StateToIndexInfo:
    offsets: NDArray[np.int64]
    shift: int
    prefix_bits: int
    range_size: int


class LoweredSymmetries:
    masks: NDArray[np.uint64]
    shifts: NDArray[np.uint64]
    characters_re: NDArray[np.float64]
    characters_im: NDArray[np.float64]
    networks: list[BenesNetwork]

    def __init__(self, representation: list[tuple[Permutation, Rational]]):
        if len(representation) == 0:
            self.masks = None
            self.shifts = None
            self.characters_re = None
            self.characters_im = None
            self.networks = None
        else:
            with measure_time() as t:
                self.networks = [permutation_to_benes_network(p) for p, _ in representation]
            msg = f"permutation_to_benes_network took {t()} seconds; {t() / len(representation)} per permutation"
            logger.debug(msg)

            self.shifts = np.asarray(self.networks[0].shifts, dtype=np.uint64)
            for i, b in enumerate(self.networks[1:], 1):
                if not np.array_equal(b.shifts, self.shifts):
                    msg = f"incompatible symmetries: {representation[0][0]} and {representation[i][0]}"
                    raise ValueError(msg)
            self.masks = np.vstack([np.asarray(b.masks, dtype=np.uint64) for b in self.networks])
            characters = [sympy.exp(-2 * sympy.pi * sympy.I * r) for _, r in representation]
            self.characters_re = np.asarray([sympy.re(c) for c in characters], dtype=np.float64)
            self.characters_im = np.asarray([sympy.im(c) for c in characters], dtype=np.float64)


class LoweredOperator:
    info: BasisInfo
    terms: PauliLoweredTerms
    kernel: CompiledKernel

    def __init__(self, info: BasisInfo, terms: list[PauliNonbranchingTerm], symm: LoweredSymmetries | None, state_to_index_info: StateToIndexInfo | None, verbose: bool = False):
        self.info = info
        self.terms = PauliLoweredTerms(terms)
        # print(self.terms)
        # self.diag_kernel = diag_matrix_kernel(self.terms) if self.terms.s_diag.size > 0 else None
        self.kernel = off_diag_matrix_kernel(info=self.info, terms=self.terms, symm=symm, state_to_index_info=state_to_index_info, verbose=verbose)


        # diag_data = _ls.ffi.new("ls_diag_terms *")
        # if self.diag_kernel is None:
        #     diag_data.kernel = _ls.ffi.NULL
        # else:
        #     diag_data.kernel = self.diag_kernel.ffi_fun_ptr
        # diag_data.number_terms = self.terms.s_diag.size
        # diag_data.v_re = _ls.ffi.from_buffer("const double *", self.terms.v_re_diag)
        # diag_data.v_im = _ls.ffi.from_buffer("const double *", self.terms.v_im_diag)
        # self.diag_data = diag_data

        # off_diag_data = _ls.ffi.new("ls_off_diag_terms *")
        # if self.kernel is None:
        #     off_diag_data.kernel = _ls.ffi.NULL
        # else:
        #     off_diag_data.kernel = self.kernel.ffi_fun_ptr
        # off_diag_data.number_terms = self.terms.s_2d.shape[0]
        # off_diag_data.number_reduced = self.terms.s_2d.shape[1]
        # off_diag_data.v_re = _ls.ffi.from_buffer("const double *", self.terms.v_re_2d)
        # off_diag_data.v_im = _ls.ffi.from_buffer("const double *", self.terms.v_im_2d)
        # off_diag_data.x = _ls.ffi.from_buffer("const uint64_t *", self.terms.mask)
        # self.data = off_diag_data

    def apply(self, states: NDArray, norms: NDArray, x: NDArray, out: NDArray):
        has_diag = len(self.terms.s_diag) > 0
        has_off_diag = self.terms.s_2d.shape[0] > 0

        if not has_diag and not has_off_diag:
            out[:] = 0
            return

        dtype = np.dtype("float32")
        states = np.asarray(states, dtype=np.uint64, order="C")
        norms = np.asarray(norms, dtype=np.uint16, order="C")
        x = np.asarray(x, dtype=dtype, order="F")
        out = np.asarray(out, dtype=dtype, order="F")

        args = [states.view(np.int64), norms, x]
        if has_off_diag:
            args += [self.terms.v_re_2d, self.terms.v_im_2d]
        if has_diag:
            args += [self.terms.v_re_diag, self.terms.v_im_diag]
        args += [states.view(np.int64), norms, x, out]

        self.kernel.callable(*args)

        # _ls.lib.ls_matrix_apply_c128(
        #     self.diag_data,
        #     self.off_diag_data,
        #     states.size,
        #     _ls.ffi.from_buffer("const uint64_t *", states),
        #     x.shape[1],
        #     _ls.ffi.from_buffer("const _complex128 *", x),
        #     _ls.ffi.from_buffer("_complex128 *", out),
        # )
        # states = np.asarray(states, dtype=np.uint64, order="C")
        # buf_re = np.zeros((self.off_diag_data.number_terms, states.size), dtype=np.float64)
        # buf_im = np.zeros((self.off_diag_data.number_terms, states.size), dtype=np.float64)
        # self.off_diag_kernel.callable(states.view(np.int64), self.terms.v_re_2d, self.terms.v_im_2d, buf_re, buf_im)
        # beta = states.reshape(-1, 1) ^ self.terms.x.reshape(1, -1)
        # out[:] = np.einsum("jb,bjk->bk", buf_re + 1j * buf_im, x[beta])


class KernelCompiler:
    temp_dir: str
    ffi: any
    target: any
    cc: str
    _runtime: any

    def __init__(self, temp_dir=None):
        if temp_dir is None:
            self.temp_dir = tempfile.mkdtemp(prefix="lattice-symmetries-cache")
        logger.trace(f"'{self.temp_dir}' will be used for compiling kernels.")
        self.ffi = cffi.FFI()
        self.ffi.cdef(_KERNEL_DEFINITIONS)
        self.target = hl.get_jit_target_from_environment()
        self.cc = "cc"
        self._runtime = None

    def _link(self, library_file, *object_files):
        subprocess.run([self.cc, "-shared", "-o", library_file] + list(object_files), check=True)

    def _generate_and_load_runtime(self):
        object_file = os.path.join(self.temp_dir, "runtime.o")
        library_file = os.path.splitext(object_file)[0] + ".so"
        hl.compile_standalone_runtime(os.path.join(self.temp_dir, "runtime.o"), target=self.target)
        self._link(library_file, object_file)
        self._runtime = self.ffi.dlopen(library_file, self.ffi.RTLD_NOW | self.ffi.RTLD_GLOBAL)

    def link_and_load(self, builder: Callable, func_name: str):
        self.ffi.init_once(self._generate_and_load_runtime, "runtime")

        tick = time.perf_counter()
        func, params, keep_alive = builder()
        tock = time.perf_counter()
        logger.trace(f"Prepared {func_name} in {tock - tick} seconds.")

        tick = time.perf_counter()
        callable = func.compile_to_callable(params, target=self.target)
        tock = time.perf_counter()
        logger.trace(f"Compiled {func_name} to a Halide::Callable in {tock - tick} seconds.")

        tick = time.perf_counter()
        _, object_file = tempfile.mkstemp(suffix=".o", dir=self.temp_dir)
        func.compile_to_object(object_file, params, fn_name=func_name, target=self.target)
        tock = time.perf_counter()
        logger.trace(f"Compiled {func_name} to {object_file} in {tock - tick} seconds.")

        tick = time.perf_counter()
        library_file = os.path.splitext(object_file)[0] + ".so"
        self._link(library_file, object_file)
        tock = time.perf_counter()
        logger.trace(f"Linked {object_file} to {library_file} in {tock - tick} seconds.")

        ffi_lib = self.ffi.dlopen(library_file, self.ffi.RTLD_NOW | self.ffi.RTLD_LOCAL)
        ffi_fun_ptr = self.ffi.gc(getattr(ffi_lib, func_name), lambda p: self.ffi.dlclose(ffi_lib))

        return CompiledKernel(callable=callable, ffi_fun_ptr=ffi_fun_ptr)


COMPILER = KernelCompiler()


@functools.lru_cache(maxsize=16)
def binomials_buffer(size: int = 65) -> hl.Buffer:
    binomial_coeffs = [comb(*t, exact=True) for t in itertools.product(range(size), range(size))]
    # Halide does not seem to creating Buffers from NumPy arrays of uint64, so we use int64 instead...
    binomial_coeffs = np.asarray(binomial_coeffs, dtype=np.int64).reshape(size, size)
    binomial_coeffs = np.ascontiguousarray(binomial_coeffs.T)
    return hl.Buffer(binomial_coeffs)


"""
def _build_fixed_hamming_state_to_index_kernel(
    number_sites: int,
    hamming_weight: int,
    vector_size: int = 8,
    unroll: bool = True,
    target=hl.get_jit_target_from_environment(),
    verbose: bool = False,
):
    alpha = hl.ImageParam(hl.Int(64), 1, "alpha")
    mask = hl.ImageParam(hl.Int(64), 1, "mask")
    out = hl.Func("out")
    i_inner = hl.Var("i_inner")
    i_outer = hl.Var("i_outer")
    i = hl.Var("i")
    j = hl.Var("j")
    fused = hl.Var("fused")
    keep_alive = dict()

    if hamming_weight == 0 or hamming_weight == number_sites:
        out[i, j] = hl.cast(hl.Int(64), 0)
        out.fuse(i, j, fused)
        out.split(fused, i_outer, i_inner, vector_size, hl.TailStrategy.GuardWithIf)
        out.vectorize(i_inner)
    else:
        binomials = binomials_buffer()
        k = hl.RDom([hl.Range(0, hamming_weight)], "k")

        temp = hl.Func("temp")
        temp[i, j] = (hl.cast(hl.Int(64), 0), hl.cast(hl.UInt(64), alpha[i] ^ mask[j]))
        index, state = temp[i, j][0], temp[i, j][1]
        n = hl.cast(hl.Int(32), hl.count_trailing_zeros(state))
        temp[i, j] = (index + binomials[n, k + 1], state & (state - 1))
        out[i, j] = temp[i, j][0]

        out.split(i, i_outer, i_inner, vector_size, hl.TailStrategy.GuardWithIf)
        out.reorder(i_inner, i_outer, j)
        out.vectorize(i_inner)

        temp.split(i, i_outer, i_inner, vector_size, hl.TailStrategy.GuardWithIf)
        temp.vectorize(i_inner)
        temp.update(0).split(i, i_outer, i_inner, vector_size, hl.TailStrategy.GuardWithIf)
        temp.update(0).reorder(i_inner, k, i_outer, j)
        temp.update(0).vectorize(i_inner)
        # if unroll:
        #     temp.update(0).unroll(k)

        temp.compute_at(out, i_outer)
        temp.store_at(out, i_outer)
        temp.bound_storage(i, vector_size)
        temp.bound_storage(j, 1)
        temp.store_in(hl.MemoryType.Register)

        # NOTE: it's important to keep binomials alive until we construct a callable from out.
        keep_alive["binomials"] = binomials

    alpha.dim(0).set_min(0).set_stride(1)
    mask.dim(0).set_min(0).set_stride(1)
    out.output_buffer().dim(0).set_min(0).set_stride(1).set_extent(alpha.dim(0).extent())
    out.output_buffer().dim(1).set_min(0).set_stride(alpha.dim(0).extent()).set_extent(
        mask.dim(0).extent()
    )

    # A dummy argument
    basis_states = hl.ImageParam(hl.Int(64), 1, "basis_states")

    if verbose:
        out.print_loop_nest()
        out.compile_to_lowered_stmt(
            "fixed_hamming.stmt.html",
            [alpha, mask, basis_states],
            fmt=hl.StmtOutputFormat.HTML,
            target=target,
        )
    return out, [alpha, mask, basis_states], keep_alive
"""

"""
def _build_xored_state_to_index_kernel(info: BasisInfo, *args, **kwargs):
    if info.number_bits > 64:
        raise NotImplementedError("xored_state_to_index not supported on large systems")

    if info.has_permutation_symmetries:
        raise NotImplementedError("binary search kernel not yet implemented")
    elif info.hamming_weight is not None:
        return _build_fixed_hamming_state_to_index_kernel(
            info.number_bits, info.hamming_weight, *args, **kwargs
        )
    else:
        # Just an identity function
        alpha = hl.ImageParam(hl.Int(64), 1, "alpha")
        mask = hl.ImageParam(hl.Int(64), 1, "mask")
        out = hl.Func("out")
        i = hl.Var("i")
        j = hl.Var("j")
        i_inner = hl.Var("i_inner")
        i_outer = hl.Var("i_outer")
        keep_alive = dict()

        out[i, j] = alpha[i] ^ mask[j]
        vector_size = 16
        out.split(i, i_outer, i_inner, vector_size, hl.TailStrategy.GuardWithIf)
        out.vectorize(i_inner)

        # A dummy argument
        basis_states = hl.ImageParam(hl.Int(64), 1, "basis_states")

        alpha.dim(0).set_min(0).set_stride(1)
        mask.dim(0).set_min(0).set_stride(1)
        out.output_buffer().dim(0).set_min(0).set_stride(1).set_extent(alpha.dim(0).extent())
        out.output_buffer().dim(1).set_min(0).set_stride(alpha.dim(0).extent()).set_extent(
            mask.dim(0).extent()
        )
        return out, [alpha, mask, basis_states], keep_alive


def xored_state_to_index_kernel(info: BasisInfo, *args, **kwargs) -> CompiledKernel:
    builder = lambda: _build_xored_state_to_index_kernel(info, *args, **kwargs)
    return COMPILER.link_and_load(builder, "xored_state_to_index_kernel")
"""


def _bit_permute_step_64(x, m, d):
    y = ((x >> d) ^ x) & m
    return (x ^ y) ^ (y << d)


def _make_permuted(alphas, masks, shifts):
    batch_idx = hl.Var("batch_idx")
    group_idx = hl.Var("group_idx")
    permuted = hl.Func("permuted")

    p = alphas[batch_idx]
    for i in range(len(shifts)):
        p = _bit_permute_step_64(p, masks[i, group_idx], shifts[i])

    permuted[batch_idx, group_idx] = p
    return permuted


def _build_symmetric_is_representative_kernel(
    info: BasisInfo,
    target=hl.get_jit_target_from_environment(),
    vector_size=8,
    verbose: bool = False,
):
    symm = LoweredSymmetries(generate_representation(info.symmetries))
    masks = hl.Buffer(symm.masks.view(np.int64), name="masks")
    characters_re = hl.Buffer(symm.characters_re, name="characters")
    alphas = hl.ImageParam(hl.Int(64), 1, name="alphas")
    inversion_mask = 2**info.number_bits - 1 if info.spin_inversion is not None else None
    keep_alive = dict(masks=masks, characters_re=characters_re)

    number_masks = symm.masks.shape[0]
    # depth = symm.masks.shape[1]

    y = _make_permuted(alphas, masks, symm.shifts)
    (batch_idx, group_idx) = y.args()

    temp = hl.Func("temp")
    r_group = hl.RDom([hl.Range(1, number_masks - 1)], "r_group")

    init_n = hl.cast(hl.UInt(16), 1)
    if info.spin_inversion is not None:
        inverted = alphas[batch_idx] ^ inversion_mask
        is_greater = inverted > alphas[batch_idx]
        init_n = hl.cast(hl.UInt(16), is_greater)

    temp[batch_idx] = init_n
    current_n = temp[batch_idx]
    is_greater = y[batch_idx, r_group] > alphas[batch_idx]
    is_equal = y[batch_idx, r_group] == alphas[batch_idx]
    is_trivial = characters_re[r_group] == hl.cast(hl.Float(64), 1)
    next_n = hl.cast(hl.UInt(16), is_greater | (is_equal & is_trivial)) * (
        current_n + hl.cast(hl.UInt(16), is_equal)
    )
    if info.spin_inversion is not None:
        inverted = y[batch_idx, r_group] ^ inversion_mask
        is_greater = inverted > alphas[batch_idx]
        is_equal = inverted == alphas[batch_idx]
        is_trivial = characters_re[r_group] == hl.cast(hl.Float(64), info.spin_inversion)
        next_n = hl.cast(hl.UInt(16), is_greater | (is_equal & is_trivial)) * (
            next_n + hl.cast(hl.UInt(16), is_equal)
        )
    r_group.where(current_n > 0)
    temp[batch_idx] = next_n

    norm = hl.Func("norm")
    norm[batch_idx] = temp[batch_idx]

    outer = hl.Var("outer")
    inner = hl.Var("inner")
    norm.split(batch_idx, outer, inner, vector_size, hl.TailStrategy.GuardWithIf)
    norm.vectorize(inner)

    temp.split(batch_idx, outer, inner, vector_size, hl.TailStrategy.GuardWithIf)
    temp.vectorize(inner)
    temp.update(0).split(batch_idx, outer, inner, vector_size, hl.TailStrategy.GuardWithIf)
    temp.update(0).reorder(inner, r_group, outer)
    temp.update(0).vectorize(inner)
    temp.compute_at(norm, outer).store_at(norm, outer)
    temp.bound_storage(batch_idx, vector_size).store_in(hl.MemoryType.Register)

    if verbose:
        norm.print_loop_nest()
        target = target.with_feature(hl.TargetFeature.NoBoundsQuery)
        target = target.with_feature(hl.TargetFeature.NoAsserts)
        target = target.with_feature(hl.TargetFeature.AVX512_Zen4)
        norm.compile_to_lowered_stmt(
            "is_representative.stmt.html",
            [alphas],
            fmt=hl.StmtOutputFormat.HTML,
            target=target,
        )
    return norm, [alphas], keep_alive


def _build_is_representative_kernel(info: BasisInfo, *args, **kwargs):
    if info.number_bits > 64:
        raise NotImplementedError("is_representative not yet implemented for large systems")

    x = hl.ImageParam(hl.Int(64), 1, "x")
    norm = hl.Func("norm")
    i = hl.Var("i_inner")
    keep_alive = dict()

    if len(info.symmetries) > 0:
        return _build_symmetric_is_representative_kernel(info, *args, **kwargs)

    if info.hamming_weight is None and info.spin_inversion is None and len(info.symmetries) == 0:
        # Identity function
        norm[i] = hl.cast(hl.UInt(16), 1)
    else:
        raise NotImplementedError("😭")

    x.dim(0).set_min(0).set_stride(1)
    norm.output_buffer().dim(0).set_min(0).set_stride(1).set_extent(x.dim(0).extent())

    return norm, [x], keep_alive


def is_representative_kernel(info: BasisInfo, *args, **kwargs) -> CompiledKernel:
    builder = lambda: _build_is_representative_kernel(info, *args, **kwargs)
    return COMPILER.link_and_load(builder, "is_representative_kernel")





def generate_offset_ranges(representatives: np.ndarray, number_bits: int, shift: int) -> tuple[np.ndarray, int]:
    """
    Generate overlapping ranges for binary search based on most significant bits.
    
    Args:
        representatives: Sorted array of uint64 bit strings
        number_bits: Number of most significant bits to use for range splitting
        shift: Number of bits to shift right (usually total_bits - number_bits)
    
    Returns:
        offsets: Array of starting positions for each range
        range_size: Maximum size of any range
    """
    assert 0 <= number_bits < 64, "invalid number_bits"
    assert shift < 64, "invalid shift"
    if number_bits == 0:
        return np.array([0, len(representatives)], dtype=np.intp), len(representatives)
    
    # Create small array of boundary values we're searching for
    num_prefixes = 1 << number_bits
    # Create shifted boundary values: [0, 1<<shift, 2<<shift, ...]
    boundaries = (np.arange(num_prefixes, dtype=np.uint64) << shift)
    # Find positions where each prefix starts
    offsets = np.searchsorted(representatives, boundaries, side='left')
    offsets = np.append(offsets, len(representatives))
    sizes = np.diff(offsets)
    max_range_size = np.max(sizes)
    # Normalize ranges to have equal size
    number_states = len(representatives)
    mask = offsets[:-1] > (number_states - max_range_size)
    offsets[:-1][mask] = number_states - max_range_size
    return offsets, max_range_size






"""
def _build_state_to_index_binary_search_kernel(
    *,
    prefix_bits: int,
    offsets: NDArray[np.int64],
    range_size: int,
    shift: int,
    target=hl.get_jit_target_from_environment(),
    verbose: bool = False,
):
    alpha = hl.ImageParam(hl.Int(64), 1, "alpha")
    representatives = hl.ImageParam(hl.Int(64), 1, "representatives")
    offsets_buf = hl.Buffer(offsets, name="offsets")
    keep_alive = dict(offsets_buf=offsets_buf)

    batch_idx = hl.Var("batch_idx")
    needle = alpha[batch_idx]
    alpha_key = hl.cast(hl.Int(32), (needle >> shift) & (2**prefix_bits - 1))
    alpha_key = hl.unsafe_promise_clamped(alpha_key, 0, offsets_buf.dim(0).extent() - 1)
    base = hl.cast(hl.Int(32), offsets_buf[alpha_key])
    size = range_size
    number_memory_accesses = 0
    while size > 1:
        half = size // 2
        size -= half
        k = hl.unsafe_promise_clamped(base + half, 0, representatives.dim(0).extent() - 1)
        base = hl.select(representatives[k] < needle, k, base)
        number_memory_accesses += 1
    base = hl.unsafe_promise_clamped(base, 0, representatives.dim(0).extent() - 1)
    base = hl.select(representatives[base] < needle, base + 1, base)
    number_memory_accesses += 1
    base = hl.unsafe_promise_clamped(base, 0, representatives.dim(0).extent() - 1)

    result = hl.Func("result")
    result[batch_idx] = hl.select(representatives[base] == needle, base, -1)
    number_memory_accesses += 1

    if verbose:
        print(f"number_memory_accesses: {number_memory_accesses}; log2(range_size): {np.log2(range_size)}")
        result.print_loop_nest()
        target = target.with_feature(hl.TargetFeature.NoBoundsQuery)
        target = target.with_feature(hl.TargetFeature.NoAsserts)
        target = target.with_feature(hl.TargetFeature.AVX512_Zen4)
        result.compile_to_lowered_stmt(
            "state_to_index_binary_search.stmt.html",
            [alpha, representatives],
            fmt=hl.StmtOutputFormat.HTML,
            target=target,
        )
    return result, [alpha, representatives], keep_alive


def state_to_index_kernel(representatives: NDArray[np.uint64], prefix_bits: int = 17, **kwargs) -> CompiledKernel:
    total_bits = int(representatives.max()).bit_length()
    shift = max(0, total_bits - prefix_bits)
    offsets, range_size = generate_offset_ranges(representatives, prefix_bits, shift)

    builder = lambda: _build_state_to_index_binary_search_kernel(
        prefix_bits=prefix_bits, offsets=offsets, range_size=range_size, shift=shift, **kwargs
    )
    return COMPILER.link_and_load(builder, "state_to_index_kernel")
"""


"""
def _build_diag_matrix_kernel(
    terms: PauliLoweredTerms,
    real_only=False,
    vector_size=8,
    target=hl.get_jit_target_from_environment(),
):
    alpha = hl.ImageParam(hl.Int(64), 1, "alpha")
    v_re_buf = hl.ImageParam(hl.Float(64), 1, "v_re")
    if not real_only:
        v_im_buf = hl.ImageParam(hl.Float(64), 1, "v_im")
    s_buf = hl.Buffer(terms.s_diag.view(np.int64), name="s")

    i = hl.Var("i")
    k = hl.RDom([hl.Range(0, terms.s_diag.size)], "k")
    i_inner = hl.Var("i_inner")
    i_outer = hl.Var("i_outer")

    temp = hl.Func("temp")
    if real_only:
        temp[i] = hl.cast(hl.Float(64), 0)
    else:
        temp[i] = (hl.cast(hl.Float(64), 0), hl.cast(hl.Float(64), 0))
    sign = hl.popcount(alpha[i] & s_buf[k]) << 63
    re = hl.reinterpret(hl.Float(64), hl.reinterpret(hl.Int(64), v_re_buf[k]) ^ sign)
    if real_only:
        temp[i] = temp[i] + re
    else:
        im = hl.reinterpret(hl.Float(64), hl.reinterpret(hl.Int(64), v_im_buf[k]) ^ sign)
        temp[i] = (temp[i][0] + re, temp[i][1] + im)

    out = hl.Func("out")
    out[i] = temp[i]

    count = alpha.dim(0).extent()
    alpha.dim(0).set_min(0).set_stride(1)
    out.output_buffers()[0].dim(0).set_min(0).set_stride(1).set_extent(count)
    if not real_only:
        out.output_buffers()[1].dim(0).set_stride(1)

    out.split(i, i_outer, i_inner, vector_size, hl.TailStrategy.GuardWithIf)
    out.reorder(i_inner, i_outer)
    out.vectorize(i_inner)

    temp.split(i, i_outer, i_inner, vector_size, hl.TailStrategy.RoundUp)
    temp.reorder(i_inner, i_outer)
    temp.vectorize(i_inner)
    temp.update(0).split(i, i_outer, i_inner, vector_size, hl.TailStrategy.GuardWithIf)
    temp.update(0).reorder(i_inner, k, i_outer)
    temp.update(0).vectorize(i_inner)
    temp.update(0).unroll(k)

    temp.compute_at(out, i_outer)
    temp.store_at(out, i_outer)
    temp.bound_storage(i, vector_size)
    temp.store_in(hl.MemoryType.Register)

    # out.print_loop_nest()
    # out.compile_to_lowered_stmt(
    #     "diag_matrix.stmt.html",
    #     [alpha, v_re_buf, v_im_buf],
    #     fmt=hl.StmtOutputFormat.HTML,
    #     target=target,
    # )
    return out, [alpha, v_re_buf, v_im_buf], dict()


def diag_matrix_kernel(*args, real_only=False, **kwargs) -> CompiledKernel:
    builder = lambda: _build_diag_matrix_kernel(*args, real_only=real_only, **kwargs)
    fn = "diag_matrix_real_kernel" if real_only else "diag_matrix_complex_kernel"
    return COMPILER.link_and_load(builder, fn)
"""








def _build_gather(info: BasisInfo, symm: LoweredSymmetries | None, terms: PauliLoweredTerms, *, init_norms, init_coeffs, norms, X, diag_coeff, off_diag_coeff, representative, state_to_index, keep_alive: dict, real_only: bool, dtype: hl.Type):
    def _real(x):
        return x if real_only else x[0]

    gather = hl.Func("gather")
    (state_idx,) = diag_coeff.args()
    init_coeff = _real(diag_coeff[state_idx]) * init_coeffs[state_idx]
    gather[state_idx] = hl.cast(dtype, init_coeff)

    if info.has_permutation_symmetries:
        characters_re = hl.Buffer(symm.characters_re, name="characters_re")
        # characters_im = hl.Buffer(symm.characters.view(np.int64), name="characters_im")
        keep_alive["characters_re"] = characters_re
        # keep_alive["characters_im"] = characters_im
        number_characters = symm.characters_re.size

    number_terms = terms.mask.shape[0]
    number_states = X.dim(0).extent()
    if terms.has_off_diag:
        r_gather = hl.RDom([hl.Range(0, number_terms)], "r_gather")
        index = state_to_index[state_idx, r_gather]
        if index.size() > 1:
            index = index[0]
        out_of_bounds = index < 0

        if info.is_state_index_identity:
            index = hl.unsafe_promise_clamped(hl.cast(hl.Int(32), index), 0, number_states - 1)
        else:
            index = hl.clamp(hl.cast(hl.Int(32), index), 0, number_states - 1)

        coeff = _real(off_diag_coeff[state_idx, r_gather])
        if info.has_permutation_symmetries:
            group_idx = representative[state_idx, r_gather][1]
            group_idx = hl.unsafe_promise_clamped(group_idx, 0, number_characters - 1)
            coeff *= characters_re[group_idx] * hl.sqrt(hl.cast(hl.Float(64), norms[index]) / hl.cast(hl.Float(64), init_norms[state_idx]))

        coeff = hl.cast(dtype, coeff)
        if info.is_state_index_identity:
            coeff = coeff * X[index]
        else:
            coeff = hl.select(out_of_bounds, hl.cast(dtype, 0), coeff * X[index])
        gather[state_idx] = gather[state_idx] + coeff
    return gather

@dataclasses.dataclass
class DynamicArgs:
    alpha0: hl.ImageParam
    norms0: hl.ImageParam
    coeff0: hl.ImageParam
    v_2d_re: hl.ImageParam | None
    v_2d_im: hl.ImageParam | None
    v_diag_re: hl.ImageParam | None
    v_diag_im: hl.ImageParam | None
    basis_states: hl.ImageParam
    norms: hl.ImageParam
    X: hl.ImageParam
    def __init__(self, dtype=hl.Float(32), real_only=False):
        self.alpha0 = hl.ImageParam(hl.Int(64),  1, "alpha0")
        self.norms0 = hl.ImageParam(hl.UInt(16), 1, "norms0")
        self.coeff0 = hl.ImageParam(dtype,       1, "coeff0")
        self.v_2d_re = hl.ImageParam(hl.Float(64), 2, "v_2d_re")
        self.v_2d_im = hl.ImageParam(hl.Float(64), 2, "v_2d_im") if not real_only else None
        self.v_diag_re = hl.ImageParam(hl.Float(64), 1, "v_diag_re")
        self.v_diag_im = hl.ImageParam(hl.Float(64), 1, "v_diag_im") if not real_only else None
        self.basis_states = hl.ImageParam(hl.Int(64), 1, "basis_states")
        self.norms = hl.ImageParam(hl.UInt(16), 1, "norms")
        self.X = hl.ImageParam(dtype, 1, "X")

@dataclasses.dataclass
class StaticArgs:
    s_2d: hl.Buffer
    s_diag: hl.Buffer
    xors: hl.Buffer # t.x
    masks: hl.Buffer # Benes masks
    binomials: hl.Buffer

def _build_off_diag_coeff(terms, ctx, static, real_only: bool):
    f2i = lambda x: hl.reinterpret(hl.Int(64), x)
    i2f = lambda x: hl.reinterpret(hl.Float(64), x)
    bi, ti, re, im = hl.Var("bi"), hl.Var("ti"), 0, 0
    for r in range(terms.s_2d.shape[1]):
        sign = hl.popcount(ctx.alpha0[bi] & static.s_2d[r, ti]) << 63
        re += i2f(f2i(ctx.v_2d_re[r, ti]) ^ sign)
        if not real_only: im += i2f(f2i(ctx.v_2d_im[r, ti]) ^ sign)
    off_diag_coeff = hl.Func("off_diag_coeff"); off_diag_coeff[bi, ti] = re if real_only else (re, im)
    return off_diag_coeff

def _build_diag_coeff(terms, ctx, static, real_only: bool):
    bi, dr, z = hl.Var("bi"), hl.RDom([hl.Range(0, terms.s_diag.size)], "dr"), hl.cast(hl.Float(64), 0)
    diag_coeff = hl.Func("diag_coeff"); diag_coeff[bi] = z if real_only else (z, z)
    sign = hl.popcount(ctx.alpha0[bi] & static.s_diag[dr]) << 63
    re = i2f(f2i(ctx.v_diag_re[r]) ^ sign)
    if not real_only: im = i2f(f2i(ctx.v_diag_im[r]) ^ sign)
    acc = diag_coeff[bi]; diag_coeff[bi] = re + acc if real_only else (re + acc[0], im + acc[1])
    return diag_coeff

def _build_state_info(info, symm, terms, static):
    representative = hl.Func("representative")
    bi, ti = hl.Var("bi"), hl.Var("ti")
    if not info.has_permutation_symmetries:
        representative[bi, ti] = alpha[state_idx] ^ static.xors[ti]
        return representative
    raise NotImplementedError("😭")

def _build_state_to_index(info, state_to_index_info, ctx, static, representative):
    state_to_index = hl.Func("state_to_index")
    bi, ti = representative.args()
    if info.is_state_index_identity: state_to_index[bi, ti] = representative[bi, ti]; return state_to_index
    if not info.has_permutation_symmetries:
        assert info.hamming_weight is not None
        hr = hl.RDom([hl.Range(0, info.hamming_weight)], "hr")
        state_to_index[bi, ti] = (hl.u64(0), hl.u64(representative[bi, ti]))
        i, s = state_to_index[bi, ti][0], state_to_index[si, ti][1]
        tr = hl.i32(hl.select(s == 0, 0, hl.count_trailing_zeros(s)))
        state_to_index[bi, ti] = (i + static.binomials[tr, hr + 1], s & (s - 1))
        return state_to_index
    raise NotImplementedError("😭")
    # else:
    #     offsets = hl.Buffer(state_to_index_info.offsets, name="offsets")
    #     keep_alive["offsets"] = offsets

    #     number_states = basis_states.dim(0).extent()
    #     needle = representative[state_idx, term_idx][0]
    #     alpha_key = hl.cast(hl.Int(32), (needle >> state_to_index_info.shift) & (2**state_to_index_info.prefix_bits - 1))
    #     alpha_key = hl.unsafe_promise_clamped(alpha_key, 0, state_to_index_info.offsets.size - 1)
    #     base = hl.cast(hl.Int(32), offsets[alpha_key])
    #     size = state_to_index_info.range_size
    #     while size > 1:
    #         half = size // 2
    #         size -= half
    #         k = hl.unsafe_promise_clamped(base + half, 0, number_states - 1)
    #         base = hl.select(basis_states[k] < needle, k, base)
    #     base = hl.unsafe_promise_clamped(base, 0, number_states - 1)
    #     base = hl.select(basis_states[base] < needle, base + 1, base)
    #     base = hl.unsafe_promise_clamped(base, 0, number_states - 1)
    #     state_to_index[state_idx, term_idx] = hl.select(basis_states[base] == needle, base, -1)


def _build_off_diag_matrix_kernel(
    info: BasisInfo,
    state_to_index_info: StateToIndexInfo | None,
    symm: LoweredSymmetries | None,
    terms: PauliLoweredTerms,
    *,
    vector_size=16,
    task_size=1024,
    target=hl.get_jit_target_from_environment(),
    dtype=hl.Float(32),
    real_only=True,
    verbose=False,
):
    GuardWithIf = hl.TailStrategy.GuardWithIf
    RoundUp = hl.TailStrategy.RoundUp
    has_diag = terms.has_diag
    has_off_diag = terms.has_off_diag

    # Initial states and coefficients
    init_alpha_buf = hl.ImageParam(hl.Int(64), 1, "init_alpha")
    init_norms_buf = hl.ImageParam(hl.UInt(16), 1, "init_norms")
    init_coeffs_buf = hl.ImageParam(dtype, 1, "init_coeffs")
    # (Potentially) time-dependent matrix elements
    v_re_buf = hl.ImageParam(hl.Float(64), 2, "v_re")
    v_im_buf = hl.ImageParam(hl.Float(64), 2, "v_im")
    v_diag_re_buf = hl.ImageParam(hl.Float(64), 1, "v_diag_re")
    v_diag_im_buf = hl.ImageParam(hl.Float(64), 1, "v_diag_im")
    # Basis states, norms
    basis_states_buf = hl.ImageParam(hl.Int(64), 1, "basis_states")
    norms_buf = hl.ImageParam(hl.UInt(16), 1, "norms")
    # State vector
    X_buf = hl.ImageParam(dtype, 1, "X")

    s_buf = hl.Buffer(terms.s_2d.view(np.int64), name="s") if has_off_diag else None
    mask_buf = hl.Buffer(terms.mask.view(np.int64), name="mask") if has_off_diag else None
    s_diag_buf = hl.Buffer(terms.s_diag.view(np.int64), name="s_diag") if has_diag else None
    keep_alive = dict(s_buf=s_buf, mask_buf=mask_buf, s_diag_buf=s_diag_buf)

    # Compute off-diagonal matrix elements
    if terms.has_off_diag: off_diag_coeff = _build_off_diag_coeff(terms, ctx, static, real_only)
    if terms.has_diag: diag_coeff = _build_diag_coeff(terms, ctx, static, real_only)
    representative = _build_state_info(basis_info, symm, terms, static)
    state_to_index = _build_state_to_index(basis_info, state_to_index_info, ctx, static, representative)
    # Compute final state vector
    gather = _build_gather(
        info,
        symm,
        terms,
        init_norms=init_norms_buf,
        init_coeffs=init_coeffs_buf,
        norms=norms_buf,
        X=X_buf,
        diag_coeff=diag_coeff,
        off_diag_coeff=off_diag_coeff,
        representative=representative,
        state_to_index=state_to_index,
        keep_alive=keep_alive,
        real_only=real_only,
        dtype=dtype,
    )
    Y = hl.Func("Y")
    (state_idx,) = gather.args()
    Y[state_idx] = gather[state_idx]

    if not target.has_gpu_feature():
        parallel = hl.Var("parallel")
        outer = hl.Var("outer")
        inner = hl.Var("inner")
        # task_size = 1024
        Y.split(state_idx, outer, inner, vector_size, GuardWithIf)
        # Y.split(outer, parallel, outer, task_size, GuardWithIf)
        # Y.reorder(inner, outer, parallel)
        # Y.vectorize(inner)

        if has_off_diag:
            (r_gather,) = gather.rvars()
            gather.split(state_idx, outer, inner, vector_size, GuardWithIf)
            gather.reorder(inner, outer).vectorize(inner)
            gather.update(0).split(state_idx, outer, inner, vector_size, GuardWithIf)
            gather.update(0).reorder(inner, r_gather, outer).vectorize(inner)
            gather.compute_at(Y, outer).store_at(Y, outer)
            gather.bound_storage(state_idx, vector_size).store_in(hl.MemoryType.Register)

            if off_diag_coeff.has_update_definition():
                off_diag_coeff.split(state_idx, outer, inner, vector_size, RoundUp)
                off_diag_coeff.reorder(inner, outer).vectorize(inner)
                (_, term_idx) = off_diag_coeff.args()
                (r_off_diag,) = off_diag_coeff.rvars()
                off_diag_coeff.update(0).split(state_idx, outer, inner, vector_size, GuardWithIf)
                off_diag_coeff.update(0).reorder(inner, r_off_diag, term_idx, outer).vectorize(inner).unroll(r_off_diag)
                off_diag_coeff.compute_at(gather, r_gather).store_at(gather, r_gather)
                off_diag_coeff.bound_storage(state_idx, vector_size).bound_storage(term_idx, 1)
                off_diag_coeff.store_in(hl.MemoryType.Register)

            if info.is_state_index_identity:
                state_to_index.compute_inline()
            elif not info.has_permutation_symmetries:
                (r_binomial,) = state_to_index.rvars()
                (state_idx, term_idx) = state_to_index.args()
                state_to_index.split(state_idx, outer, inner, vector_size, GuardWithIf)
                state_to_index.split(outer, parallel, outer, task_size, GuardWithIf)
                state_to_index.vectorize(inner)
                state_to_index.update(0).split(state_idx, outer, inner, vector_size, GuardWithIf)
                state_to_index.update(0).split(outer, parallel, outer, task_size, GuardWithIf)
                state_to_index.update(0).reorder(inner, r_binomial, term_idx, outer, parallel).vectorize(inner)
                state_to_index.compute_at(gather, r_gather).store_at(gather, r_gather)
                state_to_index.bound_storage(state_idx, vector_size).bound_storage(term_idx, 1)
                state_to_index.store_in(hl.MemoryType.Register)
        else:
            gather.compute_inline()
            off_diag_coeff.compute_inline()

        if has_diag:
            diag_coeff.split(state_idx, outer, inner, vector_size, RoundUp)
            diag_coeff.reorder(inner, outer).vectorize(inner)
            diag_coeff.update(0).split(state_idx, outer, inner, vector_size, GuardWithIf)
            (r_diag,) = diag_coeff.rvars()
            diag_coeff.update(0).reorder(inner, r_diag, outer).vectorize(inner)
            where = gather if has_off_diag else Y
            diag_coeff.compute_at(where, outer).store_at(where, outer)
            diag_coeff.bound_storage(state_idx, vector_size).store_in(hl.MemoryType.Register)
            # if has_off_diag:
            #     diag_coeff.compute_with(off_diag_coeff, inner)
            #     # diag_coeff.update(0).compute_with(off_diag_coeff.update(0), outer)
        else:
            diag_coeff.compute_inline()
        
        # Y.split(outer, parallel, outer, task_size, GuardWithIf)
        # Y.parallel(parallel)
        Y.vectorize(inner)

    count = init_alpha_buf.dim(0).extent()
    init_alpha_buf.dim(0).set_min(0).set_stride(1)
    init_coeffs_buf.dim(0).set_min(0).set_stride(1).set_extent(count)
    if has_off_diag:
        number_reduced = terms.s_2d.shape[1]
        number_terms = terms.s_2d.shape[0]
        v_re_buf.dim(0).set_min(0).set_stride(1).set_extent(number_reduced)
        v_re_buf.dim(1).set_min(0).set_stride(number_reduced).set_extent(number_terms)
        v_im_buf.dim(0).set_min(0).set_stride(1).set_extent(number_reduced)
        v_im_buf.dim(1).set_min(0).set_stride(number_reduced).set_extent(number_terms)
    v_diag_re_buf.dim(0).set_min(0).set_stride(1)
    v_diag_im_buf.dim(0).set_min(0).set_stride(1)
    basis_states_buf.dim(0).set_min(0).set_stride(1)
    X_buf.dim(0).set_min(0).set_stride(1)
    Y.output_buffer().dim(0).set_min(0).set_stride(1).set_extent(count)

    args = (init_alpha_buf, init_norms_buf, init_coeffs_buf)
    if has_off_diag:
        args += (v_re_buf, v_im_buf)
    if has_diag:
        args += (v_diag_re_buf, v_diag_im_buf)
    args += (basis_states_buf, norms_buf, X_buf,)

    if verbose:
        Y.print_loop_nest()
        target = target.with_feature(hl.TargetFeature.NoBoundsQuery)
        target = target.with_feature(hl.TargetFeature.NoAsserts)
        # .with_feature(hl.TargetFeature.AVX512_Zen4)
        Y.compile_to_lowered_stmt(
            "off_diag_matrix.stmt.html",
            args,
            fmt=hl.StmtOutputFormat.HTML,
            target=target,
        )
    return Y, args, keep_alive


def off_diag_matrix_kernel(*args, **kwargs) -> CompiledKernel:
    builder = lambda: _build_off_diag_matrix_kernel(*args, **kwargs)
    return COMPILER.link_and_load(builder, "off_diag_matrix_kernel")


def create_halide_buffer_view(arr, ffi=COMPILER.ffi):
    dtype = arr.dtype

    dim = ffi.new("halide_dimension_t[]", arr.ndim)
    for i in range(arr.ndim):
        dim[arr.ndim - 1 - i].min = 0
        dim[arr.ndim - 1 - i].extent = arr.shape[i]
        dim[arr.ndim - 1 - i].stride = arr.strides[i] // dtype.itemsize
        dim[arr.ndim - 1 - i].flags = 0

    buf = ffi.new("halide_buffer_t*")
    buf.device = 0
    buf.device_interface = ffi.NULL
    buf.host = ffi.from_buffer("uint8_t*", arr)
    buf.flags = 0
    buf.type.code = {"i": 0, "u": 1, "f": 2}[dtype.kind]
    buf.type.bits = 8 * dtype.itemsize
    buf.type.lanes = 1
    buf.dimensions = arr.ndim
    buf.dim = dim
    buf.padding = ffi.NULL

    return buf, (dim,)
