import dataclasses
import functools
import itertools
import os
import tempfile
import subprocess
import time
from typing import Callable
from lattice_symmetries import _ls

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

  int diag_matrix_complex_kernel(struct halide_buffer_t *, struct halide_buffer_t *, struct halide_buffer_t *);
  int diag_matrix_real_kernel(struct halide_buffer_t *, struct halide_buffer_t *);

  int off_diag_matrix_kernel(struct halide_buffer_t *, struct halide_buffer_t *, struct halide_buffer_t *, struct halide_buffer_t *, struct halide_buffer_t *);
"""


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


@dataclasses.dataclass(frozen=True)
class PauliNonbranchingTerm:
    v: sympy.Expr
    x: int
    s: int

    def act_on_ket(self, ket: int) -> tuple[sympy.Expr, int]:
        sign = 1 - 2 * ((ket & self.s).bit_count() % 2)
        coeff = sign * self.v
        beta = ket ^ self.x
        return coeff, beta

    def act_on_bra(self, bra: int) -> tuple[sympy.Expr, int]:
        coeff, beta = self.act_on_ket(bra)
        return coeff.conjugate(), beta

    @property
    def is_diagonal(self) -> bool:
        return self.x == 0

    def __mul__(self, other):
        if isinstance(other, PauliNonbranchingTerm):
            x = self.x ^ other.x
            s = self.s ^ other.s
            v = self.v * other.v
            return PauliNonbranchingTerm(v=v, x=x, s=s)
        return NotImplemented


@dataclasses.dataclass(init=False)
class PauliLoweredTerms:
    has_diag: bool
    has_off_diag: bool
    s_diag: NDArray[np.uint64]
    v_re_diag: NDArray[np.float64]
    v_im_diag: NDArray[np.float64]
    s_2d: NDArray[np.uint64]
    v_re_2d: NDArray[np.float64]
    v_im_2d: NDArray[np.float64]
    mask: NDArray[np.uint64]

    def __init__(self, terms: list[PauliNonbranchingTerm]):
        terms = sorted(terms, key=lambda x: (x.x, x.s))
        v = np.asarray([complex(t.v) for t in terms], dtype=np.complex128)
        x = np.asarray([t.x for t in terms], dtype=np.uint64)
        s = np.asarray([t.s for t in terms], dtype=np.uint64)

        number_diag = np.sum(x == 0)
        self.has_diag = number_diag != 0
        self.s_diag = s[:number_diag]
        self.v_re_diag = np.ascontiguousarray(v[:number_diag].real)
        self.v_im_diag = np.ascontiguousarray(v[:number_diag].imag)

        if number_diag < len(terms):
            self.has_off_diag = True
            unique_xs, counts = np.unique(x[number_diag:], return_counts=True)
            number_terms = counts.size
            number_reduced = np.max(counts)
            self.s_2d = np.zeros((number_terms, number_reduced), dtype=np.uint64)
            self.v_re_2d = np.zeros((number_terms, number_reduced), dtype=np.float64)
            self.v_im_2d = np.zeros((number_terms, number_reduced), dtype=np.float64)
            offsets = number_diag + np.pad(np.cumsum(counts), ((1, 0),))
            for i in range(number_terms):
                self.s_2d[i, : counts[i]] = s[offsets[i] : offsets[i] + counts[i]]
                self.v_re_2d[i, : counts[i]] = v[offsets[i] : offsets[i] + counts[i]].real
                self.v_im_2d[i, : counts[i]] = v[offsets[i] : offsets[i] + counts[i]].imag
            self.mask = unique_xs
        else:
            self.has_off_diag = False
            self.s_2d = np.zeros((0, 0), dtype=np.uint64)
            self.v_re_2d = np.zeros((0, 0), dtype=np.float64)
            self.v_im_2d = np.zeros((0, 0), dtype=np.float64)
            self.mask = np.zeros(0, dtype=np.uint64)


@dataclasses.dataclass(init=False)
class LoweredOperator:
    terms: PauliLoweredTerms
    kernel: CompiledKernel
    data: _ls.ffi.CData

    def __init__(self, terms: list[PauliNonbranchingTerm]):
        self.terms = PauliLoweredTerms(terms)
        # print(self.terms)
        # self.diag_kernel = diag_matrix_kernel(self.terms) if self.terms.s_diag.size > 0 else None
        self.kernel = off_diag_matrix_kernel(
            self.terms
        )  # if self.terms.s_2d.shape[0] > 0 else None

        # diag_data = _ls.ffi.new("ls_diag_terms *")
        # if self.diag_kernel is None:
        #     diag_data.kernel = _ls.ffi.NULL
        # else:
        #     diag_data.kernel = self.diag_kernel.ffi_fun_ptr
        # diag_data.number_terms = self.terms.s_diag.size
        # diag_data.v_re = _ls.ffi.from_buffer("const double *", self.terms.v_re_diag)
        # diag_data.v_im = _ls.ffi.from_buffer("const double *", self.terms.v_im_diag)
        # self.diag_data = diag_data

        off_diag_data = _ls.ffi.new("ls_off_diag_terms *")
        if self.kernel is None:
            off_diag_data.kernel = _ls.ffi.NULL
        else:
            off_diag_data.kernel = self.kernel.ffi_fun_ptr
        off_diag_data.number_terms = self.terms.s_2d.shape[0]
        off_diag_data.number_reduced = self.terms.s_2d.shape[1]
        off_diag_data.v_re = _ls.ffi.from_buffer("const double *", self.terms.v_re_2d)
        off_diag_data.v_im = _ls.ffi.from_buffer("const double *", self.terms.v_im_2d)
        off_diag_data.x = _ls.ffi.from_buffer("const uint64_t *", self.terms.mask)
        self.data = off_diag_data

    def apply(self, states: NDArray, x: NDArray, out: NDArray):
        has_diag = len(self.terms.s_diag) > 0
        has_off_diag = self.terms.s_2d.shape[0] > 0

        if not has_diag and not has_off_diag:
            out[:] = 0
            return

        dtype = np.dtype("float32")
        states = np.asarray(states, dtype=np.uint64, order="C")
        x = np.asarray(x, dtype=dtype, order="F")
        out = np.asarray(out, dtype=dtype, order="F")

        args = [states.view(np.int64)]
        if has_off_diag:
            args += [self.terms.v_re_2d, self.terms.v_im_2d]
        if has_diag:
            args += [self.terms.v_re_diag, self.terms.v_im_diag]
        args += [x, out]
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
        logger.debug(f"'{self.temp_dir}' will be used for compiling kernels.")
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
        logger.debug(f"Prepared {func_name} in {tock - tick} seconds.")

        tick = time.perf_counter()
        callable = func.compile_to_callable(params, target=self.target)
        tock = time.perf_counter()
        logger.debug(f"Compiled {func_name} to a Halide::Callable in {tock - tick} seconds.")

        tick = time.perf_counter()
        _, object_file = tempfile.mkstemp(suffix=".o", dir=self.temp_dir)
        func.compile_to_object(object_file, params, fn_name=func_name, target=self.target)
        tock = time.perf_counter()
        logger.debug(f"Compiled {func_name} to {object_file} in {tock - tick} seconds.")

        tick = time.perf_counter()
        library_file = os.path.splitext(object_file)[0] + ".so"
        self._link(library_file, object_file)
        tock = time.perf_counter()
        logger.debug(f"Linked {object_file} to {library_file} in {tock - tick} seconds.")

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


def _build_is_representative_kernel(info: BasisInfo):
    if info.number_bits > 64:
        raise NotImplementedError("is_representative not yet implemented for large systems")

    x = hl.ImageParam(hl.Int(64), 1, "x")
    norm = hl.Func("norm")
    i = hl.Var("i_inner")
    keep_alive = dict()

    if info.hamming_weight is None and info.spin_inversion is None and len(info.symmetries) == 0:
        # Identity function
        norm[i] = hl.cast(hl.UInt(16), 1)
    else:
        raise NotImplementedError("😭")

    x.dim(0).set_min(0).set_stride(1)
    norm.output_buffer().dim(0).set_min(0).set_stride(1).set_extent(x.dim(0).extent())

    return norm, [x], keep_alive


def is_representative_kernel(info: BasisInfo) -> CompiledKernel:
    builder = lambda: _build_is_representative_kernel(info)
    return COMPILER.link_and_load(builder, "is_representative_kernel")


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


def _build_off_diag_coeff(terms, *, alpha, v_re, v_im, s, real_only: bool):
    # Compute off-diagonal matrix elements
    off_diag_coeff = hl.Func("off_diag_coeff")
    state_idx = hl.Var("state_idx")
    term_idx = hl.Var("term_idx")
    if real_only:
        off_diag_coeff[state_idx, term_idx] = hl.cast(hl.Float(64), 0)
    else:
        off_diag_coeff[state_idx, term_idx] = (hl.cast(hl.Float(64), 0), hl.cast(hl.Float(64), 0))

    if terms.has_off_diag:
        number_reduced = terms.s_2d.shape[1]
        r = hl.RDom([hl.Range(0, number_reduced)], "r_off_diag")
        sign = hl.popcount(alpha[state_idx] & s[r, term_idx]) << 63
        re = hl.reinterpret(hl.Float(64), hl.reinterpret(hl.Int(64), v_re[r, term_idx]) ^ sign)
        im = hl.reinterpret(hl.Float(64), hl.reinterpret(hl.Int(64), v_im[r, term_idx]) ^ sign)
        if real_only:
            re += off_diag_coeff[state_idx, term_idx]
            off_diag_coeff[state_idx, term_idx] = re
        else:
            re += off_diag_coeff[state_idx, term_idx][0]
            im += off_diag_coeff[state_idx, term_idx][1]
            off_diag_coeff[state_idx, term_idx] = (re, im)
    return off_diag_coeff


def _build_diag_coeff(terms, *, alpha, v_re, v_im, s, real_only: bool):
    diag_coeff = hl.Func("diag_coeff")
    state_idx = hl.Var("state_idx")
    if real_only:
        diag_coeff[state_idx] = hl.cast(hl.Float(64), 0)
    else:
        diag_coeff[state_idx] = (hl.cast(hl.Float(64), 0), hl.cast(hl.Float(64), 0))

    if terms.has_diag:
        r = hl.RDom([hl.Range(0, terms.s_diag.size)], "r_diag")
        sign = hl.popcount(alpha[state_idx] & s[r]) << 63
        re = hl.reinterpret(hl.Float(64), hl.reinterpret(hl.Int(64), v_re[r]) ^ sign)
        im = hl.reinterpret(hl.Float(64), hl.reinterpret(hl.Int(64), v_im[r]) ^ sign)
        if real_only:
            diag_coeff[state_idx] = re + diag_coeff[state_idx]
        else:
            diag_coeff[state_idx] = (re + diag_coeff[state_idx][0], im + diag_coeff[state_idx][1])
    return diag_coeff


def _build_off_diag_matrix_kernel(
    terms: PauliLoweredTerms,
    *,
    spin_inversion=None,
    spin_inversion_mask=None,
    vector_size=8,
    target=hl.get_jit_target_from_environment(),
    dtype=hl.Float(32),
    real_only=True,
    verbose=False,
):
    GuardWithIf = hl.TailStrategy.GuardWithIf
    RoundUp = hl.TailStrategy.RoundUp
    has_diag = terms.has_diag
    has_off_diag = terms.has_off_diag

    alpha_buf = hl.ImageParam(hl.Int(64), 1, "alpha")
    v_re_buf = hl.ImageParam(hl.Float(64), 2, "v_re")
    v_im_buf = hl.ImageParam(hl.Float(64), 2, "v_im")
    v_diag_re_buf = hl.ImageParam(hl.Float(64), 1, "v_diag_re")
    v_diag_im_buf = hl.ImageParam(hl.Float(64), 1, "v_diag_im")
    X_buf = hl.ImageParam(dtype, 1, "X")

    s_buf = hl.Buffer(terms.s_2d.view(np.int64), name="s") if has_off_diag else None
    mask_buf = hl.Buffer(terms.mask.view(np.int64), name="mask") if has_off_diag else None
    s_diag_buf = hl.Buffer(terms.s_diag.view(np.int64), name="s_diag") if has_diag else None

    # Compute off-diagonal matrix elements
    off_diag_coeff = _build_off_diag_coeff(
        terms,
        alpha=alpha_buf,
        v_re=v_re_buf,
        v_im=v_im_buf,
        s=s_buf,
        real_only=real_only,
    )
    # coeff_temp = hl.Func("coeff_temp")
    # term_idx = hl.Var("term_idx")
    # if real_only:
    #     coeff_temp[state_idx, term_idx] = hl.cast(hl.Float(64), 0)
    # else:
    #     coeff_temp[state_idx, term_idx] = (hl.cast(hl.Float(64), 0), hl.cast(hl.Float(64), 0))
    # if has_off_diag:
    #     r1 = hl.RDom([hl.Range(0, number_reduced)], "r1")
    #     sign = hl.popcount(alpha_buf[state_idx] & s_buf[r1, term_idx]) << 63
    #     re = hl.reinterpret(hl.Float(64), hl.reinterpret(hl.Int(64), v_re_buf[r1, term_idx]) ^ sign)
    #     im = hl.reinterpret(hl.Float(64), hl.reinterpret(hl.Int(64), v_im_buf[r1, term_idx]) ^ sign)
    #     if real_only:
    #         re += coeff_temp[state_idx, term_idx]
    #         coeff_temp[state_idx, term_idx] = re
    #     else:
    #         re += coeff_temp[state_idx, term_idx][0]
    #         im += coeff_temp[state_idx, term_idx][1]
    #         coeff_temp[state_idx, term_idx] = (re, im)

    # Compute diagonal matrix elements
    diag_coeff = _build_diag_coeff(
        terms,
        alpha=alpha_buf,
        v_re=v_diag_re_buf,
        v_im=v_diag_im_buf,
        s=s_diag_buf,
        real_only=real_only,
    )

    temp = hl.Func("temp")
    (state_idx,) = diag_coeff.args()
    print(state_idx)
    number_states = X_buf.dim(0).extent()
    if has_diag:
        index = hl.cast(hl.Int(32), alpha_buf[state_idx])
        index = hl.unsafe_promise_clamped(index, 0, number_states - 1)
        if real_only:
            init_coeff = diag_coeff[state_idx] * X_buf[index]
        else:
            init_coeff = diag_coeff[state_idx][0] * X_buf[index]
    else:
        init_coeff = 0
    temp[state_idx] = hl.cast(dtype, init_coeff)

    # if has_diag:
    #     r_diag = hl.RDom([hl.Range(0, terms.s_diag.size)], "r_diag")
    #     diag_coeff = hl.Func("diag_coeff")
    #     if real_only:
    #         diag_coeff[state_idx] = hl.cast(hl.Float(64), 0)
    #     else:
    #         diag_coeff[state_idx] = (hl.cast(hl.Float(64), 0), hl.cast(hl.Float(64), 0))
    #     sign = hl.popcount(alpha_buf[state_idx] & s_diag_buf[r_diag]) << 63
    #     re = hl.reinterpret(hl.Float(64), hl.reinterpret(hl.Int(64), v_diag_re_buf[r_diag]) ^ sign)
    #     im = hl.reinterpret(hl.Float(64), hl.reinterpret(hl.Int(64), v_diag_im_buf[r_diag]) ^ sign)
    #     if real_only:
    #         diag_coeff[state_idx] = re + diag_coeff[state_idx]
    #     else:
    #         diag_coeff[state_idx] = (re + diag_coeff[state_idx][0], im + diag_coeff[state_idx][1])
    #     index = hl.cast(hl.Int(32), alpha_buf[state_idx])
    #     index = hl.unsafe_promise_clamped(index, 0, number_states - 1)
    #     if real_only:
    #         init_coeff = diag_coeff[state_idx] * X_buf[index]
    #     else:
    #         init_coeff = diag_coeff[state_idx][0] * X_buf[index]
    # else:
    #     init_coeff = 0

    if has_off_diag:
        number_terms = len(terms.mask)
        r_gather = hl.RDom([hl.Range(0, number_terms)], "r_gather")
        if real_only:
            coeff = hl.cast(dtype, off_diag_coeff[state_idx, r_gather])
        else:
            coeff = hl.cast(dtype, off_diag_coeff[state_idx, r_gather][0])
        beta = alpha_buf[state_idx] ^ mask_buf[r_gather]
        # if spin_inversion is not None:
        #     beta = hl.min(beta, beta ^ spin_inversion_mask)
        #     if spin_inversion == -1:
        #         flip = (beta ^ spin_inversion_mask) < beta
        #         coeff = hl.select(flip, -coeff, coeff)
        index = hl.cast(hl.Int(32), beta)
        index = hl.unsafe_promise_clamped(index, 0, number_states - 1)
        temp[state_idx] = temp[state_idx] + coeff * X_buf[index]

    Y = hl.Func("Y")
    Y[state_idx] = temp[state_idx]

    if not target.has_gpu_feature():
        outer = hl.Var("outer")
        inner = hl.Var("inner")
        Y.split(state_idx, outer, inner, vector_size, GuardWithIf)
        Y.vectorize(inner).parallel(outer)

        if has_off_diag:
            temp.split(state_idx, outer, inner, vector_size, GuardWithIf)
            temp.vectorize(inner)
            temp.update(0).split(state_idx, outer, inner, vector_size, GuardWithIf)
            temp.update(0).reorder(inner, r_gather, outer).vectorize(inner)
            temp.compute_at(Y, outer).store_at(Y, outer)
            temp.bound_storage(state_idx, vector_size).store_in(hl.MemoryType.Register)
            off_diag_coeff.split(state_idx, outer, inner, vector_size, RoundUp)
            off_diag_coeff.vectorize(inner)
            off_diag_coeff.update(0).split(state_idx, outer, inner, vector_size, GuardWithIf)
            (_, term_idx) = off_diag_coeff.args()
            (r_off_diag,) = off_diag_coeff.rvars()
            off_diag_coeff.update(0).reorder(inner, r_off_diag, term_idx, outer).vectorize(inner)
            off_diag_coeff.compute_at(temp, r_gather).store_at(temp, r_gather)
            off_diag_coeff.bound_storage(state_idx, vector_size).bound_storage(term_idx, 1)
            off_diag_coeff.store_in(hl.MemoryType.Register)

        if has_diag:
            diag_coeff.split(state_idx, outer, inner, vector_size, RoundUp)
            diag_coeff.vectorize(inner)
            diag_coeff.update(0).split(state_idx, outer, inner, vector_size, GuardWithIf)
            (r_diag,) = diag_coeff.rvars()
            diag_coeff.update(0).reorder(inner, r_diag, outer).vectorize(inner)
            where = temp if has_off_diag else Y
            diag_coeff.compute_at(where, outer).store_at(where, outer)
            diag_coeff.bound_storage(state_idx, vector_size).store_in(hl.MemoryType.Register)

    count = alpha_buf.dim(0).extent()
    alpha_buf.dim(0).set_min(0).set_stride(1)
    if has_off_diag:
        number_reduced = terms.s_2d.shape[1]
        v_re_buf.dim(0).set_min(0).set_stride(1).set_extent(number_reduced)
        v_re_buf.dim(1).set_min(0).set_stride(number_reduced).set_extent(number_terms)
        v_im_buf.dim(0).set_min(0).set_stride(1).set_extent(number_reduced)
        v_im_buf.dim(1).set_min(0).set_stride(number_reduced).set_extent(number_terms)
    if has_diag:
        v_diag_re_buf.dim(0).set_min(0).set_stride(1)
        v_diag_im_buf.dim(0).set_min(0).set_stride(1)
    X_buf.dim(0).set_min(0).set_stride(1)
    Y.output_buffer().dim(0).set_min(0).set_stride(1).set_extent(count)

    args = (alpha_buf,)
    if has_off_diag:
        args += (v_re_buf, v_im_buf)
    if has_diag:
        args += (v_diag_re_buf, v_diag_im_buf)
    args += (X_buf,)

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
    return Y, args, dict()


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
