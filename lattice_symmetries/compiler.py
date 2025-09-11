import cffi, numpy as np, os, pathlib, subprocess, sympy, tempfile, time, threadpoolctl, lattice_symmetries as ls
from dataclasses import dataclass, field
from loguru import logger
from sympy import S, Rational
from sympy.combinatorics import Permutation

FOLDER = pathlib.Path(__file__).parent.resolve()

class KernelCompiler:
    temp: str; ffi: any; cc: str
    def __init__(self, temp_dir=None):
        self.temp = temp_dir or tempfile.mkdtemp(prefix="lattice-symmetries-cache")
        logger.trace(f"'{self.temp}' will be used for compiling kernels.")
        _ = np.zeros(10) # dummy
        self.ffi = cffi.FFI()
        with open(FOLDER / "declarations.h", "r") as f: self.ffi.cdef(f.read())
        # C compiler. We want to support specifying CC='zig cc' hence the additional split
        self.cc = (os.getenv("CC", default="cc"),)
        if "zig" in self.cc[0]: self.cc = self.cc[0].split(" ")
        is_clang = "clang" in self._version()
        opt = lambda p, x: x if p else []
        # Optimization flags
        self.flags = ["-O3" if is_clang else "-Ofast", "-ftree-vectorize", "-ffast-math", "-DNDEBUG"] \
                   + opt(not is_clang, ["-fschedule-insns", "-fschedule-insns2"])
        # Architecture; TODO: optionally add -mcpu=... -msimd128 for WASM
        self.flags += ["-march=native", "-mtune=native"]
        # Warnings
        self.flags += ["-Wall", "-Wextra", "-W", "-Wno-comment", "-Wno-unused-parameter", "-Wno-psabi"] \
                    + opt(is_clang, ["-Wno-nan-infinity-disabled"])
        # OpenMP
        is_gomp = "libgomp" in [i["prefix"] for i in threadpoolctl.threadpool_info() if i["user_api"] == "openmp"]
        self.flags += ["-fopenmp=libgomp" if is_clang and is_gomp else "-fopenmp"]
        # Library
        self.flags += ["-fPIC", "-ffreestanding", "-I", os.getenv("LS_SIMDE_PATH", str(FOLDER))]
        M = os.getenv("LS_M")
        if M is not None:
            assert M in ["1", "2", "3"]; logger.trace(f"Kernels will be compiled for M={M}.")
            self.flags += [f"-DM={M}"]
    def _version(self): return subprocess.run([*self.cc, "--version"], check=True, capture_output=True, text=True).stdout
    def compile(self, *srcs):
        _, out = tempfile.mkstemp(suffix=".so", dir=self.temp)
        args = [*self.cc, *self.flags, "-shared", "-o", out, *map(str, srcs)]
        tick = time.perf_counter(); subprocess.run(args, check=True); tock = time.perf_counter()
        logger.trace(f"Compiled in {tock - tick} seconds. Command was '{' '.join(args)}'")
        return self.ffi.dlopen(out, self.ffi.RTLD_NOW | self.ffi.RTLD_LOCAL)

COMPILER = KernelCompiler()

@dataclass(frozen=True)
class K:
    diag64: any; off_diag64: any; norm64: any
    state_to_index: any; state_info: any
    matvec: any; has_float16: any
    enumerate_states: any; copy_finalize: any
    candidates_simple: any; candidates_hamming: any

def build_kernels():
    lib = COMPILER.compile(FOLDER / "matvec.c")
    return K(
        lib.diag64, lib.off_diag64, lib.norm64,
        lib.state_to_index, lib.state_info,
        lib.matvec, lib.has_float16,
        lib.enumerate_states, lib.copy_finalize,
        lib.candidates_simple, lib.candidates_hamming
    )

KERNELS = build_kernels()

 
@dataclass(frozen=True)
class BasisInfo:
    bits: int; hamming: int | None = None; inversion: int | None = None
    symmetries: list[tuple[Permutation, Rational]] = field(default_factory=list)
    def __post_init__(self):
        if not isinstance(self.bits, int): object.__setattr__(self, "bits", int(self.bits))
        if len(self.symmetries) > 0:
            object.__setattr__(self, "symmetries", ls.generate_representation(self.symmetries))
        for p, _ in self.symmetries: assert len(p.array_form) == self.bits
        if self.bits <= 1: object.__setattr__(self, "symmetries", [])
    has_ps = property(lambda self: len(self.symmetries) > 0)
    is_s2i_id = property(lambda self: self.hamming is None and not self.has_ps)
    @property
    def min_and_max_state_estimate(self) -> tuple[int, int]:
        l = 2**self.hamming - 1 if self.hamming else 0 
        if self.hamming is None:
            # If spin inversion is not None, leave the most significant bit as 0
            r = 2**self.bits - 1 if self.inversion is None else 2 ** (self.bits - 1) - 1
        else:
            assert self.hamming <= self.bits
            r = l << (self.bits - self.hamming) if self.inversion is None \
                else l << (self.bits - 1 - self.hamming)
        return l, r
    


NULL = COMPILER.ffi.NULL
def b_f64(arr): return COMPILER.ffi.from_buffer("f64*", arr, require_writable=True)
def b_u16(arr): return COMPILER.ffi.from_buffer("u16*", arr, require_writable=True)
def b_u64(arr): return COMPILER.ffi.from_buffer("u64*", arr, require_writable=True)
def b_i64(arr): return COMPILER.ffi.from_buffer("i64*", arr, require_writable=True)
def cb_f64(arr): return COMPILER.ffi.from_buffer("const f64*", arr, require_writable=False)
def cb_u8(arr): return COMPILER.ffi.from_buffer("const u8*", arr, require_writable=False)
def cb_u16(arr): return COMPILER.ffi.from_buffer("const u16*", arr, require_writable=False)
def cb_u32(arr): return COMPILER.ffi.from_buffer("const u32*", arr, require_writable=False)
def cb_u64(arr): return COMPILER.ffi.from_buffer("const u64*", arr, require_writable=False)
def cb_i32(arr): return COMPILER.ffi.from_buffer("const i32*", arr, require_writable=False)
def cb_i64(arr): return COMPILER.ffi.from_buffer("const i64*", arr, require_writable=False)
def b_void(arr): return COMPILER.ffi.from_buffer("void*", arr, require_writable=True)
def cb_void(arr): return COMPILER.ffi.from_buffer("const void*", arr, require_writable=False)

@dataclass(frozen=True)
class Ctx: p: any = NULL; keep_alive: any = None

@dataclass
class Term:
    n_s0: int; n_s1: int; n_s2: int; n_sX: int
    s1: any; s20: any; s21: any; sX: any
    v_re: any; v_im: any
    def __init__(self, v, s):
        cnt = np.bitwise_count(s); i = np.argsort(cnt, stable=True); s, v = s[i], v[i]
        self.n_s0, self.n_s1, self.n_s2, self.n_sX = \
            np.sum(cnt == 0), np.sum(cnt == 1), np.sum(cnt == 2), np.sum(cnt > 2)
        k = self.n_s0; self.s1 = s[k:k + self.n_s1].copy(); k += self.n_s1
        c = np.unpackbits(s[k:k + self.n_s2].view(np.uint8).reshape(-1, 8, 1),
            axis=-1, bitorder="little").reshape(-1, 64).nonzero()[1].astype(np.uint64)
        assert len(c) == 2 * self.n_s2; self.s20, self.s21 = 1 << c[::2], 1 << c[1::2]
        self.sX = s[k + self.n_s2:].copy()
        self.v_re, self.v_im = map(np.ascontiguousarray, (v.real.copy(), v.imag.copy()))
    def resize(self, n):
        assert self.n_s0 + self.n_s1 + self.n_s2 + self.n_sX <= n
        for arr in (self.s1, self.s20, self.s21, self.sX, self.v_re, self.v_im):
            arr.resize(n, refcheck=False)
def _stack_terms(ts, mask):
    if len(ts) == 0: return Ctx(COMPILER.ffi.new("oc_t *"))
    n_t, stride = len(ts), max((t.n_s0 + t.n_s1 + t.n_s2 + t.n_sX for t in ts), default=0)
    for t in ts: t.resize(stride)
    a = lambda s: [getattr(t, s) for t in ts]
    n_s0 = np.array([t.n_s0 for t in ts], dtype=np.int32)
    n_s1 = np.array([t.n_s1 for t in ts], dtype=np.int32)
    n_s2 = np.array([t.n_s2 for t in ts], dtype=np.int32)
    n_sX = np.array([t.n_sX for t in ts], dtype=np.int32)
    s1 = np.stack([t.s1 for t in ts])
    s20 = np.stack([t.s20 for t in ts])
    s21 = np.stack([t.s21 for t in ts])
    sX = np.stack([t.sX for t in ts])
    v_re, v_im = np.stack(a("v_re")), np.stack(a("v_im"))
    assert mask.dtype == np.uint64 and mask.size == n_t
    for arr in (s1, s20, s21, sX): assert arr.dtype == np.uint64
    for arr in (v_re, v_im): assert arr.dtype == np.float64
    for arr in (s1, s20, s21, sX, v_re, v_im):
        assert arr.shape == (n_t, stride) and arr.flags["C_CONTIGUOUS"]
    p = COMPILER.ffi.new("oc_t *")
    p.v_re, p.v_im = cb_f64(v_re), cb_f64(v_im)
    p.s1, p.s20, p.s21, p.sX, p.mask = cb_u64(s1), cb_u64(s20), cb_u64(s21), cb_u64(sX), cb_u64(mask)
    p.n_s0, p.n_s1, p.n_s2, p.n_sX = cb_i32(n_s0), cb_i32(n_s1), cb_i32(n_s2), cb_i32(n_sX)
    p.n_t, p.stride = n_t, stride
    return Ctx(p, (n_s0, n_s1, n_s2, n_sX, v_re, v_im, s1, s20, s21, sX, mask))
def _oc_ctx_t(terms):
    if len(terms) == 0: return Ctx(COMPILER.ffi.new("oc_t *"))
    # terms are sorted by x
    v = np.asarray([complex(t.v) for t in terms], dtype=np.complex128)
    x = np.asarray([t.x for t in terms], dtype=np.uint64)
    s = np.asarray([t.s for t in terms], dtype=np.uint64)
    xs, ns = np.unique(x, return_counts=True)
    os = np.pad(np.cumsum(ns), ((1, 0),))
    ts = [Term(v[o:o + n], s[o:o + n]) for o, n in zip(os, ns)]
    return _stack_terms(ts, xs)
def oc_ctx_t(terms):
    nd = sum(int(t.x == 0) for t in terms)
    return _oc_ctx_t(terms[:nd]), _oc_ctx_t(terms[nd:])

def _lower_symmetries(symmetries):
    assert len(symmetries) > 0, "no symmetries"
    assert len(symmetries[0][0].array_form) > 1, "need at least 2 bits"
    nets = [ls.perm2benes(p) for p, _ in symmetries]
    shifts = np.asarray(nets[0].shifts, dtype=np.uint32)
    masks = np.vstack([np.asarray(b.masks, dtype=np.uint64) for b in nets])
    # i = ~np.all(masks == 0, axis=0)
    # if not np.any(i): i[0] = True
    # masks, shifts = np.ascontiguousarray(masks[:, i]), shifts[i]
    chis = [sympy.exp(-2 * sympy.pi * sympy.I * r) for _, r in symmetries]
    is1 = np.asarray([1 if chi == S.One else -1 if chi == -S.One else 0
        for chi in chis], dtype=np.int64)
    chi_re = np.asarray([sympy.re(c) for c in chis], dtype=np.float64)
    chi_im = - np.asarray([sympy.im(c) for c in chis], dtype=np.float64) # TODO: check me!!
    return masks, shifts, is1, chi_re, chi_im
def bs_ctx_t(info):
    if not info.has_ps: return Ctx()
    assert info.hamming is None
    symmetries = info.symmetries
    masks, shifts, is1, chi_re, chi_im = _lower_symmetries(symmetries)
    flags = np.zeros((masks.shape[0], 3), dtype=np.uint8)
    flags[:, 0] = info.inversion is not None
    flags[:, 1] = is1 == 1
    flags[:, 2] = ((info.inversion == 1) & (is1 == 1)) | ((info.inversion == -1) & (is1 == -1))

    p = COMPILER.ffi.new("bs_ctx_t *")
    p.masks, p.shifts, p.flags = cb_u64(masks), cb_u32(shifts), cb_u8(flags)
    p.inversion_mask = 2**info.bits - 1
    p.chi_re, p.chi_im = cb_f64(chi_re), cb_f64(chi_im)
    p.n_m, p.n_r = masks.shape
    return Ctx(p, (masks, shifts, flags, chi_re, chi_im))

def _offset_ranges(reps, bits: int, shift: int):
    assert 0 <= bits < 64, "invalid number_bits"
    assert shift < 64, "invalid shift"
    if bits == 0: return np.array([0, len(reps)], dtype=np.int64), len(reps)
    boundaries = (np.arange(1 << bits, dtype=np.uint64) << shift)
    offsets = np.searchsorted(reps, boundaries, side='left')
    offsets = np.append(offsets, len(reps))
    size = np.max(np.diff(offsets))
    # Normalize ranges to have equal size
    offsets[:-1] = np.minimum(offsets[:-1], len(reps) - size)
    return offsets, size
def search_ctx_t(info, reps=None, norms=None, prefix_bits: int = 20):
    if reps is None and norms is None and info.is_s2i_id: return Ctx()
    prefix_bits = max(0, min(info.bits, prefix_bits))
    shift = info.bits - prefix_bits
    offsets, size = _offset_ranges(reps, prefix_bits, shift)
    p = COMPILER.ffi.new("search_ctx_t *")
    p.reps, p.norm, p.offsets = cb_u64(reps), cb_u16(norms), cb_i64(offsets)
    p.range_size, p.shift, p.mask = size, shift, 2**prefix_bits - 1
    n, steps = size, 0
    while n > 1: n -= n // 2; steps += 1
    p.steps = steps
    return Ctx(p, (reps, norms, offsets))

def enumerate_states(info, ctx=None):
    l, r = info.min_and_max_state_estimate
    if info.is_s2i_id:
        states = np.arange(l, r + 1, dtype=np.uint64)
        norms = np.ones(states.size, dtype=np.uint16)
    elif info.has_ps and info.hamming is None:
        if ctx is None: ctx = bs_ctx_t(info)
        chunk_size = max(1024, (r - l + 1) // (128 * os.cpu_count()))
        starts = np.arange(l, r + 1, chunk_size)
        sizes = np.append(np.diff(starts), [r - starts[-1] + 1])
        starts = starts - 1
        total_size = COMPILER.ffi.new("i64 *")
        chunks = KERNELS.enumerate_states(starts.size, cb_i64(sizes), cb_u64(starts),
            KERNELS.candidates_simple, ctx.p, total_size)
        if chunks == NULL: raise MemoryError("enumerate_states kernel failed to allocate memory")
        states, norms = np.empty(total_size[0], dtype=np.uint64), np.empty(total_size[0], dtype=np.uint16)
        KERNELS.copy_finalize(starts.size, chunks, b_u64(states), b_u16(norms))
    else:
        raise NotImplementedError()
    states.flags.writeable, norms.flags.writeable = False, False
    return states, norms
def state_to_index(states, ctx, out=None):
    states = np.asarray(states, dtype=np.uint64, order="C")
    assert ctx.p != NULL, "don't invoke compiler.state_to_index if is_s2i_id==True"
    assert states.ndim ==1, "expected a one-dimensional array"
    if out is None: out = np.empty(states.size, dtype=np.int64)
    else: assert out.ndim == 1 and out.dtype == np.int64 \
            and out.size == states.size and out.flags["C_CONTIGUOUS"]
    if ctx.p.range_size > 0: KERNELS.state_to_index(states.size, cb_u64(states), ctx.p, b_i64(out))
    else: out[...] = -1
    return out
def state_info(states, ctx, rep=None, idx=None):
    states = np.asarray(states, dtype=np.uint64, order="C")
    assert ctx.p != NULL, "don't invoke compiler.state_info if has_ps==False"
    assert states.ndim ==1, "expected a one-dimensional array"
    if rep is None: rep = np.empty(states.size, dtype=np.uint64)
    else: assert rep.ndim == 1 and rep.dtype == np.uint64 \
            and rep.size == states.size and rep.flags["C_CONTIGUOUS"]
    if idx is None: idx = np.empty(states.size, dtype=np.int64)
    else: assert idx.ndim == 1 and idx.dtype == np.uint64 \
            and idx.size == states.size and idx.flags["C_CONTIGUOUS"]
    KERNELS.state_info(states.size, cb_u64(states), ctx.p, b_u64(rep), b_i64(idx))
    return rep, idx

def _pad(alpha, norm, x):
    n = alpha.size; p = ((0, 64 - n),)
    if x.dtype == np.float16:
        # TODO: FIXME: we read float32 internally, so we need an extra element at the end
        pass
    if n < 64: return np.pad(alpha, p, mode="edge"), np.pad(norm, p), np.pad(x, p)
    else: return alpha, norm, x
def _suffix(dtype): return dict(float64="f64", complex128="c128")[dtype.name]
def _tc(dtype): return dict(float64=0, float32=1, complex128=3, complex64=4)[dtype.name]
    
@dataclass(frozen=True)
class Matvec:
    diag_ctx: any; off_diag_ctx: any; bs_ctx: any; search_ctx: any
    def __call__(self, alpha0, norm0, x0, x=None, out=None):
        alpha0 = np.asarray(alpha0, dtype=np.uint64, order="C")
        norm0 = np.asarray(norm0, dtype=np.uint16, order="C")
        x0 = np.asarray(x0, order="C")
        x = x0 if x is None else np.asarray(x, order="C")
        dtype, n0 = x0.dtype, alpha0.size; n = max(n0, 64)
        assert x0.size == n0 and norm0.size == n0
        if self.search_ctx.p != NULL: assert x.size == self.search_ctx.keep_alive[0].size
        if out is None: out = np.zeros(n, dtype=dtype)
        else: assert out.shape == (n,) and out.dtype == dtype and out.flags["C_CONTIGUOUS"]
        alpha0, norm0, x0 = _pad(alpha0, norm0, x0)
        KERNELS.matvec(_tc(dtype), n, cb_u64(alpha0), cb_u16(norm0), cb_void(x0), cb_void(x), b_void(out),
            self.diag_ctx.p, self.off_diag_ctx.p, self.bs_ctx.p, self.search_ctx.p)
        return out[:n0]
    def _diag64(self, alpha, x): # NOTE: for testing only
        alpha, x = map(np.ascontiguousarray, (alpha, x))
        alpha0, _, x0 = _pad(alpha, x, x)
        dtype, n = x.dtype, alpha.size
        out = np.zeros(64, dtype=dtype)
        KERNELS.diag64(_tc(dtype), cb_u64(alpha0), cb_void(x0), b_void(out), self.diag_ctx.p)
        return out[:min(n, 64)]
    def _off_diag64(self, alpha, norm, x):
        alpha, norm, x = map(np.ascontiguousarray, (alpha, norm, x))
        alpha0, norm0, x0 = _pad(alpha, norm, x)
        dtype, n = x.dtype, alpha.size
        out = np.zeros(64, dtype=dtype)
        KERNELS.off_diag64(_tc(dtype), cb_u64(alpha0), cb_u16(norm0), cb_void(x), b_void(out), self.off_diag_ctx.p, self.bs_ctx.p, self.search_ctx.p)
        return out[:min(n, 64)]


# import jax, jax.numpy as jnp, scipy
# from functools import partial
# from jax.experimental import pallas as pl
# 
# full = jnp.arange(1024)
# u1 = full[jax.lax.population_count(full) == 5]
# 
# binom = jnp.array([[scipy.special.comb(n, k, exact=True) for k in range(33)] for n in range(33)])
# binom = binom.transpose((1, 0))
# 
# def reverse_bits(v):
#     v = ((v >> 1) & 0x55555555) | ((v & 0x55555555) << 1)
#     v = ((v >> 2) & 0x33333333) | ((v & 0x33333333) << 2)
#     v = ((v >> 4) & 0x0F0F0F0F) | ((v & 0x0F0F0F0F) << 4)
#     v = ((v >> 8) & 0x00FF00FF) | ((v & 0x00FF00FF) << 8)
#     v = ( v >> 16             ) | ( v               << 16)
#     return v
# 
# @partial(jax.jit, static_argnums=(1,))
# def s2i(x, h):
#     # x = reverse_bits(x)
#     i = jnp.zeros_like(x)
#     for k in range(h):
#         n = jax.lax.clz(x)
#         i += binom[k + 1, n]
#         x &= ~(1 << (31 - n))
#     return i
# 
# display("{:010b}".format(u1[10]), 10, u1[10])
# # print(s2i.trace(u1[jnp.array([146, 10, 39])], 5).lower().as_text())
# xs = reverse_bits(u1[jax.random.choice(jax.random.PRNGKey(5432), len(u1), shape=(10_000,))])
# print(s2i(reverse_bits(u1[jnp.array([146, 10, 39])]), 5))
# 
# 
# _ = s2i(xs, h=5)
# %timeit s2i(xs, h=5).block_until_ready()
