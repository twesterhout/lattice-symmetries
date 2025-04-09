import cffi, time, numpy as np, os, igraph as ig, halide as hl, tempfile, subprocess, weakref, scipy.sparse.linalg, sympy
import lattice_symmetries as ls
from dataclasses import dataclass
from sympy import S, simplify, Rational
from sympy.combinatorics import Permutation
from loguru import logger

# def _constraints(buf, shape):
#     s = 1
#     for i in range(buf.dimensions()):
#         d = buf.dim(i).set_min(0).set_stride(s)
#         if shape[i] is not None: d.set_extent(shape[i])
#         s *= shape[i] if shape[i] is not None else buf.dim(i).extent()


# from halide import Var, RVar, Func, Buffer, ImageParam, Int, UInt, Float

# class DynamicArgs:
#     def __init__(self, dtype):
#         self.alpha0 = ImageParam(Int(64),  1, "alpha0")
#         self.norms0 = ImageParam(UInt(16), 1, "norms0")
#         self.coeff0 = ImageParam(dtype,       1, "coeff0")
#         self.v_diag_re = ImageParam(Float(64), 1, "v_diag_re")
#         self.v_diag_im = ImageParam(Float(64), 1, "v_diag_im")
#         self.v_2d_re = ImageParam(Float(64), 2, "v_2d_re")
#         self.v_2d_im = ImageParam(Float(64), 2, "v_2d_im")
#         self.basis_states = ImageParam(Int(64), 1, "basis_states")
#         self.norms = ImageParam(UInt(16), 1, "norms")
#         self.X = ImageParam(dtype, 1, "X")
#         n, n_s = self.alpha0.dim(0).extent(), self.X.dim(0).extent()
#         (n_t, n_r), (n_d,) = lowered.s_2d.shape, lowered.s_diag.shape
#         _constraints(self.alpha0, [None])
#         _constraints(self.norms0, [n])
#         _constraints(self.coeff0, [n])
#         _constraints(self.v_diag_re, [n_d])
#         _constraints(self.v_diag_im, [n_d])
#         _constraints(self.v_2d_re, [n_r, n_t])
#         _constraints(self.v_2d_im, [n_r, n_t])
#         _constraints(self.basis_states, [n])
#         _constraints(self.norms, [n])
#         _constraints(self.X, [None])

#     def args(self, lowered):
#         t = (self.alpha0,)
#         if lowered.has_diag:
#             t += (self.v_diag_re,)
#             if lowered.has_imag: t += (self.v_diag_im,)
#         if lowered.has_off_diag:
#             t += (self.v_2d_re,)
#             if lowered.has_imag: t += (self.v_2d_im,)
#         t += (self.X,)
#         return t

# class StaticArgs:
#     def __init__(self, lowered):
#         self.s_2d = hl.Buffer(lowered.s_2d.view(np.int64)) if lowered.has_off_diag else None
#         self.s_diag = hl.Buffer(lowered.s_diag.view(np.int64)) if lowered.has_diag else None
#         self.xors = hl.Buffer(lowered.mask.view(np.int64)) if lowered.has_off_diag else None
# 
# _i64 = lambda x: hl.reinterpret(hl.Int(64), x)
# _u64 = lambda x: hl.reinterpret(hl.UInt(64), x)
# _f64 = lambda x: hl.reinterpret(hl.Float(64), x)
# parity = lambda x: _u64(hl.popcount(hl.u64(x))) << hl.u64(63)

# def build_kernel(lowered, dtype=Float(32)):
#     has_imag = lowered.has_imag
#     ctx, static = DynamicArgs(dtype), StaticArgs(lowered)
#     n, n_s = ctx.alpha0.dim(0).extent(), ctx.X.dim(0).extent()
#     (n_t, n_r), (n_d,) = lowered.s_2d.shape, lowered.s_diag.shape
#     _clamp = lambda x: hl.unsafe_promise_clamped(hl.i32(x), 0, n_s - 1)
#     # Computing off diagonal coefficients
#     bi, ti, re, im = Var("bi"), Var("ti"), 0, 0
#     for r in range(n_r):
#         sign_mask = parity(ctx.alpha0[bi] & static.s_2d[r, ti])
#         re += _f64(_u64(ctx.v_2d_re[r, ti]) ^ sign_mask)
#         if has_imag: im += _f64(_u64(ctx.v_2d_im[r, ti]) ^ sign_mask)
#     ndc = Func("ndc"); ndc[ti, bi] = (re, im) if has_imag else re
#     # Computing diagonal coefficients
#     dc = Func("dc"); r = hl.RDom([hl.Range(0, n_d)], "r_diag")
#     dc[bi] = (hl.f64(0), hl.f64(0)) if has_imag else hl.f64(0)
#     sign_mask = parity(ctx.alpha0[bi] & static.s_diag[r])
#     re = _f64(_u64(ctx.v_diag_re[r]) ^ sign_mask)
#     if has_imag: im = _f64(_u64(ctx.v_diag_im[r]) ^ sign_mask)
#     dc[bi] = (dc[bi][0] + re, dc[bi][1] + im) if has_imag else dc[bi] + re
#     # Computing betas
#     beta = Func("beta"); beta[bi, ti] = ctx.alpha0[bi] ^ static.xors[ti]
#     # Gathering the coefficients
#     gather = Func("gather"); gather[bi] = dc[bi] * hl.f64(ctx.X[_clamp(ctx.alpha0[bi])])
#     r = hl.RDom([hl.Range(0, n_t)], "r_gather")
#     gather[bi] += ndc[r, bi] * hl.f64(ctx.X[_clamp(beta[bi, r])])
#     # Storing the results
#     Y = Func("Y"); Y[bi] = hl.cast(dtype, gather[bi])
#     _constraints(Y.output_buffer(), [n])

#     if True:
#         boo, bo = hl.Var("boo"), hl.Var("bo")
#         bii = hl.Var("bii")
#         n_t, n_v = 32, 8

#         # gather.compute_at(Y, bi)
#         # dc.split(bi, bo, bi, n_t * n_v, hl.TailStrategy.RoundUp) # .vectorize(bi)
#         # dc.split(bi, bi, bii, n_v, hl.TailStrategy.Predicate) # .vectorize(bi)
#         # dc.update(0).split(bi, bo, bi, n_t * n_v, hl.TailStrategy.Predicate)
#         # dc.reorder(bi, bo) # .vectorize(bi)
#         # dc.update(0).reorder(dc.rvars()[0], bi, bo) # .vectorize(bi)
#         Y.split(bi, bo, bi, n_t * n_v, hl.TailStrategy.Predicate)# .gpu_lanes(bi)
#         Y.split(bi, bi, bii, n_v, hl.TailStrategy.Predicate)# .gpu_lanes(bi)
#         Y.split(bo, boo, bo, 2**10, hl.TailStrategy.GuardWithIf)
#         Y.gpu_blocks(boo, bo).gpu_threads(bi).unroll(bii) # .gpu_lanes(bii)

#         dc.split(bi, bi, bii, n_v, hl.TailStrategy.Predicate) # .reorder(bii, *dc.rvars())
#         dc.update(0).split(bi, bi, bii, n_v, hl.TailStrategy.Predicate) # .reorder(bii, *dc.rvars())
#         dc.unroll(bii).update(0).reorder(bii, *dc.rvars()).unroll(bii)
#         dc.compute_at(Y, bi).store_in(hl.MemoryType.GPUShared)

#         gather.split(bi, bi, bii, n_v, hl.TailStrategy.Predicate) # .reorder(bii, *dc.rvars())
#         gather.update(0).split(bi, bi, bii, n_v, hl.TailStrategy.Predicate) # .reorder(bii, *dc.rvars())
#         gather.unroll(bii).update(0).reorder(bii, *gather.rvars()).unroll(bii)
#         # 
#         # gather.split(bi, bo, bi, n_t * n_v, hl.TailStrategy.Predicate) # .vectorize(bi)
#         # # gather.split(bi, bi, bii, n_v, hl.TailStrategy.Predicate) # .vectorize(bi)
#         # gather.update(0).split(bi, bo, bi, n_t * n_v, hl.TailStrategy.Predicate) # .vectorize(bi)
#         # # gather.update(0).split(bi, bi, bii, n_v, hl.TailStrategy.Predicate) # .vectorize(bi)
#         # # gather.split(bi, bo, bi, n_v) # .vectorize(bi)
#         # gather.update(0).reorder(gather.rvars()[0], bi, bo) # .vectorize(bi)
#         gather.compute_at(Y, bi).store_in(hl.MemoryType.GPUShared)
#         # Y.split(bo, boo, bo, 128).gpu(bo, boo).gpu_threads(bo).gpu_blocks(boo) # .gpu_lanes(bi)
#         # gather.compute_at(Y, by)
#     if False:
#         bo = hl.Var("bo")
#         rf = hl.Var("rf")
#         n_v = 16
#         (r,) = dc.rvars()
#         dc.split(bi, bo, bi, n_v, hl.TailStrategy.GuardWithIf).vectorize(bi)
#         dc.update(0).split(bi, bo, bi, n_v, hl.TailStrategy.GuardWithIf).reorder(bi, r, bo).vectorize(bi)
#         dc.compute_at(Y, bo)

#         Y.split(bi, bo, bi, n_v, hl.TailStrategy.GuardWithIf).vectorize(bi)
#         gather.split(bi, bo, bi, n_v, hl.TailStrategy.GuardWithIf).vectorize(bi)
#         (r,) = gather.rvars()
#         gather.update(0).split(bi, bo, bi, n_v, hl.TailStrategy.GuardWithIf).reorder(bi, r, bo).vectorize(bi)
#         # intm = gather.update(0).rfactor([(r, rf)])
#         # intm.vectorize(rf, n_v, hl.TailStrategy.GuardWithIf).vectorize(bi)
#         # intm.update(0).vectorize(rf, n_v, hl.TailStrategy.GuardWithIf).vectorize(bi)
#         gather.compute_at(Y, bo)
#         Y.parallel(bo, 2048)

#     return Y, ctx.args(lowered)
rng = np.random.default_rng(5)
p, k = Permutation(np.roll(np.arange(3), shift=-1)), 0
i = ls.BasisInfo(bits=3, symmetries=[(p, Rational(k, p.order()))])
kernels = ls.compiler.build_kernels()

h = ls.heisenberg(ig.Graph.Ring(3, circular=True), h=0.25)
terms = ls.expression.pauli2nbts(simplify(h.raw))

_states, _norms = ls.compiler.EnumerateStates(i, kernels)()
print(_states)
print(_norms)

bs_ctx = ls.compiler.bs_ctx_t(i)
search_ctx = ls.compiler.search_ctx_t(i, _states, prefix_bits=1)
_, oc_ctx = ls.compiler.oc_ctx_t(terms)


alpha = np.zeros(64, dtype=np.uint64); alpha[:_states.size] = _states
norm = np.zeros(64, dtype=np.uint16); norm[:_norms.size] = _norms
x = rng.random(64, dtype=np.float64)
print(x[:_norms.size])
out = np.zeros(64, dtype=np.float64)
oc_ctx.p.alpha0, oc_ctx.p.norm0 = ls.compiler.cb_u64(alpha), ls.compiler.cb_u16(norm)
oc_ctx.p.X, oc_ctx.p.norm = ls.compiler.cb_f64(x), ls.compiler.cb_u16(norm)

x_ref = np.zeros(2**3, dtype=np.float64)
x_ref[alpha] = x
y_ref = h.to_dense() @ x_ref - np.diag(np.diag(h.to_dense())) @ x_ref
print(y_ref)

kernels.off_diag(0, oc_ctx.p, bs_ctx.p, search_ctx.p, ls.compiler.b_f64(out))
print(out[:_states.size])
exit(0)



xs = np.arange(8, dtype=np.uint64)
rep, idx = np.zeros(xs.size, dtype=np.uint64), np.zeros(xs.size, dtype=np.int64)
kernels.state_info(xs.size, ls.compiler.cb_u64(xs), ctx.p, ls.compiler.b_u64(rep), ls.compiler.b_i64(idx))
assert list(zip(rep.tolist(), idx.tolist())) == [
    (0, 0), # 000
    (1, 0), # 001
    (1, 1), # 010
    (3, 0), # 011
    (1, 2), # 100
    (3, 2), # 101
    (3, 1), # 110
    (7, 0), # 111
]

# 0000
i = ls.BasisInfo(4)
reps = np.array([1, 3, 7, 12, 15], dtype=np.uint64)
search_ctx = ls.compiler.search_ctx_t(i, reps, prefix_bits=5)

needle = np.array([12, 95043285, 2, 1], dtype=np.int64)
out = np.zeros(2, dtype=np.int64)
kernels.state_to_index(out.size, ls.compiler.cb_u64(needle), search_ctx.p, ls.compiler.b_i64(out))
print(out)



k = 4
symmetries = [(ls.Permutation(np.roll(np.arange(k), -1)), Rational(0))]
symmetries = ls.generate_representation(symmetries)
masks, shifts, is1 = lower_symmetries(symmetries)
assert np.all(is1 == -1)
kernels = ls.compiler.build_kernels()
f = ls.compiler.build_enumerate_states()
ffi = ls.COMPILER.ffi

p = ffi.new("bs_ctx_t *")
p.masks = ffi.from_buffer("const uint64_t*", masks, require_writable=False)
p.shifts = ffi.from_buffer("const uint32_t*", shifts, require_writable=False)
p.is1 = ffi.from_buffer("const int64_t*", is1, require_writable=False)
p.n_m = masks.shape[0]
p.n_r = masks.shape[1]

# void* enumerate_states(i64 nc, i64 *os, u64 *xs, void *candidates64, void *norms64, void const* ctx) {
offsets = np.array([0, 2**k], dtype=np.int64)
xs = np.array([0], dtype=np.uint64)
f.enumerate_states(1, ffi.from_buffer("int64_t const*", offsets), ffi.from_buffer("uint64_t const*", xs), f.candidates, kernels.norm, p)
exit(0)

alpha = np.arange(2**k, dtype=np.uint64)
out = np.zeros(2**k, dtype=np.float64)
n = alpha.size
if n < 64: alpha_t = np.pad(alpha, ((0, 64 - n),), mode="edge")
else: alpha_t = alpha
out = np.zeros(max(64, n), dtype=np.uint16)





kernels.norm(
    ffi.from_buffer("const uint64_t*", alpha_t, require_writable=False),
    p,
    ffi.from_buffer("uint16_t*", out, require_writable=True),
)
print(alpha)
print(out[:n])
exit(0)




  


k = 4
h = ls.heisenberg(ig.Graph.Ring(k, circular=True), h=0.25)
terms = ls.expression.pauli2nbts(simplify(h.raw))
kernels = ls.compiler.build_kernels()

alpha = np.arange(2**k, dtype=np.uint64)
x = rng.random(2**k, dtype=np.float64)
out = np.zeros(2**k, dtype=np.float64)

# print(h.to_dense() @ x)

p_diag, p_off_diag, _keep_alive = ls.compiler.get_ctxs(terms)

n = x.size
if n < 64:
    alpha_t = np.pad(alpha, ((0, 64 - n),), mode="edge")
    x_t = np.pad(x, ((0, 64 - n),), mode="constant")
    out_t = np.pad(out, ((0, 64 - n),), mode="constant")
else:
    alpha_t = alpha
    x_t = x
    out_t = out
ffi = ls.COMPILER.ffi
p_diag.alpha0 = ffi.from_buffer("const uint64_t*", alpha_t, require_writable=False)
p_diag.X = ffi.from_buffer("const double*", x_t, require_writable=False)
p_off_diag.alpha0 = ffi.from_buffer("const uint64_t*", alpha_t, require_writable=False)
p_off_diag.X = ffi.from_buffer("const double*", x_t, require_writable=False)
tick = time.perf_counter()
kernels.matvec(max(64, n), p_diag, p_off_diag, ffi.from_buffer("double*", out_t, require_writable=True))
tock = time.perf_counter(); print(tock - tick)
out[:] = out_t[:n]
print(out)


def matvec(x):
    n = x.size
    if n < 64:
        alpha_t = np.pad(alpha, ((0, 64 - n),), mode="edge")
        x_t = np.pad(x, ((0, 64 - n),), mode="constant")
    else:
        alpha_t = alpha
        x_t = x
    out_t = np.zeros(max(n, 64), dtype=np.float64)

    p_diag.alpha0 = ffi.from_buffer("const uint64_t*", alpha_t, require_writable=False)
    p_diag.X = ffi.from_buffer("const double*", x_t, require_writable=False)
    p_off_diag.alpha0 = ffi.from_buffer("const uint64_t*", alpha_t, require_writable=False)
    p_off_diag.X = ffi.from_buffer("const double*", x_t, require_writable=False)
    kernels.matvec(max(64, n), p_diag, p_off_diag, ffi.from_buffer("double*", out_t, require_writable=True))
    return out_t[:n]
    
matrix = scipy.sparse.linalg.LinearOperator(shape=(2**k, 2**k), dtype=np.float64, matvec=matvec)
eigvals, eigvecs = scipy.sparse.linalg.eigsh(matrix, k=1, which="SA")
print(eigvals)

exit(0)



# p = ffi.new("dc_t *")

kernels.off_diag(0, p, ffi.from_buffer("double*", out, require_writable=True))

# for k in range(0, out.size, 2):
#     kernels.diag(k, p, ffi.from_buffer("double*", out[k:], require_writable=True))
print(out)
del kernels
ffi.dlclose(lib)
exit(0)



Y, args = build_kernel(lowered)
# target =  hl.get_jit_target_from_environment() # .with_feature(hl.TargetFeature.NoAsserts).with_feature(hl.TargetFeature.NoBoundsQuery)
target = hl.get_host_target().with_feature(hl.TargetFeature.JIT).with_feature(hl.TargetFeature.CUDA).with_feature(hl.TargetFeature.Debug).with_feature(hl.TargetFeature.CUDACapability61) # .with_feature(hl.TargetFeature.LargeBuffers) # .with_feature(hl.TargetFeature.NoAsserts).with_feature(hl.TargetFeature.NoBoundsQuery)
# print(target)
Y.compile_to_conceptual_stmt("matvec.stmt", args, target=target)
# Y.compile_to_assembly("matvec.asm", args, target=target)
# Y.compile_to_c("matvec.c", args, target=target)
kernel = Y.compile_to_callable(args, target)
x = rng.random(2**k, dtype=np.float32)
# ref = np.diag(np.diag(h.to_dense())) @ x
# ref = h.to_dense() @ x
x = hl.Buffer(x) # x = np.zeros(m.shape[0]); x[0] = 1
# x.copy_to_device(target)
alpha = hl.Buffer(np.arange(2**k, dtype=np.int64))
# alpha.copy_to_device(target)
out = hl.Buffer(np.zeros(2**k, dtype=np.float32))
# out.copy_to_device(target)
v_re_diag = hl.Buffer(lowered.v_re_diag)
# v_re_diag.copy_to_device(target)
v_re_2d = hl.Buffer(lowered.v_re_2d)
# v_re_2d.copy_to_device(target)
tick = time.perf_counter()
kernel(alpha, v_re_diag, v_re_2d, x, out)
out.device_sync()
tock = time.perf_counter()
print(tock - tick)
# print(ref)
# print(np.asarray(out))

rng = np.random.default_rng(5)
m = h.to_dense()
x = rng.random(m.shape[0], dtype=np.float32) # x = np.zeros(m.shape[0]); x[0] = 1
ref = m @ x
alpha = np.arange(m.shape[0], dtype=np.uint64)
# a = nbts.diag(alpha) * x
# i, c = nbts.off_diag(alpha)
# b = (c * x[i]).sum(axis=0)
# np.testing.assert_allclose(ref, a + b)

out = np.zeros(m.shape[0], dtype=np.float32)
kernel(alpha.view(np.int64), lowered.v_re_diag, lowered.v_re_2d, x, out)

print(out)
print(ref)
# print(lowered.diag(alpha) * x)
# i, c = lowered.off_diag(alpha)
# print(i.T.tolist()[1])
# print(c.T.tolist()[1])
