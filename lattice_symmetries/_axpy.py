import halide as hl
import numpy as np
from numpy.typing import NDArray
from loguru import logger
import time


def build_axpy_kernel(dtype=hl.Float(64)):
    alpha_re = hl.Param(dtype, "alpha_re")
    alpha_im = hl.Param(dtype, "alpha_im")
    x = hl.ImageParam(dtype, 2, "x")
    y = hl.ImageParam(dtype, 2, "y")

    out = hl.Func("out")
    c = hl.Var("c")
    i = hl.Var("i")

    re = hl.Func("re")
    im = hl.Func("im")
    re[i] = alpha_re * x[0, i] - alpha_im * x[1, i] + y[0, i]
    im[i] = alpha_re * x[1, i] + alpha_im * x[0, i] + y[1, i]

    out[c, i] = hl.select(c == 0, re[i], im[i])
    out.bound(c, 0, 2)

    inner = hl.Var("inner")
    outer = hl.Var("outer")

    block_size = 8
    strategy = hl.TailStrategy.GuardWithIf

    out.split(i, outer, inner, block_size, strategy)
    out.reorder(inner, c, outer)
    out.vectorize(inner)
    out.unroll(c)

    out.specialize(alpha_re == 0)
    out.specialize(alpha_im == 0)

    x.dim(0).set_min(0).set_extent(2).set_stride(1)
    x.dim(1).set_min(0).set_stride(2)
    y.dim(0).set_min(0).set_extent(2).set_stride(1)
    y.dim(1).set_min(0).set_stride(2)
    out.output_buffer().dim(0).set_min(0).set_extent(2).set_stride(1)
    out.output_buffer().dim(1).set_min(0).set_extent(x.dim(1).extent()).set_stride(2)

    return out, [alpha_re, alpha_im, x, y]


_AXPY_KERNEL = None
logger.debug("Initializing ...")


def axpy(alpha: complex, x: NDArray[np.complex128], y: NDArray[np.complex128]):
    """
    A low-level primitive to compute `y += alpha * x` inplace.

    Normally, one would reach out to BLAS operations exposed in SciPy, but
    SciPy is typically compiled with 32-bit BLAS. This means that functions in
    SciPy don't work for arrays larger than 2^31.
    """
    if x.dtype != y.dtype or x.dtype != np.dtype("complex128"):
        raise ValueError("axpy currently only supports complex128")
    if x.shape != y.shape:
        raise ValueError(f"shape mismatch: {x.shape} != {y.shape}")
    if not (
        x.flags.c_contiguous
        and y.flags.c_contiguous
        or x.flags.f_contiguous
        and y.flags.f_contiguous
    ):
        raise ValueError("x and y must be contiguous arrays")

    global _AXPY_KERNEL
    if _AXPY_KERNEL is None:
        logger.debug("Compiling axpy kernel ...")
        tick = time.perf_counter()
        func, params = build_axpy_kernel()
        _AXPY_KERNEL = func.compile_to_callable(params, target=hl.get_jit_target_from_environment())
        tock = time.perf_counter()
        logger.debug(f"Successfully compiled axpy kernel in {tock - tick} seconds.")

    size = x.size
    alpha = complex(alpha)

    offset = 0
    max_allowed_size = np.iinfo(np.dtype("int32")).max
    while size > 0:
        n = min(size, max_allowed_size)
        print(x[offset : offset + n], y[offset : offset + n])
        _AXPY_KERNEL(
            alpha.real,
            alpha.imag,
            x[offset : offset + n].view(np.float64).reshape(-1, 2),
            y[offset : offset + n].view(np.float64).reshape(-1, 2),
            y[offset : offset + n].view(np.float64).reshape(-1, 2),
        )
        offset += n
        size -= n
