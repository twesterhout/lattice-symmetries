import dataclasses
import functools
import itertools
import os
import tempfile
import subprocess
import time
from typing import Callable

import cffi
import halide as hl
import numpy as np
from loguru import logger
from scipy.special import comb

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

  int fixed_hamming_state_to_index(struct halide_buffer_t *, struct halide_buffer_t *);
"""


@dataclasses.dataclass(frozen=True)
class CompiledKernel:
    ffi_fun_ptr: any
    callable: any


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
    return hl.Buffer(np.asarray(binomial_coeffs, dtype=np.int64).reshape(size, size))


def _build_fixed_hamming_state_to_index_kernel(
    number_sites: int, hamming_weight: int, vector_size: int = 8, unroll: bool = True
):
    x = hl.ImageParam(hl.Int(64), 1, "x")
    out = hl.Func("out")
    i_inner = hl.Var("i_inner")
    i_outer = hl.Var("i_outer")

    if hamming_weight == 0 or hamming_weight == number_sites:
        batch_idx = hl.Var("batch_idx")
        out[batch_idx] = hl.cast(hl.Int(64), 0)

        out.split(batch_idx, i_outer, i_inner, vector_size, hl.TailStrategy.GuardWithIf)
        out.vectorize(i_inner)
    else:
        binomials = binomials_buffer()
        batch_idx = hl.Var("batch_idx")
        k = hl.RDom([hl.Range(0, hamming_weight)], "k")

        temp = hl.Func("temp")

        temp[batch_idx] = (hl.cast(hl.Int(64), 0), hl.cast(hl.UInt(64), x[batch_idx]))
        t = temp[batch_idx]
        index = t[0]
        state = t[1]
        n = hl.cast(hl.Int(32), hl.count_trailing_zeros(state))
        temp[batch_idx] = (index + binomials[k + 1, n], state & (state - 1))
        out[batch_idx] = temp[batch_idx][0]

        out.split(batch_idx, i_outer, i_inner, vector_size, hl.TailStrategy.GuardWithIf)
        out.vectorize(i_inner)

        temp.compute_at(out, i_inner).store_in(hl.MemoryType.Register)
        if unroll:
            temp.update(0).unroll(k)

    x.dim(0).set_min(0).set_stride(1)
    out.output_buffer().dim(0).set_min(0).set_stride(1).set_extent(x.dim(0).extent())

    # NOTE: it's important to keep binomials alive until we construct a callable from out.
    return out, [x], dict(binomials=binomials)


def fixed_hamming_state_to_index_kernel(number_sites: int, hamming_weight: int) -> CompiledKernel:
    builder = lambda: _build_fixed_hamming_state_to_index_kernel(number_sites, hamming_weight)
    return COMPILER.link_and_load(builder, "fixed_hamming_state_to_index")


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
