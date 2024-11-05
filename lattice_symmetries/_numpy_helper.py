import numpy as np
from lattice_symmetries._ls import lib, ffi


@ffi.def_extern()
def ls_alloc_numpy_array_1d(size, out):
    arr = np.empty((size,), dtype=np.uint8)
    handle = ffi.new_handle(arr)

    out.data = ffi.from_buffer("uint8_t*", arr, require_writable=True)
    out.handle = handle

    # Make sure handle stays alive when the variable goes out of scope
    lib.ls_PyObject_incref(out.handle)
