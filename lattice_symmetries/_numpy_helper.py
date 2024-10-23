import numpy as np
from loguru import logger
from lattice_symmetries._ls import lib, ffi


ls_alloc_numpy_array_1d_handle = None


@ffi.def_extern()
def ls_alloc_numpy_array_1d(size, out):
    global ls_alloc_numpy_array_1d_handle
    if ls_alloc_numpy_array_1d_handle is not None:
        msg = "'ls_alloc_numpy_array_1d_handle' has not been cleaned up properly. This is likely a bug. The code will crash ..."
        logger.error(msg)

    arr = np.empty((size,), dtype=np.uint8)
    handle = ffi.new_handle(arr)
    out.data = ffi.from_buffer("uint8_t*", arr, require_writable=True)
    out.handle = handle

    # Save the handle such that it's kept alive. This is __very__ important
    ls_alloc_numpy_array_1d_handle = handle
