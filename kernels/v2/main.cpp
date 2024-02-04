#include "lattice_symmetries_types.h"
#include <Halide.h>

using namespace Halide;

int main() {

    Runtime::Buffer<uint64_t, 2> binomials(16, 16);

    try {
        auto f = mk_fixed_hamming_state_to_index_kernel(16, 8, binomials.raw_buffer());
    } catch (CompileError &e) {
        fprintf(stderr, "%s\n", e.what());
        return 1;
    } catch (RuntimeError &e) {
        fprintf(stderr, "%s\n", e.what());
        return 1;
    }

    return 0;
}
