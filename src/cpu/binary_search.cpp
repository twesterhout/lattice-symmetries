#include <lattice_symmetries/lattice_symmetries.h>
#include <algorithm>

extern "C" LATTICE_SYMMETRIES_EXPORT uint64_t ls_binary_search(uint64_t const* const data,
                                                               uint64_t const        size,
                                                               uint64_t const        key)
{
    auto const* p = std::lower_bound(data, data + size, key);
    if (p == data + size || *p != key) { return size; }
    return static_cast<uint64_t>(p - data);
}
