#include "network.hpp"
#include "bits.hpp"
#include "cpu/benes_forward_512.hpp"
#include "cpu/benes_forward_64.hpp"
#include <algorithm>
#include <numeric>

namespace lattice_symmetries {

small_network_t::small_network_t(fat_benes_network_t const& fat) noexcept
    : masks{}
    , deltas{}
    , depth{static_cast<uint16_t>(fat.masks.size())}
    , width{static_cast<uint16_t>(fat.size)}
{
    LATTICE_SYMMETRIES_CHECK(fat.size <= 64U, "permutation too long, use big_network_t instead");
    LATTICE_SYMMETRIES_CHECK(fat.masks.size() == fat.deltas.size(), "invalid fat_benes_network_t");
    LATTICE_SYMMETRIES_CHECK(fat.masks.size() <= max_depth, "fat_benes_network_t is too deep");
    std::transform(std::begin(fat.masks), std::end(fat.masks), std::begin(masks),
                   [](auto const& m) { return m.words[0]; });
    std::copy(std::begin(fat.deltas), std::end(fat.deltas), std::begin(deltas));
    std::fill(std::next(std::begin(masks), depth), std::end(masks), uint64_t{0});
    std::fill(std::next(std::begin(deltas), depth), std::end(deltas), uint64_t{0});
}

small_network_t::small_network_t(uint16_t _depth, uint16_t _width) noexcept
    : masks{}, deltas{}, depth{_depth}, width{_width}
{
    LATTICE_SYMMETRIES_CHECK(_depth <= max_depth, "network too deep");
    LATTICE_SYMMETRIES_CHECK(_width <= 64, "network too wide");
    auto i = 0U;
    for (auto d = 1U; d < _width / 2; ++i, d *= 2) {
        LATTICE_SYMMETRIES_ASSERT(i < max_depth, "");
        deltas[i] = d;
        masks[i]  = 0U;
    }
    for (auto j = i - 1; j-- > 0; ++i) {
        LATTICE_SYMMETRIES_ASSERT(i < max_depth, "");
        deltas[i] = deltas[j];
        masks[i]  = 0;
    }
    LATTICE_SYMMETRIES_CHECK(i == depth, "depth and width are inconsistent");
    for (; i < max_depth; ++i) {
        deltas[i] = 0;
        masks[i]  = 0;
    }
}

// cppcheck-suppress unusedFunction
auto small_network_t::make_fake(uint16_t depth, uint16_t width) noexcept -> small_network_t
{
    return small_network_t{depth, width};
}

small_network_t::operator fat_benes_network_t() const
{
    std::vector<ls_bits512> new_masks(depth);
    std::vector<unsigned>   new_deltas(depth);
    for (auto i = 0U; i < depth; ++i) {
        new_masks[i]  = widen(masks[i]);
        new_deltas[i] = deltas[i];
    }
    return fat_benes_network_t{std::move(new_masks), std::move(new_deltas), width};
}

big_network_t::big_network_t(fat_benes_network_t const& fat) noexcept
    : masks{}
    , deltas{}
    , depth{static_cast<uint16_t>(fat.masks.size())}
    , width{static_cast<uint16_t>(fat.size)}
{
    LATTICE_SYMMETRIES_CHECK(fat.size <= 512U, "permutation too long");
    LATTICE_SYMMETRIES_CHECK(fat.masks.size() == fat.deltas.size(), "invalid fat_benes_network_t");
    LATTICE_SYMMETRIES_CHECK(fat.masks.size() <= max_depth, "fat_benes_network_t is too deep");
    std::copy(std::begin(fat.masks), std::end(fat.masks), std::begin(masks));
    std::copy(std::begin(fat.deltas), std::end(fat.deltas), std::begin(deltas));
    std::for_each(std::next(std::begin(masks), depth), std::end(masks),
                  [](auto& m) { set_zero(m); });
    std::fill(std::next(std::begin(deltas), depth), std::end(deltas), uint64_t{0});
}

big_network_t::operator fat_benes_network_t() const
{
    std::vector<ls_bits512> new_masks(depth);
    std::vector<unsigned>   new_deltas(depth);
    for (auto i = 0U; i < depth; ++i) {
        new_masks[i]  = masks[i];
        new_deltas[i] = deltas[i];
    }
    return fat_benes_network_t{std::move(new_masks), std::move(new_deltas), width};
}

auto small_network_t::operator()(uint64_t bits) const noexcept -> uint64_t
{
    benes_forward_64(bits, *this);
    return bits;
}

auto big_network_t::operator()(ls_bits512& bits) const noexcept -> void
{
    benes_forward_512(bits, *this);
}

batched_small_network_t::batched_small_network_t(
    std::array<small_network_t const*, batch_size> const& networks) noexcept
    : masks{}, deltas{}, depth{}, width{}
{
    // Make sure that it is safe to access members
    for (auto const* network : networks) {
        LATTICE_SYMMETRIES_CHECK(network != nullptr, "");
    }
    // Determine depth and width
    depth = networks[0]->depth;
    width = networks[0]->width;
    for (auto const* network : networks) {
        LATTICE_SYMMETRIES_CHECK(network->depth == depth, "");
        LATTICE_SYMMETRIES_CHECK(network->width == width, "");
    }
    // Determine masks
    for (auto i = 0U; i < depth; ++i) {
        std::transform(std::begin(networks), std::end(networks), std::data(masks[i]),
                       [i](auto const* network) { return network->masks[i]; });
    }
    for (auto i = depth; i < max_depth; ++i) {
        std::fill(std::begin(masks[i]), std::end(masks[i]), uint64_t{0});
    }
    // Determine deltas
    for (auto i = 0U; i < depth; ++i) {
        deltas[i] = networks[0]->deltas[i];
    }
    for (auto i = depth; i < max_depth; ++i) {
        deltas[i] = 0U;
    }
    std::for_each(std::next(std::begin(networks)), std::end(networks), [this](auto const* network) {
        for (auto i = 0U; i < max_depth; ++i) {
            LATTICE_SYMMETRIES_CHECK(network->deltas[i] == deltas[i], "");
        }
    });
}

auto batched_small_network_t::operator()(uint64_t bits[batch_size]) const noexcept -> void
{
    benes_forward_64(bits, *this);
}

} // namespace lattice_symmetries
