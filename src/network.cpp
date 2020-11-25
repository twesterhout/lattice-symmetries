#include "network.hpp"
#include "kernels.hpp"
#include <algorithm>
#include <iostream>
#include <numeric>

namespace lattice_symmetries {

namespace {
    // We set 64-bit mask to
    //
    //     0000....0011....1111
    //               ~~~~~~~~~~
    //                   n
    auto init_flip_mask(bits64& mask, unsigned const n) noexcept -> void
    {
        // Play nice and do not shift by 64 bits
        // NOLINTNEXTLINE: 64 is the number of bits in bits64
        mask = n == 0U ? uint64_t{0} : ((~uint64_t{0}) >> (64U - n));
    }

    auto init_flip_mask(bits512& mask, unsigned n) noexcept -> void
    {
        auto i = 0U;
        // NOLINTNEXTLINE: 64 is the number of bits in uint64_t
        for (; n >= 64U; ++i, n -= 64) {
            mask.words[i] = ~uint64_t{0};
        }
        if (n != 0U) {
            init_flip_mask(mask.words[i], n);
            ++i;
        }
        for (; i < std::size(mask.words); ++i) {
            mask.words[i] = uint64_t{0};
        }
    }
} // namespace

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
                   [](auto& m) { return m.words[0]; });
    std::copy(std::begin(fat.deltas), std::end(fat.deltas), std::begin(deltas));
    std::fill(std::next(std::begin(masks), depth), std::end(masks), uint64_t{0});
    std::fill(std::next(std::begin(deltas), depth), std::end(deltas), uint64_t{0});
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

constexpr auto bit_permute_step(uint64_t const x, uint64_t const m, unsigned const d) noexcept
    -> uint64_t
{
    auto const y = (x ^ (x >> d)) & m;
    return x ^ y ^ (y << d);
}

auto small_network_t::operator()(uint64_t bits) const noexcept -> uint64_t
{
    for (auto i = 0U; i < depth; ++i) {
        bits = bit_permute_step(bits, masks[i], deltas[i]);
    }
    return bits;
}

auto big_network_t::operator()(bits512& bits) const noexcept -> void
{
    benes_forward(bits, static_cast<bits512 const*>(masks), depth,
                  static_cast<uint16_t const*>(deltas));
}

inline auto next_pow_of_2(uint64_t const x) noexcept -> uint64_t
{
    return x <= 1U ? uint64_t{1}
                   // NOLINTNEXTLINE: 64 is the number of bits in uint64_t, not a magic constant
                   : uint64_t{1} << static_cast<unsigned>(64 - __builtin_clzl(x - 1U));
}

template <class Network, class Int>
auto permutation_helper(Network const& network, tcb::span<Int> const bits) noexcept -> void
{
    LATTICE_SYMMETRIES_CHECK(bits.size() == network.width, "bits has wrong size");

    auto workspace = std::vector<Int>(next_pow_of_2(bits.size()));
    auto _i        = std::copy(std::begin(bits), std::end(bits), std::begin(workspace));
    std::iota(_i, std::end(workspace), bits.size());

    for (auto i = 0U; i < network.depth; ++i) {
        auto const& mask  = network.masks[i];
        auto const  delta = network.deltas[i];
        for (auto j = 0U; j < workspace.size(); ++j) {
            if (test_bit(mask, j)) {
                LATTICE_SYMMETRIES_CHECK(j + delta < workspace.size(), "");
                std::swap(workspace[j], workspace[j + delta]);
            }
        }
    }
    std::copy(std::begin(workspace),
              std::next(std::begin(workspace), static_cast<ptrdiff_t>(bits.size())),
              std::begin(bits));
}

auto small_network_t::operator()(tcb::span<unsigned> const bits) const noexcept -> void
{
    permutation_helper(*this, bits);
}

auto big_network_t::operator()(tcb::span<unsigned> const bits) const noexcept -> void
{
    permutation_helper(*this, bits);
}

template <class Network>
auto reconstruct_permutation_helper(Network const& network) -> std::vector<uint16_t>
{
    std::vector<unsigned> permutation(network.width);
    std::iota(std::begin(permutation), std::end(permutation), 0U);
    network(permutation);
    return std::vector<uint16_t>{std::begin(permutation), std::end(permutation)};
}
auto reconstruct_permutation(small_network_t const& network) -> std::vector<uint16_t>
{
    return reconstruct_permutation_helper(network);
}
auto reconstruct_permutation(big_network_t const& network) -> std::vector<uint16_t>
{
    return reconstruct_permutation_helper(network);
}

template <class Network>
auto compose_helper(Network const& x, Network const& y) -> outcome::result<Network>
{
    if (x.width != y.width) { return LS_INCOMPATIBLE_SYMMETRIES; }
    std::vector<unsigned> permutation(x.width);
    std::iota(std::begin(permutation), std::end(permutation), 0U);
    y(permutation);
    x(permutation);
    auto r = compile(permutation);
    if (!r) { return r.as_failure(); }
    return Network{std::move(r).value()};
}

auto compose(small_network_t const& x, small_network_t const& y) -> outcome::result<small_network_t>
{
    return compose_helper(x, y);
}

// auto compose(big_network_t const& x, big_network_t const& y) -> small_network_t;

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

#if 0
batched_small_network_t::batched_small_network_t(small_network_t const& small) noexcept
    : masks{}, deltas{}, flip{}, flip_mask{small.flip_mask}, depth{small.depth}
{
    for (auto i = 0; i < depth; ++i) {
        for (auto j = 0U; j < batch_size; ++j) {
            masks[i][j] = small.masks[i];
        }
        deltas[i] = small.deltas[i];
    }
    for (auto i = depth; i < max_depth; ++i) {
        for (auto j = 0U; j < batch_size; ++j) {
            masks[i][j] = 0;
        }
        deltas[i] = 0;
    }
    for (auto j = 0U; j < batch_size; ++j) {
        flip[j] = small.flip;
    }
}
#endif

auto batched_small_network_t::operator()(uint64_t bits[batch_size]) const noexcept -> void
{
    // NOLINTNEXTLINE: we do want implicit decays of arrays to pointers here
    benes_forward(bits, masks, depth, deltas);
}

} // namespace lattice_symmetries
