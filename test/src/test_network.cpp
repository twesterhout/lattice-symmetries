#include "network.hpp"
#include <catch2/catch.hpp>
#include <iostream>
#include <random>

using namespace lattice_symmetries;
// std::random_device device;
std::mt19937 generator{0x2340987123};

TEST_CASE("constructs small networks", "[network]")
{
    auto fat = compile(std::vector<unsigned>{3, 0, 1, 2});
    REQUIRE(fat.has_value());
    auto network = small_network_t{fat.value(), false};
    REQUIRE(network.depth == fat.value().masks.size());

    fat = compile(std::vector<unsigned>{0});
    REQUIRE(fat.has_value());
    REQUIRE(fat.value().masks.empty());
}

constexpr auto apply_symmetry_slow(uint64_t const x, tcb::span<unsigned const> const permutation,
                                   bool const flip) noexcept -> uint64_t
{
    auto r = uint64_t{0};
    for (auto i = uint64_t{0}; i < permutation.size(); ++i) {
        auto b = test_bit(x, permutation[i]);
        if (flip) { b = !b; }
        r |= static_cast<uint64_t>(b) << i;
    }
    return r;
}

auto apply_symmetry_slow(bits512 const& x, tcb::span<unsigned const> const permutation,
                         bool const flip) noexcept -> bits512
{
    bits512 r;
    set_zero(r);
    for (auto i = uint64_t{0}; i < permutation.size(); ++i) {
        auto b = test_bit(x, permutation[i]);
        if (flip) { b = !b; }
        if (b) { set_bit(r, i); }
    }
    return r;
}

auto random_bits(bits64& x, unsigned count) -> void
{
    LATTICE_SYMMETRIES_CHECK(count <= 64, "");
    if (count == 0U) { x = 0U; }
    else {
        x = std::uniform_int_distribution{0UL, ~0UL >> (64U - count)}(generator);
    }
}

auto random_bits(bits512& x, unsigned count) -> void
{
    LATTICE_SYMMETRIES_CHECK(count <= 512, "");
    set_zero(x);
    auto i = 0U;
    for (; count >= 64U; ++i, count -= 64U) {
        x.words[i] = std::uniform_int_distribution{}(generator);
    }
    if (count > 0U) { random_bits(x.words[i], count); }
}

TEST_CASE("small_network_t correctly permutes up to 6 bits", "[network]")
{
    for (auto size = 2U; size <= 6U; ++size) {
        std::vector<unsigned> source(size);
        std::vector<unsigned> target(size);
        std::iota(std::begin(target), std::end(target), 0U);
        do {
            auto const fat = compile(target);
            REQUIRE(fat.has_value());
            auto const n1 = small_network_t{fat.value(), false};
            auto const n2 = small_network_t{fat.value(), true};

            std::iota(std::begin(source), std::end(source), 0U);
            n1(source);
            REQUIRE(source == target);

            for (auto x = uint64_t{0}; x < (1U << target.size()); ++x) {
                REQUIRE(n1(x) == apply_symmetry_slow(x, target, false));
                REQUIRE(n2(x) == apply_symmetry_slow(x, target, true));
            }
        } while (std::next_permutation(std::begin(target), std::end(target)));
    }
}

TEST_CASE("small_network_t correctly permutes up to 64 bits", "[network]")
{
    for (auto size = 7U; size <= 64U; ++size) {
        std::vector<unsigned> source(size);
        std::vector<unsigned> target(size);
        std::iota(std::begin(target), std::end(target), 0U);
        for (auto i = 0U; i < 20U; ++i) {
            std::shuffle(std::begin(target), std::end(target), generator);
            auto const fat = compile(target);
            REQUIRE(fat.has_value());
            auto const n1 = small_network_t{fat.value(), false};
            auto const n2 = small_network_t{fat.value(), true};

            std::iota(std::begin(source), std::end(source), 0U);
            n1(source);
            REQUIRE(source == target);

            for (auto j = 0U; j < 20U; ++j) {
                auto const x =
                    std::uniform_int_distribution<uint64_t>{0UL, ~0UL >> (64 - size)}(generator);
                REQUIRE(n1(x) == apply_symmetry_slow(x, target, false));
                REQUIRE(n2(x) == apply_symmetry_slow(x, target, true));
            }
        }
    }
}

auto operator<<(std::ostream& out, bits512 const& xs) -> std::ostream&
{
    out << '[';
    for (auto i = 0U; i < 7U; ++i) {
        out << xs.words[i] << ", ";
    }
    out << xs.words[7] << ']';
    return out;
}

TEST_CASE("big_network_t correctly permutes bits", "[network]")
{
    for (auto size = 1U; size <= 512U; ++size) {
        std::vector<unsigned> source(size);
        std::vector<unsigned> target(size);
        std::iota(std::begin(target), std::end(target), 0U);
        for (auto i = 0U; i < 10U; ++i) {
            std::shuffle(std::begin(target), std::end(target), generator);
            auto const fat = compile(target);
            REQUIRE(fat.has_value());
            auto const n1 = big_network_t{fat.value(), false};
            auto const n2 = big_network_t{fat.value(), true};

            std::iota(std::begin(source), std::end(source), 0U);
            n1(source);
            REQUIRE(source == target);

            for (auto j = 0U; j < 20U; ++j) {
                bits512 x, y;
                random_bits(x, size);
                y = x;
                n1(y);
                REQUIRE(y == apply_symmetry_slow(x, target, false));
                y = x;
                n2(y);
                REQUIRE(y == apply_symmetry_slow(x, target, true));
            }
        }
    }
}
