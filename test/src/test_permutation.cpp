#include "permutation.hpp"
#include <catch2/catch.hpp>
#include <iostream>

using namespace lattice_symmetries;

template <class T> auto operator==(std::vector<T> const& xs, std::vector<T> const& ys) -> bool
{
    if (xs.size() != ys.size()) { return false; }
    return std::equal(std::begin(xs), std::end(xs), std::begin(ys));
}

auto to_vector(bits64 const xs, unsigned const size, std::vector<int>& out) -> void
{
    LATTICE_SYMMETRIES_CHECK(size <= 64U, "index out of bounds");
    for (auto i = 0; i < size; ++i) {
        out.push_back((xs >> i) & 1);
    }
}
auto to_vector(bits64 const xs, unsigned const size = 64U) -> std::vector<int>
{
    std::vector<int> out;
    to_vector(xs, size, out);
    return out;
}
auto to_vector(bits512 const& xs, unsigned size = 512U) -> std::vector<int>
{
    LATTICE_SYMMETRIES_CHECK(size <= 512U, "index out of bounds");
    std::vector<int> out;
    for (auto const w : xs.words) {
        auto const n = std::min(64U, size);
        to_vector(w, n, out);
        size -= n;
        if (size == 0U) { break; }
    }
    return out;
}

template <class T1, class T2>
auto operator<<(std::ostream& out, std::pair<T1, T2> const& xs) -> std::ostream&
{
    return out << '(' << xs.first << ", " << xs.second << ')';
}

auto operator<<(std::ostream& out, bits512 const& x) -> std::ostream&
{
    for (auto i = 0; i < 8; ++i) {
        out << static_cast<int>(test_bit(x, i)) << ",";
    }
    return out;
}
// auto operator<<(std::ostream& out, side_t const x) -> std::ostream&
// {
//     switch (x) {
//     case side_t::left: out << "left"; break;
//     case side_t::right: out << "right"; break;
//     }
//     return out;
// }
// auto operator<<(std::ostream& out, node_t const& x) -> std::ostream&
// {
//     return out << "node{" << x.first << ", " << x.second << ", " << x.side << '}';
// }

template <class T> auto operator<<(std::ostream& out, std::vector<T> const& xs) -> std::ostream&
{
    auto const write = [&out](auto const& x) {
        if constexpr (std::is_same_v<T, uint8_t>) { out << static_cast<int>(x); }
        else {
            out << x;
        }
    };
    out << '[';
    if (!xs.empty()) {
        write(xs[0]);
        for (auto i = 1U; i < xs.size(); ++i) {
            out << ", ";
            write(xs[i]);
        }
    }
    out << ']';
    return out;
}

TEST_CASE("checks whether a vector forms a permutation", "[permutation]")
{
    REQUIRE(is_permutation({}) == true);
    REQUIRE(is_permutation({0}) == true);
    REQUIRE(is_permutation({5}) == false);
    REQUIRE(is_permutation({3, 1, 2, 0}) == true);
    REQUIRE(is_permutation({3, 1, 2, 0, 2}) == false);
}

// TEST_CASE("finds all nonoverlapping pairs", "[permutation]")
// {
//     using P = std::vector<std::pair<int, int>>;
//     using S = std::vector<int>;
//     REQUIRE(nonoverlapping_pairs(4, 1) == std::make_tuple(P{{0, 1}, {2, 3}}, S{}));
//     REQUIRE(nonoverlapping_pairs(4, 2) == std::make_tuple(P{{0, 2}, {1, 3}}, S{}));
//     REQUIRE(nonoverlapping_pairs(6, 1) == std::make_tuple(P{{0, 1}, {2, 3}, {4, 5}}, S{}));
//     // REQUIRE(nonoverlapping_pairs(6, 3) == R{{0, 3}, {1, 4}, {2, 5}});
//     // REQUIRE(nonoverlapping_pairs(8, 1) == R{{0, 1}, {2, 3}, {4, 5}, {6, 7}});
//     // REQUIRE(nonoverlapping_pairs(8, 2) == R{{0, 2}, {1, 3}, {4, 6}, {5, 7}});
//     // REQUIRE(nonoverlapping_pairs(8, 4) == R{{0, 4}, {1, 5}, {2, 6}, {3, 7}});
//     // std::cout << nonoverlapping_pairs(64, 7) << '\n';
// }

// TEST_CASE("construct Θ mask for δ-swap", "[permutation]")
// {
//     REQUIRE((to_vector(mask_for_delta_swap({1, 5, 3, 9, 2, 7}, {1, 5, 9, 3, 2, 7}, 1), 6)
//              == std::vector{0, 0, 1, 0, 0, 0}));
// }

// TEST_CASE("correctly solves cycles", "[permutation]")
// {
//     std::vector<node_t> cycle = {
//         {0, 1, side_t::left},  // ab
//         {6, 0, side_t::right}, // EF
//         {6, 7, side_t::left},  // gh
//         {5, 7, side_t::right}, // GH
//         {4, 5, side_t::left},  // ef
//         {4, 1, side_t::right}, // CD
//     };
//     REQUIRE((solve_cycle(cycle) == std::vector<uint8_t>{0, 1, 1, 1, 0, 0}));
//     REQUIRE((solve_cycle(cycle, /*start_with=*/1) == std::vector<uint8_t>{1, 0, 0, 0, 1, 1}));
//     cycle = {
//         {2, 3, side_t::left},
//         {3, 2, side_t::right},
//     };
//     REQUIRE((solve_cycle(cycle) == std::vector<uint8_t>{0, 1}));
//     REQUIRE((solve_cycle(cycle, /*start_with=*/1) == std::vector<uint8_t>{1, 0}));
// }

TEST_CASE("correctly finds cycles", "[permutation]")
{
    // auto find_cycles(std::vector<int> const& source, std::vector<int> const& target, int delta)
    auto const source = std::vector<int>{0, 1, 2, 3, 4, 5, 6, 7};
    auto const target = std::vector<int>{3, 2, 4, 1, 6, 0, 5, 7};

    // auto cycles        = find_cycles(source, target, 1);
    // auto [left, right] = solve_cycles(source, target, cycles);
    // std::cout << left << '\n';
    // std::cout << right << '\n';

    auto network = compile(target);
    std::cout << network->fwd_masks << '\n';
    std::cout << network->bwd_masks << '\n';
    std::cout << (*network)({0, 1, 2, 3, 4, 5, 6, 7}) << '\n';

    network->optimize();
    std::cout << network->fwd_masks << '\n';
    std::cout << network->bwd_masks << '\n';
    std::cout << (*network)({0, 1, 2, 3, 4, 5, 6, 7}) << '\n';

    REQUIRE(compile(target, 2) == std::nullopt);
    REQUIRE(compile({2, 1, 0, 3, 4, 5, 6, 7}, 2).has_value());
    // for_each_pair(8, 1, [](auto const i, auto const j) { std::cout << i << ", " << j << '\n'; });
    // auto cycle    = build_cycle(node_t{0, 1, side_t::left}, source, target, pairs);
    // auto expected = std::vector<node_t>{
    //     {0, 1, side_t::left},  // ab
    //     {6, 0, side_t::right}, // EF
    //     {6, 7, side_t::left},  // gh
    //     {5, 7, side_t::right}, // GH
    //     {4, 5, side_t::left},  // ef
    //     {4, 1, side_t::right}, // CD
    // };
    // REQUIRE(cycle == expected);
}
