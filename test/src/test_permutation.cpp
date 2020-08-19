#include "permutation.hpp"
#include <catch2/catch.hpp>

using namespace lattice_symmetries;

TEST_CASE("checks whether a vector forms a permutation", "[permutation]")
{
    using V = std::vector<unsigned>;
    REQUIRE(is_permutation(V{}) == true);
    REQUIRE(is_permutation(V{0}) == true);
    REQUIRE(is_permutation(V{0, 0, 0}) == false);
    REQUIRE(is_permutation(V{5}) == false);
    REQUIRE(is_permutation(V{3, 1, 2, 0}) == true);
    REQUIRE(is_permutation(V{3, 1, 2, 0, 2}) == false);
}

TEST_CASE("builds standard Benes networks", "[permutation]")
{
    auto target  = std::vector<unsigned>{3, 2, 4, 1, 6, 0, 5, 7};
    auto network = compile(target);
    REQUIRE(network.has_value());
    REQUIRE(network.value().size == target.size());

    // clang-format off
    auto source = std::vector<unsigned>{
          0,  1,  2,  3,  4,  5,  6,  7,
          8,  9, 10, 11, 12, 13, 14, 15,
         16, 17, 18, 19, 20, 21, 22, 23,
         24, 25, 26, 27, 28, 29, 30, 31,
         32, 33, 34, 35, 36, 37, 38, 39,
         40, 41, 42, 43, 44, 45, 46, 47,
         48, 49, 50, 51, 52, 53, 54, 55,
         56, 57, 58, 59, 60, 61, 62, 63
    };
    target = std::vector<unsigned>{
         0,  8, 16, 24, 32, 40, 48, 56,
         1,  9, 17, 25, 33, 41, 49, 57,
         2, 10, 18, 26, 34, 42, 50, 58,
         3, 11, 19, 27, 35, 43, 51, 59,
         4, 12, 20, 28, 36, 44, 52, 60,
         5, 13, 21, 29, 37, 45, 53, 61,
         6, 14, 22, 30, 38, 46, 54, 62,
         7, 15, 23, 31, 39, 47, 55, 63
    };
    // clang-format on
    network = compile(target);
    REQUIRE(network.has_value());
    REQUIRE(network.value().size == target.size());

    target  = std::vector<unsigned>{1, 2, 3, 4, 5, 6, 0};
    network = compile(target);
    REQUIRE(network.has_value());
}
