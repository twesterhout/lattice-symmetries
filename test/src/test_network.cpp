#include "network.hpp"
#include <catch2/catch.hpp>
#include <iostream>

using namespace lattice_symmetries;

TEST_CASE("constructs small networks", "[network]")
{
    auto fat = compile(std::vector<unsigned>{3, 0, 1, 2});
    REQUIRE(fat.has_value());
    auto network = small_network_t{fat.value(), false};
    REQUIRE(network.depth == fat.value().masks.size());
}
