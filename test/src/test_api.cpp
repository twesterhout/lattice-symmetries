#include "lattice_symmetries/lattice_symmetries.h"
#include <catch2/catch.hpp>
#include <complex>
#include <numeric>

TEST_CASE("constructs symmetries", "[api]")
{
    {
        unsigned const permutation[] = {1, 2, 3, 4, 5, 6, 0};
        ls_symmetry*   self          = nullptr;
        ls_error_code  status =
            ls_create_symmetry(&self, std::size(permutation), permutation, false, 0);
        REQUIRE(status == LS_SUCCESS);
        REQUIRE(self != nullptr);
        REQUIRE(ls_get_sector(self) == 0);
        REQUIRE(ls_get_flip(self) == false);
        REQUIRE(ls_get_periodicity(self) == 7);
        REQUIRE(ls_get_phase(self) == 0.0);
        std::complex<double> eigenvalue;
        ls_get_eigenvalue(self, &eigenvalue);
        REQUIRE(eigenvalue == 1.0);
        ls_destroy_symmetry(self);
    }

    {
        unsigned const permutation[] = {1, 2, 3, 0};
        ls_symmetry*   self          = nullptr;
        ls_error_code  status =
            ls_create_symmetry(&self, std::size(permutation), permutation, true, 3);
        REQUIRE(status == LS_SUCCESS);
        REQUIRE(self != nullptr);
        REQUIRE(ls_get_sector(self) == 3);
        REQUIRE(ls_get_flip(self) == true);
        REQUIRE(ls_get_periodicity(self) == 4);
        REQUIRE(ls_get_phase(self) == 0.75);
        ls_destroy_symmetry(self);
    }

    {
        unsigned const permutation[] = {1, 2, 0};
        ls_symmetry*   self          = nullptr;
        ls_error_code  status =
            ls_create_symmetry(&self, std::size(permutation), permutation, true, 3);
        REQUIRE(status == LS_SUCCESS);
        REQUIRE(self != nullptr);
        REQUIRE(ls_get_sector(self) == 3);
        REQUIRE(ls_get_flip(self) == true);
        REQUIRE(ls_get_periodicity(self) == 6);
        REQUIRE(ls_get_phase(self) == 0.5);
        ls_destroy_symmetry(self);
    }

    {
        unsigned      permutation[] = {1, 2, 2, 0};
        ls_symmetry*  self          = nullptr;
        ls_error_code status =
            ls_create_symmetry(&self, std::size(permutation), permutation, true, 3);
        REQUIRE(status == LS_NOT_A_PERMUTATION);
    }

    {
        std::vector<unsigned> permutation(1000);
        std::iota(std::begin(permutation), std::end(permutation), 0U);
        ls_symmetry*  self = nullptr;
        ls_error_code status =
            ls_create_symmetry(&self, std::size(permutation), permutation.data(), false, 0);
        REQUIRE(status == LS_PERMUTATION_TOO_LONG);
    }
}

TEST_CASE("constructs symmetry groups", "[api]")
{
    {
        unsigned const permutation[] = {1, 2, 3, 0};
        ls_symmetry*   symmetry      = nullptr;
        ls_group*      group         = nullptr;
        ls_error_code  status =
            ls_create_symmetry(&symmetry, std::size(permutation), permutation, false, 0);
        REQUIRE(status == LS_SUCCESS);

        ls_symmetry const* generators[] = {symmetry};
        status = ls_create_group(&group, std::size(generators), generators);
        REQUIRE(status == LS_SUCCESS);
        REQUIRE(ls_get_group_size(group) == 4);
        ls_destroy_group(group);
        ls_destroy_symmetry(symmetry);
    }
}
