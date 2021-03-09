#include "bits.hpp"
#include "cpu/search_sorted.hpp"
#include "lattice_symmetries/lattice_symmetries.h"
#include <bitset>
#include <catch2/catch.hpp>
#include <complex>
#include <cstdio>
#include <iostream>
#include <memory>
#include <numeric>

TEST_CASE("obtains CPU capabilities", "[api]")
{
    if (ls_has_avx2()) { REQUIRE(ls_has_avx() == true); }
    if (ls_has_avx()) { REQUIRE(ls_has_sse4() == true); }
    static_cast<void>(ls_has_sse4());
}

TEST_CASE("toggles debug logging", "[api]")
{
    REQUIRE(ls_is_logging_enabled() == false);
    ls_enable_logging();
    REQUIRE(ls_is_logging_enabled() == true);
    ls_disable_logging();
    REQUIRE(ls_is_logging_enabled() == false);
}

auto operator<<(std::ostream& out, ls_bits512 const& x) -> std::ostream&
{
    out << '[' << std::bitset<64>(x.words[0]);
    for (auto i = 1U; i < 8; ++i) {
        out << ", " << std::bitset<64>(x.words[i]);
    }
    out << ']';
}

template <class... Args> auto make_symmetry(Args&&... args)
{
    ls_symmetry*  self   = nullptr;
    ls_error_code status = ls_create_symmetry(&self, std::forward<Args>(args)...);
    REQUIRE(status == LS_SUCCESS);
    REQUIRE(self != nullptr);
    return std::unique_ptr<ls_symmetry, void (*)(ls_symmetry*)>{self, &ls_destroy_symmetry};
}

auto check_permutation(std::function<void(ls_bits512*)> permute,
                       std::initializer_list<uint64_t>  before,
                       std::initializer_list<uint64_t>  after)
{
    REQUIRE(before.size() <= 8);
    REQUIRE(before.size() == after.size());
    ls_bits512 bits;
    lattice_symmetries::set_zero(bits);
    std::copy(before.begin(), before.end(), std::begin(bits.words));
    permute(&bits);
    ls_bits512 expected;
    lattice_symmetries::set_zero(expected);
    std::copy(after.begin(), after.end(), std::begin(expected.words));
    REQUIRE(bits == expected);
}

auto check_permutation(ls_symmetry const* symmetry, std::initializer_list<uint64_t> before,
                       std::initializer_list<uint64_t> after)
{
    check_permutation([symmetry](auto* x) { ls_apply_symmetry(symmetry, x); }, before, after);
}

auto get_eigenvalue(ls_symmetry const* symmetry) -> std::complex<double>
{
    std::complex<double> eigenvalue;
    ls_get_eigenvalue(symmetry, &eigenvalue);
    return eigenvalue;
}

TEST_CASE("constructs symmetries", "[api]")
{
    {
        unsigned const permutation[] = {1, 2, 3, 4, 5, 6, 0};
        auto const     self          = make_symmetry(std::size(permutation), permutation, 0);
        REQUIRE(ls_get_sector(self.get()) == 0);
        REQUIRE(ls_get_periodicity(self.get()) == 7);
        REQUIRE(ls_get_phase(self.get()) == 0.0);
        REQUIRE(ls_symmetry_get_number_spins(self.get()) == 7);
        REQUIRE(get_eigenvalue(self.get()) == 1.0);
        check_permutation(self.get(), {0b1010110}, {0b0101011});
    }

    {
        unsigned const permutation[] = {1, 2, 3, 0};
        auto const     self          = make_symmetry(std::size(permutation), permutation, 3);
        REQUIRE(ls_get_sector(self.get()) == 3);
        REQUIRE(ls_get_periodicity(self.get()) == 4);
        REQUIRE(ls_get_phase(self.get()) == 0.75);
        REQUIRE(ls_symmetry_get_number_spins(self.get()) == 4);
        REQUIRE(get_eigenvalue(self.get()) == std::complex{0.0, 1.0});
        check_permutation(self.get(), {0b0110}, {0b0011});
        check_permutation(self.get(), {0b1111}, {0b1111});
    }

    {
        unsigned const permutation[] = {3, 2, 1, 0};
        auto const     self          = make_symmetry(std::size(permutation), permutation, 1);
        REQUIRE(ls_get_sector(self.get()) == 1);
        REQUIRE(ls_get_periodicity(self.get()) == 2);
        REQUIRE(ls_get_phase(self.get()) == 0.5);
        REQUIRE(ls_symmetry_get_number_spins(self.get()) == 4);
        REQUIRE(get_eigenvalue(self.get()) == -1.0);
        check_permutation(self.get(), {0b0010}, {0b0100});
        check_permutation(self.get(), {0b1000}, {0b0001});
    }

    {
        std::vector<unsigned> permutation(100);
        for (auto i = uint64_t{0}; i < permutation.size(); ++i) {
            permutation[i] = (i + 3) % 100;
        }
        auto const self = make_symmetry(permutation.size(), permutation.data(), 1);
        REQUIRE(ls_get_sector(self.get()) == 1);
        REQUIRE(ls_get_periodicity(self.get()) == 100);
        REQUIRE(ls_get_phase(self.get()) == 0.01);
        REQUIRE(ls_symmetry_get_number_spins(self.get()) == 100);
        REQUIRE(get_eigenvalue(self.get()).real() == Approx(0.998026728428272));
        REQUIRE(get_eigenvalue(self.get()).imag() == Approx(-0.062790519529313));
        check_permutation(self.get(),
                          {0b0100000000000000000000000000000100000000000000000000000000000001UL,
                           0b100000000000000000000000000000000000UL},
                          {0b0000100000000000000000000000000000100000000000000000000000000000UL,
                           0b001100000000000000000000000000000000UL});
        check_permutation(self.get(),
                          {0b1100110100111110100001110110011111000110101100001110100101000110,
                           0b100110110111001010101011001010101010},
                          {0b0101100110100111110100001110110011111000110101100001110100101000,
                           0b110100110110111001010101011001010101});
    }

    {
        unsigned      permutation[] = {1, 2, 2, 0};
        ls_symmetry*  self          = nullptr;
        ls_error_code status = ls_create_symmetry(&self, std::size(permutation), permutation, 3);
        REQUIRE(status == LS_INVALID_PERMUTATION);
    }

    {
        std::vector<unsigned> permutation(1000);
        std::iota(std::begin(permutation), std::end(permutation), 0U);
        ls_symmetry*  self = nullptr;
        ls_error_code status =
            ls_create_symmetry(&self, std::size(permutation), permutation.data(), 0);
        REQUIRE(status == LS_PERMUTATION_TOO_LONG);
    }
}

template <class Deleter = decltype(&ls_destroy_symmetry)>
auto make_group(std::initializer_list<std::unique_ptr<ls_symmetry, Deleter>> generators)
{
    std::vector<ls_symmetry const*> ptrs(generators.size());
    std::transform(std::begin(generators), std::end(generators), std::begin(ptrs),
                   [](auto const& x) { return x.get(); });
    ls_group*     self   = nullptr;
    ls_error_code status = ls_create_group(&self, ptrs.size(), ptrs.data());
    REQUIRE(status == LS_SUCCESS);
    REQUIRE(self != nullptr);
    return std::unique_ptr<ls_group, void (*)(ls_group*)>{self, &ls_destroy_group};
}

TEST_CASE("constructs symmetry groups", "[api]")
{
    {
        auto const group = make_group({});
        REQUIRE(ls_get_group_size(group.get()) == 0);
        REQUIRE(ls_group_get_symmetries(group.get()) == nullptr);
        REQUIRE(ls_group_get_number_spins(group.get()) == -1);
    }

    {
        unsigned const permutation[] = {1, 2, 3, 0};
        auto           symmetry      = make_symmetry(std::size(permutation), permutation, 0);
        auto const     group         = make_group({std::move(symmetry)});
        REQUIRE(ls_get_group_size(group.get()) == 4);
        REQUIRE(ls_group_get_symmetries(group.get()) != nullptr);
        REQUIRE(ls_group_get_number_spins(group.get()) == 4);
    }

    {
        unsigned              number_spins = 70U;
        std::vector<unsigned> T;
        std::vector<unsigned> P;
        for (auto i = 0U; i < number_spins; ++i) {
            T.push_back((i + 1U) % number_spins);
            P.push_back(number_spins - 1U - i);
        }

        auto       translation = make_symmetry(T.size(), T.data(), number_spins / 2);
        auto       parity      = make_symmetry(P.size(), P.data(), 1);
        auto const group       = make_group({std::move(translation), std::move(parity)});
        REQUIRE(ls_get_group_size(group.get()) == 2 * number_spins);
        REQUIRE(ls_group_get_symmetries(group.get()) != nullptr);
        REQUIRE(ls_group_get_number_spins(group.get()) == 70);
    }

    {
        ls_group*     group  = nullptr;
        ls_error_code status = ls_create_trivial_group(&group, 25U);
        REQUIRE(status == LS_SUCCESS);
        REQUIRE(group != nullptr);
        REQUIRE(ls_get_group_size(group) == 1);
        REQUIRE(ls_group_get_symmetries(group) != nullptr);
        REQUIRE(ls_group_get_number_spins(group) == 25);
        ls_destroy_group(group);
    }
}

TEST_CASE("constructs basis", "[api]")
{
    {
        unsigned const permutation[] = {1, 2, 3, 0};
        ls_symmetry*   symmetry      = nullptr;
        ls_error_code  status =
            ls_create_symmetry(&symmetry, std::size(permutation), permutation, 0);
        REQUIRE(status == LS_SUCCESS);

        ls_symmetry const* generators[] = {symmetry};
        ls_group*          group        = nullptr;
        status = ls_create_group(&group, std::size(generators), generators);
        REQUIRE(status == LS_SUCCESS);
        ls_destroy_symmetry(symmetry);

        ls_spin_basis* basis = nullptr;
        status               = ls_create_spin_basis(&basis, group, 4, -1, 0);
        REQUIRE(status == LS_SUCCESS);
        ls_destroy_group(group);
        REQUIRE(ls_get_number_spins(basis) == 4);
        REQUIRE(ls_get_hamming_weight(basis) == -1);
        REQUIRE(ls_get_spin_inversion(basis) == 0);
        REQUIRE(ls_has_symmetries(basis) == true);

        status = ls_build(basis);
        REQUIRE(status == LS_SUCCESS);

        uint64_t count;
        status = ls_get_number_states(basis, &count);
        REQUIRE(status == LS_SUCCESS);
        REQUIRE(count == 6); // 1 + 1 + 1 + 2 + 1

        ls_states* states = nullptr;
        status            = ls_get_states(&states, basis);
        REQUIRE(status == LS_SUCCESS);

        auto const begin = ls_states_get_data(states);
        auto const end   = begin + count;
        for (auto i = 0U; i < count; ++i) {
            uint64_t index;
            status = ls_get_index(basis, begin[i], &index);
            REQUIRE(status == LS_SUCCESS);
            REQUIRE(index == i);
        }

        ls_destroy_states(states);
        ls_destroy_spin_basis(basis);
    }

    {
        ls_group*     group  = nullptr;
        ls_error_code status = ls_create_group(&group, 0, nullptr);
        REQUIRE(status == LS_SUCCESS);

        ls_spin_basis* basis = nullptr;
        status               = ls_create_spin_basis(&basis, group, 4, 2, 1);
        REQUIRE(status == LS_SUCCESS);

        ls_destroy_group(group);

        status = ls_build(basis);
        REQUIRE(status == LS_SUCCESS);

        uint64_t count;
        status = ls_get_number_states(basis, &count);
        REQUIRE(status == LS_SUCCESS);
        REQUIRE(count == 3);

        ls_destroy_spin_basis(basis);
    }
}

TEST_CASE("finds correct states", "[api]")
{
    {
        ls_error_code status;
        // clang-format off
        unsigned const sites[] =
            { 0,  1,  2,  3,  4,  5,
              6,  7,  8,  9, 10, 11,
             12, 13, 14, 15, 16, 17,
             18, 19, 20, 21, 22, 23};
        unsigned const T_x[] =
            { 1,  2,  3,  4,  5,  0,
              7,  8,  9, 10, 11,  6,
             13, 14, 15, 16, 17, 12,
             19, 20, 21, 22, 23, 18};
        unsigned const T_y[] =
            { 6,  7,  8,  9, 10, 11,
             12, 13, 14, 15, 16, 17,
             18, 19, 20, 21, 22, 23,
              0,  1,  2,  3,  4,  5};
        unsigned const P_x[] =
            {  5,  4,  3,  2,  1,  0,
              11, 10,  9,  8,  7,  6,
              17, 16, 15, 14, 13, 12,
              23, 22, 21, 20, 19, 18};
        unsigned const P_y[] =
            {18, 19, 20, 21, 22, 23,
             12, 13, 14, 15, 16, 17,
              6,  7,  8,  9, 10, 11,
              0,  1,  2,  3,  4,  5};
        // clang-format on
        ls_symmetry* T_x_symmetry = nullptr;
        status                    = ls_create_symmetry(&T_x_symmetry, std::size(T_x), T_x, 0);
        REQUIRE(status == LS_SUCCESS);

        ls_symmetry* T_y_symmetry = nullptr;
        status                    = ls_create_symmetry(&T_y_symmetry, std::size(T_y), T_y, 0);
        REQUIRE(status == LS_SUCCESS);

        ls_symmetry* P_x_symmetry = nullptr;
        status                    = ls_create_symmetry(&P_x_symmetry, std::size(P_x), P_x, 0);
        REQUIRE(status == LS_SUCCESS);

        ls_symmetry* P_y_symmetry = nullptr;
        status                    = ls_create_symmetry(&P_y_symmetry, std::size(P_y), P_y, 0);
        REQUIRE(status == LS_SUCCESS);

        ls_symmetry const* generators[] = {T_x_symmetry, T_y_symmetry, P_x_symmetry, P_y_symmetry};
        ls_group*          group        = nullptr;
        status = ls_create_group(&group, std::size(generators), generators);
        REQUIRE(status == LS_SUCCESS);
        ls_destroy_symmetry(T_x_symmetry);
        ls_destroy_symmetry(T_y_symmetry);
        ls_destroy_symmetry(P_x_symmetry);
        ls_destroy_symmetry(P_y_symmetry);

        ls_spin_basis* basis = nullptr;
        status               = ls_create_spin_basis(&basis, group, 24, 12, 1);
        REQUIRE(status == LS_SUCCESS);
        ls_destroy_group(group);

        status = ls_build(basis);
        REQUIRE(status == LS_SUCCESS);

        uint64_t count;
        status = ls_get_number_states(basis, &count);
        REQUIRE(status == LS_SUCCESS);
        REQUIRE(count == 15578U); // Computed using QuSpin with the same symmetries

        ls_states* states = nullptr;
        status            = ls_get_states(&states, basis);
        REQUIRE(status == LS_SUCCESS);

        auto const begin = ls_states_get_data(states);
        auto const end   = begin + count;
        for (auto i = 0UL; i < count; ++i) {
            uint64_t index;
            status = ls_get_index(basis, begin[i], &index);
            REQUIRE(status == LS_SUCCESS);
            REQUIRE(index == i);

            uint64_t const other_index = lattice_symmetries::search_sorted(begin, count, begin[i]);
            REQUIRE(other_index == i);
        }

        ls_destroy_states(states);
        ls_destroy_spin_basis(basis);
    }
}
