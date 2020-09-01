#include "operator.hpp"
#include "basis.hpp"
#include "bits.hpp"
#include "error_handling.hpp"
#include "lattice_symmetries/lattice_symmetries.h"
#include <omp.h>
#include <algorithm>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <numeric>
#include <span.hpp>
#include <variant>
#include <vector>

#include <cstdio>

namespace lattice_symmetries {

namespace {
    template <unsigned N> auto transpose(std::complex<double> (&matrix)[N][N]) noexcept -> void
    {
        // for (auto i = 0U; i < N; ++i) {
        //     for (auto j = 0U; j < N; ++j) {
        //         std::printf("%f+%fi, ", matrix[i][j].real(), matrix[i][j].imag());
        //     }
        //     std::printf("\n");
        // }
        for (auto i = 0U; i < N; ++i) {
            for (auto j = 0U; j < i; ++j) {
                std::swap(matrix[i][j], matrix[j][i]);
            }
        }
    }
} // namespace

template <unsigned NumberSpins> struct interaction_t {
    struct free_fn_t {
        auto operator()(void* p) const noexcept -> void { std::free(p); }
    };

    static constexpr unsigned Dim = 1U << NumberSpins;

    struct alignas(64) matrix_t {
        std::complex<double> payload[Dim][Dim];

        explicit matrix_t(std::complex<double> const* data) noexcept
        {
            std::memcpy(payload, data, Dim * Dim * sizeof(std::complex<double>));
        }

        explicit matrix_t(std::complex<double> const (&data)[Dim][Dim]) noexcept
            : matrix_t{&data[0][0]}
        {}
    };
    using sites_t = std::vector<std::array<uint16_t, NumberSpins>>;

    std::unique_ptr<matrix_t> matrix;
    sites_t                   sites;

    interaction_t(std::complex<double> const*                        _matrix,
                  tcb::span<std::array<uint16_t, NumberSpins> const> _sites)
        : matrix{std::make_unique<matrix_t>(_matrix)}, sites{std::begin(_sites), std::end(_sites)}
    {
        transpose(matrix->payload);
    }

    interaction_t(interaction_t&&) noexcept = default;
    interaction_t(interaction_t const& other)
        : matrix{std::make_unique<matrix_t>(other.matrix->payload)}, sites{other.sites}
    {}
};

namespace {
    template <class Bits, std::size_t N>
    constexpr auto gather_bits(Bits const& bits, std::array<uint16_t, N> const& indices) noexcept
        -> unsigned
    {
        // ============================= IMPORTANT ==============================
        // The order is REALLY important here. This is done to adhere to the
        // definition of cronecker product, i.e. that
        //
        //     kron(A, B) =  A00 B     A01 B     A02 B  ...
        //                   A10 B     A11 B     A12 B  ...
        //                    .
        //                    .
        //                    .
        //
        // In other words, if you change it to
        //    r |= test_bit(bits, i);
        //    r <<= 1U;
        // shit will break in really difficult to track ways...
        auto r = 0U;
        for (auto const i : indices) {
            LATTICE_SYMMETRIES_ASSERT(i < 64, "index out of bounds");
            r <<= 1U;
            r |= test_bit(bits, i);
        }
        return r;
    }

    template <class Bits, std::size_t N>
    auto scatter_bits(Bits& bits, unsigned r, std::array<uint16_t, N> const& indices) -> void
    {
        // std::printf("scatter_bits(%zu, %u, {", bits, r);
        // for (auto i : indices) {
        //     std::printf("%u, ", i);
        // }
        // std::printf("})\n");
        for (auto i = N; i-- > 0;) {
            LATTICE_SYMMETRIES_ASSERT(indices[i] < 64, "index out of bounds");
            // std::printf("set_bit_to(%zu, %u, %u)\n", bits, indices[i], r & 1U);
            set_bit_to(bits, indices[i], r & 1U);
            // std::printf("  -> %zu\n", bits);
            r >>= 1U;
        }
        // std::printf("  -> %zu\n", bits);
    }

    template <class Bits, class OffDiag> struct interaction_apply_fn_t {
        Bits const&           x;
        std::complex<double>& diagonal;
        OffDiag               off_diag;

        static constexpr bool is_noexcept = noexcept(std::declval<OffDiag const&>()(
            std::declval<Bits const&>(), std::declval<std::complex<double> const&>()));

        template <unsigned N>
        auto operator()(interaction_t<N> const& self) const noexcept(is_noexcept) -> void
        {
            // std::printf("interaction_t<%u>::operator()(%zu, %f+%fi, off_diag)\n", N, x,
            //             diagonal.real(), diagonal.imag());
            // for (auto const& data : self.matrix->payload) {
            //     std::printf("[%f, %f, %f, %f]\n", data[0].real(), data[1].real(), data[2].real(),
            //                 data[3].real());
            // }
            for (auto const edge : self.sites) {
                auto const k = gather_bits(x, edge);
                // std::printf("k = %u, edge = {%u, %u}\n", k, edge[0], edge[1]);
                auto const& data = self.matrix->payload[k];
                // std::printf("data = {%f, %f, %f, %f}\n", data[0].real(), data[1].real(),
                //             data[2].real(), data[3].real());
                for (auto n = 0U; n < std::size(data); ++n) {
                    if (data[n] != 0.0) {
                        if (n == k) {
                            // std::printf("diagonal += %f\n", data[n].real());
                            diagonal += data[n];
                        }
                        else {
                            auto y = x;
                            scatter_bits(y, n, edge);
                            // std::printf("off_diag(%zu, %f)\n", y, data[n].real());
                            off_diag(y, data[n]);
                        }
                    }
                }
            }
        }
    };

    template <unsigned N>
    constexpr auto max_index(interaction_t<N> const& interaction) noexcept -> unsigned
    {
        return std::accumulate(std::begin(interaction.sites), std::end(interaction.sites), 0U,
                               [](auto const max, auto const& sites) {
                                   return std::max<unsigned>(
                                       max, *std::max_element(std::begin(sites), std::end(sites)));
                               });
    }

} // namespace

} // namespace lattice_symmetries

using namespace lattice_symmetries;

struct ls_interaction {
    std::variant<interaction_t<1>, interaction_t<2>, interaction_t<3>, interaction_t<4>> payload;

    template <class T, class Arg, class... Args>
    ls_interaction(std::in_place_type_t<T> tag, Arg&& arg, Args&&... args)
        : payload{tag, std::forward<Arg>(arg), std::forward<Args>(args)...}
    {}
};

namespace lattice_symmetries {
namespace {
    template <class Bits, class OffDiag>
    auto apply(ls_interaction const& interaction, Bits const& spin, std::complex<double>& diagonal,
               OffDiag off_diag) noexcept(interaction_apply_fn_t<Bits, OffDiag>::is_noexcept)
        -> void
    {
        std::visit(interaction_apply_fn_t<Bits, OffDiag>{spin, diagonal, std::move(off_diag)},
                   interaction.payload);
    }

    constexpr auto max_index(ls_interaction const& interaction) noexcept -> unsigned
    {
        return std::visit([](auto const& x) { return max_index(x); }, interaction.payload);
    }

    auto max_index(tcb::span<ls_interaction const* const> interactions) noexcept -> unsigned
    {
        return std::accumulate(std::begin(interactions), std::end(interactions), 0U,
                               [](auto const max, auto const* interaction) {
                                   return std::max(max, max_index(*interaction));
                               });
    }
} // namespace
} // namespace lattice_symmetries

struct ls_operator {
    struct basis_deleter_fn_t {
        auto operator()(ls_spin_basis* p) const noexcept -> void { ls_destroy_spin_basis(p); }
    };
    using basis_ptr_t = std::unique_ptr<ls_spin_basis, basis_deleter_fn_t>;

    basis_ptr_t                 basis;
    std::vector<ls_interaction> terms;
    bool                        is_real;

    ls_operator(ls_spin_basis const* _basis, tcb::span<ls_interaction const* const> _terms)
        : basis{ls_copy_spin_basis(_basis)}, terms{}
    {
        terms.reserve(_terms.size());
        std::transform(
            std::begin(_terms), std::end(_terms), std::back_inserter(terms),
            [](auto const* x) noexcept -> auto const& { return *x; });
        is_real = lattice_symmetries::is_real(*basis)
                  && std::all_of(std::begin(terms), std::end(terms),
                                 [](auto const& x) { return ls_interaction_is_real(&x); });
    }
};

extern "C" ls_error_code ls_create_interaction1(ls_interaction** ptr, void const* matrix_2x2,
                                                unsigned const number_nodes, uint16_t const* nodes)
{
    auto p = std::make_unique<ls_interaction>(
        std::in_place_type_t<interaction_t<1>>{},
        reinterpret_cast<std::complex<double> const*>(matrix_2x2),
        tcb::span<std::array<uint16_t, 1> const>{
            reinterpret_cast<std::array<uint16_t, 1> const*>(nodes), number_nodes});
    *ptr = p.release();
    return LS_SUCCESS;
}

extern "C" ls_error_code ls_create_interaction2(ls_interaction** ptr, void const* matrix_4x4,
                                                unsigned const number_edges,
                                                uint16_t const (*edges)[2])
{
    auto p = std::make_unique<ls_interaction>(
        std::in_place_type_t<interaction_t<2>>{},
        reinterpret_cast<std::complex<double> const*>(matrix_4x4),
        tcb::span<std::array<uint16_t, 2> const>{
            reinterpret_cast<std::array<uint16_t, 2> const*>(edges), number_edges});
    *ptr = p.release();
    return LS_SUCCESS;
}

extern "C" ls_error_code ls_create_interaction3(ls_interaction** ptr, void const* matrix_8x8,
                                                unsigned const number_triangles,
                                                uint16_t const (*triangles)[3])
{
    auto p = std::make_unique<ls_interaction>(
        std::in_place_type_t<interaction_t<3>>{},
        reinterpret_cast<std::complex<double> const*>(matrix_8x8),
        tcb::span<std::array<uint16_t, 3> const>{
            reinterpret_cast<std::array<uint16_t, 3> const*>(triangles), number_triangles});
    *ptr = p.release();
    return LS_SUCCESS;
}

extern "C" ls_error_code ls_create_interaction4(ls_interaction** ptr, void const* matrix_16x16,
                                                unsigned const number_plaquettes,
                                                uint16_t const (*plaquettes)[4])
{
    auto p = std::make_unique<ls_interaction>(
        std::in_place_type_t<interaction_t<4>>{},
        reinterpret_cast<std::complex<double> const*>(matrix_16x16),
        tcb::span<std::array<uint16_t, 4> const>{
            reinterpret_cast<std::array<uint16_t, 4> const*>(plaquettes), number_plaquettes});
    *ptr = p.release();
    return LS_SUCCESS;
}

extern "C" void ls_destroy_interaction(ls_interaction* interaction)
{
    std::default_delete<ls_interaction>{}(interaction);
}

extern "C" bool ls_interaction_is_real(ls_interaction const* interaction)
{
    return std::visit(
        [](auto const& x) noexcept {
            for (auto const& column : x.matrix->payload) {
                for (auto const element : column) {
                    if (element.imag() != 0.0) { return false; }
                }
            }
            return true;
        },
        interaction->payload);
}

#if 0
extern "C" bool ls_operator_is_real(ls_operator const* op)
{
    return is_real(*op->basis)
           && std::all_of(std::begin(op->terms), std::end(op->terms),
                          [](auto const& x) { return ls_interaction_is_real(&x); });
}
#else
extern "C" bool ls_operator_is_real(ls_operator const* op) { return op->is_real; }
#endif

namespace lattice_symmetries {

template <class T, class = void> struct is_complex : std::false_type {};
template <class T>
struct is_complex<std::complex<T>, std::enable_if_t<std::is_floating_point<T>::value>>
    : std::true_type {};
template <class T> inline constexpr bool is_complex_v = is_complex<T>::value;

constexpr auto to_bits(bits512 const& x) noexcept -> uint64_t const* { return x.words; }
constexpr auto to_bits(bits512& x) noexcept -> uint64_t* { return x.words; }
constexpr auto to_bits(bits64 const& x) noexcept -> uint64_t const* { return &x; }
constexpr auto to_bits(bits64& x) noexcept -> uint64_t* { return &x; }

template <class Bits, class Callback>
auto apply_helper(ls_operator const& op, Bits const& spin, Callback callback) noexcept(
    noexcept(std::declval<Callback&>()(std::declval<Bits const&>(),
                                       std::declval<std::complex<double> const&>())))
    -> outcome::result<void>
{
    auto                 repr = spin;
    std::complex<double> eigenvalue;
    double               norm;

    ls_get_state_info(op.basis.get(), to_bits(spin), to_bits(repr), &eigenvalue, &norm);
    if (norm == 0.0) {
        // std::printf("%zu invalid state!\n", spin);
        return LS_INVALID_STATE;
    }
    auto const old_norm = norm;
    auto       diagonal = std::complex<double>{0.0, 0.0};
    auto const off_diag = [&](Bits const& x, std::complex<double> const& c) {
        ls_get_state_info(op.basis.get(), to_bits(x), to_bits(repr), &eigenvalue, &norm);
        if (norm > 0.0) { callback(repr, c * norm / old_norm * eigenvalue); }
    };
    for (auto const& term : op.terms) {
        apply(term, spin, diagonal, std::cref(off_diag));
    }
    callback(spin, diagonal);
    return outcome::success();
}

template <class T>
auto apply_helper(ls_operator const& op, uint64_t size, T const* x, T* y) noexcept
    -> outcome::result<void>
{
    OUTCOME_TRY(states, [&op]() -> outcome::result<tcb::span<uint64_t const>> {
        ls_states* states = nullptr;
        auto const status = ls_get_states(&states, op.basis.get());
        if (status != LS_SUCCESS) { return status; }
        return tcb::span<uint64_t const>{ls_states_get_data(states), ls_states_get_size(states)};
    }());
    if constexpr (!is_complex_v<T>) {
        if (!op.is_real) { return LS_OPERATOR_IS_COMPLEX; }
    }
    if (size != states.size()) { return LS_DIMENSION_MISMATCH; }

    auto const chunk_size =
        std::max(500UL, states.size() / (100UL * static_cast<unsigned>(omp_get_max_threads())));
    using acc_t = std::conditional_t<is_complex_v<T>, std::complex<double>, double>;
#pragma omp parallel for schedule(dynamic, chunk_size) default(none)                               \
    firstprivate(x, y, chunk_size, states) shared(op)
    for (auto i = uint64_t{0}; i < states.size(); ++i) {
        auto acc    = acc_t{0.0};
        auto status = apply_helper(
            op, states[i], [&acc, &op, x](auto const spin, auto const& coeff) noexcept {
                // std::printf("%zu, %f + %fi\n", spin, coeff.real(), coeff.imag());
                uint64_t   index;
                auto const _status = ls_get_index(op.basis.get(), to_bits(spin), &index);
                LATTICE_SYMMETRIES_ASSERT(_status == LS_SUCCESS, "");
                if constexpr (is_complex_v<T>) {
                    acc += std::conj(coeff) * static_cast<acc_t>(x[index]);
                }
                else {
                    acc += coeff.real() * static_cast<acc_t>(x[index]);
                }
            });
        // if (!status) { std::printf("status = %i\n", status.error()); }
        LATTICE_SYMMETRIES_ASSERT(status.has_value(), "");
        y[i] = static_cast<T>(acc);
    }
    return LS_SUCCESS;
}

template <class T>
auto apply(ls_operator const* op, uint64_t size, T const* x, T* y) noexcept -> ls_error_code
{
    auto r = apply_helper(*op, size, x, y);
    if (!r) {
        if (r.error().category() == get_error_category()) {
            return static_cast<ls_error_code>(r.error().value());
        }
        return LS_SYSTEM_ERROR;
    }
    return LS_SUCCESS;
}

} // namespace lattice_symmetries

extern "C" ls_error_code ls_create_operator(ls_operator** ptr, ls_spin_basis const* basis,
                                            unsigned const              number_terms,
                                            ls_interaction const* const terms[])
{
    auto const _terms                = tcb::span<ls_interaction const* const>{terms, number_terms};
    auto       expected_number_spins = 1U + max_index(_terms);
    if (expected_number_spins > ls_get_number_spins(basis)) { return LS_INVALID_NUMBER_SPINS; }
    auto p = std::make_unique<ls_operator>(basis, _terms);
    *ptr   = p.release();
    return LS_SUCCESS;
}

extern "C" void ls_destroy_operator(ls_operator* op) { std::default_delete<ls_operator>{}(op); }

extern "C" ls_error_code ls_operator_matvec_f32(ls_operator const* op, uint64_t size,
                                                float const* x, float* y)
{
    return apply(op, size, x, y);
}

extern "C" ls_error_code ls_operator_matvec_f64(ls_operator const* op, uint64_t size,
                                                double const* x, double* y)
{
    return apply(op, size, x, y);
}

extern "C" ls_error_code ls_operator_matvec_c64(ls_operator const* op, uint64_t size, void const* x,
                                                void* y)
{
    using C = std::complex<float>;
    return apply(op, size, static_cast<C const*>(x), static_cast<C*>(y));
}

extern "C" ls_error_code ls_operator_matvec_c128(ls_operator const* op, uint64_t size,
                                                 void const* x, void* y)
{
    using C = std::complex<double>;
    return apply(op, size, static_cast<C const*>(x), static_cast<C*>(y));
}
