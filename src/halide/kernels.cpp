#include "kernels.hpp"
#include <HalideBuffer.h>
#include <HalideRuntime.h>
#include <cstdlib>
#include <cstring>

#define IS_X86_64 true // false
//
// Architecture-specific kernels
#if IS_X86_64
#    include "ls_internal_is_representative_general_kernel_64_avx.h"
#    include "ls_internal_is_representative_general_kernel_64_avx2.h"
#    include "ls_internal_is_representative_general_kernel_64_sse41.h"
#    include "ls_internal_state_info_general_kernel_64_avx.h"
#    include "ls_internal_state_info_general_kernel_64_avx2.h"
#    include "ls_internal_state_info_general_kernel_64_sse41.h"

#    include "ls_internal_is_representative_symmetric_kernel_64_avx.h"
#    include "ls_internal_is_representative_symmetric_kernel_64_avx2.h"
#    include "ls_internal_is_representative_symmetric_kernel_64_sse41.h"
#    include "ls_internal_state_info_symmetric_kernel_64_avx.h"
#    include "ls_internal_state_info_symmetric_kernel_64_avx2.h"
#    include "ls_internal_state_info_symmetric_kernel_64_sse41.h"

#    include "ls_internal_is_representative_antisymmetric_kernel_64_avx.h"
#    include "ls_internal_is_representative_antisymmetric_kernel_64_avx2.h"
#    include "ls_internal_is_representative_antisymmetric_kernel_64_sse41.h"
#    include "ls_internal_state_info_antisymmetric_kernel_64_avx.h"
#    include "ls_internal_state_info_antisymmetric_kernel_64_avx2.h"
#    include "ls_internal_state_info_antisymmetric_kernel_64_sse41.h"
#endif
// Generic implementation
#include "ls_internal_is_representative_antisymmetric_kernel_64.h"
#include "ls_internal_is_representative_general_kernel_64.h"
#include "ls_internal_is_representative_symmetric_kernel_64.h"
#include "ls_internal_state_info_antisymmetric_kernel_64.h"
#include "ls_internal_state_info_general_kernel_64.h"
#include "ls_internal_state_info_symmetric_kernel_64.h"

namespace lattice_symmetries {

typedef int (*ls_internal_state_info_general_kernel_t)(
    struct halide_buffer_t* _x_buffer, uint64_t _flip_mask, struct halide_buffer_t* _masks_buffer,
    struct halide_buffer_t* _eigvals_re_buffer, struct halide_buffer_t* _eigvals_im_buffer,
    struct halide_buffer_t* _shifts_buffer, struct halide_buffer_t* _representative_buffer,
    struct halide_buffer_t* _character_buffer, struct halide_buffer_t* _norm_buffer);

struct halide_kernels_list_t {
    ls_internal_state_info_general_kernel_t general;
    ls_internal_state_info_general_kernel_t symmetric;
    ls_internal_state_info_general_kernel_t antisymmetric;
};

enum class proc_arch { generic, sse4_1, avx, avx2 };

inline auto current_architecture() -> proc_arch
{
    if (auto const* arch = std::getenv("LATTICE_SYMMETRIES_ARCH")) {
        if (std::strcmp(arch, "generic") == 0) {
            LATTICE_SYMMETRIES_LOG_DEBUG("%s\n", "Kernels will use no special instructions...");
            return proc_arch::generic;
        }
        if (std::strcmp(arch, "sse4_1") == 0) {
            LATTICE_SYMMETRIES_LOG_DEBUG("%s\n", "Kernels will use SSE4.1 instructions...");
            return proc_arch::sse4_1;
        }
        if (std::strcmp(arch, "avx") == 0) {
            LATTICE_SYMMETRIES_LOG_DEBUG("%s\n", "Kernels will use AVX instructions...");
            return proc_arch::avx;
        }
        if (std::strcmp(arch, "avx2") == 0) {
            LATTICE_SYMMETRIES_LOG_DEBUG("%s\n", "Kernels will use AVX2 instructions...");
            return proc_arch::avx2;
        }
    }

    __builtin_cpu_init();
    if (__builtin_cpu_supports("avx2") > 0) {
        LATTICE_SYMMETRIES_LOG_DEBUG("%s\n", "Kernels will use AVX2 instructions...");
        return proc_arch::avx2;
    }
    if (__builtin_cpu_supports("avx") > 0) {
        LATTICE_SYMMETRIES_LOG_DEBUG("%s\n", "Kernels will use AVX instructions...");
        return proc_arch::avx;
    }
    if (__builtin_cpu_supports("sse4.1") > 0) {
        LATTICE_SYMMETRIES_LOG_DEBUG("%s\n", "Kernels will use SSE4.1 instructions...");
        return proc_arch::sse4_1;
    }
    LATTICE_SYMMETRIES_LOG_DEBUG("%s\n", "Kernels will use no special instructions...");
    return proc_arch::generic;
}

#if IS_X86_64
#    define KERNEL_FOR_ARCH(kernel_name, current_arch)                                             \
        [](proc_arch _arch923) {                                                                   \
            switch (_arch923) {                                                                    \
            case proc_arch::generic: return kernel_name;                                           \
            case proc_arch::sse4_1: return kernel_name##_sse41;                                    \
            case proc_arch::avx: return kernel_name##_avx;                                         \
            case proc_arch::avx2: return kernel_name##_avx2;                                       \
            default: LATTICE_SYMMETRIES_CHECK(false, "Unsupported architecture");                  \
            }                                                                                      \
        }(current_arch)
#else
#    define KERNEL_FOR_ARCH(kernel_name, current_arch)                                             \
        [](proc_arch _arch923) {                                                                   \
            switch (_arch923) {                                                                    \
            case proc_arch::generic: return kernel_name;                                           \
            default: LATTICE_SYMMETRIES_CHECK(false, "Unsupported architecture");                  \
            }                                                                                      \
        }(current_arch)
#endif

namespace {
    halide_kernels_list_t init_halide_kernels()
    {
        __builtin_cpu_init();
#if IS_X86_64
        if (__builtin_cpu_supports("avx2") > 0) {
            return {&ls_internal_state_info_general_kernel_64_avx2,
                    &ls_internal_state_info_symmetric_kernel_64_avx2,
                    &ls_internal_state_info_antisymmetric_kernel_64_avx2};
        }
        if (__builtin_cpu_supports("avx") > 0) {
            return {&ls_internal_state_info_general_kernel_64_avx,
                    &ls_internal_state_info_symmetric_kernel_64_avx,
                    &ls_internal_state_info_antisymmetric_kernel_64_avx};
        }
        if (__builtin_cpu_supports("sse4.1") > 0) {
            return {&ls_internal_state_info_general_kernel_64_sse41,
                    &ls_internal_state_info_symmetric_kernel_64_sse41,
                    &ls_internal_state_info_antisymmetric_kernel_64_sse41};
        }
#endif
        return {&ls_internal_state_info_general_kernel_64,
                &ls_internal_state_info_symmetric_kernel_64,
                &ls_internal_state_info_antisymmetric_kernel_64};
    }
} // namespace

struct halide_kernel_state {
    mutable halide_buffer_t _masks;
    mutable halide_buffer_t _eigvals_re;
    mutable halide_buffer_t _eigvals_im;
    mutable halide_buffer_t _shifts;
    halide_dimension_t      _masks_dims[2];
    halide_dimension_t      _shifts_dim;
    uint64_t                _flip_mask;

  private:
    static auto get_flip_mask_64(unsigned const n) noexcept -> uint64_t
    {
        return n == 0U ? uint64_t{0} : ((~uint64_t{0}) >> (64U - n));
    }

  public:
    explicit halide_kernel_state(ls_flat_spin_basis const& basis) noexcept
        : _masks{
            /*device=*/0,
            /*device_interface=*/nullptr,
            /*host=*/reinterpret_cast<uint8_t*>(basis.group.masks.get()),
            /*flags=*/0,
            /*type=*/halide_type_t{halide_type_uint, 64, 1},
            /*dimensions=*/2,
            /*dim=*/&(_masks_dims[0]),
            /*padding=*/nullptr,
        },
        _eigvals_re{
            /*device=*/0,
            /*device_interface=*/nullptr,
            /*host=*/reinterpret_cast<uint8_t*>(basis.group.eigenvalues_real.get()),
            /*flags=*/0,
            /*type=*/halide_type_t{halide_type_float, 64, 1},
            /*dimensions=*/1,
            /*dim=*/&(_masks_dims[1]),
            /*padding=*/nullptr,
        },
        _eigvals_im{
            /*device=*/0,
            /*device_interface=*/nullptr,
            /*host=*/reinterpret_cast<uint8_t*>(basis.group.eigenvalues_imag.get()),
            /*flags=*/0,
            /*type=*/halide_type_t{halide_type_float, 64, 1},
            /*dimensions=*/1,
            /*dim=*/&(_masks_dims[1]),
            /*padding=*/nullptr,
        },
        _shifts{
            /*device=*/0,
            /*device_interface=*/nullptr,
            /*host=*/reinterpret_cast<uint8_t*>(basis.group.shifts.get()),
            /*flags=*/0,
            /*type=*/halide_type_t{halide_type_uint, 64, 1},
            /*dimensions=*/1,
            /*dim=*/&_shifts_dim,
            /*padding=*/nullptr,
        },
        _masks_dims{
            halide_dimension_t{
                /*min=*/0,
                /*extent=*/static_cast<int32_t>(basis.group.shape[0]),
                /*stride=*/static_cast<int32_t>(basis.group.shape[1]),
                /*flags=*/0
            },
            halide_dimension_t{
                /*min=*/0,
                /*extent=*/static_cast<int32_t>(basis.group.shape[1]), /*stride=*/1, /*flags=*/0},
        },
        _shifts_dim{
            /*min=*/0,
            /*extent=*/static_cast<int32_t>(basis.group.shape[0]),
            /*stride=*/1,
            /*flags=*/0
        },
        _flip_mask{get_flip_mask_64(basis.number_spins)}
    {}

    halide_kernel_state(halide_kernel_state const& other) noexcept
        : _masks{other._masks}
        , _eigvals_re{other._eigvals_re}
        , _eigvals_im{other._eigvals_im}
        , _shifts{other._shifts}
        , _masks_dims{other._masks_dims[0], other._masks_dims[1]}
        , _shifts_dim{other._shifts_dim}
        , _flip_mask{other._flip_mask}
    {
        _masks.dim      = &(_masks_dims[0]);
        _eigvals_re.dim = &(_masks_dims[1]);
        _eigvals_im.dim = &(_masks_dims[1]);
        _shifts.dim     = &_shifts_dim;
    }

    auto operator=(halide_kernel_state const&) -> halide_kernel_state& = delete;
    auto operator=(halide_kernel_state&&) -> halide_kernel_state& = delete;
};

struct is_representative_no_symmetries {
    auto operator()(uint64_t const count, void const* /*x*/, uint8_t* is_repr,
                    double* norm) const noexcept -> void
    {
        std::fill(is_repr, is_repr + count, uint8_t{1});
        std::fill(norm, norm + count, 1.0);
    }
};

struct halide_is_representative_kernel {
    typedef int (*ls_internal_is_representative_kernel_t)(
        struct halide_buffer_t* /*_x_buffer*/, uint64_t /*_flip_mask*/,
        struct halide_buffer_t* /*_masks_buffer*/, struct halide_buffer_t* /*_shifts_buffer*/,
        struct halide_buffer_t* /*_eigvals_re_buffer*/,
        struct halide_buffer_t* /*_is_representative_buffer*/,
        struct halide_buffer_t* /*_norm_buffer*/);

    halide_kernel_state                    state;
    ls_internal_is_representative_kernel_t kernel;

    explicit halide_is_representative_kernel(ls_flat_spin_basis const& basis) noexcept
        : state{basis}, kernel{}
    {
        LATTICE_SYMMETRIES_LOG_DEBUG("%s\n", "Constructing halide_is_representative_kernel...");
        auto const arch = current_architecture();
        switch (basis.spin_inversion) {
        case 0:
            kernel = KERNEL_FOR_ARCH(ls_internal_is_representative_general_kernel_64, arch);
            break;
        case 1:
            kernel = KERNEL_FOR_ARCH(ls_internal_is_representative_symmetric_kernel_64, arch);
            break;
        case -1:
            kernel = KERNEL_FOR_ARCH(ls_internal_is_representative_antisymmetric_kernel_64, arch);
            break;
        default: LATTICE_SYMMETRIES_CHECK(false, "invalid spin_inversion");
        } // end switch
    }

    auto operator()(uint64_t const count, void const* x, uint8_t* is_repr,
                    double* norm) const noexcept -> void
    {
        halide_dimension_t batch_dim{/*min=*/0, /*extent=*/static_cast<int32_t>(count),
                                     /*stride=*/1, /*flags=*/0};
        halide_buffer_t    x_buf{
            /*device=*/0,
            /*device_interface=*/nullptr,
            /*host=*/static_cast<uint8_t*>(const_cast<void*>(x)),
            /*flags=*/0,
            /*type=*/halide_type_t{halide_type_uint, 64, 1},
            /*dimensions=*/1,
            /*dim=*/&batch_dim,
            /*padding=*/nullptr,
        };
        halide_buffer_t is_repr_buf{
            /*device=*/0,
            /*device_interface=*/nullptr,
            /*host=*/is_repr,
            /*flags=*/0,
            /*type=*/halide_type_t{halide_type_uint, 8, 1},
            /*dimensions=*/1,
            /*dim=*/&batch_dim,
            /*padding=*/nullptr,
        };
        halide_buffer_t norm_buf{
            /*device=*/0,
            /*device_interface=*/nullptr,
            /*host=*/reinterpret_cast<uint8_t*>(norm),
            /*flags=*/0,
            /*type=*/halide_type_t{halide_type_float, 64, 1},
            /*dimensions=*/1,
            /*dim=*/&batch_dim,
            /*padding=*/nullptr,
        };
        (*kernel)(&x_buf, state._flip_mask, &state._masks, &state._eigvals_re, &state._shifts,
                  &is_repr_buf, &norm_buf);
    }
};

struct state_info_no_symmetries {
    unsigned const number_words;

    auto operator()(uint64_t const count, void const* x, void* repr,
                    std::complex<double>* character, double* norm) const noexcept -> void
    {
        std::memcpy(repr, x, count * number_words * sizeof(uint64_t));
        std::fill(character, character + count, std::complex<double>{1.0});
        std::fill(norm, norm + count, 1.0);
    }
};

struct halide_state_info_kernel {
    typedef int (*ls_internal_state_info_kernel_t)(
        struct halide_buffer_t* /*_x_buffer*/, uint64_t /*_flip_mask*/,
        struct halide_buffer_t* /*_masks_buffer*/, struct halide_buffer_t* /*_eigvals_re_buffer*/,
        struct halide_buffer_t* /*_eigvals_im_buffer*/, struct halide_buffer_t* /*_shifts_buffer*/,
        struct halide_buffer_t* /*_representative_buffer*/,
        struct halide_buffer_t* /*_character_buffer*/, struct halide_buffer_t* /*_norm_buffer*/);

    halide_kernel_state             state;
    ls_internal_state_info_kernel_t kernel;

    explicit halide_state_info_kernel(ls_flat_spin_basis const& basis) noexcept
        : state{basis}, kernel{}
    {
        auto const arch = current_architecture();
        switch (basis.spin_inversion) {
        case 0: kernel = KERNEL_FOR_ARCH(ls_internal_state_info_general_kernel_64, arch); break;
        case 1: kernel = KERNEL_FOR_ARCH(ls_internal_state_info_symmetric_kernel_64, arch); break;
        case -1:
            kernel = KERNEL_FOR_ARCH(ls_internal_state_info_antisymmetric_kernel_64, arch);
            break;
        default: LATTICE_SYMMETRIES_CHECK(false, "invalid spin_inversion");
        } // end switch
    }

    auto operator()(uint64_t const count, void const* x, void* repr,
                    std::complex<double>* character, double* norm) const noexcept -> void
    {
        halide_dimension_t batch_dim{/*min=*/0, /*extent=*/static_cast<int32_t>(count),
                                     /*stride=*/1, /*flags=*/0};
        halide_dimension_t character_dims[2] = {
            halide_dimension_t{/*min=*/0, /*extent=*/static_cast<int32_t>(count), /*stride=*/2,
                               /*flags=*/0},
            halide_dimension_t{/*min=*/0, /*extent=*/2, /*stride=*/1,
                               /*flags=*/0}};
        halide_buffer_t x_buf{
            /*device=*/0,
            /*device_interface=*/nullptr,
            /*host=*/static_cast<uint8_t*>(const_cast<void*>(x)),
            /*flags=*/0,
            /*type=*/halide_type_t{halide_type_uint, 64, 1},
            /*dimensions=*/1,
            /*dim=*/&batch_dim,
            /*padding=*/nullptr,
        };
        halide_buffer_t repr_buf{
            /*device=*/0,
            /*device_interface=*/nullptr,
            /*host=*/static_cast<uint8_t*>(repr),
            /*flags=*/0,
            /*type=*/halide_type_t{halide_type_uint, 64, 1},
            /*dimensions=*/1,
            /*dim=*/&batch_dim,
            /*padding=*/nullptr,
        };
        halide_buffer_t character_buf{
            /*device=*/0,
            /*device_interface=*/nullptr,
            /*host=*/reinterpret_cast<uint8_t*>(character),
            /*flags=*/0,
            /*type=*/halide_type_t{halide_type_float, 64, 1},
            /*dimensions=*/2,
            /*dim=*/character_dims,
            /*padding=*/nullptr,
        };
        halide_buffer_t norm_buf{
            /*device=*/0,
            /*device_interface=*/nullptr,
            /*host=*/reinterpret_cast<uint8_t*>(norm),
            /*flags=*/0,
            /*type=*/halide_type_t{halide_type_float, 64, 1},
            /*dimensions=*/1,
            /*dim=*/&batch_dim,
            /*padding=*/nullptr,
        };
        LATTICE_SYMMETRIES_LOG_DEBUG("Extent: %i\n", state._shifts.dim[0].extent);
        (*kernel)(&x_buf, state._flip_mask, &state._masks, &state._eigvals_re, &state._eigvals_im,
                  &state._shifts, &repr_buf, &character_buf, &norm_buf);
    }
};

auto make_state_info_kernel(ls_flat_spin_basis const& basis) noexcept -> state_info_kernel_type
{
    if (basis.group.shape[0] > 0) { // Has at least one symmetry
        return halide_state_info_kernel{basis};
    }
    auto const number_words = (basis.number_spins + 64 - 1) / 64;
    return state_info_no_symmetries{number_words};
}

auto make_is_representative_kernel(ls_flat_spin_basis const& basis) noexcept
    -> is_representative_kernel_type
{
    if (basis.group.shape[0] > 0) { // Has at least one symmetry
        return halide_is_representative_kernel{basis};
    }
    return is_representative_no_symmetries{};
}

} // namespace lattice_symmetries
