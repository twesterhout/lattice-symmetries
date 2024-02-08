#include "lattice_symmetries_types.h"
#include <Halide.h>
#include <callback.h>
#include <cstdint>

using namespace Halide;

namespace lattice_symmetries {

auto mk_fixed_hamming_state_to_index_kernel(int const number_sites, int const hamming_weight,
                                            halide_buffer_t const *_binomials) -> Halide::Callable {
    ImageParam _x{halide_type_of<uint64_t>(), 1, "x"};
    Func out{"out"};

    if (hamming_weight == 0 || hamming_weight == number_sites) {
        Var batch_idx{"batch_idx"};
        out(batch_idx) = cast<int64_t>(0);
    } else {
        auto binomials = Buffer{Runtime::Buffer<uint64_t, 2>{*_binomials}.copy(), "binomials"};

        Var batch_idx{"batch_idx"};
        RDom k{0, min(Expr{hamming_weight}, binomials.dim(0).extent() - 1), "k"};

        Func temp{"temp"};

        temp(batch_idx) = Tuple{cast<int64_t>(0), _x(batch_idx)};
        auto index = temp(batch_idx)[0];
        auto state = temp(batch_idx)[1];
        auto n = cast<int32_t>(count_trailing_zeros(state));
        // min(cast<int32_t>(count_trailing_zeros(state)), binomials.dim(1).extent() - 1);
        temp(batch_idx) = {index + binomials(k + 1, n), state & (state - 1)};

        out(batch_idx) = temp(batch_idx)[0];

        Var i_inner{"i_inner"}, i_outer{"i_outer"};

        auto *const vector_size_str = getenv("LS_S2I_BLOCK_SIZE");
        auto vector_size = 8;
        if (vector_size_str != NULL) {
            vector_size = std::stoi(vector_size_str);
        }

        auto *const unroll_str = getenv("LS_S2I_UNROLL");
        auto unroll = true;
        if (unroll_str != NULL) {
            unroll = std::stoi(unroll_str) != 0;
        }

        out.split(batch_idx, i_outer, i_inner, vector_size, TailStrategy::GuardWithIf)
            .vectorize(i_inner);
        temp.compute_at(out, i_inner).store_in(MemoryType::Register);
        if (unroll) {
            temp.update(0).unroll(k);
        }
    }

    _x.dim(0).set_min(0).set_stride(1);
    out.output_buffer().dim(0).set_min(0).set_stride(1).set_extent(_x.dim(0).extent());

    // out.print_loop_nest();
    auto target = get_jit_target_from_environment();
    target.set_features({Target::Feature::NoAsserts, Target::Feature::NoBoundsQuery});
    // out.compile_to_lowered_stmt("state_to_index.html", {_x}, Halide::HTML, target);
    // out.compile_to_assembly("state_to_index.S", {_x}, "fixed_hamming_state_to_index", target);
    return out.compile_to_callable({_x}, target);
}

auto invoke_fixed_hamming_state_to_index_kernel(void *data, va_alist alist) noexcept -> void {
    auto *ctx = static_cast<std::tuple<Halide::Callable, void *> *>(data);
    auto const &fn = std::get<0>(*ctx);
    va_start_void(alist);
    auto *basis_states = va_arg_ptr(alist, halide_buffer_t *);
    auto *indices = va_arg_ptr(alist, halide_buffer_t *);
    try {
        fn(basis_states, indices);
    } catch (RuntimeError &e) {
        fprintf(stderr, "Halide::RuntimeError: %s\n", e.what());
        abort();
    } catch (InternalError &e) {
        fprintf(stderr, "Halide::InternalError: %s\n", e.what());
        abort();
    }
    va_return_void(alist);
}

auto bit_permute_step_64(Expr x, Expr m, Expr d) -> Expr {
    auto y = ((x >> d) ^ x) & m;
    return (x ^ y) ^ (y << d);
}

auto mk_is_representative_kernel(halide_buffer_t const *_masks, halide_buffer_t const *_eigvals_re,
                                 halide_buffer_t const *_shifts, int spin_inversion)
    -> Halide::Callable {
    auto masks = Halide::Buffer{Halide::Runtime::Buffer<uint64_t, 2>{*_masks}.copy(), "masks"};
    auto eigvals_re =
        Halide::Buffer{Halide::Runtime::Buffer<double, 1>{*_eigvals_re}.copy(), "eigvals_re"};
    auto shifts = Halide::Buffer{Halide::Runtime::Buffer<uint64_t, 1>{*_shifts}.copy(), "shifts"};
    // auto const flip_mask = 0;

    ImageParam _x{halide_type_of<uint64_t>(), 1, "x"};
    // ImageParam _is_representative{"is_representative", 1};
    // ImageParam _norm{"norm", 1};

    auto count = _x.dim(0).extent();
    auto depth = masks.dim(0).extent();
    auto number_masks = masks.dim(1).extent();
    auto const chunk_size = 4; // natural_vector_size(type_of<uint64_t>());
    auto number_chunks = number_masks / chunk_size;
    auto number_rest = number_masks - number_chunks * chunk_size;

    Var i{"i"};
    Var j_outer{"j_outer"};
    Var j_inner{"j_inner"};

    // Transform input
    RDom k{0, depth, "k"};
    Func y{"y"};
    y(i, j_outer, j_inner) = _x(i);
    y(i, j_outer, j_inner) = bit_permute_step_64(
        y(i, j_outer, j_inner), masks(k, j_outer * chunk_size + j_inner), shifts(k));
    Func temp{"temp"};
    if (spin_inversion == 0) {
        temp(i, j_outer, j_inner) = Tuple{cast<double>(y(i, j_outer, j_inner) == _x(i)) *
                                              eigvals_re(j_outer * chunk_size + j_inner),
                                          y(i, j_outer, j_inner) >= _x(i)};
        // } else if (spin_inversion == 1) {
        //   auto y_flipped = y(i, j_outer, j_inner) ^ flip_mask;
        //   temp(i, j_outer, j_inner) =
        //       Tuple{(cast<double>(y(i, j_outer, j_inner) == _x(i)) +
        //              cast<double>(y_flipped == _x(i))) *
        //                 eigvals_re(j_outer * chunk_size + j_inner),
        //             y(i, j_outer, j_inner) >= _x(i) && y_flipped >= _x(i)};
        // } else if (spin_inversion == -1) {
        //   auto y_flipped = y(i, j_outer, j_inner) ^ flip_mask;
        //   temp(i, j_outer, j_inner) =
        //       Tuple{(cast<double>(y(i, j_outer, j_inner) == _x(i)) -
        //              cast<double>(y_flipped == _x(i))) *
        //                 eigvals_re(j_outer * chunk_size + j_inner),
        //             y(i, j_outer, j_inner) >= _x(i) && y_flipped >= _x(i)};
    } else {
        throw std::runtime_error{"invalid spin_inversion"};
    }

    // Horizontal reduction
    RDom m_lane{0, chunk_size, "m_lane"};
    Func lane_reduction{"lane_reduction"};
    lane_reduction(i, j_outer) = Tuple{cast<double>(0), cast<bool>(true)};
    lane_reduction(i, j_outer) =
        Tuple{lane_reduction(i, j_outer)[0] + temp(i, j_outer, m_lane)[0],
              lane_reduction(i, j_outer)[1] && temp(i, j_outer, m_lane)[1]};

    // Vertical reduction
    RDom m_batch{0, number_chunks, "m_batch"};
    Func batch_reduction{"batch_reduction"};
    batch_reduction(i) = Tuple{cast<double>(0), cast<bool>(true)};
    m_batch.where(batch_reduction(i)[1]);
    batch_reduction(i) = Tuple{batch_reduction(i)[0] + lane_reduction(i, m_batch)[0],
                               batch_reduction(i)[1] && lane_reduction(i, m_batch)[1]};

    // Compute reduction over remaining elements
    Func scalar_reduction{"scalar_reduction"};
    RDom m_tail{0, number_rest, "m_tail"};
    scalar_reduction(i) = Tuple{batch_reduction(i)[0], batch_reduction(i)[1]};
    // m_tail.where(batch_reduction(i)[1]);
    scalar_reduction(i) = Tuple{scalar_reduction(i)[0] + temp(i, number_chunks, m_tail)[0],
                                scalar_reduction(i)[1] && temp(i, number_chunks, m_tail)[1]};

    // Save results
    auto const cutoff = Expr{1024 * std::numeric_limits<double>::epsilon()};
    auto const n = scalar_reduction(i)[0];
    auto const rounded = select(n < cutoff, cast<double>(0), n);

    Func norm{"norm"};
    norm(i) = select(n > cutoff && scalar_reduction(i)[1], n, 0);

    // Shapes & strides
    _x.dim(0).set_min(0).set_stride(1);

    // Schedule
    temp.in(lane_reduction).vectorize(j_inner);
    temp.in(lane_reduction).compute_at(lane_reduction, j_outer);
    temp.in(scalar_reduction).compute_at(scalar_reduction, i);
    // lane_reduction.compute_at(batch_reduction, m_batch);
    scalar_reduction.compute_at(norm, i);

    // norm.print_loop_nest();
    return norm.compile_to_callable({_x});
}

auto invoke_is_representative_kernel(void *data, va_alist alist) noexcept -> void {
    auto *ctx = static_cast<std::tuple<Halide::Callable, void *> *>(data);
    auto const &fn = std::get<0>(*ctx);
    va_start_void(alist);
    auto *basis_states = va_arg_ptr(alist, halide_buffer_t *);
    auto *norms = va_arg_ptr(alist, halide_buffer_t *);
    try {
        fn(basis_states, norms);
    } catch (InternalError &e) {
        fprintf(stderr, "Halide::InternalError: %s\n", e.what());
        abort();
    }
    va_return_void(alist);
}

auto mk_state_info_kernel(halide_buffer_t const *_masks, halide_buffer_t const *_shifts,
                          int spin_inversion) -> Halide::Callable {
    auto masks = Halide::Buffer{Halide::Runtime::Buffer<uint64_t, 2>{*_masks}.copy(), "masks"};
    auto shifts = Halide::Buffer{Halide::Runtime::Buffer<uint64_t, 1>{*_shifts}.copy(), "shifts"};
    // auto const flip_mask = 0;

    ImageParam _x{halide_type_of<uint64_t>(), 1, "basis_states"};
    // ImageParam _repr{halide_type_of<uint64_t>(), 1, "representatives"};
    // ImageParam _indices{halide_type_of<int32_t>(), 1, "indices"};

    auto depth = masks.dim(0).extent();
    auto number_masks = masks.dim(1).extent();
    auto const chunk_size = 2 * 4;
    auto number_chunks = number_masks / chunk_size;
    auto number_rest = number_masks - number_chunks * chunk_size;

    auto const reduction_step_impl = [](Expr /*i*/, Tuple current, Expr y, Expr j) -> Tuple {
        auto current_r = current[0];
        auto current_j = current[1];

        auto is_less = y < current_r;
        auto next_r = select(is_less, y, current_r);
        auto next_j = select(is_less, j, current_j);

        return {next_r, next_j};
    };

    auto const reduction_step = [spin_inversion, reduction_step_impl](Expr i, Tuple current, Expr y,
                                                                      Expr j) -> Tuple {
        if (spin_inversion == 0) {
            return reduction_step_impl(i, current, y, j);
        }
        // if (spin_inversion == 1) {
        //     auto temp = reduction_step_impl(i, current, y, j);
        //     return reduction_step_impl(i, temp, y ^ flip_mask, c_re, c_im);
        // }
        // if (spin_inversion == -1) {
        //     auto temp = reduction_step_impl(i, current, y, c_re, c_im);
        //     return reduction_step_impl(i, temp, y ^ flip_mask, -c_re, -c_im);
        // }
        throw std::runtime_error{"invalid spin_inversion"};
    };

    auto const combine_step = [](Tuple current, Tuple other) -> Tuple {
        auto current_r = current[0];
        auto current_j = current[1];
        auto r = other[0];
        auto j = other[1];

        auto is_less = r < current_r;
        auto next_r = select(is_less, r, current_r);
        auto next_j = select(is_less, j, current_j);

        return {next_r, next_j};
    };

    Var i{"i"};

    // Apply masks to generate transformed spin configurations
    RDom k{0, depth, "k"};
    Func y_batched{"y_batched"};
    Var j_outer{"j_outer"};
    Var j_inner{"j_inner"};
    y_batched(i, j_outer, j_inner) = _x(i);
    y_batched(i, j_outer, j_inner) = bit_permute_step_64(
        y_batched(i, j_outer, j_inner), masks(k, j_outer * chunk_size + j_inner), shifts(k));

    Func y_scalar{"y_scalar"};
    Var j_tail{"j_tail"};
    y_scalar(i, j_tail) = _x(i);
    y_scalar(i, j_tail) = bit_permute_step_64(y_scalar(i, j_tail), masks(k, j_tail), shifts(k));

    // Compute vectorized reduction
    RDom m_main{0, number_chunks, "m_main"};
    Func custom_reduction{"custom_reduction"};
    custom_reduction(i, j_inner) = Tuple{_x(i), cast<int32_t>(-1)};
    custom_reduction(i, j_inner) =
        reduction_step(i, custom_reduction(i, j_inner), y_batched(i, m_main, j_inner),
                       m_main * chunk_size + j_inner);

    // Reduce across one chunk
    RDom m_lane{1, chunk_size - 1, "m_lane"};
    Func reduction_lane{"reduction_lane"};
    reduction_lane(i) = Tuple{custom_reduction(i, 0)[0], custom_reduction(i, 0)[1]};
    reduction_lane(i) = combine_step(reduction_lane(i), custom_reduction(i, m_lane));

    // Compute reduction over remaining elements
    Func reduction_scalar{"reduction_scalar"};
    RDom m_tail{number_chunks * chunk_size, number_rest, "m_tail"};
    reduction_scalar(i) = Tuple{reduction_lane(i)[0], reduction_lane(i)[1]};
    reduction_scalar(i) = reduction_step(i, reduction_scalar(i), y_scalar(i, m_tail), m_tail);

    // Schedule
    y_batched.compute_at(custom_reduction, m_main);

    custom_reduction.compute_at(reduction_scalar, i);
    custom_reduction.vectorize(j_inner);
    custom_reduction.update(0).vectorize(j_inner);

    return reduction_scalar.compile_to_callable({_x});
}

auto invoke_state_info_kernel(void *data, va_alist alist) noexcept -> void {
    auto *ctx = static_cast<std::tuple<Halide::Callable, void *> *>(data);
    auto const &fn = std::get<0>(*ctx);
    va_start_void(alist);
    auto *basis_states = va_arg_ptr(alist, halide_buffer_t *);
    auto *representatives = va_arg_ptr(alist, halide_buffer_t *);
    auto *indices = va_arg_ptr(alist, halide_buffer_t *);
    try {
        fn(basis_states, representatives, indices);
    } catch (InternalError &e) {
        fprintf(stderr, "Halide::InternalError: %s\n", e.what());
        abort();
    }
    va_return_void(alist);
}

} // namespace lattice_symmetries

extern "C" ls_hs_state_to_index_kernel_type
mk_fixed_hamming_state_to_index_kernel(int const number_sites, int const hamming_weight,
                                       halide_buffer_t const *binomials) {
    try {
        // fprintf(stderr, "mk_fixed_hamming_state_to_index_kernel\n");
        auto fn = lattice_symmetries::mk_fixed_hamming_state_to_index_kernel(
            number_sites, hamming_weight, binomials);
        auto closure =
            alloc_callback(&lattice_symmetries::invoke_fixed_hamming_state_to_index_kernel,
                           new std::tuple<Halide::Callable, void *>{fn, nullptr});
        return reinterpret_cast<ls_hs_state_to_index_kernel_type>(closure);
    } catch (CompileError &e) {
        fprintf(stderr, "Halide::CompileError: %s\n", e.what());
        abort();
    } catch (InternalError &e) {
        fprintf(stderr, "Halide::InternalError: %s\n", e.what());
        abort();
    } catch (RuntimeError &e) {
        fprintf(stderr, "Halide::RuntimeError: %s\n", e.what());
        abort();
    }
}

extern "C" ls_hs_is_representative_kernel_type_v2
mk_is_representative_kernel(halide_buffer_t const *masks, halide_buffer_t const *eigvals_re,
                            halide_buffer_t const *shifts, int spin_inversion) {
    try {
        auto fn = lattice_symmetries::mk_is_representative_kernel(masks, eigvals_re, shifts,
                                                                  spin_inversion);
        auto closure = alloc_callback(&lattice_symmetries::invoke_is_representative_kernel,
                                      new std::tuple<Halide::Callable, void *>{fn, nullptr});
        return reinterpret_cast<ls_hs_is_representative_kernel_type_v2>(closure);
    } catch (CompileError &e) {
        fprintf(stderr, "Halide::CompileError: %s\n", e.what());
        abort();
    } catch (InternalError &e) {
        fprintf(stderr, "Halide::InternalError: %s\n", e.what());
        abort();
    } catch (RuntimeError &e) {
        fprintf(stderr, "Halide::RuntimeError: %s\n", e.what());
        abort();
    }
}

extern "C" ls_hs_state_info_kernel_type_v2 mk_state_info_kernel(halide_buffer_t const *masks,
                                                                halide_buffer_t const *shifts,
                                                                int spin_inversion) {
    try {
        auto fn = lattice_symmetries::mk_state_info_kernel(masks, shifts, spin_inversion);
        auto closure = alloc_callback(&lattice_symmetries::invoke_state_info_kernel,
                                      new std::tuple<Halide::Callable, void *>{fn, nullptr});
        return reinterpret_cast<ls_hs_state_info_kernel_type_v2>(closure);
    } catch (CompileError &e) {
        fprintf(stderr, "Halide::CompileError: %s\n", e.what());
        abort();
    } catch (InternalError &e) {
        fprintf(stderr, "Halide::InternalError: %s\n", e.what());
        abort();
    } catch (RuntimeError &e) {
        fprintf(stderr, "Halide::RuntimeError: %s\n", e.what());
        abort();
    }
}
