#include "intrinsics.h"
#include "declarations.h"

typedef struct Vix2 { Vi a; Vi b; } Vix2;
typedef struct Vdx2 { Vd a; Vd b; } Vdx2;

static void Pz(Vz z) { Pd(z.re); Pd(z.im); }
#define unpack_pair(ret0, ret1, ...) do { typeof(__VA_ARGS__) const _t = (__VA_ARGS__); ret0 = _t.a; ret1 = _t.b; } while (0)

// INTERNAL Vd pq2pd(Vi const x) { x = or(x, d2i(Sd(0x0010000000000000))); return subd(i2d(x), Sd(0x0010000000000000)); }
// INTERNAL Vd mask_gather_norm(u16 const *p, Vi const i, Vi const m) { Vi v = OP(mask_i64gather_epi64, Zi, (i64 const*)p, i, m, 2); return pq2pd(AND(v, Si(0xFFFF))); }
// INTERNAL Vd gather_norm(u16 const *p, Vi const i) { return pq2pd(AND(OP(i64gather_epi64, (i64 const*)p, i, 2), Si(0xFFFF))); }

// INTERNAL Vi m1(Vi const x, Vi const m) { return NOT(EQ(AND(x, m), Zi)); }
INTERNAL Vi m2(Vi const x, Vi const m_0, Vi const m_1) { return xor(m1(x, m_0), m1(x, m_1)); }
INTERNAL Vi mX(Vi const x, Vi const m) { return popcnt(and(x, m)); }
INTERNAL Vd signedd(Vd const v, Vi m) { m = shl(m, 63); return i2d(xor(d2i(v), m)); }
INTERNAL Vz signedz(Vz const v, Vi m) { m = shl(m, 63); return (Vz){i2d(xor(d2i(v.re), m)), i2d(xor(d2i(v.im), m))}; }

INTERNAL Vd zerod(void) { return Zd; }
INTERNAL Vz zeroz(void) { return (Vz){Zd, Zd}; }
INTERNAL Vd broadcast2d(f64 const *re, f64 const *im, i32 const k) { return Sd(re[k]); }
INTERNAL Vz broadcast2z(f64 const *re, f64 const *im, i32 const k) { return (Vz){Sd(re[k]), Sd(im[k])}; }
// INTERNAL Vd addd(Vd const a, Vd const b) { return ADD(a, b); }
// INTERNAL Vz addz(Vz const a, Vz const b) { return (Vz){ADD(a.re, b.re), ADD(a.im, b.im)}; }
// INTERNAL Vd muld(Vd const a, Vd const b) { return OP(mul_pd, a, b); }
// INTERNAL Vz mulz(Vz const a, Vz const b) { return (Vz){OP(fnmadd_pd, a.im, b.im, muld(a.re, b.re)), OP(fmadd_pd, a.im, b.re, mulpd(a.re, b.im))}; }
// INTERNAL Vd scaled(Vd const a, Vd const b) { return muld(a, b); }
// INTERNAL Vz scalez(Vd const a, Vz const b) { return (Vz){scaled(a, b.re), scaled(a, b.im)}; }
// INTERNAL Vd loadd(f64 const *p) { return Rd(p); }
// INTERNAL Vz loadz(c128 const *p) {
//     Vd const a = loadd((f64 const*)p), b = loadd((f64 const*)p + N);
//     Vd const c = OP(unpacklo_pd, a, b), d = OP(unpackhi_pd, a, b);
//     return (Vz){OP(permute4x64_pd, c, 0b11011000), OP(permute4x64_pd, d, 0b11011000)};
// }
// INTERNAL void stored(f64 *p, Vd z) { Wd(p, z); }
// INTERNAL void storez(c128 *p, Vz z) {
//     Vd const a = OP(permute4x64_pd, z.re, 0b11011000), b = OP(permute4x64_pd, z.im, 0b11011000);
//     Vd const c = OP(unpacklo_pd, a, b), d = OP(unpackhi_pd, a, b);
//     stored((f64*)p, c);
//     stored((f64*)p + N, d);
// }
// INTERNAL Vd gatherd(f64 const *p, Vi const i) { return OP(i64gather_pd, p, i, 8); }
// INTERNAL Vz gatherz(c128 const *p, Vi const i) { Vi const i2 = SHL(i, 1); return (Vz){gatherd((f64 const*)p, i2), gatherd((f64 const*)p + 1, i2)}; }
// INTERNAL Vd mask_gatherd(f64 const *p, Vi const i, Vi const m) { return OP(mask_i64gather_pd, zerod(), p, i, i2d(m), 8); }
// INTERNAL Vz mask_gatherz(c128 const *p, Vi const i, Vi const m) { Vi const i2 = SHL(i, 1); return (Vz){mask_gatherd((f64 const*)p, i2, m), mask_gatherd((f64 const*)p + 1, i2, m)}; }
// INTERNAL Vd gather2d(f64 const *re, f64 const *im, Vi const i) { return gatherd(re, i); }
// INTERNAL Vz gather2z(f64 const *re, f64 const *im, Vi const i) { return (Vz){gatherd(re, i), gatherd(im, i)}; }

#define coeff_template(s) \
    INTERNAL V##s coeff##s(Vi x, i32 const n_s0, i32 const n_s1, i32 const n_s2, i32 const n_sX, \
            f64 const *v_re, f64 const *v_im, u64 const *s1, u64 const *s20, u64 const *s21, u64 const *sX) { \
        i32 r = 0; Vi m; V##s acc; \
        if (n_s0 > 0) { acc = broadcast2##s(v_re, v_im, r); ++r; } else { acc = zero##s(); } \
        for (i32 k = 0; k < n_s1; ++k, ++r) { m = m1(x, Si(s1[k])); acc = add##s(acc, signed##s(broadcast2##s(v_re, v_im, r), m)); } \
        for (i32 k = 0; k < n_s2; ++k, ++r) { m = m2(x, Si(s20[k]), Si(s21[k])); acc = add##s(acc, signed##s(broadcast2##s(v_re, v_im, r), m)); } \
        for (i32 k = 0; k < n_sX; ++k, ++r) { m = mX(x, Si(sX[k])); acc = add##s(acc, signed##s(broadcast2##s(v_re, v_im, r), m)); } \
        return acc; \
    }
coeff_template(d)
coeff_template(z)
#undef coeff_template

#define diagN_template(s, t) \
    INTERNAL void diag##s##N(i64 const bi, u64 const *alpha0, void const *x0, void *out, oc_t const *ctx) { \
        Vi const alpha = Ri(alpha0 + bi); V##s const x = R##s((t const*)x0 + bi); \
        V##s const acc = ctx->n_t == 0 ? zero##s() : coeff##s(alpha, ctx->n_s0[0], ctx->n_s1[0], \
            ctx->n_s2[0], ctx->n_sX[0], ctx->v_re, ctx->v_im, ctx->s1, ctx->s20, ctx->s21, ctx->sX); \
        V##s const y = mul##s(acc, x); \
        W##s(out, y); \
    }
diagN_template(d, f64)
diagN_template(z, c128)
#undef diagN_template

#define diag64_template(s, t) \
    void diag64_##t(i64 const bi, u64 const *alpha0, void const *x0, void *out, const oc_t *ctx) { \
        for (i32 k = 0; k < 64; k += N) { diag##s##N(bi + k, alpha0, x0, (t*)out + k, ctx); } \
    }
diag64_template(d, f64)
diag64_template(z, c128)
#undef diag64_template

INTERNAL Vi pstep(Vi const x, Vi const m, u32 const d) { assume(0 <= d && d < 64); Vi const y = and(xor(shr(x, d), x), m); return xor(xor(x, y), shl(y, d)); }
// Assumes that n > 0
INTERNAL Vi permute(Vi x, u64 const *masks, u32 const *shifts, i32 const n) { i32 i = 0; do { x = pstep(x, Si(masks[i]), shifts[i]); ++i; } while (i < n); return x; }
INTERNAL Vix2 reprN(Vi x, bs_ctx_t const *ctx) {
    // the first row of ctx->masks is always the identity permutation
    i64 k = 1; u64 const *masks = ctx->masks + k * ctx->n_r; Vi r = x, i = Zi;
    for (Vi vk = Si(k), one = Si(1); k < ctx->n_m; ++k, vk = addq(vk, one), masks += ctx->n_r) {
        Vi const y = permute(x, masks, ctx->shifts, ctx->n_r), p = gt(r, y);
        r = select(p, y, r);  i = select(p, vk, i);
    }
    return (Vix2){r, i};
}
INTERNAL Vix2 searchN(Vi needle, search_ctx_t const *ctx) {
    i64 n = ctx->range_size; Vi base = gatherq(ctx->offsets, and(shr(needle, ctx->shift), Si(ctx->mask))); Vi vsz = Si(n);
    while (n > 1) {
        Vi const h = shr(vsz, 1), k = addq(base, h), y = gatherq(ctx->reps, k);
        n -= n / 2; vsz = subq(vsz, h);
        base = select(gt(needle, y), k, base);
    }
    Vi y = gatherq(ctx->reps, base);
    base = select(gt(needle, y), addq(base, Si(1)), base);
    y = gatherq(ctx->reps, base);
    return (Vix2){eqi(needle, y), base};
}

#define off_diagN_template(s, t) \
    INTERNAL void off_diag##s##N(i64 const bi, u64 const *alpha0, u16 const *norm0, void const *X, void *out,\
            oc_t const *ctx, const bs_ctx_t *bs_ctx, const search_ctx_t *search_ctx) { \
        Vi const alpha = Ri(alpha0 + bi); Vd n0; V##s acc_outer = zero##s(); \
        if (bs_ctx != NULL) { n0 = load_norm(norm0 + bi); } \
        f64 const *v_re = ctx->v_re, *v_im = ctx->v_im; u64 const *s1 = ctx->s1, *s20 = ctx->s20, *s21 = ctx->s21, *sX = ctx->sX; \
        for (i32 ti = 0; ti < ctx->n_t; ++ti) { \
            V##s x, chi; Vd norm; Vi const beta0 = xor(alpha, Si(ctx->mask[ti])); \
            if (bs_ctx != NULL) { \
                Vi rep, gid; unpack_pair(rep, gid, reprN(beta0, bs_ctx)); \
                Vi msk, idx; unpack_pair(msk, idx, searchN(rep, search_ctx)); \
                norm = mask_gather_norm(search_ctx->norm, idx, msk); \
                x = mask_gather##s((t const*)X, idx, msk); \
                chi = gather2##s(bs_ctx->chi_re, bs_ctx->chi_im, gid); \
            } \
            else { \
                x = gather##s((t const*)X, beta0); \
            } \
            V##s const acc = coeff##s(alpha, ctx->n_s0[ti], ctx->n_s1[ti], ctx->n_s2[ti], ctx->n_sX[ti], \
                v_re, v_im, s1, s20, s21, sX); \
            if (bs_ctx != NULL) { \
                Vd const c = sqrtd(divd(norm, n0)); \
                x = scale##s(c, mul##s(chi, x)); \
            } \
            acc_outer = add##s(acc_outer, mul##s(acc, x)); \
            s1 += ctx->n_s1[ti]; s20 += ctx->n_s2[ti]; s21 += ctx->n_s2[ti]; sX += ctx->n_sX[ti]; \
            i32 const k = ctx->n_s0[ti] + ctx->n_s1[ti] + ctx->n_s2[ti] + ctx->n_sX[ti]; \
            v_re += k; v_im += k; \
        } \
        W##s(out, add##s(R##s(out), acc_outer)); \
    }
off_diagN_template(d, f64)
off_diagN_template(z, c128)
#undef off_diagN_template

#define off_diag64_template(s, t) \
    void off_diag64_##t(i64 const bi, u64 const *alpha0, u16 const *norm0, void const *X, void *out, \
            oc_t const *ctx, bs_ctx_t const *bs_ctx, search_ctx_t const* search_ctx) { \
        for (i32 k = 0; k < 64; k += N) { off_diag##s##N(bi + k, alpha0, norm0, X, (t*)out + k, ctx, bs_ctx, search_ctx); } \
    }
off_diag64_template(d, f64)
off_diag64_template(z, c128)
#undef off_diag64_template

#define matvec_inner_template(t) \
    void matvec_inner_##t(i64 const i, u64 const *alpha0, u16 const *norm0, void const *X0, void const *X, void *out, \
            oc_t const *diag, oc_t const *off_diag, bs_ctx_t const *bs, search_ctx_t const *search) { \
        diag64_##t(i, alpha0, X0, (t*)out + i, diag); off_diag64_##t(i, alpha0, norm0, X, (t*)out + i, off_diag, bs, search); \
    }
matvec_inner_template(f64)
matvec_inner_template(c128)
#undef matvec_inner_template

typedef void (*matvec_internal_t)(i64, u64 const *, u16 const *, void const *, void const *, void *,
                                  oc_t const *, oc_t const *, bs_ctx_t const *, search_ctx_t const *);

INTERNAL void matvec(i64 const n, u64 const *alpha0, u16 const *norm0, void const *X0, void const *X, void *out,
        oc_t const *diag_ctx, oc_t const *off_diag_ctx, bs_ctx_t const *bs_ctx, search_ctx_t const *search_ctx, matvec_internal_t inner) {
    if (n < 64) { return; }
    
    i64 const n_b = n / 64, n_r = n % 64;
#pragma omp parallel for schedule(dynamic, 256) default(none) \
        firstprivate(n_b, alpha0, norm0, X0, X, out, diag_ctx, off_diag_ctx, bs_ctx, search_ctx, inner)
    for (i64 bi = 0; bi < n_b; ++bi) { inner(64 * bi, alpha0, norm0, X0, X, out, diag_ctx, off_diag_ctx, bs_ctx, search_ctx); }
    if (n_r != 0) { inner(n - 64, alpha0, norm0, X0, X, out, diag_ctx, off_diag_ctx, bs_ctx, search_ctx); }
}

void matvec_f64(i64 const n, u64 const *alpha0, u16 const *norm0, void const *X0, void const *X, void *out,
        oc_t const *diag_ctx, oc_t const *off_diag_ctx, bs_ctx_t const *bs_ctx, search_ctx_t const *search_ctx) {
    matvec(n, alpha0, norm0, X0, X, out, diag_ctx, off_diag_ctx, bs_ctx, search_ctx, matvec_inner_f64);
}
void matvec_c128(i64 const n, u64 const *alpha0, u16 const *norm0, void const *X0, void const *X, void *out,
        oc_t const *diag_ctx, oc_t const *off_diag_ctx, bs_ctx_t const *bs_ctx, search_ctx_t const *search_ctx) {
    matvec(n, alpha0, norm0, X0, X, out, diag_ctx, off_diag_ctx, bs_ctx, search_ctx, matvec_inner_c128);
}

// #define INNER(i, diag, off_diag, bs, search, out) \
//     do { diag64_f64(i, diag, out + (i)); off_diag64_f64(i, off_diag, bs, search, out + (i)); } while (0)
// void matvec_f64(int64_t const n, const oc_t *diag_ctx, const oc_t *off_diag_ctx,
//         const bs_ctx_t *bs_ctx, const search_ctx_t *search_ctx, void *out) {
//     if (n < 64) { return; }

//     int64_t const n_b = n / 64, n_r = n % 64;
// #pragma omp parallel for schedule(dynamic, 256)
//     for (int64_t bi = 0; bi < n_b; ++bi) { INNER(64 * bi, diag_ctx, off_diag_ctx, bs_ctx, search_ctx, out); }
//     if (n_r != 0) { INNER(n - 64, diag_ctx, off_diag_ctx, bs_ctx, search_ctx, out); }
// }
// #undef INNER

// #define SIGNED(v, m) i2d(XOR(d2i(v), SHL(m, 63)))
// #define M1(x, m) NOT(EQ(AND(x, m), Zi))
// #define M2(x, m0, m1) XOR(M1(x, m0), M1(x, m1))
// #define MX(x, m) POPCNT(AND(x, m))
// #define U(acc, v, m) acc = ADD(acc, SIGNED(Sd(v[r]), m))
// #define gcoeff(x, C0T, C0F, C1) \
//     do { \
//         i32 r = 0; Vi m; \
//         if (n_s0[ti] > 0) { C0T; ++r; } else { C0F; } \
//         for (i32 k = 0; k < n_s1[ti]; ++k, ++r) { m = M1(x, Si(s1[k])); C1; } \
//         for (i32 k = 0; k < n_s2[ti]; ++k, ++r) { m = M2(x, Si(s20[k]), Si(s21[k])); C1; } \
//         for (i32 k = 0; k < n_sX[ti]; ++k, ++r) { m = MX(x, Si(sX[k])); C1; } \
//     } while (0)
// #define unpack_oc_ctx_t(ctx) \
//     f64 const *v_re = ctx->v_re, *v_im = ctx->v_im; \
//     i32 const *n_s0 = ctx->n_s0, *n_s1 = ctx->n_s1, *n_s2 = ctx->n_s2, *n_sX = ctx->n_sX; \
//     u64 const *s1 = ctx->s1, *s20 = ctx->s20, *s21 = ctx->s21, *sX = ctx->sX
    
// static HEDLEY_ALWAYS_INLINE void diagN_f64(i64 const bi, const oc_t *ctx, double *out) {
//     unpack_oc_ctx_t(ctx); Vi const alpha = Ri(ctx->alpha0 + bi); i32 const ti = 0; Vd acc_re;
//     gcoeff(alpha, acc_re = Sd(v_re[r]), acc_re = Zd, U(acc_re, v_re, m));
//     Wd(out, MUL(acc_re, Rd(ctx->X0 + bi)));
// }
// static HEDLEY_ALWAYS_INLINE void diagN_c128(i64 const bi, const oc_t *ctx, _Complex double *out) {
//     unpack_oc_ctx_t(ctx); Vi const alpha = Ri(ctx->alpha0 + bi); i32 const ti = 0; Vd acc_re, acc_im;
//     gcoeff(alpha, (acc_re = Sd(ctx->v_re[r]), acc_im = Sd(ctx->v_im[r])),
//            (acc_re = Zd, acc_im = Zd), (U(acc_re, v_re, m), U(acc_im, v_im, m)));
//     Vdx2 const _t = cplx_mul(acc_re, acc_im, (_Complex double const*)ctx->X0 + bi);
//     Wd((double*)out, OP(unpacklo_pd, _t.a, _t.b)), Wd((double*)out + N, OP(unpackhi_pd, _t.a, _t.b));
// }

// #define diag64(t, s) \
//     void diag64_##s(i64 const bi, const oc_t *ctx, void *out) { \
//         for (i32 k = 0; k < 64; k += N) { diagN_##s(bi + k, ctx, (t*)out + k); } \
//     }
// diag64(double, f64)
// diag64(_Complex double, c128)
// #undef diag64

// #define UPDATE(v) acc = ADD(acc, SIGNED(Sd((v)[r]), m))
// #define COEFF(v_re, n_s0, n_s1, s1, n_s2, s20, s21, n_sX, sX) \
//     do { int32_t r = 0; Vi m; \
//     if ((n_s0) > 0) { acc = Sd((v_re)[r]); ++r; } else { acc = Zd; } \
//     for (int32_t k = 0; k < (n_s1); ++k, ++r) { m = M1(x, Si((s1)[k])); UPDATE(v_re); } \
//     for (int32_t k = 0; k < (n_s2); ++k, ++r) { m = M2(x, Si((s20)[k]), Si((s21)[k])); UPDATE(v_re); } \
//     for (int32_t k = 0; k < (n_sX); ++k, ++r) { m = MX(x, Si((sX)[k])); UPDATE(v_re); } \
//     } while(0)



#define with_in(tmp_x, x, t, n, ...) \
    t tmp_x[N]; __builtin_memset(tmp_x, 0, N * sizeof(t)); __builtin_memcpy(tmp_x, x, n * sizeof(t)); __VA_ARGS__
#define with_out(tmp_x, x, t, n, ...) \
    t tmp_x[N]; __builtin_memset(tmp_x, 0, N * sizeof(t)); __VA_ARGS__; __builtin_memcpy(x, tmp_x, n * sizeof(t))

#define INNER(i) \
    do { Vix2 const _t = reprN(Ri(xs + (i)), ctx); Wi(rep + (i), _t.a); Wi(idx + (i), _t.b); } while (0)
void state_info(int64_t const n, uint64_t const *xs, bs_ctx_t const *ctx, uint64_t *rep, int64_t *idx) {
    if (n <= 0) { return; }
    if (n < N) {
        with_in(tmp_xs, xs, uint64_t, n, with_out(tmp_rep, rep, uint64_t, n, with_out(tmp_idx, idx, int64_t, n,
            state_info(N, tmp_xs, ctx, tmp_rep, tmp_idx))));
        return;
    }

    int64_t const n_b = n / N, n_r = n % N;
#pragma omp parallel for default(none) firstprivate(n_b, xs, ctx, rep, idx)
    for (int64_t bi = 0; bi < n_b; ++bi) { INNER(N * bi); }
    if (n_r != 0) { INNER(n - N); }
}
#undef INNER

#define INNER(i) \
    do { Vix2 _t = searchN(Ri(xs + (i)), ctx); Wi(out + (i), select(_t.a, _t.b, Si(-1))); } while (0)
void state_to_index(int64_t const n, uint64_t const *xs, search_ctx_t const *ctx, int64_t *out) {
    if (n <= 0) { return; }
    if (n < N) {
        with_in(tmp_xs, xs, uint64_t, n, with_out(tmp_out, out, int64_t, n,
            state_to_index(N, tmp_xs, ctx, tmp_out)));
        return;
    }

    int64_t const n_b = n / N, n_r = n % N;
#pragma omp parallel for default(none) firstprivate(n_b, xs, ctx, out)
    for (int64_t bi = 0; bi < n_b; ++bi) { INNER(N * bi); }
    if (n_r != 0) { INNER(n - N); }
}
#undef INNER




// #define gather_c128(a, b, base, index) \
//     do { \
//         Vi const i2 = SHL(index, 1); \
//         a = OP(i64gather_pd, Zd, (f64 const*)base, i2, 8); \
//         b = OP(i64gather_pd, Zd, (f64 const*)base + 1, i2, 8); \
//     } while (0)
// #define mask_gather_c128(a, b, base, index, mask) \
//     do { \
//         Vi const i2 = SHL(index, 1); \
//         a = OP(mask_i64gather_pd, Zd, (f64 const*)base, i2, i2d(mask), 8); \
//         b = OP(mask_i64gather_pd, Zd, (f64 const*)base + 1, i2, i2d(mask), 8); \
//     } while (0)

// #define gather_chi chi_re = gatherqpd(bs_ctx->chi_re, gid)
// #define gather_xT x_re = OP(mask_i64gather_pd, Zd, (f64 const*)ctx->X, idx, i2d(msk), 8)
// #define gather_xF x_re = gatherqpd(ctx->X, beta0)
// #define scale_x x_re = mulpd(c, mulpd(chi_re, x_re))

// #define gather_chi chi_re = gatherqpd(bs_ctx->chi_re, gid), chi_im = gatherqpd(bs_ctx->chi_im, gid)
// #define gather_xT mask_gather_c128(x_re, x_im, ctx->X, idx, msk)
// #define gather_xF gather_c128(x_re, x_im, ctx->X, idx)
// #define scale_x do { unpack_pair(x_re, x_im, mulpz(chi_re, chi_im, x_re, x_im)); x_re = mulpd(c, x_re); x_im = mulpd(c, x_im); } while (0)

// static HEDLEY_ALWAYS_INLINE void off_diag_one(int64_t const bi, const oc_t *ctx,
//         const bs_ctx_t *bs_ctx, const search_ctx_t *search_ctx, double *out) {
//     Vi const alpha = Ri(ctx->alpha0 + bi); Vd acc_outer = Zd, norm0;
//     if (bs_ctx != NULL) { norm0 = gatherqw(ctx->norm0 + bi, OP(set_epi64x, 3, 2, 1, 0)); }
//     unpack_oc_ctx_t(ctx);
//     for (int32_t ti = 0; ti < ctx->n_t; ++ti) {
//         Vd acc_re, acc_im, x_re, x_im, chi_re, chi_im, norm;
//         Vi beta0 = XOR(alpha, Si(ctx->mask[ti]));
//         if (bs_ctx != NULL) {
//             Vi rep, gid; unpack_repr(rep, gid, repr(beta0, bs_ctx));
//             Vi msk, idx; unpack_search(msk, idx, search(rep, search_ctx));
//             norm = mask_gatherqw(ctx->norm, idx, msk); gather_xT; gather_chi;
//         }
//         else {
//             gather_xF;
//         }
//         gcoeff(alpha, acc_re = Sd(v_re[r]), acc_re = Zd, U(acc_re, v_re, m));
//         if (bs_ctx != NULL) { Vd const c = sqrtpd(divpd(norm, norm0)); scale_x; }

//         // COEFF(v_re, ctx->n_s0[ti], ctx->n_s1[ti], s1, ctx->n_s2[ti], s20, s21, ctx->n_sX[ti], sX);
//         // printf("acc="); Pd(acc);
//         acc_outer = FMA(acc_re, x_re, acc_outer);
//         s1 += ctx->n_s1[ti]; s20 += ctx->n_s2[ti]; s21 += ctx->n_s2[ti]; sX += ctx->n_sX[ti];
//         v_re += ctx->n_s0[ti] + ctx->n_s1[ti] + ctx->n_s2[ti] + ctx->n_sX[ti];
//     }
//     Wd(out, ADD(Rd(out), acc_outer));
// }





// HEDLEY_NEVER_INLINE void off_diag64(int64_t const bi, const oc_t *ctx,
//         const bs_ctx_t *bs_ctx, const search_ctx_t *search_ctx, double *out) {
//     for (int k = 0; k < 64; k += N) { off_diag_one(bi + k, ctx, bs_ctx, search_ctx, out + k); }
// }


// #define AP_F(v) flag |= OP(movemask_pd, i2d(GT(x, v)))
// #define AP_E(v) norm = addq(norm, SHR(EQ(x, v), 63))
INTERNAL u32 ap_f(Vi const x, Vi const v) { return I(B,movemask_pd,i2d(gt(x, v))); }
INTERNAL Vi ap_e(Vi const norm, Vi const x, Vi const v) { return addq(norm, shr(eqi(x, v), 63)); }
static HEDLEY_ALWAYS_INLINE Vi norm(Vi x, bs_ctx_t const* ctx) {
    int k = 1; unsigned flag = 0; Vi norm = Si(1), m = Si(ctx->inversion_mask);
    uint8_t const* flags = ctx->flags + k * 3; uint64_t const* masks = ctx->masks + k * ctx->n_r;
    for (; k < ctx->n_m; ++k, masks += ctx->n_r, flags += 3) {
        uint8_t const use_f2 = flags[0], use_e1 = flags[1], use_e2 = flags[2];
        Vi y = permute(x, masks, ctx->shifts, ctx->n_r), y2;
        flag |= ap_f(x, y);
        // if (use_f2) { y2 = XOR(y, m); AP_F(y2); }
        if (use_e1) { norm = ap_e(norm, x, y); }
        else { flag |= I(B,movemask_pd,i2d(eqi(x, y))); }
        if (flag == 0b1111) { return Zi; }
        // if (use_e2) { AP_E(y2); }
    }
    Vi const c = I(B,set_epi64x,0b1000,0b100,0b10,0b1);
    Vi const p = eqi(and(Si(flag), c), c);
    return andnot(p, norm);
}
// #undef AP_F
// #undef AP_E

#define REORDER(x) simde_mm256_castpd_si256(simde_mm256_permute4x64_pd(x, _MM_SHUFFLE(3, 1, 2, 0)))
static HEDLEY_ALWAYS_INLINE simde__m256i pack64to32(simde__m256i a, simde__m256i b)
{
    simde__m256 const combined = simde_mm256_shuffle_ps(
        simde_mm256_castsi256_ps(a), simde_mm256_castsi256_ps(b), _MM_SHUFFLE(2, 0, 2, 0));
    return REORDER(simde_mm256_castps_pd(combined));
}
static HEDLEY_ALWAYS_INLINE simde__m256i pack32to16(simde__m256i a, simde__m256i b)
{
    a = simde_mm256_and_si256(a, simde_mm256_set1_epi32(0xFFFF));
    b = simde_mm256_and_si256(b, simde_mm256_set1_epi32(0xFFFF));
    simde__m256i packed = simde_mm256_packus_epi32(a, b);
    return REORDER(simde_mm256_castsi256_pd(packed));
}
#undef REORDER
static Vi pack64to16(Vi a, Vi b, Vi c, Vi d) { return pack32to16(pack64to32(a, b), pack64to32(c, d)); }

void norm64(uint64_t const *alpha, bs_ctx_t const *ctx, uint16_t *out) {
    for (int k = 0; k < 64; k += 16) {
        Vi const n0 = norm(Ri(alpha + k + 0), ctx),
                 n1 = norm(Ri(alpha + k + 4), ctx),
                 n2 = norm(Ri(alpha + k + 8), ctx),
                 n3 = norm(Ri(alpha + k + 12), ctx);
        Wi((i64 *)(out + k), pack64to16(n0, n1, n2, n3));
    }
}



