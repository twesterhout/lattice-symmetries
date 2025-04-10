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
// INTERNAL Vi m2(Vi const x, Vi const m_0, Vi const m_1) { return xor(m1(x, m_0), m1(x, m_1)); }

Dd(zero,Zd,void);Dz(zero,(Z2(Zd,Zd)),void)
Dd(bcast2,Sd(re[k]),c(f64)*re,c(f64)*im,c(i32)k);Dz(bcast2,(Z2(Sd(re[k]),Sd(im[k]))),c(f64)*re,c(f64)*im,c(i32)k)

// Initializer for coeff##t
#define Ci(t,u,n) i32 r=0;Vi m[u];V##t acc[u];$(n>0,_(_L(_b,u,acc[_b]=bcast2##t(v_re,v_im,r));++r)){_L(_b,u,acc[_b]=zero##t())}
// Inner loop for coeff##t
#define Ck(t,u,b,x...) _L(k,b,_L(_b,u,m[_b]=(x));_L(_b,u,acc[_b]=add##t(acc[_b],signed##t(bcast2##t(v_re,v_im,r),m[_b])));++r)
static inline void coeffd1xN(c(Vi)x[],Vd o[],c(i32)n_s0,c(i32)n_s1,c(i32)n_s2,c(i32)n_sX,c(f64)*v_re,c(f64)*v_im,c(u64)*s1,c(u64)*s20,c(u64)*s21,c(u64)*sX){Ci(d,1,n_s0)Ck(d,1,n_s1,m1(x[_b],Si(s1[k])))Ck(d,1,n_s2,m2(x[_b],Si(s20[k]),Si(s21[k])))Ck(d,1,n_sX,mX(x[_b],Si(sX[k])))_L(_b,1,o[_b]=acc[_b])}
static inline void coeffz1xN(c(Vi)x[],Vz o[],c(i32)n_s0,c(i32)n_s1,c(i32)n_s2,c(i32)n_sX,c(f64)*v_re,c(f64)*v_im,c(u64)*s1,c(u64)*s20,c(u64)*s21,c(u64)*sX){Ci(z,1,n_s0)Ck(z,1,n_s1,m1(x[_b],Si(s1[k])))Ck(z,1,n_s2,m2(x[_b],Si(s20[k]),Si(s21[k])))Ck(z,1,n_sX,mX(x[_b],Si(sX[k])))_L(_b,1,o[_b]=acc[_b])}
#undef Ci
#undef Ck

static inline Vd coeffd(c(Vi)x,c(i32)n_s0,c(i32)n_s1,c(i32)n_s2,c(i32)n_sX,c(f64)*v_re,c(f64)*v_im,c(u64)*s1,c(u64)*s20,c(u64)*s21,c(u64)*sX){Vi xarr[1]={x}; Vd oarr[1];coeffd1xN(xarr,oarr,n_s0,n_s1,n_s2,n_sX,v_re,v_im,s1,s20,s21,sX);return oarr[0];}
static inline Vz coeffz(c(Vi)x,c(i32)n_s0,c(i32)n_s1,c(i32)n_s2,c(i32)n_sX,c(f64)*v_re,c(f64)*v_im,c(u64)*s1,c(u64)*s20,c(u64)*s21,c(u64)*sX){Vi xarr[1]={x}; Vz oarr[1];coeffz1xN(xarr,oarr,n_s0,n_s1,n_s2,n_sX,v_re,v_im,s1,s20,s21,sX);return oarr[0];}

// #define coeff_template(s) \
//     INTERNAL V##s coeff##s(Vi x, i32 const n_s0, i32 const n_s1, i32 const n_s2, i32 const n_sX, \
//             f64 const *v_re, f64 const *v_im, u64 const *s1, u64 const *s20, u64 const *s21, u64 const *sX) { \
//         i32 r = 0; Vi m; V##s acc; \
//         if (n_s0 > 0) { acc = bcast2##s(v_re, v_im, r); ++r; } else { acc = zero##s(); } \
//         for (i32 k = 0; k < n_s1; ++k, ++r) { m = m1(x, Si(s1[k])); acc = add##s(acc, signed##s(bcast2##s(v_re, v_im, r), m)); } \
//         for (i32 k = 0; k < n_s2; ++k, ++r) { m = m2(x, Si(s20[k]), Si(s21[k])); acc = add##s(acc, signed##s(bcast2##s(v_re, v_im, r), m)); } \
//         for (i32 k = 0; k < n_sX; ++k, ++r) { m = mX(x, Si(sX[k])); acc = add##s(acc, signed##s(bcast2##s(v_re, v_im, r), m)); } \
//         return acc; \
//     }
// coeff_template(d)
// coeff_template(z)
// #undef coeff_template

static inline void diagdN(c(i64)bi,c(u64)*alpha0,c(void)*x0,void*out,c(oc_t)*ctx){Vi alpha[1];_L(_b,1,alpha[_b]=Ri(alpha0+_b*N+bi));Vd acc[1];$(ctx->n_t>0,coeffd1xN(alpha,acc,ctx->n_s0[0],ctx->n_s1[0],ctx->n_s2[0],ctx->n_sX[0],ctx->v_re,ctx->v_im,ctx->s1,ctx->s20,ctx->s21,ctx->sX)){_L(_b,1,acc[_b]=zerod())};_L(_b,1,Wd((f64*)out+_b*N,muld(acc[_b],Rd((f64*)x0+_b*N+bi))))}
static inline void diagzN(c(i64)bi,c(u64)*alpha0,c(void)*x0,void*out,c(oc_t)*ctx){Vi alpha[1];_L(_b,1,alpha[_b]=Ri(alpha0+_b*N+bi));Vz acc[1];$(ctx->n_t>0,coeffz1xN(alpha,acc,ctx->n_s0[0],ctx->n_s1[0],ctx->n_s2[0],ctx->n_sX[0],ctx->v_re,ctx->v_im,ctx->s1,ctx->s20,ctx->s21,ctx->sX)){_L(_b,1,acc[_b]=zeroz())};_L(_b,1,Wz((c128*)out+_b*N,mulz(acc[_b],Rz((c128*)x0+_b*N+bi))))}

// #define diagN_template(s, t) \
//     INTERNAL void diag##s##N(i64 const bi, u64 const *alpha0, void const *x0, void *out, oc_t const *ctx) { \
//         Vi const alpha[1] = {Ri(alpha0 + bi)}; V##s const x = R##s((t const*)x0 + bi); \
//         V##s acc[1] = {zero##s()}; \
//         if (ctx->n_t > 0) { coeff##s##1xN(alpha, acc, ctx->n_s0[0], ctx->n_s1[0], \
//             ctx->n_s2[0], ctx->n_sX[0], ctx->v_re, ctx->v_im, ctx->s1, ctx->s20, ctx->s21, ctx->sX); } \
//         V##s const y = mul##s(acc[0], x); \
//         W##s(out, y); \
//     }
// diagN_template(d, f64)
// diagN_template(z, c128)
// #undef diagN_template

#define diag64_template(s, t) \
    void diag64_##t(i64 const bi, u64 const *alpha0, void const *x0, void *out, const oc_t *ctx) { \
        for (i32 k = 0; k < 64; k += N) { diag##s##N(bi + k, alpha0, x0, (t*)out + k, ctx); } \
    }
diag64_template(d, f64)
diag64_template(z, c128)
#undef diag64_template

INTERNAL Vi pstep(Vi const x, Vi const m, u32 const d) { assume(0 <= d && d < 64); Vi const y = and(xor(shr(x, d), x), m); return xor(xor(x, y), shl(y, d)); }

// [ 1  2  4  8 16  8  4  2  1]
HEDLEY_NEVER_INLINE Vi permute64(Vi x, u64 const *masks) {
    x = pstep(x, Si(masks[0]), 1);
    x = pstep(x, Si(masks[1]), 2);
    x = pstep(x, Si(masks[2]), 4);
    x = pstep(x, Si(masks[3]), 8);
    x = pstep(x, Si(masks[4]), 16);
    x = pstep(x, Si(masks[5]), 32);
    x = pstep(x, Si(masks[6]), 16);
    x = pstep(x, Si(masks[7]), 8);
    x = pstep(x, Si(masks[8]), 4);
    x = pstep(x, Si(masks[9]), 2);
    x = pstep(x, Si(masks[10]), 1);
    return x;
}

// Assumes that n > 0
INTERNAL Vi permute(Vi x, u64 const *masks, u32 const *shifts, i32 const n) { i32 i = 0; do { x = pstep(x, Si(masks[i]), shifts[i]); ++i; } while (i < n); return x; }

INTERNAL Vix2 reprN(Vi x, bs_ctx_t const *ctx) {
    // the first row of ctx->masks is always the identity permutation
    i64 k = 1; u64 const *masks = ctx->masks + k * ctx->n_r; Vi r = x, i = Zi;
    for (Vi vk = Si(k), one = Si(1); k < ctx->n_m; ++k, vk = addq(vk, one), masks += ctx->n_r) {
        Vi const y = permute64(x, masks); // permute(x, masks, ctx->shifts, ctx->n_r);
        typeof(gt(r, y)) const p = gt(r, y);
        r = select(p, y, r);  i = select(p, vk, i);
    }
    return (Vix2){r, i};
}

// HEDLEY_NEVER_INLINE Vix2 searchN(Vi _needle, search_ctx_t const *ctx) {
//     i64 needle[N]; Wi(needle, _needle);
//     i64 msk[N]; i64 idx[N];
//     for (i32 k = 0; k < N; ++k) {
//         i64 base = ctx->offsets[(needle[k] >> ctx->shift) & ctx->mask];
//         i64 n = ctx->range_size;
//         while (n > 1) {
//             i64 const h = n / 2;
//             base += (needle[k] > ctx->reps[base + h]) ? h : 0;
//             n -= h;
//         }
//         base += (needle[k] > ctx->reps[base]) ? 1 : 0;
//         msk[k] = (ctx->reps[base] == needle[k]) ? -1 : 0;
//         idx[k] = base;
//     }
//     return (Vix2){Ri(msk), Ri(idx)};
// }

HEDLEY_NEVER_INLINE Vix2 searchN(Vi needle, search_ctx_t const *ctx) {
    i64 n = ctx->range_size; Vi base = gatherq(ctx->offsets, and(shr(needle, ctx->shift), Si(ctx->mask))); Vi vsz = Si(n);
    while (n > 1) {
        Vi const h = shr(vsz, 1), k = addq(base, h), y = gatherq((i64 const*)ctx->reps, k);
        n -= n / 2; vsz = subq(vsz, h);
        base = select(gt(needle, y), k, base);
    }
    Vi y = gatherq((i64 const*)ctx->reps, base);
    base = select(gt(needle, y), addq(base, Si(1)), base);
    y = gatherq((i64 const*)ctx->reps, base);
    return (Vix2){eqi(needle, y), base};
}

INTERNAL f64 c2d(f64 a, f64 b) { return a; }
INTERNAL c128 c2z(f64 a, f64 b) { return __builtin_complex(a, b); }

// void repr_search(Vi a0, Vi needle0, bs_ctx_t const *bs_ctx, search_ctx_t const *search_ctx,
//         i64 *out_rep, i64 *out_gid, i64 *out_idx) {
//     Vi rep, gid; unpack_pair(rep, gid, reprN(a0, bs_ctx));
//     Wi(out_rep, rep); Wi(out_gid, gid);

//     Vi msk, idx; unpack_pair(msk, idx, searchN(needle0, search_ctx));
//     Wi(out_idx, select(msk, idx, Si(-1)));
// }

#define P(x,m,s) x=pstep(x,m,s)
#define L(i,s) m=Si(masks[i]); P(b0,m,s)
HEDLEY_NEVER_INLINE void repr_search(Vi a0, Vi needle0, bs_ctx_t const *bs_ctx, search_ctx_t const *search_ctx,
        u64 *out_rep, i64 *out_gid, i64 *out_idx) {

    // Setup for search
    i64 n = search_ctx->range_size;
    Vi base = and(shr(needle0, search_ctx->shift), Si(search_ctx->mask));
    base = gatherq(search_ctx->offsets, base);
    Vi values;
    // Setup for repr
    u64 const *masks = bs_ctx->masks + bs_ctx->n_r; Vi rep0=a0; Vi gid0=Zi;
    Vi b0;
    for (i64 k = 1; k < bs_ctx->n_m;) {
        if (n > 1) {
            _Alignas(32) i64 buf[N]; Wi(buf, base);
            for (i32 _k=0;_k<N;++_k) {
                simde_mm_prefetch(search_ctx->reps + buf[_k] + n/2, _MM_HINT_T0);
            }
        }
        Vi m; b0=a0; L(0,1); L(1,2); L(2,4); L(3,8); L(4,16); L(5,32); L(6,16); L(7,8); L(8,4); L(9,2); L(10,1);
        if (n > 1) {
            values = gatherq((i64 const*)search_ctx->reps, addq(base, Si(n / 2)));
        }
        M8 const p0=gt(rep0,b0);
        rep0=select(p0,b0,rep0);
        gid0=select(p0,Si(k),gid0);
        ++k; masks += bs_ctx->n_r;
        if (n > 1) {
            i64 const h = n / 2;
            base = select(gt(needle0, values), addq(base, Si(h)), base);
            n -= h;
        }
    }
    
    while (n > 1) {
        i64 const h = n / 2;
        values = gatherq((i64 const*)search_ctx->reps, addq(base, Si(h)));
        base = select(gt(needle0, values), addq(base, Si(h)), base);
        n -= h;
    }
    values = gatherq((i64 const*)search_ctx->reps, base);
    base = select(gt(needle0, values), addq(base, Si(1)), base);
    values = gatherq((i64 const*)search_ctx->reps, base);

    Wi(out_rep,rep0);Wi(out_gid,gid0);
    Wi(out_idx,select(eqq(needle0, values), base, Si(-1)));
}
#undef P
#undef L

INTERNAL void prefetchq(u64 const* p, Vi idx) {
    _Alignas(Vi) i64 buf[N]; Wi(buf,idx); for(i32 k=0;k<N;++k){simde_mm_prefetch(p+buf[k],_MM_HINT_T0);}
}

#define P(x,m,s) x=pstep(x,m,s)
#define L2(i,s) m=Si(masks[i]); P(b0,m,s);P(b1,m,s)
#define L4(i,s) m=Si(masks[i]); P(b0,m,s);P(b1,m,s);P(b2,m,s);P(b3,m,s)
#define E_(...) __VA_ARGS__
#define LS(L) E_(L(0,1);L(1,2);L(2,4);L(3,8);L(4,16);L(5,32);L(6,16);L(7,8);L(8,4);L(9,2);L(10,1))
// the first row of ctx->masks is always the identity permutation
HEDLEY_NEVER_INLINE void repr2xN(Vi a0, Vi a1, bs_ctx_t const *bs_ctx, u64 *rep, i64 *idx) {
    u64 const *masks = bs_ctx->masks + bs_ctx->n_r; Vi rep0=a0,rep1=a1; Vi gid0=Zi,gid1=Zi; Vi b0,b1;
    for (i64 k = 1; k < bs_ctx->n_m; ++k, masks += bs_ctx->n_r) {
        Vi m; b0=a0,b1=a1; LS(L2);M8 const p0=gt(rep0,b0),p1=gt(rep1,b1);
        rep0=select(p0,b0,rep0),rep1=select(p1,b1,rep1);
        Vi const vk=Si(k); gid0=select(p0,vk,gid0),gid1=select(p1,vk,gid1);
    }
    Wi(rep,rep0),Wi(rep+N,rep1);
    Wi(idx,gid0),Wi(idx+N,gid1);
}


#define Ln1 u64 const *masks=bs_ctx->masks+bs_ctx->n_r; Vi rep0=a0,rep1=a1,rep2=a2,rep3=a3; Vi gid0=Zi,gid1=Zi,gid2=Zi,gid3=Zi; Vi b0,b1,b2,b3
#define Ln2 Vi m; b0=a0,b1=a1,b2=a2,b3=a3; LS(L4)
#define Ln3 M8 const p0=gt(rep0,b0),p1=gt(rep1,b1),p2=gt(rep2,b2),p3=gt(rep3,b3)
#define Ln4 rep0=select(p0,b0,rep0),rep1=select(p1,b1,rep1),rep2=select(p2,b2,rep2),rep3=select(p3,b3,rep3)
#define Ln5 Vi const vk=Si(k); gid0=select(p0,vk,gid0),gid1=select(p1,vk,gid1),gid2=select(p2,vk,gid2),gid3=select(p3,vk,gid3)
#define Ln6 Wi(_rep,rep0),Wi(_rep+N,rep1),Wi(_rep+2*N,rep2),Wi(_rep+3*N,rep3); Wi(_gid,gid0),Wi(_gid+N,gid1),Wi(_gid+2*N,gid2),Wi(_gid+3*N,gid3)
void repr4xN(Vi a0, Vi a1, Vi a2, Vi a3, bs_ctx_t const *bs_ctx, u64 *_rep, i64 *_gid) {
    Ln1;for(i64 k=1; k<bs_ctx->n_m; ++k,masks+=bs_ctx->n_r){Ln2;Ln3;Ln4;Ln5;}Ln6;
}

void repr_search4xN(Vi a0, Vi a1, Vi a2, Vi a3, Vi needle0, Vi needle1, Vi needle2, Vi needle3,
        bs_ctx_t const *bs_ctx, search_ctx_t const *search_ctx, u64 *_rep, i64 *_gid, i64 *_idx) {
    // Setup for search
    i64 n = search_ctx->range_size; Vi values0,values1,values2,values3;
    Vi base0 = and(shr(needle0, search_ctx->shift), Si(search_ctx->mask)),
       base1 = and(shr(needle1, search_ctx->shift), Si(search_ctx->mask)),
       base2 = and(shr(needle2, search_ctx->shift), Si(search_ctx->mask)),
       base3 = and(shr(needle3, search_ctx->shift), Si(search_ctx->mask));
    base0 = gatherq(search_ctx->offsets, base0);
    base1 = gatherq(search_ctx->offsets, base1);
    base2 = gatherq(search_ctx->offsets, base2);
    base3 = gatherq(search_ctx->offsets, base3);
    // Setup for repr
    Ln1;
    for (i64 k = 1; k < bs_ctx->n_m;) {
        if (n > 1) {
            prefetchq(search_ctx->reps+n/2, base0);
            prefetchq(search_ctx->reps+n/2, base1);
            prefetchq(search_ctx->reps+n/2, base2);
            prefetchq(search_ctx->reps+n/2, base3);
        }
        Ln2;
        if (n > 1) {
            values0 = gatherq((i64 const*)search_ctx->reps+n/2, base0);
            values1 = gatherq((i64 const*)search_ctx->reps+n/2, base1);
            values2 = gatherq((i64 const*)search_ctx->reps+n/2, base2);
            values3 = gatherq((i64 const*)search_ctx->reps+n/2, base3);
        }
        Ln3; Ln4; Ln5; ++k,masks+=bs_ctx->n_r;
        if (n > 1) {
            i64 const h = n / 2;
            base0 = select(gt(needle0, values0), addq(base0, Si(h)), base0);
            base1 = select(gt(needle1, values1), addq(base1, Si(h)), base1);
            base2 = select(gt(needle2, values2), addq(base2, Si(h)), base2);
            base3 = select(gt(needle3, values3), addq(base3, Si(h)), base3);
            n -= h;
        }
    }

    while (n > 1) {
        i64 const h = n / 2;
        values0 = gatherq((i64 const*)search_ctx->reps+h, base0);
        values1 = gatherq((i64 const*)search_ctx->reps+h, base1);
        values2 = gatherq((i64 const*)search_ctx->reps+h, base2);
        values3 = gatherq((i64 const*)search_ctx->reps+h, base3);
        base0 = select(gt(needle0, values0), addq(base0, Si(h)), base0);
        base1 = select(gt(needle1, values1), addq(base1, Si(h)), base1);
        base2 = select(gt(needle2, values2), addq(base2, Si(h)), base2);
        base3 = select(gt(needle3, values3), addq(base3, Si(h)), base3);
        n -= h;
    }
    values0 = gatherq((i64 const*)search_ctx->reps, base0);
    values1 = gatherq((i64 const*)search_ctx->reps, base1);
    values2 = gatherq((i64 const*)search_ctx->reps, base2);
    values3 = gatherq((i64 const*)search_ctx->reps, base3);
    base0 = select(gt(needle0, values0), addq(base0, Si(1)), base0);
    base1 = select(gt(needle1, values1), addq(base1, Si(1)), base1);
    base2 = select(gt(needle2, values2), addq(base2, Si(1)), base2);
    base3 = select(gt(needle3, values3), addq(base3, Si(1)), base3);
    values0 = gatherq((i64 const*)search_ctx->reps, base0);
    values1 = gatherq((i64 const*)search_ctx->reps, base1);
    values2 = gatherq((i64 const*)search_ctx->reps, base2);
    values3 = gatherq((i64 const*)search_ctx->reps, base3);

    Ln6;
    Wi(_idx,select(eqq(needle0, values0), base0, Si(-1)));
    Wi(_idx+N,select(eqq(needle1, values1), base1, Si(-1)));
    Wi(_idx+2*N,select(eqq(needle2, values2), base2, Si(-1)));
    Wi(_idx+3*N,select(eqq(needle3, values3), base3, Si(-1)));
}

#define P(x,m,s) x=pstep(x,m,s)
#define L(i,s) m=Si(masks[i]); P(b0,m,s); P(b1,m,s)
// HEDLEY_NEVER_INLINE void repr2xN(Vi a0, Vi a1, bs_ctx_t const *bs_ctx, u64 *rep, i64 *idx) {
//     // the first row of ctx->masks is always the identity permutation
//     u64 const *masks = bs_ctx->masks + bs_ctx->n_r; Vi rep0=a0,rep1=a1; Vi gid0=Zi,gid1=Zi; Vi b0,b1;
//     for (i64 k = 1; k < bs_ctx->n_m; ++k, masks += bs_ctx->n_r) {
//         Vi m; b0=a0,b1=a1; L(0,1); L(1,2); L(2,4); L(3,8); L(4,16); L(5,32); L(6,16); L(7,8); L(8,4); L(9,2); L(10,1);
//         M8 const p0=gt(rep0,b0),p1=gt(rep1,b1);
//         rep0=select(p0,b0,rep0),rep1=select(p1,b1,rep1);
//         gid0=select(p0,Si(k),gid0),gid1=select(p1,Si(k),gid1);
//     }
//     Wi(rep, rep0); Wi(rep + N, rep1);
//     Wi(idx, gid0); Wi(idx + N, gid1);
// }


HEDLEY_NEVER_INLINE void repr_search2xN(Vi a0, Vi a1, Vi needle0, Vi needle1, bs_ctx_t const *bs_ctx, search_ctx_t const *search_ctx,
        u64 *out_rep, i64 *out_gid, i64 *out_idx) {

    // Setup for search
    i64 n = search_ctx->range_size;
    Vi base0 = and(shr(needle0, search_ctx->shift), Si(search_ctx->mask)),
       base1 = and(shr(needle1, search_ctx->shift), Si(search_ctx->mask));
    base0 = gatherq(search_ctx->offsets, base0);
    base1 = gatherq(search_ctx->offsets, base1);

    Vi values0, values1;
    // Setup for repr
    u64 const *masks = bs_ctx->masks + bs_ctx->n_r; Vi rep0=a0,rep1=a1; Vi gid0=Zi,gid1=Zi;
    Vi b0,b1;
    for (i64 k = 1; k < bs_ctx->n_m;) {
        if (n > 1) {
            _Alignas(32) i64 buf[2*N]; Wi(buf, base0), Wi(buf + N, base1);
            for (i32 _k=0;_k<2*N;++_k) {
                simde_mm_prefetch(search_ctx->reps + buf[_k] + n/2, _MM_HINT_T0);
            }
        }
        Vi m; b0=a0,b1=a1;
        L(0,1); L(1,2); L(2,4); L(3,8); L(4,16); L(5,32);
        if (n > 1) {
            values0 = gatherq((i64 const*)search_ctx->reps, addq(base0, Si(n / 2)));
        }
        L(6,16); L(7,8); L(8,4); L(9,2); L(10,1);
        if (n > 1) {
            values1 = gatherq((i64 const*)search_ctx->reps, addq(base1, Si(n / 2)));
        }

        M8 const p0=gt(rep0,b0),p1=gt(rep1,b1);
        rep0=select(p0,b0,rep0),rep1=select(p1,b1,rep1);
        gid0=select(p0,Si(k),gid0),gid1=select(p1,Si(k),gid1);
        ++k; masks += bs_ctx->n_r;
        if (n > 1) {
            i64 const h = n / 2;
            base0 = select(gt(needle0, values0), addq(base0, Si(h)), base0);
            base1 = select(gt(needle1, values1), addq(base1, Si(h)), base1);
            n -= h;
        }
    }
    
    while (n > 1) {
        i64 const h = n / 2;
        values0 = gatherq((i64 const*)search_ctx->reps, addq(base0, Si(h)));
        values1 = gatherq((i64 const*)search_ctx->reps, addq(base1, Si(h)));
        base0 = select(gt(needle0, values0), addq(base0, Si(h)), base0);
        base1 = select(gt(needle1, values1), addq(base1, Si(h)), base1);
        n -= h;
    }
    values0 = gatherq((i64 const*)search_ctx->reps, base0);
    values1 = gatherq((i64 const*)search_ctx->reps, base1);
    base0 = select(gt(needle0, values0), addq(base0, Si(1)), base0);
    base1 = select(gt(needle1, values1), addq(base1, Si(1)), base1);
    values0 = gatherq((i64 const*)search_ctx->reps, base0);
    values1 = gatherq((i64 const*)search_ctx->reps, base1);

    Wi(out_rep,rep0), Wi(out_rep + N,rep1);
    Wi(out_gid,gid0), Wi(out_gid + N,gid1);
    Wi(out_idx,select(eqq(needle0, values0), base0, Si(-1)));
    Wi(out_idx + N,select(eqq(needle1, values1), base1, Si(-1)));
}
#undef P
#undef L

#define Wd_coeff(z,b,t) Wd(acc+(z)*N,coeffd(as[b],ctx->n_s0[t],ctx->n_s1[t],ctx->n_s2[t],ctx->n_sX[t],v_re,v_im,s1,s20,s21,sX))
void off_diagd8xN(u64 const *_alpha, u16 const *_norm, f64 const *_X,
        f64 *out, oc_t const *ctx, const bs_ctx_t *bs_ctx, const search_ctx_t *search_ctx) {
    if (ctx->n_t == 0) { return; }
    f64 const *v_re = ctx->v_re, *v_im = ctx->v_im; u64 const *s1 = ctx->s1, *s20 = ctx->s20, *s21 = ctx->s21, *sX = ctx->sX;
    _Alignas(Vi) u64 rep[4*N]; _Alignas(Vi) i64 gid[4*N]; _Alignas(Vi) i64 idx[4*N]; _Alignas(Vi) f64 acc[4*N]; _Alignas(Vi) f64 chi[4*N];
    _Alignas(Vi) f64 norm0[8*N]; for(i32 k=0;k<8*N;++k){norm0[k]=_norm[k];}
    f64 outer[8*N]; for(i32 k=0;k<8;++k){Wd(outer+k*N,Zd);}
    Vi bs[4]; Vi as[8]; for(i32 k=0;k<8;++k){as[k]=Ri(_alpha+k*N);}

    i32 j = 0;
    if (bs_ctx != NULL) {
        for (i32 k=0;k<4;++k){bs[k]=xor(as[k],Si(ctx->mask[0]));}
        repr4xN(bs[0],bs[1],bs[2],bs[3],bs_ctx,rep,gid);
        j += 4;
    }

    for (; j < 8 * ctx->n_t; j += 4) {
        i32 const ti = j / 8, bi = j % 8, oti = (j - 4) / 8, obi = (j - 4) % 8;
        for (i32 k=0;k<4;++k){bs[k]=xor(as[bi+k],Si(ctx->mask[ti]));}
        if (bs_ctx != NULL) {
            for(i32 k=0;k<4*N;++k){chi[k]=bs_ctx->chi_re[gid[k]];}
            repr_search4xN(bs[0],bs[1],bs[2],bs[3],Ri(rep),Ri(rep+N),Ri(rep+2*N),Ri(rep+3*N),bs_ctx,search_ctx,rep,gid,idx);
            for(i32 k=0;k<4*N;++k){i32 const i=idx[k];if(i>=0){simde_mm_prefetch(search_ctx->norm+i,_MM_HINT_ET0);simde_mm_prefetch(_X+i,_MM_HINT_ET0);}}
            Wd_coeff(0,obi,oti),Wd_coeff(1,obi+1,oti),Wd_coeff(2,obi+2,oti),Wd_coeff(3,obi+3,oti);
            for(i32 k=0;k<4*N;++k){
                if(idx[k]>=0){
                    f64 const norm=search_ctx->norm[idx[k]], x=_X[idx[k]];
                    outer[obi*N+k]+=acc[k]*chi[k]*__builtin_sqrt(norm/norm0[obi*N+k])*x;
                }
            }
            if(obi!=0){
                s1+=ctx->n_s1[oti], s20+=ctx->n_s2[oti], s21+=ctx->n_s2[oti], sX+=ctx->n_sX[oti];
                i32 const k=ctx->n_s0[oti]+ctx->n_s1[oti]+ctx->n_s2[oti]+ctx->n_sX[oti]; v_re += k, v_im += k;
            }
        }
        else {
            for(i32 k=0;k<4;++k){Wi(idx+k*N,bs[k]);}
            for(i32 k=0;k<4*N;++k){simde_mm_prefetch(_X+idx[k],_MM_HINT_T0);}
            Wd(acc, coeffd(as[bi], ctx->n_s0[ti], ctx->n_s1[ti], ctx->n_s2[ti], ctx->n_sX[ti], v_re, v_im, s1, s20, s21, sX));
            Wd(acc+N, coeffd(as[bi+1], ctx->n_s0[ti], ctx->n_s1[ti], ctx->n_s2[ti], ctx->n_sX[ti], v_re, v_im, s1, s20, s21, sX));
            Wd(acc+2*N, coeffd(as[bi+2], ctx->n_s0[ti], ctx->n_s1[ti], ctx->n_s2[ti], ctx->n_sX[ti], v_re, v_im, s1, s20, s21, sX));
            Wd(acc+3*N, coeffd(as[bi+3], ctx->n_s0[ti], ctx->n_s1[ti], ctx->n_s2[ti], ctx->n_sX[ti], v_re, v_im, s1, s20, s21, sX));
            for(i32 k=0;k<4*N;++k){f64 const x=_X[idx[k]];outer[bi*N+k]+=acc[k]*x;}
            if(bi!=0){
                s1+=ctx->n_s1[ti], s20+=ctx->n_s2[ti], s21+=ctx->n_s2[ti], sX+=ctx->n_sX[ti];
                i32 const k=ctx->n_s0[ti]+ctx->n_s1[ti]+ctx->n_s2[ti]+ctx->n_sX[ti]; v_re += k, v_im += k;
            }
        }
    }

    if (bs_ctx != NULL) {
        i32 const oti = ctx->n_t - 1;
        { Vi _msk,_idx;unpack_pair(_msk,_idx,searchN(Ri(rep),search_ctx));Wi(idx,select(_msk,_idx,Si(-1))); }
        { Vi _msk,_idx;unpack_pair(_msk,_idx,searchN(Ri(rep+N),search_ctx));Wi(idx+N,select(_msk,_idx,Si(-1))); }
        { Vi _msk,_idx;unpack_pair(_msk,_idx,searchN(Ri(rep+2*N),search_ctx));Wi(idx+2*N,select(_msk,_idx,Si(-1))); }
        { Vi _msk,_idx;unpack_pair(_msk,_idx,searchN(Ri(rep+3*N),search_ctx));Wi(idx+3*N,select(_msk,_idx,Si(-1))); }
        Wd(acc, coeffd(as[4], ctx->n_s0[oti], ctx->n_s1[oti], ctx->n_s2[oti], ctx->n_sX[oti], v_re, v_im, s1, s20, s21, sX));
        Wd(acc+N, coeffd(as[5], ctx->n_s0[oti], ctx->n_s1[oti], ctx->n_s2[oti], ctx->n_sX[oti], v_re, v_im, s1, s20, s21, sX));
        Wd(acc+2*N, coeffd(as[6], ctx->n_s0[oti], ctx->n_s1[oti], ctx->n_s2[oti], ctx->n_sX[oti], v_re, v_im, s1, s20, s21, sX));
        Wd(acc+3*N, coeffd(as[7], ctx->n_s0[oti], ctx->n_s1[oti], ctx->n_s2[oti], ctx->n_sX[oti], v_re, v_im, s1, s20, s21, sX));
        for(i32 k=0;k<4*N;++k){
            if(idx[k]>=0){
                f64 const chi=bs_ctx->chi_re[gid[k]], norm=search_ctx->norm[idx[k]], x=_X[idx[k]];
                outer[4*N+k]+=acc[k]*chi*__builtin_sqrt(norm/norm0[4*N+k])*x;
            }
        }
    }

    for(i32 k=0;k<8;++k){Wd(out+k*N,addd(Rd(out+k*N),Rd(outer+k*N)));}
}

INTERNAL void off_diagd2xN(u64 const *_alpha, u16 const *_norm, f64 const *_X,
        f64 *out, oc_t const *ctx, const bs_ctx_t *bs_ctx, const search_ctx_t *search_ctx) {
    if (ctx->n_t == 0) { return; }

    // Vi const a0 = Ri(_alpha), a1 = Ri(_alpha + N), a2 = Ri(_alpha + 2 * N), a3 = Ri(_alpha + 3 * N);
    f64 const *v_re = ctx->v_re, *v_im = ctx->v_im; u64 const *s1 = ctx->s1, *s20 = ctx->s20, *s21 = ctx->s21, *sX = ctx->sX;

    _Alignas(Vi) u64 rep[2*N];
    _Alignas(Vi) i64 gid[2*N];
    _Alignas(Vi) i64 idx[2*N];
    _Alignas(Vi) f64 acc[2*N];
    _Alignas(Vi) f64 chi[2*N];
    f64 outer[4*N]; Wd(outer, Zd); Wd(outer+N, Zd); Wd(outer+2*N, Zd); Wd(outer+3*N, Zd);

    // linearized index j
    // j % 4 -> [beta0 beta1] or [beta2 beta3]
    // j // 4 -> ti

    Vi as[4] = {Ri(_alpha), Ri(_alpha + N), Ri(_alpha + 2 * N), Ri(_alpha + 3 * N)};
    Vi bs[2];

    i32 j = 0;
    if (bs_ctx != NULL) {
        bs[0] = xor(as[0], Si(ctx->mask[0])), bs[1] = xor(as[1], Si(ctx->mask[0]));
        repr2xN(bs[0], bs[1], bs_ctx, rep, gid);
        j += 2;
    }

    for (; j < 4 * ctx->n_t; j += 2) {
        i32 const ti = j / 4, bi = j % 4, oti = (j - 2) / 4, obi = (j - 2) % 4;
        bs[0] = xor(as[bi], Si(ctx->mask[ti])), bs[1] = xor(as[bi+1], Si(ctx->mask[ti]));
        if (bs_ctx != NULL) {
            for(i32 k=0;k<2*N;++k){chi[k]=bs_ctx->chi_re[gid[k]];}
            repr_search2xN(bs[0], bs[1], Ri(rep), Ri(rep+N), bs_ctx, search_ctx, rep, gid, idx);
            for(i32 k=0;k<2*N;++k){i32 const i=idx[k];if(i>=0){simde_mm_prefetch(search_ctx->norm+i,_MM_HINT_T0);simde_mm_prefetch(_X+i,_MM_HINT_T0);}}
            Wd(acc, coeffd(as[obi], ctx->n_s0[oti], ctx->n_s1[oti], ctx->n_s2[oti], ctx->n_sX[oti], v_re, v_im, s1, s20, s21, sX));
            Wd(acc+N, coeffd(as[obi+1], ctx->n_s0[oti], ctx->n_s1[oti], ctx->n_s2[oti], ctx->n_sX[oti], v_re, v_im, s1, s20, s21, sX));
            for(i32 k=0;k<2*N;++k){
                if(idx[k]>=0){
                    f64 const n0=_norm[obi*N+k], norm=search_ctx->norm[idx[k]], x=_X[idx[k]];
                    outer[obi*N+k]+=acc[k]*chi[k]*__builtin_sqrt(norm/n0)*x;
                }
            }
            if(obi!=0){
                s1+=ctx->n_s1[oti], s20+=ctx->n_s2[oti], s21+=ctx->n_s2[oti], sX+=ctx->n_sX[oti];
                i32 const k=ctx->n_s0[oti]+ctx->n_s1[oti]+ctx->n_s2[oti]+ctx->n_sX[oti]; v_re += k, v_im += k;
            }
        }
        else {
            Wi(idx,bs[0]), Wi(idx+N,bs[1]);
            for(i32 k=0;k<2*N;++k){simde_mm_prefetch(_X+idx[k],_MM_HINT_T0);}
            Wd(acc, coeffd(as[bi], ctx->n_s0[ti], ctx->n_s1[ti], ctx->n_s2[ti], ctx->n_sX[ti], v_re, v_im, s1, s20, s21, sX));
            Wd(acc+N, coeffd(as[bi+1], ctx->n_s0[ti], ctx->n_s1[ti], ctx->n_s2[ti], ctx->n_sX[ti], v_re, v_im, s1, s20, s21, sX));
            for(i32 k=0;k<2*N;++k){f64 const x=_X[idx[k]];outer[bi*N+k]+=acc[k]*x;}
            if(bi!=0){
                s1+=ctx->n_s1[ti], s20+=ctx->n_s2[ti], s21+=ctx->n_s2[ti], sX+=ctx->n_sX[ti];
                i32 const k=ctx->n_s0[ti]+ctx->n_s1[ti]+ctx->n_s2[ti]+ctx->n_sX[ti]; v_re += k, v_im += k;
            }
        }
    }

    if (bs_ctx != NULL) {
        i32 const oti = ctx->n_t - 1;
        { Vi _msk,_idx;unpack_pair(_msk,_idx,searchN(Ri(rep),search_ctx));Wi(idx,select(_msk,_idx,Si(-1))); }
        { Vi _msk,_idx;unpack_pair(_msk,_idx,searchN(Ri(rep+N),search_ctx));Wi(idx+N,select(_msk,_idx,Si(-1))); }
        Wd(acc, coeffd(as[2], ctx->n_s0[oti], ctx->n_s1[oti], ctx->n_s2[oti], ctx->n_sX[oti], v_re, v_im, s1, s20, s21, sX));
        Wd(acc+N, coeffd(as[3], ctx->n_s0[oti], ctx->n_s1[oti], ctx->n_s2[oti], ctx->n_sX[oti], v_re, v_im, s1, s20, s21, sX));
        for(i32 k=0;k<2*N;++k){
            if(idx[k]>=0){
                f64 const chi=bs_ctx->chi_re[gid[k]], n0=_norm[2*N+k], norm=search_ctx->norm[idx[k]], x=_X[idx[k]];
                outer[2*N+k]+=acc[k]*chi*__builtin_sqrt(norm/n0)*x;
            }
        }
    }

    for(i32 k=0;k<4;++k){Wd(out+k*N,addd(Rd(out+k*N),Rd(outer+k*N)));}
    
    // for (i32 ti = 0; ti < ctx->n_t; ++ti) {
    //     if (bs_ctx != NULL) {
    //          beta2 = xor(a2, Si(ctx->mask[ti])),
    //          beta3 = xor(a3, Si(ctx->mask[ti]));
    //         // { Vi _rep, _gid; unpack_pair(_rep, _gid, reprN(beta0, bs_ctx)); Wi(rep,_rep); Wi(gid,_gid); }
    //         // { Vi _rep, _gid; unpack_pair(_rep, _gid, reprN(beta1, bs_ctx)); Wi(rep + N,_rep); Wi(gid + N,_gid); }
    //         repr2xN(beta0, beta1, bs_ctx, rep, gid);
    //         repr_search2xN(beta2, beta3, Ri(rep + 0 * N), Ri(rep + 1 * N), bs_ctx, search_ctx,
    //                        rep + 2 * N, gid + 2 * N, idx + 0 * N);
    //         // repr_search(beta2, Ri(rep + 1 * N), bs_ctx, search_ctx, rep + 2 * N, gid + 2 * N, idx + 1 * N);
    //         // repr_search(beta3, Ri(rep + 2 * N), bs_ctx, search_ctx, rep + 3 * N, gid + 3 * N, idx + 2 * N);
    //         { Vi _msk, _idx; unpack_pair(_msk, _idx, searchN(Ri(rep + 2 * N), search_ctx)); Wi(idx + 2 * N, select(_msk, _idx, Si(-1))); }
    //         { Vi _msk, _idx; unpack_pair(_msk, _idx, searchN(Ri(rep + 3 * N), search_ctx)); Wi(idx + 3 * N, select(_msk, _idx, Si(-1))); }
    //         // Pi(Ri(rep)); Pi(Ri(rep + N));
    //     }
    //     else {
    //         Wi(idx, xor(a0, Si(ctx->mask[ti])));
    //         Wi(idx + N, xor(a1, Si(ctx->mask[ti])));
    //         Wi(idx + 2 * N, xor(a2, Si(ctx->mask[ti])));
    //         Wi(idx + 3 * N, xor(a3, Si(ctx->mask[ti])));
    //     }

    //     for(i32 k=0;k<4*N;++k) {
    //         if (idx[k] >= 0) {
    //             if (bs_ctx != NULL) { simde_mm_prefetch(search_ctx->norm + idx[k], _MM_HINT_T0); }
    //             simde_mm_prefetch((f64 const*)_X + idx[k], _MM_HINT_T0);
    //         }
    //     }
    //     Wd(acc, coeffd(a0, ctx->n_s0[ti], ctx->n_s1[ti], ctx->n_s2[ti], ctx->n_sX[ti], v_re, v_im, s1, s20, s21, sX));
    //     Wd(acc + N, coeffd(a1, ctx->n_s0[ti], ctx->n_s1[ti], ctx->n_s2[ti], ctx->n_sX[ti], v_re, v_im, s1, s20, s21, sX));
    //     Wd(acc + 2 * N, coeffd(a2, ctx->n_s0[ti], ctx->n_s1[ti], ctx->n_s2[ti], ctx->n_sX[ti], v_re, v_im, s1, s20, s21, sX));
    //     Wd(acc + 3 * N, coeffd(a3, ctx->n_s0[ti], ctx->n_s1[ti], ctx->n_s2[ti], ctx->n_sX[ti], v_re, v_im, s1, s20, s21, sX));
    //     if (bs_ctx != NULL) {
    //         for(i32 k=0;k<4*N;++k) {
    //             if (idx[k] >= 0) {
    //                 f64 const chi = bs_ctx->chi_re[gid[k]];
    //                 f64 const n0 = _norm[k];
    //                 f64 const norm = search_ctx->norm[idx[k]];
    //                 f64 const x = _X[idx[k]];
    //                 outer[k] += (acc[k] * chi * __builtin_sqrt(norm / n0)) * x;
    //             }
    //         }
    //     }
    //     else {
    //         for(i32 k=0;k<4*N;++k) {
    //             f64 const x = _X[idx[k]];
    //             outer[k] += acc[k] * x;
    //         }
    //     }
    //     s1 += ctx->n_s1[ti]; s20 += ctx->n_s2[ti]; s21 += ctx->n_s2[ti]; sX += ctx->n_sX[ti];
    //     i32 const k = ctx->n_s0[ti] + ctx->n_s1[ti] + ctx->n_s2[ti] + ctx->n_sX[ti];
    //     v_re += k; v_im += k;
    // }
    // for (i32 k=0;k<4;++k) { Wd(out + k * N, addd(Rd(out + k * N), Rd(outer + k * N))); }
}

#define off_diagN_template(s, t) \
    INTERNAL void off_diag##s##N(u64 const *alpha0, u16 const *norm0, void const *X, \
            void *out, oc_t const *ctx, const bs_ctx_t *bs_ctx, const search_ctx_t *search_ctx) { \
        Vi const alpha = Ri(alpha0); Vd n0; V##s acc_outer = zero##s(); \
        if (bs_ctx != NULL) { n0 = load_norm(norm0); } \
        f64 const *v_re = ctx->v_re, *v_im = ctx->v_im; u64 const *s1 = ctx->s1, *s20 = ctx->s20, *s21 = ctx->s21, *sX = ctx->sX; \
        _Alignas(32) t chi[N]; _Alignas(32) f64 norm[N]; \
        _Alignas(32) i64 gid[N]; _Alignas(32) t _x[N]; \
        _Alignas(32) i64 msk[N]; _Alignas(32) i64 idx[N]; \
        for (i32 ti = 0; ti < ctx->n_t; ++ti) { \
            Vi const beta0 = xor(alpha, Si(ctx->mask[ti])); V##s x; \
            if (bs_ctx != NULL) { \
                Vi rep, _gid; unpack_pair(rep, _gid, reprN(beta0, bs_ctx)); \
                Wi(gid, _gid); \
                Vi _msk, _idx; unpack_pair(_msk, _idx, searchN(rep, search_ctx)); \
                Wi(msk, _msk); Wi(idx, _idx); \
                for(i32 k=0;k<N;++k){ \
                    if (msk[k] != 0) { \
                        simde_mm_prefetch(search_ctx->norm + idx[k], _MM_HINT_T0); \
                        simde_mm_prefetch((t const*)X + idx[k], _MM_HINT_T0); \
                    } \
                    chi[k] = c2##s(bs_ctx->chi_re[gid[k]], bs_ctx->chi_im[gid[k]]); \
                } \
                /*norm = mask_gather_norm(search_ctx->norm, idx, msk);*/ \
                /*x = mask_gather##s((t const*)X, idx, msk);*/ \
                /*chi = gather2##s(bs_ctx->chi_re, bs_ctx->chi_im, gid);*/ \
            } \
            else { \
                x = gather##s((t const*)X, beta0); \
            } \
            V##s const acc = coeff##s(alpha, ctx->n_s0[ti], ctx->n_s1[ti], ctx->n_s2[ti], ctx->n_sX[ti], \
                v_re, v_im, s1, s20, s21, sX); \
            if (bs_ctx != NULL) { \
                for(i32 k=0;k<N;++k){ \
                    if (msk[k] != 0) { \
                        norm[k] = search_ctx->norm[idx[k]]; \
                        _x[k] = ((t const*)X)[idx[k]]; \
                    } \
                    else { norm[k] = 0; _x[k] = 0; } \
                } \
                Vd const c = sqrtd(divd(Rd(norm), n0)); \
                x = scale##s(c, mul##s(R##s(chi), R##s(_x))); \
            } \
            acc_outer = add##s(acc_outer, mul##s(acc, x)); \
            s1 += ctx->n_s1[ti]; s20 += ctx->n_s2[ti]; s21 += ctx->n_s2[ti]; sX += ctx->n_sX[ti]; \
            i32 const k = ctx->n_s0[ti] + ctx->n_s1[ti] + ctx->n_s2[ti] + ctx->n_sX[ti]; \
            v_re += k; v_im += k; \
        } \
        W##s(out, add##s(R##s(out), acc_outer)); \
    }
// off_diagN_template(d, f64)
off_diagN_template(z, c128)
#undef off_diagN_template

void off_diag64_f64(u64 const *alpha0, u16 const *norm0, void const *X, void *out,
        oc_t const *ctx, bs_ctx_t const *bs_ctx, search_ctx_t const* search_ctx) {
    for (i32 k = 0; k < 64; k += 8*N) { off_diagd8xN(alpha0 + k, norm0 + k, X, (f64*)out + k, ctx, bs_ctx, search_ctx); }
}
#define off_diag64_template(s, t) \
    void off_diag64_##t(u64 const *alpha0, u16 const *norm0, void const *X, void *out, \
            oc_t const *ctx, bs_ctx_t const *bs_ctx, search_ctx_t const* search_ctx) { \
        for (i32 k = 0; k < 64; k += N) { off_diag##s##N(alpha0 + k, norm0 + k, X, (t*)out + k, ctx, bs_ctx, search_ctx); } \
    }
// off_diag64_template(d, f64)
off_diag64_template(z, c128)
#undef off_diag64_template

#define matvec_inner_template(t) \
    void matvec_inner_##t(i64 const i, u64 const *alpha0, u16 const *norm0, void const *X0, void const *X, void *out, \
            oc_t const *diag, oc_t const *off_diag, bs_ctx_t const *bs, search_ctx_t const *search) { \
        diag64_##t(i, alpha0, X0, (t*)out + i, diag); off_diag64_##t(alpha0 + i, norm0 + i, X, (t*)out + i, off_diag, bs, search); \
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



