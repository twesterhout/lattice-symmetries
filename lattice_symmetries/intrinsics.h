#pragma once

// #include <stdio.h>
// #include <simde/x86/avx2.h>
#include <simde/x86/avx512.h>
// #include <simde/x86/fma.h>
// #include <simde/x86/avx512.h>

#define INTERNAL static HEDLEY_ALWAYS_INLINE

typedef int8_t i8; typedef int16_t i16; typedef int32_t i32; typedef int64_t i64;
typedef uint8_t u8; typedef uint16_t u16; typedef uint32_t u32; typedef uint64_t u64;
typedef float f32; typedef double f64; typedef float _Complex c64; typedef double _Complex c128;


#if defined(__clang__)
#define assume(cond) __builtin_assume(cond)
#else
#define assume(cond) __attribute__((__assume__(cond)))
#endif
#define M 2

// always operate on B bits of data

// #define M 1
// 1: AVX512    | 512
// 2: AVX2+FMA  | 256
// 3: SSE       | 
// 4: NEON      | 128
// 5: WASM      |
#define I_(b,f,...) simde_mm##b##_##f(__VA_ARGS__)
#define I(b,f,...) I_(b,f,__VA_ARGS__)

#if M == 1
#define B 512
#define _S(f) f##_si512
#elif M == 2
#define B 256
#define _S(f) f##_si256
#else
#define B 128
#define _S(f) f##_si128
#endif
#define N (B / 64)

#define T__(b) typedef simde__m##b##i Vi; typedef simde__m##b##d Vd
#define T(b) T__(b)
T(B);
#undef T__
#undef T
typedef struct Vz { Vd re; Vd im; } Vz;
#if M == 1
    typedef simde__mmask8 M8;
#else
    typedef Vi M8;
#endif

#define Zi I(B,_S(setzero))
#define Zd I(B,setzero_pd)
#if M == 2
#  define Si(x) I(B,set1_epi64x,x)
#else
#  define Si(x) I(B,set1_epi64,x)
#endif
#define shl(x,imm) I(B,slli_epi64,x,imm)
#define shr(x,imm) I(B,srli_epi64,x,imm)

#define F2_(t,f,i) INTERNAL t f(t const a, t const b) { return I(B,i,a,b); }
#define F2(t,f,i) F2_(t,f,i)

F2(Vi,and,_S(and))
F2(Vi,xor,_S(xor))
F2(Vi,andnot,_S(andnot))
F2(Vi,or,_S(or))
F2(Vi,addq,add_epi64)
F2(Vi,subq,sub_epi64)
INTERNAL Vi popcnt(Vi const x) { return I(B,popcnt_epi64,x); }
INTERNAL Vi Ri(void const *p) { return I(B,loadu_epi64,p); }
INTERNAL void Wi(void *p, Vi const x) { I(B,_S(storeu),p,x); }

F2(Vd,addd,add_pd)
F2(Vd,subd,sub_pd)
F2(Vd,muld,mul_pd)
F2(Vd,divd,div_pd)
INTERNAL Vd fmad(Vd const a, Vd const b, Vd const c) { return I(B,fmadd_pd,a,b,c); }
INTERNAL Vd fnmad(Vd const a, Vd const b, Vd const c) { return I(B,fnmadd_pd,a,b,c); }
INTERNAL Vd scaled(Vd const a, Vd const b) { return muld(a, b); }
INTERNAL Vd sqrtd(Vd const x) { return I(B,sqrt_pd,x); }
INTERNAL Vd Rd(f64 const *p) { return I(B,loadu_pd,p); }
INTERNAL void Wd(f64 *p, Vd const x) { I(B,storeu_pd,p,x); }
INTERNAL Vd Sd(f64 const x) { return I(B,set1_pd,x); }

INTERNAL Vz addz(Vz const a, Vz const b) { return (Vz){addd(a.re, b.re), addd(a.im, b.im)}; }
INTERNAL Vz mulz(Vz const a, Vz const b) { return (Vz){fnmad(a.im, b.im, muld(a.re, b.re)), fmad(a.im, b.re, muld(a.re, b.im))}; }
INTERNAL Vz scalez(Vd const a, Vz const b) { return (Vz){scaled(a, b.re), scaled(a, b.im)}; }
INTERNAL Vz Rz(c128 const *p) {
// #if M == 2
//     Vd const a = Rd((f64 const*)p), b = Rd((f64 const*)p + N);
//     Vd const c = OP(unpacklo_pd, a, b), d = OP(unpackhi_pd, a, b);
//     return (Vz){OP(permute4x64_pd, c, 0b11011000), OP(permute4x64_pd, d, 0b11011000)};
// #else
    double re[N]; double im[N]; for (i32 i = 0; i < N; ++i) { double _Complex const z = p[i]; re[i] = __real__ z; im[i] = __imag__ z; }
    return (Vz){Rd(re), Rd(im)};
// #endif
}
INTERNAL void Wz(c128 *p, Vz const z) {
//    Vd const a = OP(permute4x64_pd, z.re, 0b11011000), b = OP(permute4x64_pd, z.im, 0b11011000);
//    Vd const c = OP(unpacklo_pd, a, b), d = OP(unpackhi_pd, a, b);
//    stored((f64*)p, c);
//    stored((f64*)p + N, d);
    double re[N]; double im[N]; Wd(re, z.re), Wd(im, z.im);
    for (i32 i = 0; i < N; ++i) { p[i] = __builtin_complex(re[i], im[i]); }
}

#if M == 1
#  define d2i(x) I(B,castpd_si512,x)
#  define i2d(x) I(B,castsi512_pd,x)
#elif M == 2
#  define d2i(x) I(B,castpd_si256,x)
#  define i2d(x) I(B,castsi256_pd,x)
#endif

#if M == 1
#  define _gatherq(p,i,s) I(B,i64gather_epi64,i,p,s)
#  define _mask_gatherq(p,i,m,s) I(B,mask_i64gather_epi64,Zi,m,i,p,s)
#  define gatherd(p,i) I(B,i64gather_pd,i,p,8)
#  define mask_gatherd(p,i,m) I(B,mask_i64gather_pd,Zd,m,i,p,8)
#elif M == 2
#  define _gatherq(p,i,s) I(B,i64gather_epi64,p,i,s)
#  define _mask_gatherq(p,i,m,s) I(B,mask_i64gather_epi64,Zi,p,i,m,s)
#  define gatherd(p,i) I(B,i64gather_pd,p,i,8)
#  define mask_gatherd(p,i,m) I(B,mask_i64gather_pd,Zd,p,i,i2d(m),8)
#endif
#define mask_gatherz(p,i,m) (Vz){mask_gatherd((f64 const*)p, shl(i, 1), m), mask_gatherd((f64 const*)p + 1, shl(i, 1), m)}
INTERNAL Vi gatherq(i64 const *p, Vi const i) { return _gatherq(p, i, 8); }
INTERNAL Vz gatherz(c128 const *p, Vi const i) { Vi const i2 = shl(i, 1); return (Vz){gatherd((f64 const*)p, i2), gatherd((f64 const*)p + 1, i2)}; }
INTERNAL Vd gather2d(f64 const *re, f64 const *im, Vi const i) { return gatherd(re, i); }
INTERNAL Vz gather2z(f64 const *re, f64 const *im, Vi const i) { return (Vz){gatherd(re, i), gatherd(im, i)}; }
INTERNAL Vd pq2pd(Vi x) { x = or(x, d2i(Sd(0x0010000000000000))); return subd(i2d(x), Sd(0x0010000000000000)); }
// INTERNAL Vd gather_norm(u16 const *p, Vi const i) { return pq2pd(and(_gatherq((i64 const*)p, i, 2), Si(0xFFFF))); }
INTERNAL Vd load_norm(u16 const *p) {
#if M == 1
    return I(B,set_pd,p[7],p[6],p[5],p[4],p[3],p[2],p[1],p[0]);
#elif M == 2
    return I(B,set_pd,p[3],p[2],p[1],p[0]);
#endif
}
#define mask_gather_norm(p,i,m) pq2pd(and(_mask_gatherq((i64 const*)p, i, m, 2), Si(0xFFFF)))

#if M == 1
#  define select(s,a,b) I(B,mask_blend_epi8,s,b,a)
#elif M == 2
#  define select(s,a,b) I(B,blendv_epi8,b,a,s)
#endif

INTERNAL M8 eqq(Vi const a, Vi const b) {
#if M == 1
    return I(B,cmp_epi64_mask,a,b,0);
#elif M == 2
    return I(B,cmpeq_epi64,a,b);
#endif
}

#if M == 1
#  define eqi(a,b) I(B,movm_epi64,eqq(a,b))
#elif M == 2
#  define eqi(a,b) eqq(a,b)
#endif

#define NOT_(a) xor(a, eqi(Zi, Zi))
#define NOT(a) NOT_(a)
INTERNAL Vi m1(Vi const x, Vi const m) {
#if M == 1
    return I(B,movm_epi64,I(B,test_epi64_mask,x,m));
#elif M == 2
    return xor(eqi(and(x, m), Zi), eqi(Zi, Zi));
#endif
}
INTERNAL Vi m2(Vi const x, Vi const m_0, Vi const m_1) { return xor(m1(x, m_0), m1(x, m_1)); }
INTERNAL Vi mX(Vi const x, Vi const m) { return popcnt(and(x, m)); }
INTERNAL Vd signedd(Vd const v, Vi m) { m = shl(m, 63); return i2d(xor(d2i(v), m)); }
INTERNAL Vz signedz(Vz const v, Vi m) { m = shl(m, 63); return (Vz){i2d(xor(d2i(v.re), m)), i2d(xor(d2i(v.im), m))}; }

#if M == 1
#  define gt(a,b) I(B,cmp_epi64_mask,b,a,1)
#elif M == 2
#  define gt(a,b) I(B,cmpgt_epi64,a,b)
#endif

#define PX(x, t, w, f) \
    do { \
        t temp[B / (8 * sizeof(t))]; I(B,w,temp,x); \
        printf("["); for (u64 k = 0; k < B / (8 * sizeof(t)); ++k) { printf(f ",", temp[k]); } printf("]\n"); \
    } while(0)
#define Pw(x) PX(x, int32_t, storeu_epi32, "%i")
#define Pi(x) PX(x, int64_t, storeu_epi64, "%zi")
#define Ps(x) PX(x, float, storeu_ps, "%e")
#define Pd(x) PX(x, double, storeu_pd, "%e")


#define _(z) ({z;})
#define $(b,z) if(b){z;}else
#define D(t,g,k,x...) static HEDLEY_ALWAYS_INLINE t g(x){return _(k);}
#define De(t,g,k,x...) t g(x){return _(k);}
#define Dd(g,k,x...) D(Vd,g##d,k,x)
#define Dz(g,k,x...) D(Vz,g##z,k,x)
#define Z2(x,y) (Vz){x,y}
#define f64c f64 const
#define c(t) t const
// For loop over a counter v of length n
#define _L(v,n,x...) for(i32 v=0;v<(n);++v){x;}



// ===================================================================================

#if 0
#include <stdio.h>
#include <simde/x86/sse.h>
#include <simde/x86/avx2.h>
#include <simde/x86/fma.h>
#include <simde/x86/avx512.h>

#define INTERNAL static HEDLEY_ALWAYS_INLINE


// always operate on B bits of data

#define M 1
// 1: AVX512    | 512
// 2: AVX2+FMA  | 256
// 3: SSE       | 
// 4: NEON      | 128
// 5: WASM      |
#define I_(b,f,t,...) simde_mm##b##_##f##_##t(__VA_ARGS__)
#define I(b,f,t,...) I_(b,f,t,__VA_ARGS__)

#if M == 1
#define B 512
#define Bt si512
#elif M == 2
#define B 256
#define Bt si256
#else
#define B 128
#define Bt si128
#endif

INTERNAL Vi and(Vi const a, Vi const b) { return I(B,and,Bt); }
INTERNAL Vi xor(Vi const a, Vi const b) { return I(B,xor,Bt); }

#if 0
// #define N 8
// #define Vi simde__m512i
// #define Vd simde__m512d
// #define OP(s, ...) simde_mm512_##s(__VA_ARGS__)
// #define Zi OP(setzero_si512)
// #define Si(x) OP(set1_epi64, x)
// #define AND(a, b) OP(and_si512, a, b)
// #define XOR(a, b) OP(xor_si512, a, b)
// #define d2i(x) OP(castpd_si512, x)
// #define i2d(x) OP(castsi512_pd, x)
// #define GATHER(x, i) OP(i64gather_pd, i, x, 8)

// #if defined(SIMDE_X86_AVX512F_NATIVE)
// #define EQ(a, b) _mm512_movm_epi64(_mm512_cmpeq_epi64_mask(a, b))
// #else
// static HEDLEY_ALWAYS_INLINE simde__m512i EQ(simde__m512i a, simde__m512i b) {
//     simde__m512i_private r_, a_ = simde__m512i_to_private(a), b_ = simde__m512i_to_private(b);
//     r_.m256i[0] = simde_mm256_cmpeq_epi64(a_.m256i[0], b_.m256i[0]);
//     r_.m256i[1] = simde_mm256_cmpeq_epi64(a_.m256i[1], b_.m256i[1]);
//     return simde__m512i_from_private(r_);
// }
// #endif

#else

#define N 4
#define Vi simde__m256i
#define Vd simde__m256d
#define OP(s, ...) simde_mm256_##s(__VA_ARGS__)
#define Zi OP(setzero_si256)
#define Si(x) OP(set1_epi64x, x)
#define AND(a, b) OP(and_si256, a, b)
#define OR(a, b) OP(or_si256, a, b)
#define ANDNOT(a, b) OP(andnot_si256, a, b)
#define XOR(a, b) OP(xor_si256, a, b)
#define d2i(x) OP(castpd_si256, x)
#define i2d(x) OP(castsi256_pd, x)
#define EQ(a, b) OP(cmpeq_epi64, a, b)
#define GT(a, b) OP(cmpgt_epi64, a, b)
#define GATHER(x, i) OP(i64gather_pd, x, i, 8)

#define gatherqpd(x, i) OP(i64gather_pd, x, i, 8)
#define gatherqq(x, i) OP(i64gather_epi64, x, i, 8)
#define subq(a, b) OP(sub_epi64, a, b)
#define mulpd(a, b) OP(mul_pd, a, b)
#define sqrtpd(x) OP(sqrt_pd, x)

#endif

#define Zd OP(setzero_pd)
#define Ri(p) OP(loadu_epi64, p)
#define Rd(p) OP(loadu_pd, p)
#define Wi(p, x) OP(storeu_epi64, p, x)
#define Wd(p, x) OP(storeu_pd, p, x)
#define Sd(x) OP(set1_pd, x)

#define LT(a, b) GT(b, a)
#define NOT(a) XOR(a, EQ(Zi, Zi))
#define SHL(x, imm) OP(slli_epi64, x, imm)
#define SHR(x, imm) OP(srli_epi64, x, imm)
#define ADD(a, b) OP(add_pd, a, b)
#define MUL(a, b) OP(mul_pd, a, b)
#define divpd(a, b) OP(div_pd, a, b)
#define FMA(a, b, c) OP(fmadd_pd, a, b, c)
#define POPCNT(x) OP(popcnt_epi64, x)
#define SELECT(s, a, b) OP(blendv_epi8, b, a, s)

#define addq(a, b) OP(add_epi64, a, b)

#define ASSUME(cond) __attribute__((__assume__(cond)))

#define PX(x, t, w, f) \
    do { \
        t temp[B / (8 * sizeof(t))]; OP(w, temp, x); \
        printf("["); for (int k = 0; k < B / (8 * sizeof(t)); ++k) { printf(f ",", temp[k]); } printf("]\n"); \
    } while(0)
#define Pw(x) PX(x, int32_t, storeu_epi32, "%i")
#define Pi(x) PX(x, int64_t, storeu_epi64, "%zi")
#define Ps(x) PX(x, float, storeu_ps, "%e")
#define Pd(x) PX(x, double, storeu_pd, "%e")
#endif
