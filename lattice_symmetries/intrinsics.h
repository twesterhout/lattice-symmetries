#pragma once
 
// M
// 1: AVX512    | 512
// 2: AVX2+FMA  | 256
// 3: SSE       | 
//    NEON      | 128
//    WASM      |
#if !defined(M) // Architecture
#   if defined(__AVX512BW__)
#       define M 1
#   elif defined(__AVX2__)
#       define M 2
#   else
#       define M 3
#   endif
#endif

#if M == 1 || M == 3
#include <immintrin.h>
#include <simde/x86/avx512.h>
#endif

typedef float float32_t; // Need it for some weird reason at M == 3 ...

#if M == 1
#  define A(f1,f2,f3) f1
#elif M == 2
#  define A(f1,f2,f3) f2
#else
#  define A(f1,f2,f3) f3
#endif

#if defined(__clang__)
// Gather f32
#if !__has_builtin(__builtin_ia32_gatherdiv4sf256)
#  define __builtin_ia32_gatherdiv4sf256 __builtin_ia32_gatherq_ps256
#endif
// Gather f64
#if !__has_builtin(__builtin_ia32_gatherdiv4df)
#  define __builtin_ia32_gatherdiv4df __builtin_ia32_gatherq_pd256
#endif
// Gather i32
#if !__has_builtin(__builtin_ia32_gatherdiv4si256)
#  define __builtin_ia32_gatherdiv4si256 __builtin_ia32_gatherq_d256
#endif
// Gather i64
#if !__has_builtin(__builtin_ia32_gatherdiv4d)
#  define __builtin_ia32_gatherdiv4di __builtin_ia32_gatherq_q256
#endif
#endif


// Debugging
#define PX(x, t, w, f) \
    do { \
        t temp[B / (8 * sizeof(t))]; I(w,temp,x); \
        printf("["); for (u64 k = 0; k < B / (8 * sizeof(t)); ++k) { printf(f ",", temp[k]); } printf("]\n"); \
    } while(0)
#if M == 1 || M == 2
#define Pw(x) PX(x, int32_t, storeu_epi32, "%i")
#define Pi(x) PX(x, int64_t, storeu_epi64, "%zi")
#define Ps(x) PX(x, float, storeu_ps, "%e")
#define Pd(x) PX(x, double, storeu_pd, "%e")
#else
#define Pw(x) PX(x, int32_t, storeu_si128, "%i")
#define Pi(x) PX(x, int64_t, storeu_si128, "%zi")
#define Ps(x) PX(x, float, storeu_ps, "%e")
#define Pd(x) PX(x, double, storeu_pd, "%e")
#endif

typedef char i8; typedef short i16; typedef int i32; typedef long long i64;
typedef unsigned char u8; typedef unsigned short u16; typedef unsigned u32; typedef unsigned long long u64;
typedef float f32; typedef double f64; typedef float _Complex c64; typedef double _Complex c128;

#define _(z...) ({z;})
#define D(t,g,k,x...) __attribute__((__always_inline__)) static inline t g(x){return _(k);}
#define Di(t,g,k,x...) static t g(x){return _(k);}
#define De(t,g,k,x...) t g(x){return _(k);}
#define Z2(x...) (Vz){x}
#define c(t) t const
#define _L(v,n,x...) for(i32 v=0;v<(n);++v){x;}
#define VL(v,n,x...) _Pragma("omp simd") for(i32 v=0;v<(n);++v){x;}
#define U4(x...) _L(u,4,x)
#define U8(x...) _L(u,8,x)
#define I__(p,f,x...) simde_##p##_##f(x)
#define I_(p,f,x...) I__(p,f,x)
#define I2(f,x...) I_(mm,f,x)
#define I4(f,x...) I_(mm256,f,x)
#define I8(f,x...) I_(mm512,f,x)
#define O_(f) __builtin_ia32_##f
#define O(f) O_(f)

#if M == 1
SIMDE_FUNCTION_ATTRIBUTES simde__m512d simde_mm512_cvtepi32_pd (simde__m256i a) {
  #if defined(SIMDE_X86_AVX512F_NATIVE)
    return _mm512_cvtepi32_pd(a);
  #else
    simde__m512d_private r_;
    simde__m256i_private a_ = simde__m256i_to_private(a);
    r_.m256d[0] = simde_mm256_cvtepi32_pd(a_.m128i[0]);
    r_.m256d[1] = simde_mm256_cvtepi32_pd(a_.m128i[1]);
    return simde__m512d_from_private(r_);
  #endif
}
SIMDE_FUNCTION_ATTRIBUTES simde__m256i simde_mm512_cvtepi32_epi16 (simde__m512i a) {
  #if defined(SIMDE_X86_AVX512F_NATIVE)
    return _mm512_cvtepi32_epi16(a);
  #else
    simde__m512i_private const a_ = simde__m512i_to_private(a);
    simde__m256i const low = a_.m256i[0], high = a_.m256i[1];
    simde__m256i const mask  = simde_mm256_set1_epi32(0x0000FFFF);         // mask for low words
    simde__m256i const lowm  = simde_mm256_and_si256(low, mask);           // words of low
    simde__m256i const highm = simde_mm256_and_si256(high, mask);          // words of high
    simde__m256i const pk    = simde_mm256_packus_epi32(lowm,highm);       // unsigned pack
    return simde_mm256_permute4x64_epi64(pk, 0xD8);                        // put in right place
  #endif
}
SIMDE_FUNCTION_ATTRIBUTES simde__m512d simde_mm512_cvtps_pd (simde__m256 a) {
  #if defined(SIMDE_X86_AVX512F_NATIVE)
    return _mm512_cvtps_pd(a);
  #else
    simde__m512d_private r_;
    simde__m256_private const a_ = simde__m256_to_private(a);
    r_.m256d[0] = simde_mm256_cvtps_pd(a_.m128[0]);
    r_.m256d[1] = simde_mm256_cvtps_pd(a_.m128[1]);
    return simde__m512d_from_private(r_);
  #endif
}
SIMDE_FUNCTION_ATTRIBUTES simde__m256 simde_mm512_cvtpd_ps (simde__m512d a) {
  #if defined(SIMDE_X86_AVX512F_NATIVE)
    return _mm512_cvtpd_ps(a);
  #else
    simde__m512d_private const a_ = simde__m512d_to_private(a);
    simde__m128 a0 = simde_mm256_cvtpd_ps(a_.m256d[0]);
    simde__m128 a1 = simde_mm256_cvtpd_ps(a_.m256d[1]);
    return simde_mm256_insertf128_ps(simde_mm256_castps128_ps256(a0),a1,1);
  #endif
}
#endif

#define R2(x) x,x
#define R4(x) x,x,x,x
#define R8(x) x,x,x,x,x,x,x,x

#define RF2(f) f(0),f(1)
#define RF4(f) f(0),f(1),f(2),f(3)
#define RF8(f) f(0),f(1),f(2),f(3),f(4),f(5),f(6),f(7)

#if M == 1
#   define B 512
#   define _S(f) f##_si512
#   define I(f,x...) I8(f,x)
#   define Zi Z8i
#   define Zd Z8d
#   define Zf Z16f
#   define RN(x) R8(x)
#   define RFN(f) RF8(f)
#elif M == 2
#   define B 256
#   define _S(f) f##_si256
#   define I(f,x...) I4(f,x)
#   define Zi Z4i
#   define Zd Z4d
#   define Zf Z8f
#   define RN(x) R4(x)
#   define RFN(f) RF4(f)
#else
#   define B 128
#   define _S(f) f##_si128
#   define I(f,x...) I2(f,x)
#   define Zi Z2i
#   define Zd Z2d
#   define Zf Z4f
#   define RN(x) R2(x)
#   define RFN(f) RF2(f)
#endif
#define N (B / 64)
#define Va(n) __attribute__((vector_size(n),aligned(n)))
#define Vu(n) __attribute__((vector_size(n),aligned(1)))
#if defined(__clang__)
#  define shuffle1(x...) __builtin_shufflevector(x)
#else
#  define shuffle1(x...) __builtin_shuffle(x)
#endif
#define shuffle2(x...) __builtin_shufflevector(x)

#define Z8i ((V8q){0,0,0,0,0,0,0,0})
#define Z8d ((V8d){0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0})
#define Z16f ((V16f){0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f})
#define Z4i ((V4q){0,0,0,0})
#define Z4d ((V4d){0.0,0.0,0.0,0.0})
#define Z8f ((V8f){0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f})
#define Z2i ((V2q){0,0})
#define Z2d ((V2d){0.0,0.0})
#define Z4f ((V4f){0.0f,0.0f,0.0f,0.0f})

// unaligned SIMD vector types; use these only for reading and writing
typedef i64 _Vq Vu(B/8);typedef f64 _Vd Vu(B/8);typedef f32 _Vf Vu(B/8);
#define _R(t,p,i) ((c(t)*)(p))[i]
#define _W(t,p,i,z) ((t*)(p))[i]=(z)
#define Rq(p,i) _R(_Vq,p,i)
#define Wq(p,i,z) _W(_Vq,p,i,z)
#define Rd(p,i) _R(_Vd,p,i)
#define Wd(p,i,z) _W(_Vd,p,i,z)
#define Rf(p,i) _R(_Vf,p,i)
#define Wf(p,i,z) _W(_Vf,p,i,z)
// aligned SIMD vector types; use these for normal calculations & local variables
typedef i64 V8q Va(64);typedef i64 V4q Va(32); typedef i64 V2q Va(16);
typedef i32 V16i Va(64);typedef i32 V8i Va(32); typedef i32 V4i Va(16);
typedef u16 V32w Va(64);typedef u16 V16w Va(32); typedef u16 V8w Va(16);
typedef f64 V8d Va(64);typedef f64 V4d Va(32); typedef f64 V2d Va(16);
typedef f32 V16f Va(64);typedef f32 V8f Va(32); typedef f32 V4f Va(16);
// M-dependent aliases
typedef i64 Vq Va(B/8);typedef u64 Vuq Va(B/8);typedef i32 Vi Va(B/8);typedef u16 Vw Va(B/8);typedef char Vb Va(B/8);typedef f64 Vd Va(B/8);typedef f32 Vf Va(B/8);
// complex number
typedef struct Vz{Vd re;Vd im;}Vz;
// Mask type that's returned from various comparison operators and that we ues for masked gather
typedef A(simde__mmask8,Vq,Vq) M8; // AVX512 uses bitmasks; everything else just uses vectors

#if M == 1
    D(V4d,V8d_hi,I8(extractf64x4_pd,x,1),c(V8d)x)
    D(V8f,V16f_hi,I8(extractf32x8_ps,x,1),c(V16f)x)
    D(V4d,V8d_lo,I8(castpd512_pd256,x),c(V8d)x)
    D(V8f,V16f_lo,I8(castps512_ps256,x),c(V16f)x)
    D(V8d,V8f_d,I8(cvtps_pd,x),c(V8f)x)
    D(V8f,V8d_f,I8(cvtpd_ps,x),c(V8d)x)
#endif
#if M <= 2
    D(V2d,V4d_hi,(O(vextractf128_pd256)(x,1)),c(V4d)x)
    D(V4f,V8f_hi,(O(vextractf128_ps256)(x,1)),c(V8f)x)
    D(V2d,V4d_lo,shuffle2(x,x,0,1),c(V4d)x)
    D(V4f,V8f_lo,shuffle2(x,x,0,1,2,3),c(V8f)x)
    D(V4f,V4d_f,(__builtin_convertvector(x,V4f)),c(V4d)x)
    D(V4d,V4f_d,(__builtin_convertvector(x,V4d)),c(V4f)x)
    D(V4f,V2d_f,O(cvtpd2ps)(x),c(V2d)x)
#endif
#if M == 3
    D(V2d,V4f_d,((Vd){(f64)x[0],(f64)x[1]}),c(V4f)x)
    D(V4f,V2d_f,((Vf){(f32)x[0],(f32)x[1],0.0f,0.0f}),c(V2d)x)
#endif

D(Vq,Si,((Vq){RN(x)}),c(i64)x)
D(Vd,Sd,((Vd){RN(x)}),c(f64)x)
D(u32,movemask,A(m,O(movmskpd256)((Vd)m),O(movmskpd)((Vd)m)),c(M8)m)

#if M == 1
#   define hi(X) _Generic((X), V8d: V8d_hi, V4d: V4d_hi, V16f: V16f_hi, V8f: V8f_hi)(X)
#   define lo(X) _Generic((X), V8d: V8d_lo, V4d: V4d_lo, V16f: V16f_lo, V8f: V8f_lo)(X)
#   define d2f(X) _Generic((X), V8d: V8d_f, V4d: V4d_f, V2d: V2d_f)(X)
#   define f2d(X) _Generic((X), V8f: V8f_d, V4f: V4f_d)(X)
#   define select(s,a,b) I(mask_blend_epi64,s,b,a)
    // D(u32,movemask,m,c(M8)m)
#elif M == 2
#   define hi(X) _Generic((X), V4d: V4d_hi, V8f: V8f_hi)(X)
#   define lo(X) _Generic((X), V4d: V4d_lo, V8f: V8f_lo)(X)
#   define d2f(X) _Generic((X), V4d: V4d_f, V2d: V2d_f)(X)
#   define f2d(X) _Generic((X), V4f: V4f_d)(X)
#   define select(s,a,b) (Vq)O(pblendvb256)((Vb)(b),(Vb)(a),(Vb)(s))
    // D(u32,movemask,O(movmskpd256)((Vd)m),c(M8)m)
#else
#   define d2f(X) _Generic((X), V2d: V2d_f)(X)
#   define f2d(X) _Generic((X), V4f: V4f_d)(X)
#   define select(s,a,b) I(blendv_epi8,b,a,s)
    // D(u32,movemask,I(movemask_pd,(Vd)(m)),c(M8)m)
#endif

#if M == 1
#   define gt(a,b) I(cmp_epi64_mask,b,a,1)
#   define eq(a,b) I(cmp_epi64_mask,a,b,0)
#else
#   define gt(a,b) ((a)>(b))
#   define eq(a,b) ((a)==(b))
#endif

#if M == 2
#define shl(x,imm) 

#endif

#if M == 1
#  define EVEN 0,2,4,6,8,10,12,14
#  define ODD 1,3,5,7,9,11,13,15
#  define Rx4 c(Vf)b=shuffle1(a,(Vi){EVEN,ODD});r=Z2(f2d(lo(b)),f2d(hi(b)))
#  define Wx3_idx0 0,8,1,9,2,10,3,11
#  define Wx3_idx1 4,12,5,13,6,14,7,15
#  define Wx4_idx Wx3_idx0,Wx3_idx1
#  define Vhf simde__m256
#  define Vhi simde__m256i
#  define Vhw V16w
#  define Wq2w_idx0 0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30
#  define Wq2w_idx1 32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62
#  define Gmw_idx0 0,17,2,19,4,21,6,23,8,25,10,27,12,29,14,31
#  define Gmw_zero (V16w){0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
#elif M == 2
#  define EVEN 0,2,4,6
#  define ODD 1,3,5,7
#  define Rx4 c(Vf)b=shuffle1(a,(Vi){EVEN,ODD});r=Z2(f2d(lo(b)),f2d(hi(b)))
#  define Wx3_idx0 0,4,1,5
#  define Wx3_idx1 2,6,3,7
#  define Wx4_idx Wx3_idx0,Wx3_idx1
#  define Vhf simde__m128
#  define Vhi simde__m128i
#  define Vhw V8w
#  define Wq2w_idx0 0,2,4,6,8,10,12,14
#  define Wq2w_idx1 16,18,20,22,24,26,28,30
#  define Gmw_idx0 0,9,2,11,4,13,6,15
#  define Gmw_zero (V8w){0,0,0,0,0,0,0,0}
#else
#  define EVEN 0,2
#  define ODD 1,3
#  define Rx4 r=Z2((Vd){(f64)a[0],(f64)a[2]},(Vd){(f64)a[1],(f64)a[3]})
#  define Wx3_idx0 0,2
#  define Wx3_idx1 1,3
#  define Wx4_idx 0,4,1,5
#  define Vhi simde__m128i
#  define Vhw V8w
#  define Wq2w_idx0 0,2,4,6
#  define Wq2w_idx1 8,10,12,14
#  define Gmw_idx0 0,9,2,11,-1,-1,-1,-1
#  define Gmw_zero (V8w){0,0,0,0,0,0,0,0}
#endif

#if M == 1
#  define Gq(p,i) I(i64gather_epi64,i,p,8)
#  define _gthd(p,i,s) I(i64gather_pd,i,p,s)
#  define _gthf(p,i,s) I(i64gather_ps,i,p,s)
#  define _mgthd(p,i,m,s) I(mask_i64gather_pd,Zd,m,i,p,s)
#  define _mgthf(p,i,m,s) I(mask_i64gather_ps,I4(setzero_ps),m,i,p,s)
#  define _mgthi(p,i,m,s) I(mask_i64gather_epi32,I4(setzero_si256),m,i,p,s)
#elif M == 2
#  define Gq(p,i) O(gatherdiv4di)(Zi,(c(i64)*)(p),i,((Vq){~0,~0,~0,~0}),8)
#  define maski64_i32(m) shuffle2((Vi)m,(Vi)m,0,2,4,6)
#  define _mgthd(p,i,m,s) O(gatherdiv4df)(Zd,p,i,(Vd)(m),s)
#  define _mgthf(p,i,m,s) O(gatherdiv4sf256)(Z4f,p,i,(V4f)(maski64_i32(m)),s)
#  define _mgthi(p,i,m,s) O(gatherdiv4si256)(((V4i){0,0,0,0}),p,i,(V4i)(maski64_i32(m)),s)
#  define _gthd(p,i,s) _mgthd(p,i,Zd==Zd,s)
#  define _gthf(p,i,s) O(gatherdiv4sf256)(Z4f,p,i,Z4f==Z4f,s)
#else
#   define maski64_i32(m) shuffle2((Vi)m,(Vi)m,0,2,4,6)
    D(Vq,Gq,((Vq){p[i[0]],p[i[1]]}),c(u64)*p,c(Vq)i)
    D(Vd,_gthd,_(c(u8)*p0=(c(u8)*)p+i[0]*s,*p1=(c(u8)*)p+i[1]*s;(Vd){*(c(f64)*)p0,*(c(f64)*)p1}),c(f64)*p,c(Vq)i,c(i32)s)
    D(Vf,_gthf,_(c(u8)*p0=(c(u8)*)p+i[0]*s,*p1=(c(u8)*)p+i[1]*s;(Vf){*(c(f32)*)p0,*(c(f32)*)p1,0.0f,0.0f}),c(f32)*p,c(Vq)i,c(i32)s)
    D(Vd,_mgthd,_(c(u8)*p0=(c(u8)*)p+i[0]*s,*p1=(c(u8)*)p+i[1]*s;(Vd){m[0]?*(c(f64)*)p0:0.0,m[1]?*(c(f64)*)p1:0.0}),c(f64)*p,c(Vq)i,c(M8)m,c(i32)s)
    D(Vf,_mgthf,_(c(u8)*p0=(c(u8)*)p+i[0]*s,*p1=(c(u8)*)p+i[1]*s;(Vf){m[0]?*(c(f32)*)p0:0.0f,m[1]?*(c(f32)*)p1:0.0f,0.0f,0.0f}),c(f32)*p,c(Vq)i,c(M8)m,c(i32)s)
    D(Vi,_mgthi,_(c(u8)*p0=(c(u8)*)p+i[0]*s,*p1=(c(u8)*)p+i[1]*s;(Vi){m[0]?*(c(i32)*)p0:0,m[1]?*(c(i32)*)p1:0,0,0}),c(i32)*p,c(Vq)i,c(M8)m,c(i32)s)
#endif

// type code (i32)
// 0: f64
// 1: f32
// 2: f16
// 3: c128
// 4: c64
// 5: c32

// Rx---read a vector of numbers from p and convert it to Vz. The type is specified by t.
D(Vz,Rx,_(Vz r;switch(t){
    /*f64*/case 0: r=Z2(Rd(p,0),Zd);break;
    /*f32*/case 1: {c(f32)*_p=p;VL(u,N,r.re[u]=(f64)_p[u]);r.im=Zd;break;}
    /*f16*/case 2: r=Z2(Zd,Zd);break;
   /*c128*/case 3: {c(Vd)a=Rd(p,0),b=Rd(p,1);r=Z2(shuffle2(a,b,EVEN),shuffle2(a,b,ODD));break;}
    /*c64*/case 4: {c(Vf)a=Rf(p,0);Rx4;break;}
   /*c32*/default: r=Z2(Zd,Zd);break;
};r),c(i32)t,c(void)*p)

// Wx---convert Vz to type t and write to p.
D(void,Wx,_(switch(t){
    /*f64*/case 0: {Wd(p,0,z.re);break;}
    /*f32*/case 1: {f32*_p=p;VL(u,N,_p[u]=(f32)z.re[u]);break;}
    /*f16*/case 2: break;
   /*c128*/case 3: {Wd(p,0,shuffle2(z.re,z.im,Wx3_idx0));Wd(p,1,shuffle2(z.re,z.im,Wx3_idx1));break;}
    /*c64*/case 4: {Wf(p,0,shuffle2(d2f(z.re),d2f(z.im),Wx4_idx));break;}
   /*c32*/default: break;
}),c(i32)t,void*p,c(Vz)z)

D(Vz,Gmx,_(Vz r;switch(t){
    /*f64*/case 0: r=Z2(_mgthd((c(f64)*)p,i,m,8),Zd);break;
    /*f32*/case 1: r=Z2(f2d(_mgthf((c(f32)*)p,i,m,4)),Zd);break;
    /*f16*/case 2: r=Z2(Zd,Zd);break;
   /*c128*/case 3: r=Z2(_mgthd((c(f64)*)p,i<<1,m,8),_mgthd((c(f64)*)p+1,i<<1,m,8));break;
    /*c64*/case 4: {c(Vf)a=(Vf)_mgthd((c(f64)*)p,i,m,8);Rx4;break;}
   /*c32*/default: r=Z2(Zd,Zd);break;
};r),c(i32)t,c(void)*p,c(Vq)i,c(M8)m)

D(Vz,Gx,_(Vz r;switch(t){
     /*f64*/case 0: r=Z2(_gthd((c(f64)*)p,i,8),Zd);break;
     /*f32*/case 1: r=Z2(f2d(_gthf((c(f32)*)p,i,4)),Zd);break;
     /*f16*/case 2: r=Z2(Zd,Zd);break;
    /*c128*/case 3: r=Z2(_gthd((c(f64)*)p,i<<1,8),_gthd((c(f64)*)p+1,i<<1,8));break;
     /*c64*/case 4: {c(Vf)a=(Vf)_gthd((c(f64)*)p,i,8);Rx4;break;}
    /*c32*/default: r=Z2(Zd,Zd);break;
};r),c(i32)t,c(void)*p,c(Vq)i)

// Rw2d---read N u16 numbers from p and convert them to f64, returning Vd
D(Vd,Rw2d,_(Vd r;VL(u,N,r[u]=p[u]);r),c(u16)*p)

// Wq2w---convert N i64 numbers (i.e., Vq) to u16 and write them to p. Receives 4 Vq as input and writes 1 Vw to p.
D(void,Wq2w,_(*(Vw*)p=shuffle2((Vw)shuffle2((Vi)x[0],(Vi)x[1],Wq2w_idx0),(Vw)shuffle2((Vi)x[2],(Vi)x[3],Wq2w_idx0),Wq2w_idx0,Wq2w_idx1);return),u16*p,c(Vq)x[static 4])

#if M == 1 || M == 3
#  define epi32_pd(x) I(cvtepi32_pd,x)
#else
#  define epi32_pd(x) __builtin_convertvector(x,Vd)
#endif

D(Vd,Gmw,_(
    __auto_type x=(A(Vhi,V4i,Vhi))shuffle2((Vhw)_mgthi((c(i32)*)p,i,m,2),Gmw_zero,Gmw_idx0);
    A(I(cvtepi32_pd,x),epi32_pd(x),I(cvtepi32_pd,x))
),c(u16)*p,c(Vq)i,c(M8)m)

#define _F_popcnt(i) __builtin_popcountll(x[i])
D(Vq,popcnt,_((Vq){RFN(_F_popcnt)}),c(Vq)x)
#undef _F_popcnt

#if M == 1 || M == 3
    D(Vd,sqrtd,I(sqrt_pd,x),c(Vd)x)
#else
#   define sqrtd(x) O(sqrtpd256)(x)
#endif

D(Vz,mulz,Z2(a.re*b.re-a.im*b.im,a.im*b.re+a.re*b.im),c(Vz)a,c(Vz)b)
#if M == 1
    D(Vq,m1,I(movm_epi64,I(test_epi64_mask,x,m)),c(Vq)x,c(Vq)m)
#else
    D(Vq,m1,(x&m)!=(Vq)Zi,c(Vq)x,c(Vq)m)
#endif
D(Vq,m2,m1(x,m_0)^m1(x,m_1),c(Vq)x,c(Vq)m_0,c(Vq)m_1)
D(Vq,mX,popcnt(x&m),c(Vq)x,c(Vq)m)
D(Vz,bcast2,(Z2(Sd(re[k]),Sd(im[k]))),c(f64)*re,c(f64)*im,c(i32)k)
D(Vz,flipsign,_(m=m<<63;Z2((Vd)((Vq)v.re^m),(Vd)((Vq)v.im^m))),c(Vz)v,Vq m)
// TODO: The loop is ugly, but the compilers vectorize it
#define prefetch(p...) __builtin_prefetch(p)
D(void,prefetchq,_(_L(k,N,prefetch(p+idx[k],0,3))),c(u64)*p,c(Vq)idx)
// TODO: Should get rid of these ...
// static void prefetchd4xN(Vi idx[4],c(u16)*n,c(f64)*x){_L(k,4*N,c(i32)i=((c(i64)*)idx)[k];simde_mm_prefetch(n+i,_MM_HINT_ET0);simde_mm_prefetch(x+i,_MM_HINT_ET0))}
// static void prefetchz4xN(Vi idx[4],c(u16)*n,c(c128)*x){_L(k,4*N,c(i32)i=((c(i64)*)idx)[k];simde_mm_prefetch(n+i,_MM_HINT_ET0);simde_mm_prefetch(x+i,_MM_HINT_ET0))}
