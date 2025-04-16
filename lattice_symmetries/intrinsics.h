#pragma once

#include <simde/x86/avx512.h>
#include <simde/x86/f16c.h>

// M
// 1: AVX512    | 512
// 2: AVX2+FMA  | 256
// 3: SSE       | 
//    NEON      | 128
//    WASM      |
#if !defined(M) // Architecture
#   if defined(SIMDE_X86_AVX512F_NATIVE)
#       define M 1
#   elif defined(SIMDE_X86_AVX2_NATIVE)
#       define M 2
#   else
#       define M 3
#   endif
#endif

#if !defined(USE_F16)
#  define USE_F16 SIMDE_FLOAT16_API == SIMDE_FLOAT16_API_FLOAT16
#endif

// type code (i32)
// 0: f64
// 1: f32
// 2: f16
// 3: c128
// 4: c64
// 5: c32

#define INTERNAL static HEDLEY_ALWAYS_INLINE

typedef int8_t i8; typedef int16_t i16; typedef int32_t i32; typedef int64_t i64;
typedef uint8_t u8; typedef uint16_t u16; typedef uint32_t u32; typedef uint64_t u64;
typedef float f32; typedef double f64;
typedef float _Complex c64; typedef double _Complex c128;
#if USE_F16
typedef _Float16 f16; typedef _Float16 _Complex c32;
#endif

#define _(z...) ({z;})
#define $(b,z) if(b){z;}else
#define D(t,g,k,x...) static HEDLEY_ALWAYS_INLINE t g(x){return _(k);}
#define Di(t,g,k,x...) static t g(x){return _(k);}
#define De(t,g,k,x...) t g(x){return _(k);}
#define Dd(g,k,x...) D(Vd,g##d,k,x)
#define Dz(g,k,x...) D(Vz,g##z,k,x)
#define D2_(t,f,i) D(t,f,I(i,a,b),c(t)a,c(t)b)
#define D2(t,f,i) D2_(t,f,i)
#define Z2(x,y) (Vz){x,y}
#define f64c f64 const
#define c(t) t const
// For loop over a counter v of length n
#define _L(v,n,x...) for(i32 v=0;v<(n);++v){x;}
#define U4(x...) _L(u,4,x)
#define U8(x...) _L(u,8,x)

// always operate on B bits of data
#define I__(p,f,x...) simde_##p##_##f(x)
#define I_(p,f,x...) I__(p,f,x)
#define I2(f,x...) I_(mm,f,x)
#define I4(f,x...) I_(mm256,f,x)
#define I8(f,x...) I_(mm512,f,x)
// #define I_impl(b,f,x...) I_(mm##b,f,x)
// #define I(b,f,x...) I_impl(b,f,x)

#if M == 1
#   define B 512
#   define _S(f) f##_si512
#   define I(f,x...) I8(f,x)
#   define Ih(f,x...) I4(f,x)
#   define i2d(x) I(castsi512_pd,x)
#   define i2s(x) I(castsi512_ps,x)
    typedef simde__m512i Vi; typedef simde__m512d Vd; typedef simde__m512 Vs; typedef simde__mmask8 M8;
#elif M == 2
#   define B 256
#   define _S(f) f##_si256
#   define I(f,x...) I4(f,x)
#   define Ih(f,x...) I2(f,x)
#   define i2d(x) I(castsi256_pd,x)
#   define i2s(x) I(castsi256_ps,x)
    typedef simde__m256i Vi; typedef simde__m256d Vd; typedef simde__m256 Vs; typedef Vi M8;
#else
#   define B 128
#   define _S(f) f##_si128
#   define I(f,x...) I2(f,x)
#   define i2d(x) I(castsi128_pd,x)
#   define i2s(x) I(castsi128_ps,x)
    typedef simde__m128i Vi; typedef simde__m128d Vd; typedef simde__m128 Vs; typedef Vi M8;
#endif
#define d2f(x) I(castpd_ps,x)
#define f2d(x) I(castps_pd,x)
D(Vi,d2i,I(_S(castpd),x),c(Vd)x)D(Vi,s2i,I(_S(castps),x),c(Vs)x)

#if M == 1
#   define eqi(a,b) I(movm_epi64,eqq(a,b))
#   define gt(a,b) I(cmp_epi64_mask,b,a,1)
#   define select(s,a,b) I(mask_blend_epi64,s,b,a)
    D(M8,eqq,I(cmp_epi64_mask,a,b,0),c(Vi)a,c(Vi)b)
    D(Vi,Si,I(set1_epi64,x),c(i64)x)
#else
#   define eqi(a,b) eqq(a,b)
#   define gt(a,b) I(cmpgt_epi64,a,b)
#   define select(s,a,b) I(blendv_epi8,b,a,s)
    D(M8,eqq,I(cmpeq_epi64,a,b),c(Vi)a,c(Vi)b)
    D(Vi,Si,I(set1_epi64x,x),c(i64)x)
#endif

#define N (B / 64)
#define Zi I(_S(setzero))
#define Zd I(setzero_pd)
#define Zf I(setzero_ps)
#define shl(x,imm) I(slli_epi64,x,imm)
#define shr(x,imm) I(srli_epi64,x,imm)

#if M == 1 // Conversion functions that are missing from SIMDe
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

#   if USE_F16
    SIMDE_FUNCTION_ATTRIBUTES simde__m512d simde_mm512_cvtph_pd (simde__m128i a) {
      #if defined(SIMDE_X86_AVX512FP16_NATIVE)
        return _mm512_cvtph_pd((__m128h)a);
      #else
        simde__m256d a0 = simde_mm256_cvtps_pd(simde_mm_cvtph_ps(a));
        simde__m256d a1 = simde_mm256_cvtps_pd(simde_mm_cvtph_ps(simde_mm_cvtsi64_si128(simde_mm_extract_epi64(a, 1))));
        return simde_mm512_insertf64x4(simde_mm512_castpd256_pd512(a0), a1, 1);
      #endif
    }
    SIMDE_FUNCTION_ATTRIBUTES simde__m128i simde_mm512_cvtpd_ph (simde__m512d a) {
      #if defined(SIMDE_X86_AVX512FP16_NATIVE)
        return (__m128i)_mm512_cvtpd_ph(a);
      #else
        simde__m512d_private const a_ = simde__m512d_to_private(a);
        simde__m128i a0 = simde_mm_cvtps_ph(simde_mm256_cvtpd_ps(a_.m256d[0]), SIMDE_MM_FROUND_TO_NEAREST_INT);
        simde__m128i a1 = simde_mm_cvtps_ph(simde_mm256_cvtpd_ps(a_.m256d[1]), SIMDE_MM_FROUND_TO_NEAREST_INT);
        return simde_mm_insert_epi64(a0, simde_mm_cvtsi128_si64(a1), 1);
      #endif
    }
#   endif
#endif

// Vi
D(Vi,Ri,I(loadu_epi64,p),c(void)*p)D(void,Wi,I(_S(storeu),p,x),void*p,c(Vi)x)D2(Vi,and,_S(and))D2(Vi,xor,_S(xor))D2(Vi,andnot,_S(andnot))D2(Vi,or,_S(or))D2(Vi,addq,add_epi64)D2(Vi,subq,sub_epi64)D(Vi,popcnt,I(popcnt_epi64,x),c(Vi)x)
// Vd
Dd(zero,Zd,void)D(Vd,Rd,I(loadu_pd,p),c(f64)*p)D(Vd,Sd,I(set1_pd,x),c(f64)x)D(void,Wd,I(storeu_pd,p,x),f64*p,c(Vd)x)D2(Vd,addd,add_pd)D2(Vd,subd,sub_pd)D2(Vd,muld,mul_pd)D2(Vd,divd,div_pd)D(Vd,fmad,I(fmadd_pd,a,b,c),c(Vd)a,c(Vd)b,c(Vd)c)D(Vd,fnmad,I(fnmadd_pd,a,b,c),c(Vd)a,c(Vd)b,c(Vd)c)D(Vd,sqrtd,I(sqrt_pd,x),c(Vd)x)
// Vz
typedef struct Vz { Vd re; Vd im; } Vz;
Dz(zero,(Z2(Zd,Zd)),void)D(Vz,addz,Z2(addd(a.re, b.re),addd(a.im, b.im)),c(Vz)a,c(Vz)b)D(Vz,mulz,Z2(fnmad(a.im,b.im,muld(a.re,b.re)),fmad(a.im,b.re,muld(a.re,b.im))),c(Vz)a,c(Vz)b)
#if M == 1
    // D(Vz,Rz,_(
    //     c(Vd)a=Rd((c(f64)*)p),b=Rd((c(f64)*)p+N);_Alignas(Vi)c(i64)idx[N]={0,2,4,6,1,3,5,7};
    //     c(Vd)c=I(permutexvar_pd,I(load_epi64,idx),I(unpacklo_pd,a,b)),d=I(permutexvar_pd,I(load_epi64,idx),I(unpackhi_pd,a,b));
    //     Z2(c,d)), c(c128)*p)
    // TODO: boy is this one ugly... is there a better way?
    // De(void,Wz,_(
    //     _Alignas(Vi)c(i64)idx1[N]={0b0000,0b0001,0b1000,0b1001,0b0010,0b0011,0b1010,0b1011};
    //     _Alignas(Vi)c(i64)idx2[N]={0b0100,0b0101,0b1100,0b1101,0b0110,0b0111,0b1110,0b1111};
    //     c(Vd)a=I(unpacklo_pd,z.re,z.im),b=I(unpackhi_pd,z.re,z.im);
    //     c(Vd)c=I(permutex2var_pd,a,I(load_epi64,idx1),b),d=I(permutex2var_pd,a,I(load_epi64,idx2),b);
    //     Wd((f64*)p,c);Wd((f64*)p+N,d)),c(c128)*p,c(Vz)z)
#elif M == 2
    // D(Vz,Rz,_(c(Vd)a=Rd((c(f64)*)p),b=Rd((c(f64)*)p+N),c=I(permute2f128_pd,a,b,0b00110001),d=I(permute2f128_pd,a,b,0b00100000);Z2(I(unpacklo_pd,d,c),I(unpackhi_pd,d,c))),c(c128)*p)
    // D(void,Wz,_(c(Vd)a=I4(unpacklo_pd,z.re,z.im),b=I4(unpackhi_pd,z.re,z.im),c=I(permute2f128_pd,a,b,0b00110001),d=I(permute2f128_pd,a,b,0b00100000);Wd((f64*)p,d),Wd((f64*)p+N,c)),c128*p,c(Vz)z)
#else
    // D(Vz,Rz,_(c(Vd)a=Rd((c(f64)*)p),b=Rd((c(f64)*)p+N);Z2(I(unpacklo_pd,a,b),I(unpackhi_pd,a,b))),c(c128)*p)
    // D(void,Wz,_(Wd((f64*)p,I(unpacklo_pd,z.re,z.im)),Wd((f64*)p+N,I(unpackhi_pd,z.re,z.im))),c128*p,c(Vz)z)
#endif
// u16
#if M == 1
    // Why is SIMDe missing these functions...?
    D(Vd,Rw2d,I(cvtepi32_pd,I4(cvtepu16_epi32,I2(loadu_si128,(c(void)*)p))),c(u16)*p)
    D(void,Wq2w,_(
        c(Vi) y0 = I(inserti64x4,I(castsi256_si512,I(cvtepi64_epi32,x[0])),I(cvtepi64_epi32,x[1]),1);
        c(Vi) y1 = I(inserti64x4,I(castsi256_si512,I(cvtepi64_epi32,x[2])),I(cvtepi64_epi32,x[3]),1);
        c(Vi) z = I(inserti64x4,I(castsi256_si512, I(cvtepi32_epi16, y0)),I(cvtepi32_epi16, y1), 1);
        Wi((void*)p,z)),u16*p,c(Vi)x[static 4])
#elif M == 2
    D(Vd,Rw2d,I4(cvtepi32_pd,I2(cvtepu16_epi32,I2(loadu_epi64,(c(u64)*)p))),c(u16)*p)
    D(void,Wq2w,_(
        Vi y0=s2i(I(shuffle_ps,i2s(x[0]),i2s(x[1]),0b10001000)),y1=s2i(I(shuffle_ps,i2s(x[2]),i2s(x[3]),0b10001000));
        // TODO: For some reason GCC can't see through blend_epi16 and generates vpblenvb instead of vpblendw... Invoking _mm256_blend_epi16 directly seems to help
        // We only call Wq2w inside norm64, so it's probably not that important.
        y0=I(permute4x64_epi64,I(blend_epi16,y0,Zi,0b10101010),0b11011000);
        y1=I(permute4x64_epi64,I(blend_epi16,y1,Zi,0b10101010),0b11011000);
        Vi r=I(permute4x64_epi64,I(packus_epi32,y0,y1),0b11011000);
        Wi((i64*)p,r)),u16*p,c(Vi)x[static 4])
#else
    D(Vd,Rw2d,I2(set_pd,(double)p[1],(double)p[0]),c(u16)*p)
    D(void,Wq2w,_(c(Vi)y0=s2i(I(shuffle_ps,i2s(x[0]),i2s(x[1]),0b10001000)),y1=s2i(I(shuffle_ps,i2s(x[2]),i2s(x[3]),0b10001000));Wi((i64*)p,I(packus_epi32,I(blend_epi16,y0,Zi,0b10101010),I(blend_epi16,y1,Zi,0b10101010)))),u16*p,c(Vi)x[static 4])
#endif

// #if M == 1
// #define ps2pd(x) I(cvtps_pd,x)
// #define ph2pd(x) I(cvtph_pd,x)
// #define pd2ps(x) I(cvtpd_ps,x)
// #elif M == 2
// #define ps2pd(x) I(cvtps_pd,x)
// #define ph2pd(x) ps2pd(I(cvtph_ps,x))
// #define pd2ps(x) I(cvtpd_ps,x)
// #else
// // TODO: ps2pd, ph2pd
// #endif

// read
#define Rd8(p) I8(loadu_pd,p)
#define Rd4(p) I4(loadu_pd,p)
#define Rd2(p) I2(loadu_pd,p)
#define Rf16(p) I8(loadu_ps,p)
#define Rf8(p) I4(loadu_ps,p)
#define Rf4(p) I2(loadu_ps,p)
#define Rf2(p) I2(loadl_pi,I2(setzero_ps),(simde__m64 const*)(p))
#if USE_F16
#  define Rh8(p) I2(loadu_si128,p)
#  define Rh4(p) I2(loadu_si64,p)
#  define Rh2(p) I2(loadu_si32,p)
#endif
// write
#define Wd8(p,x) I8(storeu_pd,p,x)
#define Wd4(p,x) I4(storeu_pd,p,x)
#define Wd2(p,x) I2(storeu_pd,p,x)
#define Wf8(p,x) I4(storeu_ps,p,x)
#define Wf4(p,x) I2(storeu_ps,p,x)
#define Wf2(p,x) I2(storel_pi,(simde__m64*)(p),x)
#if USE_F16
#  define Wh8(p,x) I2(storeu_si128,p,x)
#  define Wh4(p,x) I2(storeu_si64,p,x)
#  define Wh2(p,x) I2(store_ss,(f32*)(p),(simde__m128)x)
#endif
// mask gather real

// #define mgthd8(p,i,m) I8(mask_i64gather_pd,I8(setzero_pd),m,i,p,8)
// #define mgthf8(p,i,m) I8(mask_i64gather_ps,I4(setzero_ps),m,i,p,4)
// #define mgthd4(p,i,m) I4(mask_i64gather_pd,I4(setzero_pd),p,i,I4(castsi256_pd,m),8)
// #define mgthf4(p,i,m) I4(mask_i64gather_ps,I2(setzero_ps),p,i,I2(castsi128_ps,m),4)
// #define mgthd2(p,i,m) I2(set_pd,I2(extract_epi64,m,1)?(p)[I2(extract_epi64,i,1)]:0.0,I2(cvtsi128_si64,m)?(p)[I2(cvtsi128_si64,i)]:0.0)
// #define mgthf2(p,i,m) I2(set_ps,0.0f,0.0f,I2(extract_epi64,m,1)?(p)[I2(extract_epi64,i,1)]:0.0f,I2(cvtsi128_si64,m)?(p)[I2(cvtsi128_si64,i)]:0.0f)
// TODO: FIXME: uses scalar code; re-write properly
// Di(simde__m128i,mgthh8,_(
//   _Alignas(simde__m512i) i64 _i[8];I8(store_epi64,_i,i);_Alignas(simde__m128i) _Float16 out[8];_L(u,8,out[u]=((m >> u)&0x1)?p[_i[u]]:(_Float16)0)
//   Rh8(out)),c(_Float16)*p,c(simde__m512i)i,c(simde__mmask8)m)
// Di(simde__m128i,mgthh4,_(
//   _Alignas(simde__m256i) i64 _i[4];I4(storeu_epi64,_i,i);_Alignas(simde__m256i) i64 _m[4];I4(storeu_epi64,_m,m);
//   _Alignas(simde__m128i) _Float16 out[4];_L(u,4,out[u]=_m[u]?p[_i[u]]:(_Float16)0)
//   Rh4(out)),c(_Float16)*p,c(simde__m256i)i,c(simde__m256i)m)
// Di(simde__m128i,mgthh2,_(
//   _Alignas(simde__m128i) _Float16 out[2];
//   out[0]=I2(cvtsi128_si64,m)?p[I2(cvtsi128_si64,i)]:(_Float16)0;out[1]=I2(extract_epi64,m,1)?p[I2(extract_epi64,i,1)]:(_Float16)0;
//   Rh2(out)),c(_Float16)*p,c(simde__m128i)i,c(simde__m128i)m)
// (mask) gather complex
#if M == 1
// D(Vd,mgthrx,_(Vd r;switch(t){
//     case 1:{c(f32)*_p=p;r=I8(cvtps_pd,mgthf8(_p,i,m));break;}
//     case 2:{c(f16)*_p=p;_Alignas(Vi)i64 _i[8];_Alignas(Vd)f64 out[8];I8(storeu_epi64,_i,i);_L(u,8,out[u]=((m>>u)&1)?(f64)_p[_i[u]]:0.0f);r=Rd8(out);break;}
//     default:{c(f64)*_p=p;r=mgthd8(_p,i,m);break;}
// };r),c(i32)t,c(void)*p,c(Vi)i,c(M8)m)
// D(Vz,mgthzx,_(Vz r;switch(t){
//     case 1:{simde__m512 a=I8(castpd_ps,mgthd8((c(f64)*)p,i,m));r=cvtf2N_Vz(a);break;}
//     case 2:{c(f16)*_p=p;_Alignas(Vi)i64 _i[8];_Alignas(Vd)f64 re[8],im[8];I8(storeu_epi64,_i,i);_L(u,8,re[u]=((m>>u)&1)?(f64)_p[2*_i[u]]:0.0;im[u]=((m>>u)&1)?(f64)_p[2*_i[u]+1]:0.0);r=Z2(Rd8(re),Rd8(im));break;}
//     default:{c(f64)*_p=p;r=Z2(mgthd8(_p,shl(i,1),m),mgthd8(_p+1,shl(i,1),m));break;}
// };r),c(i32)t,c(void)*p,c(Vi)i,c(M8)m)
#elif M == 2
// D(Vd,mgthrx,_(Vd r;switch(t){
//     case 1:{c(f32)*_p=p;r=I4(cvtps_pd,mgthf4(_p,i,maski64_i32(m)));break;}
//     case 2:{c(f16)*_p=p;_Alignas(Vi)i64 _i[4],_m[4];_Alignas(Vd)f64 out[4];I4(storeu_epi64,_i,i);I4(storeu_epi64,_m,m);_L(u,4,out[u]=_m[u]?(f64)_p[_i[u]]:0.0f);r=Rd4(out);break;}
//     default:{c(f64)*_p=p;r=mgthd4(_p,i,m);break;}
// };r),c(i32)t,c(void)*p,c(Vi)i,c(M8)m)
// D(Vd,gthrx,_(Vd r;switch(t){
//     case 1:{c(f32)*_p=p;r=I4(cvtps_pd,gthf4(_p,i));break;}
//     case 2:{c(f16)*_p=p;_Alignas(Vi)i64 _i[4],_m[4];_Alignas(Vd)f64 out[4];I4(storeu_epi64,_i,i);I4(storeu_epi64,_m,m);_L(u,4,out[u]=_m[u]?(f64)_p[_i[u]]:0.0f);r=Rd4(out);break;}
//     default:{c(f64)*_p=p;r=mgthd4(_p,i,m);break;}
// };r),c(i32)t,c(void)*p,c(Vi)i,c(M8)m)

// D(Vz,mgthzx,_(Vz r;switch(t){
//     case 1:{simde__m256 a=I4(castpd_ps,mgthd4((c(f64)*)p,i,m));r=cvtf2N_Vz(a);break;}
//     case 2:{c(f16)*_p=p;_Alignas(Vi)i64 _i[4],_m[4];_Alignas(Vd)f64 re[4],im[4];I4(storeu_epi64,_i,i);I4(storeu_epi64,_m,m);_L(u,4,re[u]=_m[u]?(f64)_p[2*_i[u]]:0.0;im[u]=_m[u]?(f64)_p[2*_i[u]+1]:0.0);r=Z2(Rd4(re),Rd4(im));break;}
//     default:{c(f64)*_p=p;r=Z2(mgthd4(_p,shl(i,1),m),mgthd4(_p+1,shl(i,1),m));break;}
// };r),c(i32)t,c(void)*p,c(Vi)i,c(M8)m)
#else
// D(Vd,mgthrx,_(Vd r;switch(t){
//     case 1:{c(f32)*_p=p;r=I2(cvtps_pd,mgthf2(_p,i,m));break;}
//     case 2:{c(f16)*_p=p;r=I2(cvtps_pd,I2(cvtph_ps,mgthh2(_p,i,m)));break;}
//     default:{c(f64)*_p=p;r=mgthd2(_p,i,m);break;}
// };r),c(i32)t,c(void)*p,c(Vi)i,c(M8)m)
// D(Vz,mgthzx,_(Vz r;switch(t){
//     case 1:{simde__m128 a=I2(castpd_ps,mgthd2((c(f64)*)p,i,m));r=cvtf4_Vz(a);break;}
//     case 2:{c(f16)*_p=p;_Alignas(Vi)i64 _i[2],_m[2];_Alignas(Vd)f64 re[2],im[2];I2(store_si128,(void*)_i,i);I2(store_si128,(void*)_m,m);_L(u,2,re[u]=_m[u]?(f64)_p[2*_i[u]]:0.0;im[u]=_m[u]?(f64)_p[2*_i[u]+1]:0.0);r=Z2(Rd2(re),Rd2(im));break;}
//     default:{c(f64)*_p=p;r=Z2(mgthd2(_p,shl(i,1),m),mgthd2(_p+1,shl(i,1),m));break;}
// };r),c(i32)t,c(void)*p,c(Vi)i,c(M8)m)
#endif

// read real
// #if M == 1
// D(Vd,Rdx,_(Vd r;switch(t){
//     case 1:{r=I8(cvtps_pd,Rf8((c(f32)*)p));break;}
//     case 2:{simde__m128i a=Rh8(p);r=I8(cvtps_pd,I4(cvtph_ps,a));break;}
//     default:{r=Rd8(p);break;}
// };r),c(i32)t,c(void)*p)
// #elif M == 2
// D(Vd,Rdx,_(Vd r;switch(t){
//     case 1:{r=I4(cvtps_pd,Rf4((c(f32)*)p));break;}
//     case 2:{simde__m128i a=Rh4(p);r=I4(cvtps_pd,I2(cvtph_ps,a));break;}
//     default:{r=Rd4(p);break;}
// };r),c(i32)t,c(void)*p)
// #else
// D(Vd,Rdx,_(Vd r;switch(t){
//     case 1:{r=I2(cvtps_pd,Rf2((c(f32)*)p));break;}
//     case 2:{c(f16)*_p=p;_Alignas(Vd)f64 out[2]={(f64)_p[0],(f64)_p[1]};r=Rd2(out);break;}
//     default:{r=Rd2(p);break;}
// };r),c(i32)t,c(void)*p)
// #endif

// write real
// #if M == 1
// D(void,Wdx,_(switch(t){
//     case 1:{Wf8(p,I8(cvtpd_ps,z));break;}
//     case 2:{Wh8(p,I4(cvtps_ph,I8(cvtpd_ps,z),SIMDE_MM_FROUND_TO_NEAREST_INT));break;}
//     default:{Wd8(p,z);break;}
// }),c(i32)t,void*p,c(Vd)z)
// #elif M == 2
// D(void,Wdx,_(switch(t){
//     case 1:{Wf4(p,I4(cvtpd_ps,z));break;}
//     case 2:{Wh4(p,I2(cvtps_ph,I4(cvtpd_ps,z),SIMDE_MM_FROUND_TO_NEAREST_INT));break;}
//     default:{Wd4(p,z);break;}
// }),c(i32)t,void*p,c(Vd)z)
// #else
// D(void,Wdx,_(switch(t){
//     case 1:{Wf2(p,I2(cvtpd_ps,z));break;}
//     case 2:{Wh2(p,I2(cvtps_ph,I2(cvtpd_ps,z),SIMDE_MM_FROUND_TO_NEAREST_INT));break;}
//     default:{Wd2(p,z);break;}
// }),c(i32)t,void*p,c(Vd)z)
// #endif

// read complex
// #if M == 1
// D(Vz,Rzx,_(Vz r;switch(t){
//     case 1:{simde__m512 a=Rf16((c(f32)*)p);r=cvtf2N_Vz(a);break;}
//     case 2:{simde__m256i a=I4(loadu_si256,p);_Alignas(Vi)c(i16)idx[16]={0,2,4,6,8,10,12,14,1,3,5,7,9,11,13,15};simde__m256i b=I4(permutexvar_epi16,I4(loadu_epi16,idx),a);r=Z2(I8(cvtps_pd,I4(cvtph_ps,I4(castsi256_si128,b))),I8(cvtps_pd,I4(cvtph_ps,I4(extracti128_si256,b,1))));break;}
//     default:{c(f64)*_p=p;c(Vd)a=Rd8(_p),b=Rd8(_p+8);_Alignas(Vi)c(i64)idx[8]={0,2,4,6,1,3,5,7};c(Vd)c=I8(permutexvar_pd,I8(load_epi64,idx),I8(unpacklo_pd,a,b)),d=I8(permutexvar_pd,I8(load_epi64,idx),I8(unpackhi_pd,a,b));r=Z2(c,d);break;}
// };r),c(i32)t,c(void)*p)
// #elif M == 2
// D(Vz,Rzx,_(Vz r;switch(t){
//     case 1:{simde__m256 a=Rf8((c(f32)*)p);r=cvtf2N_Vz(a);break;}
//     case 2:{simde__m128i a=Rh8((c(f16)*)p),b=I2(shuffle_epi8,a,I2(set_epi8,15,14,11,10,7,6,3,2,13,12,9,8,5,4,1,0));simde__m256 c=I4(cvtph_ps,b);r=Z2(I4(cvtps_pd,I4(castps256_ps128,c)),I4(cvtps_pd,I4(extractf128_ps,c,1)));break;}
//     default:{c(f64)*_p=p;c(Vd)a=Rd4(_p),b=Rd4(_p+4),c=I4(permute2f128_pd,a,b,0b00110001),d=I4(permute2f128_pd,a,b,0b00100000);r=Z2(I(unpacklo_pd,d,c),I(unpackhi_pd,d,c));break;}
// };r),c(i32)t,c(void)*p)
// #else
// D(Vz,Rzx,_(Vz r;switch(t){
//     case 1:{simde__m128 a=Rf4((c(f32)*)p);r=Z2(I2(cvtps_pd,I2(shuffle_ps,a,I2(setzero_ps),0b00001000)),I2(cvtps_pd,I2(shuffle_ps,a,I2(setzero_ps),0b00001101)));break;}
//     case 2:{c(f16)*_p=p;r=Z2(I2(set_pd,(f64)_p[2],(f64)_p[0]),I2(set_pd,(f64)_p[3],(f64)_p[1]));break;}
//     default:{c(f64)*_p=p;c(Vd)a=Rd2(_p),b=Rd2(_p+2);r=Z2(I2(unpacklo_pd,a,b),I2(unpackhi_pd,a,b));break;}
// };r),c(i32)t,c(void)*p)
// #endif

// write complex
// #if M == 1
// D(void,Wzx,_(switch(t){
//     case 1:{simde__m256 a=I8(cvtpd_ps,z.re),b=I8(cvtpd_ps,z.im),c=I4(unpacklo_ps,a,b),d=I4(unpackhi_ps,a,b);f32*_p=p;Wf8(_p,I4(permute2f128_ps,c,d,0b00100000));Wf8(_p+8,I4(permute2f128_ps,c,d,0b00110001)); break;}
//     case 2:{simde__m128i a=I8(cvtpd_ph,z.re),b=I8(cvtpd_ph,z.im),c=I2(unpacklo_epi16,a,b),d=I2(unpackhi_epi16,a,b);f16*_p=p;Wh8(_p,c);Wh8(_p+8,d);break;}
//     default:{
//         _Alignas(Vi)c(i64)idx1[N]={0b0000,0b0001,0b1000,0b1001,0b0010,0b0011,0b1010,0b1011};
//         _Alignas(Vi)c(i64)idx2[N]={0b0100,0b0101,0b1100,0b1101,0b0110,0b0111,0b1110,0b1111};
//         c(Vd)a=I(unpacklo_pd,z.re,z.im),b=I(unpackhi_pd,z.re,z.im),c=I(permutex2var_pd,a,I(load_epi64,idx1),b),d=I(permutex2var_pd,a,I(load_epi64,idx2),b);
//         f64*_p=p;Wd8(p,c);Wd8(_p+8,d);break;}
// }),c(i32)t,void*p,c(Vz)z)
// #elif M == 2
// D(void,Wzx,_(switch(t){
//     case 1:{simde__m128 a=I4(cvtpd_ps,z.re),b=I4(cvtpd_ps,z.im);f32*_p=p;Wf4(_p,I2(unpacklo_ps,a,b));Wf4(_p+4,I2(unpackhi_ps,a,b));break;} // f32
//     case 2:{_Alignas(Vz)f64 tmp[8];Wd4(tmp,z.re);Wd4(tmp+4,z.im);f16*_p=p;
//         _p[0]=(_Float16)tmp[0];_p[1]=(_Float16)tmp[4];_p[2]=(_Float16)tmp[1];_p[3]=(_Float16)tmp[5];
//         _p[4]=(_Float16)tmp[2];_p[5]=(_Float16)tmp[6];_p[6]=(_Float16)tmp[3];_p[7]=(_Float16)tmp[7];break;} // f16
//     default:{c(Vd)a=I4(unpacklo_pd,z.re,z.im),b=I4(unpackhi_pd,z.re,z.im),c=I(permute2f128_pd,a,b,0b00110001),d=I(permute2f128_pd,a,b,0b00100000);f64*_p=p;Wd4(_p,d);Wd4(_p+4,c);break;} // f64
// }),c(i32)t,void*p,c(Vz)z)
// #else
// D(void,Wzx,_(switch(t){
//     case 1:{simde__m128 a=I2(cvtpd_ps,z.re),b=I2(cvtpd_ps,z.im);f32*_p=p;Wf4(_p,I2(unpacklo_ps,a,b));break;} // f32
//     case 2:{_Alignas(Vd)f64 tmp[4];Wd2(tmp,z.re);Wd2(tmp+2,z.im);f16*_p=p;_p[0]=(_Float16)tmp[0];_p[1]=(_Float16)tmp[2];_p[2]=(_Float16)tmp[1];_p[3]=(_Float16)tmp[3];break;} // f16
//     default:{f64*_p=p;Wd2(_p,I2(unpacklo_pd,z.re,z.im));Wd2(_p+2,I2(unpackhi_pd,z.re,z.im));break;} // f64
// }),c(i32)t,void*p,c(Vz)z)
// #endif




#if M == 1
#  define _gthd(p,i,s) I(i64gather_pd,i,p,s)
#  define _gthf(p,i,s) I(i64gather_ps,i,p,s)
#  define _mgthd(p,i,m,s) I(mask_i64gather_pd,Zd,m,i,p,s)
#  define _mgthf(p,i,m,s) I(mask_i64gather_ps,I4(setzero_ps),m,i,p,s)
#  define _M8_to_i(m) m
#else
#  define _gthd(p,i,s) I(i64gather_pd,p,i,s)
#  define _gthf(p,i,s) I(i64gather_ps,p,i,s)
#  define _mgthd(p,i,m,s) I(mask_i64gather_pd,Zd,p,i,i2d(m),s)
#  define _mgthf(p,i,m,s) I(mask_i64gather_ps,I2(setzero_ps),p,i,I2(castsi128_ps,maski64_i32(m)),s)
#  define _M8_to_i(m) I(movemask_pd,i2d(m))
#endif

#if M == 1
#  define cvtf2N_Vz(a) _(_Alignas(Vi)c(i32)idx[16]={0,2,4,6,8,10,12,14,1,3,5,7,9,11,13,15};simde__m512 b=I8(permutexvar_ps,I8(load_epi32,idx),a);r=Z2(I8(cvtps_pd,I8(castps512_ps256,b)),I8(cvtps_pd,I8(extractf32x8_ps,b,1))))
#elif M == 2
#  define maski64_i32(m) I4(castsi256_si128,I4(permute4x64_epi64,I4(shuffle_epi32,m,0x88),0xd8))
#  define cvtf2N_Vz(a) _(simde__m256 b=I4(shuffle_ps,a,a,0b11011000),c=I4(castpd_ps,I4(permute4x64_pd,I4(castps_pd,b),0b11011000));Z2(I4(cvtps_pd,I4(castps256_ps128,c)),I4(cvtps_pd,I4(extractf128_ps,c,1))))
#else
#  define cvtf2N_Vz(a) _(Z2(I2(cvtps_pd,I2(shuffle_ps,a,I2(setzero_ps),0x08)),I2(cvtps_pd,I2(shuffle_ps,a,I2(setzero_ps),0x0d))))
#endif

#if USE_F16
#  if M == 1
#    define Rh(p) Rh8(p)
#    define Wh(p,z) Wh8(p,z)
#  elif M == 2
#    define Rh(p) Rh4(p)
#    define Wh(p,z) Wh4(p,z)
#  else
#    define Rh(p) Rh2(p)
#    define Wh(p,z) Wh2(p,z)
#  endif
#endif

De(Vz,Rx,_(Vz r;switch(t){
     /*f64*/case 0: r=Z2(Rd(p),Zd);break;
     /*f32*/case 1: r=Z2(I(cvtps_pd,Ih(loadu_ps,p)),Zd);break;
#if M == 1
    /*c128*/case 3: {c(f64)*_p=p;c(Vd)a=Rd(_p),b=Rd(_p+N);_Alignas(Vi)c(i64)idx[8]={0,2,4,6,1,3,5,7};c(Vd)c=I(permutexvar_pd,I(load_epi64,idx),I(unpacklo_pd,a,b)),d=I(permutexvar_pd,I(load_epi64,idx),I(unpackhi_pd,a,b));r=Z2(c,d);break;}
#elif M == 2
    /*c128*/case 3: {c(f64)*_p=p;c(Vd)a=Rd(_p),b=Rd(_p+N),c=I4(permute2f128_pd,a,b,0x31),d=I4(permute2f128_pd,a,b,0x20);r=Z2(I(unpacklo_pd,d,c),I(unpackhi_pd,d,c));break;}
#else
    /*c128*/case 3: {c(f64)*_p=p;c(Vd)a=Rd(_p),b=Rd(_p+N);r=Z2(I2(unpacklo_pd,a,b),I2(unpackhi_pd,a,b));break;}
#endif
     /*c64*/case 4: r=cvtf2N_Vz(I(loadu_ps,(c(f32)*)p));break;
#if USE_F16
     /*f16*/case 2: r=Z2(I(cvtps_pd,Ih(cvtph_ps,Rh(p))),Zd);break;
    /*c32*/default: {c(f16)*_p=p;_Alignas(Vd)f64 re[N],im[N];_L(u,N,re[u]=_p[2*u],im[u]=_p[2*u+1]);r=Z2(Rd(re),Rd(im));break;}
#else
     /*f16*/case 2: r=Z2(Zd,Zd);break;
    /*c32*/default: r=Z2(Zd,Zd);break;
#endif
};r),c(i32)t,c(void)*p)

De(void,Wx,_(switch(t){
     /*f64*/case 0: I(storeu_pd,p,z.re);break;
     /*f32*/case 1: Ih(storeu_ps,p,I(cvtpd_ps,z.re));break;
#if M == 1
    /*c128*/case 3: {_Alignas(Vi)c(i64)idx1[N]={0b0000,0b0001,0b1000,0b1001,0b0010,0b0011,0b1010,0b1011},idx2[N]={0b0100,0b0101,0b1100,0b1101,0b0110,0b0111,0b1110,0b1111};
        c(Vd)a=I(unpacklo_pd,z.re,z.im),b=I(unpackhi_pd,z.re,z.im),c=I(permutex2var_pd,a,I(load_epi64,idx1),b),d=I(permutex2var_pd,a,I(load_epi64,idx2),b);
        f64*_p=p;Wd(_p,c);Wd(_p+N,d);break;}
#elif M == 2
    /*c128*/case 3: {c(Vd)a=I(unpacklo_pd,z.re,z.im),b=I(unpackhi_pd,z.re,z.im);f64*_p=p;Wd(_p,I(permute2f128_pd,a,b,0x20));Wd(_p+N,I(permute2f128_pd,a,b,0x31));break;}
#else
    /*c128*/case 3: {c(Vd)a=I(unpacklo_pd,z.re,z.im),b=I(unpackhi_pd,z.re,z.im);f64*_p=p;Wd(_p,a);Wd(_p+N,b);break;}
#endif
     /*c64*/case 4: {f32*_p=p;typeof(I(cvtpd_ps,z.re))a=I(cvtpd_ps,z.re),b=I(cvtpd_ps,z.im),c=Ih(unpacklo_ps,a,b),d=Ih(unpackhi_ps,a,b);
#if M == 1
        I4(storeu_ps,_p,I4(permute2f128_ps,c,d,0x20));I4(storeu_ps,_p+N,I4(permute2f128_ps,c,d,0x31));
#elif M == 2
        I2(storeu_ps,_p,c);I2(storeu_ps,_p+N,d);
#else
        I2(storeu_ps,_p,c);(void)d;
#endif
        break;}
#if USE_F16
     /*f16*/case 2: Wh(p,Ih(cvtps_ph,I(cvtpd_ps,z.re),SIMDE_MM_FROUND_TO_NEAREST_INT));break;
    /*c32*/default: {f16*_p=p;_Alignas(Vd)f64 re[N],im[N];Wd(re,z.re);Wd(im,z.im);_L(u,N,_p[2*u]=(f16)re[u],_p[2*u+1]=(f16)im[u]);break;}
#else
     /*f16*/case 2: break;
    /*c32*/default: break;
#endif
}),c(i32)t,void*p,c(Vz)z)

De(Vz,Gmx,_(Vz r;switch(t){
     /*f64*/case 0: r=Z2(_mgthd((c(f64)*)p,i,m,8),Zd);break;
     /*f32*/case 1: r=Z2(I(cvtps_pd,_mgthf((c(f32)*)p,i,m,4)),Zd);break;
    /*c128*/case 3: r=Z2(_mgthd((c(f64)*)p,shl(i,1),m,8),_mgthd((c(f64)*)p+1,shl(i,1),m,8));break;
     /*c64*/case 4: r=cvtf2N_Vz(d2f(_mgthd((c(f64)*)p,i,m,8)));break;
#if USE_F16
     /*f16*/case 2: {c(f16)*_p=p;_Alignas(Vi)i64 _i[N];_Alignas(Vd)f64 re[N];c(u32)_m=_M8_to_i(m);Wi(_i,i);_L(u,N,re[u]=((_m>>u)&1)?_p[_i[u]]:0.0);r=Z2(Rd(re),Zd);break;}
    /*c32*/default: {c(f16)*_p=p;_Alignas(Vi)i64 _i[N];_Alignas(Vd)f64 re[N],im[N];c(u32)_m=_M8_to_i(m);Wi(_i,i);_L(u,N,re[u]=((_m>>u)&1)?_p[2*_i[u]]:0.0,im[u]=((_m>>u)&1)?_p[2*_i[u]+1]:0.0);r=Z2(Rd(re),Rd(im));break;}
#else
     /*f16*/case 2: r=Z2(Zd,Zd);break;
    /*c32*/default: r=Z2(Zd,Zd);break;
#endif
};r),c(i32)t,c(void)*p,c(Vi)i,c(M8)m)

De(Vz,Gx,_(Vz r; switch(t){
     /*f64*/case 0: r=Z2(_gthd((c(f64)*)p,i,8),Zd);break;
     /*f32*/case 1: r=Z2(I(cvtps_pd,_gthf((c(f32)*)p,i,4)),Zd);break;
    /*c128*/case 3: r=Z2(_gthd((c(f64)*)p,shl(i,1),8),_gthd((c(f64)*)p+1,shl(i,1),8));break;
     /*c64*/case 4: r=cvtf2N_Vz(d2f(_gthd((c(f64)*)p,i,8)));break;
#if USE_F16
     /*f16*/case 2: {c(f16)*_p=p;_Alignas(Vi)i64 _i[N];_Alignas(Vd)f64 re[N];Wi(_i,i);_L(u,N,re[u]=_p[_i[u]]);r=Z2(Rd(re),Zd);break;}
    /*c32*/default: {c(f16)*_p=p;_Alignas(Vi)i64 _i[N];_Alignas(Vd)f64 re[N],im[N];Wi(_i,i);_L(u,N,re[u]=_p[2*_i[u]],im[u]=_p[2*_i[u]+1]);r=Z2(Rd(re),Rd(im));break;}
#else
     /*f16*/case 2: r=Z2(Zd,Zd); break;
    /*c32*/default: r=Z2(Zd,Zd); break;
#endif
};r),c(i32)t,c(void)*p,c(Vi)i)


#if M == 1
    D(Vi,gthq,I(i64gather_epi64,i,p,8),c(u64)*p,c(Vi)i)
    D(Vd,gthd,I(i64gather_pd,i,p,8),c(f64)*p,c(Vi)i)D(Vd,mgthd,I8(mask_i64gather_pd,Zd,m,i,p,8),c(f64)*p,c(Vi)i,c(M8)m)
    // D(Vd,mgthps2pd,I(cvtps_pd,I(mask_i64gather_ps,I4(setzero_ps),m,i,p,4)),c(f32)*p,c(Vi)i,c(M8)m)
    D(Vd,mgthw2d,_(I(cvtepi32_pd,I4(blend_epi16,I(mask_i64gather_epi32,I4(setzero_si256),m,i,(c(i32)*)p,2),I4(setzero_si256),0b10101010))),c(u16)*p,c(Vi)i,c(M8)m)
#elif M == 2
    D(Vi,gthq,I4(i64gather_epi64,p,i,8),c(u64)*p,c(Vi)i)
    D(Vd,gthd,I4(i64gather_pd,p,i,8),c(f64)*p,c(Vi)i)D(Vd,mgthd,I4(mask_i64gather_pd,Zd,p,i,i2d(m),8),c(f64)*p,c(Vi)i,c(Vi)m)
    D(Vd,mgthw2d,_(
        simde__m128i m0=I4(castsi256_si128,m),m1=I4(extractf128_si256,m,1);
        simde__m128i _m=I2(castps_si128,I2(shuffle_ps,I2(castsi128_ps,m0),I2(castsi128_ps,m1),0b10001000));
        I4(cvtepi32_pd,I2(blend_epi16,I4(mask_i64gather_epi32,I2(setzero_si128),(c(i32)*)p,i,_m,2),I2(setzero_si128),0b10101010))
    ),c(u16)*p,c(Vi)i,c(Vi)m)
#else
    D(Vi,gthq,I2(set_epi64x,p[I2(extract_epi64,i,1)],p[I2(cvtsi128_si64,i)]),c(u64)*p,c(Vi)i)
    D(Vd,gthd,I2(set_pd,p[I2(extract_epi64,i,1)],p[I2(cvtsi128_si64,i)]),c(f64)*p,c(Vi)i)D(Vd,mgthd,I2(set_pd,I2(extract_epi64,m,1)?p[I2(extract_epi64,i,1)]:0.0,I2(cvtsi128_si64,m)?p[I2(cvtsi128_si64,i)]:0.0),c(f64)*p,c(Vi)i,c(Vi)m)
    D(Vd,mgthw2d,I2(set_pd,I2(extract_epi64,m,1)?p[I2(extract_epi64,i,1)]:0.0,I2(cvtsi128_si64,m)?p[I2(cvtsi128_si64,i)]:0.0),c(u16)*p,c(Vi)i,c(Vi)m)
#endif
// D(Vz,gthz,Z2(gthd((c(f64)*)p,shl(i,1)),gthd((c(f64)*)p+1,shl(i,1))),c(c128)*p,c(Vi)i)D(Vz,mgthz,Z2(mgthd((c(f64)*)p,shl(i,1),m),mgthd((c(f64)*)p+1,shl(i,1),m)),c(c128)*p,c(Vi)i,c(M8)m)
// D(Vz,mgthzs2zd,Z2(mgthps2pd((c(f32)*)p,shl(i,1),m),mgthps2pd((c(f32)*)p+1,shl(i,1),m)),c(c64)*p,c(Vi)i,c(M8)m)

#if M == 1
    D(Vi,m1,I(movm_epi64,I(test_epi64_mask,x,m)),c(Vi)x,c(Vi)m)
#else
    D(Vi,m1,xor(eqi(and(x, m), Zi), eqi(Zi, Zi)),c(Vi)x,c(Vi)m)
#endif
D(Vi,m2,xor(m1(x,m_0),m1(x,m_1)),c(Vi)x,c(Vi)m_0,c(Vi)m_1)D(Vi,mX,popcnt(and(x,m)),c(Vi)x,c(Vi)m)
D(Vd,bcast2d,Sd(re[k]),c(f64)*re,c(f64)*im,c(i32)k)D(Vz,bcast2z,(Z2(Sd(re[k]),Sd(im[k]))),c(f64)*re,c(f64)*im,c(i32)k)
D(Vd,signedd,i2d(xor(d2i(v),shl(m,63))),c(Vd)v,c(Vi)m)D(Vz,signedz,_(m=shl(m,63);Z2(i2d(xor(d2i(v.re),m)),i2d(xor(d2i(v.im),m)))),c(Vz)v,Vi m)

// TODO: The loop is ugly, but the compilers vectorize it
Di(void,prefetchq,_(_Alignas(Vi)i64 buf[N];Wi(buf,idx);_L(k,N,I2(prefetch,p+buf[k],_MM_HINT_T0))),c(u64)*p,c(Vi)idx)

// TODO: Should get rid of these ...
// static void prefetchd4xN(Vi idx[4],c(u16)*n,c(f64)*x){_L(k,4*N,c(i32)i=((c(i64)*)idx)[k];simde_mm_prefetch(n+i,_MM_HINT_ET0);simde_mm_prefetch(x+i,_MM_HINT_ET0))}
// static void prefetchz4xN(Vi idx[4],c(u16)*n,c(c128)*x){_L(k,4*N,c(i32)i=((c(i64)*)idx)[k];simde_mm_prefetch(n+i,_MM_HINT_ET0);simde_mm_prefetch(x+i,_MM_HINT_ET0))}

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

