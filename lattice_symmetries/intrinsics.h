#pragma once

#include <simde/x86/avx512.h>

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

#define INTERNAL static HEDLEY_ALWAYS_INLINE

typedef int8_t i8; typedef int16_t i16; typedef int32_t i32; typedef int64_t i64;
typedef uint8_t u8; typedef uint16_t u16; typedef uint32_t u32; typedef uint64_t u64;
typedef float f32; typedef double f64; typedef float _Complex c64; typedef double _Complex c128;

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
#   define i2d(x) I(castsi512_pd,x)
#   define i2s(x) I(castsi512_ps,x)
    typedef simde__m512i Vi; typedef simde__m512d Vd; typedef simde__m512 Vs; typedef simde__mmask8 M8;
#elif M == 2
#   define B 256
#   define _S(f) f##_si256
#   define I(f,x...) I4(f,x)
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
#define shl(x,imm) I(slli_epi64,x,imm)
#define shr(x,imm) I(srli_epi64,x,imm)

// Vi
D(Vi,Ri,I(loadu_epi64,p),c(void)*p)D(void,Wi,I(_S(storeu),p,x),void*p,c(Vi)x)D2(Vi,and,_S(and))D2(Vi,xor,_S(xor))D2(Vi,andnot,_S(andnot))D2(Vi,or,_S(or))D2(Vi,addq,add_epi64)D2(Vi,subq,sub_epi64)D(Vi,popcnt,I(popcnt_epi64,x),c(Vi)x)
// Vd
Dd(zero,Zd,void)D(Vd,Rd,I(loadu_pd,p),c(f64)*p)D(Vd,Sd,I(set1_pd,x),c(f64)x)D(void,Wd,I(storeu_pd,p,x),f64*p,c(Vd)x)D2(Vd,addd,add_pd)D2(Vd,subd,sub_pd)D2(Vd,muld,mul_pd)D2(Vd,divd,div_pd)D(Vd,fmad,I(fmadd_pd,a,b,c),c(Vd)a,c(Vd)b,c(Vd)c)D(Vd,fnmad,I(fnmadd_pd,a,b,c),c(Vd)a,c(Vd)b,c(Vd)c)D(Vd,sqrtd,I(sqrt_pd,x),c(Vd)x)
// Vz
typedef struct Vz { Vd re; Vd im; } Vz;
Dz(zero,(Z2(Zd,Zd)),void)D(Vz,addz,Z2(addd(a.re, b.re),addd(a.im, b.im)),c(Vz)a,c(Vz)b)D(Vz,mulz,Z2(fnmad(a.im,b.im,muld(a.re,b.re)),fmad(a.im,b.re,muld(a.re,b.im))),c(Vz)a,c(Vz)b)
#if M == 1
    D(Vz,Rz,_(
        c(Vd)a=Rd((c(f64)*)p),b=Rd((c(f64)*)p+N);_Alignas(Vi)c(i64)idx[N]={0,2,4,6,1,3,5,7};
        c(Vd)c=I(permutexvar_pd,I(load_epi64,idx),I(unpacklo_pd,a,b)),d=I(permutexvar_pd,I(load_epi64,idx),I(unpackhi_pd,a,b));
        Z2(c,d)), c(c128)*p)
    // TODO: boy is this one ugly... is there a better way?
    De(void,Wz,_(
        _Alignas(Vi)c(i64)idx1[N]={0b0000,0b0001,0b1000,0b1001,0b0010,0b0011,0b1010,0b1011};
        _Alignas(Vi)c(i64)idx2[N]={0b0100,0b0101,0b1100,0b1101,0b0110,0b0111,0b1110,0b1111};
        c(Vd)a=I(unpacklo_pd,z.re,z.im),b=I(unpackhi_pd,z.re,z.im);
        c(Vd)c=I(permutex2var_pd,a,I(load_epi64,idx1),b),d=I(permutex2var_pd,a,I(load_epi64,idx2),b);
        Wd((f64*)p,c);Wd((f64*)p+N,d)),c(c128)*p,c(Vz)z)
#elif M == 2
    D(Vz,Rz,_(c(Vd)a=Rd((c(f64)*)p),b=Rd((c(f64)*)p+N),c=I(permute2f128_pd,a,b,0b00110001),d=I(permute2f128_pd,a,b,0b00100000);Z2(I(unpacklo_pd,d,c),I(unpackhi_pd,d,c))),c(c128)*p)
    D(void,Wz,_(c(Vd)a=I4(unpacklo_pd,z.re,z.im),b=I4(unpackhi_pd,z.re,z.im),c=I(permute2f128_pd,a,b,0b00110001),d=I(permute2f128_pd,a,b,0b00100000);Wd((f64*)p,d),Wd((f64*)p+N,c)),c128*p,c(Vz)z)
#else
    D(Vz,Rz,_(c(Vd)a=Rd((c(f64)*)p),b=Rd((c(f64)*)p+N);Z2(I(unpacklo_pd,a,b),I(unpackhi_pd,a,b))),c(c128)*p)
    D(void,Wz,_(Wd((f64*)p,I(unpacklo_pd,z.re,z.im)),Wd((f64*)p+N,I(unpackhi_pd,z.re,z.im))),c128*p,c(Vz)z)
#endif
// u16
#if M == 1
    // Why is SIMDe missing these functions...?
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

#if M == 1
    D(Vi,gthq,I(i64gather_epi64,i,p,8),c(u64)*p,c(Vi)i)
    D(Vd,gthd,I(i64gather_pd,i,p,8),c(f64)*p,c(Vi)i)D(Vd,mgthd,I8(mask_i64gather_pd,Zd,m,i,p,8),c(f64)*p,c(Vi)i,c(M8)m)
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
D(Vz,gthz,Z2(gthd((c(f64)*)p,shl(i,1)),gthd((c(f64)*)p+1,shl(i,1))),c(c128)*p,c(Vi)i)D(Vz,mgthz,Z2(mgthd((c(f64)*)p,shl(i,1),m),mgthd((c(f64)*)p+1,shl(i,1),m)),c(c128)*p,c(Vi)i,c(M8)m)

#if M == 1
    D(Vi,m1,I(movm_epi64,I(test_epi64_mask,x,m)),c(Vi)x,c(Vi)m)
#else
    D(Vi,m1,xor(eqi(and(x, m), Zi), eqi(Zi, Zi)),c(Vi)x,c(Vi)m)
#endif
D(Vi,m2,xor(m1(x,m_0),m1(x,m_1)),c(Vi)x,c(Vi)m_0,c(Vi)m_1)D(Vi,mX,popcnt(and(x,m)),c(Vi)x,c(Vi)m)
Dd(bcast2,Sd(re[k]),c(f64)*re,c(f64)*im,c(i32)k)
Dz(bcast2,(Z2(Sd(re[k]),Sd(im[k]))),c(f64)*re,c(f64)*im,c(i32)k)
D(Vd,signedd,i2d(xor(d2i(v),shl(m,63))),c(Vd)v,c(Vi)m)
D(Vz,signedz,_(m=shl(m,63);Z2(i2d(xor(d2i(v.re),m)),i2d(xor(d2i(v.im),m)))),c(Vz)v,Vi m)

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

