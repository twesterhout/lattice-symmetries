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

#if M == 1
#  include <immintrin.h>
#  include <simde/x86/avx512.h>
#elif M == 3
   // We rely on SIMDe to implement WASM intrinsics using SSE or NEON
   #include <simde/wasm/simd128.h>
#endif

typedef float float32_t; // Need it for some weird reason at M == 3 ...

#if M == 1 // 512
#  define A(f1,f2,f3) f1
#elif M == 2 // 256
#  define A(f1,f2,f3) f2
#else // 128
#  define A(f1,f2,f3) f3
#endif

#if defined(__clang__)
#  define GC(g,c) clang
#else
#  define GC(g,c) g
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
#define PX(x, t, w, f) _(printf("[");for(u64 k=0;k<B/(8*sizeof(t));++k){printf(f ",",(x)[k]);}printf("]\n"))
#define Pw(x) PX(x,i32,"%i")
#define Pi(x) PX(x,i64,"%zi")
#define Ps(x) PX(x,f32,"%e")
#define Pd(x) PX(x,f64,"%e")

typedef char i8; typedef short i16; typedef int i32; typedef long long i64;
typedef unsigned char u8; typedef unsigned short u16; typedef unsigned u32; typedef unsigned long long u64;
typedef float f32; typedef double f64; typedef float _Complex c64; typedef double _Complex c128;

#define E_(x...) x
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
// #define I2(f,x...) I_(mm,f,x)
// #define I4(f,x...) I_(mm256,f,x)
#define I8(f,x...) I_(mm512,f,x)
#define O_(f) A(__builtin_ia32_##f,__builtin_ia32_##f,simde_wasm_##f)
#define O(f) O_(f)

#define R2(x) x,x
#define R4(x) x,x,x,x
#define R8(x) x,x,x,x,x,x,x,x
#define R16(x) x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x
#define RF2(f) f(0),f(1)
#define RF4(f) f(0),f(1),f(2),f(3)
#define RF8(f) f(0),f(1),f(2),f(3),f(4),f(5),f(6),f(7)
#define B A(512,256,128)
#define N (B/64)
#define RN(x) A(R8,R4,R2)(x)
#define RFN(f) A(RF8,RF4,RF2)(f)
#define Zi ((Vq){RN(0)})
#define Zd ((Vd){RN(0.0)})
#define Zf ((Vd){RN(0.0f),RN(0.0f)})
#define Z4f ((V4f){R4(0.0f)})
#define Va(n) __attribute__((vector_size(n),aligned(n)))
#define Vu(n) __attribute__((vector_size(n),aligned(1)))
#define cvt(x,t) __builtin_convertvector(x,t)
#define shfl1(x...) GC(__builtin_shuffle(x),__builtin_shufflevector(x))

#if defined(__clang__)
#  define shfl2(a,b,c...) __builtin_shufflevector(a,b,c)
#else // GCCs older than 12 don't support __builtin_shufflevector
#  define _F16(a,b,c,d,e,f,g,h,i,j,k,l,m,n,p,q,...) a,b,c,d,e,f,g,h,i,j,k,l,m,n,p,q
#  define _F8(a,b,c,d,e,f,g,h,...) a,b,c,d,e,f,g,h
#  define _F4(a,b,c,d,...) a,b,c,d
#  define _F2(a,b,...) a,b
#  define _cast_to_int(t,x...) _Generic((t), V8d: (V8q){_F8(x)}, V4d: (V4q){_F4(x)}, V2d: (V2q){_F2(x)}, \
                                             V16f: (V16i){_F16(x)}, V8f: (V8i){_F8(x)}, V4f: (V4i){_F4(x)})
#  define shfl2(a,b,c...) __builtin_shuffle(a,b,_cast_to_int(a,c,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1))
#endif
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
// complex numbers
typedef struct Vz{Vd re;Vd im;}Vz;
// Mask type that's returned from various comparison operators and that we ues for masked gather
typedef A(simde__mmask8,Vq,Vq) M8; // AVX512 uses bitmasks; everything else just uses vectors

// interleave/deinterleave for c128
#define unpack_Vz_idx0 A(E_(0,2,4,6,8,10,12,14),E_(0,2,4,6),E_(0,2))
#define unpack_Vz_idx1 A(E_(1,3,5,7,9,11,13,15),E_(1,3,5,7),E_(1,3))
#define pack_Vz_idx0 A(E_(0,8, 1,9, 2,10,3,11),E_(0,4,1,5),E_(0,2))
#define pack_Vz_idx1 A(E_(4,12,5,13,6,14,7,15),E_(2,6,3,7),E_(1,3))
D(Vz,unpack_Vz,_(Z2(shfl2(a,b,unpack_Vz_idx0),shfl2(a,b,unpack_Vz_idx1))),c(Vd)a,c(Vd)b)
D(Vz,pack_Vz,_(Z2(shfl2(a,b,pack_Vz_idx0),shfl2(a,b,pack_Vz_idx1))),c(Vd)a,c(Vd)b)

// interleave/deinterleave for c64
#if M <= 2
#  if defined(__clang__)
#    define lo_Vf(x) A(E_(shfl2(x,x,0,1,2,3,4,5,6,7)),E_(shfl2(x,x,0,1,2,3)),)
#    define hi_Vf(x) A(E_(shfl2(x,x,8,9,10,11,12,13,14,15)),E_(shfl2(x,x,4,5,6,7)),)
#    define from_lohi_Vf(lo,hi) A(E_(shfl2(lo,hi,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)),E_(shfl2(lo,hi,0,1,2,3,4,5,6,7)),)
#  else
#    define lo_Vf(x) A(I8(castps512_ps256,x),O(ps_ps256)(x),)
#    define hi_Vf(x) A(I8(extractf32x8_ps,x,1),O(vextractf128_ps256)(x,1),)
#    define from_lohi_Vf(lo,hi) A(E_(I8(insertf32x8,I8(castps256_ps512,lo),hi,1)),E_(O(vinsertf128_ps256)(O(ps256_ps)(lo),hi,1)),)
#  endif
#  define unpack_Vf_idx A(E_(0,2,4,6,8,10,12,14,1,3,5,7,9,11,13,15),E_(0,2,4,6,1,3,5,7),)
#  define pack_Vf_idx A(E_(0,8,1,9,2,10,3,11,4,12,5,13,6,14,7,15),E_(0,4,1,5,2,6,3,7),)
   D(Vz,unpack_Vf,_(c(Vf)b=shfl1(a,(Vi){unpack_Vf_idx});Z2(cvt(lo_Vf(b),Vd),cvt(hi_Vf(b),Vd))),c(Vf)a)
   D(Vf,pack_Vf,_(shfl1(from_lohi_Vf(cvt(a,A(V8f,V4f,)),cvt(b,A(V8f,V4f,))),(Vi){pack_Vf_idx})),c(Vd)a,c(Vd)b)
#else
   D(Vz,unpack_Vf,_(Z2((Vd){(f64)a[0],(f64)a[2]},(Vd){(f64)a[1],(f64)a[3]})),c(Vf)a)
   D(Vf,pack_Vf,_((Vf){(f32)a[0],(f32)b[0],(f32)a[1],(f32)b[1]}),c(Vd)a,c(Vd)b)
#endif

D(Vq,Si,((Vq){RN(x)}),c(i64)x)
D(Vd,Sd,((Vd){RN(x)}),c(f64)x)
De(u32,movemask,A(m,O(movmskpd256)((Vd)m),O(i64x2_bitmask)((simde_v128_t)m)),c(M8)m)
De(Vd,sqrtd,A(I8(sqrt_pd,x),O(sqrtpd256)(x),(Vd)O(f64x2_sqrt)((simde_v128_t)x)),c(Vd)x)
De(Vq,selectq,A(
  I8(mask_blend_epi64,s,a,b),
  E_((Vq)O(pblendvb256)((Vb)(a),(Vb)(b),(Vb)(s))),
  E_((Vq)O(v128_bitselect)((simde_v128_t)b,(simde_v128_t)a,O(i8x16_shr((simde_v128_t)s,7))))
),c(M8)s,c(Vq)b,c(Vq)a)
#define gt(a,b) A(I8(cmp_epi64_mask,b,a,1),((a)>(b)),((a)>(b)))
#define eq(a,b) A(I8(cmp_epi64_mask,b,a,0),((a)==(b)),((a)==(b)))
#define prefetch(p...) __builtin_prefetch(p)

#if M == 1
#  define Gq(p,i) I8(i64gather_epi64,i,p,8)
#  define _gthd(p,i,s) I8(i64gather_pd,i,p,s)
#  define _gthf(p,i,s) I8(i64gather_ps,i,p,s)
#  define _mgthd(p,i,m,s) I8(mask_i64gather_pd,Zd,m,i,p,s)
#  define _mgthf(p,i,m,s) I8(mask_i64gather_ps,I4(setzero_ps),m,i,p,s)
#  define _mgthi(p,i,m,s) I8(mask_i64gather_epi32,I4(setzero_si256),m,i,p,s)
#elif M == 2
#  define Gq(p,i) O(gatherdiv4di)(Zi,(c(i64)*)(p),i,((Vq){~0,~0,~0,~0}),8)
// #  define maski64_i32(m) shuffle2((Vi)m,(Vi)m,0,2,4,6)
#  define maski64_i32(m) __builtin_convertvector(m&0xFFFFFFFF,V4i)
#  define _mgthd(p,i,m,s) O(gatherdiv4df)(Zd,p,i,(Vd)(m),s)
#  define _mgthf(p,i,m,s) cvt(O(gatherdiv4sf256)(Z4f,p,i,(V4f)(maski64_i32(m)),s),Vd)
#  define _mgthi(p,i,m,s) O(gatherdiv4si256)(((V4i){0,0,0,0}),p,i,(V4i)(maski64_i32(m)),s)
#  define _gthd(p,i,s) _mgthd(p,i,Zd==Zd,s)
#  define _gthf(p,i,s) cvt(O(gatherdiv4sf256)(Z4f,p,i,Z4f==Z4f,s),Vd)
#else
// #   define maski64_i32(m) shuffle2((Vi)m,(Vi)m,0,2,4,6)
    D(Vq,Gq,((Vq){p[i[0]],p[i[1]]}),c(u64)*p,c(Vq)i)
    D(Vd,_gthd,_(c(u8)*p0=(c(u8)*)p+i[0]*s,*p1=(c(u8)*)p+i[1]*s;(Vd){*(c(f64)*)p0,*(c(f64)*)p1}),c(f64)*p,c(Vq)i,c(i32)s)
    D(Vd,_gthf,_(c(u8)*p0=(c(u8)*)p+i[0]*s,*p1=(c(u8)*)p+i[1]*s;(Vd){(f64)*(c(f32)*)p0,(f64)*(c(f32)*)p1}),c(f32)*p,c(Vq)i,c(i32)s)
    D(Vd,_mgthd,_(c(u8)*p0=(c(u8)*)p+i[0]*s,*p1=(c(u8)*)p+i[1]*s;(Vd){m[0]?*(c(f64)*)p0:0.0,m[1]?*(c(f64)*)p1:0.0}),c(f64)*p,c(Vq)i,c(M8)m,c(i32)s)
    D(Vd,_mgthf,_(c(u8)*p0=(c(u8)*)p+i[0]*s,*p1=(c(u8)*)p+i[1]*s;(Vd){m[0]?(f64)*(c(f32)*)p0:0.0,m[1]?(f64)*(c(f32)*)p1:0.0}),c(f32)*p,c(Vq)i,c(M8)m,c(i32)s)
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
   /*c128*/case 3: {c(Vd)a=Rd(p,0),b=Rd(p,1);r=unpack_Vz(a,b);break;}
    /*c64*/case 4: r=unpack_Vf(Rf(p,0));break;
   /*c32*/default: r=Z2(Zd,Zd);break;
};r),c(i32)t,c(void)*p)

// Wx---convert Vz to type t and write to p.
D(void,Wx,_(switch(t){
    /*f64*/case 0: Wd(p,0,z.re);break;
    /*f32*/case 1: {f32*_p=p;VL(u,N,_p[u]=(f32)z.re[u]);break;}
    /*f16*/case 2: break;
   /*c128*/case 3: {c(Vz)q=pack_Vz(z.re,z.im);Wd(p,0,q.re);Wd(p,1,q.im);break;}
    /*c64*/case 4: Wf(p,0,pack_Vf(z.re,z.im));break;
   /*c32*/default: break;
}),c(i32)t,void*p,c(Vz)z)

D(Vz,Gmx,_(Vz r;switch(t){
    /*f64*/case 0: r=Z2(_mgthd((c(f64)*)p,i,m,8),Zd);break;
    /*f32*/case 1: r=Z2(_mgthf((c(f32)*)p,i,m,4),Zd);break;
    /*f16*/case 2: r=Z2(Zd,Zd);break;
   /*c128*/case 3: r=Z2(_mgthd((c(f64)*)p,i<<1,m,8),_mgthd((c(f64)*)p+1,i<<1,m,8));break;
    /*c64*/case 4: r=unpack_Vf((Vf)_mgthd((c(f64)*)p,i,m,8));break;
   /*c32*/default: r=Z2(Zd,Zd);break;
};r),c(i32)t,c(void)*p,c(Vq)i,c(M8)m)

D(Vz,Gx,_(Vz r;switch(t){
     /*f64*/case 0: r=Z2(_gthd((c(f64)*)p,i,8),Zd);break;
     /*f32*/case 1: r=Z2(_gthf((c(f32)*)p,i,4),Zd);break;
     /*f16*/case 2: r=Z2(Zd,Zd);break;
    /*c128*/case 3: r=Z2(_gthd((c(f64)*)p,i<<1,8),_gthd((c(f64)*)p+1,i<<1,8));break;
     /*c64*/case 4: r=unpack_Vf((Vf)_gthd((c(f64)*)p,i,8));break;
    /*c32*/default: r=Z2(Zd,Zd);break;
};r),c(i32)t,c(void)*p,c(Vq)i)

// Rw2d---read N u16 numbers from p and convert them to f64, returning Vd
D(Vd,Rw2d,_(Vd r;VL(u,N,r[u]=p[u]);r),c(u16)*p)

// Wq2w---convert N i64 numbers (i.e., Vq) to u16 and write them to p. Receives 4 Vq as input and writes 1 Vw to p.
D(void,Wq2w,_(_L(k,4,VL(u,N,p[u+k*N]=(u16)x[k][u]))),u16*p,c(Vq)x[static 4])

De(Vd,Gmw,_(c(u32)_m=movemask(m); Vd r;VL(u,N,r[u]=((_m>>u)&1)?(f64)p[i[u]]:0.0);r),c(u16)*p,c(Vq)i,c(M8)m)
// __builtin_convertvector((A(V8i,V4i,V4i))shuffle2((Vhw)_mgthi((c(i32)*)p,i,m,2),Gmw_zero,Gmw_idx0),Vd)

#define _F_popcnt(i) __builtin_popcountll(x[i])
D(Vq,popcnt,_((Vq){RFN(_F_popcnt)}),c(Vq)x)
#undef _F_popcnt


D(Vz,mulz,Z2(a.re*b.re-a.im*b.im,a.im*b.re+a.re*b.im),c(Vz)a,c(Vz)b)
D(Vq,m1,A(I8(movm_epi64,I8(test_epi64_mask,x,m)),(x&m)!=Zi,(x&m)!=Zi),c(Vq)x,c(Vq)m)
D(Vq,m2,m1(x,m_0)^m1(x,m_1),c(Vq)x,c(Vq)m_0,c(Vq)m_1)
D(Vq,mX,popcnt(x&m),c(Vq)x,c(Vq)m)
D(Vz,bcast2,(Z2(Sd(re[k]),Sd(im[k]))),c(f64)*re,c(f64)*im,c(i32)k)
D(Vz,flipsign,_(m=m<<63;Z2((Vd)((Vq)v.re^m),(Vd)((Vq)v.im^m))),c(Vz)v,Vq m)

// TODO: The loop is ugly, but the compilers vectorize it
D(void,prefetchq,_(_L(k,N,prefetch(p+idx[k],0,3))),c(u64)*p,c(Vq)idx)

// TODO: Should get rid of these ...
// static void prefetchd4xN(Vi idx[4],c(u16)*n,c(f64)*x){_L(k,4*N,c(i32)i=((c(i64)*)idx)[k];simde_mm_prefetch(n+i,_MM_HINT_ET0);simde_mm_prefetch(x+i,_MM_HINT_ET0))}
// static void prefetchz4xN(Vi idx[4],c(u16)*n,c(c128)*x){_L(k,4*N,c(i32)i=((c(i64)*)idx)[k];simde_mm_prefetch(n+i,_MM_HINT_ET0);simde_mm_prefetch(x+i,_MM_HINT_ET0))}
