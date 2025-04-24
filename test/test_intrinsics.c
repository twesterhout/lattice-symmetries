#include "intrinsics.h"
#include <assert.h>
#include <stdio.h>

#define C(a,b) __builtin_complex(a,b)

int main() {

  // Rx
  { f64  const as[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    f32  const bs[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    c(c128)* cs = (c(c128)*)as;
    c(c64)* ds = (c(c64)*)bs;
    {c(Vz)z=Rx(0,as); _L(u,N,assert(z.re[u] == u+1)) _L(u,N,assert(z.im[u] == 0)) }
    {c(Vz)z=Rx(1,bs); _L(u,N,assert(z.re[u] == u+1)) _L(u,N,assert(z.im[u] == 0)) }
    {c(Vz)z=Rx(3,cs); _L(u,N,assert(z.re[u] == 2*u+1)) _L(u,N,assert(z.im[u] == 2*u+2)) }
    {c(Vz)z=Rx(4,ds); _L(u,N,assert(z.re[u] == 2*u+1)) _L(u,N,assert(z.im[u] == 2*u+2)) }
  }

  // Wx
  {
    c(f64) as[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    c(Vz) z = Z2(*(c(Vd)*)as,*(c(Vd)*)(as+N));
    { f64 out[N]; Wx(0,out,z); _L(u,N,assert(out[u] == u+1)) }
    { f32 out[N]; Wx(1,out,z); _L(u,N,assert(out[u] == u+1)) }
    { c128 out[N]; Wx(3,out,z); _L(u,N,assert(__real__ out[u] == u+1)) _L(u,N,assert(__imag__ out[u] == u+1+N)) }
    { c64 out[N]; Wx(4,out,z); _L(u,N,assert(__real__ out[u] == u+1)) _L(u,N,assert(__imag__ out[u] == u+1+N)) }
  }

  // Gx
  {
    c(f64) as[20] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
    c(f32) bs[20] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
    c(c128)* cs = (c(c128)*)as;
    c(c64)* ds = (c(c64)*)bs;

    c(i64) _idx[8] = {0, 5, 3, 1, 8, 8, 2, 7};
    c(Vq) idx = *(c(Vq)*)_idx;

    c(f64) rs[8] = {1, 6, 4, 2, 9, 9, 3, 8};
    c(c128) rsz[8] = {C(1.0,2.0), C(11.0,12.0), C(7.0,8.0), C(3.0,4.0), C(17.0,18.0), C(17.0,18.0), C(5.0,6.0), C(15.0,16.0)};
    { c(Vz)z = Gx(0,as,idx); _L(u,N,assert(z.re[u] == rs[u])) _L(u,N,assert(z.im[u] == 0)) }
    { c(Vz)z = Gx(1,bs,idx); _L(u,N,assert(z.re[u] == rs[u])) _L(u,N,assert(z.im[u] == 0)) }
    { c(Vz)z = Gx(3,cs,idx); _L(u,N,assert(z.re[u] == __real__ rsz[u])) _L(u,N,assert(z.im[u] == __imag__ rsz[u])) }
    { c(Vz)z = Gx(4,ds,idx); _L(u,N,assert(z.re[u] == __real__ rsz[u])) _L(u,N,assert(z.im[u] == __imag__ rsz[u])) }
  }

  // Gmx
  {
    c(f64) as[20] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
    c(f32) bs[20] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
    c(c128)* cs = (c(c128)*)as;
    c(c64)* ds = (c(c64)*)bs;

    c(i64) _idx[8] = {0, 5, 3, 1, 8, 8, 2, 7};
    c(i64) _msk[8] = {0, -1, -1, 0, 0, -1, -1, -1};
    c(Vq) idx = *(c(Vq)*)_idx;
#if M == 1
    u32 msk = 0b11100110;
#else
    c(Vq) msk = *(c(Vq)*)_msk;
#endif

    c(f64) rs[8] = {0, 6, 4, 0, 0, 9, 3, 8};
    c(c128) rsz[8] = {C(0.0,0.0), C(11.0,12.0), C(7.0,8.0), C(0.0,0.0), C(0.0,0.0), C(17.0,18.0), C(5.0,6.0), C(15.0,16.0)};
    { c(Vz)z = Gmx(0,as,idx,msk); _L(u,N,assert(z.re[u] == rs[u])) _L(u,N,assert(z.im[u] == 0)) }
    { c(Vz)z = Gmx(1,bs,idx,msk); _L(u,N,assert(z.re[u] == rs[u])) _L(u,N,assert(z.im[u] == 0)) }
    { c(Vz)z = Gmx(3,cs,idx,msk); _L(u,N,assert(z.re[u] == __real__ rsz[u])) _L(u,N,assert(z.im[u] == __imag__ rsz[u])) }
    { c(Vz)z = Gmx(4,ds,idx,msk); _L(u,N,assert(z.re[u] == __real__ rsz[u])) _L(u,N,assert(z.im[u] == __imag__ rsz[u])) }
  }

  // Gmw
  {
    c(u16) as[20] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};

    c(i64) _idx[8] = {0, 5, 3, 1, 8, 8, 2, 7};
    c(i64) _msk[8] = {-1, -1, -1, 0, 0, -1, -1, -1};
    c(Vq) idx = *(c(Vq)*)_idx;
#if M == 1
    u32 msk = 0b11100111;
#else
    c(Vq) msk = *(c(Vq)*)_msk;
#endif

    c(f64) rs[8] = {1, 6, 4, 0, 0, 9, 3, 8};
    { c(Vd)z = Gmw(as,idx,msk); _L(u,N,fprintf(stderr, "%f ", z[u])); fprintf(stderr, "\n"); _L(u,N,assert(z[u] == rs[u])) }
  }


  // Rw2d
  {
    c(u16) as[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    c(Vd) z = Rw2d(as); _L(u,N,assert(z[u] == (f64)as[u]))
  }

  // Wq2w
  {
    i64 as[4*N]; u16 out[4*N]; c(Vq)*p=(c(Vq)*)as; _L(u,4*N,as[u]=u+1)
    Wq2w(out,p); _L(u,4*N,assert(out[u]==u+1))
  }

#if 0
  double const ds[8] = {1, 2, 3, 4, 5, 6, 7, 8};
  { simde__m512d_private a = simde__m512d_to_private(Rd8(ds)); _L(u,8,assert(a.f64[u] == ds[u] && "Rd8")); }
  { simde__m256d_private a = simde__m256d_to_private(Rd4(ds)); _L(u,4,assert(a.f64[u] == ds[u] && "Rd4")); }
  { simde__m128d_private a = simde__m128d_to_private(Rd2(ds)); _L(u,2,assert(a.f64[u] == ds[u] && "Rd2")); }
  float const fs[8] = {1, 2, 3, 4, 5, 6, 7, 8};
  { simde__m256_private a = simde__m256_to_private(Rf8(fs)); _L(u,8,assert(a.f32[u] == fs[u] && "Rf8")); }
  { simde__m128_private a = simde__m128_to_private(Rf4(fs)); _L(u,4,assert(a.f32[u] == fs[u] && "Rf4")); }
  { simde__m128_private a = simde__m128_to_private(Rf2(fs)); _L(u,2,assert(a.f32[u] == fs[u] && "Rf2")); }
  _Float16 const hs[8] = {1, 2, 3, 4, 5, 6, 7, 8};
  { simde__m128i_private a = simde__m128i_to_private(Rh8(hs)); _L(u,8,assert(a.f16[u] == fs[u] && "Rh8")); }
  { simde__m128i_private a = simde__m128i_to_private(Rh4(hs)); _L(u,4,assert(a.f16[u] == fs[u] && "Rh4")); }
  { simde__m128i_private a = simde__m128i_to_private(Rh2(hs)); _L(u,2,assert(a.f16[u] == fs[u] && "Rh2")); }
  // Writing
  { double ds[8]; Wd8(ds, simde_mm512_set_pd(8, 7, 6, 5, 4, 3, 2, 1)); _L(u,8,assert(ds[u] == (double)(u + 1) && "Wd8")); }
  { double ds[4]; Wd4(ds, simde_mm256_set_pd(4, 3, 2, 1)); _L(u,4,assert(ds[u] == (double)(u + 1) && "Wd4")); }
  { double ds[2]; Wd2(ds, simde_mm_set_pd(2, 1)); _L(u,2,assert(ds[u] == (double)(u + 1) && "Wd2")); }
  { float fs[8]; Wf8(fs, simde_mm256_set_ps(8, 7, 6, 5, 4, 3, 2, 1)); _L(u,8,assert(fs[u] == (float)(u + 1) && "Wf8")); }
  { float fs[4]; Wf4(fs, simde_mm_set_ps(4, 3, 2, 1)); _L(u,4,assert(fs[u] == (float)(u + 1) && "Wf4")); }
  { float fs[4]; __builtin_memset(fs, 0, 4 * sizeof(float)); Wf2(fs, simde_mm_set_ps(4, 3, 2, 1));
    _L(u,2,assert(fs[u] == (float)(u + 1) && "Wf2")) _L(u,2,assert(fs[2 + u] == 0.0f && "Wf2")) }
  { _Float16 fs[8]; simde__m128i_private x; _L(u,8, x.f16[u]=(_Float16)(u + 1)); Wh8(fs, simde__m128i_from_private(x));
    _L(u,8,assert(fs[u] == (_Float16)(u + 1) && "Wh8")); }
  { _Float16 fs[8]; __builtin_memset(fs, 0, 8 * sizeof(_Float16)); simde__m128i_private x; _L(u,8, x.f16[u]=(_Float16)(u + 1)); Wh4(fs, simde__m128i_from_private(x));
    _L(u,4,assert(fs[u] == (_Float16)(u + 1) && "Wh4")) _L(u,4,assert(fs[4 + u] == (_Float16)0 && "Wh4")); }
  { _Float16 fs[8]; __builtin_memset(fs, 0, 8 * sizeof(_Float16)); simde__m128i_private x; _L(u,8, x.f16[u]=(_Float16)(u + 1)); Wh2(fs, simde__m128i_from_private(x));
    _L(u,2,assert(fs[u] == (_Float16)(u + 1) && "Wh2")) _L(u,6,assert(fs[2 + u] == (_Float16)0 && "Wh2")) }
  // Gathering
  { float const fs[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    simde__m256_private a = simde__m256_to_private(mgthf8(fs, simde_mm512_set_epi64(0, 5, 3, 1, 8, 8, 2, 7), 0b10101111));
    assert(a.f32[0] == 8 && a.f32[1] == 3 && a.f32[2] == 9 && a.f32[3] == 9 && a.f32[4] == 0 && a.f32[5] == 4 && a.f32[6] == 0 && a.f32[7] == 1 && "mgthf8");
  }
  { _Float16 const fs[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    simde__m128i_private a = simde__m128i_to_private(mgthh8(fs, simde_mm512_set_epi64(0, 5, 3, 1, 8, 8, 2, 7), 0b11011001));
    assert(a.f16[0] == 8 && a.f16[1] == 0 && a.f16[2] == 0 && a.f16[3] == 9 && a.f16[4] == 2 && a.f16[5] == 0 && a.f16[6] == 6 && a.f16[7] == 1 && "mgthh8");
  }
  { _Float16 const fs[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    simde__m128i_private a = simde__m128i_to_private(mgthh4(fs, simde_mm256_set_epi64x(0, 5, 3, 1), simde_mm256_set_epi64x(-1, -1, 0, 0)));
    assert(a.f16[0] == 0 && a.f16[1] == 0 && a.f16[2] == 6 && a.f16[3] == 1 && "mgthh4");
  }
  { _Float16 const fs[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    simde__m128i_private a = simde__m128i_to_private(mgthh2(fs, simde_mm_set_epi64x(5, 3), simde_mm_set_epi64x(-1, -1)));
    assert(a.f16[0] == 4 && a.f16[1] == 6 && "mgthh2");
  }
#endif

#if 0
  // Gather real numbers
#if M == 1
  { double as[10]; float bs[10]; _Float16 cs[10];
    _L(u,10,as[u]=u+1) _L(u,10,bs[u]=u+1) _L(u,10,cs[u]=u+1)
    simde__m512i i = simde_mm512_set_epi64(0, 5, 3, 1, 8, 8, 2, 7);
    simde__mmask8 m = 0b10101111;
    { double out[8]; Wd8(out, mgthrx(0,as,i,m));
      assert(out[0] == 8 && out[1] == 3 && out[2] == 9 && out[3] == 9 && out[4] == 0 && out[5] == 4 && out[6] == 0 && out[7] == 1); }
    { double out[8]; Wd8(out, mgthrx(1,bs,i,m));
      assert(out[0] == 8 && out[1] == 3 && out[2] == 9 && out[3] == 9 && out[4] == 0 && out[5] == 4 && out[6] == 0 && out[7] == 1); }
    { double out[8]; Wd8(out, mgthrx(2,cs,i,m));
      assert(out[0] == 8 && out[1] == 3 && out[2] == 9 && out[3] == 9 && out[4] == 0 && out[5] == 4 && out[6] == 0 && out[7] == 1); }
  }
#elif M == 2
  { double as[10]; float bs[10]; _Float16 cs[10];
    _L(u,10,as[u]=u+1) _L(u,10,bs[u]=u+1) _L(u,10,cs[u]=u+1)
    simde__m256i i = simde_mm256_set_epi64x(0, 5, 3, 1);
    simde__m256i m = simde_mm256_set_epi64x(-1, -1, 0, 0);
    { double out[4]; Wd4(out, mgthrx(0,as,i,m)); assert(out[0] == 0 && out[1] == 0 && out[2] == 6 && out[3] == 1); }
    { double out[4]; Wd4(out, mgthrx(1,bs,i,m)); assert(out[0] == 0 && out[1] == 0 && out[2] == 6 && out[3] == 1); }
    { double out[4]; Wd4(out, mgthrx(2,cs,i,m)); assert(out[0] == 0 && out[1] == 0 && out[2] == 6 && out[3] == 1); }
  }
#elif M == 3
  { double as[10]; float bs[10]; _Float16 cs[10];
    _L(u,10,as[u]=u+1) _L(u,10,bs[u]=u+1) _L(u,10,cs[u]=u+1)
    simde__m128i i = simde_mm_set_epi64x(5, 3);
    simde__m128i m = simde_mm_set_epi64x(-1, -1);
    { double out[2]; Wd2(out, mgthrx(0,as,i,m)); assert(out[0] == 4 && out[1] == 6); }
    { double out[2]; Wd2(out, mgthrx(1,bs,i,m)); assert(out[0] == 4 && out[1] == 6); }
    { double out[2]; Wd2(out, mgthrx(2,cs,i,m)); assert(out[0] == 4 && out[1] == 6); }
  }
#endif

  // Gather complex numbers
#if M == 1
  { c128 as[10]; c64 bs[10]; c32 cs[10];
    _L(u,10,as[u]=__builtin_complex((f64)(2*u+1),(f64)(2*u+2)))
    _L(u,10,bs[u]=__builtin_complex((f32)(2*u+1),(f32)(2*u+2)))
    _L(u,10,cs[u]=__builtin_complex((f16)(2*u+1),(f16)(2*u+2)))
    simde__m512i i = simde_mm512_set_epi64(0, 5, 3, 1, 8, 8, 2, 7);
    simde__mmask8 m = 0b10101111;
    { double re[8],im[8]; Vz z=mgthzx(0,as,i,m); Wd8(re,z.re); Wd8(im,z.im);
      assert(re[0] == 15 && re[1] == 5 && re[2] == 17 && re[3] == 17 && re[4] == 0 && re[5] == 7 && re[6] == 0 && re[7] == 1);
      assert(im[0] == 16 && im[1] == 6 && im[2] == 18 && im[3] == 18 && im[4] == 0 && im[5] == 8 && im[6] == 0 && im[7] == 2); }
    { double re[8],im[8]; Vz z=mgthzx(1,bs,i,m); Wd8(re,z.re); Wd8(im,z.im);
      assert(re[0] == 15 && re[1] == 5 && re[2] == 17 && re[3] == 17 && re[4] == 0 && re[5] == 7 && re[6] == 0 && re[7] == 1);
      assert(im[0] == 16 && im[1] == 6 && im[2] == 18 && im[3] == 18 && im[4] == 0 && im[5] == 8 && im[6] == 0 && im[7] == 2); }
    { double re[8],im[8]; Vz z=mgthzx(2,cs,i,m); Wd8(re,z.re); Wd8(im,z.im);
      assert(re[0] == 15 && re[1] == 5 && re[2] == 17 && re[3] == 17 && re[4] == 0 && re[5] == 7 && re[6] == 0 && re[7] == 1);
      assert(im[0] == 16 && im[1] == 6 && im[2] == 18 && im[3] == 18 && im[4] == 0 && im[5] == 8 && im[6] == 0 && im[7] == 2); }
  }
#elif M == 2
  { c128 as[10]; c64 bs[10]; c32 cs[10];
    _L(u,10,as[u]=__builtin_complex((f64)(2*u+1),(f64)(2*u+2)))
    _L(u,10,bs[u]=__builtin_complex((f32)(2*u+1),(f32)(2*u+2)))
    _L(u,10,cs[u]=__builtin_complex((f16)(2*u+1),(f16)(2*u+2)))
    simde__m256i i = simde_mm256_set_epi64x(0, 5, 3, 1);
    simde__m256i m = simde_mm256_set_epi64x(-1, -1, 0, 0);
    { double re[4],im[4]; Vz z=mgthzx(0,as,i,m); Wd4(re,z.re); Wd4(im,z.im);
      assert(re[0] == 0 && re[1] == 0 && re[2] == 11 && re[3] == 1); assert(im[0] == 0 && im[1] == 0 && im[2] == 12 && im[3] == 2); }
    { double re[4],im[4]; Vz z=mgthzx(1,bs,i,m); Wd4(re,z.re); Wd4(im,z.im);
      assert(re[0] == 0 && re[1] == 0 && re[2] == 11 && re[3] == 1); assert(im[0] == 0 && im[1] == 0 && im[2] == 12 && im[3] == 2); }
    { double re[4],im[4]; Vz z=mgthzx(2,cs,i,m); Wd4(re,z.re); Wd4(im,z.im);
      assert(re[0] == 0 && re[1] == 0 && re[2] == 11 && re[3] == 1); assert(im[0] == 0 && im[1] == 0 && im[2] == 12 && im[3] == 2); }
  }
#else
  { c128 as[10]; c64 bs[10]; c32 cs[10];
    _L(u,10,as[u]=__builtin_complex((f64)(2*u+1),(f64)(2*u+2)))
    _L(u,10,bs[u]=__builtin_complex((f32)(2*u+1),(f32)(2*u+2)))
    _L(u,10,cs[u]=__builtin_complex((f16)(2*u+1),(f16)(2*u+2)))
    simde__m128i i = simde_mm_set_epi64x(5, 3);
    simde__m128i m = simde_mm_set_epi64x(-1, -1);
    { double re[2],im[2]; Vz z=mgthzx(0,as,i,m); Wd2(re,z.re); Wd2(im,z.im);
      assert(re[0] == 7 && re[1] == 11); assert(im[0] == 8 && im[1] == 12); }
    { double re[2],im[2]; Vz z=mgthzx(1,bs,i,m); Wd2(re,z.re); Wd2(im,z.im);
      assert(re[0] == 7 && re[1] == 11); assert(im[0] == 8 && im[1] == 12); }
    { double re[2],im[2]; Vz z=mgthzx(2,cs,i,m); Wd2(re,z.re); Wd2(im,z.im);
      assert(re[0] == 7 && re[1] == 11); assert(im[0] == 8 && im[1] == 12); }
  }
#endif


  // Writing complex numbers
#if M == 1
  { double _Complex as[8]; float _Complex bs[8]; _Float16 _Complex cs[8];
    Wzx(0,as,Z2(simde_mm512_set_pd(15, 13, 11, 9, 7, 5, 3, 1),simde_mm512_set_pd(16, 14, 12, 10, 8, 6, 4, 2)));
    _L(u,8,assert(__real__ as[u]==2*u+1)) _L(u,8,assert(__imag__ as[u]==2*u+2))
    Wzx(1,bs,Z2(simde_mm512_set_pd(15, 13, 11, 9, 7, 5, 3, 1),simde_mm512_set_pd(16, 14, 12, 10, 8, 6, 4, 2)));
    _L(u,8,assert(__real__ bs[u]==2*u+1)) _L(u,8,assert(__imag__ bs[u]==2*u+2))
    Wzx(2,cs,Z2(simde_mm512_set_pd(15, 13, 11, 9, 7, 5, 3, 1),simde_mm512_set_pd(16, 14, 12, 10, 8, 6, 4, 2)));
    _L(u,8,assert(__real__ cs[u]==2*u+1)) _L(u,8,assert(__imag__ cs[u]==2*u+2))
  }
#elif M == 2
  { double _Complex as[4]; float _Complex bs[4]; _Float16 _Complex cs[4];
    Wzx(0,as,Z2(simde_mm256_set_pd(7, 5, 3, 1),simde_mm256_set_pd(8, 6, 4, 2)));
    assert(__real__ as[0] == 1.0); assert(__imag__ as[0] == 2.0); assert(__real__ as[1] == 3.0); assert(__imag__ as[1] == 4.0);
    assert(__real__ as[2] == 5.0); assert(__imag__ as[2] == 6.0); assert(__real__ as[3] == 7.0); assert(__imag__ as[3] == 8.0);
    Wzx(1,bs,Z2(simde_mm256_set_pd(7, 5, 3, 1),simde_mm256_set_pd(8, 6, 4, 2)));
    // fprintf(stderr, "%f %f %f %f\n", __real__ bs[0], __imag__ bs[0], __real__ bs[1], __imag__ bs[1]);
    assert(__real__ bs[0] == 1.0f); assert(__imag__ bs[0] == 2.0f); assert(__real__ bs[1] == 3.0f); assert(__imag__ bs[1] == 4.0f);
    assert(__real__ bs[2] == 5.0f); assert(__imag__ bs[2] == 6.0f); assert(__real__ bs[3] == 7.0f); assert(__imag__ bs[3] == 8.0f);
    Wzx(2,cs,Z2(simde_mm256_set_pd(7, 5, 3, 1),simde_mm256_set_pd(8, 6, 4, 2)));
    assert(__real__ cs[0] == 1.0f); assert(__imag__ cs[0] == 2.0f); assert(__real__ cs[1] == 3.0f); assert(__imag__ cs[1] == 4.0f);
    assert(__real__ cs[2] == 5.0f); assert(__imag__ cs[2] == 6.0f); assert(__real__ cs[3] == 7.0f); assert(__imag__ cs[3] == 8.0f);
  }
#elif M == 3
  { double _Complex as[2]; float _Complex bs[2]; _Float16 _Complex cs[2];
    Wzx(0,as,Z2(simde_mm_set_pd(3, 1),simde_mm_set_pd(4, 2)));
    assert(__real__ as[0] == 1.0); assert(__imag__ as[0] == 2.0); assert(__real__ as[1] == 3.0); assert(__imag__ as[1] == 4.0);
    Wzx(1,bs,Z2(simde_mm_set_pd(3, 1),simde_mm_set_pd(4, 2)));
    // fprintf(stderr, "%f %f %f %f\n", __real__ bs[0], __imag__ bs[0], __real__ bs[1], __imag__ bs[1]);
    assert(__real__ bs[0] == 1.0f); assert(__imag__ bs[0] == 2.0f); assert(__real__ bs[1] == 3.0f); assert(__imag__ bs[1] == 4.0f);
    Wzx(2,cs,Z2(simde_mm_set_pd(3, 1),simde_mm_set_pd(4, 2)));
    assert(__real__ cs[0] == (_Float16)1); assert(__imag__ cs[0] == (_Float16)2); assert(__real__ cs[1] == (_Float16)3); assert(__imag__ cs[1] == (_Float16)4);
  }
#endif

  // Reading complex numbers
#if M == 1
  { double _Complex as[8]; float _Complex bs[8]; _Float16 _Complex cs[8];
    _L(u,8,as[u]=__builtin_complex((double)(2*u+1),(double)(2*u+2)))
    _L(u,8,bs[u]=__builtin_complex((float)(2*u+1),(float)(2*u+2)))
    _L(u,8,cs[u]=__builtin_complex((_Float16)(2*u+1),(_Float16)(2*u+2)))
    { Vz z = Rzx(0,as); double out[16]; Wd8(out,z.re);Wd8(out+8,z.im); _L(u,8,assert(out[u]==__real__ as[u])) _L(u,8,assert(out[8+u]==__imag__ as[u])) }
    { Vz z = Rzx(1,bs); double out[16]; Wd8(out,z.re);Wd8(out+8,z.im); _L(u,8,assert(out[u]==__real__ bs[u])) _L(u,8,assert(out[8+u]==__imag__ bs[u])) }
    { Vz z = Rzx(2,cs); double out[16]; Wd8(out,z.re);Wd8(out+8,z.im); _L(u,8,assert(out[u]==__real__ cs[u])) _L(u,8,assert(out[8+u]==__imag__ cs[u])) }
  }
#elif M == 2
  { double _Complex as[4]; float _Complex bs[4]; _Float16 _Complex cs[4];
    _L(u,4,as[u]=__builtin_complex((double)(2*u+1),(double)(2*u+2)))
    _L(u,4,bs[u]=__builtin_complex((float)(2*u+1),(float)(2*u+2)))
    _L(u,4,cs[u]=__builtin_complex((_Float16)(2*u+1),(_Float16)(2*u+2)))
    { Vz z = Rzx(0,as); double out[8]; Wd4(out,z.re);Wd4(out+4,z.im); _L(u,4,assert(out[u]==__real__ as[u])) _L(u,4,assert(out[4+u]==__imag__ as[u])) }
    { Vz z = Rzx(1,bs); double out[8]; Wd4(out,z.re);Wd4(out+4,z.im); _L(u,4,assert(out[u]==__real__ bs[u])) _L(u,4,assert(out[4+u]==__imag__ bs[u])) }
    { Vz z = Rzx(2,cs); double out[8]; Wd4(out,z.re);Wd4(out+4,z.im); _L(u,4,assert(out[u]==__real__ cs[u])) _L(u,4,assert(out[4+u]==__imag__ cs[u])) }
  }
#elif M == 3
  { double _Complex as[2]; float _Complex bs[2]; _Float16 _Complex cs[2];
    _L(u,2,as[u]=__builtin_complex((double)(2*u+1),(double)(2*u+2)))
    _L(u,2,bs[u]=__builtin_complex((float)(2*u+1),(float)(2*u+2)))
    _L(u,2,cs[u]=__builtin_complex((_Float16)(2*u+1),(_Float16)(2*u+2)))
    { Vz z = Rzx(0,as); double out[4]; Wd2(out,z.re);Wd2(out+2,z.im); _L(u,2,assert(out[u]==__real__ as[u])) _L(u,2,assert(out[2+u]==__imag__ as[u])) }
    { Vz z = Rzx(1,bs); double out[4]; Wd2(out,z.re);Wd2(out+2,z.im); _L(u,2,assert(out[u]==__real__ bs[u])) _L(u,2,assert(out[2+u]==__imag__ bs[u])) }
    { Vz z = Rzx(2,cs); double out[4]; Wd2(out,z.re);Wd2(out+2,z.im); _L(u,2,assert(out[u]==__real__ cs[u])) _L(u,2,assert(out[2+u]==__imag__ cs[u])) }
  }
#endif
#endif

  return 0;  
}
