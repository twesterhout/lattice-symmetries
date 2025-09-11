#include "intrinsics.h"
#include "declarations.h"

De(int,has_float16,0)

// Byte increment for different data types
D(i32,stride,_(c(i32)sizes[6]={sizeof(f64),sizeof(f32),/*sizeof(f16)*/2,sizeof(c128),sizeof(c64),/*sizeof(c32)*/4};sizes[t]),c(i32)t)

// Matrix elements.
//
// acc = 0
// for r in range(n_r):
//   sign = 1 if popcount(x & s[r]) % 2 == 0 else -1
//   acc += sign * v[r]
#define Ck(u,b,x...) _L(k,b,_L(_b,u,m[_b]=(x));_L(_b,u,c(Vz)z=flipsign(bcast2(c->v_re,c->v_im,r),m[_b]);acc[_b].re+=z.re;acc[_b].im+=z.im);++r)
#define Dcoeff(u) \
    D(void,coeffz##u##xN,_( \
        i32 r=i*c->stride;Vq m[u];Vz acc[u]; \
        c(Vz)acc0=(c->n_s0[i]>0)?bcast2(c->v_re,c->v_im,r++):Z2(Zd,Zd);_L(_b,u,acc[_b]=acc0) \
        Ck(u,c->n_s1[i],m1(x[_b],Si(c->s1[i*c->stride+k]))) \
        Ck(u,c->n_s2[i],m2(x[_b],Si(c->s20[i*c->stride+k]),Si(c->s21[i*c->stride+k]))) \
        Ck(u,c->n_sX[i],mX(x[_b],Si(c->sX[i*c->stride+k]))) \
        _L(_b,u,o[_b]=acc[_b]) \
    ),c(Vq)x[u],Vz o[u],c(i32)i,c(oc_t)*c)
Dcoeff(1)Dcoeff(4)
#undef Dcoeff
#undef Ck

// Diagonal coefficients
D(void,diag1xN,_(c(Vq)alpha=Rq(alpha0,0);Vz acc=Z2(Zd,Zd),x=Rx(t,x0);if(ctx->n_t>0){coeffz1xN(&alpha,&acc,0,ctx);}Wx(t,out,mulz(acc,x))),c(i32)t,c(u64)*alpha0,c(void)*x0,void*out,c(oc_t)*ctx)
De(void,diag64,_(c(i32)inc=stride(t);_L(k,64/N,diag1xN(t,alpha0+k*N,(c(u8)*)x0+k*inc*N,(u8*)out+k*inc*N,ctx))),c(i32)t,c(u64)*alpha0,c(void)*x0,void*out,c(oc_t)*ctx)

// Representatives
D(Vuq,pstep,_(c(Vuq)y=((x>>d)^x)&m;(x^y)^(y<<d)),c(Vuq)x,c(Vuq)m,c(u32)d)
// the first row of ctx->masks is always the identity permutation
#define Rsetup(u,a,ctx) c(u64)*masks=ctx->masks+ctx->n_r;Vq b[u],rep[u],gid[u];_L(_u,u,rep[_u]=a[_u],gid[_u]=Zi)
#define Rs_(u,i,s) m=Si(masks[i]);_L(_u,u,b[_u]=(Vq)pstep((Vuq)b[_u],(Vuq)m,s))
#define Rpermute(u) _(_L(_u,u,b[_u]=a[_u]);Vq m;Rs_(u,0,1);Rs_(u,1,2);Rs_(u,2,4);Rs_(u,3,8);Rs_(u,4,16);Rs_(u,5,32);Rs_(u,6,16);Rs_(u,7,8);Rs_(u,8,4);Rs_(u,9,2);Rs_(u,10,1))
#define Rupdate(u,k) _(M8 p[u];c(Vq)vk=Si(k);_L(_u,u,p[_u]=gt(rep[_u],b[_u]))_L(_u,u,rep[_u]=selectq(p[_u],b[_u],rep[_u]),gid[_u]=selectq(p[_u],vk,gid[_u])))
D(void,repr4xN,_(Rsetup(4,a,ctx);for(i64 k=1;k<ctx->n_m;++k,masks+=ctx->n_r){Rpermute(4);Rupdate(4,k);}_L(_u,4,_rep[_u]=rep[_u],_gid[_u]=gid[_u])),c(Vq)a[static 4],c(bs_ctx_t)*ctx,Vq _rep[static 4],Vq _gid[static 4])

// Binary search
#define Ssetup(u,x,ctx) i64 n=ctx->range_size;Vq j[u],v[u];_L(_u,u,j[_u]=Gq((c(u64)*)ctx->offsets,(x[_u]>>ctx->shift)&Si(ctx->mask)))
#define Sgather(u,h,ctx) _L(_u,u,v[_u]=Gq(ctx->reps+h,j[_u]))
#define Supdate(u,h,x) _L(_u,u,j[_u]=selectq(gt(x[_u],v[_u]),j[_u]+Si(h),j[_u]))
D(void,search4xN,_(Ssetup(4,x,ctx);while(n>1){c(i64)h=n/2;Sgather(4,h,ctx);n-=h;Supdate(4,h,x)}Sgather(4,0,ctx);Supdate(4,1,x);Sgather(4,0,ctx);_L(_u,4,m[_u]=eq(x[_u],v[_u]);i[_u]=j[_u])),c(Vq)x[4],M8 m[4],Vq i[4],c(search_ctx_t)*ctx)

// Fused repr and search---the bottleneck for simulations with symmetries
void repr_search4xN(c(Vq)a[4],c(Vq)r[4],c(bs_ctx_t)*bs_ctx,c(search_ctx_t)*search_ctx,Vq _rep[4],Vq _gid[4],M8 _msk[4],Vq _idx[4]){
    Ssetup(4,r,search_ctx);Rsetup(4,a,bs_ctx);
    c(i64)cnt=bs_ctx->n_m>search_ctx->steps?search_ctx->steps:bs_ctx->n_m;
    i64 k = 1;
    for(;k<cnt;++k,masks+=bs_ctx->n_r){
        _L(_u,4,prefetchq(search_ctx->reps+n/2,j[_u]))
        Rpermute(4);
        Sgather(4,n/2,search_ctx)
        Rupdate(4,k);
        c(i64)h=n/2;Supdate(4,h,r);n -= h;
    }
    for(;k<bs_ctx->n_m;++k,masks+=bs_ctx->n_r){Rpermute(4);Rupdate(4,k);}
    while(n>1){c(i64)h=n/2;Sgather(4,h,search_ctx);n-=h;Supdate(4,h,r)}
    Sgather(4,0,search_ctx);Supdate(4,1,r);Sgather(4,0,search_ctx)
    _L(_u,4,_msk[_u]=eq(r[_u],v[_u]);_idx[_u]=j[_u]) // do this first in case the next line writes to r
    _L(_u,4,_rep[_u]=rep[_u];_gid[_u]=gid[_u])
}

// Off-diagonal part with symmetries
De(void,beta4xN,_(U4(bs[u]=as[i+u]^Si(c->mask[j]))),c(Vq)as[static 4],Vq bs[static 4],c(i32)i,c(i32)j,c(oc_t)*c)
// D(void,_chid4xN,_(U4(chi[u]=gthd(c->chi_re,gid[u]))),c(Vq)gid[static 4],Vd chi[static 4],c(bs_ctx_t)*c)
De(void,chi4xN,_(U4(chi[u].re=_gthd(c->chi_re,gid[u],8),chi[u].im=_gthd(c->chi_im,gid[u],8))),c(Vq)gid[static 4],Vz chi[static 4],c(bs_ctx_t)*c)

// D(void,_uptd4xN,_(U4(
//     Vd x=mgthrx(0,_x,idx[u],msk[u]);
//     Vd n=mgthw2d(_n,idx[u],msk[u]);
//     Vd n0=Rw2d(_n0+(i+u)*N);
//     Vd c=muld(sqrtd(divd(n,n0)),muld(chi[u],acc[u]));
//     outer[i+u]=fmad(c,x,outer[i+u]);
// )),c(i32)t,c(Vq)idx[static 4],c(M8)msk[static 4],c(Vd)acc[static 4],c(Vd)chi[static 4],Vd outer[static 4],c(i32)i,c(f64)*_x,c(u16)*_n,c(u16)*_n0)
De(void,upt4xN,_(U4(
    Vz x=Gmx(t,_x,idx[u],msk[u]);
    Vd n=Gmw(_n,idx[u],msk[u]);
    Vd n0=Rw2d(_n0+(i+u)*N);
    Vd c1=sqrtd(n/n0);
    Vz c2=mulz(chi[u],acc[u]);
    Vz c3=mulz(x,c2);
    outer[i+u].re+=c1*c3.re;
    outer[i+u].im+=c1*c3.im;
)),c(i32)t,c(Vq)idx[static 4],c(M8)msk[static 4],c(Vz)acc[static 4],c(Vz)chi[static 4],Vz outer[static 4],c(i32)i,c(c128)*_x,c(u16)*_n,c(u16)*_n0)
De(void,odT8xN,_(
    Vq as[8],bs[4],rep[4],gid[4],idx[4];M8 msk[4];Vz chi[4],acc[4],outer[8];
    c(i32)inc=stride(t);U8(outer[u]=Rx(t,(c(u8)*)out+u*inc*N))U8(as[u]=Rq(_alpha,u))
    beta4xN(as,bs,0,0,ctx);repr4xN(bs,bs_ctx,rep,gid);
    for(i32 j=0;j<8*ctx->n_t-4;j+=4){
        beta4xN(as,bs,(j+4)%8,(j+4)/8,ctx);chi4xN(gid,chi,bs_ctx);
        repr_search4xN(bs,rep,bs_ctx,search_ctx,rep,gid,msk,idx);
        _L(u,4*N,prefetch(search_ctx->norm+idx[u/4][u%4],0,3);prefetch((c(u8)*)_X+inc*idx[u/4][u%4]))
        coeffz4xN(as+j%8,acc,j/8,ctx);upt4xN(t,idx,msk,acc,chi,outer,j%8,_X,search_ctx->norm,_norm);
    }
    search4xN(rep,msk,idx,search_ctx);chi4xN(gid,chi,bs_ctx);
    _L(u,4*N,prefetch(search_ctx->norm+idx[u/4][u%4],0,3);prefetch((c(u8)*)_X+inc*idx[u/4][u%4]))
    coeffz4xN(as+4,acc,ctx->n_t-1,ctx);upt4xN(t,idx,msk,acc,chi,outer,4,_X,search_ctx->norm,_norm);
    U8(Wx(t,(u8*)out+u*inc*N,outer[u]))
),c(i32)t,c(u64)*_alpha,c(u16)*_norm,c(void)*_X,void*out,c(oc_t)*ctx,c(bs_ctx_t)*bs_ctx,c(search_ctx_t)*search_ctx)
De(void,odF4xN,_(
    Vq as[4],bs[4];Vz acc[4],outer[4];
    c(i32)inc=stride(t);U4(outer[u]=Rx(t,(c(u8)*)out+u*inc*N))U4(as[u]=Rq(_alpha,u))
    _L(ti,ctx->n_t,
        beta4xN(as,bs,0,ti,ctx);/*prefetch##t##4xN(idx,search_ctx->norm,_X);*/
        coeffz4xN(as,acc,ti,ctx);U4(c(Vz)z=mulz(acc[u],Gx(t,_X,bs[u]));outer[u].re+=z.re;outer[u].im+=z.im)
    )
    U4(Wx(t,(u8*)out+u*inc*N,outer[u]))
),c(i32)t,c(u64)*_alpha,c(u16)*_norm,c(void)*_X,void*out,c(oc_t)*ctx,c(bs_ctx_t)*bs_ctx,c(search_ctx_t)*search_ctx)

De(void,off_diag64,_(
    if(ctx->n_t<=0)return;
    c(i32)inc=stride(t);c(typeof(&odT8xN))fn=bs_ctx!=0?odT8xN:odF4xN;c(i32)step=bs_ctx!=0?8*N:4*N;
    for(i32 k=0;k<64;k+=step){fn(t,alpha0+k,norm0+k,X,(u8*)out+k*inc,ctx,bs_ctx,search_ctx);}
),c(i32)t,c(u64)*alpha0,c(u16)*norm0,c(void)*X,void*out,c(oc_t)*ctx,c(bs_ctx_t)*bs_ctx,c(search_ctx_t)*search_ctx)

Di(void,matvec_inner,_(
   c(i32)inc=stride(t);
   diag64(t,alpha0+i,(c(u8)*)X0+i*inc,(u8*)out+i*inc,diag);
   off_diag64(t,alpha0+i,norm0+i,X,(u8*)out+i*inc,off_diag,bs,search)
),c(i32)t,c(i64)i,c(u64)*alpha0,c(u16)*norm0,c(void)*X0,c(void)*X,void*out,c(oc_t)*diag,c(oc_t)*off_diag,c(bs_ctx_t)*bs,c(search_ctx_t)*search)

void matvec(i32 const t, i64 const n, u64 const *alpha0, u16 const *norm0, void const *X0, void const *X, void *out,
        oc_t const *diag_ctx, oc_t const *off_diag_ctx, bs_ctx_t const *bs_ctx, search_ctx_t const *search_ctx) {
    if (n < 64) { return; }
    
    i64 const n_b = n / 64, n_r = n % 64;
#pragma omp parallel for schedule(dynamic, 256) default(none) \
        firstprivate(t, n_b, alpha0, norm0, X0, X, out, diag_ctx, off_diag_ctx, bs_ctx, search_ctx)
    for (i64 bi = 0; bi < n_b; ++bi) { matvec_inner(t, 64 * bi, alpha0, norm0, X0, X, out, diag_ctx, off_diag_ctx, bs_ctx, search_ctx); }
    if (n_r != 0) { matvec_inner(t, n - 64, alpha0, norm0, X0, X, out, diag_ctx, off_diag_ctx, bs_ctx, search_ctx); }
}

#define with_in(tmp_x, x, t, n, ...) t tmp_x[4*N]; __builtin_memset(tmp_x, 0, 4 * N * sizeof(t)); __builtin_memcpy(tmp_x, x, n * sizeof(t)); __VA_ARGS__
#define with_out(tmp_x, x, t, n, ...) t tmp_x[4*N]; __builtin_memset(tmp_x, 0, 4 * N * sizeof(t)); __VA_ARGS__; __builtin_memcpy(x, tmp_x, n * sizeof(t))

#define INNER(_k) {c(i64)k=(_k);Vq a[4],r[4],i[4];_L(_u,4,a[_u]=Rq(xs+k,_u));repr4xN(a,ctx,r,i);_L(_u,4,Wq(rep+k,_u,r[_u]),Wq(gid+k,_u,i[_u]))}
void state_info(c(i64)n,c(u64)*xs,c(bs_ctx_t)*ctx,u64*rep,i64*gid){
    if(n<=0){return;}
    if(n<4*N){with_in(tmp_xs,xs,u64,n,with_out(tmp_rep,rep,u64,n,with_out(tmp_gid,gid,i64,n,state_info(4*N,tmp_xs,ctx,tmp_rep,tmp_gid))));return;}

    c(i64)n_b=n/(4*N),n_r=n%(4*N);
#pragma omp parallel for default(none) firstprivate(n_b,xs,ctx,rep,gid)
    for(i64 bi=0;bi<n_b;++bi){INNER(4*N*bi);}
    if(n_r!=0){INNER(n-4*N);}
}
#undef INNER

#define INNER(_k) {c(i64)k=(_k);Vq a[4],idx[4];M8 msk[4];_L(_u,4,a[_u]=Rq(xs+k,_u));search4xN(a,msk,idx,ctx);_L(_u,4,Wq(out+k,_u,selectq(msk[_u],idx[_u],Si(-1))))}
void state_to_index(c(i64)n, c(u64)*xs, search_ctx_t const *ctx, i64 *out) {
    if(n<=0){return;}
    if(n<4*N){with_in(tmp_xs,xs,u64,n,with_out(tmp_out,out,i64,n,state_to_index(4*N,tmp_xs,ctx,tmp_out)));return;}

    c(i64)n_b=n/(4*N),n_r=n%(4*N);
#pragma omp parallel for default(none) firstprivate(n_b, xs, ctx, out)
    for(i64 bi=0;bi<n_b;++bi){INNER(4*N*bi);}
    if(n_r!=0){INNER(n-4*N);}
}
#undef INNER

// Vq permute(Vq x, u64 const *masks, u32 const *shifts, i32 const n) { i32 i = 0; do { x = pstep(x, Si(masks[i]), shifts[i]); ++i; } while (i < n); return x; }


D(u32,ap_f,movemask(gt(x,v)),c(Vq)x,c(Vq)v)
D(Vq,ap_e,norm+((x==v)&0x1),c(Vq)norm,c(Vq)x,c(Vq)v)
static inline Vq norm(Vq x, bs_ctx_t const* ctx) {
    int k = 1; unsigned flag = 0; Vq norm = Si(1)/*, m = Si(ctx->inversion_mask)*/;
    u8 const* flags = ctx->flags + k * 3; u64 const* masks = ctx->masks + k * ctx->n_r;
    Vq a[1]={x};Vq b[1];
    for (; k < ctx->n_m; ++k, masks += ctx->n_r, flags += 3) {
        u8 const /*use_f2 = flags[0], */use_e1 = flags[1]/*, use_e2 = flags[2]*/;
        Rpermute(1); c(Vq)y=b[0]; // permute(x, masks, ctx->shifts, ctx->n_r);
        flag |= ap_f(x, y);
        // if (use_f2) { y2 = XOR(y, m); AP_F(y2); }
        if (use_e1) { norm = ap_e(norm, x, y); }
        else { flag |= movemask(eq(x,y)); }
        if (flag == (1<<N)-1) { return Zi; }
        // if (use_e2) { AP_E(y2); }
    }
    c(Vq)c=A(((Vq){1,1<<1,1<<2,1<<3,1<<4,1<<5,1<<6,1<<7}),
             ((Vq){1,1<<1,1<<2,1<<3}),
             ((Vq){1,1<<1}));
    Vq const p = (Si(flag)&c)==c;
    return (~p)&norm;
}

void norm64(u64 const *alpha, bs_ctx_t const *ctx, u16 *out) {
    for(i32 k=0;k<64;k+=4*N) { Vq n[4];U4(n[u]=norm(Rq(alpha+k,u),ctx)); Wq2w(out+k,n); }
}


extern void*realloc(void*,unsigned long);
extern void free(void*);
#if __APPLE__
// MacOS doesn't support aligned_alloc until recent SDKs; We use aligned_alloc for performance only
void*aligned_alloc(unsigned long alignment,unsigned long size){return malloc(size);}
#else
extern void*aligned_alloc(unsigned long,unsigned long);
#endif

typedef struct chunk_t{u64*xs;u16*ns;i64 cp,sz;;int ec;char padding[28];}chunk_t;
_Static_assert(sizeof(chunk_t) == 64, "wrong padding");

typedef void (*candidates64_fn)(u64, u64 *);
De(void,candidates_simple,_(++x0;_L(k,64,out[k]=x0++)),u64 x0, u64*out)
De(void,candidates_hamming,_(_L(k,64,c(u64)t=x0|(x0-1);x0=(t+1)|(((~t&(t+1))-1)>>(__builtin_ctzll(x0)+1));out[k]=x0)),u64 x0, u64*out)

Di(i64,up,((x+n-1)/n)*n,c(i64)x,c(i64)n)
Di(void,reset,_(f(c->xs);f(c->ns);c->cp=0;c->sz=0;return),chunk_t*c)
Di(void,expand,_(
  if(c->cp<=0){c->cp=K;c->xs=a(u64,c->cp);c->ns=a(u16,c->cp);}
  else{c->cp=up(c->cp+c->cp/4,K);c->xs=r(c->xs,u64,c->cp);c->ns=r(c->ns,u16,c->cp);}
  if(c->xs==NULL||c->ns==NULL){reset(c);c->ec=-1;}
),chunk_t*c)
Di(void,ap,_(
  if(c->sz>=c->cp){expand(c);}
  if(c->ec!=0)return;
  c->xs[c->sz]=alpha;c->ns[c->sz]=norm;++c->sz;return
),chunk_t*c,c(u64)alpha,c(u64)norm)
#define ITER(b) \
    candidates64(x0, xs); norm64(xs, ctx, ns); x0 = xs[(b) - 1]; \
    for (i64 k = 0; k < (b); ++k) { if (ns[k] != 0) { ap(&c, xs[k], ns[k]); } }
Di(chunk_t,one,_(
  chunk_t c=(chunk_t){.xs=NULL,.ns=NULL,.cp=0,.sz=0,.ec=0};u64*xs=a(u64,64);u16*ns=a(u16,64);i64 i=0;
  for(;i<n-64;i+=64){ITER(64);if(c.ec!=0){break;}}
  if(c.ec==0&&i<n){ITER(n-i);}
  if(c.ec!=0){reset(&c);}
  f(xs);f(ns); c
),c(i64)n,u64 x0,c(candidates64_fn)candidates64,c(void)*ctx)
#undef ITER

void* enumerate_states(c(i64)nc,i64*sizes,u64*starts,void*candidates64,c(void)*ctx,i64*total_size){
  chunk_t*cs=a(chunk_t,nc);if(cs==NULL){return NULL;}
  int ec = 0; // atomic
  i64 sz = 0; // atomic
#pragma omp parallel for schedule(dynamic,1) \
    default(none) firstprivate(nc,cs,sizes,starts,candidates64,ctx) shared(ec,sz)
  for(i64 k=0;k<nc;++k){
    if(__atomic_load_n(&ec,__ATOMIC_RELAXED)!=0){continue;} // skip if error
    cs[k]=one(sizes[k],starts[k],(candidates64_fn)candidates64,ctx);
    if(cs[k].ec!=0){__atomic_store_n(&ec,cs[k].ec,__ATOMIC_RELAXED);}
    else{__atomic_fetch_add(&sz,cs[k].sz,__ATOMIC_RELAXED);}
  }
  if(ec!=0){for(i64 k=0;k<nc;++k){reset(cs+k);}f(cs);return NULL;}
  *total_size=sz;return cs;
}
void copy_finalize(c(i64)nc,void*chunks,u64*states,u16 *norms){
  chunk_t*cs=chunks;
  for (i64 k=0;k<nc;++k){
    __builtin_memcpy(states,cs[k].xs,cs[k].sz*sizeof(u64));
    __builtin_memcpy(norms,cs[k].ns,cs[k].sz*sizeof(u16));
    states+=cs[k].sz;norms+=cs[k].sz;reset(cs + k);
  }
  f(cs);
}
