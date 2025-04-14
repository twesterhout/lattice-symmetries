#include "intrinsics.h"
#include "declarations.h"

// Matrix elements
#define Ci(t,u) i32 r=i*c->stride;Vi m[u];V##t acc[u];$(c->n_s0[i]>0,_(_L(_b,u,acc[_b]=bcast2##t(c->v_re,c->v_im,r));++r)){_L(_b,u,acc[_b]=zero##t())}
#define Ck(t,u,b,x...) _L(k,b,_L(_b,u,m[_b]=(x));_L(_b,u,acc[_b]=add##t(acc[_b],signed##t(bcast2##t(c->v_re,c->v_im,r),m[_b])));++r)
#define Dcoeff(t,u) D(void,coeff##t##u##xN,_(Ci(t,u)Ck(t,u,c->n_s1[i],m1(x[_b],Si(c->s1[i*c->stride+k])))Ck(t,u,c->n_s2[i],m2(x[_b],Si(c->s20[i*c->stride+k]),Si(c->s21[i*c->stride+k])))Ck(t,u,c->n_sX[i],mX(x[_b],Si(c->sX[i*c->stride+k])))_L(_b,u,o[_b]=acc[_b])),c(Vi)x[u],V##t o[u],c(i32)i,c(oc_t)*c)
Dcoeff(d,1)Dcoeff(d,4)Dcoeff(z,1)Dcoeff(z,4)
#undef Dcoeff
#undef Ci
#undef Ck

// Diagonal coefficients
#define Ddiag(t,s,u) D(void,diag##t##u##xN,_(Vi alpha[u];_L(_b,u,alpha[_b]=Ri(alpha0+_b*N));V##t acc[u];$(ctx->n_t>0,coeff##t##u##xN(alpha,acc,0,ctx))_L(_b,u,acc[_b]=zero##t());_L(_b,u,W##t((s*)out+_b*N,mul##t(acc[_b],R##t((s*)x0+_b*N))))),c(u64)*alpha0,c(s)*x0,s*out,c(oc_t)*ctx)
#define Ddiag64(t,s) De(void,diag64_##s, _(_L(k,64/N,diag##t##1xN(alpha0+k*N,(c(s)*)x0+k*N,(s*)out+k*N,ctx))),c(u64)*alpha0,c(void)*x0,void*out,c(oc_t)*ctx)
Ddiag(d,f64,1)Ddiag(z,c128,1)Ddiag64(d,f64)Ddiag64(z,c128)
#undef Ddiag
#undef Ddiag64

// Representatives
D(Vi,pstep,_(c(Vi)y=and(xor(shr(x,d),x),m);xor(xor(x,y),shl(y,d))),c(Vi)x,c(Vi)m,c(u32)d)
// the first row of ctx->masks is always the identity permutation
#define Rsetup(u,a,ctx) c(u64)*masks=ctx->masks+ctx->n_r;Vi b[u],rep[u],gid[u];_L(_u,u,rep[_u]=a[_u],gid[_u]=Zi)
#define Rs_(u,i,s) m=Si(masks[i]);_L(_u,u,b[_u]=pstep(b[_u],m,s))
#define Rpermute(u) _(_L(_u,u,b[_u]=a[_u]);Vi m;Rs_(u,0,1);Rs_(u,1,2);Rs_(u,2,4);Rs_(u,3,8);Rs_(u,4,16);Rs_(u,5,32);Rs_(u,6,16);Rs_(u,7,8);Rs_(u,8,4);Rs_(u,9,2);Rs_(u,10,1))
#define Rupdate(u,k) _(M8 p[u];c(Vi)vk=Si(k);_L(_u,u,p[_u]=gt(rep[_u],b[_u]))_L(_u,u,rep[_u]=select(p[_u],b[_u],rep[_u]),gid[_u]=select(p[_u],vk,gid[_u])))
D(void,repr4xN,_(Rsetup(4,a,ctx);for(i64 k=1;k<ctx->n_m;++k,masks+=ctx->n_r){Rpermute(4);Rupdate(4,k);}_L(_u,4,_rep[_u]=rep[_u],_gid[_u]=gid[_u])),c(Vi)a[4],c(bs_ctx_t)*ctx,Vi _rep[4],Vi _gid[4])

// Binary search
#define Ssetup(u,x,ctx) i64 n=ctx->range_size;Vi j[u],v[u];_L(_u,u,j[_u]=gthq((c(u64)*)ctx->offsets,and(shr(x[_u],ctx->shift),Si(ctx->mask))))
#define Sgather(u,h,ctx) _L(_u,u,v[_u]=gthq(ctx->reps+h,j[_u]))
#define Supdate(u,h,x) _L(_u,u,j[_u]=select(gt(x[_u],v[_u]),addq(j[_u],Si(h)),j[_u]))
D(void,search4xN,_(Ssetup(4,x,ctx);while(n>1){c(i64)h=n/2;Sgather(4,h,ctx);n-=h;Supdate(4,h,x)}Sgather(4,0,ctx);Supdate(4,1,x);Sgather(4,0,ctx);_L(_u,4,m[_u]=eqq(x[_u],v[_u]),i[_u]=j[_u])),c(Vi)x[4],M8 m[4],Vi i[4],c(search_ctx_t)*ctx)

// Fused repr and search---the bottleneck for simulations with symmetries
static void repr_search4xN(c(Vi)a[4],c(Vi)r[4],c(bs_ctx_t)*bs_ctx,c(search_ctx_t)*search_ctx,Vi _rep[4],Vi _gid[4],M8 _msk[4],Vi _idx[4]){
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
    _L(_u,4,_msk[_u]=eqq(r[_u],v[_u]),_idx[_u]=j[_u]) // do this first in case the next line writes to r
    _L(_u,4,_rep[_u]=rep[_u],_gid[_u]=gid[_u])
}

// Off-diagonal part with symmetries
D(void,_beta4xN,_(U4(bs[u]=xor(as[i+u],Si(c->mask[j])))),c(Vi)as[static 4],Vi bs[static 4],c(i32)i,c(i32)j,c(oc_t)*c)
D(void,_chid4xN,_(U4(chi[u]=gthd(c->chi_re,gid[u]))),c(Vi)gid[static 4],Vd chi[static 4],c(bs_ctx_t)*c)
D(void,_chiz4xN,_(U4(chi[u].re=gthd(c->chi_re,gid[u]),chi[u].im=gthd(c->chi_im,gid[u]))),c(Vi)gid[static 4],Vz chi[static 4],c(bs_ctx_t)*c)
D(void,_uptd4xN,_(U4(
    Vd x=mgthd(_x,idx[u],msk[u]);
    Vd n=mgthw2d(_n,idx[u],msk[u]);
    Vd n0=Rw2d(_n0+(i+u)*N);
    Vd c=muld(sqrtd(divd(n,n0)),muld(chi[u],acc[u]));
    outer[i+u]=fmad(c,x,outer[i+u]);
)),c(Vi)idx[static 4],c(M8)msk[static 4],c(Vd)acc[static 4],c(Vd)chi[static 4],Vd outer[static 4],c(i32)i,c(f64)*_x,c(u16)*_n,c(u16)*_n0)
D(void,_uptz4xN,_(U4(
    Vz x=mgthz(_x,idx[u],msk[u]);
    Vd n=mgthw2d(_n,idx[u],msk[u]);
    Vd n0=Rw2d(_n0+(i+u)*N);
    Vd c1=sqrtd(divd(n,n0));
    Vz c2=mulz(chi[u],acc[u]);
    Vz c3=mulz(x,c2);
    outer[i+u].re=fmad(c1,c3.re,outer[i+u].re);
    outer[i+u].im=fmad(c1,c3.im,outer[i+u].im);
)),c(Vi)idx[static 4],c(M8)msk[static 4],c(Vz)acc[static 4],c(Vz)chi[static 4],Vz outer[static 4],c(i32)i,c(c128)*_x,c(u16)*_n,c(u16)*_n0)
#define DodT8xN(t,s) \
    static void odT##t##8xN(c(u64)*_alpha,c(u16)*_norm,c(s)*_X,s*out,c(oc_t)*ctx,c(bs_ctx_t)*bs_ctx,c(search_ctx_t)*search_ctx){ \
        Vi as[8],bs[4],rep[4],gid[4],idx[4];M8 msk[4];V##t chi[4],acc[4],outer[8]; \
        U8(outer[u]=zero##t()) \
        U8(as[u]=Ri(_alpha+u*N)) \
        _beta4xN(as,bs,0,0,ctx); \
        repr4xN(bs,bs_ctx,rep,gid); \
        for(i32 j=0;j<8*ctx->n_t-4;j+=4){ \
            _beta4xN(as,bs,(j+4)%8,(j+4)/8,ctx); \
            _chi##t##4xN(gid,chi,bs_ctx); \
            repr_search4xN(bs,rep,bs_ctx,search_ctx,rep,gid,msk,idx); \
            /*prefetch##t##4xN(idx,search_ctx->norm,_X);*/ \
            c(i32)ti=j/8,bi=j%8; \
            coeff##t##4xN(as+bi,acc,ti,ctx); \
            _upt##t##4xN(idx,msk,acc,chi,outer,bi,_X,search_ctx->norm,_norm); \
        } \
        search4xN(rep,msk,idx,search_ctx); \
        _chi##t##4xN(gid,chi,bs_ctx); \
        /*prefetch##t##4xN(idx,search_ctx->norm,_X);*/ \
        c(i32)ti=ctx->n_t-1,bi=4; \
        coeff##t##4xN(as+bi,acc,ti,ctx); \
        _upt##t##4xN(idx,msk,acc,chi,outer,bi,_X,search_ctx->norm,_norm); \
        U8(W##t(out+u*N,add##t(R##t(out+u*N),outer[u]))) \
    }
DodT8xN(d,f64)DodT8xN(z,c128)

// Off-diagonal part without symmetries
#define DodF4xN(t,s) \
    static void odF##t##4xN(c(u64)*_alpha,c(u16)*_norm,c(s)*_X,s*out,c(oc_t)*ctx,c(bs_ctx_t)*bs_ctx,c(search_ctx_t)*search_ctx){ \
        Vi as[4],bs[4];V##t acc[4],outer[4]; \
        U4(outer[u]=zero##t())U4(as[u]=Ri(_alpha+u*N)) \
        _L(ti,ctx->n_t, \
            _beta4xN(as,bs,0,ti,ctx); \
            /*prefetch##t##4xN(idx,search_ctx->norm,_X);*/ \
            coeff##t##4xN(as,acc,ti,ctx); \
            U4(outer[u]=add##t(mul##t(acc[u],gth##t(_X,bs[u])),outer[u])) \
        ) \
        U4(W##t(out+u*N,add##t(R##t(out+u*N),outer[u]))) \
    }
DodF4xN(d,f64)DodF4xN(z,c128)

#define Doff_diag64(t,s) De(void,off_diag64_##s,_( \
    if(ctx->n_t<=0)return;i32 k=0;while(k<64){ \
        if(bs_ctx!=NULL){odT##t##8xN(alpha0+k,norm0+k,X,(s*)out+k,ctx,bs_ctx,search_ctx);k+=8*N;} \
        else{odF##t##4xN(alpha0+k,norm0+k,X,(s*)out+k,ctx,bs_ctx,search_ctx);k+=4*N;} \
    }),c(u64)*alpha0,c(u16)*norm0,c(void)*X,void*out,c(oc_t)*ctx,c(bs_ctx_t)*bs_ctx,c(search_ctx_t)*search_ctx)
Doff_diag64(d,f64)Doff_diag64(z,c128)
#undef Doff_diag64

// Matvec for 64 elements
#define matvec_inner_template(t) \
    void matvec_inner_##t(i64 const i, u64 const *alpha0, u16 const *norm0, void const *X0, void const *X, void *out, \
            oc_t const *diag, oc_t const *off_diag, bs_ctx_t const *bs, search_ctx_t const *search) { \
        diag64_##t(alpha0 + i, (t const*)X0 + i, (t*)out + i, diag); off_diag64_##t(alpha0 + i, norm0 + i, X, (t*)out + i, off_diag, bs, search); \
    }
matvec_inner_template(f64)matvec_inner_template(c128)
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


#define with_in(tmp_x, x, t, n, ...) t tmp_x[4*N]; __builtin_memset(tmp_x, 0, 4 * N * sizeof(t)); __builtin_memcpy(tmp_x, x, n * sizeof(t)); __VA_ARGS__
#define with_out(tmp_x, x, t, n, ...) t tmp_x[4*N]; __builtin_memset(tmp_x, 0, 4 * N * sizeof(t)); __VA_ARGS__; __builtin_memcpy(x, tmp_x, n * sizeof(t))

#define INNER(_k) {c(i64)k=(_k);Vi a[4],r[4],i[4];_L(_u,4,a[_u]=Ri(xs+k+_u*N));repr4xN(a,ctx,r,i);_L(_u,4,Wi(rep+k+_u*N,r[_u]),Wi(gid+k+_u*N,i[_u]))}
void state_info(c(i64)n,c(u64)*xs,c(bs_ctx_t)*ctx,u64*rep,i64*gid){
    if(n<=0){return;}
    if(n<4*N){with_in(tmp_xs,xs,u64,n,with_out(tmp_rep,rep,u64,n,with_out(tmp_gid,gid,i64,n,state_info(4*N,tmp_xs,ctx,tmp_rep,tmp_gid))));return;}

    c(i64)n_b=n/(4*N),n_r=n%(4*N);
#pragma omp parallel for default(none) firstprivate(n_b,xs,ctx,rep,gid)
    for(i64 bi=0;bi<n_b;++bi){INNER(4*N*bi);}
    if(n_r!=0){INNER(n-4*N);}
}
#undef INNER

#define INNER(_k) {c(i64)k=(_k);Vi a[4],idx[4];M8 msk[4];_L(_u,4,a[_u]=Ri(xs+k+_u*N));search4xN(a,msk,idx,ctx);_L(_u,4,Wi(out+k+_u*N,select(msk[_u],idx[_u],Si(-1))))}
void state_to_index(int64_t const n, uint64_t const *xs, search_ctx_t const *ctx, int64_t *out) {
    if(n<=0){return;}
    if(n<4*N){with_in(tmp_xs,xs,u64,n,with_out(tmp_out,out,i64,n,state_to_index(4*N,tmp_xs,ctx,tmp_out)));return;}

    c(i64)n_b=n/(4*N),n_r=n%(4*N);
#pragma omp parallel for default(none) firstprivate(n_b, xs, ctx, out)
    for(i64 bi=0;bi<n_b;++bi){INNER(4*N*bi);}
    if(n_r!=0){INNER(n-4*N);}
}
#undef INNER

Vi permute(Vi x, u64 const *masks, u32 const *shifts, i32 const n) { i32 i = 0; do { x = pstep(x, Si(masks[i]), shifts[i]); ++i; } while (i < n); return x; }

#if M == 1
D(u32,movemask,m,c(M8)m)
#else
D(u32,movemask,I(movemask_pd,i2d(m)),c(M8)m)
#endif

INTERNAL u32 ap_f(Vi const x, Vi const v) { return movemask(gt(x, v)); }
INTERNAL Vi ap_e(Vi const norm, Vi const x, Vi const v) { return addq(norm, shr(eqi(x, v), 63)); }
static HEDLEY_ALWAYS_INLINE Vi norm(Vi x, bs_ctx_t const* ctx) {
    int k = 1; unsigned flag = 0; Vi norm = Si(1)/*, m = Si(ctx->inversion_mask)*/;
    uint8_t const* flags = ctx->flags + k * 3; uint64_t const* masks = ctx->masks + k * ctx->n_r;
    for (; k < ctx->n_m; ++k, masks += ctx->n_r, flags += 3) {
        uint8_t const /*use_f2 = flags[0], */use_e1 = flags[1]/*, use_e2 = flags[2]*/;
        Vi y = permute(x, masks, ctx->shifts, ctx->n_r);
        flag |= ap_f(x, y);
        // if (use_f2) { y2 = XOR(y, m); AP_F(y2); }
        if (use_e1) { norm = ap_e(norm, x, y); }
        else { flag |= movemask(eqq(x, y)); }
        if (flag == (1<<N)-1) { return Zi; }
        // if (use_e2) { AP_E(y2); }
    }
#if M == 1
    Vi const c = I(set_epi64,1<<7,1<<6,1<<5,1<<4,1<<3,1<<2,1<<1,1);
#elif M == 2
    Vi const c = I(set_epi64x,0b1000,0b100,0b10,0b1);
#else
    Vi const c = I(set_epi64x,0b10,0b1);
#endif
    Vi const p = eqi(and(Si(flag), c), c);
    return andnot(p, norm);
}

void norm64(uint64_t const *alpha, bs_ctx_t const *ctx, uint16_t *out) {
    for(i32 k=0;k<64;k+=4*N) { Vi n[4];U4(n[u]=norm(Ri(alpha+k+u*N),ctx)); Wq2w(out+k,n); }
}
