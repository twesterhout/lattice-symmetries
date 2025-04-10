#include "intrinsics.h"
#include "declarations.h"

Dd(zero,Zd,void);Dz(zero,(Z2(Zd,Zd)),void)
Dd(bcast2,Sd(re[k]),c(f64)*re,c(f64)*im,c(i32)k);Dz(bcast2,(Z2(Sd(re[k]),Sd(im[k]))),c(f64)*re,c(f64)*im,c(i32)k)

// Initializer for coeff##t
#define Ci(t,u,n) i32 r=0;Vi m[u];V##t acc[u];$(n>0,_(_L(_b,u,acc[_b]=bcast2##t(v_re,v_im,r));++r)){_L(_b,u,acc[_b]=zero##t())}
// Inner loop for coeff##t
#define Ck(t,u,b,x...) _L(k,b,_L(_b,u,m[_b]=(x));_L(_b,u,acc[_b]=add##t(acc[_b],signed##t(bcast2##t(v_re,v_im,r),m[_b])));++r)
static inline void coeffd1xN(c(Vi)x[1],Vd o[1],c(i32)n_s0,c(i32)n_s1,c(i32)n_s2,c(i32)n_sX,c(f64)*v_re,c(f64)*v_im,c(u64)*s1,c(u64)*s20,c(u64)*s21,c(u64)*sX){Ci(d,1,n_s0)Ck(d,1,n_s1,m1(x[_b],Si(s1[k])))Ck(d,1,n_s2,m2(x[_b],Si(s20[k]),Si(s21[k])))Ck(d,1,n_sX,mX(x[_b],Si(sX[k])))_L(_b,1,o[_b]=acc[_b])}
static inline void coeffd4xN(c(Vi)x[4],Vd o[4],c(i32)n_s0,c(i32)n_s1,c(i32)n_s2,c(i32)n_sX,c(f64)*v_re,c(f64)*v_im,c(u64)*s1,c(u64)*s20,c(u64)*s21,c(u64)*sX){Ci(d,4,n_s0)Ck(d,4,n_s1,m1(x[_b],Si(s1[k])))Ck(d,4,n_s2,m2(x[_b],Si(s20[k]),Si(s21[k])))Ck(d,4,n_sX,mX(x[_b],Si(sX[k])))_L(_b,4,o[_b]=acc[_b])}
static inline void coeffz1xN(c(Vi)x[1],Vz o[1],c(i32)n_s0,c(i32)n_s1,c(i32)n_s2,c(i32)n_sX,c(f64)*v_re,c(f64)*v_im,c(u64)*s1,c(u64)*s20,c(u64)*s21,c(u64)*sX){Ci(z,1,n_s0)Ck(z,1,n_s1,m1(x[_b],Si(s1[k])))Ck(z,1,n_s2,m2(x[_b],Si(s20[k]),Si(s21[k])))Ck(z,1,n_sX,mX(x[_b],Si(sX[k])))_L(_b,1,o[_b]=acc[_b])}
static inline void coeffz4xN(c(Vi)x[4],Vz o[4],c(i32)n_s0,c(i32)n_s1,c(i32)n_s2,c(i32)n_sX,c(f64)*v_re,c(f64)*v_im,c(u64)*s1,c(u64)*s20,c(u64)*s21,c(u64)*sX){Ci(z,4,n_s0)Ck(z,4,n_s1,m1(x[_b],Si(s1[k])))Ck(z,4,n_s2,m2(x[_b],Si(s20[k]),Si(s21[k])))Ck(z,4,n_sX,mX(x[_b],Si(sX[k])))_L(_b,4,o[_b]=acc[_b])}
#undef Ci
#undef Ck

static inline Vd coeffd(c(Vi)x,c(i32)n_s0,c(i32)n_s1,c(i32)n_s2,c(i32)n_sX,c(f64)*v_re,c(f64)*v_im,c(u64)*s1,c(u64)*s20,c(u64)*s21,c(u64)*sX){Vi xarr[1]={x}; Vd oarr[1];coeffd1xN(xarr,oarr,n_s0,n_s1,n_s2,n_sX,v_re,v_im,s1,s20,s21,sX);return oarr[0];}
static inline Vz coeffz(c(Vi)x,c(i32)n_s0,c(i32)n_s1,c(i32)n_s2,c(i32)n_sX,c(f64)*v_re,c(f64)*v_im,c(u64)*s1,c(u64)*s20,c(u64)*s21,c(u64)*sX){Vi xarr[1]={x}; Vz oarr[1];coeffz1xN(xarr,oarr,n_s0,n_s1,n_s2,n_sX,v_re,v_im,s1,s20,s21,sX);return oarr[0];}

static inline void diagd1xN(c(u64)*alpha0,c(void)*x0,void*out,c(oc_t)*ctx){Vi alpha[1];_L(_b,1,alpha[_b]=Ri(alpha0+_b*N));Vd acc[1];$(ctx->n_t>0,coeffd1xN(alpha,acc,ctx->n_s0[0],ctx->n_s1[0],ctx->n_s2[0],ctx->n_sX[0],ctx->v_re,ctx->v_im,ctx->s1,ctx->s20,ctx->s21,ctx->sX)){_L(_b,1,acc[_b]=zerod())};_L(_b,1,Wd((f64*)out+_b*N,muld(acc[_b],Rd((f64*)x0+_b*N))))}
static inline void diagz1xN(c(u64)*alpha0,c(void)*x0,void*out,c(oc_t)*ctx){Vi alpha[1];_L(_b,1,alpha[_b]=Ri(alpha0+_b*N));Vz acc[1];$(ctx->n_t>0,coeffz1xN(alpha,acc,ctx->n_s0[0],ctx->n_s1[0],ctx->n_s2[0],ctx->n_sX[0],ctx->v_re,ctx->v_im,ctx->s1,ctx->s20,ctx->s21,ctx->sX)){_L(_b,1,acc[_b]=zeroz())};_L(_b,1,Wz((c128*)out+_b*N,mulz(acc[_b],Rz((c128*)x0+_b*N))))}
De(void,diag64_f64, _(_L(k,64/N,diagd1xN(alpha0+k*N,( c(f64)*)x0+k*N, (f64*)out+k*N,ctx))),c(u64)*alpha0,c(void)*x0,void*out,c(oc_t)*ctx)
De(void,diag64_c128,_(_L(k,64/N,diagz1xN(alpha0+k*N,(c(c128)*)x0+k*N,(c128*)out+k*N,ctx))),c(u64)*alpha0,c(void)*x0,void*out,c(oc_t)*ctx)

D(Vi,pstep,_(assume(d<64);c(Vi)y=and(xor(shr(x,d),x),m);xor(xor(x,y),shl(y,d))),c(Vi)x,c(Vi)m,c(u32)d)

// the first row of ctx->masks is always the identity permutation
#define Rsetup(u,a,ctx) c(u64)*masks=ctx->masks+ctx->n_r;Vi b[u],rep[u],gid[u];_L(_u,u,rep[_u]=a[_u],gid[_u]=Zi)
#define Rs_(u,i,s) m=Si(masks[i]);_L(_u,u,b[_u]=pstep(b[_u],m,s))
#define Rpermute(u) _(_L(_u,u,b[_u]=a[_u]);Vi m;Rs_(u,0,1);Rs_(u,1,2);Rs_(u,2,4);Rs_(u,3,8);Rs_(u,4,16);Rs_(u,5,32);Rs_(u,6,16);Rs_(u,7,8);Rs_(u,8,4);Rs_(u,9,2);Rs_(u,10,1))
#define Rupdate(u,k) _(M8 p[u];c(Vi)vk=Si(k);_L(_u,u,p[_u]=gt(rep[_u],b[_u]))_L(_u,u,rep[_u]=select(p[_u],b[_u],rep[_u]),gid[_u]=select(p[_u],vk,gid[_u])))

static void repr4xN(c(Vi)a[4],c(bs_ctx_t)*ctx,Vi _rep[4],Vi _gid[4]){Rsetup(4,a,ctx);for(i64 k=1;k<ctx->n_m;++k,masks+=ctx->n_r){Rpermute(4);Rupdate(4,k);}_L(_u,4,_rep[_u]=rep[_u],_gid[_u]=gid[_u])}

#define Ssetup(u,x,ctx) i64 n=ctx->range_size;Vi j[u],v[u];_L(_u,u,j[_u]=gatherq(ctx->offsets,and(shr(x[_u],ctx->shift),Si(ctx->mask))))
#define Sgather(u,h,ctx) _L(_u,u,v[_u]=gatherq((c(i64)*)ctx->reps+h,j[_u]))
#define Supdate(u,h,x) _L(_u,u,j[_u]=select(gt(x[_u],v[_u]),addq(j[_u],Si(h)),j[_u]))

static void search4xN(c(Vi)x[4],M8 m[4],Vi i[4],c(search_ctx_t)*ctx){Ssetup(4,x,ctx);while(n>1){c(i64)h=n/2;Sgather(4,h,ctx);n-=h;Supdate(4,h,x)}Sgather(4,0,ctx);Supdate(4,1,x);Sgather(4,0,ctx);_L(_u,4,m[_u]=eqq(x[_u],v[_u]),i[_u]=j[_u])}

// TODO: ugly !!
INTERNAL void prefetchq(u64 const* p, Vi idx) { _Alignas(Vi) i64 buf[N]; Wi(buf,idx); for(i32 k=0;k<N;++k){simde_mm_prefetch(p+buf[k],_MM_HINT_T0);} }

static void repr_search4xN(c(Vi)a[4],c(Vi)r[4],c(bs_ctx_t)*bs_ctx,c(search_ctx_t)*search_ctx,Vi _rep[4],Vi _gid[4],M8 _msk[4],Vi _idx[4]){
    Ssetup(4,r,search_ctx);Rsetup(4,a,bs_ctx);
    for (i64 k = 1; k < bs_ctx->n_m;++k,masks+=bs_ctx->n_r) {
        $(n>1,_L(_u,4,prefetchq(search_ctx->reps+n/2,j[_u]))){}
        Rpermute(4);
        $(n>1,Sgather(4,n/2,search_ctx)){}
        Rupdate(4,k);
        $(n>1,c(i64)h=n/2;Supdate(4,h,r);n -= h){}
    }
    while(n>1){c(i64)h=n/2;Sgather(4,h,search_ctx);n-=h;Supdate(4,h,r)}
    Sgather(4,0,search_ctx);Supdate(4,1,r);Sgather(4,0,search_ctx)
    _L(_u,4,_msk[_u]=eqq(r[_u],v[_u]),_idx[_u]=j[_u]) // do this first in case the next line writes to r
    _L(_u,4,_rep[_u]=rep[_u],_gid[_u]=gid[_u])
}

static void w2d8xN(c(u16)*src,Vd _dst[4]){f64*dst=(f64*)_dst;_L(_u,8*N,dst[_u]=src[_u])}
// TODO: this loop is ugly ...
static void prefetchd4xN(Vi idx[4],c(u16)*n,c(f64)*x){_L(k,4*N,c(i32)i=((c(i64)*)idx)[k];simde_mm_prefetch(n+i,_MM_HINT_ET0);simde_mm_prefetch(x+i,_MM_HINT_ET0))}
static void prefetchz4xN(Vi idx[4],c(u16)*n,c(c128)*x){_L(k,4*N,c(i32)i=((c(i64)*)idx)[k];simde_mm_prefetch(n+i,_MM_HINT_ET0);simde_mm_prefetch(x+i,_MM_HINT_ET0))}


static Vd updatedF(Vd acc,c(Vd)x,c(Vd)outer){return fmad(acc,x,outer);}
static Vz updatezF(Vz acc,c(Vz)x,c(Vz)outer){return addz(mulz(acc,x),outer);}
static Vd updatedT(Vd acc,c(Vd)x,c(Vd)n,c(Vd)n0,c(Vd)chi,c(Vd)outer){return updatedF(muld(scaled(sqrtd(divd(n,n0)),chi),acc),x,outer);}
static Vz updatezT(Vz acc,c(Vz)x,c(Vd)n,c(Vd)n0,c(Vz)chi,c(Vz)outer){return updatezF(mulz(scalez(sqrtd(divd(n,n0)),chi),acc),x,outer);}

#define init(t) if(ctx->n_t==0){return;}c(f64)*v_re=ctx->v_re,*v_im=ctx->v_im;c(u64)*s1=ctx->s1,*s20=ctx->s20,*s21=ctx->s21,*sX=ctx->sX;Vi as[8],bs[4],rep[4],gid[4],msk[4],idx[4];Vd norm0[8],norm[4];V##t x[4],chi[4],acc[4],outer[8];w2d8xN(_norm,norm0);_L(_u,8,outer[_u]=zero##t())_L(_u,8,as[_u]=Ri(_alpha+_u*N))
#define increment(t) {s1+=ctx->n_s1[t],s20+=ctx->n_s2[t],s21+=ctx->n_s2[t],sX+=ctx->n_sX[t];c(i32)k=ctx->n_s0[t]+ctx->n_s1[t]+ctx->n_s2[t]+ctx->n_sX[t];v_re+=k,v_im+=k;}
static void off_diagd8xN(c(u64)*_alpha,c(u16)*_norm,c(f64)*_X,f64*out,c(oc_t)*ctx,c(bs_ctx_t)*bs_ctx,c(search_ctx_t)*search_ctx){
    init(d);_L(_u,4,bs[_u]=xor(as[_u],Si(ctx->mask[0])))$(bs_ctx!=NULL,repr4xN(bs,bs_ctx,rep,gid)){}
    for(i32 j=0;j<8*ctx->n_t-4;j+=4){
        c(i32)ti=j/8,bi=j%8;$(bs_ctx==NULL,_L(_u,4,x[_u]=gatherd(_X,bs[_u]))){}_L(_u,4,bs[_u]=xor(as[(j+4)%8+_u],Si(ctx->mask[(j+4)/8])))
        $(bs_ctx!=NULL,_L(_u,4,chi[_u]=gatherd(bs_ctx->chi_re, gid[_u]));repr_search4xN(bs,rep,bs_ctx,search_ctx,rep,gid,msk,idx);prefetchd4xN(idx,search_ctx->norm,_X)){}
        coeffd4xN(as+bi,acc,ctx->n_s0[ti],ctx->n_s1[ti],ctx->n_s2[ti],ctx->n_sX[ti],v_re,v_im,s1,s20,s21,sX);
        $(bs_ctx!=NULL,_L(_u,4,norm[_u]=mask_gather_norm(search_ctx->norm,idx[_u],msk[_u]),x[_u]=mask_gatherd(_X,idx[_u],msk[_u]))_L(_u,4,outer[bi+_u]=updatedT(acc[_u],x[_u],norm[_u],norm0[bi+_u],chi[_u],outer[bi+_u])))_L(_u,4,outer[bi+_u]=updatedF(acc[_u],x[_u],outer[bi+_u]))
        if(bi!=0){increment(ti);}
    }
    c(i32)ti=ctx->n_t-1;$(bs_ctx!=NULL,search4xN(rep,msk,idx,search_ctx);_L(_u,4,chi[_u]=gatherd(bs_ctx->chi_re, gid[_u]))_L(_u,4,norm[_u]=mask_gather_norm(search_ctx->norm,idx[_u],msk[_u]),x[_u]=mask_gatherd(_X,idx[_u],msk[_u])))_L(_u,4,x[_u]=gatherd(_X,bs[_u]))
    coeffd4xN(as+4,acc,ctx->n_s0[ti],ctx->n_s1[ti],ctx->n_s2[ti],ctx->n_sX[ti],v_re,v_im,s1,s20,s21,sX);
    $(bs_ctx!=NULL,_L(_u,4,outer[4+_u]=updatedT(acc[_u],x[_u],norm[_u],norm0[4+_u],chi[_u],outer[4+_u])))_L(_u,4,outer[4+_u]=updatedF(acc[_u],x[_u],outer[4+_u]))
    _L(_u,8,Wd(out+_u*N,addd(Rd(out+_u*N),outer[_u])))
}
static void off_diagz8xN(c(u64)*_alpha,c(u16)*_norm,c(c128)*_X,c128*out,c(oc_t)*ctx,c(bs_ctx_t)*bs_ctx,c(search_ctx_t)*search_ctx){
    init(z);_L(_u,4,bs[_u]=xor(as[_u],Si(ctx->mask[0])))$(bs_ctx!=NULL,repr4xN(bs,bs_ctx,rep,gid)){}
    for(i32 j=0;j<8*ctx->n_t-4;j+=4){
        c(i32)ti=j/8,bi=j%8;$(bs_ctx==NULL,_L(_u,4,x[_u]=gatherz(_X,bs[_u]))){}_L(_u,4,bs[_u]=xor(as[(j+4)%8+_u],Si(ctx->mask[(j+4)/8])))
        $(bs_ctx!=NULL,_L(_u,4,chi[_u]=gather2z(bs_ctx->chi_re,bs_ctx->chi_im,gid[_u]));repr_search4xN(bs,rep,bs_ctx,search_ctx,rep,gid,msk,idx);prefetchz4xN(idx,search_ctx->norm,_X)){}
        coeffz4xN(as+bi,acc,ctx->n_s0[ti],ctx->n_s1[ti],ctx->n_s2[ti],ctx->n_sX[ti],v_re,v_im,s1,s20,s21,sX);
        $(bs_ctx!=NULL,_L(_u,4,norm[_u]=mask_gather_norm(search_ctx->norm,idx[_u],msk[_u]),x[_u]=mask_gatherz(_X,idx[_u],msk[_u]))_L(_u,4,outer[bi+_u]=updatezT(acc[_u],x[_u],norm[_u],norm0[bi+_u],chi[_u],outer[bi+_u])))_L(_u,4,outer[bi+_u]=updatezF(acc[_u],x[_u],outer[bi+_u]))
        if(bi!=0){increment(ti);}
    }
    c(i32)ti=ctx->n_t-1;$(bs_ctx!=NULL,search4xN(rep,msk,idx,search_ctx);_L(_u,4,chi[_u]=gather2z(bs_ctx->chi_re,bs_ctx->chi_im,gid[_u]))_L(_u,4,norm[_u]=mask_gather_norm(search_ctx->norm,idx[_u],msk[_u]),x[_u]=mask_gatherz(_X,idx[_u],msk[_u])))_L(_u,4,x[_u]=gatherz(_X,bs[_u]))
    coeffz4xN(as+4,acc,ctx->n_s0[ti],ctx->n_s1[ti],ctx->n_s2[ti],ctx->n_sX[ti],v_re,v_im,s1,s20,s21,sX);
    $(bs_ctx!=NULL,_L(_u,4,outer[4+_u]=updatezT(acc[_u],x[_u],norm[_u],norm0[4+_u],chi[_u],outer[4+_u])))_L(_u,4,outer[4+_u]=updatezF(acc[_u],x[_u],outer[4+_u]))
    _L(_u,8,Wz(out+_u*N,addz(Rz(out+_u*N),outer[_u])))
}
#undef init
#undef increment

void off_diag64_f64(u64 const *alpha0, u16 const *norm0, void const *X, void *out,
        oc_t const *ctx, bs_ctx_t const *bs_ctx, search_ctx_t const* search_ctx) {
    for (i32 k = 0; k < 64; k += 8*N) { off_diagd8xN(alpha0 + k, norm0 + k, X, (f64*)out + k, ctx, bs_ctx, search_ctx); }
}
void off_diag64_c128(u64 const *alpha0, u16 const *norm0, void const *X, void *out,
        oc_t const *ctx, bs_ctx_t const *bs_ctx, search_ctx_t const* search_ctx) {
    for (i32 k = 0; k < 64; k += 8*N) { off_diagz8xN(alpha0 + k, norm0 + k, X, (c128*)out + k, ctx, bs_ctx, search_ctx); }
}

#define matvec_inner_template(t) \
    void matvec_inner_##t(i64 const i, u64 const *alpha0, u16 const *norm0, void const *X0, void const *X, void *out, \
            oc_t const *diag, oc_t const *off_diag, bs_ctx_t const *bs, search_ctx_t const *search) { \
        diag64_##t(alpha0 + i, (t const*)X0 + i, (t*)out + i, diag); off_diag64_##t(alpha0 + i, norm0 + i, X, (t*)out + i, off_diag, bs, search); \
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

// #define AP_F(v) flag |= OP(movemask_pd, i2d(GT(x, v)))
// #define AP_E(v) norm = addq(norm, SHR(EQ(x, v), 63))
Vi permute(Vi x, u64 const *masks, u32 const *shifts, i32 const n) { i32 i = 0; do { x = pstep(x, Si(masks[i]), shifts[i]); ++i; } while (i < n); return x; }
INTERNAL u32 ap_f(Vi const x, Vi const v) { return I(B,movemask_pd,i2d(gt(x, v))); }
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
