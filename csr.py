import cffi, os, subprocess, tempfile


def _compile(unroll, temp_dir=None):
    temp = temp_dir or tempfile.mkdtemp(prefix="matvec-cache")
    ffi = cffi.FFI()
    ffi.cdef("typedef double _Complex t; typedef long long i64;void spmv(i64,i64 const*,i64 const*,t const*,i64,t const*,t*);")
    cc = os.getenv("CC", default="cc")
    flags = ["-fPIC", "-fopenmp", "-O3", "-ftree-vectorize", "-ffreestanding", "-fno-math-errno", "-ffast-math", "-DNDEBUG",
        "-march=native", "-mtune=native", "-Wall", "-Wextra", "-W", "-Wno-comment", "-Wno-unused-parameter", "-Wno-psabi"]
    with open(os.path.join(temp, "matvec.c"), "w") as f:
        f.write(f"""
        #define c(z) z const
        #define B {unroll}
        #define P _Pragma("omp for")
        #define U _Pragma("GCC unroll({unroll})")
        #define L(v,b,e,x...) for(i64 v=(b);v<(e);++v){{x;}}
        typedef double _Complex t; typedef long long i64;
        void spmv(c(i64)n,c(i64)*ii,c(i64)*jj,c(t)*dd,c(i64)s,c(t)*x,t*y){{
          P L(i,0,n,_Alignas(64)t c[B];U L(b,0,B,c[b]=0)L(j,ii[i],ii[i+1],U L(b,0,B,c[b]+=dd[j]*x[jj[j]*s+b]));U L(b,0,B,y[s*i+b]=c[b]))
        }}
        """)
    out = os.path.join(temp, "libmatvec.so")
    args = [cc, *flags, "-shared", "-o", out, os.path.join(temp, "matvec.c")]
    subprocess.run(args, check=True)
    return ffi.dlopen(out, ffi.RTLD_NOW | ffi.RTLD_LOCAL)

LIB = compile(unroll=2)

def batched_matvec(m, x, out=None):
    if not isinstance(m, scipy.sparse.csr_matrix): m = scipy.sparse.csr_matrix(m)
    if x.ndim < 2: x = x.reshape(-1, 1)
    if out is None: out = np.zeros_like(x)
    assert m.shape[1] == x.shape[0] and m.shape[0] == out.shape[0] and x.shape[1] == out.shape[1]
    assert m.dtype == x.dtype and x.dtype == out.dtype
    assert x.flags["C_CONTIGUOUS"] and out.flags["C_CONTIGUOUS"]
    LIB.spmv(m.shape[0], m.indptr, m.indices, m.data, x.shape[1], )
    pass
