// #include <immintrin.h>
// #include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h> // sysconf

#define MIN(a, b) ((a <= b) ? a : b)
#define PK i64 const P = page_size(); i64 const K = P / sizeof(u64)
#define a(t, n) ((t*)aligned_alloc(P, n * sizeof(t)))
#define r(p, t, n) ((t*)realloc(p, n * sizeof(t)))
#define f(p) do { if ((p) != NULL) free(p); } while (0)
typedef int64_t i64; typedef uint64_t u64; typedef uint16_t u16;
typedef void (*norm64_fn)(u64 const*, void const*, u16 *);
typedef void (*candidates64_fn)(u64, u64 *);
typedef void (*finish_alloc_fn)(u64);
typedef struct chunk_t { u64 *xs; u16 *ns; i64 cp; i64 sz; int ec; char padding[28]; } chunk_t;
_Static_assert(sizeof(chunk_t) == 64, "wrong padding");

static i64 page_size() { return sysconf(_SC_PAGESIZE); }
static i64 up(i64 const x, i64 const n) { return ((x + n - 1) / n) * n; }
static void reset(chunk_t* c) { f(c->xs); f(c->ns); c->cp = 0; c->sz = 0; }
void expand(chunk_t* c) {
  PK;
  if (c->cp <= 0) { c->cp = K; c->xs = a(u64, c->cp); c->ns = a(u16, c->cp); }
  else { c->cp = up(c->cp + c->cp / 4, K); c->xs = r(c->xs, u64, c->cp); c->ns = r(c->ns, u16, c->cp); }
  if (c->xs == NULL || c->ns == NULL) { reset(c); c->ec = -1; }
}
static void ap(chunk_t* c, u64 const alpha, u64 const norm) {
  if (c->sz >= c->cp) { expand(c); };
  if (c->ec != 0) return;
  c->xs[c->sz] = alpha; c->ns[c->sz] = norm; ++c->sz;
}

void candidates_simple(u64 x0, u64 *out) { ++x0; for (int k = 0; k < 64; ++k, ++x0) { out[k] = x0; } }

#define ITER(b) \
    candidates64(x0, xs); norm64(xs, ctx, ns); x0 = xs[(b) - 1]; \
    for (i64 k = 0; k < (b); ++k) { if (ns[k] != 0) { ap(&c, xs[k], ns[k]); } }
chunk_t one(i64 const n, u64 x0, candidates64_fn const candidates64, norm64_fn const norm64, void const* ctx) {
  chunk_t c = (chunk_t){.xs = NULL, .ns = NULL, .cp = 0, .sz = 0, .ec = 0};
  PK; u64 *xs = a(u64, 64); u16 *ns = a(u16, 64); i64 i = 0;
  for (; i < n - 64; i += 64) { ITER(64); if (c.ec != 0) { break; } }
  if (c.ec == 0 && i < n) { ITER(n - i); }
  if (c.ec != 0) { reset(&c); }
  f(xs); f(ns);
  return c;
}
#undef ITER

void* enumerate_states(i64 nc, i64 *sizes, u64 *starts, void *candidates64, void *norms64, void const* ctx, i64 *total_size) {
  PK; chunk_t *cs = a(chunk_t, nc); if (cs == NULL) { return NULL; }
  _Atomic int ec = 0;
  _Atomic i64 sz = 0;
#pragma omp parallel for schedule(dynamic, 1) \
    default(none) firstprivate(nc, cs, sizes, starts, candidates64, norms64, ctx) shared(ec, sz)
  for (i64 k = 0; k < nc; ++k) {
    if (__atomic_load_n(&ec, __ATOMIC_RELAXED) != 0) { continue; }
    cs[k] = one(sizes[k], starts[k], (candidates64_fn)candidates64, (norm64_fn)norms64, ctx);
    if (cs[k].ec != 0) { __atomic_store_n(&ec, cs[k].ec, __ATOMIC_RELAXED); }
    else { __atomic_fetch_add(&sz, cs[k].sz, __ATOMIC_RELAXED); }
  }
  if (ec != 0) { for (i64 k = 0; k < nc; ++k) { reset(cs + k); }; f(cs); return NULL; }
  *total_size = sz; return cs;
}
void copy_finalize(i64 nc, void* chunks, u64 *states, u16 *norms) {
  chunk_t *cs = chunks;
  for (i64 k = 0; k < nc; ++k) {
    __builtin_memcpy(states, cs[k].xs, cs[k].sz * sizeof(u64));
    __builtin_memcpy(norms, cs[k].ns, cs[k].sz * sizeof(u16));
    states += cs[k].sz; norms += cs[k].sz; reset(cs + k); 
  }
  f(cs);
}
