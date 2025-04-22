// This file should be included AFTER intrinsics.h
typedef char i8; typedef short i16; typedef int i32; typedef long long i64;
typedef unsigned char u8; typedef unsigned short u16; typedef unsigned u32; typedef unsigned long long u64;
typedef float f32; typedef double f64; typedef float _Complex c64; typedef double _Complex c128;

typedef struct oc_t {
    f64 const *v_re, *v_im;
    u64 const *s1, *s20, *s21, *sX;
    u64 const *mask;
    i32 const *n_s0, *n_s1, *n_s2, *n_sX;
    i32 n_t, stride;
} oc_t;

typedef struct bs_ctx_t {
    u64 const *masks;
    u32 const *shifts;
    u8 const *flags; // [use_f2, use_e1, use_e2]
    u64 const inversion_mask; // spin inversion mask
    f64 const *chi_re, *chi_im; // characters
    i32 n_m, n_r; // number group elements, number shifts
    char padding[8];
} bs_ctx_t;

typedef struct search_ctx_t {
    u64 const *reps;
    u16 const *norm;
    i64 const *offsets;
    i64 range_size;
    u64 mask;
    i32 shift;
    i32 steps;
    char padding[16];
} search_ctx_t;

int has_float16(void);
void norm64(u64 const *, bs_ctx_t const *, u16 *);
void diag64(i32, u64 const *, void const *, void *, const oc_t *);
void off_diag64(i32, u64 const *, u16 const *, void const *, void *, const oc_t *, const bs_ctx_t *, const search_ctx_t *);

void matvec(i32, i64, u64 const *, u16 const *, void const *, void const *, void *, const oc_t *, const oc_t *, const bs_ctx_t *, const search_ctx_t *);

void state_to_index(i64, u64 const *, search_ctx_t const *, i64 *);
void state_info(i64, u64 const *, bs_ctx_t const *, u64 *, i64 *);

void candidates_simple(u64, u64 *);
void* enumerate_states(i64, i64 *, u64 *, void *, void *, void const *, i64 *);
void copy_finalize(i64, void *, u64 *, u16 *);
