// This file should be included AFTER intrinsics.h
typedef int8_t i8; typedef int16_t i16; typedef int32_t i32; typedef int64_t i64;
typedef uint8_t u8; typedef uint16_t u16; typedef uint32_t u32; typedef uint64_t u64;
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
void norm64(uint64_t const *, bs_ctx_t const *, uint16_t *);
void diag64(i32, u64 const *, void const *, void *, const oc_t *);
void off_diag64(i32, u64 const *, u16 const *, void const *, void *, const oc_t *, const bs_ctx_t *, const search_ctx_t *);

void matvec(i32, i64, u64 const *, u16 const *, void const *, void const *, void *, const oc_t *, const oc_t *, const bs_ctx_t *, const search_ctx_t *);

void state_to_index(int64_t, uint64_t const *, search_ctx_t const *, int64_t *);
void state_info(int64_t, uint64_t const *, bs_ctx_t const *, uint64_t *, int64_t *);

void candidates_simple(uint64_t, uint64_t *);
void* enumerate_states(int64_t, int64_t *, uint64_t *, void *, void *, void const *, int64_t *);
void copy_finalize(int64_t, void *, uint64_t *, uint16_t *);
