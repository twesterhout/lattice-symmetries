#pragma once

#if defined(__cplusplus)
#include <cstdint>
#else
#include <stdint.h>
#endif

#if defined(__cplusplus)
extern "C" {
#endif

void ls_hs_init(void);
void ls_hs_exit(void);
int ls_hs_foo(void);

#ifdef __cplusplus
} // extern "C"
#endif
