#pragma once

#include <stdint.h>
#include <stddef.h>

void ls_chpl_init(void);
void ls_chpl_finalize(void);
void ls_internal_qsort_int64(int64_t* xs, size_t count, uint64_t const* keys);
void ls_internal_func_and_parent(void** func, void** parent);
