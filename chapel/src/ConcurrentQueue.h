#pragma once

#include <stdbool.h>
#include <stdint.h>

typedef struct ConcurrentQueue ConcurrentQueue;

ConcurrentQueue* ls_chpl_create_ConcurrentQueue(uint32_t capacity, uint64_t nullElement);
void ls_chpl_destroy_ConcurrentQueue(ConcurrentQueue* queue);
bool ls_chpl_ConcurrentQueue_try_push(ConcurrentQueue* queue, uint64_t value);
void ls_chpl_ConcurrentQueue_push(ConcurrentQueue* queue, uint64_t value);
bool ls_chpl_ConcurrentQueue_try_pop(ConcurrentQueue* queue, uint64_t* value);
void ls_chpl_ConcurrentQueue_pop(ConcurrentQueue* queue, uint64_t* value);
