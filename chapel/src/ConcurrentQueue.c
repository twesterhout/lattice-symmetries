#include "ConcurrentQueue.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdatomic.h>
#include <assert.h>

static uint32_t nextPowerOfTwo(uint32_t const x) {
  assert(x > 0 && "expected a positive number");
  return x == 1 ? 1 : 1 << (32 - __builtin_clz(x - 1));
}

struct ConcurrentQueue {
    uint64_t _Atomic head;
    uint64_t _Atomic tail;
    uint64_t _Atomic* buffer;
    uint64_t capacity;
    uint64_t mask;
    uint64_t nullElement;
};

ConcurrentQueue* ls_chpl_create_ConcurrentQueue(uint32_t capacity, uint64_t nullElement) {
    capacity = nextPowerOfTwo(capacity);

    ConcurrentQueue* ptr = malloc(sizeof(ConcurrentQueue));
    atomic_store(&ptr->head, 0);
    atomic_store(&ptr->tail, 0);
    ptr->buffer = malloc(capacity * sizeof(uint64_t _Atomic));
    for (uint64_t i = 0; i < capacity; ++i) {
        atomic_store(&(ptr->buffer[i]), nullElement);
    }
    ptr->capacity = capacity;
    ptr->mask = capacity - 1;
    ptr->nullElement = nullElement;
    return ptr;
}

void ls_chpl_destroy_ConcurrentQueue(ConcurrentQueue* queue) {
    free(queue->buffer);
    free(queue);
}

int64_t ls_chpl_ConcurrentQueue_size(ConcurrentQueue const* queue) {
    uint64_t const localHead = atomic_load(&queue->head);
    uint64_t const localTail = atomic_load(&queue->tail);
    return localHead - localTail;
}

bool ls_chpl_ConcurrentQueue_try_push(ConcurrentQueue* queue, uint64_t const value) {
    uint64_t localHead;
    do {
        localHead = atomic_load(&queue->head);
        uint64_t const localTail = atomic_load(&queue->tail);
        if (localHead - localTail >= queue->capacity) { return false; }
    } while (atomic_load(&(queue->buffer[localHead & queue->mask])) != queue->nullElement
             || !atomic_compare_exchange_strong(&queue->head, &localHead, localHead + 1));

    uint64_t const old = atomic_exchange(&(queue->buffer[localHead & queue->mask]), value);
    return true;
}

void ls_chpl_ConcurrentQueue_push(ConcurrentQueue* queue, uint64_t const value) {
    while (!ls_chpl_ConcurrentQueue_try_push(queue, value)) {
        // nothing
    }
}

bool ls_chpl_ConcurrentQueue_try_pop(ConcurrentQueue* queue, uint64_t* value) {
    uint64_t localTail;
    do {
        localTail = atomic_load(&queue->tail);
        uint64_t const localHead = atomic_load(&queue->head);
        if (localTail >= localHead) { return false; }

    } while (atomic_load(&(queue->buffer[localTail & queue->mask])) == queue->nullElement
             || !atomic_compare_exchange_strong(&queue->tail, &localTail, localTail + 1));

    localTail = localTail & queue->mask;
    *value = atomic_exchange(&(queue->buffer[localTail]), queue->nullElement);
    return true;
}

void ls_chpl_ConcurrentQueue_pop(ConcurrentQueue* queue, uint64_t* value) {
    while (!ls_chpl_ConcurrentQueue_try_pop(queue, value)) {
        // nothing
    }
}
