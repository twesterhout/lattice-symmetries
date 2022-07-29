#include "lattice_symmetries_haskell.h"

#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#include <stdio.h>

#define RING_BUFFER_CAPACITY 32

/* Binary search based on asynchronous memory access chaining (AMAC)
 *
 */

typedef uint64_t value_type;
typedef ptrdiff_t index_type;

typedef enum stage_type { prefetch_stage, access_stage } stage_type;

typedef struct state_type {
  stage_type stage;
  int16_t prev;
  int16_t next;
  index_type low;
  index_type size;
  index_type probe;
  value_type value;
  // index_type *result;
} state_type;

typedef struct ring_buffer_type {
  index_type size;
  // index_type first;
  // index_type last;
  state_type data[RING_BUFFER_CAPACITY];
} ring_buffer_type;

static inline index_type ring_buffer_size(ring_buffer_type const *ring_buffer) {
  return ring_buffer->size;
}
static inline index_type ring_buffer_current_index(ring_buffer_type const *ring_buffer,
                                                   state_type const *state) {
  return state - ring_buffer->data;
}

static void init_ring_buffer(ring_buffer_type *ring_buffer,
                             index_type const block_size,
                             index_type const haystack_size,
                             value_type const needles[], index_type indices[]) {
  LS_CHECK(block_size <= RING_BUFFER_CAPACITY, "block size too big");
  ring_buffer->size = block_size;
  // ring_buffer->first = 0;
  // ring_buffer->last = block_size < RING_BUFFER_CAPACITY ? block_size : 0;

  for (index_type i = 0; i < block_size; ++i) {
    ring_buffer->data[i] = (state_type){.stage = prefetch_stage,
                                        .prev = (int16_t)(i - 1),
                                        .next = (int16_t)(i + 1),
                                        .low = 0,
                                        .size = haystack_size,
                                        .probe = -1,
                                        .value = needles[i]};
  }
  ring_buffer->data[0].prev = (int16_t)(ring_buffer->size - 1);
  ring_buffer->data[ring_buffer->size - 1].next = 0;
}

static inline void deinit_ring_buffer(ring_buffer_type *ring_buffer) {
  (void)ring_buffer;
  return;
}

static inline state_type* ring_buffer_first(ring_buffer_type *ring_buffer) {
  return ring_buffer->data;
}

static inline state_type* ring_buffer_next(ring_buffer_type *ring_buffer,
                                           state_type const *state) {
  // LS_CHECK(ring_buffer->size > 0, "ring buffer is empty");
  return ring_buffer->data + state->next;
  // index_type const i = ring_buffer->first;
  // --ring_buffer->size;
  // ++ring_buffer->first;
  // if (ring_buffer->first == RING_BUFFER_CAPACITY) {
  //   ring_buffer->first = 0;
  // }
  // // fprintf(stderr, "  returning data[%zi]\n", i);
  // return ring_buffer->data + i;
}

static inline void ring_buffer_delete(ring_buffer_type *ring_buffer,
                                      state_type const *state) {
  // LS_CHECK(ring_buffer->size < RING_BUFFER_CAPACITY, "ring buffer is full");
  // fprintf(stderr, "  setting data[%zi]\n", ring_buffer->last);
  --ring_buffer->size;
  ring_buffer->data[state->prev].next = state->next;
  ring_buffer->data[state->next].prev = state->prev;

  // ring_buffer->data[ring_buffer->last] = *element;
  // ++ring_buffer->size;
  // ++ring_buffer->last;
  // if (ring_buffer->last == RING_BUFFER_CAPACITY) {
  //   ring_buffer->last = 0;
  // }
}

static inline void prefetch(value_type const *data) {
  __builtin_prefetch(data, 0, 1);
  return;
}

// void ls_hs_internal_block_binary_search(index_type const block_size,
//                                         value_type const haystack[],
//                                         index_type const haystack_size,
//                                         value_type const needles[],
//                                         index_type indices[]) {
//   ring_buffer_type ring_buffer;
//   init_ring_buffer(&ring_buffer, block_size, haystack_size, needles, indices);
//
//   // fprintf(stderr, "block_size=%zi\n", block_size);
//   state_type* state = ring_buffer_first(&ring_buffer);
//   while (ring_buffer_size(&ring_buffer) != 0) {
//     // state_type state = ring_buffer_pop(&ring_buffer);
//     // fprintf(stderr,
//     //         "popped {stage=%i, low=%zi, size=%zi, probe=%zi, value=%zu}\n",
//     //         state.stage, state.low, state.size, state.probe, state.value);
//     switch (state->stage) {
//     case prefetch_stage: {
//       // fprintf(stderr, "prefetch_stage: ");
//       if (state->size == 0) {
//         // fprintf(stderr, "state.size == 0\n");
//         indices[ring_buffer_current_index(&ring_buffer, state)] = -1;
//         ring_buffer_delete(&ring_buffer, state);
//       } else {
//         index_type const half = state->size / 2;
//         state->probe = state->low + half;
//         prefetch(haystack + state->probe);
//         // fprintf(stderr, "half=%zi, probe=%zi\n", half, state.probe);
//         state->stage = access_stage;
//         // ring_buffer_push(&ring_buffer, &state);
//         // fprintf(stderr,
//         //         "pushing {stage=%i, low=%zi, size=%zi, probe=%zi,
//         //         value=%zu}\n", state.stage, state.low, state.size,
//         //         state.probe, state.value);
//       }
//       break;
//     }
//     case access_stage: {
//       // fprintf(stderr, "access_stage: ");
//       value_type const current = haystack[state->probe];
//       // fprintf(stderr, "haystack[%zi]=%zu looking for %zu  ", state.probe,
//       //         current, state.value);
//       if (current == state->value) {
//         // fprintf(stderr, "current == state.value\n");
//         indices[ring_buffer_current_index(&ring_buffer, state)] = state->probe;
//         ring_buffer_delete(&ring_buffer, state);
//         // *state.result = state.probe;
//       } else {
//         if (current < state->value) {
//           state->size -= state->probe + 1 - state->low;
//           state->low = state->probe + 1;
//           // fprintf(stderr, "current is smaller: size=%zi, low=%zi\n",
//           // state.size,
//           //         state.low);
//         } else {
//           state->size = state->probe - state->low;
//           // fprintf(stderr, "current is bigger: size=%zi\n", state.size);
//         }
//
//         state->stage = prefetch_stage;
//         // ring_buffer_push(&ring_buffer, &state);
//         // fprintf(stderr,
//         //         "pushing {stage=%i, low=%zi, size=%zi, probe=%zi,
//         //         value=%zu}\n", state.stage, state.low, state.size,
//         //         state.probe, state.value);
//       }
//       break;
//     }
//     } // end switch
//     state = ring_buffer_next(&ring_buffer, state);
//   }
//
//   deinit_ring_buffer(&ring_buffer);
// }
