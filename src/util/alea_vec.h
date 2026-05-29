// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_vec.h
 * @brief Type-safe dynamic arrays (internal)
 *
 * Macro-based dynamic arrays with explicit error handling.
 *
 * Usage:
 *   alea_node_vec_t nodes = ALEA_VEC_INIT;
 *   if (alea_vec_push(&nodes, node, alea_node_t) != 0) { handle error }
 *   alea_vec_free(&nodes);
 */

#ifndef ALEA_VEC_H
#define ALEA_VEC_H

#include "alea_types.h"
#include <stdlib.h>
#include <string.h>
#include <stdint.h>  /* SIZE_MAX */

/* ============================================================================
 * TYPE GENERATOR
 * ============================================================================ */

/**
 * Generate a type-safe dynamic array type.
 *
 * Usage:
 *   ALEA_VEC_DEFINE(alea_node_vec, alea_node_t);
 *   // Creates: alea_node_vec_t with { .data, .count, .capacity }
 */
#define ALEA_VEC_DEFINE(name, type)      \
    typedef struct {                     \
        type* data;                      \
        size_t count;                    \
        size_t capacity;                 \
    } name##_t

/* ============================================================================
 * INITIALIZATION
 * ============================================================================ */

/** Static initializer for empty vector */
#define ALEA_VEC_INIT { .data = NULL, .count = 0, .capacity = 0 }

/** Initialize vector at runtime */
#define alea_vec_init(vec) do {  \
    (vec)->data = NULL;         \
    (vec)->count = 0;           \
    (vec)->capacity = 0;        \
} while(0)

/* ============================================================================
 * CAPACITY MANAGEMENT
 * ============================================================================ */

/**
 * Ensure vector has capacity for at least `min_cap` elements.
 * Returns: 0 on success, -1 on failure.
 *
 * If vector is empty (capacity=0), allocates exactly min_cap.
 * Otherwise uses 2x growth for push operations.
 *
 * Usage:
 *   if (alea_vec_reserve(&vec, 100, elem_type) != 0) { handle error }
 */
#define alea_vec_reserve(vec, min_cap, elem_type) __extension__({            \
    int _res = 0;                                                            \
    size_t _min = (min_cap);                                                \
    if ((vec)->capacity < _min) {                                           \
        /* Largest element count whose byte size fits in size_t. Reject up  \
         * front so neither the growth loop nor the realloc multiply below  \
         * can overflow. */                                                  \
        const size_t _max_elems = SIZE_MAX / sizeof(elem_type);             \
        if (_min > _max_elems) {                                            \
            _res = -1;                                                       \
            alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,                    \
                "Vector size %zu exceeds allocation limit", _min);          \
        } else {                                                            \
            size_t _new_cap;                                                \
            if ((vec)->capacity == 0) {                                     \
                _new_cap = _min;  /* Exact allocation for reserve hints */  \
            } else {                                                        \
                _new_cap = (vec)->capacity;                                 \
                /* 2x growth, clamped so doubling cannot overflow or stall  \
                 * (a wrapped _new_cap of 0 would loop forever). */         \
                while (_new_cap < _min) {                                   \
                    if (_new_cap > _max_elems / 2) { _new_cap = _min; break; } \
                    _new_cap *= 2;                                          \
                }                                                           \
            }                                                               \
            void* _new_data = realloc((vec)->data,                          \
                                       _new_cap * sizeof(elem_type));       \
            if (!_new_data) {                                              \
                _res = -1;                                                  \
                alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,               \
                    "Failed to grow vector to %zu elements", _new_cap);     \
            } else {                                                        \
                (vec)->data = _new_data;                                    \
                (vec)->capacity = _new_cap;                                 \
            }                                                               \
        }                                                                   \
    }                                                                       \
    _res;                                                                   \
})

/* ============================================================================
 * PUSH OPERATIONS
 * ============================================================================ */

/**
 * Push element to end of vector. Grows if needed.
 * Returns: 0 on success, -1 on failure.
 *
 * Usage:
 *   if (alea_vec_push(&vec, item, elem_type) != 0) { handle error }
 */
#define alea_vec_push(vec, item, elem_type) __extension__({                  \
    int _res = 0;                                                            \
    if ((vec)->count >= (vec)->capacity) {                                  \
        _res = alea_vec_reserve((vec), (vec)->count + 1, elem_type);         \
    }                                                                       \
    if (_res == 0) {                                                         \
        (vec)->data[(vec)->count++] = (item);                               \
    }                                                                       \
    _res;                                                                   \
})

/**
 * Push and get pointer to uninitialized element.
 * Caller must initialize the element.
 * Returns: Pointer to new element, or NULL on allocation failure.
 *
 * Usage:
 *   alea_node_t* node = alea_vec_push_uninit(&vec, alea_node_t);
 *   if (!node) { handle error }
 *   memset(node, 0, sizeof(*node));
 */
#define alea_vec_push_uninit(vec, elem_type) __extension__({                 \
    elem_type* _ptr = NULL;                                                 \
    if ((vec)->count >= (vec)->capacity) {                                  \
        if (alea_vec_reserve((vec), (vec)->count + 1, elem_type) == 0) {     \
            _ptr = &(vec)->data[(vec)->count++];                            \
        }                                                                   \
    } else {                                                                \
        _ptr = &(vec)->data[(vec)->count++];                                \
    }                                                                       \
    _ptr;                                                                   \
})

/* ============================================================================
 * ACCESS OPERATIONS
 * ============================================================================ */

/** Get element at index (no bounds check) */
#define alea_vec_get(vec, idx) ((vec)->data[idx])

/** Get pointer to element at index (no bounds check) */
#define alea_vec_get_ptr(vec, idx) (&(vec)->data[idx])

/** Get element count */
#define alea_vec_count(vec) ((vec)->count)

/** Check if vector is empty */
#define alea_vec_empty(vec) ((vec)->count == 0)

/** Get last element (undefined if empty) */
#define alea_vec_last(vec) ((vec)->data[(vec)->count - 1])

/** Get pointer to last element (undefined if empty) */
#define alea_vec_last_ptr(vec) (&(vec)->data[(vec)->count - 1])

/* ============================================================================
 * MODIFICATION OPERATIONS
 * ============================================================================ */

/** Clear vector (reset count to 0, keep capacity) */
#define alea_vec_clear(vec) ((vec)->count = 0)

/** Pop last element and return it (undefined if empty) */
#define alea_vec_pop(vec) ((vec)->data[--(vec)->count])

/** Pop last element, discarding the value (no "unused value" warning) */
#define alea_vec_pop_discard(vec) ((void)--(vec)->count)

/** Set count directly (use with caution) */
#define alea_vec_set_count(vec, n) ((vec)->count = (n))

/* ============================================================================
 * CLEANUP
 * ============================================================================ */

/** Free vector memory and reset to empty state */
#define alea_vec_free(vec) do {  \
    free((vec)->data);          \
    (vec)->data = NULL;         \
    (vec)->count = 0;           \
    (vec)->capacity = 0;        \
} while(0)

/* ============================================================================
 * ITERATION
 * ============================================================================ */

/**
 * Iterate over vector elements (pointer to each element)
 *
 * Usage:
 *   alea_vec_foreach(&vec, alea_node_t, node) {
 *       printf("node type: %d\n", node->type);
 *   }
 */
#define alea_vec_foreach(vec, elem_type, var)                                \
    for (elem_type* var = (vec)->data;                                      \
         var < (vec)->data + (vec)->count;                                  \
         var++)

/**
 * Iterate with index
 *
 * Usage:
 *   alea_vec_foreach_i(&vec, alea_node_t, i, node) {
 *       printf("node[%zu] type: %d\n", i, node->type);
 *   }
 */
#define alea_vec_foreach_i(vec, elem_type, idx_var, ptr_var)                 \
    for (size_t idx_var = 0, _cont = 1;                                     \
         idx_var < (vec)->count && _cont;                                   \
         idx_var++, _cont = 1)                                              \
        for (elem_type* ptr_var = &(vec)->data[idx_var]; _cont; _cont = 0)

/* ============================================================================
 * COMMON VECTOR TYPES
 *
 * Pre-defined vector types for common element types.
 * Specific CSG types (alea_node_vec_t, etc.) should be defined
 * where the element types are fully defined.
 * ============================================================================ */

/** Vector of size_t */
ALEA_VEC_DEFINE(alea_size_vec, size_t);

/** Vector of int */
ALEA_VEC_DEFINE(alea_int_vec, int);

/** Vector of uint32_t */
ALEA_VEC_DEFINE(alea_uint32_vec, uint32_t);

/** Vector of double */
ALEA_VEC_DEFINE(alea_double_vec, double);

/** Vector of pointers (void*) */
ALEA_VEC_DEFINE(alea_ptr_vec, void*);

#endif /* ALEA_VEC_H */
