// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef TINYPAR_ATOMIC_H
#define TINYPAR_ATOMIC_H

#include <stddef.h>

#if defined(_MSC_VER) && !defined(__clang__) && \
    (defined(__STDC_NO_ATOMICS__) || !defined(__STDC_VERSION__) || __STDC_VERSION__ < 201112L)

#include <intrin.h>

typedef volatile size_t tinypar_atomic_size_t;
typedef volatile long tinypar_atomic_int_t;

static inline void tinypar_atomic_size_init(tinypar_atomic_size_t* value,
                                            size_t initial) {
    *value = initial;
}

static inline size_t tinypar_atomic_size_load(
        const tinypar_atomic_size_t* value) {
    if (sizeof(size_t) == 8) {
        return (size_t)_InterlockedCompareExchange64(
            (volatile __int64*)value, 0, 0);
    }
    return (size_t)_InterlockedCompareExchange((volatile long*)value, 0, 0);
}

static inline void tinypar_atomic_size_store(tinypar_atomic_size_t* value,
                                             size_t desired) {
    if (sizeof(size_t) == 8) {
        (void)_InterlockedExchange64(
            (volatile __int64*)value, (__int64)desired);
    } else {
        (void)_InterlockedExchange((volatile long*)value, (long)desired);
    }
}

static inline int tinypar_atomic_size_compare_exchange(
        tinypar_atomic_size_t* value, size_t* expected, size_t desired) {
    size_t observed;
    if (sizeof(size_t) == 8) {
        observed = (size_t)_InterlockedCompareExchange64(
            (volatile __int64*)value, (__int64)desired, (__int64)*expected);
    } else {
        observed = (size_t)_InterlockedCompareExchange(
            (volatile long*)value, (long)desired, (long)*expected);
    }
    if (observed == *expected) return 1;
    *expected = observed;
    return 0;
}

static inline void tinypar_atomic_int_init(tinypar_atomic_int_t* value,
                                           int initial) {
    *value = (long)initial;
}

static inline int tinypar_atomic_int_load(const tinypar_atomic_int_t* value) {
    return (int)_InterlockedCompareExchange((volatile long*)value, 0, 0);
}

static inline void tinypar_atomic_int_store(tinypar_atomic_int_t* value,
                                            int desired) {
    (void)_InterlockedExchange(value, (long)desired);
}

static inline int tinypar_atomic_int_compare_exchange(tinypar_atomic_int_t* value,
                                                       int* expected,
                                                       int desired) {
    long observed = _InterlockedCompareExchange(value, (long)desired,
                                                (long)*expected);
    if (observed == (long)*expected) return 1;
    *expected = (int)observed;
    return 0;
}

#else

#include <stdatomic.h>

typedef _Atomic size_t tinypar_atomic_size_t;
typedef _Atomic int tinypar_atomic_int_t;

static inline void tinypar_atomic_size_init(tinypar_atomic_size_t* value,
                                            size_t initial) {
    atomic_init(value, initial);
}

static inline size_t tinypar_atomic_size_load(
        const tinypar_atomic_size_t* value) {
    return atomic_load_explicit(value, memory_order_relaxed);
}

static inline void tinypar_atomic_size_store(tinypar_atomic_size_t* value,
                                             size_t desired) {
    atomic_store_explicit(value, desired, memory_order_relaxed);
}

static inline int tinypar_atomic_size_compare_exchange(
        tinypar_atomic_size_t* value, size_t* expected, size_t desired) {
    return atomic_compare_exchange_weak_explicit(
        value, expected, desired, memory_order_relaxed, memory_order_relaxed);
}

static inline void tinypar_atomic_int_init(tinypar_atomic_int_t* value,
                                           int initial) {
    atomic_init(value, initial);
}

static inline int tinypar_atomic_int_load(const tinypar_atomic_int_t* value) {
    return atomic_load_explicit(value, memory_order_acquire);
}

static inline void tinypar_atomic_int_store(tinypar_atomic_int_t* value,
                                            int desired) {
    atomic_store_explicit(value, desired, memory_order_release);
}

static inline int tinypar_atomic_int_compare_exchange(tinypar_atomic_int_t* value,
                                                       int* expected,
                                                       int desired) {
    return atomic_compare_exchange_strong_explicit(
        value, expected, desired, memory_order_acq_rel, memory_order_acquire);
}

#endif

#endif /* TINYPAR_ATOMIC_H */
