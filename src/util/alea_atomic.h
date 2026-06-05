// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_atomic.h
 * @brief Portable C11 atomics shim (GCC/Clang native, MSVC via intrinsics)
 *
 * On compilers with a working `<stdatomic.h>` (GCC, Clang, MSVC 17.8+ in C11
 * mode) this header is a thin pass-through: it includes `<stdatomic.h>` and
 * everything behaves exactly as before.
 *
 * On MSVC without native C11 atomics it maps the subset of atomics used by
 * this project onto the `_Interlocked*` intrinsics from `<intrin.h>`:
 *
 *   - `atomic_int`              -> volatile long
 *   - `atomic_uint`             -> volatile unsigned long
 *   - `atomic_uint_fast64_t`    -> volatile unsigned __int64
 *   - `_Atomic T` (size_t/int)  -> volatile T
 *   - `atomic_flag`             -> { volatile long }
 *
 * Note on memory ordering: the `_Interlocked*` intrinsics are full barriers on
 * x86/x64, which is *stronger* than the `memory_order_relaxed` requested at
 * some call sites. That is always correct (never weaker than required), at a
 * marginal performance cost on the relaxed counters.
 *
 * Only the operations actually used in this codebase are provided:
 *   load, store, init, fetch_add (+_explicit), fetch_or, fetch_and,
 *   compare_exchange_strong, flag_test_and_set, flag_clear.
 * Add more here if new call sites need them.
 */

#ifndef ALEA_ATOMIC_H
#define ALEA_ATOMIC_H

#if defined(_MSC_VER) && !defined(__clang__) && \
    (defined(__STDC_NO_ATOMICS__) || !defined(__STDC_VERSION__) || __STDC_VERSION__ < 201112L)

/* ----------------------------------------------------------------------------
 * MSVC fallback: no usable <stdatomic.h>. Map onto Interlocked intrinsics.
 * ------------------------------------------------------------------------- */

#include <intrin.h>
#include <stddef.h>

/* `_Atomic T` declarations (e.g. `_Atomic size_t`) collapse to `volatile T`.
 * `_Atomic` is not a keyword here (atomics unsupported), so defining it as a
 * macro is safe on this path only. */
#define _Atomic volatile

typedef volatile long             atomic_int;
typedef volatile unsigned long    atomic_uint;
typedef volatile unsigned __int64 atomic_uint_fast64_t;

typedef enum {
    memory_order_relaxed,
    memory_order_consume,
    memory_order_acquire,
    memory_order_release,
    memory_order_acq_rel,
    memory_order_seq_cst
} memory_order;

/* All atomic objects in this project are 32-bit (int/uint) or 64-bit
 * (size_t/uint_fast64_t). Each macro dispatches on sizeof at compile time;
 * both branches are valid C, the constant condition selects the right width. */

#define atomic_init(p, v)   ((void)(*(p) = (v)))

#define atomic_load(p) \
    (sizeof(*(p)) == 8 \
        ? (unsigned long long)_InterlockedOr64((volatile __int64 *)(p), 0) \
        : (unsigned long long)_InterlockedOr((volatile long *)(p), 0))

#define atomic_store(p, v) \
    ((void)(sizeof(*(p)) == 8 \
        ? _InterlockedExchange64((volatile __int64 *)(p), (__int64)(v)) \
        : (__int64)_InterlockedExchange((volatile long *)(p), (long)(v))))

#define atomic_fetch_add(p, v) \
    (sizeof(*(p)) == 8 \
        ? (unsigned long long)_InterlockedExchangeAdd64((volatile __int64 *)(p), (__int64)(v)) \
        : (unsigned long long)_InterlockedExchangeAdd((volatile long *)(p), (long)(v)))

#define atomic_fetch_or(p, v) \
    (sizeof(*(p)) == 8 \
        ? (unsigned long long)_InterlockedOr64((volatile __int64 *)(p), (__int64)(v)) \
        : (unsigned long long)_InterlockedOr((volatile long *)(p), (long)(v)))

#define atomic_fetch_and(p, v) \
    (sizeof(*(p)) == 8 \
        ? (unsigned long long)_InterlockedAnd64((volatile __int64 *)(p), (__int64)(v)) \
        : (unsigned long long)_InterlockedAnd((volatile long *)(p), (long)(v)))

/* memory-order-explicit variants: ordering arg ignored (intrinsics are seq_cst) */
#define atomic_fetch_add_explicit(p, v, mo) atomic_fetch_add((p), (v))
#define atomic_load_explicit(p, mo)         atomic_load((p))
#define atomic_store_explicit(p, v, mo)     atomic_store((p), (v))

static inline int alea__atomic_cas(volatile void *p, void *expected,
                                   unsigned long long desired, size_t width) {
    if (width == 8) {
        __int64 exp = *(__int64 *)expected;
        __int64 old = _InterlockedCompareExchange64((volatile __int64 *)p,
                                                    (__int64)desired, exp);
        if (old == exp) return 1;
        *(__int64 *)expected = old;
        return 0;
    } else {
        long exp = *(long *)expected;
        long old = _InterlockedCompareExchange((volatile long *)p,
                                               (long)desired, exp);
        if (old == exp) return 1;
        *(long *)expected = old;
        return 0;
    }
}

#define atomic_compare_exchange_strong(p, expected, desired) \
    alea__atomic_cas((volatile void *)(p), (void *)(expected), \
                     (unsigned long long)(desired), sizeof(*(p)))

typedef struct { volatile long _v; } atomic_flag;
#define ATOMIC_FLAG_INIT { 0 }
#define atomic_flag_test_and_set(f) (_InterlockedExchange(&(f)->_v, 1) != 0)
#define atomic_flag_clear(f)        ((void)_InterlockedExchange(&(f)->_v, 0))

#else /* GCC / Clang / MSVC with native C11 atomics */

#include <stdatomic.h>

#endif

#endif /* ALEA_ATOMIC_H */
