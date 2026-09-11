// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file nuclear_internal.h
 * @brief Internal declarations shared between source files
 *
 * Not part of the public API.
 */

#ifndef NUCLEAR_INTERNAL_H
#define NUCLEAR_INTERNAL_H

#include "alea_nucdata.h"
#include "util/alea_log.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

/* Flag the table as corrupt after an out-of-bounds XSS access. The table is a
 * genuinely mutable, heap-owned object — the const on the helper parameters is
 * only an API courtesy — so recording the error through it is well-defined.
 * The helpers still return a safe value so the in-progress decode does not read
 * out of bounds; the loader inspects this flag afterwards and fails closed. */
static inline void xss_mark_corrupt(const alea_nuc_ace_table_t* t) {
    ((alea_nuc_ace_table_t*)t)->decode_error = true;
}

static inline void xss_mark_allocation_error(const alea_nuc_ace_table_t* t) {
    ((alea_nuc_ace_table_t*)t)->allocation_error = true;
}

/* Helper: get XSS value at 1-based index, with bounds checking */
static inline double xss(const alea_nuc_ace_table_t* t, int idx) {
    if (idx < 1 || idx > t->xss_length) {
        ALEA_LOG_ERROR("XSS out of bounds: index %d, length %d", idx, t->xss_length);
        xss_mark_corrupt(t);
        return 0.0;
    }
    int i = idx - 1;
    return t->xss[i];
}

/* Helper: get XSS value as integer at 1-based index, with bounds checking */
static inline int xss_int(const alea_nuc_ace_table_t* t, int idx) {
    if (idx < 1 || idx > t->xss_length) {
        ALEA_LOG_ERROR("XSS out of bounds: index %d, length %d", idx, t->xss_length);
        xss_mark_corrupt(t);
        return 0;
    }
    int i = idx - 1;
    double value = t->xss[i];
    if (!isfinite(value) || value < INT_MIN || value > INT_MAX || value != trunc(value)) {
        ALEA_LOG_ERROR("invalid integer in XSS at index %d: %.17g", idx, value);
        xss_mark_corrupt(t);
        return 0;
    }
    return (int)value;
}

/* Check a 1-based XSS range without overflowing signed arithmetic. */
static inline bool xss_range_valid(const alea_nuc_ace_table_t* t,
                                   int start, int n) {
    if (start < 1 || n < 0 || n > t->xss_length ||
        start - 1 > t->xss_length - n) {
        ALEA_LOG_ERROR("XSS range out of bounds: start %d, n %d, length %d",
                       start, n, t->xss_length);
        xss_mark_corrupt(t);
        return false;
    }
    return true;
}

/* Resolve a positive 1-based locator relative to a 1-based XSS block. */
static inline int xss_relative_loc(const alea_nuc_ace_table_t* t,
                                   int base, int locator) {
    if (base < 1 || locator < 1 ||
        locator - 1 > t->xss_length - base) {
        ALEA_LOG_ERROR("invalid XSS locator: base %d, locator %d, length %d",
                       base, locator, t->xss_length);
        xss_mark_corrupt(t);
        return 0;
    }
    return base + locator - 1;
}

/* Copy n doubles from XSS starting at 1-based index, with bounds checking */
static inline double* xss_copy(const alea_nuc_ace_table_t* t, int start, int n) {
    if (n <= 0) return NULL;
    if (!xss_range_valid(t, start, n)) return NULL;
    int i = start - 1;
    double* arr = malloc((size_t)n * sizeof(double));
    if (arr)
        memcpy(arr, &t->xss[i], (size_t)n * sizeof(double));
    else
        xss_mark_allocation_error(t);
    return arr;
}

/* angular.c — decode */
alea_nuc_angular_dist_t* alea_nuc_decode_angular(const alea_nuc_ace_table_t* t, int n_reactions);
void alea_nuc_decode_all_angular(alea_nuc_nuclide_t* nuc);

/* energy_dist.c — decode */
alea_nuc_energy_dist_t* alea_nuc_decode_energy_dist(const alea_nuc_ace_table_t* t, int ldlw_loc);
void alea_nuc_decode_all_energy(alea_nuc_nuclide_t* nuc);

/* lookup.c — fast path for immutable grids validated during decoding */
int alea_nuc_energy_lookup_trusted(const double* energy, int n, double E,
                                   double* frac);

#endif /* NUCLEAR_INTERNAL_H */
