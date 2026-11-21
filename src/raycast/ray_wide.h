// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_RAY_WIDE_H
#define ALEA_RAY_WIDE_H
#include <math.h>

typedef struct { double hi, lo; } ray_wide_t;

static inline ray_wide_t ray_wide_madd(ray_wide_t a, ray_wide_t b,
                                            ray_wide_t c) {
    double product = a.hi * b.hi;
    double sum = product + c.hi;
    double v = sum - product;
    double error = (product - (sum - v)) + (c.hi - v);
    error += fma(a.hi, b.hi, -product) + (a.lo * b.hi + a.hi * b.lo) + c.lo;
    double hi = sum + error;
    /* General TwoSum also covers cancellation to a tiny leading part. */
    v = hi - sum;
    return (ray_wide_t){hi, (sum - (hi - v)) + (error - v)};
}

static inline ray_wide_t ray_wide_term(double k, ray_wide_t x, ray_wide_t add) {
    return ray_wide_madd((ray_wide_t){k, 0}, x, add);
}

#endif
