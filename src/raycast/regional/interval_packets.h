// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* A packet contains four bounding volumes in nine projection directions.
 * Stored float bounds have been rounded outward. Distances and all query
 * comparisons use double precision on both paths. */
#if defined(__AVX__)
#include <immintrin.h>
#endif

static unsigned packet_point(const interval_index_t* idx,const bound_packet_t* packet, const double p[BOUND_AXES]) {
#if defined(__AVX__)
    __m256d active = _mm256_castsi256_pd(_mm256_set1_epi32(-1));
    size_t pos=packet->first;
    for (unsigned k = 0; k < BOUND_AXES; k++) {
        if(k>=3 && !(packet->mask&(1u<<(k-3))))continue;
        const packet_axis_t* axis=k<3?packet->xyz+k:idx->extra_axes.data+pos++;
        __m256d q = _mm256_set1_pd(p[k]);
        __m256d lo = _mm256_cvtps_pd(_mm_loadu_ps(axis->lo));
        __m256d hi = _mm256_cvtps_pd(_mm_loadu_ps(axis->hi));
        active = _mm256_and_pd(active, _mm256_and_pd(
            _mm256_cmp_pd(lo, q, _CMP_LE_OQ), _mm256_cmp_pd(q, hi, _CMP_LE_OQ)));
        if (k == 2 && !_mm256_movemask_pd(active)) return 0;
    }
    return (unsigned)_mm256_movemask_pd(active);
#else
    unsigned active = (1u << BOUND_LANES) - 1;
    size_t pos=packet->first;
    for (unsigned k = 0; k < BOUND_AXES && active; k++) {
        if(k>=3 && !(packet->mask&(1u<<(k-3))))continue;
        const packet_axis_t* axis=k<3?packet->xyz+k:idx->extra_axes.data+pos++;
        unsigned matching = 0;
        for (unsigned i = 0; i < BOUND_LANES; i++)
            if (p[k] >= axis->lo[i] && p[k] <= axis->hi[i]) matching |= 1u << i;
        active &= matching;
    }
    return active;
#endif
}

static unsigned packet_ray(const interval_level_t* level, const bound_packet_t* packet,
                           double lo, double hi, double enter[4], double leave[4]) {
    const interval_index_t* idx=level->query->index;
#if defined(__AVX__)
    __m256d a = _mm256_set1_pd(lo), b = _mm256_set1_pd(hi);
    size_t pos=packet->first;
    for (unsigned k = 0; k < BOUND_AXES; k++) {
        if(k>=3 && !(packet->mask&(1u<<(k-3))))continue;
        const packet_axis_t* axis=k<3?packet->xyz+k:idx->extra_axes.data+pos++;
        double direction = level->ds[k];
        const float* lower = direction < 0 ? axis->hi : axis->lo;
        const float* upper = direction < 0 ? axis->lo : axis->hi;
        __m256d pad = _mm256_set1_pd(direction < 0 ? -level->err[k] : level->err[k]);
        __m256d x = _mm256_sub_pd(_mm256_cvtps_pd(_mm_loadu_ps(lower)), pad);
        __m256d y = _mm256_add_pd(_mm256_cvtps_pd(_mm_loadu_ps(upper)), pad);
        __m256d origin = _mm256_set1_pd(level->os[k]);
        if (direction == 0) {
            __m256d outside = _mm256_or_pd(_mm256_cmp_pd(origin, x, _CMP_LT_OQ),
                                            _mm256_cmp_pd(origin, y, _CMP_GT_OQ));
            a = _mm256_blendv_pd(a, _mm256_set1_pd(INFINITY), outside);
        } else {
            __m256d inverse = _mm256_set1_pd(level->iv[k]);
            a = _mm256_max_pd(_mm256_mul_pd(_mm256_sub_pd(x, origin), inverse), a);
            b = _mm256_min_pd(_mm256_mul_pd(_mm256_sub_pd(y, origin), inverse), b);
        }
        if (k == 2 && !_mm256_movemask_pd(_mm256_cmp_pd(a, b, _CMP_LE_OQ))) return 0;
    }
    _mm256_storeu_pd(enter, a);
    _mm256_storeu_pd(leave, b);
    return (unsigned)_mm256_movemask_pd(_mm256_cmp_pd(a, b, _CMP_LE_OQ));
#else
    unsigned active = 0;
    for (unsigned i = 0; i < BOUND_LANES; i++) {
        size_t pos=packet->first;
        double a = lo, b = hi;
        for (unsigned k = 0; k < BOUND_AXES && a <= b; k++) {
            if(k>=3 && !(packet->mask&(1u<<(k-3))))continue;
            const packet_axis_t* axis=k<3?packet->xyz+k:idx->extra_axes.data+pos++;
            double origin = level->os[k], direction = level->ds[k];
            double low = axis->lo[i] - level->err[k], high = axis->hi[i] + level->err[k];
            if (direction == 0) {
                if (origin < low || origin > high) a = INFINITY;
            } else {
                double x = ((direction > 0 ? low : high) - origin) * level->iv[k];
                double y = ((direction > 0 ? high : low) - origin) * level->iv[k];
                if (x > a) a = x;
                if (y < b) b = y;
            }
        }
        enter[i] = a;
        leave[i] = b;
        if (a <= b) active |= 1u << i;
    }
    return active;
#endif
}
