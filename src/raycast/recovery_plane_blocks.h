// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_RECOVERY_PLANE_BLOCKS_H
#define ALEA_RECOVERY_PLANE_BLOCKS_H

#include "ray_intersect.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

/* Private preparation scratch. The key belongs to this sort, including when
 * separate models prepare their indices concurrently. */
typedef struct {
    alea_plane_data_t plane;
    uint32_t ordinal;
    double key;
} recovery_plane_sort_t;

typedef struct {
    double center[4], radius[4], scale;
    uint32_t begin, end;
} recovery_plane_block_t;

enum { RECOVERY_PLANE_BLOCK_SIZE = 24 };

static inline double recovery_plane_component(const alea_plane_data_t* p, unsigned axis) {
    switch (axis) {
        case 0: return p->a;
        case 1: return p->b;
        case 2: return p->c;
        default: return p->d;
    }
}

static inline int recovery_plane_key_compare(const void* a, const void* b) {
    const recovery_plane_sort_t *x = a, *y = b;
    return (x->key > y->key) - (x->key < y->key);
}

/* Partition by the widest coefficient range. Scale d by the offsets in
 * this universe so that its units do not swamp differences in normals. */
static inline void recovery_plane_partition(recovery_plane_sort_t* planes,
        uint32_t begin, uint32_t end, double offset_scale,
        recovery_plane_block_t* blocks, uint32_t* block_count, uint32_t* ordinals) {
    double lo[4], hi[4];
    for (unsigned k = 0; k < 4; k++)
        lo[k] = hi[k] = recovery_plane_component(&planes[begin].plane, k);
    for (uint32_t i = begin + 1; i < end; i++) {
        for (unsigned k = 0; k < 4; k++) {
            double v = recovery_plane_component(&planes[i].plane, k);
            lo[k] = fmin(lo[k], v);
            hi[k] = fmax(hi[k], v);
        }
    }
    if (end - begin > RECOVERY_PLANE_BLOCK_SIZE) {
        double widest = -1;
        unsigned axis = 0;
        for (unsigned k = 0; k < 4; k++) {
            double width = hi[k] - lo[k];
            if (k == 3) width /= offset_scale;
            if (width > widest) { widest = width; axis = k; }
        }
        for (uint32_t i = begin; i < end; i++)
            planes[i].key = recovery_plane_component(&planes[i].plane, axis);
        qsort(planes + begin, end - begin, sizeof(*planes), recovery_plane_key_compare);
        uint32_t mid = begin + (end - begin) / 2;
        recovery_plane_partition(planes, begin, mid, offset_scale, blocks, block_count, ordinals);
        recovery_plane_partition(planes, mid, end, offset_scale, blocks, block_count, ordinals);
        return;
    }
    recovery_plane_block_t* b = &blocks[(*block_count)++];
    b->begin = begin;
    b->end = end;
    b->scale = 1;
    for (unsigned k = 0; k < 4; k++) {
        b->center[k] = 0.5 * lo[k] + 0.5 * hi[k];
        b->radius[k] = fmax(hi[k] - b->center[k], b->center[k] - lo[k]);
        b->scale = fmax(b->scale, fmax(fabs(lo[k]), fabs(hi[k])));
    }
    for (uint32_t i = begin; i < end; i++) ordinals[i] = planes[i].ordinal;
}

/* Each plane is affine along the ray. Bound the variation around the
 * representative plane at BOTH endpoints. Strict agreement in sign excludes
 * every member throughout the closed interval. Keep uncertainty in the mask. */
static inline void recovery_plane_prune_blocks(const recovery_plane_block_t* blocks,
        uint32_t count, const uint32_t* ordinals, const alea_ray_t* ray,
        double lo, double hi, uint64_t* mask) {
    if (!count || !(lo >= 0 && hi >= lo) || !isfinite(hi)) return;
    double scale = 1 + fabs(ray->ox) + fabs(ray->oy) + fabs(ray->oz)
                     + fmax(1.0, hi) * (1 + fabs(ray->dx) + fabs(ray->dy) + fabs(ray->dz));
    if (!isfinite(scale)) return;
    double x0 = ray->ox + lo * ray->dx, x1 = ray->ox + hi * ray->dx;
    double y0 = ray->oy + lo * ray->dy, y1 = ray->oy + hi * ray->dy;
    double z0 = ray->oz + lo * ray->dz, z1 = ray->oz + hi * ray->dz;
    double x = fmax(fabs(x0), fabs(x1)), y = fmax(fabs(y0), fabs(y1));
    double z = fmax(fabs(z0), fabs(z1));
    for (uint32_t i = 0; i < count; i++) {
        const recovery_plane_block_t* b = &blocks[i];
        const double *c = b->center, *r = b->radius;
        double magnitude = b->scale * scale;
        if (!(magnitude < DBL_MAX / 64)) continue;
        /* Allow for rounding the block and its evaluations, forming the
         * endpoints, and the scalar kernel's numerator, denominator and
         * division. A lower bound of one on both scales covers underflow. */
        double error = r[0]*x + r[1]*y + r[2]*z + r[3] + (128*DBL_EPSILON)*magnitude;
        double f0 = c[0]*x0 + c[1]*y0 + c[2]*z0 + c[3];
        double f1 = c[0]*x1 + c[1]*y1 + c[2]*z1 + c[3];
        if ((f0 > error && f1 > error) || (f0 < -error && f1 < -error)) {
            for (uint32_t j = b->begin; j < b->end; j++) {
                uint32_t card = ordinals[j];
                mask[card / 64] &= ~(UINT64_C(1) << (card % 64));
            }
        }
    }
}
#endif
