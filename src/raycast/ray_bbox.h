// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef RAY_BBOX_H
#define RAY_BBOX_H

#include "raycast.h"
#include "alea_types.h"
#include <math.h>

/**
 * @file ray_bbox.h
 * @brief Shared ray-AABB slab test
 *
 * Shared by hierarchy traversal, surface culling, and lattice entry/exit.
 */

static inline int ray_bbox_clip_axis(double lo, double hi, double origin,
                                    double dir, double inv_dir,
                                    double* enter, double* exit) {
    if (!(lo <= hi)) return 0;
    if (dir == 0.0) return lo <= origin && origin <= hi;

    double near_delta = (dir > 0.0 ? lo : hi) - origin;
    double far_delta = (dir > 0.0 ? hi : lo) - origin;
    /* A nonzero subnormal direction can overflow its reciprocal. */
    double near_t, far_t;
    if (isinf(inv_dir)) {
        near_t = near_delta / dir;
        far_t = far_delta / dir;
    } else {
        near_t = near_delta * inv_dir;
        far_t = far_delta * inv_dir;
    }
    if (near_t > *enter) *enter = near_t;
    if (far_t < *exit) *exit = far_t;
    return *enter <= *exit;
}

/* Intersect the three closed slabs with the caller's closed ray interval. */
static inline int ray_bbox_slab_enter_exit(const alea_ray_t* ray, const alea_bbox_t* bbox,
                                           double t_min, double t_max,
                                           double* out_enter, double* out_exit) {
    if (!(t_min <= t_max) ||
        !ray_bbox_clip_axis(bbox->min_x, bbox->max_x, ray->ox, ray->dx, ray->inv_dx,
                            &t_min, &t_max) ||
        !ray_bbox_clip_axis(bbox->min_y, bbox->max_y, ray->oy, ray->dy, ray->inv_dy,
                            &t_min, &t_max) ||
        !ray_bbox_clip_axis(bbox->min_z, bbox->max_z, ray->oz, ray->dz, ray->inv_dz,
                            &t_min, &t_max))
        return 0;
    if (out_enter) *out_enter = t_min;
    if (out_exit) *out_exit = t_max;
    return 1;
}

static inline int ray_bbox_slab(const alea_ray_t* ray, const alea_bbox_t* bbox,
                                double t_min, double t_max) {
    return ray_bbox_slab_enter_exit(ray, bbox, t_min, t_max, NULL, NULL);
}

#endif /* RAY_BBOX_H */
