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
 * Branchless slab test using IEEE 754 inf from precomputed inv_dir.
 * Used by BVH traversal, surface culling, and lattice entry/exit.
 */

/**
 * @brief Test if ray intersects an AABB within [t_min, t_max].
 * @return 1 if ray hits, 0 if miss
 */
static inline int ray_bbox_slab(const alea_ray_t* ray, const alea_bbox_t* bbox,
                                double t_min, double t_max) {
    double t1x = (bbox->min_x - ray->ox) * ray->inv_dx;
    double t2x = (bbox->max_x - ray->ox) * ray->inv_dx;
    double tmin = fmin(t1x, t2x);
    double tmax = fmax(t1x, t2x);

    double t1y = (bbox->min_y - ray->oy) * ray->inv_dy;
    double t2y = (bbox->max_y - ray->oy) * ray->inv_dy;
    tmin = fmax(tmin, fmin(t1y, t2y));
    tmax = fmin(tmax, fmax(t1y, t2y));

    double t1z = (bbox->min_z - ray->oz) * ray->inv_dz;
    double t2z = (bbox->max_z - ray->oz) * ray->inv_dz;
    tmin = fmax(tmin, fmin(t1z, t2z));
    tmax = fmin(tmax, fmax(t1z, t2z));

    return (tmin <= t_max) && (tmax >= t_min);
}

/**
 * @brief Test ray-AABB and return entry/exit t values.
 * @return 1 if ray hits (out_enter/out_exit set), 0 if miss
 */
static inline int ray_bbox_slab_enter_exit(const alea_ray_t* ray, const alea_bbox_t* bbox,
                                           double t_min, double t_max,
                                           double* out_enter, double* out_exit) {
    double t1x = (bbox->min_x - ray->ox) * ray->inv_dx;
    double t2x = (bbox->max_x - ray->ox) * ray->inv_dx;
    double tmin = fmin(t1x, t2x);
    double tmax = fmax(t1x, t2x);

    double t1y = (bbox->min_y - ray->oy) * ray->inv_dy;
    double t2y = (bbox->max_y - ray->oy) * ray->inv_dy;
    tmin = fmax(tmin, fmin(t1y, t2y));
    tmax = fmin(tmax, fmax(t1y, t2y));

    double t1z = (bbox->min_z - ray->oz) * ray->inv_dz;
    double t2z = (bbox->max_z - ray->oz) * ray->inv_dz;
    tmin = fmax(tmin, fmin(t1z, t2z));
    tmax = fmin(tmax, fmax(t1z, t2z));

    if (tmin > t_max || tmax < t_min) return 0;
    if (tmin < t_min) tmin = t_min;
    if (tmax > t_max) tmax = t_max;
    *out_enter = tmin;
    *out_exit = tmax;
    return 1;
}

#endif /* RAY_BBOX_H */
