// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file raycast_api.c
 * @brief Public API wrappers for raycast module
 *
 * These functions are thin wrappers that provide the alea_* public API
 * for raycast functionality. They are only available when linking with
 * the raycast module.
 *
 * Volume estimation functions (alea_estimate_cell_volumes, alea_remove_cells_by_volume,
 * alea_estimate_instance_volumes) are strong-symbol overrides of the weak stubs in
 * alea_module_stubs.c. The public alea_* wrappers live in alea_public_api.c.
 */

#define _USE_MATH_DEFINES
#include "alea.h"
#include "raycast.h"
#include "core/alea_system.h"
#include "core/alea_spatial.h"
#include "primitives/bbox.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util/math.h"




alea_raycast_result_t* alea_raycast_result_create(void) {
    alea_raycast_result_t* result = calloc(1, sizeof(alea_raycast_result_t));
    if (result) {
        alea_raycast_result_init(result);
    }
    return result;
}

size_t alea_raycast_segment_count(const alea_raycast_result_t* result) {
    return result ? result->segment_count : 0;
}

int alea_raycast_segment_get(const alea_raycast_result_t* result, size_t index,
                                 double* t_enter, double* t_exit,
                                 int* cell_id, int* material_id, double* density) {
    if (!result || index >= result->segment_count) return -1;
    alea_ray_segment_t* seg = &result->segments[index];
    if (t_enter) *t_enter = seg->t_enter;
    if (t_exit) *t_exit = seg->t_exit;
    if (cell_id) *cell_id = seg->cell_id;
    if (material_id) *material_id = seg->material_id;
    if (density) *density = seg->density;
    return 0;
}



/* Simple LCG for deterministic random numbers (no global state mutation). */
static uint32_t lcg_next(uint32_t* state) {
    *state = *state * 1664525u + 1013904223u;
    return *state;
}

static double lcg_double(uint32_t* state) {
    return (double)lcg_next(state) / 4294967296.0;
}

int alea_estimate_cell_volumes(const alea_system_t* sys,
                              double ox, double oy, double oz,
                              double radius, int n_rays,
                              double* volumes, double* rel_errors) {
    if (!sys || radius <= 0.0 || n_rays <= 0 || !volumes) return -1;

    /* Ensure all raycast caches are built (BVH, spatial index, adjacency) */
    alea_raycast_ensure_caches(sys);

    size_t n_cells = alea_vec_count(&sys->cells);
    if (n_cells == 0) return 0;

    double cx = ox, cy = oy, cz = oz;
    double R = radius;

    /* Zero output */
    memset(volumes, 0, n_cells * sizeof(double));
    if (rel_errors) memset(rel_errors, 0, n_cells * sizeof(double));

    /* Allocate sum-of-squares for error estimation */
    double* sum_l2 = NULL;
    if (rel_errors) {
        sum_l2 = calloc(n_cells, sizeof(double));
        if (!sum_l2) return -1;
    }

    /* Per-ray track length accumulator for a single ray.
     * Needed because one ray can cross a cell multiple times, and we need
     * the total per-ray contribution L_i for that cell to compute L_i^2. */
    double* ray_l = NULL;
    if (rel_errors) {
        ray_l = calloc(n_cells, sizeof(double));
        if (!ray_l) { free(sum_l2); return -1; }
    }

    /* Trace mu-random rays (Cauchy-Crofton):
     *   1. Pick uniform random direction u
     *   2. Pick uniform random point in disk of radius R perpendicular to u
     *   3. Ray origin = disk_point - 2R * u  (starts well outside sphere)
     *   4. Trace with t_max = 4R (covers full sphere diameter + margin)
     *
     * Volume formula: V_cell = pi * R^2 * (sum_track_length / n_rays)
     * Error:  sigma_V = pi * R^2 / sqrt(N) * sqrt(E[L^2] - E[L]^2)
     */
    alea_raycast_result_t result;
    alea_raycast_result_init(&result);
    uint32_t rng = 42;  /* deterministic seed */

    for (int ray = 0; ray < n_rays; ray++) {
        /* Random unit direction */
        double phi = 2.0 * M_PI * lcg_double(&rng);
        double cos_theta = 2.0 * lcg_double(&rng) - 1.0;
        double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
        double ux = sin_theta * cos(phi);
        double uy = sin_theta * sin(phi);
        double uz = cos_theta;

        /* Build orthonormal frame (v1, v2) perpendicular to u */
        double tx, ty, tz;
        if (fabs(ux) < 0.9) {
            tx = 1; ty = 0; tz = 0;
        } else {
            tx = 0; ty = 1; tz = 0;
        }
        /* v1 = normalize(cross(u, t)) */
        double v1x = uy * tz - uz * ty;
        double v1y = uz * tx - ux * tz;
        double v1z = ux * ty - uy * tx;
        double v1_len = sqrt(v1x * v1x + v1y * v1y + v1z * v1z);
        v1x /= v1_len; v1y /= v1_len; v1z /= v1_len;
        /* v2 = cross(u, v1) */
        double v2x = uy * v1z - uz * v1y;
        double v2y = uz * v1x - ux * v1z;
        double v2z = ux * v1y - uy * v1x;

        /* Random point in disk of radius R */
        double r_disk = R * sqrt(lcg_double(&rng));
        double angle = 2.0 * M_PI * lcg_double(&rng);
        double d1 = r_disk * cos(angle);
        double d2 = r_disk * sin(angle);

        /* Ray origin: center + disk offset - 2R along direction */
        double rox = cx + v1x * d1 + v2x * d2 - ux * 2.0 * R;
        double roy = cy + v1y * d1 + v2y * d2 - uy * 2.0 * R;
        double roz = cz + v1z * d1 + v2z * d2 - uz * 2.0 * R;

        if (ray_l) memset(ray_l, 0, n_cells * sizeof(double));

        alea_raycast_result_clear(&result);
        int rc = alea_raycast(sys, rox, roy, roz, ux, uy, uz, 4.0 * R, &result);
        if (rc != 0) {
            /* Ray contributed L=0 to every cell — still counts for variance */
            continue;
        }

        /* Accumulate track lengths per cell */
        for (size_t s = 0; s < result.segment_count; s++) {
            int seg_cell_id = result.segments[s].cell_id;
            if (seg_cell_id < 0) continue;
            double len = result.segments[s].t_exit - result.segments[s].t_enter;
            if (len <= 0) continue;

            int ci = alea_find_cell_by_id(sys, seg_cell_id);
            if (ci >= 0 && (size_t)ci < n_cells) {
                volumes[ci] += len;
                if (ray_l) ray_l[ci] += len;
            }
        }

        /* Accumulate L^2 for this ray */
        if (sum_l2) {
            for (size_t ci = 0; ci < n_cells; ci++) {
                sum_l2[ci] += ray_l[ci] * ray_l[ci];
            }
        }
    }
    alea_raycast_result_free(&result);

    /* Convert accumulated track lengths to volumes and relative errors:
     * V_cell = pi * R^2 * (sum_L / N)
     * sigma_V = pi * R^2 / sqrt(N) * sqrt(E[L^2] - E[L]^2)
     * rel_error = sigma_V / V = sqrt(E[L^2] - E[L]^2) / (E[L] * sqrt(N))  */
    double scale = M_PI * R * R / (double)n_rays;
    for (size_t i = 0; i < n_cells; i++) {
        double mean_l = volumes[i] / (double)n_rays;
        if (rel_errors) {
            if (mean_l > 0.0) {
                double mean_l2 = sum_l2[i] / (double)n_rays;
                double var_l = mean_l2 - mean_l * mean_l;
                if (var_l < 0.0) var_l = 0.0;
                rel_errors[i] = sqrt(var_l) / (mean_l * sqrt((double)n_rays));
            } else {
                rel_errors[i] = -1.0;
            }
        }
        volumes[i] *= scale;
    }

    free(sum_l2);
    free(ray_l);
    return 0;
}

int alea_remove_cells_by_volume(alea_system_t* sys,
                               const double* volumes,
                               double threshold) {
    if (!sys || !volumes) return -1;

    size_t n_cells = alea_vec_count(&sys->cells);
    if (n_cells == 0) return 0;

    size_t write = 0;
    int removed = 0;
    for (size_t i = 0; i < n_cells; i++) {
        if (volumes[i] <= threshold) {
            free(sys->cells.data[i].surface_indices);
            free(sys->cells.data[i].neighbors);
            free(sys->cells.data[i].lat_fill);
            removed++;
        } else {
            if (write != i)
                sys->cells.data[write] = sys->cells.data[i];
            write++;
        }
    }
    sys->cells.count = write;

    return removed;
}

/* ============================================================================
 * PER-INSTANCE VOLUME ESTIMATION
 *
 * Like alea_estimate_cell_volumes but resolves each ray segment to a specific
 * cell instance via the spatial index.  This distinguishes the same cell
 * appearing in multiple fill contexts.
 * ============================================================================ */

int alea_estimate_instance_volumes(const alea_system_t* sys,
                                   int n_rays,
                                   double* volumes, double* rel_errors) {
    if (!sys || n_rays <= 0 || !volumes) return -1;
    if (!sys->spatial_index || !sys->spatial_index->built) return -1;

    size_t n_instances = sys->spatial_index->instance_count;
    if (n_instances == 0) return 0;

    /* Ensure raycast caches */
    alea_raycast_ensure_caches(sys);

    /* Compute bounding sphere from spatial index global bounds */
    const alea_bbox_t* bounds = &sys->spatial_index->bounds;
    double cx = (bounds->min_x + bounds->max_x) * 0.5;
    double cy = (bounds->min_y + bounds->max_y) * 0.5;
    double cz = (bounds->min_z + bounds->max_z) * 0.5;
    double dx = bounds->max_x - bounds->min_x;
    double dy = bounds->max_y - bounds->min_y;
    double dz = bounds->max_z - bounds->min_z;
    double R = 0.5 * sqrt(dx * dx + dy * dy + dz * dz) * 1.01;

    if (R <= 0.0) return -1;

    memset(volumes, 0, n_instances * sizeof(double));
    if (rel_errors) memset(rel_errors, 0, n_instances * sizeof(double));

    /* Allocate sum-of-squares and per-ray accumulators */
    double* sum_l2 = NULL;
    double* ray_l = NULL;
    if (rel_errors) {
        sum_l2 = calloc(n_instances, sizeof(double));
        ray_l = calloc(n_instances, sizeof(double));
        if (!sum_l2 || !ray_l) { free(sum_l2); free(ray_l); return -1; }
    }

    alea_raycast_result_t result;
    alea_raycast_result_init(&result);
    uint32_t rng = 42;

    /* Pre-allocate hit buffer for spatial queries */
    size_t max_hits = 64;
    alea_spatial_hit_t* hits = malloc(max_hits * sizeof(alea_spatial_hit_t));
    if (!hits) {
        free(sum_l2); free(ray_l);
        return -1;
    }

    for (int ray = 0; ray < n_rays; ray++) {
        /* Random unit direction */
        double phi = 2.0 * M_PI * lcg_double(&rng);
        double cos_theta = 2.0 * lcg_double(&rng) - 1.0;
        double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
        double ux = sin_theta * cos(phi);
        double uy = sin_theta * sin(phi);
        double uz = cos_theta;

        /* Build orthonormal frame */
        double tx, ty, tz;
        if (fabs(ux) < 0.9) { tx = 1; ty = 0; tz = 0; }
        else                 { tx = 0; ty = 1; tz = 0; }
        double v1x = uy * tz - uz * ty;
        double v1y = uz * tx - ux * tz;
        double v1z = ux * ty - uy * tx;
        double v1_len = sqrt(v1x * v1x + v1y * v1y + v1z * v1z);
        v1x /= v1_len; v1y /= v1_len; v1z /= v1_len;
        double v2x = uy * v1z - uz * v1y;
        double v2y = uz * v1x - ux * v1z;
        double v2z = ux * v1y - uy * v1x;

        /* Random disk point */
        double r_disk = R * sqrt(lcg_double(&rng));
        double angle = 2.0 * M_PI * lcg_double(&rng);
        double d1 = r_disk * cos(angle);
        double d2 = r_disk * sin(angle);

        /* Ray origin */
        double rox = cx + v1x * d1 + v2x * d2 - ux * 2.0 * R;
        double roy = cy + v1y * d1 + v2y * d2 - uy * 2.0 * R;
        double roz = cz + v1z * d1 + v2z * d2 - uz * 2.0 * R;

        if (ray_l) memset(ray_l, 0, n_instances * sizeof(double));

        alea_raycast_result_clear(&result);
        int rc = alea_raycast(sys, rox, roy, roz, ux, uy, uz, 4.0 * R, &result);
        if (rc != 0) continue;

        /* For each segment, find the matching instance via spatial query */
        for (size_t s = 0; s < result.segment_count; s++) {
            int seg_cell_id = result.segments[s].cell_id;
            if (seg_cell_id < 0) continue;
            double len = result.segments[s].t_exit - result.segments[s].t_enter;
            if (len <= 0) continue;

            /* Query at segment midpoint */
            double t_mid = (result.segments[s].t_enter + result.segments[s].t_exit) * 0.5;
            double px = rox + t_mid * ux;
            double py = roy + t_mid * uy;
            double pz = roz + t_mid * uz;

            int n_hits = alea_spatial_query_point(sys, px, py, pz, hits, max_hits);
            if (n_hits <= 0) continue;

            /* Find deepest terminal instance matching the segment's cell ID */
            int best_idx = -1;
            int best_depth = -1;
            for (int h = 0; h < n_hits; h++) {
                if (!hits[h].is_terminal) continue;
                if (hits[h].cell_id != seg_cell_id) continue;
                if (hits[h].depth > best_depth) {
                    best_depth = hits[h].depth;
                    best_idx = (int)hits[h].instance_index;
                }
            }

            if (best_idx >= 0 && (size_t)best_idx < n_instances) {
                volumes[best_idx] += len;
                if (ray_l) ray_l[best_idx] += len;
            }
        }

        /* Accumulate L^2 for this ray */
        if (sum_l2) {
            for (size_t i = 0; i < n_instances; i++) {
                sum_l2[i] += ray_l[i] * ray_l[i];
            }
        }
    }

    free(hits);
    alea_raycast_result_free(&result);

    /* Convert to volumes and relative errors: V = pi * R^2 * sum_L / N */
    double scale = M_PI * R * R / (double)n_rays;
    for (size_t i = 0; i < n_instances; i++) {
        double mean_l = volumes[i] / (double)n_rays;
        if (rel_errors) {
            if (mean_l > 0.0) {
                double mean_l2 = sum_l2[i] / (double)n_rays;
                double var_l = mean_l2 - mean_l * mean_l;
                if (var_l < 0.0) var_l = 0.0;
                rel_errors[i] = sqrt(var_l) / (mean_l * sqrt((double)n_rays));
            } else {
                rel_errors[i] = -1.0;
            }
        }
        volumes[i] *= scale;
    }

    free(sum_l2);
    free(ray_l);
    return 0;
}
