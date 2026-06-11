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
 * alea_estimate_path_volumes) are strong-symbol overrides of the weak stubs in
 * alea_module_stubs.c. The public alea_* wrappers live in alea_public_api.c.
 */

#define _USE_MATH_DEFINES
#include "alea.h"
#include "raycast.h"
#include "bvh.h"
#include "core/alea_system.h"
#include "primitives/bbox.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util/math.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(__GNUC__) || defined(__clang__)
#define ALEA_RAYCAST_DEPRECATED_CALL_BEGIN \
    _Pragma("GCC diagnostic push") \
    _Pragma("GCC diagnostic ignored \"-Wdeprecated-declarations\"")
#define ALEA_RAYCAST_DEPRECATED_CALL_END \
    _Pragma("GCC diagnostic pop")
#else
#define ALEA_RAYCAST_DEPRECATED_CALL_BEGIN
#define ALEA_RAYCAST_DEPRECATED_CALL_END
#endif


alea_raycast_result_t* alea_raycast_result_create(void) {
    alea_raycast_result_t* result = calloc(1, sizeof(alea_raycast_result_t));
    if (result) {
        alea_raycast_result_init(result);
    }
    return result;
}

void alea_raycast_result_destroy(alea_raycast_result_t* result) {
    if (!result) return;
    alea_raycast_result_free(result);
    free(result);
}

size_t alea_raycast_segment_count(const alea_raycast_result_t* result) {
    return result ? result->segments.count : 0;
}

int alea_raycast_segment_get(const alea_raycast_result_t* result, size_t index,
                                 double* t_enter, double* t_exit,
                                 int* cell_id, int* material_id, double* density,
                                 int* enter_surface_id, int* exit_surface_id) {
    if (!result || index >= result->segments.count) return -1;
    alea_ray_segment_t* seg = &result->segments.data[index];
    if (t_enter) *t_enter = seg->t_enter;
    if (t_exit) *t_exit = seg->t_exit;
    if (cell_id) *cell_id = seg->cell_id;
    if (material_id) *material_id = seg->material_id;
    if (density) *density = seg->density;
    if (enter_surface_id) *enter_surface_id = seg->enter_surface_id;
    if (exit_surface_id) *exit_surface_id = seg->exit_surface_id;
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

/**
 * Generate a random Cauchy-Crofton ray for volume estimation.
 * Picks a uniform random direction and a random point on a disk of radius R
 * perpendicular to that direction, centered at (cx,cy,cz).
 */
static void generate_cauchy_crofton_ray(uint32_t* rng,
                                        double cx, double cy, double cz, double R,
                                        double* rox, double* roy, double* roz,
                                        double* ux, double* uy, double* uz) {
    double phi = 2.0 * M_PI * lcg_double(rng);
    double cos_theta = 2.0 * lcg_double(rng) - 1.0;
    double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
    *ux = sin_theta * cos(phi);
    *uy = sin_theta * sin(phi);
    *uz = cos_theta;

    /* Build orthonormal frame (v1, v2) perpendicular to u */
    double tx, ty, tz;
    if (fabs(*ux) < 0.9) { tx = 1; ty = 0; tz = 0; }
    else                  { tx = 0; ty = 1; tz = 0; }
    double v1x = *uy * tz - *uz * ty;
    double v1y = *uz * tx - *ux * tz;
    double v1z = *ux * ty - *uy * tx;
    double v1_len = sqrt(v1x * v1x + v1y * v1y + v1z * v1z);
    v1x /= v1_len; v1y /= v1_len; v1z /= v1_len;
    double v2x = *uy * v1z - *uz * v1y;
    double v2y = *uz * v1x - *ux * v1z;
    double v2z = *ux * v1y - *uy * v1x;

    double r_disk = R * sqrt(lcg_double(rng));
    double angle = 2.0 * M_PI * lcg_double(rng);
    double d1 = r_disk * cos(angle);
    double d2 = r_disk * sin(angle);

    *rox = cx + v1x * d1 + v2x * d2 - *ux * 2.0 * R;
    *roy = cy + v1y * d1 + v2y * d2 - *uy * 2.0 * R;
    *roz = cz + v1z * d1 + v2z * d2 - *uz * 2.0 * R;
}

/**
 * Compute volume errors from accumulated track lengths and sum-of-squares.
 */
/**
 * Compute relative errors from raw sum_L (volumes before scaling) and sum_L^2.
 * Called before volumes[] are multiplied by scale.
 */
static void compute_volume_errors(const double* volumes, const double* sum_l2,
                                  double* rel_errors, size_t count, int n_rays) {
    for (size_t i = 0; i < count; i++) {
        double mean_l = volumes[i] / (double)n_rays;
        if (mean_l > 0.0) {
            double mean_l2 = sum_l2[i] / (double)n_rays;
            double var_l = mean_l2 - mean_l * mean_l;
            if (var_l < 0.0) var_l = 0.0;
            rel_errors[i] = sqrt(var_l) / (mean_l * sqrt((double)n_rays));
        } else {
            rel_errors[i] = -1.0;
        }
    }
}

static int compute_path_bounding_sphere(alea_system_t* sys,
                                        size_t n_paths,
                                        double* cx,
                                        double* cy,
                                        double* cz,
                                        double* radius) {
    if (!sys || n_paths == 0 || !cx || !cy || !cz || !radius) return -1;

    alea_volume_path_t* paths = calloc(n_paths, sizeof(*paths));
    if (!paths) return -1;
    size_t got = alea_volume_paths_get(sys, paths, n_paths);
    if (got > n_paths) got = n_paths;

    alea_bbox_t bounds = alea_bbox_empty();
    int bounded = 0;
    for (size_t i = 0; i < got; i++) {
        int cell_index = paths[i].terminal_cell_index;
        if (cell_index < 0 || (size_t)cell_index >= alea_vec_count(&sys->cells))
            continue;
        const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
        if (cell->root_node_id == ALEA_NODE_ID_INVALID ||
            cell->root_node_id >= alea_vec_count(&sys->nodes)) {
            continue;
        }

        const alea_bbox_t local_v = alea_node_bbox_get(&sys->nodes.data[cell->root_node_id].bbox);
        const alea_bbox_t* local = &local_v;
        if (!alea_bbox_is_valid(local)) continue;
        double dx = local->max_x - local->min_x;
        double dy = local->max_y - local->min_y;
        double dz = local->max_z - local->min_z;
        if (!isfinite(dx) || !isfinite(dy) || !isfinite(dz) ||
            dx > 9e5 || dy > 9e5 || dz > 9e5) {
            continue;
        }

        alea_matrix_t transform;
        memset(&transform, 0, sizeof(transform));
        memcpy(transform.m, paths[i].world_to_local, sizeof(transform.m));
        alea_bbox_t world = alea_bbox_transform(local, &transform);
        if (!alea_bbox_is_valid(&world)) continue;

        bounds = bounded ? alea_bbox_union(&bounds, &world) : world;
        bounded++;
    }
    free(paths);

    if (bounded == 0) return -1;

    *cx = (bounds.min_x + bounds.max_x) * 0.5;
    *cy = (bounds.min_y + bounds.max_y) * 0.5;
    *cz = (bounds.min_z + bounds.max_z) * 0.5;

    double dx = bounds.max_x - bounds.min_x;
    double dy = bounds.max_y - bounds.min_y;
    double dz = bounds.max_z - bounds.min_z;
    *radius = 0.5 * sqrt(dx * dx + dy * dy + dz * dz) * 1.01;
    return *radius > 0.0 ? 0 : -1;
}

int alea_estimate_cell_volumes(alea_system_t* sys,
                              double ox, double oy, double oz,
                              double radius, int n_rays,
                              double* volumes, double* rel_errors) {
    if (!sys || radius <= 0.0 || n_rays <= 0 || !volumes) return -1;

    /* Ensure all raycast caches are built before parallel section */
    if (alea_raycast_ensure_caches(sys) != 0) return -1;

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

    int error_flag = 0;

    /* Trace mu-random rays (Cauchy-Crofton):
     *   Volume formula: V_cell = pi * R^2 * (sum_track_length / n_rays)
     *   Error:  sigma_V = pi * R^2 / sqrt(N) * sqrt(E[L^2] - E[L]^2)
     */
    #pragma omp parallel
    {
        /* Thread-local accumulators */
        double* local_vol = calloc(n_cells, sizeof(double));
        double* local_l2 = rel_errors ? calloc(n_cells, sizeof(double)) : NULL;
        double* ray_l = rel_errors ? calloc(n_cells, sizeof(double)) : NULL;
        alea_raycast_result_t result;
        alea_raycast_result_init(&result);

        if (!local_vol || (rel_errors && (!local_l2 || !ray_l))) {
            #pragma omp atomic
            error_flag |= 1;
        }

        /* Per-thread deterministic seed */
        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif
        uint32_t rng = 42 + (uint32_t)tid * 2654435761u;

        #pragma omp for schedule(dynamic, 16)
        for (int ray = 0; ray < n_rays; ray++) {
            if (error_flag) continue;

            double rox, roy, roz, ux, uy, uz;
            generate_cauchy_crofton_ray(&rng, cx, cy, cz, R,
                                        &rox, &roy, &roz, &ux, &uy, &uz);

            if (ray_l) memset(ray_l, 0, n_cells * sizeof(double));

            alea_raycast_result_clear(&result);
            int rc = alea_raycast_hier_cell_aware(sys, rox, roy, roz,
                                                  ux, uy, uz, 4.0 * R,
                                                  &result);
            if (rc != 0) continue;

            /* Accumulate track lengths per cell */
            for (size_t s = 0; s < result.segments.count; s++) {
                int seg_cell_id = result.segments.data[s].cell_id;
                if (seg_cell_id < 0) continue;
                double len = result.segments.data[s].t_exit - result.segments.data[s].t_enter;
                if (len <= 0) continue;

                int ci = alea_find_cell_by_id(sys, seg_cell_id);
                if (ci >= 0 && (size_t)ci < n_cells) {
                    local_vol[ci] += len;
                    if (ray_l) ray_l[ci] += len;
                }
            }

            /* Accumulate L^2 for this ray */
            if (local_l2) {
                for (size_t ci = 0; ci < n_cells; ci++) {
                    local_l2[ci] += ray_l[ci] * ray_l[ci];
                }
            }
        }

        /* Merge thread-local results into global arrays */
        #pragma omp critical
        {
            if (local_vol) {
                for (size_t i = 0; i < n_cells; i++)
                    volumes[i] += local_vol[i];
            }
            if (local_l2 && sum_l2) {
                for (size_t i = 0; i < n_cells; i++)
                    sum_l2[i] += local_l2[i];
            }
        }

        alea_raycast_result_free(&result);
        free(local_vol);
        free(local_l2);
        free(ray_l);
    }

    if (error_flag) {
        free(sum_l2);
        return -1;
    }

    double scale = M_PI * R * R / (double)n_rays;
    if (rel_errors) {
        compute_volume_errors(volumes, sum_l2, rel_errors, n_cells, n_rays);
    }
    for (size_t i = 0; i < n_cells; i++) {
        volumes[i] *= scale;
    }

    free(sum_l2);
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
            if (sys->on_cell_removed)
                sys->on_cell_removed(sys->cell_hook_userdata, i);
            free(sys->cells.data[i].surface_indices);
            if (!sys->neighbor_pool)
                free(sys->cells.data[i].neighbors);
            free(sys->cells.data[i].lat_fill);
            free(sys->cells.data[i].comments);
            free(sys->cells.data[i].inline_comment);
            removed++;
        } else {
            if (write != i)
                sys->cells.data[write] = sys->cells.data[i];
            write++;
        }
    }
    sys->cells.count = write;

    if (removed > 0) {
        cell_hashmap_clear(&sys->cell_index);
        for (size_t i = 0; i < write; i++) {
            cell_hashmap_put(&sys->cell_index,
                             sys->cells.data[i].mc_cell_id,
                             (int)i);
            free(sys->cells.data[i].surface_indices);
            sys->cells.data[i].surface_indices = NULL;
            sys->cells.data[i].surface_index_count = 0;
            if (!sys->neighbor_pool)
                free(sys->cells.data[i].neighbors);
            sys->cells.data[i].neighbors = NULL;
            sys->cells.data[i].neighbor_count = 0;
        }

        alea_system_invalidate_query_caches(sys, ALEA_CACHE_ALL);
    }

    return removed;
}

/* ============================================================================
 * PER-INSTANCE VOLUME ESTIMATION
 *
 * Like alea_estimate_cell_volumes but resolves each ray segment to a specific
 * cell instance via the spatial index.  This distinguishes the same cell
 * appearing in multiple fill contexts.
 * ============================================================================ */

int alea_estimate_path_volumes(alea_system_t* sys,
                               int n_rays,
                               double* volumes,
                               double* rel_errors) {
    if (!sys || n_rays <= 0 || !volumes) return -1;

    alea_error_clear();
    size_t n_paths = alea_volume_path_count(sys);
    if (n_paths == 0 && alea_error_code() != (int)ALEA_OK) return -1;
    if (n_paths == 0) return 0;

    if (alea_raycast_ensure_hier_caches(sys) != 0) return -1;

    double cx, cy, cz, R;
    if (compute_path_bounding_sphere(sys, n_paths, &cx, &cy, &cz, &R) != 0 ||
        R <= 0.0) {
        if (alea_compute_bounding_sphere(sys, 1.0, &cx, &cy, &cz, &R) != 0 ||
            R <= 0.0) {
            return -1;
        }
    }
    R *= 1.01;

    memset(volumes, 0, n_paths * sizeof(double));
    if (rel_errors) memset(rel_errors, 0, n_paths * sizeof(double));

    double* sum_l2 = NULL;
    if (rel_errors) {
        sum_l2 = calloc(n_paths, sizeof(double));
        if (!sum_l2) return -1;
    }

    int error_flag = 0;

    #pragma omp parallel
    {
        double* local_vol = calloc(n_paths, sizeof(double));
        double* local_l2 = rel_errors ? calloc(n_paths, sizeof(double)) : NULL;
        double* ray_l = rel_errors ? calloc(n_paths, sizeof(double)) : NULL;
        alea_raycast_result_t result;
        alea_raycast_result_init(&result);

        if (!local_vol || (rel_errors && (!local_l2 || !ray_l))) {
            #pragma omp atomic
            error_flag |= 1;
        }

        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif
        uint32_t rng = 42 + (uint32_t)tid * 2654435761u;

        #pragma omp for schedule(dynamic, 16)
        for (int ray = 0; ray < n_rays; ray++) {
            if (error_flag) continue;

            double rox, roy, roz, ux, uy, uz;
            generate_cauchy_crofton_ray(&rng, cx, cy, cz, R,
                                        &rox, &roy, &roz, &ux, &uy, &uz);

            if (ray_l) memset(ray_l, 0, n_paths * sizeof(double));

            alea_raycast_result_clear(&result);
            int rc = alea_raycast(sys, rox, roy, roz, ux, uy, uz, 4.0 * R, &result);
            if (rc != 0) continue;

            for (size_t s = 0; s < result.segments.count; s++) {
                int seg_cell_id = result.segments.data[s].cell_id;
                if (seg_cell_id < 0) continue;
                double len = result.segments.data[s].t_exit -
                             result.segments.data[s].t_enter;
                if (len <= 0.0) continue;

                double t_mid = (result.segments.data[s].t_enter +
                                result.segments.data[s].t_exit) * 0.5;
                double px = rox + t_mid * ux;
                double py = roy + t_mid * uy;
                double pz = roz + t_mid * uz;

                alea_volume_path_t path;
                int found = alea_volume_path_at_point(sys, px, py, pz, &path);
                if (found <= 0) continue;
                if (path.terminal_cell_id != seg_cell_id) continue;
                if (path.path_id >= n_paths) continue;

                local_vol[path.path_id] += len;
                if (ray_l) ray_l[path.path_id] += len;
            }

            if (local_l2) {
                for (size_t i = 0; i < n_paths; i++) {
                    local_l2[i] += ray_l[i] * ray_l[i];
                }
            }
        }

        #pragma omp critical
        {
            if (local_vol) {
                for (size_t i = 0; i < n_paths; i++)
                    volumes[i] += local_vol[i];
            }
            if (local_l2 && sum_l2) {
                for (size_t i = 0; i < n_paths; i++)
                    sum_l2[i] += local_l2[i];
            }
        }

        alea_raycast_result_free(&result);
        free(local_vol);
        free(local_l2);
        free(ray_l);
    }

    if (error_flag) {
        free(sum_l2);
        return -1;
    }

    double scale = M_PI * R * R / (double)n_rays;
    if (rel_errors) {
        compute_volume_errors(volumes, sum_l2, rel_errors, n_paths, n_rays);
    }
    for (size_t i = 0; i < n_paths; i++) {
        volumes[i] *= scale;
    }

    free(sum_l2);
    return 0;
}
