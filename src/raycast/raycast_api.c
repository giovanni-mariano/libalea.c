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
#include "alea_raycast.h"
#include "alea_slice.h"
#include "raycast.h"
#include "ray_intersect.h"
#include "ray_epsilon.h"
#include "bvh.h"
#include "core/alea_system.h"
#include "primitives/bbox.h"
#include "primitives/primitive_eval.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <limits.h>
#include "util/math.h"

#ifdef _OPENMP
#include <omp.h>
#endif

int alea_surface_project_along(const alea_system_t* sys, int surface_id,
                               const double point[3], const double direction[3],
                               double* parameter, double projected[3],
                               alea_primitive_type_t* primitive_type) {
    if (!sys || !point || !direction || !parameter || !projected) return -1;
    int surface_index = alea_surface_find(sys, surface_id);
    if (surface_index < 0) return -1;
    const alea_surface_entry_t* surface = &sys->surfaces.data[surface_index];
    if (surface->primitive_id >= alea_vec_count(&sys->primitives)) return -1;

    alea_primitive_data_t data;
    if (!alea_primitive_copy_data(sys, surface->primitive_id, &data)) return -1;
    alea_primitive_type_t type = sys->primitives.data[surface->primitive_id].type;
    if (primitive_type) *primitive_type = type;

    alea_ray_t forward;
    if (alea_ray_init(&forward, point[0], point[1], point[2],
                      direction[0], direction[1], direction[2]) != 0) return -1;

    /* MCNP tori can have unequal radial and axial semiwidths.  The general
     * ray/torus quartic currently models a circular minor section, so use the
     * canonical primitive evaluator for the small receipt correction here. */
    if (type == ALEA_PRIMITIVE_TORUS_X ||
        type == ALEA_PRIMITIVE_TORUS_Y ||
        type == ALEA_PRIMITIVE_TORUS_Z) {
        double t = 0.0;
        for (int iteration = 0; iteration < 12; ++iteration) {
            double x = point[0] + t * forward.dx;
            double y = point[1] + t * forward.dy;
            double z = point[2] + t * forward.dz;
            double value = alea_primitive_eval(type, &data, x, y, z);
            if (!isfinite(value)) return 0;
            if (value == 0.0) break;
            double scale = fmax(1.0, fmax(fabs(x), fmax(fabs(y), fabs(z))));
            double h = 1.0e-8 * scale;
            double plus = alea_primitive_eval(
                type, &data, x + h*forward.dx, y + h*forward.dy,
                z + h*forward.dz);
            double minus = alea_primitive_eval(
                type, &data, x - h*forward.dx, y - h*forward.dy,
                z - h*forward.dz);
            double derivative = (plus - minus) / (2.0 * h);
            if (!isfinite(derivative) || fabs(derivative) < 1.0e-14) return 0;
            double step = value / derivative;
            if (!isfinite(step)) return 0;
            t -= step;
            if (fabs(step) <= 1.0e-13 * fmax(1.0, fabs(t))) break;
        }
        *parameter = t;
        projected[0] = point[0] + t * forward.dx;
        projected[1] = point[1] + t * forward.dy;
        projected[2] = point[2] + t * forward.dz;
        double final_value = alea_primitive_eval(
            type, &data, projected[0], projected[1], projected[2]);
        if (!isfinite(final_value) || fabs(final_value) > 1.0e-9) return 0;
        return 1;
    }

    double best = INFINITY;
    for (int reverse = 0; reverse < 2; ++reverse) {
        alea_ray_t ray;
        if (alea_ray_init(&ray, point[0], point[1], point[2],
                          reverse ? -forward.dx : forward.dx,
                          reverse ? -forward.dy : forward.dy,
                          reverse ? -forward.dz : forward.dz) != 0) return -1;
        double roots[4];
        int count = ray_intersect_primitive(&ray, type, &data, roots);
        for (int index = 0; index < count; ++index) {
            double candidate = reverse ? -roots[index] : roots[index];
            if (isfinite(candidate) && fabs(candidate) < fabs(best))
                best = candidate;
        }
    }
    if (!isfinite(best)) return 0;
    *parameter = best;
    projected[0] = point[0] + best * forward.dx;
    projected[1] = point[1] + best * forward.dy;
    projected[2] = point[2] + best * forward.dz;
    return 1;
}

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

struct alea_raycast_batch_result {
    size_t ray_count;
    size_t segment_count;
    uint32_t fields;
    uint64_t* ray_offsets;
    double* t_enter;
    double* t_exit;
    int32_t* cell_ids;
    int32_t* material_ids;
    double* densities;
    int32_t* enter_surface_ids;
    int32_t* exit_surface_ids;
    uint8_t* resolution_flags;
    int32_t* projected_cell_ids;
    int32_t* projected_material_ids;
    int32_t* projected_universe_ids;
    int32_t* projected_fill_universes;
    int32_t* projected_depths;
    uint8_t* projected_is_lattice;
    uint64_t* projected_occurrence_keys;
    size_t path_entry_count;
    uint64_t* segment_path_offsets;
    int32_t* path_cell_ids;
    int32_t* path_material_ids;
    int32_t* path_universe_ids;
    int32_t* path_fill_universes;
    int32_t* path_depths;
    uint8_t* path_is_lattice;
    double* path_lattice_origins_xyz;
    uint64_t* path_occurrence_keys;
    struct {
        uint8_t valid;
        uint64_t system_id;
        uint64_t geometry_generation;
        double origin[3];
        double u_axis[3];
        double v_axis[3];
        double u_min, u_max, v_min, v_max;
        size_t row_count;
        int projected_depth;
    } fast_slice_cache;
    alea_raycast_batch_work_stats_t work_stats;
};

typedef struct {
    alea_raycast_result_t trace;
    int status;
} alea_batch_trace_tmp_t;

/* Common segment products retain only segments in a worker-owned arena.  The
 * richer per-ray trace remains necessary for projected-owner and full-path
 * products, which depend on the captured hierarchy path. */
typedef struct {
    alea_raycast_result_t trace;
    alea_ray_segment_vec_t arena;
    uint64_t* row_offsets;
    size_t row_count;
    int status;
    alea_raycast_batch_work_stats_t work_stats;
} alea_batch_segment_worker_t;


static void batch_result_free_buffers(alea_raycast_batch_result_t* result) {
    if (!result) return;
    free(result->ray_offsets);
    free(result->t_enter);
    free(result->t_exit);
    free(result->cell_ids);
    free(result->material_ids);
    free(result->densities);
    free(result->enter_surface_ids);
    free(result->exit_surface_ids);
    free(result->resolution_flags);
    free(result->projected_cell_ids);
    free(result->projected_material_ids);
    free(result->projected_universe_ids);
    free(result->projected_fill_universes);
    free(result->projected_depths);
    free(result->projected_is_lattice);
    free(result->projected_occurrence_keys);
    free(result->segment_path_offsets);
    free(result->path_cell_ids);
    free(result->path_material_ids);
    free(result->path_universe_ids);
    free(result->path_fill_universes);
    free(result->path_depths);
    free(result->path_is_lattice);
    free(result->path_lattice_origins_xyz);
    free(result->path_occurrence_keys);
    memset(result, 0, sizeof(*result));
}

static void* batch_alloc_array(size_t count, size_t element_size) {
    if (count != 0 && element_size > SIZE_MAX / count) return NULL;
    return malloc(count ? count * element_size : 1);
}

static int batch_add_output_bytes(size_t* total, size_t count, size_t element_size) {
    if (count != 0 && element_size > SIZE_MAX / count) return -1;
    size_t bytes = count * element_size;
    if (bytes > SIZE_MAX - *total) return -1;
    *total += bytes;
    return 0;
}

static uint64_t batch_trace_retained_bytes(const alea_raycast_result_t* trace) {
    size_t bytes = 0;
    const struct { size_t count, size; } buffers[] = {
        { trace->hits.capacity, sizeof(*trace->hits.data) },
        { trace->segments.capacity, sizeof(*trace->segments.data) },
        { trace->paths.capacity, sizeof(*trace->paths.data) },
        { trace->path_entries.capacity, sizeof(*trace->path_entries.data) }
    };
    for (size_t i = 0; i < sizeof(buffers) / sizeof(buffers[0]); i++) {
        if (buffers[i].count > SIZE_MAX / buffers[i].size ||
            bytes > SIZE_MAX - buffers[i].count * buffers[i].size)
            return UINT64_MAX;
        bytes += buffers[i].count * buffers[i].size;
    }
    return bytes > UINT64_MAX ? UINT64_MAX : (uint64_t)bytes;
}

static void batch_work_stats_add_staging(alea_raycast_batch_work_stats_t* stats,
                                         const alea_raycast_result_t* trace) {
    const uint64_t retained = batch_trace_retained_bytes(trace);
    if (UINT64_MAX - stats->total_result_buffer_growths <
        trace->result_buffer_growths)
        stats->total_result_buffer_growths = UINT64_MAX;
    else
        stats->total_result_buffer_growths += trace->result_buffer_growths;
    if (UINT64_MAX - stats->total_result_buffer_growth_bytes <
        trace->result_buffer_growth_bytes)
        stats->total_result_buffer_growth_bytes = UINT64_MAX;
    else
        stats->total_result_buffer_growth_bytes += trace->result_buffer_growth_bytes;
    if (UINT64_MAX - stats->peak_trace_staging_bytes < retained)
        stats->peak_trace_staging_bytes = UINT64_MAX;
    else
        stats->peak_trace_staging_bytes += retained;
}


static int batch_validate_ray_inputs(const double* origins_xyz,
                                     const double* directions_xyz,
                                     size_t ray_count) {
    if (ray_count == 0) return 0;
    if (!origins_xyz || !directions_xyz) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "batch ray origins and directions are required");
        return -1;
    }
    if (ray_count > SIZE_MAX / 3) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW, "batch ray count is too large");
        return -1;
    }
    for (size_t i = 0; i < ray_count; i++) {
        const double* o = &origins_xyz[i * 3];
        const double* d = &directions_xyz[i * 3];
        if (!isfinite(o[0]) || !isfinite(o[1]) || !isfinite(o[2]) ||
            !isfinite(d[0]) || !isfinite(d[1]) || !isfinite(d[2]) ||
            (d[0] == 0.0 && d[1] == 0.0 && d[2] == 0.0)) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "batch ray %zu has non-finite data or zero direction", i);
            return -1;
        }
    }
    return 0;
}

static int batch_validate_slice_view(const alea_slice_view_t* view) {
    const alea_slice_plane_t* plane;
    const double tolerance = 1e-8;
    if (!view) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG, "ray slice view is required");
        return -1;
    }
    plane = &view->plane;
    for (size_t i = 0; i < 3; i++) {
        if (!isfinite(plane->origin[i]) || !isfinite(plane->normal[i]) ||
            !isfinite(plane->u_axis[i]) || !isfinite(plane->v_axis[i])) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "ray slice view contains non-finite coordinates");
            return -1;
        }
    }
    if (!isfinite(view->u_min) || !isfinite(view->u_max) ||
        !isfinite(view->v_min) || !isfinite(view->v_max) ||
        view->u_min >= view->u_max || view->v_min >= view->v_max) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "ray slice viewport bounds must be finite and increasing");
        return -1;
    }
    double un = sqrt(plane->u_axis[0] * plane->u_axis[0] +
                     plane->u_axis[1] * plane->u_axis[1] +
                     plane->u_axis[2] * plane->u_axis[2]);
    double vn = sqrt(plane->v_axis[0] * plane->v_axis[0] +
                     plane->v_axis[1] * plane->v_axis[1] +
                     plane->v_axis[2] * plane->v_axis[2]);
    double nn = sqrt(plane->normal[0] * plane->normal[0] +
                     plane->normal[1] * plane->normal[1] +
                     plane->normal[2] * plane->normal[2]);
    double uv = plane->u_axis[0] * plane->v_axis[0] +
                plane->u_axis[1] * plane->v_axis[1] +
                plane->u_axis[2] * plane->v_axis[2];
    double cross[3] = {
        plane->u_axis[1] * plane->v_axis[2] - plane->u_axis[2] * plane->v_axis[1],
        plane->u_axis[2] * plane->v_axis[0] - plane->u_axis[0] * plane->v_axis[2],
        plane->u_axis[0] * plane->v_axis[1] - plane->u_axis[1] * plane->v_axis[0]
    };
    double handedness = cross[0] * plane->normal[0] +
                        cross[1] * plane->normal[1] +
                        cross[2] * plane->normal[2];
    if (fabs(un - 1.0) > tolerance || fabs(vn - 1.0) > tolerance ||
        fabs(nn - 1.0) > tolerance || fabs(uv) > tolerance ||
        fabs(fabs(handedness) - 1.0) > tolerance) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "ray slice basis must be orthonormal and match its normal");
        return -1;
    }
    return 0;
}

static void batch_set_fast_slice_cache(alea_raycast_batch_result_t* result,
                                       const alea_system_t* sys,
                                       const alea_slice_view_t* view,
                                       size_t row_count,
                                       int projected_depth) {
    if (!result || !sys || !view) return;
    result->fast_slice_cache.valid = 1;
    result->fast_slice_cache.system_id = sys->system_id;
    result->fast_slice_cache.geometry_generation =
        alea_system_geometry_generation(sys);
    memcpy(result->fast_slice_cache.origin, view->plane.origin,
           sizeof(result->fast_slice_cache.origin));
    memcpy(result->fast_slice_cache.u_axis, view->plane.u_axis,
           sizeof(result->fast_slice_cache.u_axis));
    memcpy(result->fast_slice_cache.v_axis, view->plane.v_axis,
           sizeof(result->fast_slice_cache.v_axis));
    result->fast_slice_cache.u_min = view->u_min;
    result->fast_slice_cache.u_max = view->u_max;
    result->fast_slice_cache.v_min = view->v_min;
    result->fast_slice_cache.v_max = view->v_max;
    result->fast_slice_cache.row_count = row_count;
    result->fast_slice_cache.projected_depth = projected_depth;
}

int alea_raycast_batch_result_matches_fast_slice_cache(
    const alea_raycast_batch_result_t* result,
    const alea_system_t* sys,
    const alea_slice_view_t* view,
    size_t row_count,
    const alea_raycast_batch_options_t* render_options,
    int projected_depth) {
    uint32_t required_fields = render_options ? render_options->fields : 0;
    if (!result || !sys || !view || !result->fast_slice_cache.valid) return 0;
    if (result->fast_slice_cache.system_id != sys->system_id ||
        result->fast_slice_cache.geometry_generation !=
            alea_system_geometry_generation(sys) ||
        result->fast_slice_cache.row_count != row_count ||
        result->ray_count != row_count) return 0;
    if (memcmp(result->fast_slice_cache.origin, view->plane.origin,
               sizeof(result->fast_slice_cache.origin)) != 0 ||
        memcmp(result->fast_slice_cache.u_axis, view->plane.u_axis,
               sizeof(result->fast_slice_cache.u_axis)) != 0 ||
        memcmp(result->fast_slice_cache.v_axis, view->plane.v_axis,
               sizeof(result->fast_slice_cache.v_axis)) != 0 ||
        result->fast_slice_cache.u_min != view->u_min ||
        result->fast_slice_cache.u_max != view->u_max ||
        result->fast_slice_cache.v_min != view->v_min ||
        result->fast_slice_cache.v_max != view->v_max) return 0;
    if ((result->fields & required_fields) != required_fields) return 0;
    if (projected_depth >= 0 &&
        (!(result->fields & ALEA_RAY_BATCH_PROJECTED_OWNER) ||
         result->fast_slice_cache.projected_depth != projected_depth)) return 0;
    return result->ray_offsets && result->t_enter && result->t_exit &&
           result->cell_ids;
}

int alea_raycast_batch_result_get_compact_slice_provenance(
    const alea_raycast_batch_result_t* result,
    const alea_slice_view_t* view,
    size_t row_count,
    int* out_projected_depth) {
    if (out_projected_depth) *out_projected_depth = -1;
    if (!result || !view || !result->fast_slice_cache.valid ||
        result->fast_slice_cache.row_count != row_count ||
        result->ray_count != row_count ||
        !result->ray_offsets || !result->t_enter || !result->t_exit ||
        !result->cell_ids)
        return 0;
    if (memcmp(result->fast_slice_cache.origin, view->plane.origin,
               sizeof(result->fast_slice_cache.origin)) != 0 ||
        memcmp(result->fast_slice_cache.u_axis, view->plane.u_axis,
               sizeof(result->fast_slice_cache.u_axis)) != 0 ||
        memcmp(result->fast_slice_cache.v_axis, view->plane.v_axis,
               sizeof(result->fast_slice_cache.v_axis)) != 0 ||
        result->fast_slice_cache.u_min != view->u_min ||
        result->fast_slice_cache.u_max != view->u_max ||
        result->fast_slice_cache.v_min != view->v_min ||
        result->fast_slice_cache.v_max != view->v_max)
        return 0;
    if (out_projected_depth)
        *out_projected_depth = result->fast_slice_cache.projected_depth;
    return 1;
}

void alea_raycast_batch_result_swap_internal(
    alea_raycast_batch_result_t* a,
    alea_raycast_batch_result_t* b) {
    if (!a || !b) return;
    alea_raycast_batch_result_t tmp = *a;
    *a = *b;
    *b = tmp;
}

static void batch_store_projected_owner(alea_raycast_batch_result_t* result,
                                        size_t index,
                                        const alea_raycast_result_t* trace,
                                        const alea_ray_segment_t* segment,
                                        int requested_depth) {
    const alea_ray_path_entry_t* selected = NULL;
    result->projected_cell_ids[index] = -1;
    result->projected_material_ids[index] = 0;
    result->projected_universe_ids[index] = -1;
    result->projected_fill_universes[index] = -1;
    result->projected_depths[index] = -1;
    result->projected_is_lattice[index] = 0;
    result->projected_occurrence_keys[index] = 0;

    if (segment->path_index == UINT32_MAX ||
        segment->path_index >= trace->paths.count) return;
    const alea_ray_path_t* path = &trace->paths.data[segment->path_index];
    if (path->count == 0 || (size_t)path->offset + path->count > trace->path_entries.count)
        return;

    selected = &trace->path_entries.data[path->offset + path->count - 1];
    if (requested_depth >= 0) {
        const alea_ray_path_entry_t* best = NULL;
        for (size_t i = 0; i < path->count; i++) {
            const alea_ray_path_entry_t* entry =
                &trace->path_entries.data[path->offset + i];
            if (entry->depth <= requested_depth &&
                (!best || entry->depth > best->depth)) {
                best = entry;
            }
        }
        if (best) selected = best;
    }

    result->projected_cell_ids[index] = selected->cell_id;
    result->projected_material_ids[index] = selected->material_id;
    result->projected_universe_ids[index] = selected->universe_id;
    result->projected_fill_universes[index] = selected->fill_universe;
    result->projected_depths[index] = selected->depth;
    result->projected_is_lattice[index] = selected->is_lattice;
    result->projected_occurrence_keys[index] = selected->occurrence_key;
}

static const alea_ray_path_entry_t* batch_segment_path(
    const alea_raycast_result_t* trace, const alea_ray_segment_t* segment,
    size_t* out_count) {
    *out_count = 0;
    if (segment->path_index == UINT32_MAX ||
        segment->path_index >= trace->paths.count) return NULL;
    const alea_ray_path_t* path = &trace->paths.data[segment->path_index];
    if (path->count == 0 || (size_t)path->offset + path->count > trace->path_entries.count)
        return NULL;
    *out_count = path->count;
    return &trace->path_entries.data[path->offset];
}

alea_raycast_batch_result_t* alea_raycast_batch_result_create(void) {
    return calloc(1, sizeof(alea_raycast_batch_result_t));
}

void alea_raycast_batch_result_destroy(alea_raycast_batch_result_t* result) {
    if (!result) return;
    batch_result_free_buffers(result);
    free(result);
}

size_t alea_raycast_batch_ray_count(const alea_raycast_batch_result_t* result) {
    return result ? result->ray_count : 0;
}

size_t alea_raycast_batch_segment_count(const alea_raycast_batch_result_t* result) {
    return result ? result->segment_count : 0;
}

int alea_raycast_batch_result_get_work_stats_internal(
    const alea_raycast_batch_result_t* result,
    alea_raycast_batch_work_stats_t* out_stats) {
    if (!result || !out_stats) return -1;
    *out_stats = result->work_stats;
    return 0;
}

static void batch_work_stats_accumulate(
    alea_raycast_batch_work_stats_t* stats,
    const alea_raycast_result_t* trace) {
#define BATCH_WORK_MAX(field, source) \
    do { if (stats->field < trace->source) stats->field = trace->source; } while (0)
    BATCH_WORK_MAX(max_owner_neighbor_attempts, owner_neighbor_attempts);
    BATCH_WORK_MAX(max_owner_neighbor_hits, owner_neighbor_hits);
    BATCH_WORK_MAX(max_owner_path_attempts, owner_path_attempts);
    BATCH_WORK_MAX(max_owner_path_hits, owner_path_hits);
    BATCH_WORK_MAX(max_owner_root_queries, owner_root_queries);
    BATCH_WORK_MAX(max_owner_root_hits, owner_root_hits);
    BATCH_WORK_MAX(max_owner_full_queries, owner_full_queries);
    BATCH_WORK_MAX(max_owner_full_hits, owner_full_hits);
    BATCH_WORK_MAX(max_boundary_event_enrichments,
                   boundary_event_enrichments);
    BATCH_WORK_MAX(max_path_snapshot_copies, path_snapshot_copies);
    BATCH_WORK_MAX(max_path_snapshot_entries, path_snapshot_entries);
    BATCH_WORK_MAX(max_selected_intervals_yielded,
                   selected_intervals_yielded);
    BATCH_WORK_MAX(max_result_buffer_growths, result_buffer_growths);
    BATCH_WORK_MAX(max_result_buffer_growth_bytes, result_buffer_growth_bytes);
    BATCH_WORK_MAX(max_lattice_entry_calls, lattice_entry_calls);
    BATCH_WORK_MAX(max_lattice_entry_tlas_nodes_tested,
                   lattice_entry_tlas_nodes_tested);
    BATCH_WORK_MAX(max_lattice_entry_tlas_leaves_visited,
                   lattice_entry_tlas_leaves_visited);
    BATCH_WORK_MAX(max_lattice_entry_candidates, lattice_entry_candidates);
    BATCH_WORK_MAX(max_lattice_entry_dda_steps, lattice_entry_dda_steps);
    BATCH_WORK_MAX(max_lattice_entry_no_entry_results,
                   lattice_entry_no_entry_results);
    BATCH_WORK_MAX(max_lattice_entry_future_entry_results,
                   lattice_entry_future_entry_results);
    BATCH_WORK_MAX(max_lattice_entry_already_inside_results,
                   lattice_entry_already_inside_results);
    BATCH_WORK_MAX(max_lattice_entry_ancestor_surface_tests,
                   lattice_entry_ancestor_surface_tests);
    BATCH_WORK_MAX(max_lattice_entry_ancestor_events,
                   lattice_entry_ancestor_events);
    BATCH_WORK_MAX(max_lattice_entry_canonical_rejections,
                   lattice_entry_canonical_rejections);
#undef BATCH_WORK_MAX
}

static void batch_work_stats_merge(alea_raycast_batch_work_stats_t* target,
                                   const alea_raycast_batch_work_stats_t* source) {
#define BATCH_WORK_MERGE(field) \
    do { if (target->field < source->field) target->field = source->field; } while (0)
    BATCH_WORK_MERGE(max_owner_neighbor_attempts);
    BATCH_WORK_MERGE(max_owner_neighbor_hits);
    BATCH_WORK_MERGE(max_owner_path_attempts);
    BATCH_WORK_MERGE(max_owner_path_hits);
    BATCH_WORK_MERGE(max_owner_root_queries);
    BATCH_WORK_MERGE(max_owner_root_hits);
    BATCH_WORK_MERGE(max_owner_full_queries);
    BATCH_WORK_MERGE(max_owner_full_hits);
    BATCH_WORK_MERGE(max_boundary_event_enrichments);
    BATCH_WORK_MERGE(max_path_snapshot_copies);
    BATCH_WORK_MERGE(max_path_snapshot_entries);
    BATCH_WORK_MERGE(max_selected_intervals_yielded);
    BATCH_WORK_MERGE(max_result_buffer_growths);
    BATCH_WORK_MERGE(max_result_buffer_growth_bytes);
    BATCH_WORK_MERGE(max_lattice_entry_calls);
    BATCH_WORK_MERGE(max_lattice_entry_tlas_nodes_tested);
    BATCH_WORK_MERGE(max_lattice_entry_tlas_leaves_visited);
    BATCH_WORK_MERGE(max_lattice_entry_candidates);
    BATCH_WORK_MERGE(max_lattice_entry_dda_steps);
    BATCH_WORK_MERGE(max_lattice_entry_no_entry_results);
    BATCH_WORK_MERGE(max_lattice_entry_future_entry_results);
    BATCH_WORK_MERGE(max_lattice_entry_already_inside_results);
    BATCH_WORK_MERGE(max_lattice_entry_ancestor_surface_tests);
    BATCH_WORK_MERGE(max_lattice_entry_ancestor_events);
    BATCH_WORK_MERGE(max_lattice_entry_canonical_rejections);
#undef BATCH_WORK_MERGE
}

static void batch_segment_worker_free(alea_batch_segment_worker_t* worker) {
    if (!worker) return;
    alea_raycast_result_free(&worker->trace);
    alea_vec_free(&worker->arena);
    free(worker->row_offsets);
    memset(worker, 0, sizeof(*worker));
}

static int batch_segment_worker_build(
    alea_system_t* sys, const double* origins_xyz, const double* directions_xyz,
    size_t ray_count, double scalar_t_max, const double* t_mins,
    const double* t_maxs, uint64_t max_segments,
    atomic_uint_fast64_t* live_segment_count,
    alea_batch_segment_worker_t* worker, size_t worker_index,
    size_t worker_count) {
    const size_t owned_rows = ray_count > worker_index
        ? 1 + (ray_count - 1 - worker_index) / worker_count : 0;
    worker->row_offsets = calloc(owned_rows + 1, sizeof(*worker->row_offsets));
    if (!worker->row_offsets) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate segment worker offsets");
        return -1;
    }
    worker->row_count = owned_rows;
    alea_raycast_result_init(&worker->trace);
    alea_vec_init(&worker->arena);
    for (size_t row = worker_index, local_row = 0; row < ray_count;
         row += worker_count, local_row++) {
        if (alea_interrupted()) return -1;
        alea_raycast_result_clear(&worker->trace);
        if (max_segments != 0) {
            worker->trace.segment_counter = live_segment_count;
            worker->trace.segment_limit = max_segments;
        }
        const double* o = &origins_xyz[row * 3];
        const double* d = &directions_xyz[row * 3];
        alea_ray_t ray;
        if (alea_ray_init(&ray, o[0], o[1], o[2], d[0], d[1], d[2]) != 0 ||
            alea_raycast_hier_segments_nocache(
                sys, &ray, t_maxs ? t_maxs[row] : scalar_t_max,
                &worker->trace) != 0) {
            worker->status = -1;
            return -1;
        }
        if (t_mins && t_mins[row] > 0.0) {
            const double t_min = t_mins[row];
            size_t write = 0;
            for (size_t j = 0; j < worker->trace.segments.count; j++) {
                alea_ray_segment_t segment = worker->trace.segments.data[j];
                if (segment.t_exit <= t_min + RAY_EPSILON) continue;
                if (segment.t_enter < t_min) {
                    segment.t_enter = t_min;
                    segment.enter_surface_id = -1;
                    segment.enter_hit_index = -1;
                }
                worker->trace.segments.data[write++] = segment;
            }
            worker->trace.segments.count = write;
        }
        worker->row_offsets[local_row] = worker->arena.count;
        for (size_t j = 0; j < worker->trace.segments.count; j++) {
            if (alea_vec_push(&worker->arena, worker->trace.segments.data[j],
                              alea_ray_segment_t) != 0) {
                alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                                      "failed to append segment worker arena");
                worker->status = -1;
                return -1;
            }
        }
        worker->row_offsets[local_row + 1] = worker->arena.count;
        batch_work_stats_accumulate(&worker->work_stats, &worker->trace);
    }
    return 0;
}

static int raycast_hier_batch_execute_common_arena(
    alea_system_t* sys, const double* origins_xyz, const double* directions_xyz,
    size_t ray_count, double scalar_t_max, const double* t_mins,
    const double* t_maxs, uint32_t fields, uint64_t max_segments,
    uint64_t max_output_bytes, atomic_uint_fast64_t* live_segment_count,
    alea_raycast_batch_result_t* result) {
    size_t worker_count = 1;
    alea_batch_segment_worker_t* workers = NULL;
    alea_raycast_batch_result_t next = {0};
#ifdef _OPENMP
    worker_count = (size_t)omp_get_max_threads();
#endif
    if (worker_count == 0) worker_count = 1;
    if (ray_count != 0 && worker_count > ray_count) worker_count = ray_count;
    if (worker_count > SIZE_MAX / sizeof(*workers)) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW, "batch worker storage overflows");
        return -1;
    }
    workers = calloc(worker_count, sizeof(*workers));
    if (!workers) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate segment worker arenas");
        return -1;
    }

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (size_t worker_index = 0; worker_index < worker_count; worker_index++)
        workers[worker_index].status = batch_segment_worker_build(
            sys, origins_xyz, directions_xyz, ray_count, scalar_t_max, t_mins,
            t_maxs, max_segments, live_segment_count, &workers[worker_index],
            worker_index, worker_count);

    for (size_t i = 0; i < worker_count; i++) {
        if (workers[i].status != 0 || alea_interrupted()) {
            if (alea_interrupted())
                alea_set_error_detail(ALEA_ERR_INTERRUPTED, "batch raycast interrupted");
            else if (workers[i].trace.segment_limit_exceeded)
                alea_set_error_detail(ALEA_ERR_OVERFLOW,
                                      "batch segment limit (%llu) exceeded",
                                      (unsigned long long)max_segments);
            else if (alea_error_code() == ALEA_OK)
                alea_set_error_detail(ALEA_ERR_INVALID_STATE, "batch raycast failed");
            goto cleanup;
        }
    }

    next.ray_count = ray_count;
    next.fields = fields;
    next.ray_offsets = batch_alloc_array(ray_count + 1, sizeof(*next.ray_offsets));
    if (!next.ray_offsets) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "failed to allocate batch offsets");
        goto cleanup;
    }
    next.ray_offsets[0] = 0;
    for (size_t row = 0; row < ray_count; row++) {
        const size_t worker_index = row % worker_count;
        const size_t local_row = row / worker_count;
        const alea_batch_segment_worker_t* worker = &workers[worker_index];
        const uint64_t count = worker->row_offsets[local_row + 1] -
                               worker->row_offsets[local_row];
        if (count > UINT64_MAX - next.ray_offsets[row]) {
            alea_set_error_detail(ALEA_ERR_OVERFLOW, "batch segment offsets overflow");
            goto cleanup;
        }
        next.ray_offsets[row + 1] = next.ray_offsets[row] + count;
    }
    if (next.ray_offsets[ray_count] > SIZE_MAX) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW, "batch segment count overflows size_t");
        goto cleanup;
    }
    next.segment_count = (size_t)next.ray_offsets[ray_count];
    if (max_segments != 0 && next.ray_offsets[ray_count] > max_segments) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW, "batch segment limit (%llu) exceeded",
                              (unsigned long long)max_segments);
        goto cleanup;
    }
    if (max_output_bytes != 0) {
        size_t output_bytes = 0;
        int overflow =
            batch_add_output_bytes(&output_bytes, ray_count + 1,
                                   sizeof(*next.ray_offsets)) ||
            batch_add_output_bytes(&output_bytes, next.segment_count,
                                   sizeof(*next.t_enter)) ||
            batch_add_output_bytes(&output_bytes, next.segment_count,
                                   sizeof(*next.t_exit)) ||
            batch_add_output_bytes(&output_bytes, next.segment_count,
                                   sizeof(*next.cell_ids));
        if (fields & ALEA_RAY_BATCH_MATERIAL)
            overflow |= batch_add_output_bytes(&output_bytes, next.segment_count,
                                               sizeof(*next.material_ids));
        if (fields & ALEA_RAY_BATCH_DENSITY)
            overflow |= batch_add_output_bytes(&output_bytes, next.segment_count,
                                               sizeof(*next.densities));
        if (fields & ALEA_RAY_BATCH_SURFACES) {
            overflow |= batch_add_output_bytes(&output_bytes, next.segment_count,
                                               sizeof(*next.enter_surface_ids));
            overflow |= batch_add_output_bytes(&output_bytes, next.segment_count,
                                               sizeof(*next.exit_surface_ids));
        }
        if (fields & ALEA_RAY_BATCH_RESOLUTION_FLAGS)
            overflow |= batch_add_output_bytes(&output_bytes, next.segment_count,
                                               sizeof(*next.resolution_flags));
        if (overflow || output_bytes > max_output_bytes) {
            alea_set_error_detail(ALEA_ERR_OVERFLOW,
                                  "batch compact output exceeds byte limit (%llu)",
                                  (unsigned long long)max_output_bytes);
            goto cleanup;
        }
    }
    next.t_enter = batch_alloc_array(next.segment_count, sizeof(*next.t_enter));
    next.t_exit = batch_alloc_array(next.segment_count, sizeof(*next.t_exit));
    next.cell_ids = batch_alloc_array(next.segment_count, sizeof(*next.cell_ids));
    if (fields & ALEA_RAY_BATCH_MATERIAL)
        next.material_ids = batch_alloc_array(next.segment_count, sizeof(*next.material_ids));
    if (fields & ALEA_RAY_BATCH_DENSITY)
        next.densities = batch_alloc_array(next.segment_count, sizeof(*next.densities));
    if (fields & ALEA_RAY_BATCH_SURFACES) {
        next.enter_surface_ids = batch_alloc_array(next.segment_count,
                                                    sizeof(*next.enter_surface_ids));
        next.exit_surface_ids = batch_alloc_array(next.segment_count,
                                                   sizeof(*next.exit_surface_ids));
    }
    if (fields & ALEA_RAY_BATCH_RESOLUTION_FLAGS)
        next.resolution_flags = batch_alloc_array(next.segment_count,
                                                  sizeof(*next.resolution_flags));
    if (!next.t_enter || !next.t_exit || !next.cell_ids ||
        ((fields & ALEA_RAY_BATCH_MATERIAL) && !next.material_ids) ||
        ((fields & ALEA_RAY_BATCH_DENSITY) && !next.densities) ||
        ((fields & ALEA_RAY_BATCH_SURFACES) &&
         (!next.enter_surface_ids || !next.exit_surface_ids)) ||
        ((fields & ALEA_RAY_BATCH_RESOLUTION_FLAGS) && !next.resolution_flags)) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate batch segment fields");
        goto cleanup;
    }

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (size_t row = 0; row < ray_count; row++) {
        const size_t worker_index = row % worker_count;
        const size_t local_row = row / worker_count;
        const alea_batch_segment_worker_t* worker = &workers[worker_index];
        const size_t begin = (size_t)worker->row_offsets[local_row];
        const size_t end = (size_t)worker->row_offsets[local_row + 1];
        size_t dst = (size_t)next.ray_offsets[row];
        for (size_t source = begin; source < end; source++, dst++) {
            const alea_ray_segment_t* seg = &worker->arena.data[source];
            next.t_enter[dst] = seg->t_enter;
            next.t_exit[dst] = seg->t_exit;
            next.cell_ids[dst] = (int32_t)seg->cell_id;
            if (next.material_ids) next.material_ids[dst] = (int32_t)seg->material_id;
            if (next.densities) next.densities[dst] = seg->density;
            if (next.enter_surface_ids) next.enter_surface_ids[dst] = (int32_t)seg->enter_surface_id;
            if (next.exit_surface_ids) next.exit_surface_ids[dst] = (int32_t)seg->exit_surface_id;
            if (next.resolution_flags) next.resolution_flags[dst] = seg->resolution_flags;
        }
    }
    for (size_t i = 0; i < worker_count; i++) {
        batch_work_stats_merge(&next.work_stats, &workers[i].work_stats);
        batch_segment_worker_free(&workers[i]);
    }
    free(workers);
    batch_result_free_buffers(result);
    *result = next;
    return 0;

cleanup:
    for (size_t i = 0; i < worker_count; i++) batch_segment_worker_free(&workers[i]);
    free(workers);
    batch_result_free_buffers(&next);
    return -1;
}

uint32_t alea_raycast_batch_fields(const alea_raycast_batch_result_t* result) {
    return result ? result->fields : 0;
}

const uint64_t* alea_raycast_batch_ray_offsets(const alea_raycast_batch_result_t* result) {
    return result ? result->ray_offsets : NULL;
}
const double* alea_raycast_batch_t_enter(const alea_raycast_batch_result_t* result) {
    return result ? result->t_enter : NULL;
}
const double* alea_raycast_batch_t_exit(const alea_raycast_batch_result_t* result) {
    return result ? result->t_exit : NULL;
}
const int32_t* alea_raycast_batch_cell_ids(const alea_raycast_batch_result_t* result) {
    return result ? result->cell_ids : NULL;
}
const int32_t* alea_raycast_batch_material_ids(const alea_raycast_batch_result_t* result) {
    return result ? result->material_ids : NULL;
}
const double* alea_raycast_batch_densities(const alea_raycast_batch_result_t* result) {
    return result ? result->densities : NULL;
}
const int32_t* alea_raycast_batch_enter_surface_ids(const alea_raycast_batch_result_t* result) {
    return result ? result->enter_surface_ids : NULL;
}
const int32_t* alea_raycast_batch_exit_surface_ids(const alea_raycast_batch_result_t* result) {
    return result ? result->exit_surface_ids : NULL;
}
const uint8_t* alea_raycast_batch_resolution_flags(const alea_raycast_batch_result_t* result) {
    return result ? result->resolution_flags : NULL;
}
const int32_t* alea_raycast_batch_projected_cell_ids(const alea_raycast_batch_result_t* result) {
    return result ? result->projected_cell_ids : NULL;
}
const int32_t* alea_raycast_batch_projected_material_ids(const alea_raycast_batch_result_t* result) {
    return result ? result->projected_material_ids : NULL;
}
const int32_t* alea_raycast_batch_projected_universe_ids(const alea_raycast_batch_result_t* result) {
    return result ? result->projected_universe_ids : NULL;
}
const int32_t* alea_raycast_batch_projected_fill_universes(const alea_raycast_batch_result_t* result) {
    return result ? result->projected_fill_universes : NULL;
}
const int32_t* alea_raycast_batch_projected_depths(const alea_raycast_batch_result_t* result) {
    return result ? result->projected_depths : NULL;
}
const uint8_t* alea_raycast_batch_projected_is_lattice(const alea_raycast_batch_result_t* result) {
    return result ? result->projected_is_lattice : NULL;
}
const uint64_t* alea_raycast_batch_projected_occurrence_keys(const alea_raycast_batch_result_t* result) {
    return result ? result->projected_occurrence_keys : NULL;
}
size_t alea_raycast_batch_path_entry_count(const alea_raycast_batch_result_t* result) {
    return result ? result->path_entry_count : 0;
}
const uint64_t* alea_raycast_batch_segment_path_offsets(const alea_raycast_batch_result_t* result) {
    return result ? result->segment_path_offsets : NULL;
}
const int32_t* alea_raycast_batch_path_cell_ids(const alea_raycast_batch_result_t* result) {
    return result ? result->path_cell_ids : NULL;
}
const int32_t* alea_raycast_batch_path_material_ids(const alea_raycast_batch_result_t* result) {
    return result ? result->path_material_ids : NULL;
}
const int32_t* alea_raycast_batch_path_universe_ids(const alea_raycast_batch_result_t* result) {
    return result ? result->path_universe_ids : NULL;
}
const int32_t* alea_raycast_batch_path_fill_universes(const alea_raycast_batch_result_t* result) {
    return result ? result->path_fill_universes : NULL;
}
const int32_t* alea_raycast_batch_path_depths(const alea_raycast_batch_result_t* result) {
    return result ? result->path_depths : NULL;
}
const uint8_t* alea_raycast_batch_path_is_lattice(const alea_raycast_batch_result_t* result) {
    return result ? result->path_is_lattice : NULL;
}
const double* alea_raycast_batch_path_lattice_origins_xyz(const alea_raycast_batch_result_t* result) {
    return result ? result->path_lattice_origins_xyz : NULL;
}
const uint64_t* alea_raycast_batch_path_occurrence_keys(const alea_raycast_batch_result_t* result) {
    return result ? result->path_occurrence_keys : NULL;
}

static int raycast_hier_batch_execute(
    alea_system_t* sys, const double* origins_xyz,
    const double* directions_xyz, size_t ray_count,
    double scalar_t_max, const double* t_mins, const double* t_maxs,
    const alea_raycast_batch_options_t* options,
    alea_raycast_batch_result_t* result) {
    const uint32_t known_fields = ALEA_RAY_BATCH_MATERIAL |
                                  ALEA_RAY_BATCH_DENSITY |
                                  ALEA_RAY_BATCH_SURFACES |
                                  ALEA_RAY_BATCH_RESOLUTION_FLAGS |
                                  ALEA_RAY_BATCH_PROJECTED_OWNER |
                                  ALEA_RAY_BATCH_FULL_PATHS;
    uint32_t fields = 0;
    uint64_t max_segments = 0;
    uint64_t max_path_entries = 0;
    uint64_t max_output_bytes = 0;
    atomic_uint_fast64_t live_segment_count;
    atomic_uint_fast64_t live_path_entry_count;
    alea_batch_trace_tmp_t* traces = NULL;
    alea_raycast_batch_result_t next = {0};

    if (!sys || !result) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "batch raycast requires a system and result");
        return -1;
    }
    if (!t_maxs && !isfinite(scalar_t_max)) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG, "batch t_max must be finite");
        return -1;
    }
    if (t_mins || t_maxs) {
        for (size_t i = 0; i < ray_count; i++) {
            const double t_min = t_mins ? t_mins[i] : 0.0;
            const double t_max = t_maxs ? t_maxs[i] : scalar_t_max;
            if (!isfinite(t_min) || !isfinite(t_max) || t_min < 0.0 ||
                (t_max > 0.0 && t_min > t_max)) {
                alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                      "batch ray %zu has an invalid t range", i);
                return -1;
            }
        }
    }
    if (options) {
        if (options->struct_size < sizeof(*options)) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "batch options struct_size is too small");
            return -1;
        }
        if (options->fields & ~known_fields) {
            alea_set_error_detail(ALEA_ERR_UNSUPPORTED,
                                  "batch options request unsupported fields");
            return -1;
        }
        fields = options->fields;
        max_segments = options->max_segments;
        max_path_entries = options->max_path_entries;
        max_output_bytes = options->max_output_bytes;
    }
    if (batch_validate_ray_inputs(origins_xyz, directions_xyz, ray_count) != 0)
        return -1;
    if (alea_interrupted()) {
        alea_set_error_detail(ALEA_ERR_INTERRUPTED, "batch raycast interrupted");
        return -1;
    }
    if (alea_system_prepare_query_caches(sys,
            ALEA_CACHE_HIER_SPATIAL | ALEA_CACHE_CELL_SURFACES) != 0)
        return -1;
    atomic_init(&live_segment_count, 0);
    atomic_init(&live_path_entry_count, 0);
    if ((fields & (ALEA_RAY_BATCH_PROJECTED_OWNER |
                   ALEA_RAY_BATCH_FULL_PATHS)) == 0) {
        return raycast_hier_batch_execute_common_arena(
            sys, origins_xyz, directions_xyz, ray_count, scalar_t_max,
            t_mins, t_maxs, fields, max_segments, max_output_bytes,
            &live_segment_count, result);
    }
    if (ray_count > SIZE_MAX / sizeof(*traces)) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW, "batch temporary storage overflows");
        return -1;
    }
    traces = calloc(ray_count ? ray_count : 1, sizeof(*traces));
    if (!traces) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate batch ray temporaries");
        return -1;
    }

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (size_t i = 0; i < ray_count; i++) {
        if (alea_interrupted()) {
            traces[i].status = -1;
            continue;
        }
        alea_raycast_result_init(&traces[i].trace);
        if (fields & (ALEA_RAY_BATCH_PROJECTED_OWNER | ALEA_RAY_BATCH_FULL_PATHS))
            traces[i].trace.capture_paths = 1;
        if (max_segments != 0) {
            traces[i].trace.segment_counter = &live_segment_count;
            traces[i].trace.segment_limit = max_segments;
        }
        if ((fields & ALEA_RAY_BATCH_FULL_PATHS) && max_path_entries != 0) {
            traces[i].trace.path_entry_counter = &live_path_entry_count;
            traces[i].trace.path_entry_limit = max_path_entries;
        }
        const double* o = &origins_xyz[i * 3];
        const double* d = &directions_xyz[i * 3];
        alea_ray_t ray;
        if (alea_ray_init(&ray, o[0], o[1], o[2], d[0], d[1], d[2]) != 0) {
            traces[i].status = -1;
            continue;
        }
        const double ray_t_max = t_maxs ? t_maxs[i] : scalar_t_max;
        traces[i].status = alea_raycast_hier_segments_nocache(
            sys, &ray, ray_t_max, &traces[i].trace);
        if (traces[i].status == 0 && t_mins && t_mins[i] > 0.0) {
            const double ray_t_min = t_mins[i];
            size_t write = 0;
            for (size_t j = 0; j < traces[i].trace.segments.count; j++) {
                alea_ray_segment_t segment = traces[i].trace.segments.data[j];
                if (segment.t_exit <= ray_t_min + RAY_EPSILON) continue;
                if (segment.t_enter < ray_t_min) {
                    segment.t_enter = ray_t_min;
                    segment.enter_surface_id = -1;
                    segment.enter_hit_index = -1;
                }
                traces[i].trace.segments.data[write++] = segment;
            }
            traces[i].trace.segments.count = write;
        }
    }
    for (size_t i = 0; i < ray_count; i++) {
        if (traces[i].status != 0 || alea_interrupted()) {
            if (alea_interrupted())
                alea_set_error_detail(ALEA_ERR_INTERRUPTED, "batch raycast interrupted");
            else if (traces[i].trace.segment_limit_exceeded)
                alea_set_error_detail(ALEA_ERR_OVERFLOW,
                                      "batch segment limit (%llu) exceeded",
                                      (unsigned long long)max_segments);
            else if (traces[i].trace.path_entry_limit_exceeded)
                alea_set_error_detail(ALEA_ERR_OVERFLOW,
                                      "batch path-entry limit (%llu) exceeded",
                                      (unsigned long long)max_path_entries);
            else
                alea_set_error_detail(ALEA_ERR_INVALID_STATE, "batch raycast failed");
            goto cleanup;
        }
    }

    for (size_t i = 0; i < ray_count; i++) {
        batch_work_stats_accumulate(&next.work_stats, &traces[i].trace);
        batch_work_stats_add_staging(&next.work_stats, &traces[i].trace);
    }

    next.ray_count = ray_count;
    next.fields = fields;
    next.ray_offsets = batch_alloc_array(ray_count + 1, sizeof(*next.ray_offsets));
    if (!next.ray_offsets) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "failed to allocate batch offsets");
        goto cleanup;
    }
    next.ray_offsets[0] = 0;
    for (size_t i = 0; i < ray_count; i++) {
        size_t count = traces[i].trace.segments.count;
        if (count > UINT64_MAX - next.ray_offsets[i]) {
            alea_set_error_detail(ALEA_ERR_OVERFLOW, "batch segment offsets overflow");
            goto cleanup;
        }
        next.ray_offsets[i + 1] = next.ray_offsets[i] + (uint64_t)count;
    }
    if (next.ray_offsets[ray_count] > SIZE_MAX) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW, "batch segment count overflows size_t");
        goto cleanup;
    }
    next.segment_count = (size_t)next.ray_offsets[ray_count];
    if (max_segments != 0 && next.ray_offsets[ray_count] > max_segments) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "batch segment limit (%llu) exceeded",
                              (unsigned long long)max_segments);
        goto cleanup;
    }
    if (fields & ALEA_RAY_BATCH_FULL_PATHS) {
        next.segment_path_offsets = batch_alloc_array(next.segment_count + 1,
                                                      sizeof(*next.segment_path_offsets));
        if (!next.segment_path_offsets) {
            alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                                  "failed to allocate batch path offsets");
            goto cleanup;
        }
        next.segment_path_offsets[0] = 0;
        size_t segment_index = 0;
        for (size_t i = 0; i < ray_count; i++) {
            const alea_ray_segment_vec_t* segments = &traces[i].trace.segments;
            for (size_t j = 0; j < segments->count; j++, segment_index++) {
                size_t path_count = 0;
                (void)batch_segment_path(&traces[i].trace, &segments->data[j],
                                         &path_count);
                if (path_count > UINT64_MAX - next.segment_path_offsets[segment_index]) {
                    alea_set_error_detail(ALEA_ERR_OVERFLOW,
                                          "batch path offsets overflow");
                    goto cleanup;
                }
                next.segment_path_offsets[segment_index + 1] =
                    next.segment_path_offsets[segment_index] + path_count;
            }
        }
        if (next.segment_path_offsets[next.segment_count] > SIZE_MAX) {
            alea_set_error_detail(ALEA_ERR_OVERFLOW,
                                  "batch path-entry count overflows size_t");
            goto cleanup;
        }
        next.path_entry_count = (size_t)next.segment_path_offsets[next.segment_count];
        if (max_path_entries != 0 &&
            next.segment_path_offsets[next.segment_count] > max_path_entries) {
            alea_set_error_detail(ALEA_ERR_OVERFLOW,
                                  "batch path-entry limit (%llu) exceeded",
                                  (unsigned long long)max_path_entries);
            goto cleanup;
        }
    }

    if (max_output_bytes != 0) {
        size_t output_bytes = 0;
        int overflow =
            batch_add_output_bytes(&output_bytes, ray_count + 1,
                                   sizeof(*next.ray_offsets)) ||
            batch_add_output_bytes(&output_bytes, next.segment_count,
                                   sizeof(*next.t_enter)) ||
            batch_add_output_bytes(&output_bytes, next.segment_count,
                                   sizeof(*next.t_exit)) ||
            batch_add_output_bytes(&output_bytes, next.segment_count,
                                   sizeof(*next.cell_ids));
        if (fields & ALEA_RAY_BATCH_MATERIAL)
            overflow |= batch_add_output_bytes(&output_bytes, next.segment_count,
                                               sizeof(*next.material_ids));
        if (fields & ALEA_RAY_BATCH_DENSITY)
            overflow |= batch_add_output_bytes(&output_bytes, next.segment_count,
                                               sizeof(*next.densities));
        if (fields & ALEA_RAY_BATCH_SURFACES) {
            overflow |= batch_add_output_bytes(&output_bytes, next.segment_count,
                                               sizeof(*next.enter_surface_ids));
            overflow |= batch_add_output_bytes(&output_bytes, next.segment_count,
                                               sizeof(*next.exit_surface_ids));
        }
        if (fields & ALEA_RAY_BATCH_RESOLUTION_FLAGS)
            overflow |= batch_add_output_bytes(&output_bytes, next.segment_count,
                                               sizeof(*next.resolution_flags));
        if (fields & ALEA_RAY_BATCH_PROJECTED_OWNER) {
            overflow |= batch_add_output_bytes(&output_bytes, next.segment_count,
                                               sizeof(*next.projected_cell_ids));
            overflow |= batch_add_output_bytes(&output_bytes, next.segment_count,
                                               sizeof(*next.projected_material_ids));
            overflow |= batch_add_output_bytes(&output_bytes, next.segment_count,
                                               sizeof(*next.projected_universe_ids));
            overflow |= batch_add_output_bytes(&output_bytes, next.segment_count,
                                               sizeof(*next.projected_fill_universes));
            overflow |= batch_add_output_bytes(&output_bytes, next.segment_count,
                                               sizeof(*next.projected_depths));
            overflow |= batch_add_output_bytes(&output_bytes, next.segment_count,
                                               sizeof(*next.projected_is_lattice));
            overflow |= batch_add_output_bytes(&output_bytes, next.segment_count,
                                               sizeof(*next.projected_occurrence_keys));
        }
        if (fields & ALEA_RAY_BATCH_FULL_PATHS) {
            overflow |= batch_add_output_bytes(&output_bytes, next.segment_count + 1,
                                               sizeof(*next.segment_path_offsets));
            overflow |= batch_add_output_bytes(&output_bytes, next.path_entry_count,
                                               sizeof(*next.path_cell_ids));
            overflow |= batch_add_output_bytes(&output_bytes, next.path_entry_count,
                                               sizeof(*next.path_material_ids));
            overflow |= batch_add_output_bytes(&output_bytes, next.path_entry_count,
                                               sizeof(*next.path_universe_ids));
            overflow |= batch_add_output_bytes(&output_bytes, next.path_entry_count,
                                               sizeof(*next.path_fill_universes));
            overflow |= batch_add_output_bytes(&output_bytes, next.path_entry_count,
                                               sizeof(*next.path_depths));
            overflow |= batch_add_output_bytes(&output_bytes, next.path_entry_count,
                                               sizeof(*next.path_is_lattice));
            if (next.path_entry_count > SIZE_MAX / 3 ||
                batch_add_output_bytes(&output_bytes, next.path_entry_count * 3,
                                       sizeof(*next.path_lattice_origins_xyz))) {
                overflow = 1;
            }
            overflow |= batch_add_output_bytes(&output_bytes, next.path_entry_count,
                                               sizeof(*next.path_occurrence_keys));
        }
        if (overflow || output_bytes > max_output_bytes) {
            alea_set_error_detail(ALEA_ERR_OVERFLOW,
                                  "batch compact output exceeds byte limit (%llu)",
                                  (unsigned long long)max_output_bytes);
            goto cleanup;
        }
    }

    next.t_enter = batch_alloc_array(next.segment_count, sizeof(*next.t_enter));
    next.t_exit = batch_alloc_array(next.segment_count, sizeof(*next.t_exit));
    next.cell_ids = batch_alloc_array(next.segment_count, sizeof(*next.cell_ids));
    if (fields & ALEA_RAY_BATCH_MATERIAL)
        next.material_ids = batch_alloc_array(next.segment_count, sizeof(*next.material_ids));
    if (fields & ALEA_RAY_BATCH_DENSITY)
        next.densities = batch_alloc_array(next.segment_count, sizeof(*next.densities));
    if (fields & ALEA_RAY_BATCH_SURFACES) {
        next.enter_surface_ids = batch_alloc_array(next.segment_count, sizeof(*next.enter_surface_ids));
        next.exit_surface_ids = batch_alloc_array(next.segment_count, sizeof(*next.exit_surface_ids));
    }
    if (fields & ALEA_RAY_BATCH_RESOLUTION_FLAGS)
        next.resolution_flags = batch_alloc_array(next.segment_count, sizeof(*next.resolution_flags));
    if (fields & ALEA_RAY_BATCH_PROJECTED_OWNER) {
        next.projected_cell_ids = batch_alloc_array(next.segment_count,
                                                    sizeof(*next.projected_cell_ids));
        next.projected_material_ids = batch_alloc_array(next.segment_count,
                                                        sizeof(*next.projected_material_ids));
        next.projected_universe_ids = batch_alloc_array(next.segment_count,
                                                        sizeof(*next.projected_universe_ids));
        next.projected_fill_universes = batch_alloc_array(next.segment_count,
                                                          sizeof(*next.projected_fill_universes));
        next.projected_depths = batch_alloc_array(next.segment_count,
                                                  sizeof(*next.projected_depths));
        next.projected_is_lattice = batch_alloc_array(next.segment_count,
                                                      sizeof(*next.projected_is_lattice));
        next.projected_occurrence_keys = batch_alloc_array(next.segment_count,
                                                           sizeof(*next.projected_occurrence_keys));
    }
    if (fields & ALEA_RAY_BATCH_FULL_PATHS) {
        next.path_cell_ids = batch_alloc_array(next.path_entry_count, sizeof(*next.path_cell_ids));
        next.path_material_ids = batch_alloc_array(next.path_entry_count, sizeof(*next.path_material_ids));
        next.path_universe_ids = batch_alloc_array(next.path_entry_count, sizeof(*next.path_universe_ids));
        next.path_fill_universes = batch_alloc_array(next.path_entry_count, sizeof(*next.path_fill_universes));
        next.path_depths = batch_alloc_array(next.path_entry_count, sizeof(*next.path_depths));
        next.path_is_lattice = batch_alloc_array(next.path_entry_count, sizeof(*next.path_is_lattice));
        if (next.path_entry_count > SIZE_MAX / 3) {
            alea_set_error_detail(ALEA_ERR_OVERFLOW, "batch path origins overflow");
            goto cleanup;
        }
        next.path_lattice_origins_xyz = batch_alloc_array(next.path_entry_count * 3,
                                                           sizeof(*next.path_lattice_origins_xyz));
        next.path_occurrence_keys = batch_alloc_array(next.path_entry_count,
                                                       sizeof(*next.path_occurrence_keys));
    }
    if (!next.t_enter || !next.t_exit || !next.cell_ids ||
        ((fields & ALEA_RAY_BATCH_MATERIAL) && !next.material_ids) ||
        ((fields & ALEA_RAY_BATCH_DENSITY) && !next.densities) ||
        ((fields & ALEA_RAY_BATCH_SURFACES) &&
         (!next.enter_surface_ids || !next.exit_surface_ids)) ||
        ((fields & ALEA_RAY_BATCH_RESOLUTION_FLAGS) && !next.resolution_flags) ||
        ((fields & ALEA_RAY_BATCH_PROJECTED_OWNER) &&
         (!next.projected_cell_ids || !next.projected_material_ids ||
          !next.projected_universe_ids || !next.projected_fill_universes ||
          !next.projected_depths || !next.projected_is_lattice ||
          !next.projected_occurrence_keys)) ||
        ((fields & ALEA_RAY_BATCH_FULL_PATHS) &&
         (!next.path_cell_ids || !next.path_material_ids || !next.path_universe_ids ||
          !next.path_fill_universes || !next.path_depths || !next.path_is_lattice ||
          !next.path_lattice_origins_xyz || !next.path_occurrence_keys))) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate batch segment fields");
        goto cleanup;
    }

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (size_t i = 0; i < ray_count; i++) {
        size_t dst = (size_t)next.ray_offsets[i];
        const alea_ray_segment_vec_t* segments = &traces[i].trace.segments;
        for (size_t j = 0; j < segments->count; j++) {
            const alea_ray_segment_t* seg = &segments->data[j];
            size_t k = dst + j;
            next.t_enter[k] = seg->t_enter;
            next.t_exit[k] = seg->t_exit;
            next.cell_ids[k] = (int32_t)seg->cell_id;
            if (next.material_ids) next.material_ids[k] = (int32_t)seg->material_id;
            if (next.densities) next.densities[k] = seg->density;
            if (next.enter_surface_ids) next.enter_surface_ids[k] = (int32_t)seg->enter_surface_id;
            if (next.exit_surface_ids) next.exit_surface_ids[k] = (int32_t)seg->exit_surface_id;
            if (next.resolution_flags) next.resolution_flags[k] = seg->resolution_flags;
            if (next.projected_cell_ids)
                batch_store_projected_owner(&next, k, &traces[i].trace, seg,
                                            options ? options->projected_depth : -1);
            if (next.segment_path_offsets) {
                size_t path_count = 0;
                const alea_ray_path_entry_t* path =
                    batch_segment_path(&traces[i].trace, seg, &path_count);
                size_t path_dst = (size_t)next.segment_path_offsets[k];
                for (size_t p = 0; p < path_count; p++) {
                    const alea_ray_path_entry_t* entry = &path[p];
                    size_t q = path_dst + p;
                    next.path_cell_ids[q] = entry->cell_id;
                    next.path_material_ids[q] = entry->material_id;
                    next.path_universe_ids[q] = entry->universe_id;
                    next.path_fill_universes[q] = entry->fill_universe;
                    next.path_depths[q] = entry->depth;
                    next.path_is_lattice[q] = entry->is_lattice;
                    next.path_lattice_origins_xyz[q * 3] = entry->lattice_origin[0];
                    next.path_lattice_origins_xyz[q * 3 + 1] = entry->lattice_origin[1];
                    next.path_lattice_origins_xyz[q * 3 + 2] = entry->lattice_origin[2];
                    next.path_occurrence_keys[q] = entry->occurrence_key;
                }
            }
        }
    }
    for (size_t i = 0; i < ray_count; i++) alea_raycast_result_free(&traces[i].trace);
    free(traces);
    batch_result_free_buffers(result);
    *result = next;
    return 0;

cleanup:
    for (size_t i = 0; i < ray_count; i++) alea_raycast_result_free(&traces[i].trace);
    free(traces);
    batch_result_free_buffers(&next);
    return -1;
}

int alea_raycast_hier_batch(alea_system_t* sys,
                            const double* origins_xyz,
                            const double* directions_xyz,
                            size_t ray_count,
                            double t_max,
                            const alea_raycast_batch_options_t* options,
                            alea_raycast_batch_result_t* result) {
    return raycast_hier_batch_execute(sys, origins_xyz, directions_xyz,
                                      ray_count, t_max, NULL, NULL,
                                      options, result);
}

int alea_raycast_hier_batch_query_nocache(
    alea_system_t* sys, const double* origins_xyz, const double* directions_xyz,
    size_t ray_count, const alea_ray_batch_query_t* query,
    const alea_raycast_batch_options_t* options,
    alea_raycast_batch_result_t* result) {
    if (!query) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG, "batch query descriptor is required");
        return -1;
    }
    if (query->kind != ALEA_RAY_QUERY_SEGMENTS) {
        alea_set_error_detail(ALEA_ERR_UNSUPPORTED,
                              "compact batch query kind is not implemented");
        return -1;
    }
    /* A NULL t_max array is an unbounded range for every ray. */
    return raycast_hier_batch_execute(sys, origins_xyz, directions_xyz,
                                      ray_count, 0.0, query->t_mins,
                                      query->t_maxs, options, result);
}

void alea_ray_first_visible_batch_result_init(
    alea_ray_first_visible_batch_result_t* result) {
    if (result) memset(result, 0, sizeof(*result));
}

void alea_ray_first_visible_batch_result_free(
    alea_ray_first_visible_batch_result_t* result) {
    if (!result) return;
    free(result->found);
    free(result->t);
    free(result->cell_ids);
    free(result->material_ids);
    free(result->densities);
    free(result->surface_ids);
    free(result->primitive_ids);
    free(result->resolution_flags);
    free(result->normals_xyz);
    memset(result, 0, sizeof(*result));
}

static int first_visible_batch_add_bytes(size_t* total, size_t count,
                                         size_t element_size) {
    return batch_add_output_bytes(total, count, element_size);
}

/* Fixed-output queries publish directly into preallocated SoA arrays, so the
 * executor boundary only needs to own one operation-wide region and collect
 * per-ray status.  Keep that policy here, rather than letting each semantic
 * query grow its own OpenMP wrapper.  The worker remains a packet traversal
 * body and is also used unchanged by serial builds. */
typedef int (*raycast_fixed_batch_worker_t)(
    alea_system_t* sys, const double* origins_xyz, const double* directions_xyz,
    size_t ray_count, const alea_ray_batch_query_t* query, void* result,
    int* statuses);

static int raycast_fixed_batch_execute(
    raycast_fixed_batch_worker_t worker, alea_system_t* sys,
    const double* origins_xyz, const double* directions_xyz, size_t ray_count,
    const alea_ray_batch_query_t* query, void* result, int* statuses) {
    int execute_rc = 0;
#ifdef _OPENMP
    #pragma omp parallel shared(execute_rc)
    {
        const int local_rc = worker(sys, origins_xyz, directions_xyz, ray_count,
                                    query, result, statuses);
        if (local_rc != 0) {
            #pragma omp atomic write
            execute_rc = local_rc;
        }
    }
#else
    execute_rc = worker(sys, origins_xyz, directions_xyz, ray_count, query,
                        result, statuses);
#endif
    return execute_rc;
}

static int raycast_first_visible_batch_worker(
    alea_system_t* sys, const double* origins_xyz, const double* directions_xyz,
    size_t ray_count, const alea_ray_batch_query_t* query, void* result,
    int* statuses) {
    return alea_raycast_hier_first_visible_batch_execute_nocache(
        sys, origins_xyz, directions_xyz, ray_count, query, result, statuses);
}

static int raycast_any_hit_batch_worker(
    alea_system_t* sys, const double* origins_xyz, const double* directions_xyz,
    size_t ray_count, const alea_ray_batch_query_t* query, void* result,
    int* statuses) {
    return alea_raycast_hier_any_hit_batch_execute_nocache(
        sys, origins_xyz, directions_xyz, ray_count, query, result, statuses);
}

int alea_raycast_hier_first_visible_batch_nocache(
    alea_system_t* sys, const double* origins_xyz, const double* directions_xyz,
    size_t ray_count, const alea_ray_batch_query_t* query,
    alea_ray_first_visible_batch_result_t* result) {
    const uint32_t known_fields = ALEA_RAY_QUERY_FIELD_CELL_ID |
                                  ALEA_RAY_QUERY_FIELD_MATERIAL_ID |
                                  ALEA_RAY_QUERY_FIELD_DENSITY |
                                  ALEA_RAY_QUERY_FIELD_SURFACE_ID |
                                  ALEA_RAY_QUERY_FIELD_SURFACE_NORMAL |
                                  ALEA_RAY_QUERY_FIELD_RESOLUTION_FLAGS |
                                  ALEA_RAY_QUERY_FIELD_PRIMITIVE_ID;
    alea_ray_first_visible_batch_result_t next;
    int* statuses = NULL;
    size_t output_bytes = 0;
    memset(&next, 0, sizeof(next));
    if (!sys || !query || !result) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "first-visible batch requires system, query, and result");
        return -1;
    }
    if (query->kind != ALEA_RAY_QUERY_FIRST_VISIBLE ||
        (query->fields & ~known_fields)) {
        alea_set_error_detail(ALEA_ERR_UNSUPPORTED,
                              "unsupported first-visible batch query descriptor");
        return -1;
    }
    if (batch_validate_ray_inputs(origins_xyz, directions_xyz, ray_count) != 0)
        return -1;
    for (size_t i = 0; i < ray_count; i++) {
        const double t_min = query->t_mins ? query->t_mins[i] : 0.0;
        const double t_max = query->t_maxs ? query->t_maxs[i] : 0.0;
        if (!isfinite(t_min) || !isfinite(t_max) || t_min < 0.0 ||
            (t_max > 0.0 && t_min > t_max)) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "batch ray %zu has an invalid t range", i);
            return -1;
        }
    }
    if (alea_interrupted()) {
        alea_set_error_detail(ALEA_ERR_INTERRUPTED, "batch raycast interrupted");
        return -1;
    }
    if ((sys->has_lattice ? alea_raycast_ensure_caches(sys) :
                             alea_raycast_ensure_hier_caches(sys)) != 0)
        return -1;

    int overflow =
        first_visible_batch_add_bytes(&output_bytes, ray_count, sizeof(*next.found)) ||
        first_visible_batch_add_bytes(&output_bytes, ray_count, sizeof(*next.t)) ||
        first_visible_batch_add_bytes(&output_bytes, ray_count, sizeof(*next.cell_ids)) ||
        first_visible_batch_add_bytes(&output_bytes, ray_count, sizeof(*next.material_ids));
    if (query->fields & ALEA_RAY_QUERY_FIELD_DENSITY)
        overflow |= first_visible_batch_add_bytes(&output_bytes, ray_count,
                                                  sizeof(*next.densities));
    if (query->fields & (ALEA_RAY_QUERY_FIELD_SURFACE_ID |
                         ALEA_RAY_QUERY_FIELD_PRIMITIVE_ID)) {
        overflow |= first_visible_batch_add_bytes(&output_bytes, ray_count,
                                                  sizeof(*next.surface_ids));
        overflow |= first_visible_batch_add_bytes(&output_bytes, ray_count,
                                                  sizeof(*next.primitive_ids));
    }
    if (query->fields & ALEA_RAY_QUERY_FIELD_RESOLUTION_FLAGS)
        overflow |= first_visible_batch_add_bytes(&output_bytes, ray_count,
                                                  sizeof(*next.resolution_flags));
    if (query->fields & ALEA_RAY_QUERY_FIELD_SURFACE_NORMAL) {
        if (ray_count > SIZE_MAX / 3 ||
            first_visible_batch_add_bytes(&output_bytes, ray_count * 3,
                                          sizeof(*next.normals_xyz)))
            overflow = 1;
    }
    if (overflow || (query->max_output_bytes &&
                     output_bytes > query->max_output_bytes)) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "first-visible batch output exceeds byte limit");
        return -1;
    }

    next.ray_count = ray_count;
    next.found = batch_alloc_array(ray_count, sizeof(*next.found));
    next.t = batch_alloc_array(ray_count, sizeof(*next.t));
    next.cell_ids = batch_alloc_array(ray_count, sizeof(*next.cell_ids));
    next.material_ids = batch_alloc_array(ray_count, sizeof(*next.material_ids));
    if (query->fields & ALEA_RAY_QUERY_FIELD_DENSITY)
        next.densities = batch_alloc_array(ray_count, sizeof(*next.densities));
    if (query->fields & (ALEA_RAY_QUERY_FIELD_SURFACE_ID |
                         ALEA_RAY_QUERY_FIELD_PRIMITIVE_ID)) {
        next.surface_ids = batch_alloc_array(ray_count, sizeof(*next.surface_ids));
        next.primitive_ids = batch_alloc_array(ray_count, sizeof(*next.primitive_ids));
    }
    if (query->fields & ALEA_RAY_QUERY_FIELD_RESOLUTION_FLAGS)
        next.resolution_flags = batch_alloc_array(ray_count, sizeof(*next.resolution_flags));
    if (query->fields & ALEA_RAY_QUERY_FIELD_SURFACE_NORMAL)
        next.normals_xyz = batch_alloc_array(ray_count * 3, sizeof(*next.normals_xyz));
    if (!next.found || !next.t || !next.cell_ids || !next.material_ids ||
        ((query->fields & ALEA_RAY_QUERY_FIELD_DENSITY) && !next.densities) ||
        ((query->fields & (ALEA_RAY_QUERY_FIELD_SURFACE_ID |
                           ALEA_RAY_QUERY_FIELD_PRIMITIVE_ID)) &&
         (!next.surface_ids || !next.primitive_ids)) ||
        ((query->fields & ALEA_RAY_QUERY_FIELD_RESOLUTION_FLAGS) &&
         !next.resolution_flags) ||
        ((query->fields & ALEA_RAY_QUERY_FIELD_SURFACE_NORMAL) &&
         !next.normals_xyz)) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate first-visible batch output");
        goto fail;
    }
    statuses = calloc(ray_count ? ray_count : 1, sizeof(*statuses));
    if (!statuses) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate first-visible batch statuses");
        goto fail;
    }

    const int execute_rc = raycast_fixed_batch_execute(
        raycast_first_visible_batch_worker, sys, origins_xyz, directions_xyz,
        ray_count, query, &next, statuses);
    if (execute_rc != 0 && !alea_interrupted()) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              "native first-visible packet traversal failed");
        goto fail;
    }
    for (size_t i = 0; i < ray_count; i++) {
        if (statuses[i] != 0 || alea_interrupted()) {
            alea_set_error_detail(alea_interrupted() ? ALEA_ERR_INTERRUPTED :
                                  ALEA_ERR_INVALID_STATE,
                                  "first-visible batch raycast failed");
            goto fail;
        }
    }
    free(statuses);
    alea_ray_first_visible_batch_result_free(result);
    *result = next;
    return 0;

fail:
    free(statuses);
    alea_ray_first_visible_batch_result_free(&next);
    return -1;
}

void alea_ray_any_hit_batch_result_init(alea_ray_any_hit_batch_result_t* result) {
    if (result) memset(result, 0, sizeof(*result));
}

void alea_ray_any_hit_batch_result_free(alea_ray_any_hit_batch_result_t* result) {
    if (!result) return;
    free(result->hits);
    memset(result, 0, sizeof(*result));
}

int alea_raycast_hier_any_hit_batch_nocache(
    alea_system_t* sys, const double* origins_xyz, const double* directions_xyz,
    size_t ray_count, const alea_ray_batch_query_t* query,
    alea_ray_any_hit_batch_result_t* result) {
    alea_ray_any_hit_batch_result_t next = {0};
    int* statuses = NULL;
    if (!sys || !query || !result) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "any-hit batch requires system, query, and result");
        return -1;
    }
    if (query->kind != ALEA_RAY_QUERY_ANY_HIT || query->fields != 0) {
        alea_set_error_detail(ALEA_ERR_UNSUPPORTED,
                              "any-hit batch does not materialize fields");
        return -1;
    }
    if (batch_validate_ray_inputs(origins_xyz, directions_xyz, ray_count) != 0)
        return -1;
    for (size_t i = 0; i < ray_count; i++) {
        const double t_min = query->t_mins ? query->t_mins[i] : 0.0;
        const double t_max = query->t_maxs ? query->t_maxs[i] : 0.0;
        if (!isfinite(t_min) || !isfinite(t_max) || t_min < 0.0 ||
            (t_max > 0.0 && t_min > t_max)) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "batch ray %zu has an invalid t range", i);
            return -1;
        }
    }
    if (query->max_output_bytes && ray_count > query->max_output_bytes) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "any-hit batch output exceeds byte limit");
        return -1;
    }
    if (alea_interrupted()) {
        alea_set_error_detail(ALEA_ERR_INTERRUPTED, "batch raycast interrupted");
        return -1;
    }
    if ((sys->has_lattice ? alea_raycast_ensure_caches(sys) :
                             alea_raycast_ensure_hier_caches(sys)) != 0)
        return -1;
    next.ray_count = ray_count;
    next.hits = batch_alloc_array(ray_count, sizeof(*next.hits));
    statuses = calloc(ray_count ? ray_count : 1, sizeof(*statuses));
    if (!next.hits || !statuses) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate any-hit batch output");
        goto fail;
    }

    const int execute_rc = raycast_fixed_batch_execute(
        raycast_any_hit_batch_worker, sys, origins_xyz, directions_xyz,
        ray_count, query, &next, statuses);
    if (execute_rc != 0 && !alea_interrupted()) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              "native any-hit packet traversal failed");
        goto fail;
    }
    for (size_t i = 0; i < ray_count; i++) {
        if (statuses[i] != 0 || alea_interrupted()) {
            alea_set_error_detail(alea_interrupted() ? ALEA_ERR_INTERRUPTED :
                                  ALEA_ERR_INVALID_STATE,
                                  "any-hit batch raycast failed");
            goto fail;
        }
    }
    free(statuses);
    alea_ray_any_hit_batch_result_free(result);
    *result = next;
    return 0;

fail:
    free(statuses);
    alea_ray_any_hit_batch_result_free(&next);
    return -1;
}

/* Variable-output boundary batches keep one reusable scalar trace/event pair
 * per worker.  The worker arena retains only compact events and local CSR
 * offsets; it never keeps one rich trace per input ray. */
typedef struct {
    alea_raycast_result_t trace;
    alea_ray_boundary_event_result_t events;
    alea_ray_boundary_event_vec_t arena;
    uint64_t* row_offsets;
    size_t row_count;
    uint64_t breakpoint_hits;
    uint64_t selected_segments;
} alea_batch_event_worker_t;

static void batch_event_worker_free(alea_batch_event_worker_t* worker) {
    if (!worker) return;
    alea_raycast_result_free(&worker->trace);
    alea_ray_boundary_event_result_free(&worker->events);
    alea_vec_free(&worker->arena);
    free(worker->row_offsets);
    memset(worker, 0, sizeof(*worker));
}

static int batch_event_worker_build(
    alea_system_t* sys, const double* origins_xyz, const double* directions_xyz,
    size_t ray_count, const alea_ray_batch_query_t* query,
    alea_batch_event_worker_t* worker, size_t worker_index, size_t worker_count,
    atomic_uint_fast64_t* live_event_count) {
    const size_t owned_rows = ray_count > worker_index
        ? 1 + (ray_count - 1 - worker_index) / worker_count : 0;
    worker->row_offsets = calloc(owned_rows + 1, sizeof(*worker->row_offsets));
    if (!worker->row_offsets) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate boundary-event worker offsets");
        return -1;
    }
    worker->row_count = owned_rows;
    alea_raycast_result_init(&worker->trace);
    alea_ray_boundary_event_result_init(&worker->events);
    alea_vec_init(&worker->arena);
    size_t local_row = 0;
    for (size_t row = worker_index; row < ray_count;
         row += worker_count, local_row++) {
        if (alea_interrupted()) return -1;
        alea_raycast_result_clear(&worker->trace);
        alea_ray_boundary_event_result_clear(&worker->events);
        alea_ray_t ray;
        const double* o = &origins_xyz[row * 3];
        const double* d = &directions_xyz[row * 3];
        const double t_min = query->t_mins ? query->t_mins[row] : 0.0;
        const double t_max = query->t_maxs ? query->t_maxs[row] : 0.0;
        const alea_ray_boundary_event_options_internal_t event_options = {
            .include_all_coincident_physical =
                query->include_all_coincident_physical,
            .use_hier_blas = query->use_hier_blas
        };
        if (alea_ray_init(&ray, o[0], o[1], o[2], d[0], d[1], d[2]) != 0 ||
            alea_raycast_boundary_events_with_options(
                sys, &ray, t_max, &event_options,
                &worker->trace, &worker->events) != 0)
            return -1;
        worker->breakpoint_hits += worker->events.breakpoint_hits;
        worker->selected_segments += worker->events.selected_segments;
        worker->row_offsets[local_row] = worker->arena.count;
        for (size_t event = 0; event < worker->events.events.count; event++) {
            const alea_ray_boundary_event_t value =
                worker->events.events.data[event];
            if (value.t + RAY_EPSILON < t_min) continue;
            if (query->max_events != 0 &&
                atomic_fetch_add(live_event_count, 1) >= query->max_events) {
                alea_set_error_detail(ALEA_ERR_OVERFLOW,
                                      "boundary-event batch event limit exceeded");
                return -1;
            }
            if (alea_vec_push(&worker->arena, value,
                              alea_ray_boundary_event_t) != 0) {
                alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                                      "failed to append boundary-event worker arena");
                return -1;
            }
        }
        worker->row_offsets[local_row + 1] = worker->arena.count;
    }
    return 0;
}

void alea_ray_boundary_event_batch_result_init(
    alea_ray_boundary_event_batch_result_t* result) {
    if (result) memset(result, 0, sizeof(*result));
}

void alea_ray_boundary_event_batch_result_free(
    alea_ray_boundary_event_batch_result_t* result) {
    if (!result) return;
    free(result->ray_offsets);
    free(result->t);
    free(result->kinds);
    free(result->surface_ids);
    free(result->cell_before);
    free(result->cell_after);
    free(result->material_before);
    free(result->material_after);
    free(result->resolution_flags);
    free(result->primitive_ids);
    free(result->normals_xyz);
    memset(result, 0, sizeof(*result));
}

int alea_raycast_boundary_events_batch_nocache(
    alea_system_t* sys, const double* origins_xyz, const double* directions_xyz,
    size_t ray_count, const alea_ray_batch_query_t* query,
    alea_ray_boundary_event_batch_result_t* result) {
    const uint32_t known_fields = ALEA_RAY_QUERY_FIELD_SURFACE_NORMAL |
                                  ALEA_RAY_QUERY_FIELD_PRIMITIVE_ID;
    alea_ray_boundary_event_batch_result_t next = {0};
    alea_batch_event_worker_t* workers = NULL;
    size_t worker_count = 1;
    atomic_uint_fast64_t live_event_count;
    size_t output_bytes = 0;
    if (!sys || !query || !result) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "boundary-event batch requires system, query, and result");
        return -1;
    }
    if (query->kind != ALEA_RAY_QUERY_BOUNDARY_EVENTS ||
        (query->fields & ~known_fields) || query->material_filter >= 0) {
        alea_set_error_detail(ALEA_ERR_UNSUPPORTED,
                              "unsupported boundary-event batch query descriptor");
        return -1;
    }
    if (batch_validate_ray_inputs(origins_xyz, directions_xyz, ray_count) != 0)
        return -1;
    for (size_t i = 0; i < ray_count; i++) {
        const double t_min = query->t_mins ? query->t_mins[i] : 0.0;
        const double t_max = query->t_maxs ? query->t_maxs[i] : 0.0;
        if (!isfinite(t_min) || !isfinite(t_max) || t_min < 0.0 ||
            (t_max > 0.0 && t_min > t_max)) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "batch ray %zu has an invalid t range", i);
            return -1;
        }
    }
    if (alea_interrupted()) {
        alea_set_error_detail(ALEA_ERR_INTERRUPTED, "batch raycast interrupted");
        return -1;
    }
    if ((query->use_hier_blas
            ? alea_raycast_ensure_hier_caches(sys)
            : alea_raycast_ensure_caches(sys)) != 0)
        return -1;
    atomic_init(&live_event_count, 0);
#ifdef _OPENMP
    worker_count = (size_t)omp_get_max_threads();
#endif
    if (query->max_workers && worker_count > query->max_workers)
        worker_count = query->max_workers;
    if (worker_count > ray_count && ray_count != 0) worker_count = ray_count;
    if (worker_count == 0 || worker_count > SIZE_MAX / sizeof(*workers)) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "boundary-event batch worker storage overflows");
        return -1;
    }
    workers = calloc(worker_count, sizeof(*workers));
    if (!workers) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate boundary-event batch workers");
        return -1;
    }

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (size_t worker = 0; worker < worker_count; worker++) {
        if (batch_event_worker_build(sys, origins_xyz, directions_xyz,
                                     ray_count, query, &workers[worker], worker,
                                     worker_count, &live_event_count) != 0) {
            /* The caller checks every worker after the region so no thread
             * mutates a shared result or error state while traversing. */
            workers[worker].row_count = SIZE_MAX;
        }
    }
    for (size_t worker = 0; worker < worker_count; worker++) {
        if (workers[worker].row_count == SIZE_MAX || alea_interrupted()) {
            alea_set_error_detail(alea_interrupted() ? ALEA_ERR_INTERRUPTED :
                                  ALEA_ERR_INVALID_STATE,
                                  "boundary-event batch raycast failed");
            goto fail;
        }
        next.breakpoint_hits += workers[worker].breakpoint_hits;
        next.selected_segments += workers[worker].selected_segments;
    }

    next.ray_count = ray_count;
    next.ray_offsets = batch_alloc_array(ray_count + 1, sizeof(*next.ray_offsets));
    if (!next.ray_offsets) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate boundary-event batch offsets");
        goto fail;
    }
    next.ray_offsets[0] = 0;
    for (size_t i = 0; i < ray_count; i++) {
        const size_t worker_index = i % worker_count;
        const size_t local_row = i / worker_count;
        const alea_batch_event_worker_t* worker = &workers[worker_index];
        const size_t count = (size_t)(worker->row_offsets[local_row + 1] -
                                      worker->row_offsets[local_row]);
        if (count > UINT64_MAX - next.ray_offsets[i]) {
            alea_set_error_detail(ALEA_ERR_OVERFLOW,
                                  "boundary-event batch offsets overflow");
            goto fail;
        }
        next.ray_offsets[i + 1] = next.ray_offsets[i] + count;
    }
    if (next.ray_offsets[ray_count] > SIZE_MAX) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "boundary-event batch count overflows size_t");
        goto fail;
    }
    next.event_count = (size_t)next.ray_offsets[ray_count];
    if (query->max_events && next.event_count > query->max_events) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "boundary-event batch event limit exceeded");
        goto fail;
    }
    int overflow =
        batch_add_output_bytes(&output_bytes, ray_count + 1,
                               sizeof(*next.ray_offsets)) ||
        batch_add_output_bytes(&output_bytes, next.event_count, sizeof(*next.t)) ||
        batch_add_output_bytes(&output_bytes, next.event_count, sizeof(*next.kinds)) ||
        batch_add_output_bytes(&output_bytes, next.event_count, sizeof(*next.surface_ids)) ||
        batch_add_output_bytes(&output_bytes, next.event_count, sizeof(*next.cell_before)) ||
        batch_add_output_bytes(&output_bytes, next.event_count, sizeof(*next.cell_after)) ||
        batch_add_output_bytes(&output_bytes, next.event_count, sizeof(*next.material_before)) ||
        batch_add_output_bytes(&output_bytes, next.event_count, sizeof(*next.material_after)) ||
        batch_add_output_bytes(&output_bytes, next.event_count,
                               sizeof(*next.resolution_flags));
    if (query->fields & ALEA_RAY_QUERY_FIELD_PRIMITIVE_ID)
        overflow |= batch_add_output_bytes(&output_bytes, next.event_count,
                                           sizeof(*next.primitive_ids));
    if (query->fields & ALEA_RAY_QUERY_FIELD_SURFACE_NORMAL) {
        if (next.event_count > SIZE_MAX / 3 ||
            batch_add_output_bytes(&output_bytes, next.event_count * 3,
                                   sizeof(*next.normals_xyz)))
            overflow = 1;
    }
    if (overflow || (query->max_output_bytes &&
                     output_bytes > query->max_output_bytes)) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "boundary-event batch output exceeds byte limit");
        goto fail;
    }
    next.t = batch_alloc_array(next.event_count, sizeof(*next.t));
    next.kinds = batch_alloc_array(next.event_count, sizeof(*next.kinds));
    next.surface_ids = batch_alloc_array(next.event_count, sizeof(*next.surface_ids));
    next.cell_before = batch_alloc_array(next.event_count, sizeof(*next.cell_before));
    next.cell_after = batch_alloc_array(next.event_count, sizeof(*next.cell_after));
    next.material_before = batch_alloc_array(next.event_count, sizeof(*next.material_before));
    next.material_after = batch_alloc_array(next.event_count, sizeof(*next.material_after));
    next.resolution_flags = batch_alloc_array(next.event_count,
                                              sizeof(*next.resolution_flags));
    if (query->fields & ALEA_RAY_QUERY_FIELD_PRIMITIVE_ID)
        next.primitive_ids = batch_alloc_array(next.event_count,
                                               sizeof(*next.primitive_ids));
    if (query->fields & ALEA_RAY_QUERY_FIELD_SURFACE_NORMAL)
        next.normals_xyz = batch_alloc_array(next.event_count * 3,
                                             sizeof(*next.normals_xyz));
    if (!next.t || !next.kinds || !next.surface_ids || !next.cell_before ||
        !next.cell_after || !next.material_before || !next.material_after ||
        !next.resolution_flags ||
        ((query->fields & ALEA_RAY_QUERY_FIELD_PRIMITIVE_ID) &&
         !next.primitive_ids) ||
        ((query->fields & ALEA_RAY_QUERY_FIELD_SURFACE_NORMAL) &&
         !next.normals_xyz)) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate boundary-event batch output");
        goto fail;
    }
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (size_t i = 0; i < ray_count; i++) {
        const size_t dst = (size_t)next.ray_offsets[i];
        const size_t worker_index = i % worker_count;
        const size_t local_row = i / worker_count;
        const alea_batch_event_worker_t* worker = &workers[worker_index];
        const size_t begin = (size_t)worker->row_offsets[local_row];
        const size_t end = (size_t)worker->row_offsets[local_row + 1];
        for (size_t j = begin; j < end; j++) {
            const alea_ray_boundary_event_t* event =
                &worker->arena.data[j];
            const size_t k = dst + j - begin;
            next.t[k] = event->t;
            next.kinds[k] = (uint8_t)event->kind;
            next.surface_ids[k] = event->surface_id;
            next.cell_before[k] = event->cell_before;
            next.cell_after[k] = event->cell_after;
            next.material_before[k] = event->material_before;
            next.material_after[k] = event->material_after;
            next.resolution_flags[k] = event->resolution_flags;
            if (next.primitive_ids) next.primitive_ids[k] = event->primitive_id;
            if (next.normals_xyz) {
                next.normals_xyz[k * 3] = event->nx;
                next.normals_xyz[k * 3 + 1] = event->ny;
                next.normals_xyz[k * 3 + 2] = event->nz;
            }
        }
    }
    for (size_t worker = 0; worker < worker_count; worker++)
        batch_event_worker_free(&workers[worker]);
    free(workers);
    alea_ray_boundary_event_batch_result_free(result);
    *result = next;
    return 0;

fail:
    for (size_t worker = 0; worker < worker_count; worker++)
        batch_event_worker_free(&workers[worker]);
    free(workers);
    alea_ray_boundary_event_batch_result_free(&next);
    return -1;
}

int alea_trace_ray_slice_compact(
    alea_system_t* sys,
    const alea_slice_view_t* view,
    size_t row_count,
    const alea_raycast_batch_options_t* options,
    alea_raycast_batch_result_t* result) {
    double* origins = NULL;
    double* directions = NULL;
    const alea_slice_plane_t* plane;
    double span;
    int rc;

    if (!sys || !result) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "ray slice requires a system and result");
        return -1;
    }
    if (row_count == 0) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "ray slice row_count must be positive");
        return -1;
    }
    if (batch_validate_slice_view(view) != 0) return -1;
    if (row_count > SIZE_MAX / (3 * sizeof(double))) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "ray slice row storage overflows");
        return -1;
    }

    origins = malloc(row_count * 3 * sizeof(*origins));
    directions = malloc(row_count * 3 * sizeof(*directions));
    if (!origins || !directions) {
        free(origins);
        free(directions);
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate ray slice rows");
        return -1;
    }

    plane = &view->plane;
    span = view->u_max - view->u_min;
    for (size_t row = 0; row < row_count; row++) {
        double v = view->v_min + ((double)row + 0.5) *
                   (view->v_max - view->v_min) / (double)row_count;
        origins[row * 3] = plane->origin[0] + view->u_min * plane->u_axis[0] +
                           v * plane->v_axis[0];
        origins[row * 3 + 1] = plane->origin[1] + view->u_min * plane->u_axis[1] +
                               v * plane->v_axis[1];
        origins[row * 3 + 2] = plane->origin[2] + view->u_min * plane->u_axis[2] +
                               v * plane->v_axis[2];
        directions[row * 3] = plane->u_axis[0];
        directions[row * 3 + 1] = plane->u_axis[1];
        directions[row * 3 + 2] = plane->u_axis[2];
    }

    rc = alea_raycast_hier_batch(sys, origins, directions, row_count, span,
                                 options, result);
    free(origins);
    free(directions);
    if (rc != 0) return rc;

    /* The generic batch records distance from the left viewport origin. The
     * slice contract exposes view-U coordinates, clipped to its viewport. */
    for (size_t i = 0; i < result->segment_count; i++) {
        double u_enter = result->t_enter[i] + view->u_min;
        double u_exit = result->t_exit[i] + view->u_min;
        result->t_enter[i] = fmax(view->u_min, fmin(view->u_max, u_enter));
        result->t_exit[i] = fmax(view->u_min, fmin(view->u_max, u_exit));
    }
    batch_set_fast_slice_cache(result, sys, view, row_count,
                               options ? options->projected_depth : -1);
    return 0;
}


struct alea_ray_boundary_event_query_result {
    alea_raycast_result_t scratch;
    alea_ray_boundary_event_result_t events;
    uint32_t fields;
    int includes_occurrence_provenance;
};

void alea_ray_boundary_event_options_init(
    alea_ray_boundary_event_options_t* options) {
    if (!options) return;
    memset(options, 0, sizeof(*options));
    options->struct_size = sizeof(*options);
}

alea_ray_boundary_event_query_result_t* alea_ray_boundary_event_query_result_create(void) {
    alea_ray_boundary_event_query_result_t* result = calloc(1, sizeof(*result));
    if (!result) return NULL;
    alea_raycast_result_init(&result->scratch);
    alea_ray_boundary_event_result_init(&result->events);
    return result;
}

void alea_ray_boundary_event_query_result_destroy(
    alea_ray_boundary_event_query_result_t* result) {
    if (!result) return;
    alea_raycast_result_free(&result->scratch);
    alea_ray_boundary_event_result_free(&result->events);
    free(result);
}

int alea_ray_boundary_event_query(alea_system_t* sys,
    double ox, double oy, double oz, double dx, double dy, double dz,
    const alea_ray_boundary_event_options_t* input,
    alea_ray_boundary_event_query_result_t* result) {
    alea_ray_boundary_event_options_t defaults, options;
    if (!sys || !result) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "boundary-event query requires system and result");
        return -1;
    }
    alea_ray_boundary_event_options_init(&defaults); options = defaults;
    if (input) {
        if (input->struct_size < sizeof(input->struct_size)) goto invalid;
        size_t bytes = input->struct_size < sizeof(options) ? input->struct_size : sizeof(options);
        memcpy(&options, input, bytes);
    }
    if ((options.fields & ~(ALEA_RAY_BOUNDARY_EVENT_PRIMITIVE_ID |
                            ALEA_RAY_BOUNDARY_EVENT_NORMAL)) ||
        options.t_min < 0 || (options.t_max > 0 && options.t_min > options.t_max)) goto invalid;
    alea_ray_t ray;
    if (alea_ray_init(&ray, ox, oy, oz, dx, dy, dz) != 0) goto invalid;
    if ((options.include_occurrence_provenance
            ? alea_raycast_ensure_hier_caches(sys)
            : alea_system_prepare_query_caches(sys, ALEA_CACHE_RAYCAST)) != 0) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              "boundary-event query failed to prepare raycast caches");
        goto fail;
    }
    const alea_ray_boundary_event_options_internal_t internal_options = {
        .include_all_coincident_physical = options.include_all_coincident_physical != 0,
        .max_events = options.max_events,
        .max_output_bytes = options.max_output_bytes
    };
    int event_rc = options.include_occurrence_provenance
        ? alea_raycast_selected_boundary_events_with_options_nocache(
            sys, &ray, options.t_max > 0.0 ? options.t_max : DBL_MAX,
            &internal_options, &result->scratch, &result->events)
        : alea_raycast_boundary_events_with_options(
            sys, &ray, options.t_max, &internal_options, &result->scratch,
            &result->events);
    if (event_rc != 0) {
        const char* detail = alea_error();
        if (!detail || !detail[0])
            alea_set_error_detail(
                ALEA_ERR_INVALID_STATE,
                "boundary-event query failed to materialize events");
        goto fail;
    }
    size_t write = 0;
    for (size_t i = 0; i < result->events.events.count; i++) {
        alea_ray_boundary_event_t event = result->events.events.data[i];
        if (event.t + RAY_EPSILON < options.t_min) continue;
        result->events.events.data[write++] = event;
    }
    result->events.events.count = write;
    if ((options.max_events && write > options.max_events) ||
        (options.max_output_bytes &&
         write > options.max_output_bytes / sizeof(*result->events.events.data))) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW, "boundary-event query output limit exceeded");
        goto fail;
    }
    result->fields = options.fields;
    result->includes_occurrence_provenance =
        options.include_occurrence_provenance != 0;
    return 0;
invalid:
    alea_set_error_detail(ALEA_ERR_INVALID_ARG, "invalid boundary-event query options");
fail:
    alea_raycast_result_clear(&result->scratch);
    alea_ray_boundary_event_result_clear(&result->events);
    result->fields = 0;
    result->includes_occurrence_provenance = 0;
    return -1;
}

size_t alea_ray_boundary_event_count(const alea_ray_boundary_event_query_result_t* result) {
    return result ? result->events.events.count : 0;
}

int alea_ray_boundary_event_get(const alea_ray_boundary_event_query_result_t* result,
    size_t index, double* t, int* kind, int* surface_id, int* cell_before,
    int* cell_after, int* material_before, int* material_after,
    uint32_t* resolution_flags, uint32_t* primitive_id,
    double* nx, double* ny, double* nz) {
    if (!result || index >= result->events.events.count) return -1;
    const alea_ray_boundary_event_t* event = &result->events.events.data[index];
    if (t) *t = event->t;
    if (kind) *kind = (int)event->kind;
    if (surface_id) *surface_id = event->surface_id;
    if (cell_before) *cell_before = event->cell_before;
    if (cell_after) *cell_after = event->cell_after;
    if (material_before) *material_before = event->material_before;
    if (material_after) *material_after = event->material_after;
    if (resolution_flags) *resolution_flags = event->resolution_flags;
    if (primitive_id) *primitive_id = (result->fields & ALEA_RAY_BOUNDARY_EVENT_PRIMITIVE_ID) ? event->primitive_id : UINT32_MAX;
    if (nx) *nx = (result->fields & ALEA_RAY_BOUNDARY_EVENT_NORMAL) ? event->nx : 0;
    if (ny) *ny = (result->fields & ALEA_RAY_BOUNDARY_EVENT_NORMAL) ? event->ny : 0;
    if (nz) *nz = (result->fields & ALEA_RAY_BOUNDARY_EVENT_NORMAL) ? event->nz : 0;
    return 0;
}

int alea_ray_boundary_event_provenance_get(
    const alea_ray_boundary_event_query_result_t* result, size_t index,
    alea_ray_boundary_event_provenance_t* out) {
    if (!result || !out || !result->includes_occurrence_provenance ||
        index >= result->events.events.count) return -1;
    const alea_ray_boundary_event_t* event = &result->events.events.data[index];
    memset(out, 0, sizeof(*out));
    out->flags = event->provenance_flags;
    out->active_cell_id = event->active_cell_id;
    out->active_universe_id = event->active_universe_id;
    out->active_depth = event->active_depth;
    out->active_occurrence_key = event->active_occurrence_key;
    out->active_parent_occurrence_key = event->active_parent_occurrence_key;
    out->before_occurrence_key = event->before_occurrence_key;
    out->before_parent_occurrence_key = event->before_parent_occurrence_key;
    out->after_occurrence_key = event->after_occurrence_key;
    out->after_parent_occurrence_key = event->after_parent_occurrence_key;
    memcpy(out->local_point, event->local_point, sizeof(out->local_point));
    memcpy(out->local_direction, event->local_direction,
           sizeof(out->local_direction));
    out->local_surface_count = event->local_surface_count;
    out->local_surface_complete = event->local_surface_complete;
    memcpy(out->local_surface_ids, event->local_surface_ids,
           sizeof(out->local_surface_ids));
    return 0;
}

struct alea_ray_first_visible_query_result {
    alea_raycast_result_t scratch;
    alea_ray_first_visible_result_t answer;
    uint32_t fields;
};

void alea_ray_first_visible_options_init(
    alea_ray_first_visible_options_t* options) {
    if (!options) return;
    memset(options, 0, sizeof(*options));
    options->struct_size = sizeof(*options);
    options->material_filter = -1;
}

alea_ray_first_visible_query_result_t* alea_ray_first_visible_query_result_create(void) {
    alea_ray_first_visible_query_result_t* result = calloc(1, sizeof(*result));
    if (!result) return NULL;
    alea_raycast_result_init(&result->scratch);
    result->answer.cell_id = -1;
    result->answer.surface_id = -1;
    result->answer.primitive_id = UINT32_MAX;
    return result;
}

void alea_ray_first_visible_query_result_destroy(
    alea_ray_first_visible_query_result_t* result) {
    if (!result) return;
    alea_raycast_result_free(&result->scratch);
    free(result);
}

int alea_ray_first_visible_query(alea_system_t* sys,
    double ox, double oy, double oz, double dx, double dy, double dz,
    const alea_ray_first_visible_options_t* input,
    alea_ray_first_visible_query_result_t* result) {
    alea_ray_first_visible_options_t defaults, options;
    if (!sys || !result) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "first-visible query requires system and result");
        return -1;
    }
    alea_ray_first_visible_options_init(&defaults);
    options = defaults;
    if (input) {
        if (input->struct_size < sizeof(input->struct_size)) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "first-visible options struct_size is too small");
            goto fail;
        }
        size_t bytes = input->struct_size < sizeof(options) ?
            input->struct_size : sizeof(options);
        memcpy(&options, input, bytes);
        if (bytes < offsetof(alea_ray_first_visible_options_t, material_filter) +
                        sizeof(options.material_filter))
            options.material_filter = -1;
    }
    if (options.fields & ~(ALEA_RAY_FIRST_VISIBLE_SURFACE_ID |
                           ALEA_RAY_FIRST_VISIBLE_SURFACE_NORMAL) ||
        options.t_min < 0 || (options.t_max > 0 && options.t_min > options.t_max)) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "invalid first-visible query options");
        goto fail;
    }
    alea_ray_t ray;
    if (alea_ray_init(&ray, ox, oy, oz, dx, dy, dz) != 0) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG, "first-visible ray direction is zero");
        goto fail;
    }
    if (alea_system_prepare_query_caches(sys, ALEA_CACHE_RAYCAST_HIER) != 0)
        goto fail;
    const alea_ray_query_t query = {
        .kind = ALEA_RAY_QUERY_FIRST_VISIBLE,
        .fields = (options.fields & ALEA_RAY_FIRST_VISIBLE_SURFACE_NORMAL) ?
            ALEA_RAY_QUERY_FIELD_SURFACE_NORMAL : 0,
        .t_min = options.t_min, .t_max = options.t_max,
        .material_filter = options.material_filter
    };
    alea_ray_query_output_t output;
    if (alea_raycast_query_reuse_nocache(sys, &ray, &query, &result->scratch,
                                         NULL, &output) != 0)
        goto fail;
    result->answer = output.first_visible;
    result->fields = options.fields;
    return 0;
fail:
    alea_raycast_result_clear(&result->scratch);
    memset(&result->answer, 0, sizeof(result->answer));
    result->answer.cell_id = -1; result->answer.surface_id = -1;
    result->answer.primitive_id = UINT32_MAX; result->fields = 0;
    return -1;
}

int alea_ray_first_visible_found(const alea_ray_first_visible_query_result_t* result) { return result && result->answer.found; }
double alea_ray_first_visible_t(const alea_ray_first_visible_query_result_t* result) { return result && result->answer.found ? result->answer.t : 0; }
int alea_ray_first_visible_cell_id(const alea_ray_first_visible_query_result_t* result) { return result && result->answer.found ? result->answer.cell_id : -1; }
int alea_ray_first_visible_material_id(const alea_ray_first_visible_query_result_t* result) { return result && result->answer.found ? result->answer.material_id : 0; }
double alea_ray_first_visible_density(const alea_ray_first_visible_query_result_t* result) { return result && result->answer.found ? result->answer.density : 0; }
int alea_ray_first_visible_surface_id(const alea_ray_first_visible_query_result_t* result) { return result && result->answer.found && (result->fields & ALEA_RAY_FIRST_VISIBLE_SURFACE_ID) ? result->answer.surface_id : -1; }
int alea_ray_first_visible_normal(const alea_ray_first_visible_query_result_t* result, double* nx, double* ny, double* nz) {
    if (!result || !result->answer.found || !(result->fields & ALEA_RAY_FIRST_VISIBLE_SURFACE_NORMAL)) return -1;
    if (nx) *nx = result->answer.nx;
    if (ny) *ny = result->answer.ny;
    if (nz) *nz = result->answer.nz;
    return 0;
}

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

int alea_raycast_segment_resolution_flags(const alea_raycast_result_t* result,
                                          size_t index,
                                          uint8_t* out_flags) {
    if (!result || !out_flags || index >= result->segments.count) return -1;
    *out_flags = result->segments.data[index].resolution_flags;
    return 0;
}

size_t alea_raycast_hit_count(const alea_raycast_result_t* result) {
    return result ? result->hits.count : 0;
}

int alea_raycast_hit_get(const alea_raycast_result_t* result,
                         size_t index,
                         double* out_t,
                         int* out_surface_id) {
    if (!result || !out_t || !out_surface_id || index >= result->hits.count) {
        return -1;
    }
    *out_t = result->hits.data[index].t;
    *out_surface_id = result->hits.data[index].surface_id;
    return 0;
}

void alea_raycast_result_set_path_capture(alea_raycast_result_t* result,
                                          int enabled) {
    if (result) result->capture_paths = enabled ? 1 : 0;
}

size_t alea_raycast_segment_path_count(const alea_raycast_result_t* result,
                                       size_t segment_index) {
    if (!result || segment_index >= result->segments.count) return 0;
    uint32_t path_index = result->segments.data[segment_index].path_index;
    if (path_index == UINT32_MAX || path_index >= result->paths.count) return 0;
    return result->paths.data[path_index].count;
}

int alea_raycast_segment_path_get(const alea_raycast_result_t* result,
                                  size_t segment_index,
                                  size_t path_entry_index,
                                  alea_raycast_path_entry_t* out_entry) {
    if (!result || !out_entry || segment_index >= result->segments.count) return -1;
    uint32_t path_index = result->segments.data[segment_index].path_index;
    if (path_index == UINT32_MAX || path_index >= result->paths.count) return -1;
    const alea_ray_path_t* path = &result->paths.data[path_index];
    if (path_entry_index >= path->count) return -1;
    const alea_ray_path_entry_t* src =
        &result->path_entries.data[path->offset + path_entry_index];
    out_entry->cell_index = src->cell_index;
    out_entry->cell_id = src->cell_id;
    out_entry->material_id = src->material_id;
    out_entry->universe_id = src->universe_id;
    out_entry->fill_universe = src->fill_universe;
    out_entry->depth = src->depth;
    out_entry->is_lattice = src->is_lattice;
    out_entry->lattice_origin[0] = src->lattice_origin[0];
    out_entry->lattice_origin[1] = src->lattice_origin[1];
    out_entry->lattice_origin[2] = src->lattice_origin[2];
    out_entry->occurrence_key = src->occurrence_key;
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

int alea_generate_cauchy_crofton_rays(double cx, double cy, double cz,
                                      double radius, uint32_t* rng_state,
                                      size_t ray_count, double* origins_xyz,
                                      double* directions_xyz) {
    if (!rng_state) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "Cauchy-Crofton generator requires an RNG state");
        return -1;
    }
    if (!isfinite(cx) || !isfinite(cy) || !isfinite(cz) ||
        !isfinite(radius) || radius <= 0.0) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "Cauchy-Crofton sphere must be finite with positive radius");
        return -1;
    }
    if (ray_count > SIZE_MAX / 3) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "Cauchy-Crofton ray count overflows packed XYZ buffers");
        return -1;
    }
    if (ray_count != 0 && (!origins_xyz || !directions_xyz)) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "Cauchy-Crofton output buffers are required");
        return -1;
    }
    for (size_t i = 0; i < ray_count; i++) {
        generate_cauchy_crofton_ray(rng_state, cx, cy, cz, radius,
                                    &origins_xyz[i * 3], &origins_xyz[i * 3 + 1],
                                    &origins_xyz[i * 3 + 2],
                                    &directions_xyz[i * 3], &directions_xyz[i * 3 + 1],
                                    &directions_xyz[i * 3 + 2]);
    }
    return 0;
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
        /* Volume paths expose world-to-local.  Bounding the placed local CSG
         * needs its inverse (local-to-world). */
        memcpy(transform.m, paths[i].world_to_local, sizeof(transform.m));
        transform.has_inverse = false;
        if (!alea_matrix_invert(&transform)) continue;
        memcpy(transform.m, transform.inv, sizeof(transform.m));
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

/* ============================================================================
 * Interval defect classification (see alea_raycast.h)
 * ============================================================================ */

int alea_ray_classify_intervals(alea_system_t* sys,
                                double ox, double oy, double oz,
                                double dx, double dy, double dz,
                                double t_max,
                                alea_ray_interval_finding_t* out,
                                size_t max_out) {
    if (!sys || t_max <= 0.0) return -1;

    /* The recursive owner-set query requires the universe index; the global
     * diagnostic breakpoint enumerator requires the raycast caches. */
    if (!sys->universe_index_built && alea_build_universe_index(sys) != 0)
        return -1;
    if (alea_raycast_ensure_caches(sys) != 0)
        return -1;

    alea_ray_t ray;
    if (alea_ray_init(&ray, ox, oy, oz, dx, dy, dz) != 0) return -1;

    alea_raycast_result_t breakpoints;
    alea_raycast_result_init(&breakpoints);
    const int count = alea_ray_coverage_classify_reuse_nocache(
        sys, &ray, t_max, &breakpoints, out, max_out);
    alea_raycast_result_free(&breakpoints);
    return count;
}
