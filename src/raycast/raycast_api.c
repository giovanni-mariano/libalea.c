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
#include "bvh.h"
#include "core/alea_system.h"
#include "primitives/bbox.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <limits.h>
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
};

typedef struct {
    alea_raycast_result_t trace;
    int status;
} alea_batch_trace_tmp_t;

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

int alea_raycast_hier_batch(alea_system_t* sys,
                            const double* origins_xyz,
                            const double* directions_xyz,
                            size_t ray_count,
                            double t_max,
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
    if (!isfinite(t_max)) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG, "batch t_max must be finite");
        return -1;
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
        traces[i].status = alea_raycast_hier_segments_nocache(
            sys, &ray, t_max, &traces[i].trace);
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

/* ============================================================================
 * Interval defect classification (see alea_raycast.h)
 * ============================================================================ */

#include "core/alea_universe.h"

int alea_ray_classify_intervals(alea_system_t* sys,
                                double ox, double oy, double oz,
                                double dx, double dy, double dz,
                                double t_max,
                                alea_ray_interval_finding_t* out,
                                size_t max_out) {
    if (!sys || t_max <= 0.0) return -1;

    /* The recursive owner-set query requires the universe index; the global
     * hit pipeline requires the raycast caches. */
    if (!sys->universe_index_built && alea_build_universe_index(sys) != 0)
        return -1;
    if (alea_raycast_ensure_caches(sys) != 0)
        return -1;

    /* Global pipeline: every crossing along the ray (physical surfaces,
     * fill-transformed universe content, synthetic lattice boundaries).
     * The precedence segments it also builds are ignored — only the
     * breakpoint list matters here. */
    alea_raycast_result_t r;
    alea_raycast_result_init(&r);
    if (alea_raycast(sys, ox, oy, oz, dx, dy, dz, t_max, &r) != 0) {
        alea_raycast_result_free(&r);
        return -1;
    }

    alea_ray_t ray;
    alea_ray_init(&ray, ox, oy, oz, dx, dy, dz);

    enum { CLASSIFY_MAX_OWNERS = 32 };
    alea_cell_hit_t hits[CLASSIFY_MAX_OWNERS];
    int prev_ids[CLASSIFY_MAX_OWNERS];
    int prev_n = -1;
    alea_ray_interval_finding_t cur = {0};
    int have_cur = 0;
    size_t total = 0;
    int failed = 0;

    double t_prev = 0.0;
    for (size_t i = 0; i <= r.hits.count && !failed; i++) {
        double t_curr = (i < r.hits.count) ? r.hits.data[i].t : t_max;
        if (t_curr > t_max) t_curr = t_max;

        /* Skip slivers below the merge tolerance: grazing/tangential hits. */
        if (t_curr - t_prev > 1e-9) {
            /* Interior probe at an irrational fraction (never lands on a
             * periodic internal plane even for symmetric intervals). */
            double tp = t_prev + 0.381966011250105 * (t_curr - t_prev);
            double px, py, pz;
            alea_ray_point_at(&ray, tp, &px, &py, &pz);

            /* Complete owner set: uncached recursive query — the coherence
             * caches hide co-claimants by design and must stay out. */
            int n = alea_find_all_cells_at_point_recursive(
                sys, px, py, pz, hits, CLASSIFY_MAX_OWNERS);
            if (n < 0) { failed = 1; break; }

            int ids[CLASSIFY_MAX_OWNERS];
            for (int k = 0; k < n; k++) ids[k] = hits[k].cell_id;

            int same = have_cur && n == prev_n &&
                       memcmp(ids, prev_ids, (size_t)n * sizeof(int)) == 0;
            if (same) {
                cur.t_exit = t_curr;
            } else {
                if (have_cur) {
                    if (total < max_out && out) out[total] = cur;
                    total++;
                }

                int kind = ALEA_INTERVAL_OK;
                int cell_id = -1, overlap_cell = -1, depth = -1;
                if (n == 0) {
                    kind = ALEA_INTERVAL_GAP;
                } else {
                    /* Overlap: some depth claimed by more than one cell.
                     * Report the shallowest duplicated depth and its first
                     * two claimants (hits are in DFS deck order). */
                    int dup_depth = -1;
                    for (int a = 0; a < n && dup_depth < 0; a++) {
                        for (int b = a + 1; b < n; b++) {
                            if (hits[b].depth == hits[a].depth) {
                                dup_depth = hits[a].depth;
                                cell_id = hits[a].cell_id;
                                overlap_cell = hits[b].cell_id;
                                break;
                            }
                        }
                    }
                    if (dup_depth >= 0) {
                        kind = ALEA_INTERVAL_OVERLAP;
                        depth = dup_depth;
                    } else {
                        /* Single chain: deepest hit of the first DFS chain */
                        int ti = 0;
                        while (ti + 1 < n &&
                               hits[ti + 1].depth == hits[ti].depth + 1)
                            ti++;
                        cell_id = hits[ti].cell_id;
                        depth = hits[ti].depth;
                        kind = (hits[ti].resolution_flags &
                                ALEA_RESOLVE_UNDEFINED_FILL)
                                   ? ALEA_INTERVAL_UNDEFINED_FILL
                                   : ALEA_INTERVAL_OK;
                    }
                }

                cur.t_enter = t_prev;
                cur.t_exit = t_curr;
                cur.kind = kind;
                cur.cell_id = cell_id;
                cur.overlap_cell_id = overlap_cell;
                cur.depth = depth;
                have_cur = 1;
                prev_n = n;
                memcpy(prev_ids, ids, (size_t)n * sizeof(int));
            }
        }

        t_prev = t_curr;
        if (t_prev >= t_max) break;
    }

    if (have_cur && !failed) {
        if (total < max_out && out) out[total] = cur;
        total++;
    }

    alea_raycast_result_free(&r);
    return failed ? -1 : (int)total;
}
