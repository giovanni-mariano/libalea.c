// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "cluster_internal.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define SLICE_FIELDS (ALEA_SLICE_RASTER_CELL_ID | \
    ALEA_SLICE_RASTER_MATERIAL_ID | ALEA_SLICE_RASTER_UNIVERSE_ID | \
    ALEA_SLICE_RASTER_FILL_UNIVERSE | ALEA_SLICE_RASTER_DENSITY | \
    ALEA_SLICE_RASTER_RESOLUTION_FLAGS)

static size_t row_count(size_t total, size_t rank, size_t ranks) {
    return total / ranks + (rank < total % ranks);
}

static size_t first_row(size_t total, size_t rank, size_t ranks) {
    return rank * (total / ranks) + (rank < total % ranks ? rank : total % ranks);
}

static uint64_t view_fingerprint(const alea_slice_view_t* view) {
    uint64_t hash = UINT64_C(1469598103934665603);
    const double values[] = {
        view->plane.origin[0], view->plane.origin[1], view->plane.origin[2],
        view->plane.normal[0], view->plane.normal[1], view->plane.normal[2],
        view->plane.u_axis[0], view->plane.u_axis[1], view->plane.u_axis[2],
        view->plane.v_axis[0], view->plane.v_axis[1], view->plane.v_axis[2],
        view->u_min, view->u_max, view->v_min, view->v_max
    };
    const unsigned char* bytes = (const unsigned char*)values;
    for (size_t i = 0; i < sizeof(values); ++i) {
        hash ^= bytes[i];
        hash *= UINT64_C(1099511628211);
    }
    return hash ? hash : 1;
}

static int trace_bytes(size_t rows, size_t segments, uint32_t fields,
                       size_t* result) {
    if (rows >= SIZE_MAX / sizeof(uint64_t) ||
        segments > SIZE_MAX / 48) return -1;
    size_t per_segment = 2 * sizeof(double) + sizeof(int32_t);
    if (fields & ALEA_RAY_BATCH_MATERIAL) per_segment += sizeof(int32_t);
    if (fields & ALEA_RAY_BATCH_DENSITY) per_segment += sizeof(double);
    if (fields & ALEA_RAY_BATCH_RESOLUTION_FLAGS) per_segment += sizeof(uint8_t);
    if (fields & ALEA_RAY_BATCH_PROJECTED_OWNER)
        per_segment += 5 * sizeof(int32_t) + sizeof(uint8_t) + sizeof(uint64_t);
    if (segments > (SIZE_MAX - (rows + 1) * sizeof(uint64_t)) / per_segment)
        return -1;
    *result = (rows + 1) * sizeof(uint64_t) + segments * per_segment;
    return 0;
}

alea_cluster_status_t alea_cluster_slice_raster(
        alea_cluster_t* cluster, alea_system_t* sys,
        const alea_slice_view_t* view,
        const alea_slice_raster_options_t* options,
        alea_slice_raster_t* root_output) {
    if (!cluster) return ALEA_CLUSTER_INVALID_ARGUMENT;
    alea_cluster_status_t local = sys && view &&
        isfinite(view->u_min) && isfinite(view->u_max) &&
        isfinite(view->v_min) && isfinite(view->v_max) &&
        view->u_max > view->u_min && view->v_max > view->v_min
        ? ALEA_CLUSTER_OK : ALEA_CLUSTER_INVALID_ARGUMENT;
    uint64_t wire[5] = {0};
    int depth = -1;
    if (cluster->rank == 0) {
        if (!root_output || root_output->struct_size < sizeof(*root_output) ||
            !root_output->nu || !root_output->nv ||
            root_output->nu > SIZE_MAX / root_output->nv ||
            !root_output->fields || (root_output->fields & ~SLICE_FIELDS) ||
            (options && (options->struct_size < sizeof(*options) ||
                         options->projected_depth < -1)))
            local = ALEA_CLUSTER_INVALID_ARGUMENT;
        if (local == ALEA_CLUSTER_OK) {
#define REQUIRE(field, member) \
    if ((root_output->fields & (field)) && !root_output->member) \
        local = ALEA_CLUSTER_INVALID_ARGUMENT
            REQUIRE(ALEA_SLICE_RASTER_CELL_ID, cell_ids);
            REQUIRE(ALEA_SLICE_RASTER_MATERIAL_ID, material_ids);
            REQUIRE(ALEA_SLICE_RASTER_UNIVERSE_ID, universe_ids);
            REQUIRE(ALEA_SLICE_RASTER_FILL_UNIVERSE, fill_universe_ids);
            REQUIRE(ALEA_SLICE_RASTER_DENSITY, densities);
            REQUIRE(ALEA_SLICE_RASTER_RESOLUTION_FLAGS, resolution_flags);
#undef REQUIRE
        }
        if (local == ALEA_CLUSTER_OK) {
            wire[0] = root_output->nu;
            wire[1] = root_output->nv;
            wire[2] = root_output->fields;
            wire[3] = options ? options->max_segments : 0;
            wire[4] = options ? options->max_trace_output_bytes : 0;
            depth = options ? options->projected_depth : -1;
        }
    }
    alea_cluster_status_t status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    if (alea_cluster_backend_broadcast_bytes(cluster->backend, wire,
            sizeof(wire), 0) ||
        alea_cluster_backend_broadcast_int(cluster->backend, &depth, 0)) {
        cluster->usable = 0;
        return ALEA_CLUSTER_BACKEND_ERROR;
    }
    const size_t nu = (size_t)wire[0], nv = (size_t)wire[1];
    const uint32_t fields = (uint32_t)wire[2];
    local = nu && nv && nu <= SIZE_MAX / nv &&
        wire[2] <= UINT32_MAX && !(fields & ~SLICE_FIELDS)
        ? ALEA_CLUSTER_OK : ALEA_CLUSTER_INVALID_ARGUMENT;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    uint64_t model = alea_cluster_system_fingerprint(sys);
    local = model ? ALEA_CLUSTER_OK : ALEA_CLUSTER_COMPUTE_ERROR;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    int match = alea_cluster_fingerprints_match(cluster, model);
    if (match < 0) return ALEA_CLUSTER_BACKEND_ERROR;
    if (!match) return ALEA_CLUSTER_MODEL_MISMATCH;
    match = alea_cluster_fingerprints_match(cluster, view_fingerprint(view));
    if (match < 0) return ALEA_CLUSTER_BACKEND_ERROR;
    if (!match) return ALEA_CLUSTER_INVALID_ARGUMENT;

    const size_t ranks = (size_t)cluster->size;
    const size_t rows = row_count(nv, (size_t)cluster->rank, ranks);
    const size_t begin = first_row(nv, (size_t)cluster->rank, ranks);
    local = rows > SIZE_MAX / nu ? ALEA_CLUSTER_OUTPUT_LIMIT : ALEA_CLUSTER_OK;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    const size_t pixels = rows * nu;
    alea_slice_raster_t raster;
    alea_slice_raster_init(&raster);
    raster.nu = nu;
    raster.nv = rows;
    raster.fields = fields;
#define ALLOC_FIELD(flag, member, type) \
    if (fields & (flag)) raster.member = malloc((pixels ? pixels : 1) * sizeof(type))
    ALLOC_FIELD(ALEA_SLICE_RASTER_CELL_ID, cell_ids, int32_t);
    ALLOC_FIELD(ALEA_SLICE_RASTER_MATERIAL_ID, material_ids, int32_t);
    ALLOC_FIELD(ALEA_SLICE_RASTER_UNIVERSE_ID, universe_ids, int32_t);
    ALLOC_FIELD(ALEA_SLICE_RASTER_FILL_UNIVERSE, fill_universe_ids, int32_t);
    ALLOC_FIELD(ALEA_SLICE_RASTER_DENSITY, densities, double);
    ALLOC_FIELD(ALEA_SLICE_RASTER_RESOLUTION_FLAGS, resolution_flags, uint8_t);
#undef ALLOC_FIELD
    alea_raycast_batch_result_t* compact = rows
        ? alea_raycast_batch_result_create() : NULL;
    uint64_t* segment_counts = cluster->rank == 0
        ? malloc(ranks * sizeof(uint64_t)) : NULL;
    size_t* rank_bytes = cluster->rank == 0
        ? malloc(ranks * sizeof(size_t)) : NULL;
    local = (!rows || compact) &&
        (cluster->rank != 0 || (segment_counts && rank_bytes))
        ? ALEA_CLUSTER_OK : ALEA_CLUSTER_OUT_OF_MEMORY;
#define CHECK_FIELD(flag, member) \
    if ((fields & (flag)) && !raster.member) local = ALEA_CLUSTER_OUT_OF_MEMORY
    CHECK_FIELD(ALEA_SLICE_RASTER_CELL_ID, cell_ids);
    CHECK_FIELD(ALEA_SLICE_RASTER_MATERIAL_ID, material_ids);
    CHECK_FIELD(ALEA_SLICE_RASTER_UNIVERSE_ID, universe_ids);
    CHECK_FIELD(ALEA_SLICE_RASTER_FILL_UNIVERSE, fill_universe_ids);
    CHECK_FIELD(ALEA_SLICE_RASTER_DENSITY, densities);
    CHECK_FIELD(ALEA_SLICE_RASTER_RESOLUTION_FLAGS, resolution_flags);
#undef CHECK_FIELD
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) goto cleanup;

    uint32_t trace_fields = 0;
    if ((fields & (ALEA_SLICE_RASTER_UNIVERSE_ID |
                   ALEA_SLICE_RASTER_FILL_UNIVERSE)) ||
        (depth >= 0 && (fields & (ALEA_SLICE_RASTER_CELL_ID |
                                  ALEA_SLICE_RASTER_MATERIAL_ID))))
        trace_fields |= ALEA_RAY_BATCH_PROJECTED_OWNER;
    if ((fields & ALEA_SLICE_RASTER_MATERIAL_ID) && depth < 0)
        trace_fields |= ALEA_RAY_BATCH_MATERIAL;
    if (fields & ALEA_SLICE_RASTER_DENSITY)
        trace_fields |= ALEA_RAY_BATCH_DENSITY;
    if (fields & ALEA_SLICE_RASTER_RESOLUTION_FLAGS)
        trace_fields |= ALEA_RAY_BATCH_RESOLUTION_FLAGS;
    alea_raycast_batch_options_t trace_options = {
        sizeof(trace_options), trace_fields, depth, 0, 0, 0
    };
    alea_slice_view_t local_view = *view;
    const double dv = (view->v_max - view->v_min) / (double)nv;
    local_view.v_min = view->v_min + (double)begin * dv;
    local_view.v_max = view->v_min + (double)(begin + rows) * dv;
    local = rows && alea_trace_ray_slice_compact(sys, &local_view, rows,
        &trace_options, compact) != 0
        ? ALEA_CLUSTER_COMPUTE_ERROR : ALEA_CLUSTER_OK;
    if (alea_interrupted()) local = ALEA_CLUSTER_INTERRUPTED;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) goto cleanup;
    uint64_t local_segments = rows
        ? (uint64_t)alea_raycast_batch_segment_count(compact) : 0;
    if (alea_cluster_backend_gather_u64(cluster->backend, local_segments,
                                        segment_counts, 0)) {
        cluster->usable = 0;
        status = ALEA_CLUSTER_BACKEND_ERROR;
        goto cleanup;
    }
    local = ALEA_CLUSTER_OK;
    if (cluster->rank == 0) {
        size_t segments = 0;
        for (size_t rank = 0; rank < ranks; ++rank) {
            if (segment_counts[rank] > SIZE_MAX - segments) {
                local = ALEA_CLUSTER_OUTPUT_LIMIT;
                break;
            }
            segments += (size_t)segment_counts[rank];
        }
        size_t bytes = 0;
        if (local == ALEA_CLUSTER_OK &&
            (trace_bytes(nv, segments, trace_fields, &bytes) ||
             (wire[3] && segments > wire[3]) ||
             (wire[4] && bytes > wire[4])))
            local = ALEA_CLUSTER_OUTPUT_LIMIT;
    }
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) goto cleanup;
    local = rows && alea_rasterize_ray_slice_compact(
        &local_view, compact, &raster) != 0
        ? ALEA_CLUSTER_COMPUTE_ERROR : ALEA_CLUSTER_OK;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) goto cleanup;

#define GATHER_FIELD(flag, member, type) do { \
    if (fields & (flag)) { \
        if (cluster->rank == 0) \
            for (size_t rank = 0; rank < ranks; ++rank) \
                rank_bytes[rank] = row_count(nv, rank, ranks) * nu * sizeof(type); \
        if (alea_cluster_backend_gather_bytes(cluster->backend, \
                raster.member, pixels * sizeof(type), \
                cluster->rank == 0 ? root_output->member : NULL, \
                rank_bytes, 0)) { \
            cluster->usable = 0; \
            status = ALEA_CLUSTER_BACKEND_ERROR; goto cleanup; \
        } \
    } \
} while (0)
    GATHER_FIELD(ALEA_SLICE_RASTER_CELL_ID, cell_ids, int32_t);
    GATHER_FIELD(ALEA_SLICE_RASTER_MATERIAL_ID, material_ids, int32_t);
    GATHER_FIELD(ALEA_SLICE_RASTER_UNIVERSE_ID, universe_ids, int32_t);
    GATHER_FIELD(ALEA_SLICE_RASTER_FILL_UNIVERSE, fill_universe_ids, int32_t);
    GATHER_FIELD(ALEA_SLICE_RASTER_DENSITY, densities, double);
    GATHER_FIELD(ALEA_SLICE_RASTER_RESOLUTION_FLAGS, resolution_flags, uint8_t);
#undef GATHER_FIELD
    status = ALEA_CLUSTER_OK;

cleanup:
    free(raster.cell_ids); free(raster.material_ids);
    free(raster.universe_ids); free(raster.fill_universe_ids);
    free(raster.densities); free(raster.resolution_flags);
    alea_raycast_batch_result_destroy(compact);
    free(segment_counts); free(rank_bytes);
    return status;
}

alea_cluster_status_t alea_cluster_slice_stack_stream(
        alea_cluster_t* cluster, alea_system_t* sys,
        const alea_slice_view_t* views, size_t view_count,
        const alea_slice_raster_options_t* options,
        alea_slice_raster_t* root_output,
        alea_cluster_slice_plane_callback_t callback, void* user_data) {
    if (!cluster) return ALEA_CLUSTER_INVALID_ARGUMENT;
    alea_cluster_status_t local = sys && (!view_count || views)
        ? ALEA_CLUSTER_OK : ALEA_CLUSTER_INVALID_ARGUMENT;
    if (cluster->rank == 0 && view_count && (!root_output || !callback))
        local = ALEA_CLUSTER_INVALID_ARGUMENT;
    alea_cluster_status_t status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    uint64_t count_wire = cluster->rank == 0 ? (uint64_t)view_count : 0;
    if (alea_cluster_backend_broadcast_u64(cluster->backend, &count_wire, 0)) {
        cluster->usable = 0;
        return ALEA_CLUSTER_BACKEND_ERROR;
    }
    local = count_wire == (uint64_t)view_count
        ? ALEA_CLUSTER_OK : ALEA_CLUSTER_INVALID_ARGUMENT;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    for (size_t i = 0; i < view_count; ++i) {
        status = alea_cluster_slice_raster(cluster, sys, &views[i],
            cluster->rank == 0 ? options : NULL,
            cluster->rank == 0 ? root_output : NULL);
        if (status != ALEA_CLUSTER_OK) return status;
        local = cluster->rank == 0 && callback(i, root_output, user_data)
            ? ALEA_CLUSTER_INTERRUPTED : ALEA_CLUSTER_OK;
        status = alea_cluster_agree(cluster, local);
        if (status != ALEA_CLUSTER_OK) return status;
    }
    return ALEA_CLUSTER_OK;
}
