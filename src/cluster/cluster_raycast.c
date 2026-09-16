// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "cluster_internal.h"
#include "raycast/batch_internal.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define CLUSTER_RAY_BATCH 1024u

alea_cluster_status_t alea_cluster_raycast_first_segments(
        alea_cluster_t* cluster, alea_system_t* sys,
        const double* origins_xyz, const double* directions_xyz,
        size_t ray_count, double t_max,
        unsigned char* hit, int32_t* cell_ids,
        double* t_enter, double* t_exit) {
    if (!cluster) return ALEA_CLUSTER_INVALID_ARGUMENT;
    alea_cluster_status_t local = sys ? ALEA_CLUSTER_OK
                                      : ALEA_CLUSTER_INVALID_ARGUMENT;
    if (cluster->rank == 0 &&
        ((ray_count && (!origins_xyz || !directions_xyz || !hit || !cell_ids ||
                        !t_enter || !t_exit)) || ray_count > SIZE_MAX / 3 ||
         !isfinite(t_max)))
        local = ALEA_CLUSTER_INVALID_ARGUMENT;
    alea_cluster_status_t status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;

    uint64_t wire_count = (uint64_t)ray_count;
    if (alea_cluster_backend_broadcast_u64(cluster->backend, &wire_count, 0) ||
        alea_cluster_backend_broadcast_bytes(cluster->backend, &t_max,
                                             sizeof(t_max), 0)) {
        cluster->usable = 0;
        return ALEA_CLUSTER_BACKEND_ERROR;
    }
    local = wire_count > (uint64_t)(SIZE_MAX / 3)
        ? ALEA_CLUSTER_INVALID_ARGUMENT : ALEA_CLUSTER_OK;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    const size_t total = (size_t)wire_count;

    uint64_t fingerprint = alea_cluster_system_fingerprint(sys);
    local = fingerprint ? ALEA_CLUSTER_OK : ALEA_CLUSTER_COMPUTE_ERROR;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    int matching = alea_cluster_fingerprints_match(cluster, fingerprint);
    if (matching < 0) return ALEA_CLUSTER_BACKEND_ERROR;
    if (!matching) return ALEA_CLUSTER_MODEL_MISMATCH;

    double* inputs = malloc(6 * CLUSTER_RAY_BATCH * sizeof(double));
    double* reduced = malloc(4 * CLUSTER_RAY_BATCH * sizeof(double));
    alea_raycast_batch_result_t* batch = alea_raycast_batch_result_create();
    local = inputs && reduced && batch ? ALEA_CLUSTER_OK
                                      : ALEA_CLUSTER_OUT_OF_MEMORY;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) goto done;

    for (size_t offset = 0; offset < total; ) {
        const size_t count = total - offset < CLUSTER_RAY_BATCH
            ? total - offset : CLUSTER_RAY_BATCH;
        if (cluster->rank == 0) {
            memcpy(inputs, origins_xyz + 3 * offset, 3 * count * sizeof(double));
            memcpy(inputs + 3 * count, directions_xyz + 3 * offset,
                   3 * count * sizeof(double));
        }
        if (alea_cluster_backend_broadcast_bytes(cluster->backend, inputs,
                6 * count * sizeof(double), 0)) {
            cluster->usable = 0;
            status = ALEA_CLUSTER_BACKEND_ERROR;
            goto done;
        }
        const size_t ranks = (size_t)cluster->size;
        const size_t rank = (size_t)cluster->rank;
        const size_t base = count / ranks;
        const size_t extra = count % ranks;
        const size_t local_begin = rank * base + (rank < extra ? rank : extra);
        const size_t local_count = base + (rank < extra ? 1u : 0u);
        local = alea_raycast_hier_batch(sys,
            inputs + 3 * local_begin,
            inputs + 3 * count + 3 * local_begin,
            local_count, t_max, NULL, batch) == 0
            ? ALEA_CLUSTER_OK : ALEA_CLUSTER_COMPUTE_ERROR;
        if (alea_interrupted()) local = ALEA_CLUSTER_INTERRUPTED;
        status = alea_cluster_agree(cluster, local);
        if (status != ALEA_CLUSTER_OK) goto done;

        memset(reduced, 0, 4 * count * sizeof(double));
        const uint64_t* offsets = alea_raycast_batch_ray_offsets(batch);
        const int32_t* cells = alea_raycast_batch_cell_ids(batch);
        const double* enters = alea_raycast_batch_t_enter(batch);
        const double* exits = alea_raycast_batch_t_exit(batch);
        for (size_t i = 0; i < local_count; ++i) {
            size_t src = (size_t)offsets[i];
            const size_t end = (size_t)offsets[i + 1];
            while (src < end && cells[src] < 0) ++src;
            if (src == end) continue;
            const size_t dst = 4 * (local_begin + i);
            reduced[dst] = 1.0;
            reduced[dst + 1] = (double)cells[src];
            reduced[dst + 2] = enters[src];
            reduced[dst + 3] = exits[src];
        }
        if (alea_cluster_backend_sum_doubles(cluster->backend,
                                              reduced, 4 * count)) {
            cluster->usable = 0;
            status = ALEA_CLUSTER_BACKEND_ERROR;
            goto done;
        }
        if (cluster->rank == 0) {
            for (size_t i = 0; i < count; ++i) {
                hit[offset + i] = reduced[4 * i] != 0.0;
                if (!hit[offset + i]) continue;
                cell_ids[offset + i] = (int32_t)reduced[4 * i + 1];
                t_enter[offset + i] = reduced[4 * i + 2];
                t_exit[offset + i] = reduced[4 * i + 3];
            }
        }
        offset += count;
    }
    status = ALEA_CLUSTER_OK;

done:
    alea_raycast_batch_result_destroy(batch);
    free(inputs);
    free(reduced);
    return status;
}

static size_t cluster_rank_count(size_t count, size_t rank, size_t ranks) {
    return count / ranks + (rank < count % ranks ? 1u : 0u);
}

static int cluster_add_bytes(size_t* total, size_t count, size_t element_size) {
    if (count > PTRDIFF_MAX / element_size ||
        count * element_size > SIZE_MAX - *total) return -1;
    *total += count * element_size;
    return 0;
}

static int cluster_result_bytes(size_t rays, size_t segments, size_t paths,
                                uint32_t fields, size_t* result) {
    size_t bytes = 0;
#define ADD(count, type) if (cluster_add_bytes(&bytes, (count), sizeof(type))) return -1
    if (rays == SIZE_MAX ||
        ((fields & ALEA_RAY_BATCH_FULL_PATHS) && segments == SIZE_MAX))
        return -1;
    ADD(rays + 1, uint64_t);
    ADD(segments, double);
    ADD(segments, double);
    ADD(segments, int32_t);
    if (fields & ALEA_RAY_BATCH_MATERIAL) ADD(segments, int32_t);
    if (fields & ALEA_RAY_BATCH_DENSITY) ADD(segments, double);
    if (fields & ALEA_RAY_BATCH_SURFACES) {
        ADD(segments, int32_t); ADD(segments, int32_t);
    }
    if (fields & ALEA_RAY_BATCH_RESOLUTION_FLAGS) ADD(segments, uint8_t);
    if (fields & ALEA_RAY_BATCH_PROJECTED_OWNER) {
        for (int i = 0; i < 5; ++i) { ADD(segments, int32_t); }
        ADD(segments, uint8_t); ADD(segments, uint64_t);
    }
    if (fields & ALEA_RAY_BATCH_FULL_PATHS) {
        ADD(segments + 1, uint64_t);
        for (int i = 0; i < 5; ++i) { ADD(paths, int32_t); }
        ADD(paths, uint8_t);
        if (paths > SIZE_MAX / 3) return -1;
        ADD(paths * 3, double);
        ADD(paths, uint64_t);
    }
#undef ADD
    *result = bytes;
    return 0;
}

static int cluster_grow_result(alea_raycast_batch_result_t* result,
                               size_t segments, size_t paths) {
#define GROW(member, count) do { \
    if ((count) > PTRDIFF_MAX / sizeof(*result->member)) return -1; \
    size_t bytes = (count) * sizeof(*result->member); \
    void* grown = realloc(result->member, bytes ? bytes : 1); \
    if (!grown) return -1; \
    result->member = grown; \
} while (0)
    GROW(t_enter, segments); GROW(t_exit, segments); GROW(cell_ids, segments);
    if (result->fields & ALEA_RAY_BATCH_MATERIAL) GROW(material_ids, segments);
    if (result->fields & ALEA_RAY_BATCH_DENSITY) GROW(densities, segments);
    if (result->fields & ALEA_RAY_BATCH_SURFACES) {
        GROW(enter_surface_ids, segments); GROW(exit_surface_ids, segments);
    }
    if (result->fields & ALEA_RAY_BATCH_RESOLUTION_FLAGS)
        GROW(resolution_flags, segments);
    if (result->fields & ALEA_RAY_BATCH_PROJECTED_OWNER) {
        GROW(projected_cell_ids, segments);
        GROW(projected_material_ids, segments);
        GROW(projected_universe_ids, segments);
        GROW(projected_fill_universes, segments);
        GROW(projected_depths, segments);
        GROW(projected_is_lattice, segments);
        GROW(projected_occurrence_keys, segments);
    }
    if (result->fields & ALEA_RAY_BATCH_FULL_PATHS) {
        GROW(segment_path_offsets, segments + 1);
        GROW(path_cell_ids, paths); GROW(path_material_ids, paths);
        GROW(path_universe_ids, paths); GROW(path_fill_universes, paths);
        GROW(path_depths, paths); GROW(path_is_lattice, paths);
        GROW(path_lattice_origins_xyz, paths * 3);
        GROW(path_occurrence_keys, paths);
    }
#undef GROW
    return 0;
}

static int cluster_gather_field(alea_cluster_t* cluster,
                                const void* local_data, size_t local_count,
                                void* root_data, const uint64_t* root_counts,
                                size_t element_size, size_t* rank_bytes) {
    size_t local_bytes = local_count * element_size;
    if (cluster->rank == 0) {
        for (int peer = 0; peer < cluster->size; ++peer)
            rank_bytes[peer] = (size_t)root_counts[peer] * element_size;
    }
    if (alea_cluster_backend_gather_bytes(cluster->backend, local_data,
            local_bytes, root_data, rank_bytes, 0)) {
        cluster->usable = 0;
        return -1;
    }
    return 0;
}

static alea_cluster_status_t cluster_raycast_batch_impl(
        alea_cluster_t* cluster, alea_system_t* sys,
        const double* origins_xyz, const double* directions_xyz,
        size_t ray_count, double t_max,
        const alea_raycast_batch_options_t* options,
        alea_raycast_batch_result_t* result,
        alea_cluster_ray_batch_callback_t callback, void* user_data,
        int streaming) {
    if (!cluster) return ALEA_CLUSTER_INVALID_ARGUMENT;
    const uint32_t known_fields = ALEA_RAY_BATCH_MATERIAL |
        ALEA_RAY_BATCH_DENSITY | ALEA_RAY_BATCH_SURFACES |
        ALEA_RAY_BATCH_RESOLUTION_FLAGS | ALEA_RAY_BATCH_PROJECTED_OWNER |
        ALEA_RAY_BATCH_FULL_PATHS;
    uint64_t wire[5] = {0};
    int projected_depth = -1;
    alea_cluster_status_t local = sys ? ALEA_CLUSTER_OK
                                      : ALEA_CLUSTER_INVALID_ARGUMENT;
    if (cluster->rank == 0) {
        if ((streaming ? !callback : !result) ||
            (ray_count && (!origins_xyz || !directions_xyz)) ||
            ray_count > SIZE_MAX / 3 || !isfinite(t_max) ||
            (options && (options->struct_size < sizeof(*options) ||
                         (options->fields & ~known_fields) ||
                         options->projected_depth < -1)))
            local = ALEA_CLUSTER_INVALID_ARGUMENT;
        if (local == ALEA_CLUSTER_OK) {
            wire[0] = (uint64_t)ray_count;
            wire[1] = options ? options->fields : 0;
            wire[2] = options ? options->max_segments : 0;
            wire[3] = options ? options->max_path_entries : 0;
            wire[4] = options ? options->max_output_bytes : 0;
            projected_depth = options ? options->projected_depth : -1;
        }
    }
    alea_cluster_status_t status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    if (alea_cluster_backend_broadcast_bytes(cluster->backend, wire,
            sizeof(wire), 0) ||
        alea_cluster_backend_broadcast_bytes(cluster->backend, &t_max,
            sizeof(t_max), 0) ||
        alea_cluster_backend_broadcast_int(cluster->backend, &projected_depth,
            0)) {
        cluster->usable = 0;
        return ALEA_CLUSTER_BACKEND_ERROR;
    }
    local = wire[0] > (uint64_t)(SIZE_MAX / 3) ||
            wire[1] > UINT32_MAX || !isfinite(t_max)
        ? ALEA_CLUSTER_INVALID_ARGUMENT : ALEA_CLUSTER_OK;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    const size_t total = (size_t)wire[0];
    const uint32_t fields = (uint32_t)wire[1];
    uint64_t fingerprint = alea_cluster_system_fingerprint(sys);
    local = fingerprint ? ALEA_CLUSTER_OK : ALEA_CLUSTER_COMPUTE_ERROR;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    int matching = alea_cluster_fingerprints_match(cluster, fingerprint);
    if (matching < 0) return ALEA_CLUSTER_BACKEND_ERROR;
    if (!matching) return ALEA_CLUSTER_MODEL_MISMATCH;

    alea_raycast_batch_options_t local_options = {
        sizeof(local_options), fields, projected_depth, wire[2], wire[3], wire[4]
    };
    const size_t ranks = (size_t)cluster->size;
    double* local_origins = malloc(3 * CLUSTER_RAY_BATCH * sizeof(double));
    double* local_directions = malloc(3 * CLUSTER_RAY_BATCH * sizeof(double));
    alea_raycast_batch_result_t* batch = alea_raycast_batch_result_create();
    alea_raycast_batch_result_t* staged = cluster->rank == 0
        ? alea_raycast_batch_result_create() : NULL;
    uint64_t* segment_counts = cluster->rank == 0
        ? calloc(ranks, sizeof(uint64_t)) : NULL;
    uint64_t* path_counts = cluster->rank == 0
        ? calloc(ranks, sizeof(uint64_t)) : NULL;
    uint64_t* ray_counts = cluster->rank == 0
        ? calloc(ranks, sizeof(uint64_t)) : NULL;
    size_t* rank_bytes = cluster->rank == 0
        ? calloc(ranks, sizeof(size_t)) : NULL;
    local = local_origins && local_directions && batch &&
        (cluster->rank != 0 || (staged && segment_counts && path_counts &&
                                ray_counts && rank_bytes))
        ? ALEA_CLUSTER_OK : ALEA_CLUSTER_OUT_OF_MEMORY;
    if (local == ALEA_CLUSTER_OK && cluster->rank == 0) {
        size_t initial_bytes;
        const size_t initial_rays = streaming && total > CLUSTER_RAY_BATCH
            ? CLUSTER_RAY_BATCH : total;
        if (cluster_result_bytes(initial_rays, 0, 0, fields,
                                 &initial_bytes) ||
            (wire[4] && initial_bytes > wire[4]))
            local = ALEA_CLUSTER_OUTPUT_LIMIT;
        else {
            staged->ray_count = initial_rays;
            staged->fields = fields;
            staged->ray_offsets = malloc((initial_rays + 1) * sizeof(uint64_t));
            if (!staged->ray_offsets || cluster_grow_result(staged, 0, 0))
                local = ALEA_CLUSTER_OUT_OF_MEMORY;
            else {
                staged->ray_offsets[0] = 0;
                if (staged->segment_path_offsets)
                    staged->segment_path_offsets[0] = 0;
            }
        }
    }
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) goto cleanup_batch;

    for (size_t offset = 0; offset < total; ) {
        const size_t count = total - offset < CLUSTER_RAY_BATCH
            ? total - offset : CLUSTER_RAY_BATCH;
        const size_t result_offset = streaming ? 0 : offset;
        if (streaming && cluster->rank == 0) {
            staged->ray_count = count;
            staged->segment_count = 0;
            staged->path_entry_count = 0;
            staged->ray_offsets[0] = 0;
            if (staged->segment_path_offsets)
                staged->segment_path_offsets[0] = 0;
        }
        const size_t local_count = cluster_rank_count(
            count, (size_t)cluster->rank, ranks);
        if (alea_cluster_backend_scatter_blocks(cluster->backend,
                cluster->rank == 0 ? origins_xyz + 3 * offset : NULL,
                local_origins, 3 * sizeof(double), count, 0) ||
            alea_cluster_backend_scatter_blocks(cluster->backend,
                cluster->rank == 0 ? directions_xyz + 3 * offset : NULL,
                local_directions, 3 * sizeof(double), count, 0)) {
            cluster->usable = 0;
            status = ALEA_CLUSTER_BACKEND_ERROR;
            goto cleanup_batch;
        }
        int trace_status = alea_raycast_hier_batch(
            sys, local_origins, local_directions, local_count, t_max,
            &local_options, batch);
        local = trace_status == 0 ? ALEA_CLUSTER_OK :
            alea_error_code() == ALEA_ERR_OVERFLOW
                ? ALEA_CLUSTER_OUTPUT_LIMIT : ALEA_CLUSTER_COMPUTE_ERROR;
        if (alea_interrupted()) local = ALEA_CLUSTER_INTERRUPTED;
        status = alea_cluster_agree(cluster, local);
        if (status != ALEA_CLUSTER_OK) goto cleanup_batch;

        const size_t local_segments = batch->segment_count;
        const size_t local_paths = batch->path_entry_count;
        if (alea_cluster_backend_gather_u64(cluster->backend,
                (uint64_t)local_segments, segment_counts, 0) ||
            alea_cluster_backend_gather_u64(cluster->backend,
                (uint64_t)local_paths, path_counts, 0)) {
            cluster->usable = 0;
            status = ALEA_CLUSTER_BACKEND_ERROR;
            goto cleanup_batch;
        }
        size_t old_segments = 0, old_paths = 0;
        local = ALEA_CLUSTER_OK;
        if (cluster->rank == 0) {
            old_segments = staged->segment_count;
            old_paths = staged->path_entry_count;
            size_t new_segments = old_segments;
            size_t new_paths = old_paths;
            for (size_t peer = 0; peer < ranks; ++peer) {
                ray_counts[peer] = cluster_rank_count(count, peer, ranks);
                if (segment_counts[peer] > SIZE_MAX - new_segments ||
                    path_counts[peer] > SIZE_MAX - new_paths) {
                    local = ALEA_CLUSTER_OUTPUT_LIMIT;
                    break;
                }
                new_segments += (size_t)segment_counts[peer];
                new_paths += (size_t)path_counts[peer];
            }
            if (local == ALEA_CLUSTER_OK) {
                size_t output_bytes;
                if (cluster_result_bytes(streaming ? count : total,
                                         new_segments, new_paths,
                                         fields, &output_bytes) ||
                    (wire[2] && new_segments > wire[2]) ||
                    (wire[3] && new_paths > wire[3]) ||
                    (wire[4] && output_bytes > wire[4]))
                    local = ALEA_CLUSTER_OUTPUT_LIMIT;
                else if (cluster_grow_result(staged, new_segments, new_paths))
                    local = ALEA_CLUSTER_OUT_OF_MEMORY;
                else {
                    staged->segment_count = new_segments;
                    staged->path_entry_count = new_paths;
                }
            }
        }
        status = alea_cluster_agree(cluster, local);
        if (status != ALEA_CLUSTER_OK) goto cleanup_batch;

#define GATHER(source, destination, local_items, counts, type) do { \
    if (cluster_gather_field(cluster, (source), (local_items), \
            cluster->rank == 0 ? (destination) : NULL, (counts), \
            sizeof(type), rank_bytes)) { \
        status = ALEA_CLUSTER_BACKEND_ERROR; goto cleanup_batch; \
    } \
} while (0)
#define SEGMENT_FIELD(member, type) GATHER(batch->member, \
    staged ? staged->member + old_segments : NULL, local_segments, \
    segment_counts, type)
#define PATH_FIELD(member, type) GATHER(batch->member, \
    staged ? staged->member + old_paths : NULL, local_paths, path_counts, type)
        GATHER(batch->ray_offsets,
               staged ? staged->ray_offsets + result_offset : NULL,
               local_count, ray_counts, uint64_t);
        if (cluster->rank == 0) {
            size_t row = result_offset, preceding = old_segments;
            for (size_t peer = 0; peer < ranks; ++peer) {
                for (size_t i = 0; i < ray_counts[peer]; ++i)
                    staged->ray_offsets[row++] += preceding;
                preceding += (size_t)segment_counts[peer];
            }
            staged->ray_offsets[result_offset + count] = staged->segment_count;
        }
        if (fields & ALEA_RAY_BATCH_FULL_PATHS) {
            GATHER(batch->segment_path_offsets,
                   staged ? staged->segment_path_offsets + old_segments : NULL,
                   local_segments, segment_counts, uint64_t);
            if (cluster->rank == 0) {
                size_t segment = old_segments, preceding = old_paths;
                for (size_t peer = 0; peer < ranks; ++peer) {
                    for (size_t i = 0; i < segment_counts[peer]; ++i)
                        staged->segment_path_offsets[segment++] += preceding;
                    preceding += (size_t)path_counts[peer];
                }
                staged->segment_path_offsets[staged->segment_count] =
                    staged->path_entry_count;
            }
        }
        SEGMENT_FIELD(t_enter, double);
        SEGMENT_FIELD(t_exit, double);
        SEGMENT_FIELD(cell_ids, int32_t);
        if (fields & ALEA_RAY_BATCH_MATERIAL)
            SEGMENT_FIELD(material_ids, int32_t);
        if (fields & ALEA_RAY_BATCH_DENSITY)
            SEGMENT_FIELD(densities, double);
        if (fields & ALEA_RAY_BATCH_SURFACES) {
            SEGMENT_FIELD(enter_surface_ids, int32_t);
            SEGMENT_FIELD(exit_surface_ids, int32_t);
        }
        if (fields & ALEA_RAY_BATCH_RESOLUTION_FLAGS)
            SEGMENT_FIELD(resolution_flags, uint8_t);
        if (fields & ALEA_RAY_BATCH_PROJECTED_OWNER) {
            SEGMENT_FIELD(projected_cell_ids, int32_t);
            SEGMENT_FIELD(projected_material_ids, int32_t);
            SEGMENT_FIELD(projected_universe_ids, int32_t);
            SEGMENT_FIELD(projected_fill_universes, int32_t);
            SEGMENT_FIELD(projected_depths, int32_t);
            SEGMENT_FIELD(projected_is_lattice, uint8_t);
            SEGMENT_FIELD(projected_occurrence_keys, uint64_t);
        }
        if (fields & ALEA_RAY_BATCH_FULL_PATHS) {
            PATH_FIELD(path_cell_ids, int32_t);
            PATH_FIELD(path_material_ids, int32_t);
            PATH_FIELD(path_universe_ids, int32_t);
            PATH_FIELD(path_fill_universes, int32_t);
            PATH_FIELD(path_depths, int32_t);
            PATH_FIELD(path_is_lattice, uint8_t);
            GATHER(batch->path_lattice_origins_xyz,
                   staged ? staged->path_lattice_origins_xyz + old_paths * 3
                          : NULL,
                   local_paths, path_counts, double[3]);
            PATH_FIELD(path_occurrence_keys, uint64_t);
        }
#undef PATH_FIELD
#undef SEGMENT_FIELD
#undef GATHER
        if (streaming) {
            local = cluster->rank == 0 && callback(offset, staged, user_data)
                ? ALEA_CLUSTER_INTERRUPTED : ALEA_CLUSTER_OK;
            status = alea_cluster_agree(cluster, local);
            if (status != ALEA_CLUSTER_OK) goto cleanup_batch;
        }
        offset += count;
    }
    if (!streaming && cluster->rank == 0)
        alea_raycast_batch_result_replace_internal(result, staged);
    status = ALEA_CLUSTER_OK;

cleanup_batch:
    alea_raycast_batch_result_destroy(batch);
    alea_raycast_batch_result_destroy(staged);
    free(local_origins); free(local_directions);
    free(segment_counts); free(path_counts); free(ray_counts); free(rank_bytes);
    return status;
}

alea_cluster_status_t alea_cluster_raycast_batch(
        alea_cluster_t* cluster, alea_system_t* sys,
        const double* origins_xyz, const double* directions_xyz,
        size_t ray_count, double t_max,
        const alea_raycast_batch_options_t* options,
        alea_raycast_batch_result_t* result) {
    return cluster_raycast_batch_impl(cluster, sys, origins_xyz,
        directions_xyz, ray_count, t_max, options, result, NULL, NULL, 0);
}

alea_cluster_status_t alea_cluster_raycast_batch_stream(
        alea_cluster_t* cluster, alea_system_t* sys,
        const double* origins_xyz, const double* directions_xyz,
        size_t ray_count, double t_max,
        const alea_raycast_batch_options_t* options,
        alea_cluster_ray_batch_callback_t callback, void* user_data) {
    return cluster_raycast_batch_impl(cluster, sys, origins_xyz,
        directions_xyz, ray_count, t_max, options, NULL,
        callback, user_data, 1);
}

alea_cluster_status_t alea_cluster_raycast_batch_shards(
        alea_cluster_t* cluster, alea_system_t* sys,
        const double* origins_xyz, const double* directions_xyz,
        size_t ray_count, double t_max,
        const alea_raycast_batch_options_t* options,
        alea_cluster_ray_shard_callback_t callback, void* user_data) {
    if (!cluster) return ALEA_CLUSTER_INVALID_ARGUMENT;
    const uint32_t known_fields = ALEA_RAY_BATCH_MATERIAL |
        ALEA_RAY_BATCH_DENSITY | ALEA_RAY_BATCH_SURFACES |
        ALEA_RAY_BATCH_RESOLUTION_FLAGS | ALEA_RAY_BATCH_PROJECTED_OWNER |
        ALEA_RAY_BATCH_FULL_PATHS;
    alea_cluster_status_t local = sys && callback ? ALEA_CLUSTER_OK
                                                 : ALEA_CLUSTER_INVALID_ARGUMENT;
    uint64_t wire[5] = {0};
    int projected_depth = -1;
    if (cluster->rank == 0) {
        if ((ray_count && (!origins_xyz || !directions_xyz)) ||
            ray_count > SIZE_MAX / 3 || !isfinite(t_max) ||
            (options && (options->struct_size < sizeof(*options) ||
                         (options->fields & ~known_fields) ||
                         options->projected_depth < -1)))
            local = ALEA_CLUSTER_INVALID_ARGUMENT;
        if (local == ALEA_CLUSTER_OK) {
            wire[0] = ray_count;
            wire[1] = options ? options->fields : 0;
            wire[2] = options ? options->max_segments : 0;
            wire[3] = options ? options->max_path_entries : 0;
            wire[4] = options ? options->max_output_bytes : 0;
            projected_depth = options ? options->projected_depth : -1;
        }
    }
    alea_cluster_status_t status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    if (alea_cluster_backend_broadcast_bytes(cluster->backend, wire,
            sizeof(wire), 0) ||
        alea_cluster_backend_broadcast_bytes(cluster->backend, &t_max,
            sizeof(t_max), 0) ||
        alea_cluster_backend_broadcast_int(cluster->backend,
            &projected_depth, 0)) {
        cluster->usable = 0;
        return ALEA_CLUSTER_BACKEND_ERROR;
    }
    local = wire[0] > (uint64_t)(SIZE_MAX / 3) ||
            wire[1] > UINT32_MAX || !isfinite(t_max)
        ? ALEA_CLUSTER_INVALID_ARGUMENT : ALEA_CLUSTER_OK;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    uint64_t fingerprint = alea_cluster_system_fingerprint(sys);
    local = fingerprint ? ALEA_CLUSTER_OK : ALEA_CLUSTER_COMPUTE_ERROR;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    int matching = alea_cluster_fingerprints_match(cluster, fingerprint);
    if (matching < 0) return ALEA_CLUSTER_BACKEND_ERROR;
    if (!matching) return ALEA_CLUSTER_MODEL_MISMATCH;

    alea_raycast_batch_options_t local_options = {
        sizeof(local_options), (uint32_t)wire[1], projected_depth,
        wire[2], wire[3], wire[4]
    };
    double* local_origins = malloc(3 * CLUSTER_RAY_BATCH * sizeof(double));
    double* local_directions = malloc(3 * CLUSTER_RAY_BATCH * sizeof(double));
    alea_raycast_batch_result_t* shard = alea_raycast_batch_result_create();
    local = local_origins && local_directions && shard
        ? ALEA_CLUSTER_OK : ALEA_CLUSTER_OUT_OF_MEMORY;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) goto cleanup_shards;
    const size_t total = (size_t)wire[0];
    const size_t ranks = (size_t)cluster->size;
    const size_t rank = (size_t)cluster->rank;
    for (size_t offset = 0; offset < total; ) {
        const size_t count = total - offset < CLUSTER_RAY_BATCH
            ? total - offset : CLUSTER_RAY_BATCH;
        const size_t mine = cluster_rank_count(count, rank, ranks);
        const size_t first = rank * (count / ranks) +
            (rank < count % ranks ? rank : count % ranks);
        if (alea_cluster_backend_scatter_blocks(cluster->backend,
                cluster->rank == 0 ? origins_xyz + 3 * offset : NULL,
                local_origins, 3 * sizeof(double), count, 0) ||
            alea_cluster_backend_scatter_blocks(cluster->backend,
                cluster->rank == 0 ? directions_xyz + 3 * offset : NULL,
                local_directions, 3 * sizeof(double), count, 0)) {
            cluster->usable = 0;
            status = ALEA_CLUSTER_BACKEND_ERROR;
            goto cleanup_shards;
        }
        int rc = alea_raycast_hier_batch(sys, local_origins,
            local_directions, mine, t_max, &local_options, shard);
        local = rc == 0 ? ALEA_CLUSTER_OK :
            alea_error_code() == ALEA_ERR_OVERFLOW
                ? ALEA_CLUSTER_OUTPUT_LIMIT : ALEA_CLUSTER_COMPUTE_ERROR;
        if (alea_interrupted()) local = ALEA_CLUSTER_INTERRUPTED;
        status = alea_cluster_agree(cluster, local);
        if (status != ALEA_CLUSTER_OK) goto cleanup_shards;
        local = mine && callback(offset + first, shard, user_data)
            ? ALEA_CLUSTER_INTERRUPTED : ALEA_CLUSTER_OK;
        status = alea_cluster_agree(cluster, local);
        if (status != ALEA_CLUSTER_OK) goto cleanup_shards;
        offset += count;
    }
    status = ALEA_CLUSTER_OK;

cleanup_shards:
    alea_raycast_batch_result_destroy(shard);
    free(local_origins);
    free(local_directions);
    return status;
}
