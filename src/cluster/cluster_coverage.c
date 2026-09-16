// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "cluster_internal.h"
#include "raycast/raycast.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define COVERAGE_BATCH_ROWS 256u

static size_t rank_count(size_t total, size_t rank, size_t ranks) {
    return total / ranks + (rank < total % ranks);
}

static size_t rank_first(size_t total, size_t rank, size_t ranks) {
    return rank * (total / ranks) + (rank < total % ranks ? rank : total % ranks);
}

static int coverage_bytes(size_t rows, size_t intervals, size_t owners,
                          size_t* bytes) {
    if (rows == SIZE_MAX || intervals == SIZE_MAX) return -1;
    size_t count = 0;
#define ADD(n, type) do { \
    if ((n) > (SIZE_MAX - count) / sizeof(type)) return -1; \
    count += (n) * sizeof(type); \
} while (0)
    ADD(rows + 1, size_t); ADD(rows, uint8_t); ADD(rows, double);
    ADD(intervals, double); ADD(intervals, double); ADD(intervals, uint8_t);
    ADD(intervals + 1, size_t); ADD(intervals, size_t);
    for (int i = 0; i < 6; ++i) ADD(owners, int);
    ADD(owners, uint64_t); ADD(owners, uint64_t); ADD(owners, uint8_t);
#undef ADD
    *bytes = count;
    return 0;
}

static int gather_field(alea_cluster_t* cluster, const void* local_data,
        size_t local_count, void* root_data, const uint64_t* root_counts,
        size_t element_size, size_t* rank_bytes) {
    if (cluster->rank == 0)
        for (int rank = 0; rank < cluster->size; ++rank)
            rank_bytes[rank] = (size_t)root_counts[rank] * element_size;
    if (alea_cluster_backend_gather_bytes(cluster->backend, local_data,
            local_count * element_size, root_data, rank_bytes, 0)) {
        cluster->usable = 0;
        return -1;
    }
    return 0;
}

static alea_cluster_status_t assemble_coverage_batch(
        alea_cluster_t* cluster,
        const alea_ray_coverage_slice_result_t* shard,
        size_t first_row, size_t rows, size_t local_row_count,
        const alea_ray_coverage_slice_options_t* options,
        alea_ray_coverage_slice_result_t* staged,
        uint64_t* row_counts, uint64_t* interval_counts,
        uint64_t* owner_counts, size_t* rank_bytes,
        alea_cluster_coverage_shard_callback_t callback, void* user_data) {
    uint64_t local_rows = (uint64_t)local_row_count;
    uint64_t local_intervals = local_row_count
        ? (uint64_t)shard->interval_count : 0;
    uint64_t local_owners = local_row_count
        ? (uint64_t)shard->owner_count : 0;
    if (alea_cluster_backend_gather_u64(cluster->backend, local_rows,
            row_counts, 0) ||
        alea_cluster_backend_gather_u64(cluster->backend, local_intervals,
            interval_counts, 0) ||
        alea_cluster_backend_gather_u64(cluster->backend, local_owners,
            owner_counts, 0)) {
        cluster->usable = 0;
        return ALEA_CLUSTER_BACKEND_ERROR;
    }
    alea_cluster_status_t local = ALEA_CLUSTER_OK;
    size_t intervals = 0, owners = 0;
    if (cluster->rank == 0) {
        size_t counted_rows = 0;
        for (int rank = 0; rank < cluster->size; ++rank) {
            if (row_counts[rank] > rows - counted_rows ||
                interval_counts[rank] > SIZE_MAX - intervals ||
                owner_counts[rank] > SIZE_MAX - owners) {
                local = ALEA_CLUSTER_OUTPUT_LIMIT;
                break;
            }
            counted_rows += (size_t)row_counts[rank];
            intervals += (size_t)interval_counts[rank];
            owners += (size_t)owner_counts[rank];
        }
        size_t bytes = 0;
        if (local == ALEA_CLUSTER_OK &&
            (counted_rows != rows ||
             coverage_bytes(rows, intervals, owners, &bytes) ||
             (options->max_rows && rows > options->max_rows) ||
             (options->max_intervals &&
              intervals > options->max_intervals) ||
             (options->max_owners && owners > options->max_owners) ||
             (options->max_output_bytes &&
              bytes > options->max_output_bytes)))
            local = ALEA_CLUSTER_OUTPUT_LIMIT;
        if (local == ALEA_CLUSTER_OK) {
            alea_ray_coverage_slice_result_free(staged);
            staged->row_count = rows;
            staged->interval_count = intervals;
            staged->owner_count = owners;
#define ALLOC(member, count, type) do { \
    staged->member = malloc(((count) ? (count) : 1) * sizeof(type)); \
    if (!staged->member) local = ALEA_CLUSTER_OUT_OF_MEMORY; \
} while (0)
            ALLOC(row_offsets, rows + 1, size_t);
            ALLOC(row_direction_tags, rows, uint8_t);
            ALLOC(row_transverse_coordinates, rows, double);
            ALLOC(t_enter, intervals, double); ALLOC(t_exit, intervals, double);
            ALLOC(kinds, intervals, uint8_t);
            ALLOC(owner_offsets, intervals + 1, size_t);
            ALLOC(owner_count_lower_bounds, intervals, size_t);
            ALLOC(owner_cell_ids, owners, int);
            ALLOC(owner_cell_indices, owners, int);
            ALLOC(owner_material_ids, owners, int);
            ALLOC(owner_universe_ids, owners, int);
            ALLOC(owner_fill_universes, owners, int);
            ALLOC(owner_depths, owners, int);
            ALLOC(owner_occurrence_keys, owners, uint64_t);
            ALLOC(owner_parent_occurrence_keys, owners, uint64_t);
            ALLOC(owner_resolution_flags, owners, uint8_t);
#undef ALLOC
        }
    }
    alea_cluster_status_t status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
#define GATHER(member, count, counts, type) do { \
    if (gather_field(cluster, shard->member, (count), \
            staged ? staged->member : NULL, (counts), \
            sizeof(type), rank_bytes)) return ALEA_CLUSTER_BACKEND_ERROR; \
} while (0)
    GATHER(row_offsets, local_rows, row_counts, size_t);
    if (cluster->rank == 0) {
        size_t row = 0, preceding = 0;
        for (int rank = 0; rank < cluster->size; ++rank) {
            for (size_t i = 0; i < row_counts[rank]; ++i)
                staged->row_offsets[row++] += preceding;
            preceding += (size_t)interval_counts[rank];
        }
        staged->row_offsets[rows] = intervals;
    }
    GATHER(row_direction_tags, local_rows, row_counts, uint8_t);
    GATHER(row_transverse_coordinates, local_rows, row_counts, double);
    GATHER(owner_offsets, local_intervals, interval_counts, size_t);
    if (cluster->rank == 0) {
        size_t interval = 0, preceding = 0;
        for (int rank = 0; rank < cluster->size; ++rank) {
            for (size_t i = 0; i < interval_counts[rank]; ++i)
                staged->owner_offsets[interval++] += preceding;
            preceding += (size_t)owner_counts[rank];
        }
        staged->owner_offsets[intervals] = owners;
    }
    GATHER(t_enter, local_intervals, interval_counts, double);
    GATHER(t_exit, local_intervals, interval_counts, double);
    GATHER(kinds, local_intervals, interval_counts, uint8_t);
    GATHER(owner_count_lower_bounds, local_intervals,
           interval_counts, size_t);
    GATHER(owner_cell_ids, local_owners, owner_counts, int);
    GATHER(owner_cell_indices, local_owners, owner_counts, int);
    GATHER(owner_material_ids, local_owners, owner_counts, int);
    GATHER(owner_universe_ids, local_owners, owner_counts, int);
    GATHER(owner_fill_universes, local_owners, owner_counts, int);
    GATHER(owner_depths, local_owners, owner_counts, int);
    GATHER(owner_occurrence_keys, local_owners, owner_counts, uint64_t);
    GATHER(owner_parent_occurrence_keys, local_owners,
           owner_counts, uint64_t);
    GATHER(owner_resolution_flags, local_owners, owner_counts, uint8_t);
#undef GATHER
    local = cluster->rank == 0 && callback(first_row, staged, user_data)
        ? ALEA_CLUSTER_INTERRUPTED : ALEA_CLUSTER_OK;
    return alea_cluster_agree(cluster, local);
}

static alea_cluster_status_t coverage_impl(
        alea_cluster_t* cluster, alea_system_t* sys,
        const double* origins_xyz, const double* directions_xyz,
        size_t row_count, const uint8_t* direction_tags,
        const double* transverse_coordinates,
        const alea_ray_coverage_slice_options_t* options,
        alea_cluster_coverage_shard_callback_t callback, void* user_data,
        int root_stream) {
    if (!cluster) return ALEA_CLUSTER_INVALID_ARGUMENT;
    alea_cluster_status_t local = sys &&
        (root_stream ? (cluster->rank != 0 || callback != NULL)
                     : callback != NULL)
        ? ALEA_CLUSTER_OK : ALEA_CLUSTER_INVALID_ARGUMENT;
    uint64_t count_wire = 0;
    int has_tags = 0, has_coordinates = 0;
    alea_ray_coverage_slice_options_t wire_options;
    memset(&wire_options, 0, sizeof(wire_options));
    if (cluster->rank == 0) {
        if (!options || options->struct_size < sizeof(*options) ||
            (options->flags & ~(ALEA_RAY_COVERAGE_DOMAIN |
                                ALEA_RAY_COVERAGE_REPORT_EXTERIOR)) ||
            !(options->t_max > 0.0) || !isfinite(options->t_max) ||
            options->max_refinement_depth ||
            ((options->flags & ALEA_RAY_COVERAGE_REPORT_EXTERIOR) &&
             !(options->flags & ALEA_RAY_COVERAGE_DOMAIN)) ||
            (row_count && (!origins_xyz || !directions_xyz)) ||
            row_count > SIZE_MAX / 3)
            local = ALEA_CLUSTER_INVALID_ARGUMENT;
        if (local == ALEA_CLUSTER_OK) {
            count_wire = row_count;
            has_tags = direction_tags != NULL;
            has_coordinates = transverse_coordinates != NULL;
            wire_options.struct_size = sizeof(wire_options);
            wire_options.flags = options->flags;
            wire_options.t_max = options->t_max;
            wire_options.domain_t_min = options->domain_t_min;
            wire_options.domain_t_max = options->domain_t_max;
            wire_options.max_rows = options->max_rows;
            wire_options.max_intervals = options->max_intervals;
            wire_options.max_owners = options->max_owners;
            wire_options.max_output_bytes = options->max_output_bytes;
        }
    }
    alea_cluster_status_t status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    if (alea_cluster_backend_broadcast_u64(cluster->backend, &count_wire, 0) ||
        alea_cluster_backend_broadcast_int(cluster->backend, &has_tags, 0) ||
        alea_cluster_backend_broadcast_int(cluster->backend,
            &has_coordinates, 0) ||
        alea_cluster_backend_broadcast_bytes(cluster->backend,
            &wire_options, sizeof(wire_options), 0)) {
        cluster->usable = 0;
        return ALEA_CLUSTER_BACKEND_ERROR;
    }
    local = count_wire > SIZE_MAX / 3
        ? ALEA_CLUSTER_INVALID_ARGUMENT : ALEA_CLUSTER_OK;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    uint64_t model = alea_cluster_system_fingerprint(sys);
    local = model ? ALEA_CLUSTER_OK : ALEA_CLUSTER_COMPUTE_ERROR;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    int match = alea_cluster_fingerprints_match(cluster, model);
    if (match < 0) return ALEA_CLUSTER_BACKEND_ERROR;
    if (!match) return ALEA_CLUSTER_MODEL_MISMATCH;

    double* local_origins = malloc(3 * COVERAGE_BATCH_ROWS * sizeof(double));
    double* local_directions = malloc(3 * COVERAGE_BATCH_ROWS * sizeof(double));
    uint8_t* local_tags = has_tags ? malloc(COVERAGE_BATCH_ROWS) : NULL;
    double* local_coordinates = has_coordinates
        ? malloc(COVERAGE_BATCH_ROWS * sizeof(double)) : NULL;
    alea_ray_coverage_slice_result_t* shard =
        alea_ray_coverage_slice_result_create();
    const size_t rank_total = (size_t)cluster->size;
    alea_ray_coverage_slice_result_t* staged =
        root_stream && cluster->rank == 0
            ? alea_ray_coverage_slice_result_create() : NULL;
    uint64_t* row_counts = root_stream && cluster->rank == 0
        ? calloc(rank_total, sizeof(uint64_t)) : NULL;
    uint64_t* interval_counts = root_stream && cluster->rank == 0
        ? calloc(rank_total, sizeof(uint64_t)) : NULL;
    uint64_t* owner_counts = root_stream && cluster->rank == 0
        ? calloc(rank_total, sizeof(uint64_t)) : NULL;
    size_t* rank_bytes = root_stream && cluster->rank == 0
        ? calloc(rank_total, sizeof(size_t)) : NULL;
    local = local_origins && local_directions && shard &&
        (!has_tags || local_tags) && (!has_coordinates || local_coordinates) &&
        (!root_stream || cluster->rank != 0 ||
         (staged && row_counts && interval_counts && owner_counts && rank_bytes))
        ? ALEA_CLUSTER_OK : ALEA_CLUSTER_OUT_OF_MEMORY;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) goto cleanup;

    const size_t ranks = (size_t)cluster->size;
    const size_t rank = (size_t)cluster->rank;
    const size_t total = (size_t)count_wire;
    for (size_t base = 0; base < total; ) {
        const size_t count = total - base < COVERAGE_BATCH_ROWS
            ? total - base : COVERAGE_BATCH_ROWS;
        const size_t mine = rank_count(count, rank, ranks);
        const size_t first = rank_first(count, rank, ranks);
        int transfer_error =
            alea_cluster_backend_scatter_blocks(cluster->backend,
                cluster->rank == 0 ? origins_xyz + 3 * base : NULL,
                local_origins, 3 * sizeof(double), count, 0) ||
            alea_cluster_backend_scatter_blocks(cluster->backend,
                cluster->rank == 0 ? directions_xyz + 3 * base : NULL,
                local_directions, 3 * sizeof(double), count, 0);
        if (!transfer_error && has_tags)
            transfer_error = alea_cluster_backend_scatter_blocks(
                cluster->backend,
                cluster->rank == 0 ? direction_tags + base : NULL,
                local_tags, sizeof(uint8_t), count, 0);
        if (!transfer_error && has_coordinates)
            transfer_error = alea_cluster_backend_scatter_blocks(
                cluster->backend,
                cluster->rank == 0 ? transverse_coordinates + base : NULL,
                local_coordinates, sizeof(double), count, 0);
        if (transfer_error) {
            cluster->usable = 0;
            status = ALEA_CLUSTER_BACKEND_ERROR;
            goto cleanup;
        }
        int rc = mine ? alea_ray_coverage_slice_query(sys,
            local_origins, local_directions, mine,
            local_tags, local_coordinates, &wire_options, shard) : 0;
        local = rc == 0 ? ALEA_CLUSTER_OK :
            alea_error_code() == ALEA_ERR_OVERFLOW
                ? ALEA_CLUSTER_OUTPUT_LIMIT : ALEA_CLUSTER_COMPUTE_ERROR;
        if (alea_interrupted()) local = ALEA_CLUSTER_INTERRUPTED;
        status = alea_cluster_agree(cluster, local);
        if (status != ALEA_CLUSTER_OK) goto cleanup;
        if (root_stream)
            status = assemble_coverage_batch(cluster, shard, base, count, mine,
                &wire_options, staged, row_counts, interval_counts,
                owner_counts, rank_bytes, callback, user_data);
        else {
            local = mine && callback(base + first, shard, user_data)
                ? ALEA_CLUSTER_INTERRUPTED : ALEA_CLUSTER_OK;
            status = alea_cluster_agree(cluster, local);
        }
        if (status != ALEA_CLUSTER_OK) goto cleanup;
        base += count;
    }
    status = ALEA_CLUSTER_OK;

cleanup:
    alea_ray_coverage_slice_result_destroy(shard);
    alea_ray_coverage_slice_result_destroy(staged);
    free(local_origins); free(local_directions);
    free(local_tags); free(local_coordinates);
    free(row_counts); free(interval_counts); free(owner_counts); free(rank_bytes);
    return status;
}

alea_cluster_status_t alea_cluster_coverage_shards(
        alea_cluster_t* cluster, alea_system_t* sys,
        const double* origins_xyz, const double* directions_xyz,
        size_t row_count, const uint8_t* direction_tags,
        const double* transverse_coordinates,
        const alea_ray_coverage_slice_options_t* options,
        alea_cluster_coverage_shard_callback_t callback, void* user_data) {
    return coverage_impl(cluster, sys, origins_xyz, directions_xyz,
        row_count, direction_tags, transverse_coordinates, options,
        callback, user_data, 0);
}

alea_cluster_status_t alea_cluster_coverage_stream(
        alea_cluster_t* cluster, alea_system_t* sys,
        const double* origins_xyz, const double* directions_xyz,
        size_t row_count, const uint8_t* direction_tags,
        const double* transverse_coordinates,
        const alea_ray_coverage_slice_options_t* options,
        alea_cluster_coverage_shard_callback_t callback, void* user_data) {
    return coverage_impl(cluster, sys, origins_xyz, directions_xyz,
        row_count, direction_tags, transverse_coordinates, options,
        callback, user_data, 1);
}
