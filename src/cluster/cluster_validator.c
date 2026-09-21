// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "cluster_internal.h"
#include "geo_validator/geo_validator.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define VALIDATOR_BATCH_RAYS 64u
#define VALIDATOR_RAY_ERROR_LIMIT 4096u
#define VALIDATOR_META_FIELDS 10u

static size_t rank_count(size_t total, size_t rank, size_t ranks) {
    return total / ranks + (rank < total % ranks);
}

static size_t rank_first(size_t total, size_t rank, size_t ranks) {
    return rank * (total / ranks) + (rank < total % ranks ? rank : total % ranks);
}

uint64_t alea_cluster_validator_options_fingerprint(
        const alea_geom_validator_options_t* options) {
    uint64_t hash = UINT64_C(1469598103934665603);
#define HASH(value) do { \
    const unsigned char* bytes = (const unsigned char*)&(value); \
    for (size_t i = 0; i < sizeof(value); ++i) { \
        hash ^= bytes[i]; hash *= UINT64_C(1099511628211); \
    } \
} while (0)
    HASH(options->flags); HASH(options->universe_depth);
    HASH(options->max_errors); HASH(options->max_samples_per_signature);
    HASH(options->max_samples_per_curve); HASH(options->max_crossings);
    HASH(options->max_breakpoints);
    HASH(options->sample_offset); HASH(options->t_max);
    HASH(options->seed); HASH(options->ray_count);
    for (size_t i = 0; i < 6; ++i) HASH(options->validation_bounds[i]);
#undef HASH
    return hash ? hash : 1;
}

alea_cluster_status_t alea_cluster_validate_geometry(
        alea_cluster_t* cluster, alea_system_t* sys,
        const alea_geom_validator_options_t* options,
        alea_geom_validator_result_t* root_result) {
    if (!cluster) return ALEA_CLUSTER_INVALID_ARGUMENT;
    alea_cluster_status_t local = sys ? ALEA_CLUSTER_OK
                                      : ALEA_CLUSTER_INVALID_ARGUMENT;
    if (cluster->rank == 0 && !root_result)
        local = ALEA_CLUSTER_INVALID_ARGUMENT;
    alea_cluster_status_t status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    alea_geom_validator_options_t prepared;
    double bounds[6], t_max;
    uint64_t rng;
    local = alea_validator_cluster_source_init(sys, options, &prepared,
        bounds, &t_max, &rng) == 0 ? ALEA_CLUSTER_OK
                                : ALEA_CLUSTER_COMPUTE_ERROR;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    uint64_t model = alea_cluster_system_fingerprint(sys);
    local = model ? ALEA_CLUSTER_OK : ALEA_CLUSTER_COMPUTE_ERROR;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    int match = alea_cluster_fingerprints_match(cluster, model);
    if (match < 0) return ALEA_CLUSTER_BACKEND_ERROR;
    if (!match) return ALEA_CLUSTER_MODEL_MISMATCH;
    match = alea_cluster_fingerprints_match(cluster,
        alea_cluster_validator_options_fingerprint(&prepared));
    if (match < 0) return ALEA_CLUSTER_BACKEND_ERROR;
    if (!match) return ALEA_CLUSTER_INVALID_ARGUMENT;
    int skip = cluster->rank == 0 && root_result->truncated;
    if (alea_cluster_backend_broadcast_int(cluster->backend, &skip, 0)) {
        cluster->usable = 0;
        return ALEA_CLUSTER_BACKEND_ERROR;
    }
    if (!(prepared.flags & ALEA_GEOM_VALIDATE_RAYS) || skip)
        return ALEA_CLUSTER_OK;

    const size_t ranks = (size_t)cluster->size;
    double rays[VALIDATOR_BATCH_RAYS * 6];
    uint64_t local_meta[VALIDATOR_BATCH_RAYS * VALIDATOR_META_FIELDS];
    alea_geom_validator_result_t per_ray[VALIDATOR_BATCH_RAYS];
    uint64_t* all_meta = cluster->rank == 0
        ? malloc(VALIDATOR_BATCH_RAYS * VALIDATOR_META_FIELDS *
                 sizeof(uint64_t)) : NULL;
    uint64_t* error_counts = cluster->rank == 0
        ? malloc(ranks * sizeof(uint64_t)) : NULL;
    size_t* rank_bytes = cluster->rank == 0
        ? malloc(ranks * sizeof(size_t)) : NULL;
    local = cluster->rank != 0 || (all_meta && error_counts && rank_bytes)
        ? ALEA_CLUSTER_OK : ALEA_CLUSTER_OUT_OF_MEMORY;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) goto cleanup;
    alea_geom_validator_options_t worker_options = prepared;
    worker_options.max_errors = VALIDATOR_RAY_ERROR_LIMIT;
    worker_options.max_samples_per_signature = 0;
    worker_options.max_crossings = SIZE_MAX;
    const size_t total = (size_t)prepared.ray_count;
    for (size_t base = 0; base < total; ) {
        size_t count = total - base < VALIDATOR_BATCH_RAYS
            ? total - base : VALIDATOR_BATCH_RAYS;
        if (cluster->rank == 0)
            for (size_t i = 0; i < count; ++i)
                alea_validator_cluster_next_ray(&rng, bounds,
                    rays + 6 * i, rays + 6 * i + 3);
        if (alea_cluster_backend_broadcast_bytes(cluster->backend, rays,
                count * 6 * sizeof(double), 0)) {
            cluster->usable = 0;
            status = ALEA_CLUSTER_BACKEND_ERROR;
            goto cleanup;
        }
        for (size_t i = 0; i < count; ++i)
            alea_geom_validator_result_init(&per_ray[i]);
        size_t mine = rank_count(count, (size_t)cluster->rank, ranks);
        size_t first = rank_first(count, (size_t)cluster->rank, ranks);
        local = ALEA_CLUSTER_OK;
        size_t local_errors = 0;
        for (size_t i = 0; i < mine; ++i) {
            const double* ray = rays + 6 * (first + i);
            alea_geom_validator_result_t* current = &per_ray[i];
            if (alea_validate_geometry_ray(sys, &worker_options,
                    ray[0], ray[1], ray[2], ray[3], ray[4], ray[5],
                    t_max, current) != 0) {
                local = ALEA_CLUSTER_COMPUTE_ERROR;
                break;
            }
            if (current->truncated ||
                current->error_count > SIZE_MAX - local_errors) {
                local = ALEA_CLUSTER_OUTPUT_LIMIT;
                break;
            }
            local_errors += current->error_count;
            uint64_t* meta = local_meta + i * VALIDATOR_META_FIELDS;
            meta[0] = current->error_count;
            meta[1] = current->crossings_checked;
            meta[2] = current->adjacency_hits;
            meta[3] = current->exact_queries;
            meta[4] = current->ambiguous_crossings;
            meta[5] = current->suppressed_samples;
            meta[6] = current->sample_limited_curves;
            meta[7] = current->truncated;
            meta[8] = current->incomplete_rays;
            meta[9] = current->incomplete_slice_samples;
        }
        if (alea_interrupted()) local = ALEA_CLUSTER_INTERRUPTED;
        status = alea_cluster_agree(cluster, local);
        if (status != ALEA_CLUSTER_OK) goto free_batch;
        alea_geom_error_t* local_errors_flat = malloc(
            (local_errors ? local_errors : 1) * sizeof(alea_geom_error_t));
        local = local_errors_flat ? ALEA_CLUSTER_OK
                                  : ALEA_CLUSTER_OUT_OF_MEMORY;
        status = alea_cluster_agree(cluster, local);
        if (status != ALEA_CLUSTER_OK) { free(local_errors_flat); goto free_batch; }
        size_t copied = 0;
        for (size_t i = 0; i < mine; ++i) {
            size_t n = per_ray[i].error_count;
            if (n) memcpy(local_errors_flat + copied, per_ray[i].errors,
                          n * sizeof(alea_geom_error_t));
            copied += n;
        }
        if (cluster->rank == 0)
            for (size_t rank = 0; rank < ranks; ++rank)
                rank_bytes[rank] = rank_count(count, rank, ranks) *
                    VALIDATOR_META_FIELDS * sizeof(uint64_t);
        if (alea_cluster_backend_gather_bytes(cluster->backend, local_meta,
                mine * VALIDATOR_META_FIELDS * sizeof(uint64_t),
                all_meta, rank_bytes, 0) ||
            alea_cluster_backend_gather_u64(cluster->backend,
                (uint64_t)local_errors, error_counts, 0)) {
            cluster->usable = 0;
            status = ALEA_CLUSTER_BACKEND_ERROR;
            free(local_errors_flat);
            goto free_batch;
        }
        size_t total_errors = 0;
        local = ALEA_CLUSTER_OK;
        if (cluster->rank == 0)
            for (size_t rank = 0; rank < ranks; ++rank) {
                if (error_counts[rank] > SIZE_MAX - total_errors ||
                    error_counts[rank] > SIZE_MAX / sizeof(alea_geom_error_t)) {
                    local = ALEA_CLUSTER_OUTPUT_LIMIT;
                    break;
                }
                total_errors += (size_t)error_counts[rank];
            }
        alea_geom_error_t* all_errors = NULL;
        if (cluster->rank == 0 && local == ALEA_CLUSTER_OK) {
            if (total_errors > SIZE_MAX / sizeof(alea_geom_error_t))
                local = ALEA_CLUSTER_OUTPUT_LIMIT;
            else {
                all_errors = malloc((total_errors ? total_errors : 1) *
                                    sizeof(alea_geom_error_t));
                if (!all_errors) local = ALEA_CLUSTER_OUT_OF_MEMORY;
            }
        }
        status = alea_cluster_agree(cluster, local);
        if (status != ALEA_CLUSTER_OK) {
            free(all_errors); free(local_errors_flat); goto free_batch;
        }
        if (cluster->rank == 0)
            for (size_t rank = 0; rank < ranks; ++rank)
                rank_bytes[rank] = (size_t)error_counts[rank] *
                    sizeof(alea_geom_error_t);
        if (alea_cluster_backend_gather_bytes(cluster->backend,
                local_errors_flat, local_errors * sizeof(alea_geom_error_t),
                all_errors, rank_bytes, 0)) {
            cluster->usable = 0;
            status = ALEA_CLUSTER_BACKEND_ERROR;
            free(all_errors); free(local_errors_flat); goto free_batch;
        }
        free(local_errors_flat);
        local = ALEA_CLUSTER_OK;
        if (cluster->rank == 0) {
            size_t error_offset = 0;
            for (size_t i = 0; i < count && !root_result->truncated; ++i) {
                const uint64_t* meta = all_meta + i * VALIDATOR_META_FIELDS;
                if (meta[0] > total_errors - error_offset) {
                    local = ALEA_CLUSTER_COMPUTE_ERROR;
                    break;
                }
                alea_geom_validator_result_t candidate = {0};
                candidate.errors = all_errors + error_offset;
                candidate.error_count = (size_t)meta[0];
                candidate.crossings_checked = (size_t)meta[1];
                candidate.adjacency_hits = (size_t)meta[2];
                candidate.exact_queries = (size_t)meta[3];
                candidate.ambiguous_crossings = (size_t)meta[4];
                candidate.suppressed_samples = (size_t)meta[5];
                candidate.sample_limited_curves = (size_t)meta[6];
                candidate.truncated = (int)meta[7];
                candidate.incomplete_rays = (size_t)meta[8];
                candidate.incomplete_slice_samples = (size_t)meta[9];
                const double* ray = rays + 6 * i;
                if (alea_validator_cluster_merge_one(sys, ray, ray + 3,
                        t_max, &prepared, root_result, &candidate) != 0) {
                    local = ALEA_CLUSTER_COMPUTE_ERROR;
                    break;
                }
                error_offset += candidate.error_count;
            }
        }
        free(all_errors);
        status = alea_cluster_agree(cluster, local);
        if (status != ALEA_CLUSTER_OK) goto free_batch;
        int stopped = cluster->rank == 0 && root_result->truncated;
        if (alea_cluster_backend_broadcast_int(cluster->backend,
                                               &stopped, 0)) {
            cluster->usable = 0;
            status = ALEA_CLUSTER_BACKEND_ERROR;
            goto free_batch;
        }
        for (size_t i = 0; i < count; ++i)
            alea_geom_validator_result_free(&per_ray[i]);
        if (stopped) break;
        base += count;
        continue;

free_batch:
        for (size_t i = 0; i < count; ++i)
            alea_geom_validator_result_free(&per_ray[i]);
        goto cleanup;
    }
    status = ALEA_CLUSTER_OK;

cleanup:
    free(all_meta); free(error_counts); free(rank_bytes);
    return status;
}
