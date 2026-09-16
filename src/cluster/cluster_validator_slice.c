// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "cluster_internal.h"
#include "geo_validator/geo_validator.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define VALIDATOR_BATCH_CURVES 64u
#define VALIDATOR_CURVE_ERROR_LIMIT 4096u
#define VALIDATOR_META_FIELDS 8u

static size_t rank_count(size_t total, size_t rank, size_t ranks) {
    return total / ranks + (rank < total % ranks);
}

static size_t rank_first(size_t total, size_t rank, size_t ranks) {
    return rank * (total / ranks) + (rank < total % ranks ? rank : total % ranks);
}

static uint64_t slice_fingerprint(const alea_slice_view_t* view,
                                  const alea_slice_curves_t* curves) {
    uint64_t hash = UINT64_C(1469598103934665603);
#define HASH(value) do { \
    const unsigned char* bytes = (const unsigned char*)&(value); \
    for (size_t i = 0; i < sizeof(value); ++i) { \
        hash ^= bytes[i]; hash *= UINT64_C(1099511628211); \
    } \
} while (0)
    for (size_t i = 0; i < 3; ++i) {
        HASH(view->plane.origin[i]); HASH(view->plane.normal[i]);
        HASH(view->plane.u_axis[i]); HASH(view->plane.v_axis[i]);
    }
    HASH(view->u_min); HASH(view->u_max);
    HASH(view->v_min); HASH(view->v_max);
    size_t count = alea_slice_curves_count(curves);
    HASH(count);
    for (size_t i = 0; i < count; ++i) {
        alea_curve_t curve;
        if (alea_slice_curves_get(curves, i, &curve) != 0) return 0;
        HASH(curve.type); HASH(curve.surface_id);
        HASH(curve.primitive_id); HASH(curve.t_min); HASH(curve.t_max);
        switch (curve.type) {
        case ALEA_CURVE_LINE:
        case ALEA_CURVE_LINE_SEGMENT:
        case ALEA_CURVE_RAY:
            for (size_t j = 0; j < 2; ++j) {
                HASH(curve.data.line.point[j]);
                HASH(curve.data.line.direction[j]);
            }
            break;
        case ALEA_CURVE_CIRCLE:
        case ALEA_CURVE_ARC:
            for (size_t j = 0; j < 2; ++j)
                HASH(curve.data.circle.center[j]);
            HASH(curve.data.circle.radius);
            break;
        case ALEA_CURVE_ELLIPSE:
        case ALEA_CURVE_ELLIPSE_ARC:
            for (size_t j = 0; j < 2; ++j)
                HASH(curve.data.ellipse.center[j]);
            HASH(curve.data.ellipse.semi_a);
            HASH(curve.data.ellipse.semi_b);
            HASH(curve.data.ellipse.angle);
            break;
        case ALEA_CURVE_POLYGON:
            HASH(curve.data.polygon.count);
            HASH(curve.data.polygon.closed);
            if (curve.data.polygon.count < 0 ||
                curve.data.polygon.count > 16) return 0;
            for (int j = 0; j < curve.data.polygon.count; ++j) {
                HASH(curve.data.polygon.vertices[j][0]);
                HASH(curve.data.polygon.vertices[j][1]);
            }
            break;
        case ALEA_CURVE_PARALLEL_LINES:
            for (size_t j = 0; j < 2; ++j) {
                HASH(curve.data.parallel_lines.point1[j]);
                HASH(curve.data.parallel_lines.point2[j]);
                HASH(curve.data.parallel_lines.direction[j]);
            }
            break;
        default:
            break;
        }
    }
#undef HASH
    return hash ? hash : 1;
}

alea_cluster_status_t alea_cluster_validate_slice_curves(
        alea_cluster_t* cluster, alea_system_t* sys,
        const alea_slice_view_t* view,
        const alea_slice_curves_t* curves,
        const alea_geom_validator_options_t* options,
        alea_geom_validator_result_t* root_result) {
    if (!cluster) return ALEA_CLUSTER_INVALID_ARGUMENT;
    alea_cluster_status_t local = sys && view && curves &&
        isfinite(view->u_min) && isfinite(view->u_max) &&
        isfinite(view->v_min) && isfinite(view->v_max) &&
        view->u_min < view->u_max && view->v_min < view->v_max
        ? ALEA_CLUSTER_OK : ALEA_CLUSTER_INVALID_ARGUMENT;
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
    uint64_t slice = slice_fingerprint(view, curves);
    local = slice ? ALEA_CLUSTER_OK : ALEA_CLUSTER_COMPUTE_ERROR;
    status = alea_cluster_agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    match = alea_cluster_fingerprints_match(cluster, slice);
    if (match < 0) return ALEA_CLUSTER_BACKEND_ERROR;
    if (!match) return ALEA_CLUSTER_INVALID_ARGUMENT;
    int skip = cluster->rank == 0 && root_result->truncated;
    if (alea_cluster_backend_broadcast_int(cluster->backend, &skip, 0)) {
        cluster->usable = 0;
        return ALEA_CLUSTER_BACKEND_ERROR;
    }
    if (skip)
        return ALEA_CLUSTER_OK;

    const size_t ranks = (size_t)cluster->size;
    uint64_t local_meta[VALIDATOR_BATCH_CURVES * VALIDATOR_META_FIELDS];
    alea_geom_validator_result_t per_curve[VALIDATOR_BATCH_CURVES];
    uint64_t* all_meta = cluster->rank == 0
        ? malloc(VALIDATOR_BATCH_CURVES * VALIDATOR_META_FIELDS *
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
    worker_options.max_errors = VALIDATOR_CURVE_ERROR_LIMIT;
    worker_options.max_samples_per_signature = 0;
    worker_options.max_crossings = SIZE_MAX;
    const size_t total = alea_slice_curves_count(curves);
    for (size_t base = 0; base < total; ) {
        size_t count = total - base < VALIDATOR_BATCH_CURVES
            ? total - base : VALIDATOR_BATCH_CURVES;
        for (size_t i = 0; i < count; ++i)
            alea_geom_validator_result_init(&per_curve[i]);
        size_t mine = rank_count(count, (size_t)cluster->rank, ranks);
        size_t first = rank_first(count, (size_t)cluster->rank, ranks);
        local = ALEA_CLUSTER_OK;
        size_t local_errors = 0;
        for (size_t i = 0; i < mine; ++i) {
            alea_geom_validator_result_t* current = &per_curve[i];
            if (alea_validator_cluster_slice_range(sys, view, curves,
                    &worker_options, current, base + first + i,
                    base + first + i + 1) != 0) {
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
            size_t n = per_curve[i].error_count;
            if (n) memcpy(local_errors_flat + copied, per_curve[i].errors,
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
                if (alea_validator_cluster_merge_curve_one(sys, view, curves,
                        &prepared, root_result, &candidate, base + i) != 0) {
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
            alea_geom_validator_result_free(&per_curve[i]);
        if (stopped) break;
        base += count;
        continue;

free_batch:
        for (size_t i = 0; i < count; ++i)
            alea_geom_validator_result_free(&per_curve[i]);
        goto cleanup;
    }
    status = ALEA_CLUSTER_OK;

cleanup:
    free(all_meta); free(error_counts); free(rank_bytes);
    return status;
}
