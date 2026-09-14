// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_cluster.h"
#include "cluster_internal.h"
#include "raycast/volume_internal.h"

#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct alea_cluster {
    void* backend;
    int rank;
    int size;
    int usable;
};

static int g_initialized;
static int g_context_count;

#ifndef ALEA_CLUSTER_USE_MPI
int alea_cluster_backend_initialize(int* argc, char*** argv) {
    (void)argc; (void)argv; return 0;
}
int alea_cluster_backend_finalize(void) { return 0; }
int alea_cluster_backend_create(void** state, int* rank, int* size) {
    if (!state || !rank || !size) return -1;
    *state = (void*)1; *rank = 0; *size = 1; return 0;
}
void alea_cluster_backend_destroy(void* state) { (void)state; }
const char* alea_cluster_backend_name(void) { return "local"; }
int alea_cluster_backend_agree_status(void* state, int local, int* global) {
    if (!state || !global) return -1;
    *global = local;
    return 0;
}
int alea_cluster_backend_u64_minmax(void* state, uint64_t local,
                                    uint64_t* minimum, uint64_t* maximum) {
    if (!state || !minimum || !maximum) return -1;
    *minimum = local; *maximum = local; return 0;
}
int alea_cluster_backend_sum_doubles(void* state, double* values,
                                     size_t count) {
    return state && (values || count == 0) ? 0 : -1;
}
int alea_cluster_backend_broadcast_int(void* state, int* value, int root) {
    return state && value && root == 0 ? 0 : -1;
}
#endif

static uint64_t hash_bytes(uint64_t hash, const void* data, size_t size) {
    const unsigned char* bytes = data;
    for (size_t i = 0; i < size; ++i) {
        hash ^= bytes[i];
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

#define HASH_FIELD(hash, value) hash_bytes((hash), &(value), sizeof(value))

static uint64_t options_fingerprint(
        const alea_volume_estimate_options_t* options) {
    uint64_t hash = UINT64_C(1469598103934665603);
    hash = HASH_FIELD(hash, options->max_rays);
    hash = HASH_FIELD(hash, options->seed);
    hash = HASH_FIELD(hash, options->rng_algorithm);
    hash = HASH_FIELD(hash, options->requested_workers);
    hash = HASH_FIELD(hash, options->batch_size);
    hash = HASH_FIELD(hash, options->target_rel_error);
    hash = HASH_FIELD(hash, options->use_sampling_sphere);
    hash = hash_bytes(hash, options->sampling_center,
                      sizeof(options->sampling_center));
    hash = HASH_FIELD(hash, options->sampling_radius);
    return hash;
}

static uint64_t problem_fingerprint(alea_system_t* sys,
                                    const alea_volume_problem_t* problem) {
    uint64_t hash = UINT64_C(1469598103934665603);
    hash = HASH_FIELD(hash, problem->path_count);
    hash = HASH_FIELD(hash, problem->cx);
    hash = HASH_FIELD(hash, problem->cy);
    hash = HASH_FIELD(hash, problem->cz);
    hash = HASH_FIELD(hash, problem->radius);
    if (problem->path_count == 0) return hash;
    alea_volume_path_t* paths = calloc(problem->path_count, sizeof(*paths));
    if (!paths) return 0;
    const size_t got = alea_volume_paths_get(sys, paths, problem->path_count);
    if (got != problem->path_count) { free(paths); return 0; }
    for (size_t p = 0; p < got; ++p) {
        const alea_volume_path_t* path = &paths[p];
        hash = HASH_FIELD(hash, path->path_id);
        hash = HASH_FIELD(hash, path->terminal_cell_index);
        hash = HASH_FIELD(hash, path->terminal_cell_id);
        hash = HASH_FIELD(hash, path->material_id);
        hash = HASH_FIELD(hash, path->universe_id);
        hash = HASH_FIELD(hash, path->depth);
        hash = HASH_FIELD(hash, path->ancestor_count);
        hash = HASH_FIELD(hash, path->lattice_step_count);
        for (size_t i = 0; i < path->ancestor_count; ++i) {
            hash = HASH_FIELD(hash, path->ancestor_cell_indices[i]);
            hash = HASH_FIELD(hash, path->ancestor_universe_ids[i]);
        }
        for (size_t i = 0; i < path->lattice_step_count; ++i) {
            const alea_volume_lattice_step_t* step = &path->lattice_steps[i];
            hash = HASH_FIELD(hash, step->lattice_cell_index);
            hash = HASH_FIELD(hash, step->fill_universe);
            hash = HASH_FIELD(hash, step->i);
            hash = HASH_FIELD(hash, step->j);
            hash = HASH_FIELD(hash, step->k);
            hash = HASH_FIELD(hash, step->linear_index);
        }
        hash = hash_bytes(hash, path->world_to_local,
                          sizeof(path->world_to_local));
    }
    free(paths);
    return hash;
}

static alea_cluster_status_t agree(alea_cluster_t* cluster,
                                   alea_cluster_status_t local) {
    int global = (int)local;
    if (!cluster || !cluster->usable || alea_cluster_backend_agree_status(
            cluster->backend, (int)local, &global) != 0) {
        if (cluster) cluster->usable = 0;
        return ALEA_CLUSTER_BACKEND_ERROR;
    }
    return (alea_cluster_status_t)global;
}

alea_cluster_status_t alea_cluster_initialize(int* argc, char*** argv) {
    if (g_initialized) return ALEA_CLUSTER_OK;
    if (alea_cluster_backend_initialize(argc, argv) != 0)
        return ALEA_CLUSTER_BACKEND_ERROR;
    g_initialized = 1;
    return ALEA_CLUSTER_OK;
}

alea_cluster_status_t alea_cluster_finalize(void) {
    if (!g_initialized) return ALEA_CLUSTER_INVALID_STATE;
    if (g_context_count != 0) return ALEA_CLUSTER_INVALID_STATE;
    if (alea_cluster_backend_finalize() != 0)
        return ALEA_CLUSTER_BACKEND_ERROR;
    g_initialized = 0;
    return ALEA_CLUSTER_OK;
}

alea_cluster_t* alea_cluster_create(void) {
    if (!g_initialized || g_context_count != 0) return NULL;
    alea_cluster_t* cluster = calloc(1, sizeof(*cluster));
    if (!cluster) return NULL;
    if (alea_cluster_backend_create(&cluster->backend, &cluster->rank,
                                    &cluster->size) != 0) {
        free(cluster); return NULL;
    }
    cluster->usable = 1;
    g_context_count = 1;
    return cluster;
}

void alea_cluster_destroy(alea_cluster_t* cluster) {
    if (!cluster) return;
    alea_cluster_backend_destroy(cluster->backend);
    free(cluster);
    g_context_count = 0;
}

int alea_cluster_rank(const alea_cluster_t* c) { return c ? c->rank : -1; }
int alea_cluster_size(const alea_cluster_t* c) { return c ? c->size : 0; }
int alea_cluster_is_root(const alea_cluster_t* c) {
    return c && c->rank == 0;
}
const char* alea_cluster_backend(const alea_cluster_t* c) {
    return c ? alea_cluster_backend_name() : NULL;
}

static int fingerprints_match(alea_cluster_t* cluster, uint64_t fingerprint) {
    uint64_t minimum = 0, maximum = 0;
    if (alea_cluster_backend_u64_minmax(cluster->backend, fingerprint,
                                       &minimum, &maximum) != 0) {
        cluster->usable = 0; return -1;
    }
    return minimum == maximum ? 1 : 0;
}

alea_cluster_status_t alea_cluster_estimate_volumes(
        alea_cluster_t* cluster, alea_system_t* sys,
        const alea_volume_estimate_options_t* options, double* volumes,
        double* rel_errors, alea_cluster_volume_stats_t* out_stats) {
    alea_cluster_status_t local = ALEA_CLUSTER_OK;
    if (!cluster || !cluster->usable || !sys || !options || !volumes ||
        options->max_rays == 0 ||
        options->requested_workers > (size_t)INT_MAX ||
        options->rng_algorithm != ALEA_RNG_PHILOX4X32_10 ||
        !isfinite(options->target_rel_error) ||
        options->target_rel_error < 0.0 || options->target_rel_error > 1.0)
        local = ALEA_CLUSTER_INVALID_ARGUMENT;
    if (local == ALEA_CLUSTER_OK && options->use_sampling_sphere &&
        (!isfinite(options->sampling_center[0]) ||
         !isfinite(options->sampling_center[1]) ||
         !isfinite(options->sampling_center[2]) ||
         !isfinite(options->sampling_radius) || options->sampling_radius <= 0.0))
        local = ALEA_CLUSTER_INVALID_ARGUMENT;
    alea_cluster_status_t status = agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;

    int matching = fingerprints_match(cluster, options_fingerprint(options));
    if (matching < 0) return ALEA_CLUSTER_BACKEND_ERROR;
    if (!matching) return ALEA_CLUSTER_INVALID_ARGUMENT;

    alea_volume_problem_t problem;
    local = alea_volume_problem_prepare(sys, options, &problem) == 0
        ? ALEA_CLUSTER_OK : ALEA_CLUSTER_COMPUTE_ERROR;
    status = agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    uint64_t fingerprint = problem_fingerprint(sys, &problem);
    local = fingerprint ? ALEA_CLUSTER_OK : ALEA_CLUSTER_OUT_OF_MEMORY;
    status = agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) return status;
    matching = fingerprints_match(cluster, fingerprint);
    if (matching < 0) return ALEA_CLUSTER_BACKEND_ERROR;
    if (!matching) return ALEA_CLUSTER_MODEL_MISMATCH;

    const size_t count = problem.path_count;
    double* local_l = count ? calloc(count, sizeof(*local_l)) : NULL;
    double* local_l2 = count ? calloc(count, sizeof(*local_l2)) : NULL;
    double* sum_l2 = count ? calloc(count, sizeof(*sum_l2)) : NULL;
    double* errors = rel_errors;
    if (count && !errors) errors = calloc(count, sizeof(*errors));
    local = count && (!local_l || !local_l2 || !sum_l2 || !errors)
        ? ALEA_CLUSTER_OUT_OF_MEMORY : ALEA_CLUSTER_OK;
    status = agree(cluster, local);
    if (status != ALEA_CLUSTER_OK) goto cleanup;

    if (count) memset(volumes, 0, count * sizeof(*volumes));
    if (rel_errors && count) memset(rel_errors, 0, count * sizeof(*rel_errors));
    if (out_stats) memset(out_stats, 0, sizeof(*out_stats));
    if (count == 0) goto cleanup;

    size_t batch_size = options->batch_size ? options->batch_size : 10000;
    if (batch_size > options->max_rays) batch_size = options->max_rays;
    size_t completed = 0;
    size_t local_completed = 0;
    size_t local_workers = 0;
    double maximum_error = INFINITY;
    int converged = 0, cancelled = 0, interrupted = 0;

    while (completed < options->max_rays) {
        if (alea_interrupted()) {
            local = ALEA_CLUSTER_INTERRUPTED;
        } else {
            local = ALEA_CLUSTER_OK;
        }
        status = agree(cluster, local);
        if (status == ALEA_CLUSTER_INTERRUPTED) { interrupted = 1; break; }
        if (status != ALEA_CLUSTER_OK) goto cleanup;

        size_t batch_end = completed + batch_size;
        if (batch_end < completed || batch_end > options->max_rays)
            batch_end = options->max_rays;
        const size_t rays = batch_end - completed;
        const size_t ranks = (size_t)cluster->size;
        const size_t rank = (size_t)cluster->rank;
        const size_t base = rays / ranks;
        const size_t remainder = rays % ranks;
        const size_t local_count = base + (rank < remainder ? 1u : 0u);
        const size_t local_begin = completed + rank * base +
            (rank < remainder ? rank : remainder);

        memset(local_l, 0, count * sizeof(*local_l));
        memset(local_l2, 0, count * sizeof(*local_l2));
        size_t actual_workers = 0;
        local = alea_volume_accumulate_ray_range(
            sys, &problem, local_begin, local_begin + local_count,
            options->rng_algorithm, options->seed, options->requested_workers,
            local_l, local_l2, &actual_workers) == 0
            ? ALEA_CLUSTER_OK : ALEA_CLUSTER_COMPUTE_ERROR;
        if (alea_interrupted()) local = ALEA_CLUSTER_INTERRUPTED;
        status = agree(cluster, local);
        if (status == ALEA_CLUSTER_INTERRUPTED) { interrupted = 1; break; }
        if (status != ALEA_CLUSTER_OK) goto cleanup;

        if (alea_cluster_backend_sum_doubles(cluster->backend, local_l, count) ||
            alea_cluster_backend_sum_doubles(cluster->backend, local_l2, count)) {
            cluster->usable = 0; status = ALEA_CLUSTER_BACKEND_ERROR; goto cleanup;
        }
        for (size_t i = 0; i < count; ++i) {
            volumes[i] += local_l[i];
            sum_l2[i] += local_l2[i];
        }
        completed = batch_end;
        local_completed += local_count;
        if (actual_workers > local_workers) local_workers = actual_workers;
        alea_volume_compute_errors(volumes, sum_l2, errors, count, completed);
        maximum_error = 0.0;
        for (size_t i = 0; i < count; ++i) {
            if (errors[i] < 0.0 || !isfinite(errors[i])) {
                maximum_error = INFINITY; break;
            }
            if (errors[i] > maximum_error) maximum_error = errors[i];
        }
        int stop = 0;
        if (cluster->rank == 0) {
            if (options->target_rel_error > 0.0 &&
                maximum_error <= options->target_rel_error) {
                converged = 1; stop = 1;
            }
            if (options->progress && options->progress(
                    completed, options->max_rays, maximum_error,
                    options->progress_user_data)) {
                cancelled = 1; stop = 1;
            }
        }
        if (alea_cluster_backend_broadcast_int(cluster->backend, &stop, 0) ||
            alea_cluster_backend_broadcast_int(cluster->backend, &converged, 0) ||
            alea_cluster_backend_broadcast_int(cluster->backend, &cancelled, 0)) {
            cluster->usable = 0; status = ALEA_CLUSTER_BACKEND_ERROR; goto cleanup;
        }
        if (stop) break;
    }

    if (completed) {
        const double scale = M_PI * problem.radius * problem.radius /
                             (double)completed;
        for (size_t i = 0; i < count; ++i) volumes[i] *= scale;
    }
    if (out_stats) {
        out_stats->rank_count = cluster->size;
        out_stats->local_rays_completed = local_completed;
        out_stats->local_workers = local_workers;
        out_stats->volume.rays_completed = completed;
        out_stats->volume.requested_workers = options->requested_workers;
        out_stats->volume.actual_workers = local_workers;
        out_stats->volume.batch_size = batch_size;
        out_stats->volume.rng_algorithm = options->rng_algorithm;
        out_stats->volume.rng_address_version = ALEA_RNG_ADDRESS_VERSION;
        out_stats->volume.seed = options->seed;
        out_stats->volume.maximum_relative_error = maximum_error;
        out_stats->volume.converged = converged != 0;
        out_stats->volume.cancelled = (cancelled || interrupted) != 0;
    }
    status = interrupted ? ALEA_CLUSTER_INTERRUPTED : ALEA_CLUSTER_OK;

cleanup:
    free(local_l); free(local_l2); free(sum_l2);
    if (errors != rel_errors) free(errors);
    return status;
}

const char* alea_cluster_status_string(alea_cluster_status_t status) {
    switch (status) {
        case ALEA_CLUSTER_OK: return "success";
        case ALEA_CLUSTER_INVALID_ARGUMENT: return "invalid argument";
        case ALEA_CLUSTER_INVALID_STATE: return "invalid state";
        case ALEA_CLUSTER_OUT_OF_MEMORY: return "out of memory";
        case ALEA_CLUSTER_BACKEND_ERROR: return "cluster backend error";
        case ALEA_CLUSTER_MODEL_MISMATCH: return "model mismatch between ranks";
        case ALEA_CLUSTER_COMPUTE_ERROR: return "local computation failed";
        case ALEA_CLUSTER_INTERRUPTED: return "operation interrupted";
        default: return "unknown cluster status";
    }
}
