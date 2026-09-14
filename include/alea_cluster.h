// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_CLUSTER_H
#define ALEA_CLUSTER_H

#include "alea.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct alea_cluster alea_cluster_t;

typedef enum {
    ALEA_CLUSTER_OK = 0,
    ALEA_CLUSTER_INVALID_ARGUMENT,
    ALEA_CLUSTER_INVALID_STATE,
    ALEA_CLUSTER_OUT_OF_MEMORY,
    ALEA_CLUSTER_BACKEND_ERROR,
    ALEA_CLUSTER_MODEL_MISMATCH,
    ALEA_CLUSTER_COMPUTE_ERROR,
    ALEA_CLUSTER_INTERRUPTED
} alea_cluster_status_t;

typedef struct {
    alea_volume_estimate_stats_t volume;
    int rank_count;
    size_t local_rays_completed;
    size_t local_workers;
} alea_cluster_volume_stats_t;

/**
 * Initialize the process-wide cluster runtime.
 *
 * In an MPI build this collectively attaches to MPI_COMM_WORLD, initializing
 * MPI with MPI_THREAD_FUNNELED when needed. MPI calls must be made by the
 * thread that calls this function. A local build creates a one-rank runtime.
 */
alea_cluster_status_t alea_cluster_initialize(int* argc, char*** argv);

/**
 * Finalize a runtime initialized by Alea. All cluster contexts must first be
 * destroyed. MPI is finalized only when Alea initialized it.
 */
alea_cluster_status_t alea_cluster_finalize(void);

/** Create a cluster context. V1 permits one live context per process. */
alea_cluster_t* alea_cluster_create(void);
void alea_cluster_destroy(alea_cluster_t* cluster);

int alea_cluster_rank(const alea_cluster_t* cluster);
int alea_cluster_size(const alea_cluster_t* cluster);
int alea_cluster_is_root(const alea_cluster_t* cluster);
const char* alea_cluster_backend(const alea_cluster_t* cluster);

/**
 * Collectively estimate physical volumes on every rank.
 *
 * Every rank must call this function in the same order with equivalent
 * systems and identical scalar options. Output arrays are sized to
 * alea_volume_path_count(sys), and receive the final result on every rank.
 * An explicit sampling sphere in options is used identically on every rank;
 * the caller is responsible for its coverage. The progress callback, when
 * supplied, is invoked on rank zero only.
 */
alea_cluster_status_t alea_cluster_estimate_volumes(
    alea_cluster_t* cluster,
    alea_system_t* sys,
    const alea_volume_estimate_options_t* options,
    double* volumes,
    double* rel_errors,
    alea_cluster_volume_stats_t* out_stats);

const char* alea_cluster_status_string(alea_cluster_status_t status);

#ifdef __cplusplus
}
#endif

#endif /* ALEA_CLUSTER_H */
