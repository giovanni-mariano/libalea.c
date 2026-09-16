// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_CLUSTER_BASE_H
#define ALEA_CLUSTER_BASE_H

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
    ALEA_CLUSTER_INTERRUPTED,
    ALEA_CLUSTER_IO_ERROR,
    ALEA_CLUSTER_OUTPUT_LIMIT
} alea_cluster_status_t;

/** Initialize the process-wide local or MPI runtime before creating a context.
 * MPI calls must use the initializing thread. If MPI is already initialized,
 * libalea attaches without taking ownership and requires at least
 * MPI_THREAD_FUNNELED; otherwise it initializes MPI with that thread level. */
alea_cluster_status_t alea_cluster_initialize(int* argc, char*** argv);
/** Close the libalea runtime after destroying all contexts. MPI is finalized
 * only when libalea initialized it; application-owned MPI remains active. */
alea_cluster_status_t alea_cluster_finalize(void);
/** Collectively create a context on MPI_COMM_WORLD in an MPI build, or a
 * one-rank local context otherwise. The MPI backend duplicates the world
 * communicator, so library messages do not share the caller's context.
 * One live cluster context is permitted per process. */
alea_cluster_t* alea_cluster_create(void);
/** Destroy the context collectively on its MPI group. This frees libalea's
 * duplicate communicator, but does not finalize MPI. */
void alea_cluster_destroy(alea_cluster_t* cluster);

int alea_cluster_rank(const alea_cluster_t* cluster);
int alea_cluster_size(const alea_cluster_t* cluster);
int alea_cluster_is_root(const alea_cluster_t* cluster);
const char* alea_cluster_backend(const alea_cluster_t* cluster);

/** Collectively read one binary file on root and broadcast its bytes. Only
 * root uses path; every rank receives malloc-owned, NUL-terminated data and
 * frees it with free(). Outputs are reset on error. */
alea_cluster_status_t alea_cluster_read_file(
    alea_cluster_t* cluster, const char* path, char** data, size_t* length);
/** Expand nested MCNP READ FILE= cards on root, resolving relative paths
 * against the containing file, then broadcast the text. The returned buffer
 * is malloc-owned on every rank. Missing files return IO_ERROR; malformed
 * cards or nesting beyond 32 files return COMPUTE_ERROR. */
alea_cluster_status_t alea_cluster_read_mcnp_input(
    alea_cluster_t* cluster, const char* path, char** data, size_t* length);
/** Collectively agree on a local status before the next operation. */
alea_cluster_status_t alea_cluster_agree(
    alea_cluster_t* cluster, alea_cluster_status_t local_status);

const char* alea_cluster_status_string(alea_cluster_status_t status);

#ifdef __cplusplus
}
#endif

#endif /* ALEA_CLUSTER_BASE_H */
