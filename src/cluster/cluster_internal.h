// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_CLUSTER_INTERNAL_H
#define ALEA_CLUSTER_INTERNAL_H

#include "alea_cluster.h"

struct alea_cluster {
    void* backend;
    int rank;
    int size;
    int usable;
};

uint64_t alea_cluster_system_fingerprint(const alea_system_t* sys);
int alea_cluster_fingerprints_match(alea_cluster_t* cluster,
                                    uint64_t fingerprint);
alea_cluster_status_t alea_cluster_read_path_root(
    const char* path, char** data, size_t* length);
alea_cluster_status_t alea_cluster_broadcast_owned_bytes(
    alea_cluster_t* cluster, char* root_data, size_t root_length,
    alea_cluster_status_t root_status, char** data, size_t* length);

int alea_cluster_backend_initialize(int* argc, char*** argv);
int alea_cluster_backend_finalize(void);
int alea_cluster_backend_create(void** state, int* rank, int* size);
void alea_cluster_backend_destroy(void* state);
const char* alea_cluster_backend_name(void);
int alea_cluster_backend_agree_status(void* state, int local, int* global);
int alea_cluster_backend_u64_minmax(void* state, uint64_t local,
                                    uint64_t* minimum, uint64_t* maximum);
int alea_cluster_backend_sum_doubles(void* state, double* values,
                                     size_t count);
int alea_cluster_backend_broadcast_int(void* state, int* value, int root);
int alea_cluster_backend_broadcast_u64(void* state, uint64_t* value, int root);
int alea_cluster_backend_broadcast_bytes(void* state, void* bytes,
                                         size_t count, int root);
int alea_cluster_backend_scatter_blocks(void* state, const void* root_data,
                                        void* local_data, size_t item_size,
                                        size_t item_count, int root);
int alea_cluster_backend_gather_u64(void* state, uint64_t local,
                                    uint64_t* root_values, int root);
int alea_cluster_backend_gather_bytes(void* state, const void* local_data,
                                      size_t local_bytes, void* root_data,
                                      const size_t* root_bytes_by_rank,
                                      int root);

#endif
