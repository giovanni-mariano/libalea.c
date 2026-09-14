// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_CLUSTER_INTERNAL_H
#define ALEA_CLUSTER_INTERNAL_H

#include "alea_cluster.h"

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

#endif
