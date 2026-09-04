// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_PARALLEL_H
#define ALEA_PARALLEL_H

#include <stddef.h>

typedef enum alea_parallel_schedule {
    ALEA_PARALLEL_DYNAMIC = 0,
    ALEA_PARALLEL_STATIC_BLOCK
} alea_parallel_schedule_t;

typedef int (*alea_parallel_range_fn)(void* context, size_t worker_index,
                                      size_t begin, size_t end);

typedef enum alea_parallel_status {
    ALEA_PARALLEL_OK = 0,
    ALEA_PARALLEL_INVALID_ARGUMENT,
    ALEA_PARALLEL_ALLOCATION_FAILED,
    ALEA_PARALLEL_THREAD_CREATE_FAILED,
    ALEA_PARALLEL_THREAD_JOIN_FAILED,
    ALEA_PARALLEL_CALLBACK_FAILED
} alea_parallel_status_t;

int alea_parallel_enabled(void);
size_t alea_parallel_max_workers(void);
int alea_parallel_in_region(void);

size_t alea_parallel_effective_workers(size_t item_count, size_t grain_size,
                                       size_t max_workers);

alea_parallel_status_t alea_parallel_for(
    size_t item_count, size_t grain_size, size_t max_workers,
    alea_parallel_schedule_t schedule, alea_parallel_range_fn callback,
    void* context, size_t* out_actual_workers);

const char* alea_parallel_status_string(alea_parallel_status_t status);

#endif /* ALEA_PARALLEL_H */
