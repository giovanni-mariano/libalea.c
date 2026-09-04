// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_parallel.h"

#include "tinypar.h"

int alea_parallel_enabled(void) {
    return tinypar_threading_enabled();
}

size_t alea_parallel_max_workers(void) {
    return tinypar_threading_enabled() ? tinypar_hardware_threads() : 1;
}

int alea_parallel_in_region(void) {
    return tinypar_in_parallel();
}

size_t alea_parallel_effective_workers(size_t item_count, size_t grain_size,
                                       size_t max_workers) {
    if (tinypar_in_parallel()) max_workers = 1;
    return tinypar_effective_workers(item_count, grain_size, max_workers);
}

alea_parallel_status_t alea_parallel_for(
    size_t item_count, size_t grain_size, size_t max_workers,
    alea_parallel_schedule_t schedule, alea_parallel_range_fn callback,
    void* context, size_t* out_actual_workers) {
    if (tinypar_in_parallel()) max_workers = 1;
    size_t workers = tinypar_effective_workers(
        item_count, grain_size, max_workers);
    if (out_actual_workers) *out_actual_workers = workers;

    tinypar_config_t config = {
        .item_count = item_count,
        .chunk_size = grain_size,
        .max_workers = max_workers,
        .schedule = schedule == ALEA_PARALLEL_STATIC_BLOCK
            ? TINYPAR_SCHEDULE_STATIC_BLOCK : TINYPAR_SCHEDULE_DYNAMIC
    };
    tinypar_status_t status = tinypar_parallel_for(&config, callback, context);
    switch (status) {
        case TINYPAR_OK: return ALEA_PARALLEL_OK;
        case TINYPAR_INVALID_ARGUMENT: return ALEA_PARALLEL_INVALID_ARGUMENT;
        case TINYPAR_ALLOCATION_FAILED:
            return ALEA_PARALLEL_ALLOCATION_FAILED;
        case TINYPAR_THREAD_CREATE_FAILED:
            return ALEA_PARALLEL_THREAD_CREATE_FAILED;
        case TINYPAR_THREAD_JOIN_FAILED:
            return ALEA_PARALLEL_THREAD_JOIN_FAILED;
        case TINYPAR_CALLBACK_FAILED: return ALEA_PARALLEL_CALLBACK_FAILED;
        default: return ALEA_PARALLEL_THREAD_CREATE_FAILED;
    }
}

const char* alea_parallel_status_string(alea_parallel_status_t status) {
    switch (status) {
        case ALEA_PARALLEL_OK: return tinypar_status_string(TINYPAR_OK);
        case ALEA_PARALLEL_INVALID_ARGUMENT:
            return tinypar_status_string(TINYPAR_INVALID_ARGUMENT);
        case ALEA_PARALLEL_ALLOCATION_FAILED:
            return tinypar_status_string(TINYPAR_ALLOCATION_FAILED);
        case ALEA_PARALLEL_THREAD_CREATE_FAILED:
            return tinypar_status_string(TINYPAR_THREAD_CREATE_FAILED);
        case ALEA_PARALLEL_THREAD_JOIN_FAILED:
            return tinypar_status_string(TINYPAR_THREAD_JOIN_FAILED);
        case ALEA_PARALLEL_CALLBACK_FAILED:
            return tinypar_status_string(TINYPAR_CALLBACK_FAILED);
        default: return "unknown parallel error";
    }
}
