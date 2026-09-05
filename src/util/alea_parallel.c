// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_parallel.h"

#include "tinypar.h"
#include "util/alea_atomic.h"

#include <stdlib.h>

static atomic_flag g_executor_lock = ATOMIC_FLAG_INIT;
static tinypar_executor_t* g_executor;
static int g_executor_cleanup_registered;

static void executor_lock(void) {
    while (atomic_flag_test_and_set(&g_executor_lock)) {
        /* Executor initialization happens once and does not run callbacks. */
    }
}

static void executor_unlock(void) {
    atomic_flag_clear(&g_executor_lock);
}

static alea_parallel_status_t map_status(tinypar_status_t status) {
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
        case TINYPAR_SYNCHRONIZATION_FAILED:
            return ALEA_PARALLEL_SYNCHRONIZATION_FAILED;
        default: return ALEA_PARALLEL_SYNCHRONIZATION_FAILED;
    }
}

static void destroy_process_executor(void) {
    executor_lock();
    if (g_executor) (void)tinypar_executor_destroy(&g_executor);
    executor_unlock();
}

static alea_parallel_status_t get_process_executor(
        tinypar_executor_t** executor) {
    executor_lock();
    if (!g_executor) {
        tinypar_executor_config_t config;
        tinypar_executor_config_init(&config);
        tinypar_status_t status = tinypar_executor_create(&config, &g_executor);
        if (status != TINYPAR_OK) {
            executor_unlock();
            return map_status(status);
        }
        if (!g_executor_cleanup_registered) {
            g_executor_cleanup_registered = atexit(destroy_process_executor) == 0;
        }
    }
    *executor = g_executor;
    executor_unlock();
    return ALEA_PARALLEL_OK;
}

int alea_parallel_enabled(void) {
    return tinypar_threading_enabled();
}

size_t alea_parallel_max_workers(void) {
    executor_lock();
    size_t workers = g_executor
        ? tinypar_executor_workers(g_executor) : tinypar_default_workers();
    executor_unlock();
    return workers;
}

int alea_parallel_in_region(void) {
    return tinypar_in_callback();
}

alea_parallel_status_t alea_parallel_set_default_workers(size_t workers) {
    executor_lock();
    if (g_executor) {
        executor_unlock();
        return ALEA_PARALLEL_ALREADY_INITIALIZED;
    }
    tinypar_status_t status = tinypar_set_default_workers(workers);
    executor_unlock();
    return map_status(status);
}

size_t alea_parallel_effective_workers(size_t item_count, size_t grain_size,
                                       size_t max_workers) {
    if (tinypar_in_callback()) max_workers = 1;
    else {
        size_t executor_workers = alea_parallel_max_workers();
        if (max_workers == 0 || max_workers > executor_workers)
            max_workers = executor_workers;
    }
    return tinypar_effective_workers(item_count, grain_size, max_workers);
}

alea_parallel_status_t alea_parallel_for(
    size_t item_count, size_t grain_size, size_t max_workers,
    alea_parallel_schedule_t schedule, alea_parallel_range_fn callback,
    void* context, size_t* out_actual_workers) {
    if (!callback || (item_count != 0 && grain_size == 0) ||
        (schedule != ALEA_PARALLEL_DYNAMIC &&
         schedule != ALEA_PARALLEL_STATIC_BLOCK)) {
        if (out_actual_workers) *out_actual_workers = 0;
        return ALEA_PARALLEL_INVALID_ARGUMENT;
    }
    if (item_count == 0) {
        if (out_actual_workers) *out_actual_workers = 1;
        return ALEA_PARALLEL_OK;
    }

    tinypar_config_t config = {
        .item_count = item_count,
        .chunk_size = grain_size,
        .max_workers = max_workers,
        .schedule = schedule == ALEA_PARALLEL_STATIC_BLOCK
            ? TINYPAR_SCHEDULE_STATIC_BLOCK : TINYPAR_SCHEDULE_DYNAMIC
    };

    /* Select the lexical nested path before touching the process executor.
     * This prevents a first nested call from creating an otherwise idle team. */
    if (tinypar_in_callback()) {
        config.max_workers = 1;
        if (out_actual_workers) *out_actual_workers = 1;
        return map_status(tinypar_parallel_for(&config, callback, context));
    }

    size_t available_workers = alea_parallel_max_workers();
    if (max_workers == 0 || max_workers > available_workers)
        max_workers = available_workers;
    size_t workers = tinypar_effective_workers(
        item_count, grain_size, max_workers);
    if (workers == 1) {
        config.max_workers = 1;
        if (out_actual_workers) *out_actual_workers = 1;
        return map_status(tinypar_parallel_for(&config, callback, context));
    }

    tinypar_executor_t* executor = NULL;
    alea_parallel_status_t executor_status =
        get_process_executor(&executor);
    if (executor_status != ALEA_PARALLEL_OK) {
        if (out_actual_workers) *out_actual_workers = 0;
        return executor_status;
    }

    size_t executor_workers = tinypar_executor_workers(executor);
    if (max_workers > executor_workers)
        max_workers = executor_workers;
    workers = tinypar_effective_workers(
        item_count, grain_size, max_workers);
    if (out_actual_workers) *out_actual_workers = workers;

    config.max_workers = max_workers;
    tinypar_status_t status = tinypar_executor_parallel_for(
        executor, &config, callback, context);
    return map_status(status);
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
        case ALEA_PARALLEL_SYNCHRONIZATION_FAILED:
            return tinypar_status_string(TINYPAR_SYNCHRONIZATION_FAILED);
        case ALEA_PARALLEL_ALREADY_INITIALIZED:
            return "parallel runtime already initialized";
        default: return "unknown parallel error";
    }
}
