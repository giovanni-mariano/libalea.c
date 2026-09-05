// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef TINYPAR_H
#define TINYPAR_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define TINYPAR_VERSION_MAJOR 0
#define TINYPAR_VERSION_MINOR 3
#define TINYPAR_VERSION_PATCH 0

typedef enum tinypar_status {
    TINYPAR_OK = 0,
    TINYPAR_INVALID_ARGUMENT,
    TINYPAR_ALLOCATION_FAILED,
    TINYPAR_THREAD_CREATE_FAILED,
    TINYPAR_THREAD_JOIN_FAILED,
    TINYPAR_CALLBACK_FAILED,
    TINYPAR_SYNCHRONIZATION_FAILED
} tinypar_status_t;

/**
 * Process one half-open range [begin, end). A callback may be called multiple
 * times for a worker. Returning non-zero cancels unclaimed work and makes
 * tinypar_parallel_for() return TINYPAR_CALLBACK_FAILED.
 */
typedef int (*tinypar_range_fn)(void* context, size_t worker_index,
                                size_t begin, size_t end);

typedef enum tinypar_schedule {
    /* Workers claim fixed-size chunks as they become available. */
    TINYPAR_SCHEDULE_DYNAMIC = 0,
    /* Each worker receives one deterministic, contiguous, balanced range. */
    TINYPAR_SCHEDULE_STATIC_BLOCK
} tinypar_schedule_t;

typedef struct tinypar_config {
    size_t item_count;
    /* Dynamic chunk size, or worker-grain hint for static scheduling. */
    size_t chunk_size;
    /* Zero selects the process default worker count. */
    size_t max_workers;
    tinypar_schedule_t schedule;
} tinypar_config_t;

typedef struct tinypar_executor tinypar_executor_t;

typedef struct tinypar_executor_config {
    /* Set by tinypar_executor_config_init(); reserved for ABI extension. */
    size_t struct_size;
    /* Zero selects the process default worker count. */
    size_t max_workers;
} tinypar_executor_config_t;

/** Initialize a configuration for dynamic scheduling with a chunk size of 1. */
void tinypar_config_init(tinypar_config_t* config);

/** Return non-zero when this build can create worker threads. */
int tinypar_threading_enabled(void);

/** Return non-zero inside a callback of a multi-worker invocation. */
int tinypar_in_parallel(void);

/** Return non-zero inside any TinyPar callback, including serial execution. */
int tinypar_in_callback(void);

/** Return the available logical processor count, never less than one. */
size_t tinypar_hardware_threads(void);

/** Return the process default worker count, never less than one. */
size_t tinypar_default_workers(void);

/** Set the process default worker count. Zero restores the hardware default. */
tinypar_status_t tinypar_set_default_workers(size_t workers);

/**
 * Return the number of workers that tinypar_parallel_for() would use.
 * Returns zero when chunk_size is zero for non-empty work.
 */
size_t tinypar_effective_workers(size_t item_count, size_t chunk_size,
                                 size_t max_workers);

/**
 * Run callback over [0, config->item_count). Dynamic scheduling assigns one
 * initial chunk to every effective worker, then distributes remaining chunks
 * as workers become available. Static-block scheduling assigns each worker
 * one deterministic contiguous range. The calling thread participates as
 * worker zero. A nested invocation inside a TinyPar callback executes serially.
 * Independent top-level invocations own separate state and may run
 * concurrently. On return, no worker can access callback or context.
 */
tinypar_status_t tinypar_parallel_for(const tinypar_config_t* config,
                                      tinypar_range_fn callback,
                                      void* context);

/** Initialize an executor configuration with the process worker default. */
void tinypar_executor_config_init(tinypar_executor_config_t* config);

/**
 * Create a reusable executor. The actual stable participant count can be less
 * than requested if only part of the native worker team can be created.
 */
tinypar_status_t tinypar_executor_create(
    const tinypar_executor_config_t* config,
    tinypar_executor_t** executor);

/**
 * Destroy an executor and set *executor to NULL. Calls using an executor must
 * finish before destruction begins. If worker termination is unknown,
 * destruction returns an error and retains *executor for a later attempt. A
 * cleanup-only error may be returned after safely releasing the executor.
 */
tinypar_status_t tinypar_executor_destroy(tinypar_executor_t** executor);

/** Return the stable number of participants owned by an executor. */
size_t tinypar_executor_workers(const tinypar_executor_t* executor);

/**
 * Run a synchronous parallel range on a reusable executor. Calls on the same
 * executor serialize; calls on separate executors may proceed concurrently.
 * Nested calls inside a TinyPar callback execute serially.
 */
tinypar_status_t tinypar_executor_parallel_for(
    tinypar_executor_t* executor,
    const tinypar_config_t* config,
    tinypar_range_fn callback,
    void* context);

/** Return a stable textual description of a tinypar status code. */
const char* tinypar_status_string(tinypar_status_t status);

#ifdef __cplusplus
}
#endif

#endif /* TINYPAR_H */
