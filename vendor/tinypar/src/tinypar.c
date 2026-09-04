// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "tinypar.h"
#include "tinypar_atomic.h"
#include "tinypar_platform.h"

#include <stdint.h>
#include <stdlib.h>

typedef struct tinypar_job {
    tinypar_range_fn callback;
    void* context;
    size_t item_count;
    size_t chunk_size;
    size_t chunk_count;
    size_t worker_count;
    tinypar_schedule_t schedule;
    tinypar_atomic_size_t next_chunk;
    tinypar_atomic_int_t cancelled;
    tinypar_atomic_int_t status;
} tinypar_job_t;

typedef struct tinypar_worker_arg {
    tinypar_job_t* job;
    size_t worker_index;
} tinypar_worker_arg_t;

#if defined(_MSC_VER) && !defined(__clang__)
#define TINYPAR_THREAD_LOCAL __declspec(thread)
#else
#define TINYPAR_THREAD_LOCAL _Thread_local
#endif

static TINYPAR_THREAD_LOCAL unsigned tinypar_parallel_depth;

static void tinypar_fail_job(tinypar_job_t* job) {
    int expected = TINYPAR_OK;
    (void)tinypar_atomic_int_compare_exchange(
        &job->status, &expected, TINYPAR_CALLBACK_FAILED);
    tinypar_atomic_int_store(&job->cancelled, 1);
}

static int tinypar_claim_chunk(tinypar_job_t* job, size_t* chunk) {
    size_t next = tinypar_atomic_size_load(&job->next_chunk);
    while (next < job->chunk_count) {
        if (tinypar_atomic_size_compare_exchange(
                &job->next_chunk, &next, next + 1)) {
            *chunk = next;
            return 1;
        }
    }
    return 0;
}

static void tinypar_run_dynamic(tinypar_worker_arg_t* argument) {
    tinypar_job_t* job = argument->job;

    while (tinypar_atomic_int_load(&job->cancelled) == 0) {
        size_t chunk;
        if (!tinypar_claim_chunk(job, &chunk)) return;
        size_t begin = chunk * job->chunk_size;

        size_t remaining = job->item_count - begin;
        size_t end = remaining < job->chunk_size
            ? job->item_count : begin + job->chunk_size;
        if (tinypar_atomic_int_load(&job->cancelled) != 0) return;
        if (job->callback(job->context, argument->worker_index, begin, end) != 0) {
            tinypar_fail_job(job);
            return;
        }
    }
}

static void tinypar_run_static(tinypar_worker_arg_t* argument) {
    tinypar_job_t* job = argument->job;
    size_t worker = argument->worker_index;
    size_t base = job->item_count / job->worker_count;
    size_t remainder = job->item_count % job->worker_count;
    size_t extra_before = worker < remainder ? worker : remainder;
    size_t begin = worker * base + extra_before;
    size_t end = begin + base + (worker < remainder ? 1 : 0);

    if (tinypar_atomic_int_load(&job->cancelled) == 0 && begin < end &&
        job->callback(job->context, worker, begin, end) != 0) {
        tinypar_fail_job(job);
    }
}

static void tinypar_run_worker(tinypar_worker_arg_t* argument) {
    tinypar_job_t* job = argument->job;
    int is_parallel = job->worker_count > 1;
    if (is_parallel) tinypar_parallel_depth++;
    if (job->schedule == TINYPAR_SCHEDULE_STATIC_BLOCK)
        tinypar_run_static(argument);
    else
        tinypar_run_dynamic(argument);
    if (is_parallel) tinypar_parallel_depth--;
}

#if defined(_WIN32)
static unsigned __stdcall tinypar_worker_entry(void* argument) {
    tinypar_run_worker((tinypar_worker_arg_t*)argument);
    return 0;
}
#else
static void* tinypar_worker_entry(void* argument) {
    tinypar_run_worker((tinypar_worker_arg_t*)argument);
    return NULL;
}
#endif

size_t tinypar_hardware_threads(void) {
    size_t count = tinypar_platform_hardware_threads();
    return count == 0 ? 1 : count;
}

void tinypar_config_init(tinypar_config_t* config) {
    if (!config) return;
    config->item_count = 0;
    config->chunk_size = 1;
    config->max_workers = 0;
    config->schedule = TINYPAR_SCHEDULE_DYNAMIC;
}

int tinypar_threading_enabled(void) {
    return tinypar_platform_threading_enabled();
}

int tinypar_in_parallel(void) {
    return tinypar_parallel_depth != 0;
}

size_t tinypar_effective_workers(size_t item_count, size_t chunk_size,
                                 size_t max_workers) {
    if (item_count == 0) return 1;
    if (chunk_size == 0) return 0;
    if (!tinypar_threading_enabled()) return 1;

    size_t requested = max_workers == 0 ? tinypar_hardware_threads() : max_workers;
    if (requested == 0) requested = 1;

    size_t ranges = 1 + (item_count - 1) / chunk_size;
    return requested < ranges ? requested : ranges;
}

tinypar_status_t tinypar_parallel_for(const tinypar_config_t* config,
                                      tinypar_range_fn callback,
                                      void* context) {
    if (!config || !callback) return TINYPAR_INVALID_ARGUMENT;
    if (config->item_count == 0) return TINYPAR_OK;
    if (config->chunk_size == 0) return TINYPAR_INVALID_ARGUMENT;
    if (config->schedule != TINYPAR_SCHEDULE_DYNAMIC &&
        config->schedule != TINYPAR_SCHEDULE_STATIC_BLOCK)
        return TINYPAR_INVALID_ARGUMENT;

    size_t worker_count = tinypar_effective_workers(
        config->item_count, config->chunk_size, config->max_workers);
    if (worker_count == 0) return TINYPAR_INVALID_ARGUMENT;

    tinypar_job_t job;
    job.callback = callback;
    job.context = context;
    job.item_count = config->item_count;
    job.chunk_size = config->chunk_size;
    job.chunk_count = 1 + (config->item_count - 1) / config->chunk_size;
    job.worker_count = worker_count;
    job.schedule = config->schedule;
    tinypar_atomic_size_init(&job.next_chunk, 0);
    tinypar_atomic_int_init(&job.cancelled, 0);
    tinypar_atomic_int_init(&job.status, TINYPAR_OK);

    tinypar_worker_arg_t caller = { &job, 0 };
    if (worker_count == 1) {
        tinypar_run_worker(&caller);
        return (tinypar_status_t)tinypar_atomic_int_load(&job.status);
    }

    size_t spawned_count = worker_count - 1;
    if (spawned_count > SIZE_MAX / sizeof(tinypar_thread_t) ||
        spawned_count > SIZE_MAX / sizeof(tinypar_worker_arg_t))
        return TINYPAR_ALLOCATION_FAILED;
    tinypar_thread_t* threads = calloc(spawned_count, sizeof(*threads));
    tinypar_worker_arg_t* arguments = calloc(spawned_count, sizeof(*arguments));
    if (!threads || !arguments) {
        free(arguments);
        free(threads);
        return TINYPAR_ALLOCATION_FAILED;
    }

    size_t started = 0;
    for (; started < spawned_count; started++) {
        arguments[started].job = &job;
        arguments[started].worker_index = started + 1;
        if (!tinypar_thread_start(&threads[started], tinypar_worker_entry,
                                  &arguments[started])) {
            tinypar_atomic_int_store(&job.cancelled, 1);
            break;
        }
    }

    tinypar_status_t result = TINYPAR_OK;
    if (started != spawned_count) {
        result = TINYPAR_THREAD_CREATE_FAILED;
    } else {
        tinypar_run_worker(&caller);
    }

    for (size_t i = 0; i < started; i++) {
        if (!tinypar_thread_join(&threads[i])) result = TINYPAR_THREAD_JOIN_FAILED;
    }

    free(arguments);
    free(threads);

    if (result != TINYPAR_OK) return result;
    return (tinypar_status_t)tinypar_atomic_int_load(&job.status);
}

const char* tinypar_status_string(tinypar_status_t status) {
    switch (status) {
    case TINYPAR_OK: return "success";
    case TINYPAR_INVALID_ARGUMENT: return "invalid argument";
    case TINYPAR_ALLOCATION_FAILED: return "allocation failed";
    case TINYPAR_THREAD_CREATE_FAILED: return "thread creation failed";
    case TINYPAR_THREAD_JOIN_FAILED: return "thread join failed";
    case TINYPAR_CALLBACK_FAILED: return "callback failed";
    default: return "unknown tinypar status";
    }
}
