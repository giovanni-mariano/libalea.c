// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "tinypar.h"
#include "tinypar_atomic.h"
#include "tinypar_platform.h"

#include <stdint.h>
#include <stdio.h>
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

typedef enum tinypar_start_state {
    TINYPAR_START_WAITING = 0,
    TINYPAR_START_RUN,
    TINYPAR_START_ABORT
} tinypar_start_state_t;

typedef struct tinypar_start_gate {
    tinypar_mutex_t mutex;
    tinypar_condition_t condition;
    tinypar_start_state_t state;
    size_t ready_workers;
} tinypar_start_gate_t;

typedef struct tinypar_worker_arg {
    tinypar_job_t* job;
    tinypar_start_gate_t* start_gate;
    size_t worker_index;
} tinypar_worker_arg_t;

typedef struct tinypar_executor_worker_arg {
    tinypar_executor_t* executor;
    size_t worker_index;
} tinypar_executor_worker_arg_t;

struct tinypar_executor {
    tinypar_mutex_t submission_mutex;
    tinypar_mutex_t state_mutex;
    tinypar_condition_t done_condition;
    tinypar_condition_t* worker_conditions;
    tinypar_thread_t* threads;
    tinypar_executor_worker_arg_t* arguments;
    unsigned char* terminated;
    size_t worker_count;
    size_t spawned_count;
    size_t generation;
    size_t completed_workers;
    size_t active_workers;
    size_t initialized_worker_conditions;
    size_t ready_workers;
    int shutdown;
    int cleanup_failed;
    tinypar_job_t job;
};

#if defined(_MSC_VER) && !defined(__clang__)
#define TINYPAR_THREAD_LOCAL __declspec(thread)
#else
#define TINYPAR_THREAD_LOCAL _Thread_local
#endif

static TINYPAR_THREAD_LOCAL unsigned tinypar_parallel_depth;
static TINYPAR_THREAD_LOCAL unsigned tinypar_callback_depth;
static tinypar_atomic_size_t tinypar_default_worker_limit = 0;

static void tinypar_fatal_runtime_failure(const char* reason) {
    fprintf(stderr, "tinypar: unrecoverable runtime failure: %s\n", reason);
    abort();
}

static void tinypar_require_sync(int succeeded) {
    if (!succeeded)
        tinypar_fatal_runtime_failure("native synchronization failed");
}

static int tinypar_valid_schedule(tinypar_schedule_t schedule) {
    return schedule == TINYPAR_SCHEDULE_DYNAMIC ||
           schedule == TINYPAR_SCHEDULE_STATIC_BLOCK;
}

static void tinypar_job_init(tinypar_job_t* job,
                             const tinypar_config_t* config,
                             size_t worker_count,
                             tinypar_range_fn callback,
                             void* context) {
    job->callback = callback;
    job->context = context;
    job->item_count = config->item_count;
    job->chunk_size = config->chunk_size;
    job->chunk_count = 1 + (config->item_count - 1) / config->chunk_size;
    job->worker_count = worker_count;
    job->schedule = config->schedule;
    /* Reserve one initial chunk for every participant. This retains dynamic
     * scheduling for the remaining chunks while preventing one fast worker
     * from consuming a short job before the rest of the team runs. */
    tinypar_atomic_size_init(&job->next_chunk, worker_count);
    tinypar_atomic_int_init(&job->cancelled, 0);
    tinypar_atomic_int_init(&job->status, TINYPAR_OK);
}

static void tinypar_fail_job(tinypar_job_t* job) {
    int expected = TINYPAR_OK;
    (void)tinypar_atomic_int_compare_exchange(
        &job->status, &expected, TINYPAR_CALLBACK_FAILED);
    tinypar_atomic_int_store(&job->cancelled, 1);
}

static int tinypar_claim_chunk(tinypar_job_t* job, size_t* chunk) {
    if (job->chunk_count < SIZE_MAX) {
        size_t next = tinypar_atomic_size_fetch_add(&job->next_chunk, 1);
        if (next >= job->chunk_count) return 0;
        *chunk = next;
        return 1;
    }

    /* Preserve the full size_t input domain without allowing the claim counter
     * to wrap after the final representable chunk. */
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

static int tinypar_run_dynamic_chunk(tinypar_worker_arg_t* argument,
                                     size_t chunk) {
    tinypar_job_t* job = argument->job;
    size_t begin = chunk * job->chunk_size;
    size_t remaining = job->item_count - begin;
    size_t end = remaining < job->chunk_size
        ? job->item_count : begin + job->chunk_size;
    if (tinypar_atomic_int_load(&job->cancelled) != 0) return 0;
    if (job->callback(job->context, argument->worker_index,
                      begin, end) != 0) {
        tinypar_fail_job(job);
        return 0;
    }
    return 1;
}

static void tinypar_run_dynamic(tinypar_worker_arg_t* argument) {
    tinypar_job_t* job = argument->job;

    if (argument->worker_index < job->chunk_count &&
        !tinypar_run_dynamic_chunk(argument, argument->worker_index))
        return;

    while (tinypar_atomic_int_load(&job->cancelled) == 0) {
        size_t chunk;
        if (!tinypar_claim_chunk(job, &chunk)) return;
        if (!tinypar_run_dynamic_chunk(argument, chunk)) return;
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
    tinypar_callback_depth++;
    if (job->schedule == TINYPAR_SCHEDULE_STATIC_BLOCK)
        tinypar_run_static(argument);
    else
        tinypar_run_dynamic(argument);
    tinypar_callback_depth--;
    if (is_parallel) tinypar_parallel_depth--;
}

static int tinypar_start_gate_init(tinypar_start_gate_t* gate) {
    if (!tinypar_mutex_init(&gate->mutex)) return 0;
    if (!tinypar_condition_init(&gate->condition)) {
        tinypar_require_sync(tinypar_mutex_destroy(&gate->mutex));
        return 0;
    }
    gate->state = TINYPAR_START_WAITING;
    gate->ready_workers = 0;
    return 1;
}

static void tinypar_start_gate_destroy(tinypar_start_gate_t* gate) {
    tinypar_require_sync(tinypar_condition_destroy(&gate->condition));
    tinypar_require_sync(tinypar_mutex_destroy(&gate->mutex));
}

static int tinypar_start_gate_wait(tinypar_start_gate_t* gate) {
    tinypar_require_sync(tinypar_mutex_lock(&gate->mutex));
    gate->ready_workers++;
    tinypar_require_sync(tinypar_condition_broadcast(&gate->condition));
    while (gate->state == TINYPAR_START_WAITING)
        tinypar_require_sync(
            tinypar_condition_wait(&gate->condition, &gate->mutex));
    int run = gate->state == TINYPAR_START_RUN;
    tinypar_require_sync(tinypar_mutex_unlock(&gate->mutex));
    return run;
}

static void tinypar_start_gate_publish(tinypar_start_gate_t* gate,
                                       tinypar_start_state_t state,
                                       size_t expected_ready_workers) {
    tinypar_require_sync(tinypar_mutex_lock(&gate->mutex));
    while (state == TINYPAR_START_RUN &&
           gate->ready_workers != expected_ready_workers) {
        tinypar_require_sync(
            tinypar_condition_wait(&gate->condition, &gate->mutex));
    }
    gate->state = state;
    tinypar_require_sync(tinypar_condition_broadcast(&gate->condition));
    tinypar_require_sync(tinypar_mutex_unlock(&gate->mutex));
}

#if defined(_WIN32)
static unsigned __stdcall tinypar_worker_entry(void* opaque) {
    tinypar_worker_arg_t* argument = opaque;
    tinypar_platform_worker_enter(argument->worker_index);
    if (tinypar_start_gate_wait(argument->start_gate))
        tinypar_run_worker(argument);
    return 0;
}
#else
static void* tinypar_worker_entry(void* opaque) {
    tinypar_worker_arg_t* argument = opaque;
    tinypar_platform_worker_enter(argument->worker_index);
    if (tinypar_start_gate_wait(argument->start_gate))
        tinypar_run_worker(argument);
    return NULL;
}
#endif

static void tinypar_join_one_shot(tinypar_thread_t* threads,
                                  size_t started,
                                  tinypar_status_t* result) {
    for (size_t i = 0; i < started; i++) {
        tinypar_join_result_t joined = tinypar_thread_join(&threads[i]);
        if (joined == TINYPAR_JOIN_TERMINATION_UNKNOWN)
            tinypar_fatal_runtime_failure(
                "worker termination could not be established");
        if (joined == TINYPAR_JOIN_TERMINATED_CLEANUP_FAILED)
            *result = TINYPAR_THREAD_JOIN_FAILED;
    }
}

size_t tinypar_hardware_threads(void) {
    size_t count = tinypar_platform_hardware_threads();
    return count == 0 ? 1 : count;
}

size_t tinypar_default_workers(void) {
    if (!tinypar_threading_enabled()) return 1;
    size_t configured = tinypar_atomic_size_load(&tinypar_default_worker_limit);
    return configured == 0 ? tinypar_hardware_threads() : configured;
}

tinypar_status_t tinypar_set_default_workers(size_t workers) {
    tinypar_atomic_size_store(&tinypar_default_worker_limit, workers);
    return TINYPAR_OK;
}

void tinypar_config_init(tinypar_config_t* config) {
    if (!config) return;
    config->item_count = 0;
    config->chunk_size = 1;
    config->max_workers = 0;
    config->schedule = TINYPAR_SCHEDULE_DYNAMIC;
}

void tinypar_executor_config_init(tinypar_executor_config_t* config) {
    if (!config) return;
    config->struct_size = sizeof(*config);
    config->max_workers = 0;
}

int tinypar_threading_enabled(void) {
    return tinypar_platform_threading_enabled();
}

int tinypar_in_parallel(void) {
    return tinypar_parallel_depth != 0;
}

int tinypar_in_callback(void) {
    return tinypar_callback_depth != 0;
}

size_t tinypar_effective_workers(size_t item_count, size_t chunk_size,
                                 size_t max_workers) {
    if (item_count == 0) return 1;
    if (chunk_size == 0) return 0;
    if (!tinypar_threading_enabled()) return 1;

    size_t requested = max_workers == 0
        ? tinypar_default_workers() : max_workers;
    if (requested == 0) requested = 1;
    size_t ranges = 1 + (item_count - 1) / chunk_size;
    return requested < ranges ? requested : ranges;
}

static tinypar_status_t tinypar_run_serial(const tinypar_config_t* config,
                                           tinypar_range_fn callback,
                                           void* context) {
    tinypar_job_t job;
    tinypar_job_init(&job, config, 1, callback, context);
    tinypar_worker_arg_t caller = { &job, NULL, 0 };
    tinypar_run_worker(&caller);
    return (tinypar_status_t)tinypar_atomic_int_load(&job.status);
}

tinypar_status_t tinypar_parallel_for(const tinypar_config_t* config,
                                      tinypar_range_fn callback,
                                      void* context) {
    if (!config || !callback) return TINYPAR_INVALID_ARGUMENT;
    if (config->item_count == 0) return TINYPAR_OK;
    if (config->chunk_size == 0 || !tinypar_valid_schedule(config->schedule))
        return TINYPAR_INVALID_ARGUMENT;
    if (tinypar_in_callback())
        return tinypar_run_serial(config, callback, context);

    size_t worker_count = tinypar_effective_workers(
        config->item_count, config->chunk_size, config->max_workers);
    if (worker_count == 0) return TINYPAR_INVALID_ARGUMENT;
    if (worker_count == 1)
        return tinypar_run_serial(config, callback, context);

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

    tinypar_start_gate_t start_gate;
    if (!tinypar_start_gate_init(&start_gate)) {
        free(arguments);
        free(threads);
        return TINYPAR_SYNCHRONIZATION_FAILED;
    }

    tinypar_job_t job;
    tinypar_job_init(&job, config, worker_count, callback, context);
    size_t started = 0;
    for (; started < spawned_count; started++) {
        arguments[started].job = &job;
        arguments[started].start_gate = &start_gate;
        arguments[started].worker_index = started + 1;
        if (!tinypar_thread_start(&threads[started], tinypar_worker_entry,
                                  &arguments[started]))
            break;
    }

    tinypar_status_t result = TINYPAR_OK;
    if (started != spawned_count) {
        result = TINYPAR_THREAD_CREATE_FAILED;
        tinypar_start_gate_publish(
            &start_gate, TINYPAR_START_ABORT, 0);
    } else {
        tinypar_start_gate_publish(
            &start_gate, TINYPAR_START_RUN, spawned_count);
        tinypar_worker_arg_t caller = { &job, NULL, 0 };
        tinypar_run_worker(&caller);
    }

    tinypar_join_one_shot(threads, started, &result);
    tinypar_start_gate_destroy(&start_gate);
    free(arguments);
    free(threads);

    if (result != TINYPAR_OK) return result;
    return (tinypar_status_t)tinypar_atomic_int_load(&job.status);
}

#if defined(_WIN32)
static unsigned __stdcall tinypar_executor_worker_entry(void* opaque)
#else
static void* tinypar_executor_worker_entry(void* opaque)
#endif
{
    tinypar_executor_worker_arg_t* argument = opaque;
    tinypar_executor_t* executor = argument->executor;
    tinypar_platform_worker_enter(argument->worker_index);
    size_t observed_generation = 0;

    tinypar_require_sync(tinypar_mutex_lock(&executor->state_mutex));
    executor->ready_workers++;
    tinypar_require_sync(
        tinypar_condition_signal(&executor->done_condition));
    for (;;) {
        while (!executor->shutdown &&
               executor->generation == observed_generation) {
            tinypar_require_sync(tinypar_condition_wait(
                &executor->worker_conditions[argument->worker_index - 1],
                &executor->state_mutex));
        }
        if (executor->shutdown) {
            tinypar_require_sync(tinypar_mutex_unlock(&executor->state_mutex));
#if defined(_WIN32)
            return 0;
#else
            return NULL;
#endif
        }

        observed_generation = executor->generation;
        if (argument->worker_index > executor->active_workers)
            continue;
        tinypar_job_t* job = &executor->job;
        tinypar_require_sync(tinypar_mutex_unlock(&executor->state_mutex));
        tinypar_worker_arg_t worker = {
            job, NULL, argument->worker_index
        };
        tinypar_run_worker(&worker);

        tinypar_require_sync(tinypar_mutex_lock(&executor->state_mutex));
        executor->completed_workers++;
        if (executor->completed_workers == executor->active_workers)
            tinypar_require_sync(
                tinypar_condition_signal(&executor->done_condition));
    }
}

static void tinypar_executor_destroy_sync(tinypar_executor_t* executor) {
    for (size_t i = 0; i < executor->initialized_worker_conditions; i++) {
        tinypar_require_sync(
            tinypar_condition_destroy(&executor->worker_conditions[i]));
    }
    tinypar_require_sync(
        tinypar_condition_destroy(&executor->done_condition));
    tinypar_require_sync(tinypar_mutex_destroy(&executor->state_mutex));
    tinypar_require_sync(tinypar_mutex_destroy(&executor->submission_mutex));
}

tinypar_status_t tinypar_executor_create(
        const tinypar_executor_config_t* config,
        tinypar_executor_t** output) {
    if (!output) return TINYPAR_INVALID_ARGUMENT;
    *output = NULL;
    if (config && config->struct_size != sizeof(*config))
        return TINYPAR_INVALID_ARGUMENT;

    size_t requested = config && config->max_workers != 0
        ? config->max_workers : tinypar_default_workers();
    if (!tinypar_threading_enabled()) requested = 1;
    if (requested == 0) requested = 1;

    tinypar_executor_t* executor = calloc(1, sizeof(*executor));
    if (!executor) return TINYPAR_ALLOCATION_FAILED;

    int submission_initialized = 0;
    int state_initialized = 0;
    if (!(submission_initialized =
              tinypar_mutex_init(&executor->submission_mutex)) ||
        !(state_initialized = tinypar_mutex_init(&executor->state_mutex)) ||
        !tinypar_condition_init(&executor->done_condition)) {
        if (state_initialized)
            tinypar_require_sync(tinypar_mutex_destroy(&executor->state_mutex));
        if (submission_initialized)
            tinypar_require_sync(
                tinypar_mutex_destroy(&executor->submission_mutex));
        free(executor);
        return TINYPAR_SYNCHRONIZATION_FAILED;
    }

    size_t wanted_threads = requested - 1;
    if (wanted_threads > SIZE_MAX / sizeof(*executor->threads) ||
        wanted_threads > SIZE_MAX / sizeof(*executor->arguments) ||
        wanted_threads > SIZE_MAX / sizeof(*executor->worker_conditions)) {
        tinypar_executor_destroy_sync(executor);
        free(executor);
        return TINYPAR_ALLOCATION_FAILED;
    }

    if (wanted_threads != 0) {
        executor->worker_conditions = calloc(
            wanted_threads, sizeof(*executor->worker_conditions));
        executor->threads = calloc(wanted_threads, sizeof(*executor->threads));
        executor->arguments = calloc(
            wanted_threads, sizeof(*executor->arguments));
        executor->terminated = calloc(
            wanted_threads, sizeof(*executor->terminated));
        if (!executor->worker_conditions || !executor->threads ||
            !executor->arguments ||
            !executor->terminated) {
            free(executor->terminated);
            free(executor->arguments);
            free(executor->threads);
            free(executor->worker_conditions);
            tinypar_executor_destroy_sync(executor);
            free(executor);
            return TINYPAR_ALLOCATION_FAILED;
        }
        for (; executor->initialized_worker_conditions < wanted_threads;
             executor->initialized_worker_conditions++) {
            if (!tinypar_condition_init(
                    &executor->worker_conditions[
                        executor->initialized_worker_conditions])) {
                tinypar_executor_destroy_sync(executor);
                free(executor->terminated);
                free(executor->arguments);
                free(executor->threads);
                free(executor->worker_conditions);
                free(executor);
                return TINYPAR_SYNCHRONIZATION_FAILED;
            }
        }
    }

    for (; executor->spawned_count < wanted_threads;
         executor->spawned_count++) {
        size_t index = executor->spawned_count;
        executor->arguments[index].executor = executor;
        executor->arguments[index].worker_index = index + 1;
        if (!tinypar_thread_start(&executor->threads[index],
                                  tinypar_executor_worker_entry,
                                  &executor->arguments[index]))
            break;
    }
    executor->worker_count = executor->spawned_count + 1;
    tinypar_require_sync(tinypar_mutex_lock(&executor->state_mutex));
    while (executor->ready_workers != executor->spawned_count) {
        tinypar_require_sync(tinypar_condition_wait(
            &executor->done_condition, &executor->state_mutex));
    }
    tinypar_require_sync(tinypar_mutex_unlock(&executor->state_mutex));
    *output = executor;
    return TINYPAR_OK;
}

size_t tinypar_executor_workers(const tinypar_executor_t* executor) {
    return executor ? executor->worker_count : 0;
}

void tinypar_executor_abandon_after_fork(tinypar_executor_t** executor) {
    if (executor) *executor = NULL;
}

tinypar_status_t tinypar_executor_parallel_for(
        tinypar_executor_t* executor,
        const tinypar_config_t* config,
        tinypar_range_fn callback,
        void* context) {
    if (!executor || !config || !callback)
        return TINYPAR_INVALID_ARGUMENT;
    if (config->item_count == 0) return TINYPAR_OK;
    if (config->chunk_size == 0 || !tinypar_valid_schedule(config->schedule))
        return TINYPAR_INVALID_ARGUMENT;
    if (tinypar_in_callback())
        return tinypar_run_serial(config, callback, context);

    tinypar_require_sync(tinypar_mutex_lock(&executor->submission_mutex));
    tinypar_require_sync(tinypar_mutex_lock(&executor->state_mutex));
    if (executor->shutdown) {
        tinypar_require_sync(tinypar_mutex_unlock(&executor->state_mutex));
        tinypar_require_sync(
            tinypar_mutex_unlock(&executor->submission_mutex));
        return TINYPAR_THREAD_JOIN_FAILED;
    }
    tinypar_require_sync(tinypar_mutex_unlock(&executor->state_mutex));

    size_t limit = config->max_workers == 0
        ? executor->worker_count : config->max_workers;
    if (limit > executor->worker_count) limit = executor->worker_count;
    size_t worker_count = tinypar_effective_workers(
        config->item_count, config->chunk_size, limit);
    tinypar_job_init(&executor->job, config, worker_count, callback, context);

    tinypar_worker_arg_t caller = { &executor->job, NULL, 0 };
    if (worker_count == 1) {
        tinypar_run_worker(&caller);
    } else {
        tinypar_require_sync(tinypar_mutex_lock(&executor->state_mutex));
        executor->completed_workers = 0;
        executor->active_workers = worker_count - 1;
        executor->generation++;
        for (size_t i = 0; i < executor->active_workers; i++) {
            tinypar_require_sync(
                tinypar_condition_signal(&executor->worker_conditions[i]));
        }
        tinypar_require_sync(tinypar_mutex_unlock(&executor->state_mutex));

        tinypar_run_worker(&caller);

        tinypar_require_sync(tinypar_mutex_lock(&executor->state_mutex));
        while (executor->completed_workers != executor->active_workers) {
            tinypar_require_sync(tinypar_condition_wait(
                &executor->done_condition, &executor->state_mutex));
        }
        tinypar_require_sync(tinypar_mutex_unlock(&executor->state_mutex));
    }

    tinypar_status_t result = (tinypar_status_t)
        tinypar_atomic_int_load(&executor->job.status);
    tinypar_require_sync(tinypar_mutex_unlock(&executor->submission_mutex));
    return result;
}

tinypar_status_t tinypar_executor_destroy(tinypar_executor_t** pointer) {
    if (!pointer || !*pointer) return TINYPAR_INVALID_ARGUMENT;
    tinypar_executor_t* executor = *pointer;

    tinypar_require_sync(tinypar_mutex_lock(&executor->submission_mutex));
    tinypar_require_sync(tinypar_mutex_lock(&executor->state_mutex));
    executor->shutdown = 1;
    for (size_t i = 0; i < executor->spawned_count; i++) {
        tinypar_require_sync(
            tinypar_condition_signal(&executor->worker_conditions[i]));
    }
    tinypar_require_sync(tinypar_mutex_unlock(&executor->state_mutex));

    int unknown = 0;
    for (size_t i = 0; i < executor->spawned_count; i++) {
        if (executor->terminated[i]) continue;
        tinypar_join_result_t joined = tinypar_thread_join(&executor->threads[i]);
        if (joined == TINYPAR_JOIN_TERMINATION_UNKNOWN) {
            unknown = 1;
            continue;
        }
        executor->terminated[i] = 1;
        if (joined == TINYPAR_JOIN_TERMINATED_CLEANUP_FAILED)
            executor->cleanup_failed = 1;
    }

    tinypar_require_sync(tinypar_mutex_unlock(&executor->submission_mutex));
    if (unknown) return TINYPAR_THREAD_JOIN_FAILED;

    int cleanup_failed = executor->cleanup_failed;
    tinypar_executor_destroy_sync(executor);
    free(executor->terminated);
    free(executor->arguments);
    free(executor->threads);
    free(executor->worker_conditions);
    free(executor);
    *pointer = NULL;
    return cleanup_failed ? TINYPAR_THREAD_JOIN_FAILED : TINYPAR_OK;
}

const char* tinypar_status_string(tinypar_status_t status) {
    switch (status) {
    case TINYPAR_OK: return "success";
    case TINYPAR_INVALID_ARGUMENT: return "invalid argument";
    case TINYPAR_ALLOCATION_FAILED: return "allocation failed";
    case TINYPAR_THREAD_CREATE_FAILED: return "thread creation failed";
    case TINYPAR_THREAD_JOIN_FAILED: return "thread join failed";
    case TINYPAR_SYNCHRONIZATION_FAILED: return "synchronization failed";
    case TINYPAR_CALLBACK_FAILED: return "callback failed";
    default: return "unknown tinypar status";
    }
}
