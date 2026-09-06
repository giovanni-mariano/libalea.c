// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_test.h"
#include "alea.h"
#include "util/alea_parallel.h"

#include <string.h>

#if !defined(_WIN32)
#include <pthread.h>
#include <signal.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

typedef struct count_context {
    size_t counts[4];
} count_context_t;

static int count_range(void* opaque, size_t worker_index,
                       size_t begin, size_t end) {
    count_context_t* context = opaque;
    if (worker_index >= 4) return -1;
    context->counts[worker_index] += end - begin;
    return 0;
}

static size_t count_total(const count_context_t* context) {
    return context->counts[0] + context->counts[1] +
           context->counts[2] + context->counts[3];
}

TEST(runtime_reuses_bounded_executor) {
    int configured = alea_parallel_set_threads(4);
    ASSERT(configured == 0 || configured == -1);
    size_t worker_limit = alea_parallel_enabled() ? 4u : 1u;
    size_t expected_workers = 0;

    count_context_t context = {{0}};
    for (size_t iteration = 0; iteration < 10000; iteration++) {
        size_t actual_workers = 0;
        ASSERT_EQ(alea_parallel_for(
            64, 1, 4, ALEA_PARALLEL_DYNAMIC, count_range, &context,
            &actual_workers), ALEA_PARALLEL_OK);
        if (iteration == 0) {
            ASSERT(actual_workers >= 1);
            ASSERT(actual_workers <= worker_limit);
            expected_workers = actual_workers;
        }
        ASSERT_EQ(actual_workers, expected_workers);
    }
    ASSERT_EQ(count_total(&context), 640000);
    ASSERT_EQ(alea_parallel_set_threads(2),
              alea_parallel_enabled() ? -1 : 0);
}

typedef struct nested_context {
    unsigned char visits[4][17];
    size_t inner_workers[4];
    int failed[4];
} nested_context_t;

static int mark_nested(void* opaque, size_t worker_index,
                       size_t begin, size_t end) {
    unsigned char* visits = opaque;
    if (worker_index != 0) return -1;
    for (size_t i = begin; i < end; i++) visits[i]++;
    return 0;
}

static int run_nested(void* opaque, size_t worker_index,
                      size_t begin, size_t end) {
    nested_context_t* context = opaque;
    (void)begin;
    (void)end;
    if (worker_index >= 4) return -1;
    alea_parallel_status_t status = alea_parallel_for(
        17, 1, 4, ALEA_PARALLEL_DYNAMIC, mark_nested,
        context->visits[worker_index], &context->inner_workers[worker_index]);
    if (status != ALEA_PARALLEL_OK) {
        context->failed[worker_index] = 1;
        return -1;
    }
    return 0;
}

TEST(runtime_nested_submission_is_serial) {
    nested_context_t context;
    memset(&context, 0, sizeof(context));
    size_t outer_workers = 0;
    ASSERT_EQ(alea_parallel_for(
        4, 1, 4, ALEA_PARALLEL_STATIC_BLOCK, run_nested, &context,
        &outer_workers), ALEA_PARALLEL_OK);
    for (size_t worker = 0; worker < outer_workers; worker++) {
        ASSERT_EQ(context.failed[worker], 0);
        ASSERT_EQ(context.inner_workers[worker], 1);
        for (size_t item = 0; item < 17; item++)
            ASSERT_EQ(context.visits[worker][item], 1);
    }

    memset(&context, 0, sizeof(context));
    ASSERT_EQ(alea_parallel_for(
        1, 1, 1, ALEA_PARALLEL_STATIC_BLOCK, run_nested, &context,
        &outer_workers), ALEA_PARALLEL_OK);
    ASSERT_EQ(outer_workers, 1);
    ASSERT_EQ(context.inner_workers[0], 1);
    for (size_t item = 0; item < 17; item++)
        ASSERT_EQ(context.visits[0][item], 1);
}

#if !defined(_WIN32)
typedef struct submit_context {
    count_context_t count;
    alea_parallel_status_t status;
} submit_context_t;

static void* submit_job(void* opaque) {
    submit_context_t* context = opaque;
    size_t actual_workers = 0;
    context->status = alea_parallel_for(
        10000, 17, 4, ALEA_PARALLEL_DYNAMIC, count_range, &context->count,
        &actual_workers);
    return NULL;
}

TEST(runtime_serializes_concurrent_submitters) {
    if (!alea_parallel_enabled()) SKIP("native worker backend disabled");

    submit_context_t contexts[2] = {
        { {{0}}, ALEA_PARALLEL_INVALID_ARGUMENT },
        { {{0}}, ALEA_PARALLEL_INVALID_ARGUMENT }
    };
    pthread_t threads[2];
    ASSERT_EQ(pthread_create(&threads[0], NULL, submit_job, &contexts[0]), 0);
    ASSERT_EQ(pthread_create(&threads[1], NULL, submit_job, &contexts[1]), 0);
    ASSERT_EQ(pthread_join(threads[0], NULL), 0);
    ASSERT_EQ(pthread_join(threads[1], NULL), 0);
    for (size_t i = 0; i < 2; i++) {
        ASSERT_EQ(contexts[i].status, ALEA_PARALLEL_OK);
        ASSERT_EQ(count_total(&contexts[i].count), 10000);
    }
}

TEST(runtime_recreates_executor_in_post_fork_child) {
    if (!alea_parallel_enabled()) SKIP("native worker backend disabled");

    count_context_t parent_context = {{0}};
    ASSERT_EQ(alea_parallel_for(
        64, 1, 4, ALEA_PARALLEL_DYNAMIC, count_range, &parent_context, NULL),
        ALEA_PARALLEL_OK);
    ASSERT_EQ(count_total(&parent_context), 64);

    pid_t child = fork();
    ASSERT(child >= 0);
    if (child == 0) {
        alarm(5);
        count_context_t child_context = {{0}};
        alea_parallel_status_t status = alea_parallel_for(
            64, 1, 4, ALEA_PARALLEL_DYNAMIC,
            count_range, &child_context, NULL);
        _exit(status == ALEA_PARALLEL_OK &&
              count_total(&child_context) == 64 ? 0 : 1);
    }

    int status = 0;
    ASSERT_EQ(waitpid(child, &status, 0), child);
    ASSERT(WIFEXITED(status));
    ASSERT_EQ(WEXITSTATUS(status), 0);

    memset(&parent_context, 0, sizeof(parent_context));
    ASSERT_EQ(alea_parallel_for(
        64, 1, 4, ALEA_PARALLEL_DYNAMIC, count_range, &parent_context, NULL),
        ALEA_PARALLEL_OK);
    ASSERT_EQ(count_total(&parent_context), 64);
}
#endif

TEST_MAIN()
