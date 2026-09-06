// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#if !defined(_WIN32)
#define _POSIX_C_SOURCE 200809L
#endif

#include "alea.h"
#include "alea_test.h"
#include "util/alea_parallel.h"

#include <stdlib.h>

static int count_range(void* context, size_t worker, size_t begin, size_t end) {
    (void)worker;
    size_t* visits = context;
    for (size_t item = begin; item < end; item++) visits[item]++;
    return 0;
}

TEST(parallel_environment_sets_process_default) {
#if defined(_WIN32)
    ASSERT_EQ(_putenv_s("ALEA_NUM_THREADS", "3"), 0);
#else
    ASSERT_EQ(setenv("ALEA_NUM_THREADS", "3", 1), 0);
#endif

    const int enabled = alea_parallel_enabled();
    ASSERT_EQ(alea_parallel_max_threads(), enabled ? 3 : 1);

    size_t visits[12] = {0};
    size_t actual_workers = 0;
    ASSERT_EQ(alea_parallel_for(
        12, 1, 0, ALEA_PARALLEL_DYNAMIC, count_range, visits,
        &actual_workers), ALEA_PARALLEL_OK);
    ASSERT_EQ(actual_workers, (size_t)(enabled ? 3 : 1));
    for (size_t item = 0; item < 12; item++) ASSERT_EQ(visits[item], 1);
}

TEST_MAIN()
