// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_test.h"
#include "tinypar.h"
#include "util/alea_parallel.h"

typedef struct nested_first_use_context {
    unsigned char visits[17];
    size_t inner_workers;
} nested_first_use_context_t;

static int mark_inner(void* opaque, size_t worker_index,
                      size_t begin, size_t end) {
    nested_first_use_context_t* context = opaque;
    if (worker_index != 0) return -1;
    for (size_t item = begin; item < end; item++) context->visits[item]++;
    return 0;
}

static int call_libalea_from_tinypar(void* opaque, size_t worker_index,
                                     size_t begin, size_t end) {
    nested_first_use_context_t* context = opaque;
    (void)begin;
    (void)end;
    if (worker_index != 0 || !tinypar_in_callback()) return -1;
    return alea_parallel_for(
        17, 1, 4, ALEA_PARALLEL_DYNAMIC, mark_inner, context,
        &context->inner_workers) == ALEA_PARALLEL_OK ? 0 : -1;
}

TEST(nested_first_use_does_not_create_process_executor) {
    nested_first_use_context_t context = {{0}, 0};
    tinypar_config_t outer = {
        1, 1, 1, TINYPAR_SCHEDULE_DYNAMIC
    };
    ASSERT_EQ(tinypar_parallel_for(
        &outer, call_libalea_from_tinypar, &context), TINYPAR_OK);
    ASSERT_EQ(context.inner_workers, 1);
    for (size_t item = 0; item < 17; item++)
        ASSERT_EQ(context.visits[item], 1);

    /* Configuration remains possible because the nested call did not create
     * the process executor. */
    ASSERT_EQ(alea_parallel_set_default_workers(2), ALEA_PARALLEL_OK);
}

TEST_MAIN()
