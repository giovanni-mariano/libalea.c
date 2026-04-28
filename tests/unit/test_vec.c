// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_vec.c - Unit tests for alea_vec dynamic array
 */

#include "alea_test.h"
#include "core/alea_system.h"
#include "util/alea_vec.h"

/* Define a test struct */
typedef struct {
    int id;
    double value;
} test_item_t;

ALEA_VEC_DEFINE(test_item_vec, test_item_t);

/* ------------------------------------------------------------------------- */
/* Basic Operations                                                           */
/* ------------------------------------------------------------------------- */

TEST(vec_init_and_free) {
    test_item_vec_t vec = ALEA_VEC_INIT;
    ASSERT_NULL(vec.data);
    ASSERT_EQ(vec.count, 0);
    ASSERT_EQ(vec.capacity, 0);

    alea_vec_free(&vec);
    ASSERT_NULL(vec.data);
    ASSERT_EQ(vec.count, 0);
}

TEST(vec_push_single) {
    test_item_vec_t vec = ALEA_VEC_INIT;
    test_item_t item = { .id = 42, .value = 3.14 };

    size_t index = alea_vec_count(&vec);
    int res = alea_vec_push(&vec, item, test_item_t);
    ASSERT_EQ(res, 0);
    ASSERT_EQ(index, 0);
    ASSERT_EQ(alea_vec_count(&vec), 1);
    ASSERT_EQ(vec.data[0].id, 42);
    ASSERT_NEAR(vec.data[0].value, 3.14, 0.001);

    alea_vec_free(&vec);
}

TEST(vec_push_many) {
    alea_int_vec_t vec = ALEA_VEC_INIT;

    /* Push 100 elements to trigger multiple growths */
    for (int i = 0; i < 100; i++) {
        size_t index = alea_vec_count(&vec);
        int res = alea_vec_push(&vec, i, int);
        ASSERT_EQ(res, 0);
        ASSERT_EQ(index, (size_t)i);
    }

    ASSERT_EQ(alea_vec_count(&vec), 100);
    ASSERT(vec.capacity >= 100);

    /* Verify values */
    for (int i = 0; i < 100; i++) {
        ASSERT_EQ(alea_vec_get(&vec, i), i);
    }

    alea_vec_free(&vec);
}

TEST(vec_push_uninit) {
    test_item_vec_t vec = ALEA_VEC_INIT;

    test_item_t* item = alea_vec_push_uninit(&vec, test_item_t);
    ASSERT_NOT_NULL(item);

    item->id = 99;
    item->value = 2.71;

    ASSERT_EQ(alea_vec_count(&vec), 1);
    ASSERT_EQ(vec.data[0].id, 99);

    alea_vec_free(&vec);
}

/* ------------------------------------------------------------------------- */
/* Capacity Management                                                        */
/* ------------------------------------------------------------------------- */

TEST(vec_reserve) {
    alea_int_vec_t vec = ALEA_VEC_INIT;

    int res = alea_vec_reserve(&vec, 1000, int);
    ASSERT_EQ(res, 0);
    ASSERT(vec.capacity >= 1000);
    ASSERT_EQ(vec.count, 0);

    alea_vec_free(&vec);
}

TEST(vec_clear) {
    alea_int_vec_t vec = ALEA_VEC_INIT;

    for (int i = 0; i < 10; i++) {
        alea_vec_push(&vec, i, int);
    }

    size_t old_capacity = vec.capacity;
    alea_vec_clear(&vec);

    ASSERT_EQ(alea_vec_count(&vec), 0);
    ASSERT_EQ(vec.capacity, old_capacity);
    ASSERT_NOT_NULL(vec.data);

    alea_vec_free(&vec);
}

/* ------------------------------------------------------------------------- */
/* Element Access                                                             */
/* ------------------------------------------------------------------------- */

TEST(vec_last) {
    alea_int_vec_t vec = ALEA_VEC_INIT;

    alea_vec_push(&vec, 10, int);
    alea_vec_push(&vec, 20, int);
    alea_vec_push(&vec, 30, int);

    ASSERT_EQ(alea_vec_last(&vec), 30);

    int* last_ptr = alea_vec_last_ptr(&vec);
    ASSERT_EQ(*last_ptr, 30);
    *last_ptr = 99;
    ASSERT_EQ(alea_vec_last(&vec), 99);

    alea_vec_free(&vec);
}

TEST(vec_pop) {
    alea_int_vec_t vec = ALEA_VEC_INIT;

    alea_vec_push(&vec, 1, int);
    alea_vec_push(&vec, 2, int);
    alea_vec_push(&vec, 3, int);

    ASSERT_EQ(alea_vec_pop(&vec), 3);
    ASSERT_EQ(alea_vec_count(&vec), 2);

    ASSERT_EQ(alea_vec_pop(&vec), 2);
    ASSERT_EQ(alea_vec_count(&vec), 1);

    alea_vec_free(&vec);
}

/* ------------------------------------------------------------------------- */
/* Iteration                                                                  */
/* ------------------------------------------------------------------------- */

TEST(vec_foreach) {
    alea_int_vec_t vec = ALEA_VEC_INIT;

    for (int i = 0; i < 5; i++) {
        alea_vec_push(&vec, i * 10, int);
    }

    int sum = 0;
    alea_vec_foreach(&vec, int, val) {
        sum += *val;
    }

    ASSERT_EQ(sum, 0 + 10 + 20 + 30 + 40);

    alea_vec_free(&vec);
}

TEST(vec_foreach_i) {
    alea_int_vec_t vec = ALEA_VEC_INIT;

    for (int i = 0; i < 5; i++) {
        alea_vec_push(&vec, i * 10, int);
    }

    size_t last_idx = 0;
    alea_vec_foreach_i(&vec, int, i, val) {
        ASSERT_EQ(*val, (int)(i * 10));
        last_idx = i;
    }

    ASSERT_EQ(last_idx, 4);

    alea_vec_free(&vec);
}

/* ------------------------------------------------------------------------- */
/* Queries                                                                    */
/* ------------------------------------------------------------------------- */

TEST(vec_empty) {
    alea_int_vec_t vec = ALEA_VEC_INIT;

    ASSERT_TRUE(alea_vec_empty(&vec));

    alea_vec_push(&vec, 1, int);
    ASSERT_FALSE(alea_vec_empty(&vec));

    alea_vec_clear(&vec);
    ASSERT_TRUE(alea_vec_empty(&vec));

    alea_vec_free(&vec);
}

TEST_MAIN()
