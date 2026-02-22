// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_hashmap.c - Unit tests for alea_hashmap
 */

#include "alea_test.h"
#include "util/alea_hashmap.h"
#include <limits.h>

/* Define a simple int->int hash map for testing */
static inline uint64_t test_hash_int(int key) {
    return (uint64_t)(unsigned int)key * 0x9E3779B97F4A7C15ULL;
}
static inline bool test_eq_int(int a, int b) { return a == b; }

ALEA_HASHMAP_DEFINE(test_map, int, int, test_hash_int, test_eq_int, INT_MIN);

/* ------------------------------------------------------------------------- */
/* Basic Operations                                                           */
/* ------------------------------------------------------------------------- */

TEST(hashmap_create_destroy) {
    test_map_t map = test_map_create(16);
    ASSERT_NOT_NULL(map.entries);
    ASSERT_EQ(map.count, 0);
    ASSERT(map.capacity >= 16);

    test_map_destroy(&map);
    ASSERT_NULL(map.entries);
    ASSERT_EQ(map.count, 0);
}

TEST(hashmap_put_get) {
    test_map_t map = test_map_create(16);

    ASSERT_TRUE(test_map_put(&map, 42, 100));
    ASSERT_EQ(map.count, 1);

    int* val = test_map_get(&map, 42);
    ASSERT_NOT_NULL(val);
    ASSERT_EQ(*val, 100);

    /* Not found */
    int* missing = test_map_get(&map, 99);
    ASSERT_NULL(missing);

    test_map_destroy(&map);
}

TEST(hashmap_put_overwrite) {
    test_map_t map = test_map_create(16);

    test_map_put(&map, 10, 100);
    test_map_put(&map, 10, 200);

    ASSERT_EQ(map.count, 1);  /* should not increase count */
    int* val = test_map_get(&map, 10);
    ASSERT_NOT_NULL(val);
    ASSERT_EQ(*val, 200);

    test_map_destroy(&map);
}

TEST(hashmap_remove) {
    test_map_t map = test_map_create(16);

    test_map_put(&map, 1, 10);
    test_map_put(&map, 2, 20);
    test_map_put(&map, 3, 30);
    ASSERT_EQ(map.count, 3);

    ASSERT_TRUE(test_map_remove(&map, 2));
    ASSERT_EQ(map.count, 2);
    ASSERT_NULL(test_map_get(&map, 2));

    /* Remaining keys still accessible */
    ASSERT_NOT_NULL(test_map_get(&map, 1));
    ASSERT_NOT_NULL(test_map_get(&map, 3));

    /* Remove non-existent */
    ASSERT_FALSE(test_map_remove(&map, 999));
    ASSERT_EQ(map.count, 2);

    test_map_destroy(&map);
}

TEST(hashmap_clear) {
    test_map_t map = test_map_create(16);

    for (int i = 0; i < 10; i++) {
        test_map_put(&map, i, i * 10);
    }
    ASSERT_EQ(map.count, 10);

    test_map_clear(&map);
    ASSERT_EQ(map.count, 0);

    for (int i = 0; i < 10; i++) {
        ASSERT_NULL(test_map_get(&map, i));
    }

    test_map_destroy(&map);
}

/* ------------------------------------------------------------------------- */
/* Collision and Resize                                                       */
/* ------------------------------------------------------------------------- */

TEST(hashmap_many_entries_trigger_resize) {
    test_map_t map = test_map_create(8);
    size_t initial_cap = map.capacity;

    /* Insert enough entries to trigger resize (load > 70%) */
    for (int i = 0; i < 100; i++) {
        ASSERT_TRUE(test_map_put(&map, i, i * 10));
    }

    ASSERT(map.capacity > initial_cap);
    ASSERT_EQ(map.count, 100);

    /* Verify all entries survived resize */
    for (int i = 0; i < 100; i++) {
        int* val = test_map_get(&map, i);
        ASSERT_NOT_NULL(val);
        ASSERT_EQ(*val, i * 10);
    }

    test_map_destroy(&map);
}

TEST(hashmap_remove_with_collisions) {
    /* Use a small capacity to force collisions */
    test_map_t map = test_map_create(8);

    /* Insert keys that will likely collide in 8-bucket table */
    for (int i = 0; i < 5; i++) {
        test_map_put(&map, i, i * 100);
    }

    /* Remove a middle entry */
    test_map_remove(&map, 2);

    /* All other entries still accessible */
    for (int i = 0; i < 5; i++) {
        int* val = test_map_get(&map, i);
        if (i == 2) {
            ASSERT_NULL(val);
        } else {
            ASSERT_NOT_NULL(val);
            ASSERT_EQ(*val, i * 100);
        }
    }

    test_map_destroy(&map);
}

/* ------------------------------------------------------------------------- */
/* Large Map Stress Test                                                      */
/* ------------------------------------------------------------------------- */

TEST(hashmap_large) {
    test_map_t map = test_map_create(8);

    for (int i = 0; i < 10000; i++) {
        ASSERT_TRUE(test_map_put(&map, i, i));
    }
    ASSERT_EQ(map.count, 10000);

    /* Verify all */
    for (int i = 0; i < 10000; i++) {
        int* val = test_map_get(&map, i);
        ASSERT_NOT_NULL(val);
        ASSERT_EQ(*val, i);
    }

    /* Remove half */
    for (int i = 0; i < 10000; i += 2) {
        ASSERT_TRUE(test_map_remove(&map, i));
    }
    ASSERT_EQ(map.count, 5000);

    /* Verify remaining */
    for (int i = 0; i < 10000; i++) {
        int* val = test_map_get(&map, i);
        if (i % 2 == 0) {
            ASSERT_NULL(val);
        } else {
            ASSERT_NOT_NULL(val);
            ASSERT_EQ(*val, i);
        }
    }

    test_map_destroy(&map);
}

TEST(hashmap_empty_get) {
    test_map_t map = test_map_create(16);
    ASSERT_NULL(test_map_get(&map, 42));
    test_map_destroy(&map);
}

TEST(hashmap_negative_keys) {
    test_map_t map = test_map_create(16);

    test_map_put(&map, -1, 100);
    test_map_put(&map, -100, 200);
    test_map_put(&map, 0, 300);

    int* v1 = test_map_get(&map, -1);
    int* v2 = test_map_get(&map, -100);
    int* v3 = test_map_get(&map, 0);

    ASSERT_NOT_NULL(v1);
    ASSERT_EQ(*v1, 100);
    ASSERT_NOT_NULL(v2);
    ASSERT_EQ(*v2, 200);
    ASSERT_NOT_NULL(v3);
    ASSERT_EQ(*v3, 300);

    test_map_destroy(&map);
}

TEST_MAIN()
