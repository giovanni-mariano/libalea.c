// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_arena.c - Arena allocator tests
 */

#include "alea_test.h"
#include "util/arena.h"
#include <string.h>
#include <stdint.h>

TEST(arena_alloc_basic) {
    arena_t arena;
    ASSERT_TRUE(arena_init(&arena));

    void* ptr = arena_alloc(&arena, 64);
    ASSERT_NOT_NULL(ptr);

    /* Write and read back */
    memset(ptr, 0xAB, 64);
    ASSERT_EQ(((unsigned char*)ptr)[0], 0xAB);
    ASSERT_EQ(((unsigned char*)ptr)[63], 0xAB);

    arena_free(&arena);
}

TEST(arena_alloc_multiple) {
    arena_t arena;
    ASSERT_TRUE(arena_init(&arena));

    void* p1 = arena_alloc(&arena, 100);
    void* p2 = arena_alloc(&arena, 200);
    void* p3 = arena_alloc(&arena, 300);

    ASSERT_NOT_NULL(p1);
    ASSERT_NOT_NULL(p2);
    ASSERT_NOT_NULL(p3);

    /* Pointers should be different */
    ASSERT(p1 != p2);
    ASSERT(p2 != p3);

    /* Total used should be >= 600 (with alignment padding) */
    ASSERT(arena_used(&arena) >= 600);

    arena_free(&arena);
}

TEST(arena_reset) {
    arena_t arena;
    ASSERT_TRUE(arena_init(&arena));

    arena_alloc(&arena, 1000);
    ASSERT(arena_used(&arena) >= 1000);

    arena_reset(&arena);
    /* After reset, used should be 0 */
    ASSERT_EQ(arena_used(&arena), 0);

    /* Should be able to allocate again */
    void* ptr = arena_alloc(&arena, 500);
    ASSERT_NOT_NULL(ptr);

    arena_free(&arena);
}

TEST(arena_large_alloc) {
    arena_t arena;
    ASSERT_TRUE(arena_init_with_size(&arena, 1024)); /* Small blocks */

    /* Allocate more than one block */
    void* ptr = arena_alloc(&arena, 2048);
    ASSERT_NOT_NULL(ptr);

    memset(ptr, 0, 2048); /* Should not crash */

    arena_free(&arena);
}

TEST(arena_alignment) {
    arena_t arena;
    ASSERT_TRUE(arena_init(&arena));

    /* Allocate various sizes and check alignment */
    for (int i = 0; i < 10; i++) {
        void* ptr = arena_alloc(&arena, 1 + i * 7);
        ASSERT_NOT_NULL(ptr);
        /* Should be aligned to at least sizeof(void*) */
        ASSERT_EQ((uintptr_t)ptr % sizeof(void*), 0);
    }

    arena_free(&arena);
}

TEST(arena_strdup) {
    arena_t arena;
    ASSERT_TRUE(arena_init(&arena));

    const char* original = "Hello, arena!";
    char* copy = arena_strdup(&arena, original);
    ASSERT_NOT_NULL(copy);
    ASSERT_STR_EQ(copy, original);

    /* Modify copy shouldn't affect original */
    copy[0] = 'X';
    ASSERT(original[0] == 'H');

    arena_free(&arena);
}

TEST_MAIN()
