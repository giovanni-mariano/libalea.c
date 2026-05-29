// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_str_builder.c - String builder tests
 */

#include "alea_test.h"
#include "util/str_builder.h"
#include "util/arena.h"
#include <string.h>
#include <stdio.h>

TEST(str_basic) {
    arena_t arena;
    ASSERT_TRUE(arena_init(&arena));

    str_builder_t sb;
    str_builder_init(&sb, &arena, 64);

    str_builder_puts(&sb, "Hello");
    str_builder_puts(&sb, ", ");
    str_builder_puts(&sb, "world!");
    str_builder_finish(&sb);

    ASSERT_STR_EQ(str_builder_get(&sb), "Hello, world!");
    ASSERT_FALSE(str_builder_error(&sb));

    arena_free(&arena);
}

TEST(str_printf) {
    arena_t arena;
    ASSERT_TRUE(arena_init(&arena));

    str_builder_t sb;
    str_builder_init(&sb, &arena, 64);

    str_builder_printf(&sb, "Value: %d, Float: %.2f", 42, 3.14);
    str_builder_finish(&sb);

    ASSERT_STR_EQ(str_builder_get(&sb), "Value: 42, Float: 3.14");

    arena_free(&arena);
}

TEST(str_grow) {
    arena_t arena;
    ASSERT_TRUE(arena_init(&arena));

    str_builder_t sb;
    str_builder_init(&sb, &arena, 8); /* Very small initial capacity */

    /* Append enough to force growth */
    for (int i = 0; i < 100; i++) {
        str_builder_puts(&sb, "X");
    }
    str_builder_finish(&sb);

    ASSERT_EQ(strlen(str_builder_get(&sb)), 100);
    ASSERT_FALSE(str_builder_error(&sb));

    arena_free(&arena);
}

TEST(str_empty) {
    arena_t arena;
    ASSERT_TRUE(arena_init(&arena));

    str_builder_t sb;
    str_builder_init(&sb, &arena, 64);
    str_builder_finish(&sb);

    const char* s = str_builder_get(&sb);
    ASSERT_NOT_NULL(s);
    ASSERT_EQ(strlen(s), 0);

    arena_free(&arena);
}

TEST(str_long) {
    arena_t arena;
    ASSERT_TRUE(arena_init(&arena));

    str_builder_t sb;
    str_builder_init(&sb, &arena, 64);

    /* Build a >4KB string */
    for (int i = 0; i < 5000; i++) {
        str_builder_putc(&sb, 'A');
    }
    str_builder_finish(&sb);

    ASSERT_EQ(strlen(str_builder_get(&sb)), 5000);
    ASSERT_FALSE(str_builder_error(&sb));

    arena_free(&arena);
}

TEST(str_int) {
    arena_t arena;
    ASSERT_TRUE(arena_init(&arena));

    str_builder_t sb;
    str_builder_init(&sb, &arena, 64);

    str_builder_int(&sb, 42);
    str_builder_puts(&sb, " ");
    str_builder_int(&sb, -7);
    str_builder_finish(&sb);

    ASSERT_STR_EQ(str_builder_get(&sb), "42 -7");

    arena_free(&arena);
}

TEST(str_double) {
    arena_t arena;
    ASSERT_TRUE(arena_init(&arena));

    str_builder_t sb;
    str_builder_init(&sb, &arena, 64);

    /* precision is significant digits (%g format) */
    str_builder_double(&sb, 3.14159, 4);
    str_builder_finish(&sb);

    ASSERT_STR_EQ(str_builder_get(&sb), "3.142");

    arena_free(&arena);
}

TEST(str_reset) {
    arena_t arena;
    ASSERT_TRUE(arena_init(&arena));

    str_builder_t sb;
    str_builder_init(&sb, &arena, 64);

    str_builder_puts(&sb, "first");
    str_builder_finish(&sb);
    ASSERT_STR_EQ(str_builder_get(&sb), "first");

    str_builder_reset(&sb);
    str_builder_puts(&sb, "second");
    str_builder_finish(&sb);
    ASSERT_STR_EQ(str_builder_get(&sb), "second");

    arena_free(&arena);
}

/* Regression: a stream builder initialised with zero capacity must report an
   error rather than leaving a zero-capacity buffer that loops forever (or
   overflows a malloc(0) allocation) on the next write. */
TEST(str_stream_zero_capacity_errors) {
    FILE* f = tmpfile();
    ASSERT_NOT_NULL(f);

    str_builder_t sb;
    str_builder_init_stream(&sb, 0, f);
    ASSERT_TRUE(str_builder_error(&sb));

    /* Writes are no-ops in the error state and must not hang or crash. */
    ASSERT_FALSE(str_builder_puts(&sb, "data"));

    str_builder_destroy(&sb);
    fclose(f);
}

TEST_MAIN()
