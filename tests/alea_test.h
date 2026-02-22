// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * alea_test.h - Minimal unit test framework for CSG library
 *
 * Single-header, zero dependencies, ~100 lines.
 *
 * Usage:
 *   #define ALEA_TEST_IMPLEMENTATION
 *   #include "alea_test.h"
 *
 *   TEST(my_test) {
 *       ASSERT(1 + 1 == 2);
 *       ASSERT_EQ(42, 42);
 *       ASSERT_NEAR(3.14, 3.14159, 0.01);
 *   }
 *
 *   TEST(another_test) {
 *       ASSERT_STR_EQ("hello", "hello");
 *   }
 *
 *   TEST_MAIN()
 *
 * Compile: gcc -o test test.c -I../include -L../bin -lalea -lm
 */

#ifndef ALEA_TEST_H
#define ALEA_TEST_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

/* ------------------------------------------------------------------------- */
/* Test registration                                                          */
/* ------------------------------------------------------------------------- */

typedef void (*alea_test_fn)(void);

typedef struct alea_test_entry {
    const char *name;
    alea_test_fn fn;
    struct alea_test_entry *next;
} alea_test_entry_t;

/* Global state */
extern alea_test_entry_t *alea_test_list;
extern alea_test_entry_t **alea_test_tail;
extern int alea_test_passed;
extern int alea_test_failed;
extern int alea_test_current_failed;
extern const char *alea_test_current_name;

/* Register a test */
#define TEST(name) \
    static void test_##name(void); \
    static void register_##name(void) __attribute__((constructor)); \
    static alea_test_entry_t entry_##name = { #name, test_##name, NULL }; \
    static void register_##name(void) { \
        *alea_test_tail = &entry_##name; \
        alea_test_tail = &entry_##name.next; \
    } \
    static void test_##name(void)

/* ------------------------------------------------------------------------- */
/* Assertions                                                                  */
/* ------------------------------------------------------------------------- */

#define ASSERT(cond) do { \
    if (!(cond)) { \
        printf("    FAIL: %s:%d: %s\n", __FILE__, __LINE__, #cond); \
        alea_test_current_failed = 1; \
        return; \
    } \
} while(0)

#define ASSERT_MSG(cond, msg) do { \
    if (!(cond)) { \
        printf("    FAIL: %s:%d: %s (%s)\n", __FILE__, __LINE__, #cond, msg); \
        alea_test_current_failed = 1; \
        return; \
    } \
} while(0)

#define ASSERT_EQ(a, b) do { \
    long long _a = (long long)(a), _b = (long long)(b); \
    if (_a != _b) { \
        printf("    FAIL: %s:%d: %s == %s (got %lld, expected %lld)\n", \
               __FILE__, __LINE__, #a, #b, _a, _b); \
        alea_test_current_failed = 1; \
        return; \
    } \
} while(0)

#define ASSERT_NE(a, b) do { \
    long long _a = (long long)(a), _b = (long long)(b); \
    if (_a == _b) { \
        printf("    FAIL: %s:%d: %s != %s (both are %lld)\n", \
               __FILE__, __LINE__, #a, #b, _a); \
        alea_test_current_failed = 1; \
        return; \
    } \
} while(0)

#define ASSERT_NEAR(a, b, eps) do { \
    double _a = (double)(a), _b = (double)(b), _eps = (double)(eps); \
    if (fabs(_a - _b) > _eps) { \
        printf("    FAIL: %s:%d: |%s - %s| <= %s (got %g, expected %g, diff %g)\n", \
               __FILE__, __LINE__, #a, #b, #eps, _a, _b, fabs(_a - _b)); \
        alea_test_current_failed = 1; \
        return; \
    } \
} while(0)

#define ASSERT_STR_EQ(a, b) do { \
    const char *_a = (a), *_b = (b); \
    if (strcmp(_a, _b) != 0) { \
        printf("    FAIL: %s:%d: strcmp(%s, %s) == 0\n", __FILE__, __LINE__, #a, #b); \
        printf("      got:      \"%s\"\n", _a); \
        printf("      expected: \"%s\"\n", _b); \
        alea_test_current_failed = 1; \
        return; \
    } \
} while(0)

#define ASSERT_NULL(ptr) do { \
    if ((ptr) != NULL) { \
        printf("    FAIL: %s:%d: %s == NULL (got %p)\n", \
               __FILE__, __LINE__, #ptr, (void*)(ptr)); \
        alea_test_current_failed = 1; \
        return; \
    } \
} while(0)

#define ASSERT_NOT_NULL(ptr) do { \
    if ((ptr) == NULL) { \
        printf("    FAIL: %s:%d: %s != NULL\n", __FILE__, __LINE__, #ptr); \
        alea_test_current_failed = 1; \
        return; \
    } \
} while(0)

#define ASSERT_TRUE(cond)  ASSERT(cond)
#define ASSERT_FALSE(cond) ASSERT(!(cond))

/* Skip remaining assertions in this test (test still passes) */
#define SKIP(reason) do { \
    printf("    SKIP: %s\n", reason); \
    return; \
} while(0)

/* ------------------------------------------------------------------------- */
/* Test runner                                                                */
/* ------------------------------------------------------------------------- */

#define TEST_MAIN() \
    alea_test_entry_t *alea_test_list = NULL; \
    alea_test_entry_t **alea_test_tail = &alea_test_list; \
    int alea_test_passed = 0; \
    int alea_test_failed = 0; \
    int alea_test_current_failed = 0; \
    const char *alea_test_current_name = NULL; \
    \
    int main(int argc, char **argv) { \
        (void)argc; (void)argv; \
        const char *filter = argc > 1 ? argv[1] : NULL; \
        printf("\n"); \
        for (alea_test_entry_t *t = alea_test_list; t; t = t->next) { \
            if (filter && strstr(t->name, filter) == NULL) continue; \
            alea_test_current_failed = 0; \
            alea_test_current_name = t->name; \
            printf("  %-50s ", t->name); \
            fflush(stdout); \
            t->fn(); \
            if (alea_test_current_failed) { \
                alea_test_failed++; \
            } else { \
                printf("OK\n"); \
                alea_test_passed++; \
            } \
        } \
        printf("\n----------------------------------------\n"); \
        printf("Results: %d passed, %d failed\n", alea_test_passed, alea_test_failed); \
        printf("----------------------------------------\n\n"); \
        return alea_test_failed > 0 ? 1 : 0; \
    }

#endif /* ALEA_TEST_H */
