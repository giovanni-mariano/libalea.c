// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "alea_alloc.h"
#include "compat.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static ALEA_THREAD_LOCAL long failure_countdown = -1;

void alea_alloc_fail_after(long allocations) {
    failure_countdown = allocations;
}

static int allocation_fails(void) {
    if (failure_countdown < 0) return 0;
    if (failure_countdown-- != 0) return 0;
    return 1;
}

void* alea_malloc(size_t size) {
    return allocation_fails() ? NULL : malloc(size);
}

void* alea_calloc(size_t count, size_t size) {
    if (size && count > SIZE_MAX / size) return NULL;
    return allocation_fails() ? NULL : calloc(count, size);
}

void* alea_realloc(void* pointer, size_t size) {
    return allocation_fails() ? NULL : realloc(pointer, size);
}

char* alea_alloc_strdup(const char* string) {
    if (!string) return NULL;
    size_t size = strlen(string) + 1;
    char* copy = alea_malloc(size);
    if (copy) memcpy(copy, string, size);
    return copy;
}
