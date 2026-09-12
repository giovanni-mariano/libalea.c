// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "nuclear_internal.h"

#include <stdint.h>

static bool (*allocation_should_fail)(void* context);
static void* allocation_failure_context;

void alea_nuc_set_alloc_failure(bool (*should_fail)(void*), void* context) {
    allocation_should_fail = should_fail;
    allocation_failure_context = context;
}

static bool fail_allocation(void) {
    return allocation_should_fail &&
           allocation_should_fail(allocation_failure_context);
}

void* alea_nuc_malloc(size_t size) {
    if (fail_allocation()) return NULL;
    return malloc(size);
}

void* alea_nuc_calloc(size_t count, size_t size) {
    if (size && count > SIZE_MAX / size) return NULL;
    if (fail_allocation()) return NULL;
    return calloc(count, size);
}

void* alea_nuc_realloc(void* pointer, size_t size) {
    if (fail_allocation()) return NULL;
    return realloc(pointer, size);
}

char* alea_nuc_strdup(const char* text) {
    if (!text) return NULL;
    size_t length = strlen(text);
    if (length == SIZE_MAX || fail_allocation()) return NULL;
    char* copy = malloc(length + 1);
    if (copy) memcpy(copy, text, length + 1);
    return copy;
}
