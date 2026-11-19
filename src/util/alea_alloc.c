// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "alea_alloc.h"
#include "compat.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

static ALEA_THREAD_LOCAL long failure_countdown = -1;
static ALEA_THREAD_LOCAL alea_alloc_budget_t* active_budget;

void alea_alloc_budget_begin(alea_alloc_budget_t* budget, size_t limit) {
    *budget = (alea_alloc_budget_t){.limit=limit, .previous=active_budget};
    active_budget = budget;
}

void alea_alloc_budget_end(alea_alloc_budget_t* budget) {
    assert(active_budget == budget);
    active_budget = budget->previous;
}

static void budget_failed(void) {
    for (alea_alloc_budget_t* b=active_budget; b; b=b->previous) b->failed=1;
}

static int budget_charge(size_t size) {
    for (alea_alloc_budget_t* b=active_budget; b; b=b->previous) {
        if (b->failed || size > b->limit-b->requested) {
            if (!b->failed) b->exceeded=1;
            budget_failed();
            return -1;
        }
    }
    for (alea_alloc_budget_t* b=active_budget; b; b=b->previous) b->requested+=size;
    return 0;
}

void alea_alloc_fail_after(long allocations) {
    failure_countdown = allocations;
}

static int allocation_fails(void) {
    if (failure_countdown < 0) return 0;
    if (failure_countdown-- != 0) return 0;
    return 1;
}

void* alea_malloc(size_t size) {
    if (budget_charge(size)) return NULL;
    void* p=allocation_fails() ? NULL : malloc(size);
    if (!p && size) budget_failed();
    return p;
}

void* alea_calloc(size_t count, size_t size) {
    if (size && count > SIZE_MAX / size) { budget_failed(); return NULL; }
    if (budget_charge(count*size)) return NULL;
    void* p=allocation_fails() ? NULL : calloc(count, size);
    if (!p && count && size) budget_failed();
    return p;
}

void* alea_realloc(void* pointer, size_t size) {
    if (budget_charge(size)) return NULL;
    void* p=allocation_fails() ? NULL : realloc(pointer, size);
    if (!p && size) budget_failed();
    return p;
}

char* alea_alloc_strdup(const char* string) {
    if (!string) return NULL;
    size_t size = strlen(string) + 1;
    char* copy = alea_malloc(size);
    if (copy) memcpy(copy, string, size);
    return copy;
}
