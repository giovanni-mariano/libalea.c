// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_ALLOC_H
#define ALEA_ALLOC_H

#include <stddef.h>

/* Internal allocation helpers. Failure injection is local to the calling
 * thread; -1 disables it, 0 fails the next allocation, N skips N allocations.
 * A failure automatically disables injection so cleanup and retry can proceed. */
void alea_alloc_fail_after(long allocations);
void* alea_malloc(size_t size);
void* alea_calloc(size_t count, size_t size);
void* alea_realloc(void* pointer, size_t size);
char* alea_alloc_strdup(const char* string);

/* Thread-local scopes bound cumulative requested bytes, including full realloc
 * capacities and temporary allocations. Freeing storage does not refund bytes.
 * Nested scopes charge all active budgets. These are allocation limits, not RSS
 * limits: allocator overhead, stack, stdio buffers and mappings are excluded. */
typedef struct alea_alloc_budget {
    size_t limit, requested;
    int failed, exceeded;
    struct alea_alloc_budget* previous;
} alea_alloc_budget_t;
void alea_alloc_budget_begin(alea_alloc_budget_t* budget, size_t limit);
void alea_alloc_budget_end(alea_alloc_budget_t* budget);

#endif
