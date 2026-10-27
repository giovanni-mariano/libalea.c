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

#endif
