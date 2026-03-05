// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ARENA_H
#define ARENA_H

#include <stddef.h>
#include <stdbool.h>


/**
 * @file arena.h
 * @brief Fast arena/bump allocator for temporary memory
 * 
 * Arena allocators provide fast O(1) allocation with bulk deallocation.
 * Perfect for temporary data structures that are freed all at once.
 * 
 * Features:
 * - Fast bump-pointer allocation
 * - Automatic growth with chained blocks
 * - Zero-cost bulk reset/free
 * - Cache-friendly sequential allocation
 */

#define ARENA_DEFAULT_ALIGNMENT sizeof(void*)
#define ARENA_DEFAULT_BLOCK_SIZE (1024 * 1024) // 1MB per block

/**
 * A single block in the arena's chain
 */
typedef struct arena_block_s {
    unsigned char* buffer;           // Memory buffer
    size_t total_size;              // Total size of this block
    size_t current_offset;          // Current allocation offset
    struct arena_block_s* next;     // Next block in chain
} arena_block_t;

/**
 * Main arena structure - growable chain of blocks
 */
typedef struct arena {
    arena_block_t* head;            // First block
    arena_block_t* current;         // Current allocation block
    size_t default_block_size;      // Size for new blocks
} arena_t;

/**
 * Initialize an arena with default block size
 * @param arena Arena to initialize
 * @return true on success, false on allocation failure
 */
bool arena_init(arena_t* arena);

/**
 * Initialize arena with custom block size
 * @param arena Arena to initialize
 * @param block_size Size of each block (in bytes)
 * @return true on success, false on allocation failure
 */
bool arena_init_with_size(arena_t* arena, size_t block_size);

/**
 * Free all memory associated with arena
 * @param arena Arena to free
 */
void arena_free(arena_t* arena);

/**
 * Allocate memory from arena with alignment
 * @param arena Arena to allocate from
 * @param size Number of bytes to allocate
 * @param align Alignment requirement (must be power of 2)
 * @return Pointer to allocated memory, or NULL on failure
 */
void* arena_alloc_align(arena_t* arena, size_t size, size_t align);

/**
 * Allocate memory from arena with default alignment
 * @param arena Arena to allocate from
 * @param size Number of bytes to allocate
 * @return Pointer to allocated memory, or NULL on failure
 */
static inline void* arena_alloc(arena_t* arena, size_t size) {
    return arena_alloc_align(arena, size, ARENA_DEFAULT_ALIGNMENT);
}

/**
 * Duplicate a string into arena memory
 * @param arena Arena to allocate from
 * @param str String to duplicate (NULL-terminated)
 * @return Pointer to duplicated string, or NULL on failure
 */
char* arena_strdup(arena_t* arena, const char* str);

/**
 * Reset arena to initial state (keeps allocated blocks for reuse)
 * This is much faster than free + reinit for temporary allocations
 * @param arena Arena to reset
 */
void arena_reset(arena_t* arena);

/**
 * Get total memory allocated by arena (including unused portions)
 * @param arena Arena to query
 * @return Total bytes allocated from system
 */
size_t arena_total_allocated(const arena_t* arena);

/**
 * Get memory currently used in arena
 * @param arena Arena to query
 * @return Bytes currently used
 */
size_t arena_used(const arena_t* arena);


#endif // ARENA_H