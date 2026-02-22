// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "arena.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "util/alea_log.h"

/**
 * @file arena.c
 * @brief Implementation of fast arena allocator
 * 
 */

// Helper: Create a new arena block
static arena_block_t* arena_new_block(size_t size) {
    arena_block_t* block = (arena_block_t*)malloc(sizeof(arena_block_t));
    if (!block) {
        return NULL;
    }

    block->buffer = (unsigned char*)malloc(size);
    if (!block->buffer) {
        free(block);
        return NULL;
    }
    
    block->total_size = size;
    block->current_offset = 0;
    block->next = NULL;
    return block;
}

// Helper: Align an offset forward to alignment boundary
static inline size_t align_forward(size_t ptr, size_t align) {
    return (ptr + (align - 1)) & ~(align - 1);
}

bool arena_init(arena_t* arena) {
    return arena_init_with_size(arena, ARENA_DEFAULT_BLOCK_SIZE);
}

bool arena_init_with_size(arena_t* arena, size_t block_size) {
    if (!arena || block_size == 0) {
        return false;
    }
    
    arena->default_block_size = block_size;
    arena->head = arena_new_block(block_size);
    arena->current = arena->head;
    
    return arena->head != NULL;
}

void arena_free(arena_t* arena) {
    if (!arena) {
        return;
    }
    
    arena_block_t* block = arena->head;
    while (block) {
        arena_block_t* next = block->next;
        free(block->buffer);
        free(block);
        block = next;
    }
    
    arena->head = NULL;
    arena->current = NULL;
}

void* arena_alloc_align(arena_t* arena, size_t size, size_t align) {
    if (!arena || !arena->current || size == 0) {
        return NULL;
    }

    // Try to allocate from the current block
    size_t aligned_offset = align_forward(arena->current->current_offset, align);
    
    if (aligned_offset + size <= arena->current->total_size) {
        // Allocation fits in current block
        void* ptr = arena->current->buffer + aligned_offset;
        arena->current->current_offset = aligned_offset + size;
        memset(ptr, 0, size);
        return ptr;
    }

    // Current block is full - need a new one
    // Make sure the new block is large enough for the requested size
    size_t new_block_size = arena->default_block_size;
    if (size > new_block_size) {
        // For large allocations, allocate a bigger block
        new_block_size = size * 2;
    }
    
    arena_block_t* new_block = arena_new_block(new_block_size);
    if (!new_block) {
        ALEA_LOG_ERROR("Arena failed to allocate new block of %zu bytes",
                new_block_size);
        return NULL;
    }
    
    // Link new block to chain
    arena->current->next = new_block;
    arena->current = new_block;

    // Now retry allocation in the new block (guaranteed to succeed)
    aligned_offset = align_forward(arena->current->current_offset, align);
    void* ptr = arena->current->buffer + aligned_offset;
    arena->current->current_offset = aligned_offset + size;
    memset(ptr, 0, size);
    return ptr;
}

char* arena_strdup(arena_t* arena, const char* str) {
    if (!str) {
        return NULL;
    }
    
    size_t len = strlen(str);
    char* dest = (char*)arena_alloc(arena, len + 1);
    if (dest) {
        memcpy(dest, str, len);
        dest[len] = '\0';
    }
    return dest;
}

void arena_reset(arena_t* arena) {
    if (!arena) {
        return;
    }
    
    // Reset all blocks to unused state
    arena_block_t* block = arena->head;
    while (block) {
        block->current_offset = 0;
        block = block->next;
    }
    
    // Reset current pointer to head
    arena->current = arena->head;
}

size_t arena_total_allocated(const arena_t* arena) {
    if (!arena) {
        return 0;
    }
    
    size_t total = 0;
    const arena_block_t* block = arena->head;
    while (block) {
        total += block->total_size;
        block = block->next;
    }
    return total;
}

size_t arena_used(const arena_t* arena) {
    if (!arena) {
        return 0;
    }
    
    size_t used = 0;
    const arena_block_t* block = arena->head;
    while (block) {
        used += block->current_offset;
        block = block->next;
    }
    return used;
}