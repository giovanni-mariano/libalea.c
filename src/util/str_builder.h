// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file str_builder.h
 * @brief Generic string builder with automatic growth
 *
 * A simple, safe string builder that uses arena allocation.
 * Can be used by any exporter (MCNP, OpenMC, etc.)
 */

#ifndef STR_BUILDER_H
#define STR_BUILDER_H

#include <stddef.h>
#include <stdbool.h>
#include <stdio.h>
#include "arena.h"


typedef struct {
    char* buf;              /* Output buffer */
    size_t capacity;        /* Buffer capacity */
    size_t len;             /* Current length (not including null) */
    arena_t* arena;         /* Arena for allocations (NULL in stream mode) */
    bool error;             /* Sticky error flag */
    /* Stream mode: when non-NULL the buffer is malloc'd and flushed to this
     * file instead of growing via the arena.  This keeps memory constant
     * regardless of total output size.  arena must be NULL in this mode. */
    FILE* stream;
} str_builder_t;

/**
 * @brief Initialize string builder (arena-backed, grows without bound)
 */
void str_builder_init(str_builder_t* sb, arena_t* arena, size_t initial_capacity);

/**
 * @brief Initialize string builder in streaming mode.
 *
 * Uses a fixed malloc'd buffer of @p buf_size bytes.  Whenever the buffer
 * would overflow, its contents are flushed to @p stream and the buffer is
 * reused from the start — so memory stays bounded at @p buf_size bytes
 * regardless of total output size.
 *
 * Call str_builder_flush() after the last write to drain any remaining bytes,
 * then str_builder_destroy() to free the buffer.
 */
void str_builder_init_stream(str_builder_t* sb, size_t buf_size, FILE* stream);

/**
 * @brief Flush any buffered bytes to the underlying stream (stream mode only).
 * @return true on success, false on write error.
 */
bool str_builder_flush(str_builder_t* sb);

/**
 * @brief Free the malloc'd buffer (stream mode only).  No-op for arena mode.
 */
void str_builder_destroy(str_builder_t* sb);

/**
 * @brief Reset builder for reuse (keeps buffer)
 */
void str_builder_reset(str_builder_t* sb);

/**
 * @brief Append a single character
 */
bool str_builder_putc(str_builder_t* sb, char c);

/**
 * @brief Append a null-terminated string
 */
bool str_builder_puts(str_builder_t* sb, const char* str);

/**
 * @brief Append n bytes from a buffer
 */
bool str_builder_write(str_builder_t* sb, const char* data, size_t n);

/**
 * @brief Append an integer
 */
bool str_builder_int(str_builder_t* sb, int value);

/**
 * @brief Append a double with specified precision
 */
bool str_builder_double(str_builder_t* sb, double value, int precision);

/**
 * @brief Append formatted string (like sprintf)
 */
bool str_builder_printf(str_builder_t* sb, const char* fmt, ...);

/**
 * @brief Ensure null termination and return length
 */
size_t str_builder_finish(str_builder_t* sb);

/**
 * @brief Check if builder is in error state
 */
bool str_builder_error(const str_builder_t* sb);

/**
 * @brief Get current string (null-terminated)
 */
const char* str_builder_get(str_builder_t* sb);


#endif /* STR_BUILDER_H */
