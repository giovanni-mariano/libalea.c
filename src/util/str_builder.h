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
    arena_t* arena;         /* Arena for allocations */
    bool error;             /* Sticky error flag */
} str_builder_t;

/**
 * @brief Initialize string builder
 */
void str_builder_init(str_builder_t* sb, arena_t* arena, size_t initial_capacity);

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
