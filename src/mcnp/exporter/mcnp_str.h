// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file mcnp_str.h
 * @brief MCNP string builder with 80-column formatting
 */

#ifndef MCNP_STR_H
#define MCNP_STR_H

#include <stddef.h>
#include <stdbool.h>
#include <stdio.h>
#include "util/str_builder.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief MCNP string builder
 *
 * Builds strings respecting MCNP 80-column format:
 * - First line starts at column 1
 * - Continuation lines start with 5 spaces
 * - Line breaks at spaces, colons, or parentheses when possible
 *
 * Delegates buffer management to an embedded str_builder_t.
 */
typedef struct {
    str_builder_t sb;       /* Underlying buffer management */
    int col;                /* Current column (1-based) */
    int max_col;            /* Max column (default 80) */
    int cont_indent;        /* Continuation indent (default 5) */
} mcnp_str_t;

/**
 * @brief Initialize string builder with external buffer
 */
void mcnp_str_init(mcnp_str_t* s, arena_t* arena, size_t initial_capacity,
                    int max_col, int cont_indent);

/**
 * @brief Initialize unwrapped string builder (no line breaks)
 */
void mcnp_str_init_unwrapped(mcnp_str_t* s, arena_t* arena, size_t initial_capacity);

/**
 * @brief Reset builder for reuse (keeps buffer)
 */
void mcnp_str_reset(mcnp_str_t* s);

/**
 * @brief Append a character
 */
bool mcnp_str_putc(mcnp_str_t* s, char c);

/**
 * @brief Append a token (word) that must not be split across lines
 *
 * If token won't fit on current line, breaks BEFORE it.
 */
bool mcnp_str_token(mcnp_str_t* s, const char* token);

/**
 * @brief Append a string (may contain spaces, handles word breaks)
 */
bool mcnp_str_puts(mcnp_str_t* s, const char* str);

/**
 * @brief Append formatted string (like sprintf)
 */
bool mcnp_str_printf(mcnp_str_t* s, const char* fmt, ...);

/**
 * @brief Append an integer (handles negative sign properly)
 */
bool mcnp_str_int(mcnp_str_t* s, int value);

/**
 * @brief Append a double with specified precision
 */
bool mcnp_str_double(mcnp_str_t* s, double value, int precision);

/**
 * @brief Force a line break (continuation)
 */
bool mcnp_str_newline(mcnp_str_t* s);

/**
 * @brief Write a full-line comment (c text...)
 *
 * Wraps long comments across multiple "c " lines.
 * Starts a new line if not already at column 1.
 */
bool mcnp_str_comment(mcnp_str_t* s, const char* text);

/**
 * @brief Append an inline comment ($ text...)
 *
 * Appends " $ text" at the current position on the current line.
 */
bool mcnp_str_inline_comment(mcnp_str_t* s, const char* text);

/**
 * @brief Ensure null termination and return length
 */
size_t mcnp_str_finish(mcnp_str_t* s);

/**
 * @brief Check if builder is in error state
 */
bool mcnp_str_error(const mcnp_str_t* s);

/**
 * @brief Get current string (null-terminated)
 */
const char* mcnp_str_get(mcnp_str_t* s);

/**
 * @brief Write completed string to file with newline
 */
bool mcnp_str_write(mcnp_str_t* s, FILE* out);

#ifdef __cplusplus
}
#endif

#endif /* MCNP_STR_H */
