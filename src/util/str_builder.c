// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file str_builder.c
 * @brief Generic string builder implementation
 */

#include "str_builder.h"
#include <string.h>
#include <stdarg.h>
#include <stdio.h>

static bool ensure_space(str_builder_t* sb, size_t n)
{
    if (sb->error) return false;
    if (sb->len + n + 1 <= sb->capacity) return true;
    if (!sb->arena) { sb->error = true; return false; }

    size_t new_cap = sb->capacity ? sb->capacity * 2 : 256;
    while (new_cap < sb->len + n + 1) new_cap *= 2;

    char* new_buf = (char*)arena_alloc(sb->arena, new_cap);
    if (!new_buf) { sb->error = true; return false; }

    if (sb->buf && sb->len) memcpy(new_buf, sb->buf, sb->len + 1);
    else new_buf[0] = '\0';

    sb->buf = new_buf;
    sb->capacity = new_cap;
    return true;
}

void str_builder_init(str_builder_t* sb, arena_t* arena, size_t initial_capacity)
{
    sb->arena = arena;
    sb->buf = (char*)arena_alloc(arena, initial_capacity);
    sb->capacity = sb->buf ? initial_capacity : 0;
    sb->len = 0;
    sb->error = false;
    if (sb->buf && sb->capacity > 0) sb->buf[0] = '\0';
}

void str_builder_reset(str_builder_t* sb)
{
    sb->len = 0;
    sb->error = false;
    if (sb->buf && sb->capacity > 0) sb->buf[0] = '\0';
}

bool str_builder_putc(str_builder_t* sb, char c)
{
    if (sb->error) return false;
    if (!ensure_space(sb, 1)) return false;
    sb->buf[sb->len++] = c;
    sb->buf[sb->len] = '\0';
    return true;
}

bool str_builder_puts(str_builder_t* sb, const char* str)
{
    if (sb->error) return false;
    if (!str) return true;
    size_t n = strlen(str);
    if (n == 0) return true;
    if (!ensure_space(sb, n)) return false;
    memcpy(sb->buf + sb->len, str, n);
    sb->len += n;
    sb->buf[sb->len] = '\0';
    return true;
}

bool str_builder_write(str_builder_t* sb, const char* data, size_t n)
{
    if (sb->error) return false;
    if (!data || n == 0) return true;
    if (!ensure_space(sb, n)) return false;
    memcpy(sb->buf + sb->len, data, n);
    sb->len += n;
    sb->buf[sb->len] = '\0';
    return true;
}

bool str_builder_int(str_builder_t* sb, int value)
{
    if (sb->error) return false;
    char tmp[32];
    int n = snprintf(tmp, sizeof(tmp), "%d", value);
    if (n <= 0) return false;
    return str_builder_write(sb, tmp, (size_t)n);
}

bool str_builder_double(str_builder_t* sb, double value, int precision)
{
    if (sb->error) return false;
    char tmp[64];
    int n = snprintf(tmp, sizeof(tmp), "%.*g", precision, value);
    if (n <= 0) return false;
    return str_builder_write(sb, tmp, (size_t)n);
}

bool str_builder_printf(str_builder_t* sb, const char* fmt, ...)
{
    if (sb->error) return false;

    va_list args;
    va_start(args, fmt);

    /* First, determine how much space we need */
    va_list args_copy;
    va_copy(args_copy, args);
    int needed = vsnprintf(NULL, 0, fmt, args_copy);
    va_end(args_copy);

    if (needed < 0) {
        va_end(args);
        return false;
    }

    if (!ensure_space(sb, (size_t)needed)) {
        va_end(args);
        return false;
    }

    int written = vsnprintf(sb->buf + sb->len, sb->capacity - sb->len, fmt, args);
    va_end(args);

    if (written < 0 || written > needed) return false;
    sb->len += (size_t)written;
    return true;
}

size_t str_builder_finish(str_builder_t* sb)
{
    if (sb->buf && sb->len < sb->capacity) sb->buf[sb->len] = '\0';
    return sb->len;
}

bool str_builder_error(const str_builder_t* sb)
{
    return sb->error;
}

const char* str_builder_get(str_builder_t* sb)
{
    str_builder_finish(sb);
    return sb->buf;
}
