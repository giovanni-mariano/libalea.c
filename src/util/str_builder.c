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
#include <stdlib.h>
#include <stdint.h>  /* SIZE_MAX */

/* Compute a buffer capacity >= needed using ~2x growth. `needed` must already
 * be a representable size (<= SIZE_MAX); the caller guarantees this. Starting
 * from at least 1 prevents the doubling loop from stalling on a zero capacity,
 * and the SIZE_MAX/2 guard prevents it from overflowing. */
static size_t grow_capacity(size_t current, size_t needed)
{
    size_t cap = current ? current : 1;
    while (cap < needed) {
        if (cap > SIZE_MAX / 2) return needed;  /* doubling would overflow */
        cap *= 2;
    }
    return cap;
}

static bool ensure_space(str_builder_t* sb, size_t n)
{
    if (sb->error) return false;

    /* Guard the `len + n + 1` arithmetic against size_t overflow before it is
     * used as an allocation size. */
    if (n > SIZE_MAX - 1 - sb->len) { sb->error = true; return false; }
    size_t needed = sb->len + n + 1;
    if (needed <= sb->capacity) return true;

    /* Stream mode: flush current buffer to file and reuse it. */
    if (sb->stream) {
        if (sb->len > 0) {
            size_t wrote = fwrite(sb->buf, 1, sb->len, sb->stream);
            if (wrote != sb->len) { sb->error = true; return false; }
            sb->len = 0;
            if (sb->capacity > 0) sb->buf[0] = '\0';
        }
        /* After flush the buffer is empty; if n still doesn't fit (single
         * write larger than the buffer), grow it once with realloc. */
        if (n + 1 <= sb->capacity) return true;
        size_t new_cap = grow_capacity(sb->capacity, n + 1);
        char* new_buf = (char*)realloc(sb->buf, new_cap);
        if (!new_buf) { sb->error = true; return false; }
        sb->buf = new_buf;
        sb->capacity = new_cap;
        return true;
    }

    /* Arena mode: allocate a larger block (old block stays in arena). */
    if (!sb->arena) { sb->error = true; return false; }

    size_t new_cap = grow_capacity(sb->capacity, needed);

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
    sb->stream = NULL;   /* must be NULL so ensure_space uses arena path */
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

void str_builder_init_stream(str_builder_t* sb, size_t buf_size, FILE* stream)
{
    memset(sb, 0, sizeof(*sb));
    sb->stream = stream;
    /* A zero-size buffer would make malloc(0) ambiguous and the buf[0]
     * terminator write below an overflow; reject it as an error. */
    if (buf_size == 0) { sb->error = true; return; }
    sb->buf = (char*)malloc(buf_size);
    if (sb->buf) {
        sb->capacity = buf_size;
        sb->buf[0] = '\0';
    } else {
        sb->error = true;
    }
}

bool str_builder_flush(str_builder_t* sb)
{
    if (!sb->stream || sb->len == 0) return !sb->error;
    size_t wrote = fwrite(sb->buf, 1, sb->len, sb->stream);
    if (wrote != sb->len) { sb->error = true; return false; }
    sb->len = 0;
    if (sb->capacity > 0) sb->buf[0] = '\0';
    return true;
}

void str_builder_destroy(str_builder_t* sb)
{
    if (sb->stream) {
        free(sb->buf);
        sb->buf = NULL;
        sb->capacity = 0;
        sb->len = 0;
    }
    /* Arena-backed builders: the arena itself is freed by the caller. */
}
