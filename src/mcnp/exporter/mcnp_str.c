// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file mcnp_str.c
 * @brief MCNP string builder implementation
 *
 * Key requirements:
 * - Never split "words" (numbers, tokens) across lines
 * - Can break after space or colon (:)
 * - If a token won't fit on current line, break BEFORE it
 * - Maximum 80 columns per line
 * - Continuation lines start with 5 spaces
 *
 * Buffer management is delegated to str_builder_t.
 */

#include "mcnp_str.h"
#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include <limits.h>

void mcnp_str_init(mcnp_str_t* s, arena_t* arena, size_t initial_capacity,
                    int max_col, int cont_indent) {
    str_builder_init(&s->sb, arena, initial_capacity);
    s->col = 1;
    s->max_col = max_col;
    s->cont_indent = cont_indent;
}

void mcnp_str_init_unwrapped(mcnp_str_t* s, arena_t* arena, size_t initial_capacity) {
    mcnp_str_init(s, arena, initial_capacity, INT_MAX, 0);
}

void mcnp_str_reset(mcnp_str_t* s) {
    str_builder_reset(&s->sb);
    s->col = 1;
}

/* Internal: check if character is a valid break point */
static bool is_break_char(char c) {
    return c == ' ' || c == ':';
}

bool mcnp_str_newline(mcnp_str_t* s) {
    if (s->sb.error) return false;

    /* Add newline + continuation indent */
    char cont[8];
    int cont_len = snprintf(cont, sizeof(cont), "\n%*s", s->cont_indent, "");

    if (!str_builder_write(&s->sb, cont, cont_len)) return false;
    s->col = s->cont_indent + 1;

    return true;
}

bool mcnp_str_putc(mcnp_str_t* s, char c) {
    if (s->sb.error) return false;

    /* If at/past column limit and this is a break char, break instead */
    if (s->col >= s->max_col && is_break_char(c)) {
        return mcnp_str_newline(s);
    }

    /* If past column limit and not a break char, we have a problem -
     * the token was too long. Just continue (will exceed 80 cols). */

    if (!str_builder_putc(&s->sb, c)) return false;

    if (c == '\n') {
        s->col = 1;
    } else {
        s->col++;
    }

    return true;
}

/**
 * @brief Append a token (word) that must not be split
 *
 * If token won't fit on current line, break BEFORE it.
 */
bool mcnp_str_token(mcnp_str_t* s, const char* token) {
    if (s->sb.error || !token) return false;

    size_t len = strlen(token);
    if (len == 0) return true;

    /* Check if token fits on current line */
    if (s->col + (int)len > s->max_col && s->col > s->cont_indent + 1) {
        /* Won't fit - break first */
        if (!mcnp_str_newline(s)) return false;
    }

    /* Now append the token */
    if (!str_builder_write(&s->sb, token, len)) return false;
    s->col += (int)len;

    return true;
}

/**
 * @brief Append string, treating it as sequence of tokens
 *
 * Tokens are separated by spaces or colons. Colons are kept with
 * the preceding token (so "-1:" stays together, break happens after).
 */
bool mcnp_str_puts(mcnp_str_t* s, const char* str) {
    if (s->sb.error || !str) return false;

    while (*str) {
        /* Skip any leading spaces - output them directly */
        while (*str == ' ') {
            if (!mcnp_str_putc(s, ' ')) return false;
            str++;
        }

        if (!*str) break;

        const char* token_start = str;
        const char* token_end = str;

        /* Scan token: include trailing colon if present */
        while (*token_end && *token_end != ' ') {
            token_end++;
            /* If we hit a colon, include it and stop */
            if (*(token_end - 1) == ':') {
                break;
            }
        }

        /* Extract and append token */
        if (token_end > token_start) {
            size_t token_len = token_end - token_start;

            /* Check if token fits */
            if (s->col + (int)token_len > s->max_col && s->col > s->cont_indent + 1) {
                if (!mcnp_str_newline(s)) return false;
            }

            /* Append token */
            if (!str_builder_write(&s->sb, token_start, token_len)) return false;
            s->col += (int)token_len;
        }

        str = token_end;
    }

    return true;
}

bool mcnp_str_printf(mcnp_str_t* s, const char* fmt, ...) {
    if (s->sb.error) return false;

    char tmp[256];
    va_list args;
    va_start(args, fmt);
    int len = vsnprintf(tmp, sizeof(tmp), fmt, args);
    va_end(args);

    if (len < 0 || len >= (int)sizeof(tmp)) {
        s->sb.error = true;
        return false;
    }

    return mcnp_str_puts(s, tmp);
}

bool mcnp_str_int(mcnp_str_t* s, int value) {
    char tmp[16];
    int len = snprintf(tmp, sizeof(tmp), "%d", value);
    if (len <= 0) return false;

    return mcnp_str_token(s, tmp);
}

bool mcnp_str_double(mcnp_str_t* s, double value, int precision) {
    char tmp[32];
    int len = snprintf(tmp, sizeof(tmp), "%.*g", precision, value);
    if (len <= 0) return false;

    return mcnp_str_token(s, tmp);
}

bool mcnp_str_comment(mcnp_str_t* s, const char* text) {
    if (s->sb.error) return false;

    /* Start a new line if we're not at column 1 */
    if (s->col > 1) {
        if (!str_builder_putc(&s->sb, '\n')) return false;
        s->col = 1;
    }

    if (!text || !*text) {
        /* Empty comment line */
        if (!str_builder_putc(&s->sb, 'c')) return false;
        s->col = 2;
        return true;
    }

    const char* p = text;
    while (*p) {
        /* Write "c " prefix */
        if (!str_builder_write(&s->sb, "c ", 2)) return false;
        s->col = 3;

        /* Write words until line is full */
        while (*p) {
            /* Skip leading spaces */
            while (*p == ' ') p++;
            if (!*p) break;

            /* Measure next word */
            const char* word = p;
            while (*p && *p != ' ' && *p != '\n') p++;
            size_t wlen = (size_t)(p - word);

            /* Break to new comment line if word won't fit */
            if (s->col + (int)wlen > s->max_col && s->col > 3) {
                break;
            }

            /* Space before word (except at start of line) */
            if (s->col > 3) {
                if (!str_builder_putc(&s->sb, ' ')) return false;
                s->col++;
            }

            if (!str_builder_write(&s->sb, word, wlen)) return false;
            s->col += (int)wlen;

            /* Explicit newline in text forces a new comment line */
            if (*p == '\n') { p++; break; }
        }

        /* End this comment line */
        if (*p) {
            if (!str_builder_putc(&s->sb, '\n')) return false;
            s->col = 1;
        }
    }

    return true;
}

bool mcnp_str_inline_comment(mcnp_str_t* s, const char* text) {
    if (s->sb.error) return false;
    if (!text || !*text) return true;

    if (!str_builder_write(&s->sb, " $ ", 3)) return false;
    s->col += 3;

    if (!str_builder_puts(&s->sb, text)) return false;
    s->col += (int)strlen(text);

    return true;
}

size_t mcnp_str_finish(mcnp_str_t* s) {
    return str_builder_finish(&s->sb);
}

bool mcnp_str_error(const mcnp_str_t* s) {
    return str_builder_error(&s->sb);
}

const char* mcnp_str_get(mcnp_str_t* s) {
    return str_builder_get(&s->sb);
}

bool mcnp_str_write(mcnp_str_t* s, FILE* out) {
    if (s->sb.error || !out) return false;
    str_builder_finish(&s->sb);
    fprintf(out, "%s\n", s->sb.buf);
    return true;
}
