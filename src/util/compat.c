// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file compat.c
 * @brief Portable replacements for POSIX functions
 */

#ifndef _WIN32
#define _POSIX_C_SOURCE 200809L
#endif

#include "compat.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

char* alea_strdup(const char* s) {
    if (!s) return NULL;
    size_t len = strlen(s) + 1;
    char* copy = malloc(len);
    if (copy) {
        memcpy(copy, s, len);
    }
    return copy;
}

FILE* alea_tmpfile(char* path_out) {
#ifdef _WIN32
    char tmp_dir[MAX_PATH];
    if (GetTempPathA(MAX_PATH, tmp_dir) == 0) return NULL;
    if (GetTempFileNameA(tmp_dir, "alea", 0, path_out) == 0) return NULL;
    return fopen(path_out, "w");
#else
    strcpy(path_out, "/tmp/alea_XXXXXX");
    int fd = mkstemp(path_out);
    if (fd < 0) return NULL;
    FILE* f = fdopen(fd, "w");
    if (!f) { close(fd); remove(path_out); return NULL; }
    return f;
#endif
}

int alea_strcasecmp(const char* s1, const char* s2) {
    if (!s1 || !s2) {
        if (s1 == s2) return 0;
        return s1 ? 1 : -1;
    }
    while (*s1 && *s2) {
        char c1 = (*s1 >= 'A' && *s1 <= 'Z') ? *s1 + 32 : *s1;
        char c2 = (*s2 >= 'A' && *s2 <= 'Z') ? *s2 + 32 : *s2;
        if (c1 != c2) return c1 - c2;
        s1++; s2++;
    }
    return *s1 - *s2;
}

int alea_strncasecmp(const char* s1, const char* s2, size_t n) {
    if (n == 0) return 0;
    if (!s1 || !s2) {
        if (s1 == s2) return 0;
        return s1 ? 1 : -1;
    }
    for (; n > 0; n--, s1++, s2++) {
        char c1 = (*s1 >= 'A' && *s1 <= 'Z') ? *s1 + 32 : *s1;
        char c2 = (*s2 >= 'A' && *s2 <= 'Z') ? *s2 + 32 : *s2;
        if (c1 != c2) return c1 - c2;
        if (c1 == '\0') return 0;
    }
    return 0;
}
