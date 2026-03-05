// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file compat.h
 * @brief Portable replacements for POSIX functions
 *
 * Provides cross-platform implementations of common string functions
 * that may not be available on all systems.
 */

#ifndef ALEA_COMPAT_H
#define ALEA_COMPAT_H

#include <stddef.h>
#include <stdio.h>


/**
 * @brief Portable strdup - duplicate a string
 * @param s String to duplicate
 * @return Newly allocated copy (caller must free), or NULL on failure
 */
char* alea_strdup(const char* s);

/**
 * @brief Case-insensitive string comparison
 * @param s1 First string
 * @param s2 Second string
 * @return <0 if s1<s2, 0 if equal, >0 if s1>s2
 */
int alea_strcasecmp(const char* s1, const char* s2);

/**
 * @brief Create and open a temporary file with a unique name
 * @param path_out Buffer to receive the temp file path (at least 256 bytes)
 * @return FILE* open for writing, or NULL on failure. Caller must fclose().
 */
FILE* alea_tmpfile(char* path_out);


#endif /* ALEA_COMPAT_H */
