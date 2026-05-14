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

#include <stdbool.h>
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
 * @brief Case-insensitive string comparison (bounded)
 * @param s1 First string
 * @param s2 Second string
 * @param n Maximum number of characters to compare
 * @return <0 if s1<s2, 0 if equal, >0 if s1>s2
 */
int alea_strncasecmp(const char* s1, const char* s2, size_t n);

/**
 * @brief Create and open a temporary file with a unique name
 * @param path_out Buffer to receive the temp file path (at least 256 bytes)
 * @return FILE* open for writing, or NULL on failure. Caller must fclose().
 */
FILE* alea_tmpfile(char* path_out);

/**
 * @brief Opaque directory handle for portable directory iteration
 */
typedef struct alea_dir alea_dir_t;

/**
 * @brief Open a directory for iteration
 * @param path Directory path
 * @return Handle, or NULL on failure
 */
alea_dir_t* alea_dir_open(const char* path);

/**
 * @brief Read the next entry name from an open directory
 * @param dir Handle from alea_dir_open
 * @return Entry name (valid until next call or close), or NULL when done
 */
const char* alea_dir_next(alea_dir_t* dir);

/**
 * @brief Close a directory handle
 * @param dir Handle from alea_dir_open
 */
void alea_dir_close(alea_dir_t* dir);

/**
 * @brief Approximate physical memory currently available for new allocations
 *
 * Linux: parses MemAvailable from /proc/meminfo.
 * macOS: uses host_statistics64 / HOST_VM_INFO64.
 * Windows: uses GlobalMemoryStatusEx (ullAvailPhys).
 *
 * Returns 0 if the value cannot be determined; callers should treat 0 as
 * "unknown" and skip any memory-budget enforcement rather than refusing work.
 */
size_t alea_mem_available_bytes(void);

#endif /* ALEA_COMPAT_H */
