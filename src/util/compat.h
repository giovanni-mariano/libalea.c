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
 * @brief Portable thread-local storage qualifier
 *
 * Maps to MSVC's `__declspec(thread)` and to C11 `_Thread_local` everywhere
 * else. Both require *constant* initializers, which all current uses satisfy
 * (= 0 / = {0} / NULL).
 */
#if defined(_MSC_VER) && !defined(__clang__)
  #define ALEA_THREAD_LOCAL __declspec(thread)
#else
  #define ALEA_THREAD_LOCAL _Thread_local
#endif

/**
 * @brief Monotonic wall-clock in seconds (for elapsed-time measurement)
 *
 * Windows: QueryPerformanceCounter. POSIX: gettimeofday. The absolute value is
 * unspecified; only differences between two readings are meaningful.
 */
double alea_monotonic_seconds(void);

/**
 * @brief Portable, statically-initializable blocking mutex
 *
 * Maps to a Windows SRWLOCK (ABI: a single pointer; SRWLOCK_INIT == {0}) or to
 * a POSIX pthread_mutex_t. Designed for file-scope/static locks:
 *     static alea_mutex_t m = ALEA_MUTEX_INIT;
 *     alea_mutex_lock(&m); ... alea_mutex_unlock(&m);
 */
#if defined(_WIN32)
  typedef struct { void* _ptr; } alea_mutex_t;  /* layout-compatible with SRWLOCK */
  #define ALEA_MUTEX_INIT { 0 }
  void alea_mutex_lock(alea_mutex_t* m);
  void alea_mutex_unlock(alea_mutex_t* m);
#else
  #include <pthread.h>
  typedef pthread_mutex_t alea_mutex_t;
  #define ALEA_MUTEX_INIT PTHREAD_MUTEX_INITIALIZER
  static inline void alea_mutex_lock(alea_mutex_t* m) { pthread_mutex_lock(m); }
  static inline void alea_mutex_unlock(alea_mutex_t* m) { pthread_mutex_unlock(m); }
#endif

/**
 * @brief Report whether a stdio stream refers to a terminal (honest check)
 *
 * POSIX: isatty(). Windows: _isatty(_fileno(f)). Use this for interactivity
 * decisions (e.g. whether to start a REPL). For colorization use
 * alea_file_use_color(), which additionally applies the platform color policy.
 *
 * @param f Stream to test
 * @return Non-zero if f is an interactive terminal, 0 otherwise
 */
int alea_file_is_tty(FILE* f);

/**
 * @brief Report whether colorized output should be used on a stream
 *
 * Equivalent to alea_file_is_tty() on POSIX, but returns 0 on Windows, where
 * ANSI colors are disabled by default. Kept separate from alea_file_is_tty()
 * so interactivity and color policy do not get conflated.
 *
 * @param f Stream to test (typically stdout/stderr)
 * @return Non-zero if color escapes should be emitted, 0 otherwise
 */
int alea_file_use_color(FILE* f);

/**
 * @brief Read-only memory-mapped view of a whole file
 *
 * POSIX: mmap(MAP_PRIVATE, PROT_READ). Windows: CreateFileMapping +
 * MapViewOfFile. HANDLEs are kept as void* so <windows.h> need not leak here.
 */
typedef struct {
    const char* data;  /* mapped file contents (read-only) */
    size_t size;       /* file size in bytes */
#ifdef _WIN32
    void* _file;       /* HANDLE for the open file */
    void* _mapping;    /* HANDLE for the file mapping */
#else
    int _fd;           /* open file descriptor */
#endif
} alea_mapped_file_t;

/**
 * @brief Memory-map a file read-only
 * @param path File to map
 * @param out  Handle populated on success (data/size valid until unmap)
 * @return 1 on success, 0 on failure
 */
int alea_file_map(const char* path, alea_mapped_file_t* out);

/**
 * @brief Release a mapping created by alea_file_map
 * @param mf Handle from alea_file_map
 */
void alea_file_unmap(alea_mapped_file_t* mf);

/**
 * @brief Return the file-name component of a path (no directory part)
 *
 * POSIX: strips everything up to the last '/'. Windows: strips up to the last
 * '/' or '\'. The result points into the input string (no allocation).
 *
 * @param path Path string (may be NULL, which is returned unchanged)
 * @return Pointer to the basename within @p path
 */
const char* alea_path_basename(const char* path);

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
 * @brief Portable unsigned long long parser.
 *
 * POSIX/C99: strtoull(). Older MSVC CRTs: _strtoui64().
 *
 * @param nptr Input string
 * @param endptr Receives first unparsed character when non-NULL
 * @param base Numeric base, as for strtoull()
 * @return Parsed value
 */
unsigned long long alea_strtoull(const char* nptr, char** endptr, int base);

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
