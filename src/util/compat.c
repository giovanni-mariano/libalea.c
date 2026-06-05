// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file compat.c
 * @brief Portable replacements for POSIX functions
 */

#if defined(__APPLE__)
/* _DARWIN_C_SOURCE keeps BSD types (u_int, u_char, u_short) visible alongside
 * the POSIX surface area — sys/sysctl.h needs them on macOS. */
#define _DARWIN_C_SOURCE
#elif !defined(_WIN32)
#define _POSIX_C_SOURCE 200809L
#endif

#include "compat.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifdef _WIN32
#include <windows.h>
#include <io.h>
#else
#include <unistd.h>
#include <dirent.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#endif

int alea_file_is_tty(FILE* f) {
#ifdef _WIN32
    return _isatty(_fileno(f));
#else
    int fd = -1;
    if (f == stdout) fd = 1;
    else if (f == stderr) fd = 2;
    else if (f == stdin) fd = 0;
    else return 0;
    return isatty(fd);
#endif
}

int alea_file_use_color(FILE* f) {
#ifdef _WIN32
    /* ANSI colors are disabled by default on Windows consoles. */
    (void)f;
    return 0;
#else
    return alea_file_is_tty(f);
#endif
}

int alea_file_map(const char* path, alea_mapped_file_t* mf) {
#ifdef _WIN32
    mf->_file = CreateFileA(path, GENERIC_READ, FILE_SHARE_READ, NULL,
                            OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if (mf->_file == INVALID_HANDLE_VALUE) return 0;

    DWORD sz = GetFileSize((HANDLE)mf->_file, NULL);
    if (sz == INVALID_FILE_SIZE) { CloseHandle((HANDLE)mf->_file); return 0; }
    mf->size = sz;

    mf->_mapping = CreateFileMapping((HANDLE)mf->_file, NULL, PAGE_READONLY, 0, 0, NULL);
    if (mf->_mapping == NULL) { CloseHandle((HANDLE)mf->_file); return 0; }

    mf->data = (const char*)MapViewOfFile((HANDLE)mf->_mapping, FILE_MAP_READ, 0, 0, 0);
    if (mf->data == NULL) {
        CloseHandle((HANDLE)mf->_mapping);
        CloseHandle((HANDLE)mf->_file);
        return 0;
    }
#else
    mf->_fd = open(path, O_RDONLY);
    if (mf->_fd == -1) return 0;

    struct stat s;
    if (fstat(mf->_fd, &s) == -1) { close(mf->_fd); return 0; }
    mf->size = (size_t)s.st_size;

    mf->data = (const char*)mmap(NULL, mf->size, PROT_READ, MAP_PRIVATE, mf->_fd, 0);
    if (mf->data == MAP_FAILED) { close(mf->_fd); return 0; }
#endif
    return 1;
}

void alea_file_unmap(alea_mapped_file_t* mf) {
#ifdef _WIN32
    UnmapViewOfFile(mf->data);
    CloseHandle((HANDLE)mf->_mapping);
    CloseHandle((HANDLE)mf->_file);
#else
    munmap((void*)mf->data, mf->size);
    close(mf->_fd);
#endif
}

double alea_monotonic_seconds(void) {
#ifdef _WIN32
    LARGE_INTEGER freq, counter;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&counter);
    return (double)counter.QuadPart / (double)freq.QuadPart;
#else
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;
#endif
}

#ifdef _WIN32
/* alea_mutex_t is layout-compatible with SRWLOCK (one pointer); cast through. */
void alea_mutex_lock(alea_mutex_t* m)   { AcquireSRWLockExclusive((PSRWLOCK)m); }
void alea_mutex_unlock(alea_mutex_t* m) { ReleaseSRWLockExclusive((PSRWLOCK)m); }
#endif

const char* alea_path_basename(const char* path) {
    if (!path) return path;
    const char* base = path;
    const char* slash = strrchr(path, '/');
    if (slash) base = slash + 1;
#ifdef _WIN32
    slash = strrchr(base, '\\');
    if (slash) base = slash + 1;
#endif
    return base;
}

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

unsigned long long alea_strtoull(const char* nptr, char** endptr, int base) {
#if defined(_MSC_VER) && !defined(__clang__)
    return _strtoui64(nptr, endptr, base);
#else
    return strtoull(nptr, endptr, base);
#endif
}

/* --- Portable directory iteration --- */

#ifdef _WIN32

struct alea_dir {
    HANDLE handle;
    WIN32_FIND_DATAA data;
    bool first;
};

alea_dir_t* alea_dir_open(const char* path) {
    if (!path) return NULL;

    /* Build search pattern: path\* */
    size_t len = strlen(path);
    char* pattern = malloc(len + 3);
    if (!pattern) return NULL;
    memcpy(pattern, path, len);
    /* Normalize trailing separator */
    if (len > 0 && path[len - 1] != '\\' && path[len - 1] != '/')
        pattern[len++] = '\\';
    pattern[len] = '*';
    pattern[len + 1] = '\0';

    alea_dir_t* dir = malloc(sizeof(*dir));
    if (!dir) { free(pattern); return NULL; }

    dir->handle = FindFirstFileA(pattern, &dir->data);
    free(pattern);
    if (dir->handle == INVALID_HANDLE_VALUE) { free(dir); return NULL; }

    dir->first = true;
    return dir;
}

const char* alea_dir_next(alea_dir_t* dir) {
    if (!dir) return NULL;
    if (dir->first) {
        dir->first = false;
        return dir->data.cFileName;
    }
    if (FindNextFileA(dir->handle, &dir->data))
        return dir->data.cFileName;
    return NULL;
}

void alea_dir_close(alea_dir_t* dir) {
    if (!dir) return;
    FindClose(dir->handle);
    free(dir);
}

#else /* POSIX */

struct alea_dir {
    DIR* dp;
};

alea_dir_t* alea_dir_open(const char* path) {
    if (!path) return NULL;
    DIR* dp = opendir(path);
    if (!dp) return NULL;
    alea_dir_t* dir = malloc(sizeof(*dir));
    if (!dir) { closedir(dp); return NULL; }
    dir->dp = dp;
    return dir;
}

const char* alea_dir_next(alea_dir_t* dir) {
    if (!dir) return NULL;
    struct dirent* de = readdir(dir->dp);
    return de ? de->d_name : NULL;
}

void alea_dir_close(alea_dir_t* dir) {
    if (!dir) return;
    closedir(dir->dp);
    free(dir);
}

#endif /* _WIN32 */

/* ============================================================================
 * Available physical memory
 * ============================================================================ */

#if defined(_WIN32)

size_t alea_mem_available_bytes(void) {
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    if (!GlobalMemoryStatusEx(&status)) return 0;
    return (size_t)status.ullAvailPhys;
}

#elif defined(__APPLE__)

#include <mach/mach.h>
#include <mach/vm_statistics.h>
#include <sys/sysctl.h>

size_t alea_mem_available_bytes(void) {
    mach_port_t host = mach_host_self();
    vm_size_t page_size = 0;
    if (host_page_size(host, &page_size) != KERN_SUCCESS || page_size == 0) {
        return 0;
    }
    vm_statistics64_data_t stats;
    mach_msg_type_number_t count = HOST_VM_INFO64_COUNT;
    if (host_statistics64(host, HOST_VM_INFO64,
                          (host_info64_t)&stats, &count) != KERN_SUCCESS) {
        return 0;
    }
    /* Pages that can be reclaimed without paging from disk. Mirrors what
     * Activity Monitor labels "Memory Free" plus "Inactive". */
    uint64_t pages = (uint64_t)stats.free_count +
                     (uint64_t)stats.inactive_count +
                     (uint64_t)stats.purgeable_count;
    return (size_t)(pages * (uint64_t)page_size);
}

#elif defined(__linux__)

size_t alea_mem_available_bytes(void) {
    FILE* f = fopen("/proc/meminfo", "r");
    if (!f) return 0;
    char line[256];
    size_t avail = 0;
    while (fgets(line, sizeof(line), f)) {
        unsigned long kb = 0;
        if (sscanf(line, "MemAvailable: %lu kB", &kb) == 1) {
            avail = (size_t)kb * 1024;
            break;
        }
    }
    fclose(f);
    return avail;
}

#else

size_t alea_mem_available_bytes(void) {
    return 0; /* Unknown platform: callers skip enforcement. */
}

#endif
