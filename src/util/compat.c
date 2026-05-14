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
#include <dirent.h>
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
