// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#if defined(__linux__) && !defined(_GNU_SOURCE)
#define _GNU_SOURCE
#endif

#if !defined(_WIN32) && !defined(TINYPAR_NO_THREADS)

#include "tinypar_platform.h"

#if defined(__linux__)
#include <sched.h>
#endif
#include <unistd.h>

size_t tinypar_platform_hardware_threads(void) {
#if defined(__linux__)
    cpu_set_t affinity;
    if (sched_getaffinity(0, sizeof(affinity), &affinity) == 0) {
        int available = CPU_COUNT(&affinity);
        if (available > 0) return (size_t)available;
    }
#endif
    long count = sysconf(_SC_NPROCESSORS_ONLN);
    return count > 0 ? (size_t)count : 1;
}

int tinypar_platform_threading_enabled(void) {
    return 1;
}

int tinypar_thread_start(tinypar_thread_t* thread, tinypar_thread_entry_t entry,
                         void* argument) {
    return pthread_create(thread, NULL, entry, argument) == 0;
}

int tinypar_thread_join(tinypar_thread_t* thread) {
    return pthread_join(*thread, NULL) == 0;
}

#endif
