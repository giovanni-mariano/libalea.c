// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef TINYPAR_PLATFORM_H
#define TINYPAR_PLATFORM_H

#include <stddef.h>

#if defined(TINYPAR_NO_THREADS)

typedef int tinypar_thread_t;
typedef int tinypar_mutex_t;
typedef int tinypar_condition_t;
#if defined(_WIN32)
typedef unsigned(__stdcall* tinypar_thread_entry_t)(void* argument);
#else
typedef void* (*tinypar_thread_entry_t)(void* argument);
#endif

#elif defined(_WIN32)

#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0601
#endif

#include <windows.h>

typedef HANDLE tinypar_thread_t;
typedef CRITICAL_SECTION tinypar_mutex_t;
typedef CONDITION_VARIABLE tinypar_condition_t;
typedef unsigned(__stdcall* tinypar_thread_entry_t)(void* argument);

#else

#include <pthread.h>

typedef pthread_t tinypar_thread_t;
typedef pthread_mutex_t tinypar_mutex_t;
typedef pthread_cond_t tinypar_condition_t;
typedef void* (*tinypar_thread_entry_t)(void* argument);

#endif

size_t tinypar_platform_hardware_threads(void);
int tinypar_platform_threading_enabled(void);
int tinypar_thread_start(tinypar_thread_t* thread, tinypar_thread_entry_t entry,
                         void* argument);

typedef enum tinypar_join_result {
    TINYPAR_JOIN_TERMINATED = 0,
    TINYPAR_JOIN_TERMINATED_CLEANUP_FAILED,
    TINYPAR_JOIN_TERMINATION_UNKNOWN
} tinypar_join_result_t;

tinypar_join_result_t tinypar_thread_join(tinypar_thread_t* thread);

int tinypar_mutex_init(tinypar_mutex_t* mutex);
int tinypar_mutex_destroy(tinypar_mutex_t* mutex);
int tinypar_mutex_lock(tinypar_mutex_t* mutex);
int tinypar_mutex_unlock(tinypar_mutex_t* mutex);

int tinypar_condition_init(tinypar_condition_t* condition);
int tinypar_condition_destroy(tinypar_condition_t* condition);
int tinypar_condition_wait(tinypar_condition_t* condition,
                           tinypar_mutex_t* mutex);
int tinypar_condition_signal(tinypar_condition_t* condition);
int tinypar_condition_broadcast(tinypar_condition_t* condition);

#endif /* TINYPAR_PLATFORM_H */
