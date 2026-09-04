// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef TINYPAR_PLATFORM_H
#define TINYPAR_PLATFORM_H

#include <stddef.h>

#if defined(TINYPAR_NO_THREADS)

typedef int tinypar_thread_t;
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
typedef unsigned(__stdcall* tinypar_thread_entry_t)(void* argument);

#else

#include <pthread.h>

typedef pthread_t tinypar_thread_t;
typedef void* (*tinypar_thread_entry_t)(void* argument);

#endif

size_t tinypar_platform_hardware_threads(void);
int tinypar_platform_threading_enabled(void);
int tinypar_thread_start(tinypar_thread_t* thread, tinypar_thread_entry_t entry,
                         void* argument);
int tinypar_thread_join(tinypar_thread_t* thread);

#endif /* TINYPAR_PLATFORM_H */
