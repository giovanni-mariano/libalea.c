// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#if defined(_WIN32) && !defined(TINYPAR_NO_THREADS)

#include "tinypar_platform.h"

#include <process.h>
#include <stdint.h>

size_t tinypar_platform_hardware_threads(void) {
    PROCESSOR_NUMBER processor;
    GetCurrentProcessorNumberEx(&processor);
    DWORD count = GetActiveProcessorCount(processor.Group);
    if (count != 0) return (size_t)count;

    SYSTEM_INFO info;
    GetSystemInfo(&info);
    return info.dwNumberOfProcessors == 0 ? 1 : (size_t)info.dwNumberOfProcessors;
}

int tinypar_platform_threading_enabled(void) {
    return 1;
}

int tinypar_thread_start(tinypar_thread_t* thread, tinypar_thread_entry_t entry,
                         void* argument) {
    uintptr_t handle = _beginthreadex(NULL, 0, entry, argument, 0, NULL);
    if (handle == 0) return 0;
    *thread = (HANDLE)handle;
    return 1;
}

tinypar_join_result_t tinypar_thread_join(tinypar_thread_t* thread) {
    DWORD waited = WaitForSingleObject(*thread, INFINITE);
    if (waited != WAIT_OBJECT_0) return TINYPAR_JOIN_TERMINATION_UNKNOWN;
    return CloseHandle(*thread) != 0
        ? TINYPAR_JOIN_TERMINATED
        : TINYPAR_JOIN_TERMINATED_CLEANUP_FAILED;
}

int tinypar_mutex_init(tinypar_mutex_t* mutex) {
    InitializeCriticalSection(mutex);
    return 1;
}

int tinypar_mutex_destroy(tinypar_mutex_t* mutex) {
    DeleteCriticalSection(mutex);
    return 1;
}

int tinypar_mutex_lock(tinypar_mutex_t* mutex) {
    EnterCriticalSection(mutex);
    return 1;
}

int tinypar_mutex_unlock(tinypar_mutex_t* mutex) {
    LeaveCriticalSection(mutex);
    return 1;
}

int tinypar_condition_init(tinypar_condition_t* condition) {
    InitializeConditionVariable(condition);
    return 1;
}

int tinypar_condition_destroy(tinypar_condition_t* condition) {
    (void)condition;
    return 1;
}

int tinypar_condition_wait(tinypar_condition_t* condition,
                           tinypar_mutex_t* mutex) {
    return SleepConditionVariableCS(condition, mutex, INFINITE) != 0;
}

int tinypar_condition_signal(tinypar_condition_t* condition) {
    WakeConditionVariable(condition);
    return 1;
}

int tinypar_condition_broadcast(tinypar_condition_t* condition) {
    WakeAllConditionVariable(condition);
    return 1;
}

#endif
