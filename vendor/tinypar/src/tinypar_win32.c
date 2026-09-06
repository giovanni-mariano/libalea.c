// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#if defined(_WIN32) && !defined(TINYPAR_NO_THREADS)

#include "tinypar_platform.h"

#include <process.h>
#include <stdint.h>
#include <stdlib.h>

size_t tinypar_platform_hardware_threads(void) {
    DWORD count = GetActiveProcessorCount(ALL_PROCESSOR_GROUPS);
    if (count != 0) return (size_t)count;

    SYSTEM_INFO info;
    GetSystemInfo(&info);
    return info.dwNumberOfProcessors == 0 ? 1 : (size_t)info.dwNumberOfProcessors;
}

int tinypar_platform_threading_enabled(void) {
    return 1;
}

static int tinypar_platform_group_mask(WORD group, KAFFINITY* mask) {
    DWORD bytes = 0;
    if (GetLogicalProcessorInformationEx(RelationGroup, NULL, &bytes) ||
        GetLastError() != ERROR_INSUFFICIENT_BUFFER || bytes == 0)
        return 0;

    SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX* information = malloc(bytes);
    if (!information) return 0;
    int found = 0;
    if (GetLogicalProcessorInformationEx(
            RelationGroup, information, &bytes) &&
        information->Relationship == RelationGroup &&
        group < information->Group.ActiveGroupCount) {
        *mask = information->Group.GroupInfo[group].ActiveProcessorMask;
        found = *mask != 0;
    }
    free(information);
    return found;
}

void tinypar_platform_worker_enter(size_t worker_index) {
    WORD group_count = GetActiveProcessorGroupCount();
    if (group_count <= 1) return;

    size_t ordinal = worker_index;
    for (WORD group = 0; group < group_count; group++) {
        DWORD processors = GetActiveProcessorCount(group);
        if (ordinal < processors) {
            GROUP_AFFINITY affinity;
            ZeroMemory(&affinity, sizeof(affinity));
            affinity.Group = group;
            if (!tinypar_platform_group_mask(group, &affinity.Mask)) return;
            (void)SetThreadGroupAffinity(
                GetCurrentThread(), &affinity, NULL);
            return;
        }
        ordinal -= processors;
    }
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
