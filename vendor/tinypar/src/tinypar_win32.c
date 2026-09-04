// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#if defined(_WIN32) && !defined(TINYPAR_NO_THREADS)

#include "tinypar_platform.h"

#include <process.h>
#include <stdint.h>

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

int tinypar_thread_start(tinypar_thread_t* thread, tinypar_thread_entry_t entry,
                         void* argument) {
    uintptr_t handle = _beginthreadex(NULL, 0, entry, argument, 0, NULL);
    if (handle == 0) return 0;
    *thread = (HANDLE)handle;
    return 1;
}

int tinypar_thread_join(tinypar_thread_t* thread) {
    DWORD waited = WaitForSingleObject(*thread, INFINITE);
    BOOL closed = CloseHandle(*thread);
    return waited == WAIT_OBJECT_0 && closed != 0;
}

#endif
