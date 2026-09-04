// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "tinypar_platform.h"

size_t tinypar_platform_hardware_threads(void) {
    return 1;
}

int tinypar_platform_threading_enabled(void) {
    return 0;
}

int tinypar_thread_start(tinypar_thread_t* thread, tinypar_thread_entry_t entry,
                         void* argument) {
    (void)thread;
    (void)entry;
    (void)argument;
    return 0;
}

int tinypar_thread_join(tinypar_thread_t* thread) {
    (void)thread;
    return 0;
}
