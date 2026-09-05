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

tinypar_join_result_t tinypar_thread_join(tinypar_thread_t* thread) {
    (void)thread;
    return TINYPAR_JOIN_TERMINATION_UNKNOWN;
}

int tinypar_mutex_init(tinypar_mutex_t* mutex) {
    (void)mutex;
    return 1;
}

int tinypar_mutex_destroy(tinypar_mutex_t* mutex) {
    (void)mutex;
    return 1;
}

int tinypar_mutex_lock(tinypar_mutex_t* mutex) {
    (void)mutex;
    return 1;
}

int tinypar_mutex_unlock(tinypar_mutex_t* mutex) {
    (void)mutex;
    return 1;
}

int tinypar_condition_init(tinypar_condition_t* condition) {
    (void)condition;
    return 1;
}

int tinypar_condition_destroy(tinypar_condition_t* condition) {
    (void)condition;
    return 1;
}

int tinypar_condition_wait(tinypar_condition_t* condition,
                           tinypar_mutex_t* mutex) {
    (void)condition;
    (void)mutex;
    return 0;
}

int tinypar_condition_signal(tinypar_condition_t* condition) {
    (void)condition;
    return 1;
}

int tinypar_condition_broadcast(tinypar_condition_t* condition) {
    (void)condition;
    return 1;
}
