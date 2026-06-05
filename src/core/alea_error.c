// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_error.c
 * @brief Unified error handling implementation
 *
 * Single source of truth for error state. Thread-local storage
 * ensures thread safety without locks.
 */

#include "alea_types.h"
#include "util/compat.h"
#include <stdarg.h>
#include <stdio.h>
#include <signal.h>

/* Global interrupt flag (signal-safe, cross-thread visible) */
volatile sig_atomic_t g_alea_interrupted = 0;

/* Thread-local error state */
static ALEA_THREAD_LOCAL alea_error_t g_last_error = ALEA_OK;
static ALEA_THREAD_LOCAL char g_error_message[256] = {0};

/* ============================================================================
 * ERROR CODE TO STRING
 * ============================================================================ */

const char* alea_error_string(alea_error_t error) {
    switch (error) {
        case ALEA_OK:                  return "Success";
        case ALEA_ERR_NULL_ARG:        return "NULL argument";
        case ALEA_ERR_INVALID_ID:      return "Invalid ID";
        case ALEA_ERR_INVALID_ARG:     return "Invalid argument";
        case ALEA_ERR_INVALID_STATE:   return "Invalid state";
        case ALEA_ERR_OUT_OF_MEMORY:   return "Out of memory";
        case ALEA_ERR_FILE_NOT_FOUND:  return "File not found";
        case ALEA_ERR_FILE_READ:       return "File read error";
        case ALEA_ERR_FILE_WRITE:      return "File write error";
        case ALEA_ERR_PARSE_ERROR:     return "Parse error";
        case ALEA_ERR_UNSUPPORTED:     return "Unsupported feature";
        case ALEA_ERR_UNSUPPORTED_SURFACE: return "Unsupported surface type";
        case ALEA_ERR_EXPORT_FAILED:   return "Export failed";
        case ALEA_ERR_NOT_IMPLEMENTED: return "Not implemented";
        case ALEA_ERR_NOT_FOUND:       return "Not found";
        case ALEA_ERR_EMPTY:           return "Empty";
        case ALEA_ERR_INTERRUPTED:     return "Operation interrupted";
        case ALEA_ERR_OVERFLOW:        return "Buffer overflow";
        default:                      return "Unknown error";
    }
}

/* ============================================================================
 * ERROR DETAIL (thread-local detailed messages)
 * ============================================================================ */

void alea_set_error_detail(alea_error_t code, const char* fmt, ...) {
    g_last_error = code;
    if (fmt) {
        va_list args;
        va_start(args, fmt);
        vsnprintf(g_error_message, sizeof(g_error_message), fmt, args);
        va_end(args);
    } else {
        g_error_message[0] = '\0';
    }
}

const char* alea_get_error_detail(void) {
    if (g_error_message[0]) {
        return g_error_message;
    }
    return alea_error_string(g_last_error);
}

alea_error_t alea_get_last_error(void) {
    return g_last_error;
}

void alea_clear_error_detail(void) {
    g_last_error = ALEA_OK;
    g_error_message[0] = '\0';
}

/* ============================================================================
 * LEGACY ALIASES (for existing internal code)
 *
 * These will be removed once all internal code is migrated.
 * ============================================================================ */

void alea_set_error(alea_error_t code, const char* fmt, ...) {
    g_last_error = code;
    if (fmt) {
        va_list args;
        va_start(args, fmt);
        vsnprintf(g_error_message, sizeof(g_error_message), fmt, args);
        va_end(args);
    } else {
        g_error_message[0] = '\0';
    }
}

const char* alea_get_error(void) {
    return g_error_message[0] ? g_error_message : "No error";
}

