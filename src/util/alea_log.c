// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_log.c
 * @brief Logging system implementation
 */

#include "alea_log.h"
#include <string.h>
#include <time.h>
#include <unistd.h>

// ============================================================================
// GLOBAL STATE
// ============================================================================
// NOTE: Logging configuration is intentionally global (not thread-local).
// Log level and output settings should be configured BEFORE spawning threads.
// Reading the log level is safe from multiple threads (word-sized reads are atomic
// on common architectures), but changing configuration during parallel execution
// may cause benign races (missed/extra log messages).

// Global log level (exposed for inline check in header)
alea_log_level_t g_alea_log_level = ALEA_LOG_LEVEL_WARN;

// Full configuration (set before parallelization)
static alea_log_config_t g_config = {
    .level = ALEA_LOG_LEVEL_WARN,
    .target = ALEA_LOG_TARGET_STDERR,
    .file = NULL,
    .callback = NULL,
    .callback_data = NULL,
    .show_timestamp = 0,
    .show_level = 1,
    .show_location = 0
};

// Level names for output
static const char* level_names[] = {
    "NONE",
    "ERROR",
    "WARN",
    "INFO",
    "DEBUG",
    "TRACE"
};

// ANSI color codes for terminals
static const char* level_colors[] = {
    "",           // NONE
    "\033[31m",   // ERROR - red
    "\033[33m",   // WARN - yellow
    "\033[32m",   // INFO - green
    "\033[36m",   // DEBUG - cyan
    "\033[35m"    // TRACE - magenta
};

static const char* color_reset = "\033[0m";

// ============================================================================
// CONFIGURATION FUNCTIONS
// ============================================================================

void alea_log_set_level(alea_log_level_t level) {
    g_config.level = level;
    g_alea_log_level = level;
}

alea_log_level_t alea_log_get_level(void) {
    return g_config.level;
}

void alea_log_set_file(FILE* file) {
    if (file) {
        g_config.target = ALEA_LOG_TARGET_FILE;
        g_config.file = file;
    } else {
        g_config.target = ALEA_LOG_TARGET_STDERR;
        g_config.file = NULL;
    }
}

void alea_log_set_callback(alea_log_callback_t callback, void* user_data) {
    if (callback) {
        g_config.target = ALEA_LOG_TARGET_CALLBACK;
        g_config.callback = callback;
        g_config.callback_data = user_data;
    } else {
        g_config.target = ALEA_LOG_TARGET_STDERR;
        g_config.callback = NULL;
        g_config.callback_data = NULL;
    }
}

void alea_log_show_timestamp(int enable) {
    g_config.show_timestamp = enable;
}

void alea_log_show_location(int enable) {
    g_config.show_location = enable;
}

// ============================================================================
// CORE LOGGING
// ============================================================================

// Check if output is a terminal (for color support)
static int is_tty(FILE* f) {
#ifdef _WIN32
    (void)f;
    return 0;  // Disable colors on Windows by default
#else
    int fd = -1;
    if (f == stdout) fd = 1;
    else if (f == stderr) fd = 2;
    else if (f == stdin) fd = 0;
    else return 0;
    return isatty(fd);
#endif
}

static void alea_log_writev(alea_log_level_t level, const char* file, int line,
                            const char* fmt, va_list args) {
    // Check if this level is enabled
    if (level > g_config.level || level == ALEA_LOG_LEVEL_NONE) {
        return;
    }

    // Format the user message first
    char message[2048];
    vsnprintf(message, sizeof(message), fmt, args);

    // Handle callback target
    if (g_config.target == ALEA_LOG_TARGET_CALLBACK && g_config.callback) {
        g_config.callback(level, file, line, message, g_config.callback_data);
        return;
    }

    // Determine output file
    FILE* out;
    switch (g_config.target) {
        case ALEA_LOG_TARGET_STDOUT:
            out = stdout;
            break;
        case ALEA_LOG_TARGET_FILE:
            out = g_config.file ? g_config.file : stderr;
            break;
        default:
            out = stderr;
            break;
    }

    int use_color = is_tty(out);

    // Build output line
    char buffer[4096];
    char* ptr = buffer;
    size_t remaining = sizeof(buffer);
    int n;

    // Timestamp
    if (g_config.show_timestamp) {
        time_t now = time(NULL);
        struct tm* tm_info = localtime(&now);
        n = snprintf(ptr, remaining, "[%02d:%02d:%02d] ",
                     tm_info->tm_hour, tm_info->tm_min, tm_info->tm_sec);
        if (n > 0 && (size_t)n < remaining) {
            ptr += n;
            remaining -= n;
        }
    }

    // Level
    if (g_config.show_level) {
        if (use_color) {
            n = snprintf(ptr, remaining, "%s[%-5s]%s ",
                        level_colors[level], level_names[level], color_reset);
        } else {
            n = snprintf(ptr, remaining, "[%-5s] ", level_names[level]);
        }
        if (n > 0 && (size_t)n < remaining) {
            ptr += n;
            remaining -= n;
        }
    }

    // Location (file:line) - typically only for DEBUG/TRACE
    if (g_config.show_location && file && (level >= ALEA_LOG_LEVEL_DEBUG)) {
        // Extract just the filename, not full path
        const char* basename = file;
        const char* slash = strrchr(file, '/');
        if (slash) basename = slash + 1;
#ifdef _WIN32
        slash = strrchr(basename, '\\');
        if (slash) basename = slash + 1;
#endif
        n = snprintf(ptr, remaining, "(%s:%d) ", basename, line);
        if (n > 0 && (size_t)n < remaining) {
            ptr += n;
            remaining -= n;
        }
    }

    // Message
    n = snprintf(ptr, remaining, "%s\n", message);

    // Write to output
    fputs(buffer, out);
    fflush(out);
}

void alea_log_write(alea_log_level_t level, const char* file, int line,
                   const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    alea_log_writev(level, file, line, fmt, args);
    va_end(args);
}
