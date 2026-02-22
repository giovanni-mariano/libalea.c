// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_log.h
 * @brief Logging system for Alea
 *
 * Provides configurable logging with multiple levels. Can be used globally
 * or with per-system configuration.
 *
 * Usage:
 *   ALEA_LOG_ERROR("Failed to load file: %s", filename);
 *   ALEA_LOG_WARN("Using default tolerance");
 *   ALEA_LOG_INFO("Loaded %zu cells", cell_count);
 *   ALEA_LOG_DEBUG("Node %u evaluated to %.6f", node_id, result);
 *   ALEA_LOG_TRACE("Entering function %s", __func__);
 */

#ifndef ALEA_LOG_H
#define ALEA_LOG_H

#include <stdio.h>
#include <stdarg.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Log levels (in order of severity)
 */
#ifndef ALEA_LOG_LEVEL_DEFINED
#define ALEA_LOG_LEVEL_DEFINED
typedef enum {
    ALEA_LOG_LEVEL_NONE  = 0,   // No logging
    ALEA_LOG_LEVEL_ERROR = 1,   // Errors only
    ALEA_LOG_LEVEL_WARN  = 2,   // Warnings and errors
    ALEA_LOG_LEVEL_INFO  = 3,   // Informational messages
    ALEA_LOG_LEVEL_DEBUG = 4,   // Debug output
    ALEA_LOG_LEVEL_TRACE = 5    // Detailed tracing
} alea_log_level_t;
#endif

/**
 * @brief Log output target
 */
typedef enum {
    ALEA_LOG_TARGET_STDERR = 0,  // Default: stderr
    ALEA_LOG_TARGET_STDOUT = 1,  // Standard output
    ALEA_LOG_TARGET_FILE   = 2,  // Custom file
    ALEA_LOG_TARGET_CALLBACK = 3 // User callback
} alea_log_target_t;

/**
 * @brief Callback type for custom log handlers
 *
 * @param level Log level of the message
 * @param file Source file (may be NULL if not available)
 * @param line Source line (0 if not available)
 * @param message Formatted message string
 * @param user_data User-provided context
 */
typedef void (*alea_log_callback_t)(alea_log_level_t level, const char* file,
                                    int line, const char* message, void* user_data);

/**
 * @brief Log configuration structure
 */
typedef struct {
    alea_log_level_t level;           // Current log level
    alea_log_target_t target;         // Output target
    FILE* file;                      // File handle (for ALEA_LOG_TARGET_FILE)
    alea_log_callback_t callback;     // Callback (for ALEA_LOG_TARGET_CALLBACK)
    void* callback_data;             // User data for callback
    int show_timestamp;              // Include timestamp in output
    int show_level;                  // Include level name in output
    int show_location;               // Include file:line in output
} alea_log_config_t;

// ============================================================================
// GLOBAL CONFIGURATION
// ============================================================================

/**
 * @brief Set global log level
 */
void alea_log_set_level(alea_log_level_t level);

/**
 * @brief Get current global log level
 */
alea_log_level_t alea_log_get_level(void);

/**
 * @brief Set log output to a file
 */
void alea_log_set_file(FILE* file);

/**
 * @brief Set custom log callback
 */
void alea_log_set_callback(alea_log_callback_t callback, void* user_data);

/**
 * @brief Configure timestamp display
 */
void alea_log_show_timestamp(int enable);

/**
 * @brief Configure file:line display (for DEBUG and TRACE)
 */
void alea_log_show_location(int enable);

// ============================================================================
// CORE LOGGING FUNCTIONS
// ============================================================================

/**
 * @brief Core logging function (use macros instead)
 */
void alea_log_write(alea_log_level_t level, const char* file, int line,
                   const char* fmt, ...);

/**
 * @brief Check if a log level is enabled
 */
static inline int alea_log_enabled(alea_log_level_t level) {
    extern alea_log_level_t g_alea_log_level;
    return level <= g_alea_log_level;
}

// ============================================================================
// LOGGING MACROS
// ============================================================================

// Enable/disable compile-time log levels
#ifndef ALEA_LOG_MIN_LEVEL
#define ALEA_LOG_MIN_LEVEL ALEA_LOG_LEVEL_TRACE
#endif

// Conditional logging based on compile-time minimum level
#if ALEA_LOG_MIN_LEVEL >= ALEA_LOG_LEVEL_ERROR
#define ALEA_LOG_ERROR(fmt, ...) \
    alea_log_write(ALEA_LOG_LEVEL_ERROR, __FILE__, __LINE__, fmt, ##__VA_ARGS__)
#else
#define ALEA_LOG_ERROR(fmt, ...) ((void)0)
#endif

#if ALEA_LOG_MIN_LEVEL >= ALEA_LOG_LEVEL_WARN
#define ALEA_LOG_WARN(fmt, ...) \
    alea_log_write(ALEA_LOG_LEVEL_WARN, __FILE__, __LINE__, fmt, ##__VA_ARGS__)
#else
#define ALEA_LOG_WARN(fmt, ...) ((void)0)
#endif

#if ALEA_LOG_MIN_LEVEL >= ALEA_LOG_LEVEL_INFO
#define ALEA_LOG_INFO(fmt, ...) \
    alea_log_write(ALEA_LOG_LEVEL_INFO, __FILE__, __LINE__, fmt, ##__VA_ARGS__)
#else
#define ALEA_LOG_INFO(fmt, ...) ((void)0)
#endif

#if ALEA_LOG_MIN_LEVEL >= ALEA_LOG_LEVEL_DEBUG
#define ALEA_LOG_DEBUG(fmt, ...) \
    alea_log_write(ALEA_LOG_LEVEL_DEBUG, __FILE__, __LINE__, fmt, ##__VA_ARGS__)
#else
#define ALEA_LOG_DEBUG(fmt, ...) ((void)0)
#endif

#if ALEA_LOG_MIN_LEVEL >= ALEA_LOG_LEVEL_TRACE
#define ALEA_LOG_TRACE(fmt, ...) \
    alea_log_write(ALEA_LOG_LEVEL_TRACE, __FILE__, __LINE__, fmt, ##__VA_ARGS__)
#else
#define ALEA_LOG_TRACE(fmt, ...) ((void)0)
#endif

// Conditional logging that checks runtime level first (more efficient for disabled logs)
#define ALEA_LOG_IF(level, fmt, ...) \
    do { if (alea_log_enabled(level)) alea_log_write(level, __FILE__, __LINE__, fmt, ##__VA_ARGS__); } while(0)

// Short aliases
#define LOG_E ALEA_LOG_ERROR
#define LOG_W ALEA_LOG_WARN
#define LOG_I ALEA_LOG_INFO
#define LOG_D ALEA_LOG_DEBUG
#define LOG_T ALEA_LOG_TRACE

#ifdef __cplusplus
}
#endif

#endif // ALEA_LOG_H
