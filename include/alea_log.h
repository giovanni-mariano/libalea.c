// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_log.h
 * @brief Public logging configuration API
 */

#ifndef ALEA_LOG_PUBLIC_H
#define ALEA_LOG_PUBLIC_H

#ifdef __cplusplus
extern "C" {
#endif

/** Log level: 0=none, 1=error, 2=warn, 3=info, 4=debug, 5=trace. */
#ifndef ALEA_LOG_LEVEL_DEFINED
#define ALEA_LOG_LEVEL_DEFINED
typedef enum {
    ALEA_LOG_LEVEL_NONE  = 0,
    ALEA_LOG_LEVEL_ERROR = 1,
    ALEA_LOG_LEVEL_WARN  = 2,
    ALEA_LOG_LEVEL_INFO  = 3,
    ALEA_LOG_LEVEL_DEBUG = 4,
    ALEA_LOG_LEVEL_TRACE = 5
} alea_log_level_t;
#endif

/**
 * @brief Receive a formatted libalea log message.
 *
 * The callback can be invoked by any thread performing libalea work. The
 * callback and user_data are borrowed until replaced or unregistered; callers
 * must synchronize callback state and unregister before destroying it. Calling
 * libalea logging APIs recursively from this callback is unsupported.
 */
typedef void (*alea_log_callback_t)(alea_log_level_t level,
                                    const char* file,
                                    int line,
                                    const char* message,
                                    void* user_data);

/** Set the process-wide log level. Configure before starting worker threads. */
void alea_log_set_level(alea_log_level_t level);

/** Get the process-wide log level. */
alea_log_level_t alea_log_get_level(void);

/**
 * @brief Set or clear the process-wide log callback.
 *
 * Passing NULL unregisters the current callback and restores stderr logging.
 */
void alea_log_set_callback(alea_log_callback_t callback, void* user_data);

#ifdef __cplusplus
}
#endif

#endif /* ALEA_LOG_PUBLIC_H */
