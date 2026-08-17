// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_log.c - Unit tests for logging system
 */

#include "alea_test.h"
#include "alea_log.h"
#include "util/alea_log.h"

static int callback_count;
static alea_log_level_t callback_level;
static void* callback_user_data;

static void test_log_callback(alea_log_level_t level, const char* file,
                              int line, const char* message, void* user_data) {
    (void)file;
    (void)line;
    (void)message;
    callback_count++;
    callback_level = level;
    callback_user_data = user_data;
}

TEST(log_default_level) {
    /* Default level is WARN, so INFO/DEBUG should be filtered */
    alea_log_set_level(ALEA_LOG_LEVEL_WARN);

    /* These shouldn't crash - just verify the macros work */
    ALEA_LOG_ERROR("Test error message");
    ALEA_LOG_WARN("Test warning message");
    ALEA_LOG_INFO("This should be filtered");
    ALEA_LOG_DEBUG("This should be filtered");

    /* If we got here without crashing, pass */
    ASSERT_TRUE(1);
}

TEST(log_level_change) {
    /* Test changing log levels */
    alea_log_set_level(ALEA_LOG_LEVEL_DEBUG);
    ALEA_LOG_DEBUG("Debug should now work");

    alea_log_set_level(ALEA_LOG_LEVEL_ERROR);
    ALEA_LOG_WARN("Warning should be filtered");

    /* Reset to default */
    alea_log_set_level(ALEA_LOG_LEVEL_WARN);
    ASSERT_TRUE(1);
}

TEST(log_formatted_output) {
    alea_log_set_level(ALEA_LOG_LEVEL_INFO);

    int value = 42;
    double result = 3.14159;
    ALEA_LOG_INFO("Formatted: value=%d, result=%.4f", value, result);

    /* Reset */
    alea_log_set_level(ALEA_LOG_LEVEL_WARN);
    ASSERT_TRUE(1);
}

TEST(log_location_and_timestamp) {
    alea_log_set_level(ALEA_LOG_LEVEL_DEBUG);
    alea_log_show_location(1);
    alea_log_show_timestamp(1);

    ALEA_LOG_DEBUG("Debug with location and timestamp");

    /* Reset */
    alea_log_show_location(0);
    alea_log_show_timestamp(0);
    alea_log_set_level(ALEA_LOG_LEVEL_WARN);
    ASSERT_TRUE(1);
}

TEST(log_short_aliases) {
    alea_log_set_level(ALEA_LOG_LEVEL_INFO);

    LOG_I("Using LOG_I alias");
    LOG_W("Using LOG_W alias");
    LOG_E("Using LOG_E alias");

    alea_log_set_level(ALEA_LOG_LEVEL_WARN);
    ASSERT_TRUE(1);
}

TEST(log_public_callback) {
    int user_data = 42;
    callback_count = 0;
    callback_level = ALEA_LOG_LEVEL_NONE;
    callback_user_data = NULL;

    alea_log_set_level(ALEA_LOG_LEVEL_INFO);
    alea_log_set_callback(test_log_callback, &user_data);
    ALEA_LOG_INFO("Public callback test");

    ASSERT_EQ(callback_count, 1);
    ASSERT_EQ(callback_level, ALEA_LOG_LEVEL_INFO);
    ASSERT_EQ(callback_user_data, &user_data);

    alea_log_set_callback(NULL, NULL);
    ALEA_LOG_INFO("Callback has been cleared");
    ASSERT_EQ(callback_count, 1);
    alea_log_set_level(ALEA_LOG_LEVEL_WARN);
}

TEST_MAIN()
