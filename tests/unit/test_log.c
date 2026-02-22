// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_log.c - Unit tests for logging system
 */

#include "alea_test.h"
#include "util/alea_log.h"

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

TEST_MAIN()
