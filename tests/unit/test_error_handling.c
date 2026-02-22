// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_error_handling.c - Error handling and robustness tests
 */

#include "alea_test.h"
#include "alea.h"
#include "core/alea_system.h"
#include "core/alea_eval.h"
#include <string.h>
#include <stdio.h>

/* ========================================================================= */
/* NULL argument handling                                                    */
/* ========================================================================= */

TEST(null_system_create) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    alea_destroy(sys);
}

TEST(null_system_destroy) {
    /* Should not crash */
    alea_destroy(NULL);
}

TEST(null_find_cell) {
    int result = alea_find_cell(NULL, 0, 0, 0);
    ASSERT(result < 0);
}

TEST(null_export) {
    int result = alea_export_mcnp(NULL, "null_test_tmp.mcnp");
    ASSERT(result != 0);
}

TEST(null_load_mcnp) {
    alea_system_t* sys = alea_load_mcnp(NULL);
    ASSERT_NULL(sys);
}

TEST(null_load_mcnp_string) {
    alea_system_t* sys = alea_load_mcnp_string(NULL, 0);
    ASSERT_NULL(sys);
}

/* ========================================================================= */
/* Invalid IDs                                                               */
/* ========================================================================= */

TEST(invalid_node_id) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Evaluating invalid node should not crash */
    bool inside = alea_point_inside(sys, ALEA_NODE_ID_INVALID, 0, 0, 0);
    ASSERT_FALSE(inside);

    alea_destroy(sys);
}

TEST(invalid_cell_index) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    alea_cell_info_t info;
    int rc = alea_cell_get_info(sys, 99999, &info);
    ASSERT(rc != 0);

    alea_destroy(sys);
}

TEST(invalid_surface_id) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int idx = alea_surface_find(sys, 99999);
    ASSERT(idx < 0);

    alea_destroy(sys);
}

/* ========================================================================= */
/* Parse error handling                                                      */
/* ========================================================================= */

TEST(parse_empty_string) {
    alea_system_t* sys = alea_load_mcnp_string("", 0);
    ASSERT_NULL(sys);
}

TEST(parse_title_only) {
    const char* input = "Title only\n";
    alea_system_t* sys = alea_load_mcnp_string(input, strlen(input));
    /* Should parse but have no cells */
    if (sys) {
        ASSERT_EQ(alea_cell_count(sys), 0);
        alea_destroy(sys);
    }
}

TEST(parse_bad_surface_type) {
    const char* input =
        "Test bad surface\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 XYZZY 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    alea_system_t* sys = alea_load_mcnp_string(input, strlen(input));
    /* Should either fail to parse or handle gracefully */
    if (sys) alea_destroy(sys);
    /* Pass if no crash */
}

/* ========================================================================= */
/* Resource limits                                                           */
/* ========================================================================= */

TEST(empty_system_export) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Export empty system */
    int rc = alea_export_mcnp(sys, "test_empty_tmp.mcnp");
    /* Should succeed with empty output or return error */
    remove("test_empty_tmp.mcnp");

    alea_destroy(sys);
    (void)rc; /* May succeed or fail, just don't crash */
}

TEST(empty_system_queries) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Query on empty system */
    int cell = alea_find_cell(sys, 0, 0, 0);
    ASSERT(cell < 0);

    int mat = alea_material_at(sys, 0, 0, 0);
    ASSERT_EQ(mat, 0);

    alea_destroy(sys);
}

/* ========================================================================= */
/* Error message functions                                                   */
/* ========================================================================= */

TEST(error_string) {
    const char* msg = alea_error_string(ALEA_OK);
    ASSERT_NOT_NULL(msg);

    msg = alea_error_string(ALEA_ERR_NULL_ARG);
    ASSERT_NOT_NULL(msg);

    msg = alea_error_string(ALEA_ERR_PARSE_ERROR);
    ASSERT_NOT_NULL(msg);
}

TEST(error_clear) {
    alea_set_error_detail(ALEA_ERR_NULL_ARG, "test error");
    ASSERT_EQ(alea_get_last_error(), ALEA_ERR_NULL_ARG);

    alea_clear_error_detail();
    ASSERT_EQ(alea_get_last_error(), ALEA_OK);
}

/* ========================================================================= */
/* System lifecycle                                                          */
/* ========================================================================= */

TEST(system_reset) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Add some geometry */
    alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    ASSERT_EQ(alea_surface_count(sys), 1);

    /* Reset should clear everything */
    alea_reset(sys);
    ASSERT_EQ(alea_cell_count(sys), 0);

    alea_destroy(sys);
}

TEST(system_clone) {
    const char* input =
        "Test clone\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    alea_system_t* sys = alea_load_mcnp_string(input, strlen(input));
    ASSERT_NOT_NULL(sys);
    alea_build_universe_index(sys);

    alea_system_t* clone = alea_clone(sys);
    ASSERT_NOT_NULL(clone);
    alea_build_universe_index(clone);

    /* Both should give same answers */
    ASSERT_EQ(alea_material_at(sys, 0, 0, 0),
              alea_material_at(clone, 0, 0, 0));
    ASSERT_EQ(alea_material_at(sys, 10, 0, 0),
              alea_material_at(clone, 10, 0, 0));

    alea_destroy(clone);
    alea_destroy(sys);
}

TEST_MAIN()
