// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_error_handling.c - Error handling and robustness tests
 */

#include "alea_test.h"
#include "alea.h"
#include "alea_mcnp.h"
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
    int result = mcnp_export(NULL, "null_test_tmp.mcnp");
    ASSERT(result != 0);
}

TEST(null_load_mcnp) {
    mcnp_model_t* model = mcnp_load(NULL);
    ASSERT_NULL(model);
}

TEST(null_load_mcnp_string) {
    mcnp_model_t* model = mcnp_load_string(NULL, 0);
    ASSERT_NULL(model);
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
    mcnp_model_t* model = mcnp_load_string("", 0);
    ASSERT_NULL(model);
}

TEST(parse_title_only) {
    const char* input = "Title only\n";
    mcnp_model_t* model = mcnp_load_string(input, strlen(input));
    /* Should parse but have no cells */
    if (model) {
        ASSERT_EQ(alea_cell_count(model->sys), 0);
        mcnp_model_destroy(model);
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
    mcnp_model_t* model = mcnp_load_string(input, strlen(input));
    /* Should either fail to parse or handle gracefully */
    if (model) mcnp_model_destroy(model);
    /* Pass if no crash */
}

/* ========================================================================= */
/* Resource limits                                                           */
/* ========================================================================= */

TEST(empty_system_export) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Export empty system */
    FILE* f = tmpfile();
    int rc = f ? mcnp_export_system_stream(sys, f) : -1;
    if (f) fclose(f);
    /* Should succeed with empty output or return error */

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
    mcnp_model_t* model = mcnp_load_string(input, strlen(input));
    ASSERT_NOT_NULL(model);
    alea_build_universe_index(model->sys);

    alea_system_t* clone = alea_clone(model->sys);
    ASSERT_NOT_NULL(clone);
    alea_build_universe_index(clone);

    /* Both should give same answers */
    ASSERT_EQ(alea_material_at(model->sys, 0, 0, 0),
              alea_material_at(clone, 0, 0, 0));
    ASSERT_EQ(alea_material_at(model->sys, 10, 0, 0),
              alea_material_at(clone, 10, 0, 0));

    alea_destroy(clone);
    mcnp_model_destroy(model);
}

TEST_MAIN()
