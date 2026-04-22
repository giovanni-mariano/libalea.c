// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_mcnp_parser_edge.c - MCNP parser edge case tests
 *
 * Tests parser handling of continuation lines, comments, whitespace,
 * case sensitivity, and other edge cases.
 */

#include "alea_test.h"
#include "alea.h"
#include "alea_mcnp.h"
#include "core/alea_system.h"
#include <string.h>
#include <stdio.h>

static mcnp_model_t* parse_mcnp(const char* input) {
    mcnp_model_t* model = mcnp_load_string(input, strlen(input));
    if (model) alea_build_universe_index(model->sys);
    return model;
}

/* ========================================================================= */
/* Continuation line tests                                                   */
/* ========================================================================= */

TEST(continuation_5_blanks) {
    /* 5+ leading spaces = continuation of previous card */
    const char* input =
        "Test continuation\n"
        "1 1 -1.0\n"
        "     -1\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_cell_count(sys), 2);
    ASSERT_EQ(alea_material_at(sys, 0, 0, 0), 1);
    mcnp_model_destroy(model);
}

TEST(continuation_ampersand) {
    /* & at end of line = continuation */
    const char* input =
        "Test ampersand\n"
        "1 1 -1.0 -1 &\n"
        "     IMP:N=1\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    ASSERT_EQ(alea_cell_count(model->sys), 2);
    mcnp_model_destroy(model);
}

/* ========================================================================= */
/* Comment tests                                                             */
/* ========================================================================= */

TEST(comment_c_column1) {
    /* C in column 1 = full-line comment */
    const char* input =
        "Test comments\n"
        "c This is a comment\n"
        "C This is also a comment\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "c Surface comment\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    ASSERT_EQ(alea_cell_count(model->sys), 2);
    mcnp_model_destroy(model);
}

TEST(comment_dollar) {
    /* $ mid-line comment */
    const char* input =
        "Test dollar comment\n"
        "1 1 -1.0 -1 $ This is mat 1\n"
        "2 0 1 $ void region\n"
        "\n"
        "1 SO 5.0 $ sphere of radius 5\n"
        "\n"
        "M1 92235.80c 1.0 $ uranium\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_cell_count(sys), 2);
    ASSERT_EQ(alea_material_at(sys, 0, 0, 0), 1);
    mcnp_model_destroy(model);
}

/* ========================================================================= */
/* Whitespace and formatting tests                                           */
/* ========================================================================= */

TEST(blank_line_separator) {
    /* Blank lines separate cell/surface/data sections */
    const char* input =
        "Test blank lines\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    ASSERT_EQ(alea_cell_count(model->sys), 2);
    mcnp_model_destroy(model);
}

TEST(trailing_whitespace) {
    /* Extra whitespace after values should be ignored */
    const char* input =
        "Test trailing ws   \n"
        "1 1 -1.0 -1   \n"
        "2 0 1   \n"
        "\n"
        "1 SO 5.0   \n"
        "\n"
        "M1 92235.80c 1.0   \n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    ASSERT_EQ(alea_cell_count(model->sys), 2);
    mcnp_model_destroy(model);
}

/* ========================================================================= */
/* Case sensitivity tests                                                    */
/* ========================================================================= */

TEST(case_insensitive) {
    /* MCNP is case insensitive for keywords */
    const char* input =
        "Test case insensitive\n"
        "1 1 -1.0 -1 imp:n=1\n"
        "2 0 1\n"
        "\n"
        "1 so 5.0\n"
        "\n"
        "m1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    ASSERT_EQ(alea_cell_count(model->sys), 2);
    mcnp_model_destroy(model);
}

/* ========================================================================= */
/* Large ID tests                                                            */
/* ========================================================================= */

TEST(large_cell_id) {
    const char* input =
        "Test large IDs\n"
        "99999 1 -1.0 -1\n"
        "99998 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_cell_count(sys), 2);
    alea_cell_info_t info;
    ASSERT_EQ(alea_cell_find_info(sys, 99999, &info), 0);
    ASSERT_EQ(info.material_id, 1);
    mcnp_model_destroy(model);
}

TEST(large_surface_id) {
    const char* input =
        "Test large surface ID\n"
        "1 1 -1.0 -99999\n"
        "2 0 99999\n"
        "\n"
        "99999 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    ASSERT_EQ(alea_material_at(model->sys, 0, 0, 0), 1);
    mcnp_model_destroy(model);
}

/* ========================================================================= */
/* Surface sense tests                                                       */
/* ========================================================================= */

TEST(negative_surface_sense) {
    /* -1 (negative sense = inside) vs 1 (positive sense = outside) */
    const char* input =
        "Test sense\n"
        "1 1 -1.0 -1\n"
        "2 2 -1.0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n"
        "M2 26056.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    /* Inside sphere → mat 1 */
    ASSERT_EQ(alea_material_at(sys, 0, 0, 0), 1);
    /* Outside sphere → mat 2 */
    ASSERT_EQ(alea_material_at(sys, 10, 0, 0), 2);
    mcnp_model_destroy(model);
}

/* ========================================================================= */
/* Material tests                                                            */
/* ========================================================================= */

TEST(multiple_materials) {
    const char* input =
        "Test multi-mat\n"
        "1 1 -1.0 -1\n"
        "2 2 -8.0 1 -2\n"
        "3 0 2\n"
        "\n"
        "1 SO 3.0\n"
        "2 SO 10.0\n"
        "\n"
        "M1 92235.80c 0.05 92238.80c 0.95\n"
        "M2 26056.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_material_at(sys, 0, 0, 0), 1);
    ASSERT_EQ(alea_material_at(sys, 5, 0, 0), 2);
    ASSERT_EQ(alea_material_at(sys, 15, 0, 0), 0);
    mcnp_model_destroy(model);
}

TEST(material_nuclides) {
    /* Material card with multiple nuclides */
    const char* input =
        "Test nuclides\n"
        "1 1 -10.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 0.05 92238.80c 0.95\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    ASSERT_EQ(alea_material_at(model->sys, 0, 0, 0), 1);
    mcnp_model_destroy(model);
}

/* ========================================================================= */
/* Transform card tests                                                      */
/* ========================================================================= */

TEST(transform_card_cosines) {
    /* TR1 with direction cosines (12 values) */
    const char* input =
        "Test TR cosines\n"
        "1 1 -1.0 -1 TRCL=1\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "TR1 10 0 0 1 0 0 0 1 0 0 0 1\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    /* Pure translation to (10,0,0) */
    ASSERT_EQ(alea_material_at(sys, 10, 0, 0), 1);
    ASSERT_EQ(alea_material_at(sys, 0, 0, 0), 0);
    mcnp_model_destroy(model);
}

TEST(transform_card_angles) {
    /* *TR1 with angles in degrees */
    const char* input =
        "Test *TR angles\n"
        "1 1 -1.0 -1 TRCL=1\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "*TR1 10 0 0 0 90 90 90 0 90 90 90 0\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    /* Translation to (10,0,0) with identity rotation */
    ASSERT_EQ(alea_material_at(model->sys, 10, 0, 0), 1);
    mcnp_model_destroy(model);
}

/* ========================================================================= */
/* Empty parameter tests                                                     */
/* ========================================================================= */

TEST(empty_params) {
    /* Cell with no parameters after geometry */
    const char* input =
        "Test no params\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    ASSERT_EQ(alea_cell_count(model->sys), 2);
    mcnp_model_destroy(model);
}

TEST(surface_with_transform) {
    /* Surface card with TRn reference: "1 1 PZ 5.0" */
    const char* input =
        "Test surface with TR\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 1 PZ 5.0\n"
        "\n"
        "TR1 10 0 0\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    /* PZ 5.0 translated by TR1 → plane at z=5 shifted by (10,0,0) */
    /* Actually the plane normal is along z, so the translation in x/y
       doesn't change the z-intercept. The plane is still at z=5+0=5. */
    ASSERT_EQ(alea_material_at(model->sys, 10, 0, 3), 1);
    mcnp_model_destroy(model);
}

/* ========================================================================= */
/* Parse error handling tests                                                */
/* ========================================================================= */

TEST(parse_empty_string) {
    mcnp_model_t* model = mcnp_load_string("", 0);
    ASSERT_NULL(model);
}

TEST(parse_missing_surface) {
    /* Cell references surface 99 which doesn't exist */
    const char* input =
        "Test missing surf\n"
        "1 1 -1.0 -99\n"
        "2 0 99\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = mcnp_load_string(input, strlen(input));
    ASSERT_NULL(model);
}

TEST(parse_malformed_cell_missing_density) {
    const char* input =
        "Test malformed cell\n"
        "1 1\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = mcnp_load_string(input, strlen(input));
    ASSERT_NULL(model);
}

TEST(parse_malformed_surface_missing_mnemonic) {
    const char* input =
        "Test malformed surface\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = mcnp_load_string(input, strlen(input));
    ASSERT_NULL(model);
}

TEST_MAIN()
