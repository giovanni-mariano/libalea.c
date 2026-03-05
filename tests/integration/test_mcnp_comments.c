// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_mcnp_comments.c - MCNP cell comment preservation tests
 *
 * Tests "C" comment lines before cells and inline "$" comments on the last
 * line of cell cards, through parse, conversion, export, and roundtrip.
 */

#include "alea_test.h"
#include "alea.h"
#include "alea_mcnp.h"
#include "core/alea_system.h"
#include "core/alea_export.h"
#include <string.h>
#include <stdio.h>

/* ========================================================================= */
/* Helper                                                                    */
/* ========================================================================= */

static mcnp_model_t* parse_mcnp(const char* input) {
    mcnp_model_t* model = mcnp_load_string(input, strlen(input));
    if (model) alea_build_universe_index(model->sys);
    return model;
}

/* ========================================================================= */
/* "C" comment line tests                                                    */
/* ========================================================================= */

TEST(comment_single_line) {
    const char* input =
        "Test comments\n"
        "c Fuel pin region\n"
        "1 1 -10.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";

    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);

    const alea_cell_entry_t* cell = &model->sys->cells.data[0];
    ASSERT_NOT_NULL(cell->comments);
    ASSERT(strstr(cell->comments, "Fuel pin region") != NULL);

    // Cell 2 has no comment
    const alea_cell_entry_t* cell2 = &model->sys->cells.data[1];
    ASSERT_NULL(cell2->comments);

    mcnp_model_destroy(model);
}

TEST(comment_multi_line) {
    const char* input =
        "Test multi-line comments\n"
        "c Fuel pin inner region\n"
        "c Material: UO2\n"
        "c Density: 10.0 g/cm3\n"
        "1 1 -10.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";

    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);

    const alea_cell_entry_t* cell = &model->sys->cells.data[0];
    ASSERT_NOT_NULL(cell->comments);
    ASSERT(strstr(cell->comments, "Fuel pin") != NULL);
    ASSERT(strstr(cell->comments, "Material: UO2") != NULL);
    ASSERT(strstr(cell->comments, "Density: 10.0") != NULL);

    mcnp_model_destroy(model);
}

TEST(comment_per_cell) {
    const char* input =
        "Test per-cell comments\n"
        "c First cell - fuel\n"
        "1 1 -10.0 -1\n"
        "c Second cell - clad\n"
        "2 1 -6.5 1 -2\n"
        "3 0 2\n"
        "\n"
        "1 SO 5.0\n"
        "2 SO 10.0\n"
        "\n"
        "M1 92235.80c 1.0\n";

    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);

    const alea_cell_entry_t* c1 = &model->sys->cells.data[0];
    ASSERT_NOT_NULL(c1->comments);
    ASSERT(strstr(c1->comments, "First cell") != NULL);

    const alea_cell_entry_t* c2 = &model->sys->cells.data[1];
    ASSERT_NOT_NULL(c2->comments);
    ASSERT(strstr(c2->comments, "Second cell") != NULL);

    // Cell 3 has no comment
    const alea_cell_entry_t* c3 = &model->sys->cells.data[2];
    ASSERT_NULL(c3->comments);

    mcnp_model_destroy(model);
}

TEST(comment_case_insensitive) {
    const char* input =
        "Test case insensitive\n"
        "C Upper case comment\n"
        "1 1 -10.0 -1\n"
        "c lower case comment\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";

    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);

    const alea_cell_entry_t* c1 = &model->sys->cells.data[0];
    ASSERT_NOT_NULL(c1->comments);
    ASSERT(strstr(c1->comments, "Upper case") != NULL);

    const alea_cell_entry_t* c2 = &model->sys->cells.data[1];
    ASSERT_NOT_NULL(c2->comments);
    ASSERT(strstr(c2->comments, "lower case") != NULL);

    mcnp_model_destroy(model);
}

/* ========================================================================= */
/* Inline "$" comment tests                                                  */
/* ========================================================================= */

TEST(inline_comment_single_line) {
    const char* input =
        "Test inline\n"
        "1 1 -10.0 -1 $ fuel region\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";

    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);

    const alea_cell_entry_t* cell = &model->sys->cells.data[0];
    ASSERT_NOT_NULL(cell->inline_comment);
    ASSERT_STR_EQ(cell->inline_comment, "fuel region");

    // Cell 2 has no inline comment
    const alea_cell_entry_t* cell2 = &model->sys->cells.data[1];
    ASSERT_NULL(cell2->inline_comment);

    mcnp_model_destroy(model);
}

TEST(inline_comment_continuation_last_line) {
    // Only the $ on the LAST physical line should be captured
    const char* input =
        "Test continuation\n"
        "1 1 -10.0 -1 $ first line comment\n"
        "     2 -3 $ last line comment\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "2 SO 10.0\n"
        "3 SO 15.0\n"
        "\n"
        "M1 92235.80c 1.0\n";

    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);

    const alea_cell_entry_t* cell = &model->sys->cells.data[0];
    ASSERT_NOT_NULL(cell->inline_comment);
    ASSERT_STR_EQ(cell->inline_comment, "last line comment");

    mcnp_model_destroy(model);
}

TEST(inline_comment_continuation_no_dollar_on_last) {
    // If the last continuation line has no $, inline_comment should be NULL
    const char* input =
        "Test no dollar last\n"
        "1 1 -10.0 -1 $ first line only\n"
        "     2 -3\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "2 SO 10.0\n"
        "3 SO 15.0\n"
        "\n"
        "M1 92235.80c 1.0\n";

    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);

    const alea_cell_entry_t* cell = &model->sys->cells.data[0];
    ASSERT_NULL(cell->inline_comment);

    mcnp_model_destroy(model);
}

/* ========================================================================= */
/* Both comment types together                                               */
/* ========================================================================= */

TEST(both_comments) {
    const char* input =
        "Test both\n"
        "c Shielding layer\n"
        "c Steel structure\n"
        "1 1 -8.0 -1 $ steel\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 26056.80c 1.0\n";

    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);

    const alea_cell_entry_t* cell = &model->sys->cells.data[0];
    ASSERT_NOT_NULL(cell->comments);
    ASSERT(strstr(cell->comments, "Shielding layer") != NULL);
    ASSERT(strstr(cell->comments, "Steel structure") != NULL);
    ASSERT_NOT_NULL(cell->inline_comment);
    ASSERT_STR_EQ(cell->inline_comment, "steel");

    mcnp_model_destroy(model);
}

TEST(no_comments) {
    const char* input =
        "Test no comments\n"
        "1 1 -10.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";

    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);

    const alea_cell_entry_t* cell = &model->sys->cells.data[0];
    ASSERT_NULL(cell->comments);
    ASSERT_NULL(cell->inline_comment);

    mcnp_model_destroy(model);
}

/* ========================================================================= */
/* Roundtrip tests                                                           */
/* ========================================================================= */

TEST(comment_roundtrip) {
    const char* input =
        "Test roundtrip\n"
        "c Important shielding cell\n"
        "c Material: Iron\n"
        "1 1 -8.0 -1 $ steel layer\n"
        "2 0 1\n"
        "\n"
        "1 SO 10.0\n"
        "\n"
        "M1 26056.80c 1.0\n";

    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);

    const char* tmpfile = "test_comments_rt_tmp.mcnp";
    int rc = mcnp_export(model, tmpfile);
    ASSERT_EQ(rc, 0);
    mcnp_model_destroy(model);

    // Re-parse and verify comments preserved
    mcnp_model_t* model2 = mcnp_load(tmpfile);
    ASSERT_NOT_NULL(model2);
    alea_build_universe_index(model2->sys);

    const alea_cell_entry_t* cell = &model2->sys->cells.data[0];
    ASSERT_NOT_NULL(cell->comments);
    ASSERT(strstr(cell->comments, "shielding") != NULL);
    ASSERT(strstr(cell->comments, "Iron") != NULL);
    ASSERT_NOT_NULL(cell->inline_comment);
    ASSERT_STR_EQ(cell->inline_comment, "steel layer");

    // Cell 2 should have no comments after roundtrip
    const alea_cell_entry_t* cell2 = &model2->sys->cells.data[1];
    ASSERT_NULL(cell2->comments);
    ASSERT_NULL(cell2->inline_comment);

    mcnp_model_destroy(model2);
    remove(tmpfile);
}

TEST(comment_roundtrip_multiple_cells) {
    const char* input =
        "Test multi-cell roundtrip\n"
        "c Cell 1: fuel\n"
        "1 1 -10.0 -1 $ UO2 fuel\n"
        "c Cell 2: clad\n"
        "2 1 -6.5 1 -2 $ zircaloy\n"
        "3 0 2\n"
        "\n"
        "1 SO 5.0\n"
        "2 SO 6.0\n"
        "\n"
        "M1 92235.80c 1.0\n";

    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);

    const char* tmpfile = "test_comments_multi_rt_tmp.mcnp";
    int rc = mcnp_export(model, tmpfile);
    ASSERT_EQ(rc, 0);
    mcnp_model_destroy(model);

    mcnp_model_t* model2 = mcnp_load(tmpfile);
    ASSERT_NOT_NULL(model2);
    alea_build_universe_index(model2->sys);

    const alea_cell_entry_t* c1 = &model2->sys->cells.data[0];
    ASSERT_NOT_NULL(c1->comments);
    ASSERT(strstr(c1->comments, "fuel") != NULL);
    ASSERT_NOT_NULL(c1->inline_comment);
    ASSERT_STR_EQ(c1->inline_comment, "UO2 fuel");

    const alea_cell_entry_t* c2 = &model2->sys->cells.data[1];
    ASSERT_NOT_NULL(c2->comments);
    ASSERT(strstr(c2->comments, "clad") != NULL);
    ASSERT_NOT_NULL(c2->inline_comment);
    ASSERT_STR_EQ(c2->inline_comment, "zircaloy");

    const alea_cell_entry_t* c3 = &model2->sys->cells.data[2];
    ASSERT_NULL(c3->comments);
    ASSERT_NULL(c3->inline_comment);

    mcnp_model_destroy(model2);
    remove(tmpfile);
}

TEST(comment_roundtrip_system_export) {
    // Test via mcnp_export_system (core path, no model)
    const char* input =
        "Test system export\n"
        "c Core region\n"
        "1 1 -10.0 -1 $ core\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";

    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);

    const char* tmpfile = "test_comments_sys_rt_tmp.mcnp";
    int rc = mcnp_export_system(model->sys, tmpfile);
    ASSERT_EQ(rc, 0);
    mcnp_model_destroy(model);

    mcnp_model_t* model2 = mcnp_load(tmpfile);
    ASSERT_NOT_NULL(model2);
    alea_build_universe_index(model2->sys);

    const alea_cell_entry_t* cell = &model2->sys->cells.data[0];
    ASSERT_NOT_NULL(cell->comments);
    ASSERT(strstr(cell->comments, "Core region") != NULL);
    ASSERT_NOT_NULL(cell->inline_comment);
    ASSERT_STR_EQ(cell->inline_comment, "core");

    mcnp_model_destroy(model2);
    remove(tmpfile);
}

TEST_MAIN()
