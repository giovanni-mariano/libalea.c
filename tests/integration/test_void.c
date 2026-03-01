// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_void.c - Void generation integration tests using the public alea API
 */

#include "alea_test.h"
#include "alea.h"

/* Reduce octree depth for fast tests (default 8 is too slow for CI) */
static void set_fast_void_config(alea_system_t* sys) {
    alea_config_t cfg = alea_get_config(sys);
    cfg.void_max_depth = 4;
    cfg.void_min_size = 0.5;
    alea_set_config(sys, &cfg);
}

/* Helper: create a system with a single box cell at [-1,1]^3 in universe 0 */
static alea_system_t* make_box_system(void) {
    alea_system_t* sys = alea_create();
    if (!sys) return NULL;

    set_fast_void_config(sys);

    int si = alea_box_surface(sys, 0, -1, 1, -1, 1, -1, 1);
    alea_node_id_t box = alea_halfspace(sys, si, -1);

    int m1 = alea_add_material(sys, 1);
    alea_add_cell(sys, 1, box, m1, 1.0, 0);

    alea_build_universe_index(sys);
    return sys;
}

/* ------------------------------------------------------------------------- */
/* Core correctness                                                           */
/* ------------------------------------------------------------------------- */

TEST(void_single_box) {
    alea_system_t* sys = make_box_system();
    ASSERT_NOT_NULL(sys);

    alea_bbox_t bounds = {-5, 5, -5, 5, -5, 5};
    void_result_t* vr = alea_void_generate(sys, &bounds);
    ASSERT_NOT_NULL(vr);

    /* There IS void outside the box */
    ASSERT(alea_void_count(vr) > 0);

    /* Each returned bbox should be valid and within search bounds */
    for (size_t i = 0; i < alea_void_count(vr); i++) {
        alea_bbox_t box;
        ASSERT_EQ(alea_void_get(vr, i, &box), 0);

        /* Positive volume */
        double vol = (box.max_x - box.min_x) *
                     (box.max_y - box.min_y) *
                     (box.max_z - box.min_z);
        ASSERT(vol > 0);

        /* Within search bounds (with tolerance for float rounding) */
        ASSERT(box.min_x >= bounds.min_x - 1e-6);
        ASSERT(box.max_x <= bounds.max_x + 1e-6);
        ASSERT(box.min_y >= bounds.min_y - 1e-6);
        ASSERT(box.max_y <= bounds.max_y + 1e-6);
        ASSERT(box.min_z >= bounds.min_z - 1e-6);
        ASSERT(box.max_z <= bounds.max_z + 1e-6);
    }

    alea_void_free(vr);
    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* CSG node creation                                                          */
/* ------------------------------------------------------------------------- */

TEST(void_to_node) {
    alea_system_t* sys = make_box_system();
    ASSERT_NOT_NULL(sys);

    alea_bbox_t bounds = {-5, 5, -5, 5, -5, 5};
    void_result_t* vr = alea_void_generate(sys, &bounds);
    ASSERT_NOT_NULL(vr);
    ASSERT(alea_void_count(vr) > 0);

    alea_node_id_t node = alea_void_to_node(sys, vr);
    ASSERT_NE(node, ALEA_NODE_ID_INVALID);

    alea_void_free(vr);
    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Cell registration                                                          */
/* ------------------------------------------------------------------------- */

TEST(void_add_cells) {
    alea_system_t* sys = make_box_system();
    ASSERT_NOT_NULL(sys);

    size_t before = alea_cell_count(sys);

    alea_bbox_t bounds = {-5, 5, -5, 5, -5, 5};
    void_result_t* vr = alea_void_generate(sys, &bounds);
    ASSERT_NOT_NULL(vr);

    size_t void_count = alea_void_count(vr);
    ASSERT(void_count > 0);

    int added = alea_void_add_cells(sys, vr);
    ASSERT(added > 0);
    ASSERT_EQ((size_t)added, void_count);
    ASSERT_EQ(alea_cell_count(sys), before + void_count);

    alea_void_free(vr);
    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Merge reduces count                                                        */
/* ------------------------------------------------------------------------- */

TEST(void_merge) {
    alea_system_t* sys = make_box_system();
    ASSERT_NOT_NULL(sys);

    alea_bbox_t bounds = {-5, 5, -5, 5, -5, 5};
    void_result_t* vr = alea_void_generate(sys, &bounds);
    ASSERT_NOT_NULL(vr);

    size_t count_before = alea_void_count(vr);
    ASSERT(count_before > 0);

    int merged = alea_void_merge(sys, vr);
    ASSERT(merged > 0);

    size_t count_after = alea_void_count(vr);
    ASSERT(count_after <= count_before);
    ASSERT(count_after > 0);

    alea_void_free(vr);
    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* No void when cell covers entire bounds                                     */
/* ------------------------------------------------------------------------- */

TEST(void_no_void) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    set_fast_void_config(sys);

    /* Box [-10,10]^3 fully encloses bounds [-5,5]^3 */
    int si = alea_box_surface(sys, 0, -10, 10, -10, 10, -10, 10);
    alea_node_id_t box = alea_halfspace(sys, si, -1);

    int m1 = alea_add_material(sys, 1);
    alea_add_cell(sys, 1, box, m1, 1.0, 0);
    alea_build_universe_index(sys);

    alea_bbox_t bounds = {-5, 5, -5, 5, -5, 5};
    void_result_t* vr = alea_void_generate(sys, &bounds);
    ASSERT_NOT_NULL(vr);

    ASSERT_EQ(alea_void_count(vr), 0);

    alea_void_free(vr);
    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Non-box geometry (sphere)                                                  */
/* ------------------------------------------------------------------------- */

TEST(void_sphere) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    set_fast_void_config(sys);

    int si = alea_sphere_surface(sys, 0, 0, 0, 0, 3.0);
    alea_node_id_t s = alea_halfspace(sys, si, -1);

    int m1 = alea_add_material(sys, 1);
    alea_add_cell(sys, 1, s, m1, 1.0, 0);
    alea_build_universe_index(sys);

    alea_bbox_t bounds = {-5, 5, -5, 5, -5, 5};
    void_result_t* vr = alea_void_generate(sys, &bounds);
    ASSERT_NOT_NULL(vr);

    /* Void in corners outside sphere */
    ASSERT(alea_void_count(vr) > 0);

    /* All boxes within bounds */
    for (size_t i = 0; i < alea_void_count(vr); i++) {
        alea_bbox_t box;
        ASSERT_EQ(alea_void_get(vr, i, &box), 0);
        ASSERT(box.min_x >= bounds.min_x - 1e-6);
        ASSERT(box.max_x <= bounds.max_x + 1e-6);
        ASSERT(box.min_y >= bounds.min_y - 1e-6);
        ASSERT(box.max_y <= bounds.max_y + 1e-6);
        ASSERT(box.min_z >= bounds.min_z - 1e-6);
        ASSERT(box.max_z <= bounds.max_z + 1e-6);
    }

    alea_void_free(vr);
    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Point verification: CSG correctness after void cell registration           */
/* ------------------------------------------------------------------------- */

TEST(void_point_verification) {
    alea_system_t* sys = make_box_system();
    ASSERT_NOT_NULL(sys);

    int mat_cell = alea_find_cell(sys, 0, 0, 0);
    ASSERT(mat_cell >= 0);

    /* Before void: point outside box is unresolved */
    ASSERT_EQ(alea_find_cell(sys, 3, 3, 3), -1);

    alea_bbox_t bounds = {-5, 5, -5, 5, -5, 5};
    void_result_t* vr = alea_void_generate(sys, &bounds);
    ASSERT_NOT_NULL(vr);
    ASSERT(alea_void_count(vr) > 0);

    alea_void_add_cells(sys, vr);
    alea_build_universe_index(sys);

    /* Point inside material cell is still found */
    int found = alea_find_cell(sys, 0, 0, 0);
    ASSERT(found >= 0);

    /* Point that was void is now covered by a void cell */
    int void_cell = alea_find_cell(sys, 3, 3, 3);
    ASSERT(void_cell >= 0);

    /* Verify it's a void material cell (MCNP convention: material_id=0 for void) */
    int material_id;
    ASSERT_EQ(alea_cell_get(sys, (size_t)void_cell, NULL, &material_id,
                             NULL, NULL, NULL, NULL), 0);
    ASSERT_EQ(material_id, 0);

    alea_void_free(vr);
    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Null safety                                                                */
/* ------------------------------------------------------------------------- */

TEST(void_null_safety) {
    /* NULL sys returns NULL */
    void_result_t* vr = alea_void_generate(NULL, NULL);
    ASSERT_NULL(vr);

    /* NULL bounds with valid sys uses auto-bounds (should not crash) */
    alea_system_t* sys = make_box_system();
    ASSERT_NOT_NULL(sys);

    vr = alea_void_generate(sys, NULL);
    ASSERT_NOT_NULL(vr);
    ASSERT(alea_void_count(vr) > 0);

    alea_void_free(vr);
    alea_destroy(sys);
}

TEST_MAIN()
