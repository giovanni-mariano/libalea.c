// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_core.c - Core integration tests using the public alea API
 */

#include "alea_test.h"
#include "alea.h"

/* ------------------------------------------------------------------------- */
/* System Lifecycle Tests                                                     */
/* ------------------------------------------------------------------------- */

TEST(lifecycle_create_destroy) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    alea_destroy(sys);
}

TEST(lifecycle_reset) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int si = alea_box_surface(sys, 0, -1, 1, -1, 1, -1, 1);
    alea_node_id_t box = alea_halfspace(sys, si, -1);
    ASSERT_NE(box, ALEA_NODE_ID_INVALID);

    int m1 = alea_add_material(sys, 1);
    alea_add_cell(sys, 1, box, m1, 1.0, 0);

    alea_reset(sys);
    ASSERT_EQ(alea_cell_count(sys), 0);
    ASSERT_EQ(alea_surface_count(sys), 0);

    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Primitive Creation Tests                                                   */
/* ------------------------------------------------------------------------- */

TEST(primitive_box) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int si = alea_box_surface(sys, 0, -1, 1, -1, 1, -1, 1);
    alea_node_id_t box = alea_halfspace(sys, si, -1);
    ASSERT_NE(box, ALEA_NODE_ID_INVALID);

    int m5 = alea_add_material(sys, 5);
    int cell_idx = alea_add_cell(sys, 100, box, m5, 1.0, 0);
    ASSERT(cell_idx >= 0);

    ASSERT_EQ(alea_build_universe_index(sys), 0);

    ASSERT_EQ(alea_find_cell(sys, 0, 0, 0), cell_idx);
    ASSERT_EQ(alea_find_cell(sys, 2, 0, 0), -1);
    ASSERT_EQ(alea_material_at(sys, 0, 0, 0), 5);

    alea_destroy(sys);
}

TEST(primitive_sphere) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int si = alea_sphere_surface(sys, 0, 0, 0, 0, 2.0);
    alea_node_id_t s = alea_halfspace(sys, si, -1);
    ASSERT_NE(s, ALEA_NODE_ID_INVALID);

    int m1 = alea_add_material(sys, 1);
    int cell = alea_add_cell(sys, 1, s, m1, 1.0, 0);
    ASSERT(cell >= 0);

    ASSERT_EQ(alea_build_universe_index(sys), 0);

    ASSERT_EQ(alea_find_cell(sys, 0, 0, 0), cell);
    ASSERT_EQ(alea_find_cell(sys, 5, 0, 0), -1);

    alea_destroy(sys);
}

TEST(primitive_cylinder) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int si = alea_cylinder_z_surface(sys, 0, 0, 0, 2.0);
    alea_node_id_t cyl = alea_halfspace(sys, si, -1);
    ASSERT_NE(cyl, ALEA_NODE_ID_INVALID);

    int m1 = alea_add_material(sys, 1);
    int cell = alea_add_cell(sys, 1, cyl, m1, 1.0, 0);
    ASSERT(cell >= 0);

    ASSERT_EQ(alea_build_universe_index(sys), 0);

    ASSERT_EQ(alea_find_cell(sys, 0, 0, 0), cell);
    ASSERT_EQ(alea_find_cell(sys, 1, 0, 100), cell);
    ASSERT_EQ(alea_find_cell(sys, 5, 0, 0), -1);

    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* CSG Operation Tests                                                        */
/* ------------------------------------------------------------------------- */

TEST(alea_union) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int si1 = alea_sphere_surface(sys, 0, 0, 0, 0, 2.0);
    int si2 = alea_sphere_surface(sys, 0, 3, 0, 0, 2.0);
    alea_node_id_t s1 = alea_halfspace(sys, si1, -1);
    alea_node_id_t s2 = alea_halfspace(sys, si2, -1);

    alea_node_id_t uni = alea_union(sys, s1, s2);
    ASSERT_NE(uni, ALEA_NODE_ID_INVALID);

    int m1 = alea_add_material(sys, 1);
    int cell = alea_add_cell(sys, 1, uni, m1, 1.0, 0);
    ASSERT(cell >= 0);

    ASSERT_EQ(alea_build_universe_index(sys), 0);

    ASSERT_EQ(alea_find_cell(sys, 0, 0, 0), cell);
    ASSERT_EQ(alea_find_cell(sys, 3, 0, 0), cell);
    ASSERT_EQ(alea_find_cell(sys, 10, 0, 0), -1);

    alea_destroy(sys);
}

TEST(alea_intersection) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int si1 = alea_sphere_surface(sys, 0, 0, 0, 0, 3.0);
    int si2 = alea_sphere_surface(sys, 0, 2, 0, 0, 3.0);
    alea_node_id_t s1 = alea_halfspace(sys, si1, -1);
    alea_node_id_t s2 = alea_halfspace(sys, si2, -1);

    alea_node_id_t inter = alea_intersection(sys, s1, s2);
    ASSERT_NE(inter, ALEA_NODE_ID_INVALID);

    int m1 = alea_add_material(sys, 1);
    int cell = alea_add_cell(sys, 1, inter, m1, 1.0, 0);
    ASSERT(cell >= 0);

    ASSERT_EQ(alea_build_universe_index(sys), 0);

    ASSERT_EQ(alea_find_cell(sys, 1, 0, 0), cell);
    ASSERT_EQ(alea_find_cell(sys, -2, 0, 0), -1);

    alea_destroy(sys);
}

TEST(alea_difference) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int si1 = alea_sphere_surface(sys, 0, 0, 0, 0, 5.0);
    int si2 = alea_sphere_surface(sys, 0, 0, 0, 0, 2.0);
    alea_node_id_t outer = alea_halfspace(sys, si1, -1);
    alea_node_id_t inner = alea_halfspace(sys, si2, -1);

    alea_node_id_t diff = alea_difference(sys, outer, inner);
    ASSERT_NE(diff, ALEA_NODE_ID_INVALID);

    int m1 = alea_add_material(sys, 1);
    int cell = alea_add_cell(sys, 1, diff, m1, 1.0, 0);
    ASSERT(cell >= 0);

    ASSERT_EQ(alea_build_universe_index(sys), 0);

    ASSERT_EQ(alea_find_cell(sys, 3, 0, 0), cell);
    ASSERT_EQ(alea_find_cell(sys, 0, 0, 0), -1);

    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Cell Query Tests                                                           */
/* ------------------------------------------------------------------------- */

TEST(find_all_cells) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int si = alea_sphere_surface(sys, 0, 0, 0, 0, 2.0);
    alea_node_id_t s = alea_halfspace(sys, si, -1);

    int m99 = alea_add_material(sys, 99);
    int cell_idx = alea_add_cell(sys, 42, s, m99, 1.0, 0);
    ASSERT(cell_idx >= 0);

    ASSERT_EQ(alea_build_universe_index(sys), 0);

    alea_cell_hit_t hits[16];
    int count = alea_find_all_cells(sys, 0, 0, 0, hits, 16);
    ASSERT_EQ(count, 1);
    ASSERT_EQ(hits[0].cell_id, 42);
    ASSERT_EQ(hits[0].material_id, 99);

    alea_destroy(sys);
}

TEST(overlap_detection) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int si1 = alea_box_surface(sys, 0, 0, 2, 0, 2, 0, 2);
    int si2 = alea_box_surface(sys, 0, 1, 3, 1, 3, 1, 3);
    alea_node_id_t a = alea_halfspace(sys, si1, -1);
    alea_node_id_t b = alea_halfspace(sys, si2, -1);

    int m1 = alea_add_material(sys, 1);
    ASSERT(alea_add_cell(sys, 1, a, m1, 1.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 2, b, m1, 1.0, 0) >= 0);

    int pairs[16];
    int count = alea_find_overlaps(sys, pairs, 8);
    ASSERT(count >= 1);

    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Universe and Flattening Tests                                              */
/* ------------------------------------------------------------------------- */

TEST(universe_flatten) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int si = alea_box_surface(sys, 0, -1, 1, -1, 1, -1, 1);
    alea_node_id_t box = alea_halfspace(sys, si, -1);

    int m1 = alea_add_material(sys, 1);
    ASSERT(alea_add_cell(sys, 10, box, m1, 1.0, 1) >= 0);

    ASSERT_EQ(alea_build_universe_index(sys), 0);
    ASSERT(alea_flatten(sys, 1) >= 0);
    ASSERT(alea_cell_count(sys) >= 1);

    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Export Tests                                                               */
/* ------------------------------------------------------------------------- */

TEST(export_mcnp) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int si = alea_box_surface(sys, 0, 0, 1, 0, 1, 0, 1);
    alea_node_id_t box = alea_halfspace(sys, si, -1);

    int m1 = alea_add_material(sys, 1);
    ASSERT(alea_add_cell(sys, 1, box, m1, 1.0, 0) >= 0);

    FILE* f = tmpfile();
    ASSERT_NOT_NULL(f);
    int rc = alea_export_mcnp_stream(sys, f);
    ASSERT_EQ(rc, 0);
    fclose(f);

    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Void Generation Tests                                                      */
/* ------------------------------------------------------------------------- */

TEST(void_generation) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int si = alea_box_surface(sys, 0, -1, 1, -1, 1, -1, 1);
    alea_node_id_t box = alea_halfspace(sys, si, -1);

    int m1 = alea_add_material(sys, 1);
    ASSERT(alea_add_cell(sys, 1, box, m1, 1.0, 0) >= 0);

    ASSERT_EQ(alea_build_universe_index(sys), 0);

    alea_bbox_t bounds = {-5, 5, -5, 5, -5, 5};
    void_result_t* vr = alea_void_generate(sys, &bounds);
    ASSERT_NOT_NULL(vr);

    alea_void_free(vr);
    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Stress Tests                                                               */
/* ------------------------------------------------------------------------- */

TEST(stress_many_cells) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Register materials 1..99; material 0 is void */
    int mat_idx[100];
    mat_idx[0] = ALEA_MATERIAL_VOID;
    for (int i = 1; i < 100; ++i)
        mat_idx[i] = alea_add_material(sys, i);

    for (int i = 0; i < 100; ++i) {
        int si = alea_sphere_surface(sys, 0, i * 3.0, 0, 0, 1.0);
        alea_node_id_t s = alea_halfspace(sys, si, -1);
        ASSERT_NE(s, ALEA_NODE_ID_INVALID);
        ASSERT(alea_add_cell(sys, i + 1, s, mat_idx[i], 1.0, 0) >= 0);
    }

    ASSERT_EQ(alea_cell_count(sys), 100);

    alea_destroy(sys);
}

TEST_MAIN()
