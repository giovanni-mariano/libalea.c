// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_cell_setters.c - Tests for cell property setters and cell removal
 */

#include "alea_test.h"
#include "alea.h"

/* Helper: create a system with one material (id=1) and one cell */
static alea_system_t* make_simple_system(void) {
    alea_system_t* sys = alea_create();
    if (!sys) return NULL;

    int mat = alea_add_material(sys, 1);  /* material index 0, id 1 */
    (void)mat;

    int si = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    alea_node_id_t inside = alea_halfspace(sys, si, -1);
    alea_add_cell(sys, 10, inside, 0, -2.5, 0);  /* cell id=10, mat idx=0, density=-2.5 (g/cm3) */

    return sys;
}

/* ========================================================================= */
/* alea_cell_set_material                                                    */
/* ========================================================================= */

TEST(set_material_basic) {
    alea_system_t* sys = make_simple_system();
    ASSERT_NOT_NULL(sys);

    /* Verify initial material */
    int mat_id = -1;
    alea_cell_get(sys, 0, NULL, &mat_id, NULL, NULL, NULL, NULL);
    ASSERT_EQ(mat_id, 1);

    /* Add a second material and switch to it */
    int mat2 = alea_add_material(sys, 5);  /* material index 1, id 5 */
    ASSERT(mat2 >= 0);

    int rc = alea_cell_set_material(sys, 0, mat2);
    ASSERT_EQ(rc, 0);

    alea_cell_get(sys, 0, NULL, &mat_id, NULL, NULL, NULL, NULL);
    ASSERT_EQ(mat_id, 5);

    alea_destroy(sys);
}

TEST(set_material_void) {
    alea_system_t* sys = make_simple_system();
    ASSERT_NOT_NULL(sys);

    int rc = alea_cell_set_material(sys, 0, ALEA_MATERIAL_VOID);
    ASSERT_EQ(rc, 0);

    int mat_id = -1;
    alea_cell_get(sys, 0, NULL, &mat_id, NULL, NULL, NULL, NULL);
    ASSERT_EQ(mat_id, 0);

    alea_destroy(sys);
}

TEST(set_material_out_of_range) {
    alea_system_t* sys = make_simple_system();
    ASSERT_NOT_NULL(sys);

    int rc = alea_cell_set_material(sys, 0, 999);
    ASSERT_EQ(rc, -1);

    /* Material unchanged */
    int mat_id = -1;
    alea_cell_get(sys, 0, NULL, &mat_id, NULL, NULL, NULL, NULL);
    ASSERT_EQ(mat_id, 1);

    alea_destroy(sys);
}

TEST(set_material_invalid_cell) {
    alea_system_t* sys = make_simple_system();
    ASSERT_NOT_NULL(sys);

    ASSERT_EQ(alea_cell_set_material(sys, -1, 0), -1);
    ASSERT_EQ(alea_cell_set_material(sys, 100, 0), -1);
    ASSERT_EQ(alea_cell_set_material(NULL, 0, 0), -1);

    alea_destroy(sys);
}

/* ========================================================================= */
/* alea_cell_set_density                                                     */
/* ========================================================================= */

TEST(set_density_mass) {
    alea_system_t* sys = make_simple_system();
    ASSERT_NOT_NULL(sys);

    /* Initial: -2.5 -> density=2.5, mass density */
    alea_cell_info_t info;
    alea_cell_get_info(sys, 0, &info);
    ASSERT_NEAR(info.density, 2.5, 1e-10);
    ASSERT_TRUE(info.is_mass_density);

    /* Set to atom density (positive) */
    int rc = alea_cell_set_density(sys, 0, 0.05);
    ASSERT_EQ(rc, 0);

    alea_cell_get_info(sys, 0, &info);
    ASSERT_NEAR(info.density, 0.05, 1e-10);
    ASSERT_FALSE(info.is_mass_density);

    /* Set back to mass density (negative) */
    rc = alea_cell_set_density(sys, 0, -7.8);
    ASSERT_EQ(rc, 0);

    alea_cell_get_info(sys, 0, &info);
    ASSERT_NEAR(info.density, 7.8, 1e-10);
    ASSERT_TRUE(info.is_mass_density);

    alea_destroy(sys);
}

TEST(set_density_invalid_cell) {
    alea_system_t* sys = make_simple_system();
    ASSERT_NOT_NULL(sys);

    ASSERT_EQ(alea_cell_set_density(sys, -1, 1.0), -1);
    ASSERT_EQ(alea_cell_set_density(sys, 100, 1.0), -1);
    ASSERT_EQ(alea_cell_set_density(NULL, 0, 1.0), -1);

    alea_destroy(sys);
}

/* ========================================================================= */
/* alea_cell_set_universe                                                    */
/* ========================================================================= */

TEST(set_universe_basic) {
    alea_system_t* sys = make_simple_system();
    ASSERT_NOT_NULL(sys);

    /* Initial: universe 0 */
    int univ = -1;
    alea_cell_get(sys, 0, NULL, NULL, NULL, &univ, NULL, NULL);
    ASSERT_EQ(univ, 0);

    int rc = alea_cell_set_universe(sys, 0, 5);
    ASSERT_EQ(rc, 0);

    alea_cell_get(sys, 0, NULL, NULL, NULL, &univ, NULL, NULL);
    ASSERT_EQ(univ, 5);

    alea_destroy(sys);
}

TEST(set_universe_invalid_cell) {
    alea_system_t* sys = make_simple_system();
    ASSERT_NOT_NULL(sys);

    ASSERT_EQ(alea_cell_set_universe(sys, -1, 1), -1);
    ASSERT_EQ(alea_cell_set_universe(sys, 100, 1), -1);
    ASSERT_EQ(alea_cell_set_universe(NULL, 0, 1), -1);

    alea_destroy(sys);
}

/* ========================================================================= */
/* alea_cell_remove                                                          */
/* ========================================================================= */

TEST(remove_single_cell) {
    alea_system_t* sys = make_simple_system();
    ASSERT_NOT_NULL(sys);

    ASSERT_EQ(alea_cell_count(sys), 1);

    int rc = alea_cell_remove(sys, 0);
    ASSERT_EQ(rc, 0);
    ASSERT_EQ(alea_cell_count(sys), 0);

    alea_destroy(sys);
}

TEST(remove_middle_cell) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);

    /* Create 3 cells with ids 10, 20, 30 */
    int s1 = alea_sphere_surface(sys, 0, 0, 0, 0, 1.0);
    int s2 = alea_sphere_surface(sys, 0, 10, 0, 0, 1.0);
    int s3 = alea_sphere_surface(sys, 0, 20, 0, 0, 1.0);
    alea_add_cell(sys, 10, alea_halfspace(sys, s1, -1), mat, -1.0, 0);
    alea_add_cell(sys, 20, alea_halfspace(sys, s2, -1), mat, -2.0, 0);
    alea_add_cell(sys, 30, alea_halfspace(sys, s3, -1), mat, -3.0, 0);

    ASSERT_EQ(alea_cell_count(sys), 3);

    /* Remove middle cell (id=20, index=1) */
    int rc = alea_cell_remove(sys, 1);
    ASSERT_EQ(rc, 0);
    ASSERT_EQ(alea_cell_count(sys), 2);

    /* Remaining: id=10 at index 0, id=30 at index 1 */
    int cell_id = -1;
    alea_cell_get(sys, 0, &cell_id, NULL, NULL, NULL, NULL, NULL);
    ASSERT_EQ(cell_id, 10);

    alea_cell_get(sys, 1, &cell_id, NULL, NULL, NULL, NULL, NULL);
    ASSERT_EQ(cell_id, 30);

    /* Hashmap lookup still works */
    ASSERT_EQ(alea_cell_find(sys, 10), 0);
    ASSERT_EQ(alea_cell_find(sys, 30), 1);
    ASSERT_EQ(alea_cell_find(sys, 20), -1);  /* removed */

    alea_destroy(sys);
}

TEST(remove_first_cell) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    int s1 = alea_sphere_surface(sys, 0, 0, 0, 0, 1.0);
    int s2 = alea_sphere_surface(sys, 0, 10, 0, 0, 1.0);
    alea_add_cell(sys, 10, alea_halfspace(sys, s1, -1), mat, -1.0, 0);
    alea_add_cell(sys, 20, alea_halfspace(sys, s2, -1), mat, -2.0, 0);

    int rc = alea_cell_remove(sys, 0);
    ASSERT_EQ(rc, 0);
    ASSERT_EQ(alea_cell_count(sys), 1);

    int cell_id = -1;
    alea_cell_get(sys, 0, &cell_id, NULL, NULL, NULL, NULL, NULL);
    ASSERT_EQ(cell_id, 20);

    ASSERT_EQ(alea_cell_find(sys, 20), 0);
    ASSERT_EQ(alea_cell_find(sys, 10), -1);

    alea_destroy(sys);
}

TEST(remove_last_cell) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    int s1 = alea_sphere_surface(sys, 0, 0, 0, 0, 1.0);
    int s2 = alea_sphere_surface(sys, 0, 10, 0, 0, 1.0);
    alea_add_cell(sys, 10, alea_halfspace(sys, s1, -1), mat, -1.0, 0);
    alea_add_cell(sys, 20, alea_halfspace(sys, s2, -1), mat, -2.0, 0);

    int rc = alea_cell_remove(sys, 1);
    ASSERT_EQ(rc, 0);
    ASSERT_EQ(alea_cell_count(sys), 1);

    int cell_id = -1;
    alea_cell_get(sys, 0, &cell_id, NULL, NULL, NULL, NULL, NULL);
    ASSERT_EQ(cell_id, 10);

    ASSERT_EQ(alea_cell_find(sys, 10), 0);
    ASSERT_EQ(alea_cell_find(sys, 20), -1);

    alea_destroy(sys);
}

TEST(remove_invalid_index) {
    alea_system_t* sys = make_simple_system();
    ASSERT_NOT_NULL(sys);

    ASSERT_EQ(alea_cell_remove(sys, -1), -1);
    ASSERT_EQ(alea_cell_remove(sys, 100), -1);
    ASSERT_EQ(alea_cell_remove(NULL, 0), -1);

    /* Cell still there */
    ASSERT_EQ(alea_cell_count(sys), 1);

    alea_destroy(sys);
}

TEST(remove_then_add) {
    alea_system_t* sys = make_simple_system();
    ASSERT_NOT_NULL(sys);

    /* Remove the only cell */
    alea_cell_remove(sys, 0);
    ASSERT_EQ(alea_cell_count(sys), 0);

    /* Add a new cell */
    int si = alea_sphere_surface(sys, 0, 0, 0, 0, 3.0);
    alea_node_id_t inside = alea_halfspace(sys, si, -1);
    int idx = alea_add_cell(sys, 99, inside, ALEA_MATERIAL_VOID, 0.0, 0);
    ASSERT(idx >= 0);
    ASSERT_EQ(alea_cell_count(sys), 1);

    int cell_id = -1;
    alea_cell_get(sys, 0, &cell_id, NULL, NULL, NULL, NULL, NULL);
    ASSERT_EQ(cell_id, 99);

    ASSERT_EQ(alea_cell_find(sys, 99), 0);

    alea_destroy(sys);
}

TEST_MAIN()
