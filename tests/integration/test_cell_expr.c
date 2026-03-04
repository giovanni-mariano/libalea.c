// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_cell_expr.c - Tests for alea_cell_expr() CSG expression stringifier
 */

#include "alea_test.h"
#include "alea.h"

#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------------- */
/* Single primitive                                                           */
/* ------------------------------------------------------------------------- */

TEST(expr_single_primitive) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int si = alea_sphere_surface(sys, 10, 0, 0, 0, 5.0);
    alea_node_id_t inside = alea_halfspace(sys, si, -1);

    int m = alea_add_material(sys, 1);
    int ci = alea_add_cell(sys, 1, inside, m, 1.0, 0);
    ASSERT(ci >= 0);

    char* expr = alea_cell_expr(sys, (size_t)ci, ":", " ", "#");
    ASSERT_NOT_NULL(expr);
    ASSERT_STR_EQ(expr, "-10");
    free(expr);

    alea_destroy(sys);
}

TEST(expr_positive_sense) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int si = alea_sphere_surface(sys, 10, 0, 0, 0, 5.0);
    alea_node_id_t outside = alea_halfspace(sys, si, +1);

    int m = alea_add_material(sys, 1);
    int ci = alea_add_cell(sys, 1, outside, m, 1.0, 0);
    ASSERT(ci >= 0);

    char* expr = alea_cell_expr(sys, (size_t)ci, ":", " ", "#");
    ASSERT_NOT_NULL(expr);
    ASSERT_STR_EQ(expr, "10");
    free(expr);

    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Intersection                                                               */
/* ------------------------------------------------------------------------- */

TEST(expr_intersection) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Box: -S1 & -S2  (inside sphere1 AND inside sphere2) */
    int s1 = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    int s2 = alea_sphere_surface(sys, 2, 1, 0, 0, 3.0);
    alea_node_id_t n1 = alea_halfspace(sys, s1, -1);
    alea_node_id_t n2 = alea_halfspace(sys, s2, -1);
    alea_node_id_t root = alea_intersection(sys, n1, n2);

    int m = alea_add_material(sys, 1);
    int ci = alea_add_cell(sys, 1, root, m, 1.0, 0);
    ASSERT(ci >= 0);

    char* expr = alea_cell_expr(sys, (size_t)ci, ":", " ", "#");
    ASSERT_NOT_NULL(expr);
    ASSERT_STR_EQ(expr, "-1 -2");
    free(expr);

    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Union                                                                      */
/* ------------------------------------------------------------------------- */

TEST(expr_union) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int s1 = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    int s2 = alea_sphere_surface(sys, 2, 10, 0, 0, 3.0);
    alea_node_id_t n1 = alea_halfspace(sys, s1, -1);
    alea_node_id_t n2 = alea_halfspace(sys, s2, -1);
    alea_node_id_t root = alea_union(sys, n1, n2);

    int m = alea_add_material(sys, 1);
    int ci = alea_add_cell(sys, 1, root, m, 1.0, 0);
    ASSERT(ci >= 0);

    char* mcnp = alea_cell_expr(sys, (size_t)ci, ":", " ", "#");
    ASSERT_NOT_NULL(mcnp);
    ASSERT_STR_EQ(mcnp, "-1:-2");
    free(mcnp);

    /* OpenMC style */
    char* omc = alea_cell_expr(sys, (size_t)ci, " | ", " ", "~");
    ASSERT_NOT_NULL(omc);
    ASSERT_STR_EQ(omc, "-1 | -2");
    free(omc);

    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Union inside intersection needs parens                                     */
/* ------------------------------------------------------------------------- */

TEST(expr_union_in_intersection) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* (-S1 | -S2) & -S3 */
    int s1 = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    int s2 = alea_sphere_surface(sys, 2, 10, 0, 0, 3.0);
    int s3 = alea_sphere_surface(sys, 3, 5, 0, 0, 8.0);
    alea_node_id_t n1 = alea_halfspace(sys, s1, -1);
    alea_node_id_t n2 = alea_halfspace(sys, s2, -1);
    alea_node_id_t n3 = alea_halfspace(sys, s3, -1);
    alea_node_id_t u = alea_union(sys, n1, n2);
    alea_node_id_t root = alea_intersection(sys, u, n3);

    int m = alea_add_material(sys, 1);
    int ci = alea_add_cell(sys, 1, root, m, 1.0, 0);
    ASSERT(ci >= 0);

    char* expr = alea_cell_expr(sys, (size_t)ci, ":", " ", "#");
    ASSERT_NOT_NULL(expr);
    ASSERT_STR_EQ(expr, "(-1:-2) -3");
    free(expr);

    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Difference                                                                 */
/* ------------------------------------------------------------------------- */

TEST(expr_difference) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* S1 - S2 = inside S1 minus inside S2 */
    int s1 = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    int s2 = alea_sphere_surface(sys, 2, 0, 0, 0, 2.0);
    alea_node_id_t n1 = alea_halfspace(sys, s1, -1);
    alea_node_id_t n2 = alea_halfspace(sys, s2, -1);
    alea_node_id_t root = alea_difference(sys, n1, n2);

    int m = alea_add_material(sys, 1);
    int ci = alea_add_cell(sys, 1, root, m, 1.0, 0);
    ASSERT(ci >= 0);

    char* expr = alea_cell_expr(sys, (size_t)ci, ":", " ", "#");
    ASSERT_NOT_NULL(expr);
    /* A - B where B is a primitive: "A #B" (no parens on primitive) */
    ASSERT_STR_EQ(expr, "-1 #-2");
    free(expr);

    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Complement of complex expression                                           */
/* ------------------------------------------------------------------------- */

TEST(expr_complement_complex) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* #(S1 & S2) */
    int s1 = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    int s2 = alea_sphere_surface(sys, 2, 0, 0, 0, 2.0);
    alea_node_id_t n1 = alea_halfspace(sys, s1, -1);
    alea_node_id_t n2 = alea_halfspace(sys, s2, -1);
    alea_node_id_t inter = alea_intersection(sys, n1, n2);
    alea_node_id_t root = alea_complement(sys, inter);

    int m = alea_add_material(sys, 1);
    int ci = alea_add_cell(sys, 1, root, m, 1.0, 0);
    ASSERT(ci >= 0);

    char* expr = alea_cell_expr(sys, (size_t)ci, ":", " ", "#");
    ASSERT_NOT_NULL(expr);
    ASSERT_STR_EQ(expr, "#(-1 -2)");
    free(expr);

    /* OpenMC style */
    char* omc = alea_cell_expr(sys, (size_t)ci, " | ", " ", "~");
    ASSERT_NOT_NULL(omc);
    ASSERT_STR_EQ(omc, "~(-1 -2)");
    free(omc);

    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Complement of single primitive (flips sign)                                */
/* ------------------------------------------------------------------------- */

TEST(expr_complement_primitive) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int si = alea_sphere_surface(sys, 5, 0, 0, 0, 3.0);
    alea_node_id_t inside = alea_halfspace(sys, si, -1);
    alea_node_id_t root = alea_complement(sys, inside);

    int m = alea_add_material(sys, 1);
    int ci = alea_add_cell(sys, 1, root, m, 1.0, 0);
    ASSERT(ci >= 0);

    /* complement of -5 should give 5 */
    char* expr = alea_cell_expr(sys, (size_t)ci, ":", " ", "#");
    ASSERT_NOT_NULL(expr);
    ASSERT_STR_EQ(expr, "5");
    free(expr);

    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Error cases                                                                */
/* ------------------------------------------------------------------------- */

TEST(expr_null_system) {
    char* expr = alea_cell_expr(NULL, 0, ":", " ", "#");
    ASSERT_NULL(expr);
}

TEST(expr_invalid_index) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    char* expr = alea_cell_expr(sys, 999, ":", " ", "#");
    ASSERT_NULL(expr);

    alea_destroy(sys);
}

TEST_MAIN()
