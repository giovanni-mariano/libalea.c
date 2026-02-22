// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_ops.c - CSG tree construction and evaluation edge cases
 */

#include "alea_test.h"
#include "alea.h"
#include "core/alea_system.h"
#include "core/alea_ops.h"
#include "core/alea_eval.h"
#include "core/alea_simplify.h"
#include <string.h>

/* ========================================================================= */
/* Basic tree operations                                                     */
/* ========================================================================= */

TEST(ops_single_primitive) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int si = alea_sphere_surface(sys, 0, 0, 0, 0, 5.0);
    alea_node_id_t node = alea_halfspace(sys, si, -1);
    ASSERT(node != ALEA_NODE_ID_INVALID);

    /* Inside sphere */
    ASSERT_TRUE(alea_contains_point(sys, node, 0, 0, 0));
    /* Outside sphere */
    ASSERT_FALSE(alea_contains_point(sys, node, 10, 0, 0));

    alea_destroy(sys);
}

TEST(ops_deep_nesting) {
    /* Build a deeply nested intersection chain: S1 ∩ S2 ∩ ... ∩ Sn */
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Create 100 planes as an intersection chain */
    int depth = 100;
    int s0 = alea_plane_surface(sys, 0, 0, 0, 1, -100);
    alea_node_id_t current = alea_halfspace(sys, s0, -1);  /* z < 100 */
    for (int i = 1; i < depth; i++) {
        int si = alea_plane_surface(sys, 0, 0, 0, 1, -(100 + i));
        alea_node_id_t next = alea_halfspace(sys, si, -1);
        current = alea_intersection(sys, current, next);
    }
    ASSERT(current != ALEA_NODE_ID_INVALID);

    /* Point at z=50 should be inside all (z < 100, z < 101, ...) */
    ASSERT_TRUE(alea_contains_point(sys, current, 0, 0, 50));

    alea_destroy(sys);
}

TEST(ops_wide_union) {
    /* Union of many spheres */
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int n = 50;
    alea_node_id_t nodes[50];
    for (int i = 0; i < n; i++) {
        int si = alea_sphere_surface(sys, 0, (double)(i * 10), 0, 0, 3.0);
        nodes[i] = alea_halfspace(sys, si, -1);
    }
    alea_node_id_t root = alea_union_n(sys, nodes, n);
    ASSERT(root != ALEA_NODE_ID_INVALID);

    /* Center of first sphere */
    ASSERT_TRUE(alea_contains_point(sys, root, 0, 0, 0));
    /* Center of last sphere */
    ASSERT_TRUE(alea_contains_point(sys, root, 490, 0, 0));
    /* Between spheres (gap) */
    ASSERT_FALSE(alea_contains_point(sys, root, 5, 0, 0));

    alea_destroy(sys);
}

TEST(ops_complement_twice) {
    /* NOT(NOT(A)) = A */
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int si = alea_sphere_surface(sys, 0, 0, 0, 0, 5.0);
    alea_node_id_t sphere = alea_halfspace(sys, si, -1);
    alea_node_id_t not1 = alea_complement(sys, sphere);
    alea_node_id_t not2 = alea_complement(sys, not1);

    /* Double complement should behave same as original */
    ASSERT_TRUE(alea_contains_point(sys, not2, 0, 0, 0));
    ASSERT_FALSE(alea_contains_point(sys, not2, 10, 0, 0));

    alea_destroy(sys);
}

TEST(ops_difference_eval) {
    /* A - B: point in A but not B → inside */
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int si1 = alea_sphere_surface(sys, 0, 0, 0, 0, 5.0);
    int si2 = alea_sphere_surface(sys, 0, 0, 0, 0, 2.0);
    alea_node_id_t big = alea_halfspace(sys, si1, -1);
    alea_node_id_t small = alea_halfspace(sys, si2, -1);
    alea_node_id_t diff = alea_difference(sys, big, small);

    /* Between radii: inside */
    ASSERT_TRUE(alea_contains_point(sys, diff, 3, 0, 0));
    /* Inside small sphere: outside (subtracted) */
    ASSERT_FALSE(alea_contains_point(sys, diff, 0, 0, 0));
    /* Outside big sphere: outside */
    ASSERT_FALSE(alea_contains_point(sys, diff, 10, 0, 0));

    alea_destroy(sys);
}

TEST(ops_difference_eval_neg) {
    /* A - B: point in both A and B → outside */
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int si1 = alea_sphere_surface(sys, 0, 0, 0, 0, 5.0);
    int si2 = alea_sphere_surface(sys, 0, 0, 0, 0, 5.0);
    alea_node_id_t a = alea_halfspace(sys, si1, -1);
    alea_node_id_t b = alea_halfspace(sys, si2, -1);
    alea_node_id_t diff = alea_difference(sys, a, b);

    /* A - A = empty everywhere */
    ASSERT_FALSE(alea_contains_point(sys, diff, 0, 0, 0));
    ASSERT_FALSE(alea_contains_point(sys, diff, 3, 0, 0));

    alea_destroy(sys);
}

TEST(ops_empty_intersection) {
    /* A ∩ NOT(A) should be empty everywhere */
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int si = alea_sphere_surface(sys, 0, 0, 0, 0, 5.0);
    alea_node_id_t sphere = alea_halfspace(sys, si, -1);
    alea_node_id_t not_sphere = alea_complement(sys, sphere);
    alea_node_id_t empty = alea_intersection(sys, sphere, not_sphere);

    ASSERT_FALSE(alea_contains_point(sys, empty, 0, 0, 0));
    ASSERT_FALSE(alea_contains_point(sys, empty, 3, 0, 0));
    ASSERT_FALSE(alea_contains_point(sys, empty, 10, 0, 0));

    alea_destroy(sys);
}

TEST(ops_universe_complement) {
    /* A ∪ NOT(A) should be full everywhere */
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int si = alea_sphere_surface(sys, 0, 0, 0, 0, 5.0);
    alea_node_id_t sphere = alea_halfspace(sys, si, -1);
    alea_node_id_t not_sphere = alea_complement(sys, sphere);
    alea_node_id_t full = alea_union(sys, sphere, not_sphere);

    ASSERT_TRUE(alea_contains_point(sys, full, 0, 0, 0));
    ASSERT_TRUE(alea_contains_point(sys, full, 10, 0, 0));
    ASSERT_TRUE(alea_contains_point(sys, full, -100, 50, 30));

    alea_destroy(sys);
}

TEST(ops_demorgan_union) {
    /* NOT(A ∪ B) should equal NOT(A) ∩ NOT(B) after NNF */
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int si1 = alea_sphere_surface(sys, 0, -3, 0, 0, 3.0);
    int si2 = alea_sphere_surface(sys, 0, 3, 0, 0, 3.0);
    alea_node_id_t a = alea_halfspace(sys, si1, -1);
    alea_node_id_t b = alea_halfspace(sys, si2, -1);

    /* NOT(A ∪ B) */
    alea_node_id_t ab_union = alea_union(sys, a, b);
    alea_node_id_t not_ab = alea_complement(sys, ab_union);

    /* NOT(A) ∩ NOT(B) */
    alea_node_id_t not_a = alea_complement(sys, a);
    alea_node_id_t not_b = alea_complement(sys, b);
    alea_node_id_t demorgan = alea_intersection(sys, not_a, not_b);

    /* Should agree everywhere */
    double test_pts[][3] = {{0,0,0}, {-3,0,0}, {3,0,0}, {10,0,0}, {0,10,0}};
    for (int i = 0; i < 5; i++) {
        bool r1 = alea_contains_point(sys, not_ab, test_pts[i][0], test_pts[i][1], test_pts[i][2]);
        bool r2 = alea_contains_point(sys, demorgan, test_pts[i][0], test_pts[i][1], test_pts[i][2]);
        ASSERT_EQ(r1, r2);
    }

    alea_destroy(sys);
}

TEST_MAIN()
