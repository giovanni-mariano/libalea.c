// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_simplify.c - Unit tests for CSG tree simplification
 */

#include "alea_test.h"
#include "alea.h"
#include "core/alea_system.h"
#include "core/alea_simplify.h"
#include "core/alea_ops.h"
#include "core/alea_eval.h"
#include "primitives/bbox.h"

/* ------------------------------------------------------------------------- */
/* Test Utilities                                                             */
/* ------------------------------------------------------------------------- */

/* Compare evaluation at multiple points */
static int trees_equivalent(alea_system_t* sys, alea_node_id_t a, alea_node_id_t b) {
    double test_points[][3] = {
        {0, 0, 0},
        {1, 0, 0}, {-1, 0, 0},
        {0, 1, 0}, {0, -1, 0},
        {0, 0, 1}, {0, 0, -1},
        {0.5, 0.5, 0.5}, {-0.5, -0.5, -0.5},
        {2, 2, 2}, {-2, -2, -2},
        {0.1, 0.2, 0.3}, {-0.1, -0.2, -0.3},
    };
    int n_points = sizeof(test_points) / sizeof(test_points[0]);

    for (int i = 0; i < n_points; i++) {
        double x = test_points[i][0];
        double y = test_points[i][1];
        double z = test_points[i][2];

        double val_a = alea_evaluate_point(sys, a, x, y, z);
        double val_b = alea_evaluate_point(sys, b, x, y, z);

        int inside_a = (val_a <= 0.0);
        int inside_b = (val_b <= 0.0);

        if (inside_a != inside_b) {
            return 0;
        }
    }
    return 1;
}

/* Create a sphere surface with both sense nodes */
static void create_sphere(alea_system_t* sys, int id, double cx, double cy, double cz, double r,
                          alea_node_id_t* pos, alea_node_id_t* neg) {
    alea_primitive_data_t data;
    data.sphere.center_x = cx;
    data.sphere.center_y = cy;
    data.sphere.center_z = cz;
    data.sphere.radius = r;

    int8_t inverted = 0;
    alea_primitive_id_t prim_id = alea_get_or_create_primitive(sys, ALEA_PRIMITIVE_SPHERE, &data, &inverted);

    *pos = alea_add_primitive_node(sys, prim_id, +1, inverted, id);
    *neg = alea_add_primitive_node(sys, prim_id, -1, inverted, id);

    alea_surface_entry_t* surf = alea_vec_push_uninit(&sys->surfaces, alea_surface_entry_t);
    if (surf) {
        memset(surf, 0, sizeof(*surf));
        surf->mc_surface_id = id;
        surf->primitive_id = prim_id;
        surf->pos_node = *pos;
        surf->neg_node = *neg;
    }
}

/* ------------------------------------------------------------------------- */
/* NNF (Negation Normal Form) Tests                                          */
/* ------------------------------------------------------------------------- */

TEST(nnf_primitive_unchanged) {
    alea_system_t* sys = alea_create();
    alea_node_id_t pos, neg;
    create_sphere(sys, 1, 0, 0, 0, 5.0, &pos, &neg);

    alea_node_id_t result = alea_tree_to_nnf(sys, pos);
    ASSERT_EQ(result, pos);

    alea_destroy(sys);
}

TEST(nnf_complement_to_opposite_sense) {
    alea_system_t* sys = alea_create();
    alea_node_id_t pos, neg;
    create_sphere(sys, 1, 0, 0, 0, 5.0, &pos, &neg);

    /* NOT(+S) should become -S */
    alea_node_id_t comp = alea_create_complement(sys, pos);
    alea_node_id_t result = alea_tree_to_nnf(sys, comp);

    ASSERT_EQ(result, neg);
    ASSERT_TRUE(trees_equivalent(sys, comp, result));

    alea_destroy(sys);
}

TEST(nnf_double_complement) {
    alea_system_t* sys = alea_create();
    alea_node_id_t pos, neg;
    create_sphere(sys, 1, 0, 0, 0, 5.0, &pos, &neg);

    /* NOT(NOT(+S)) should become +S */
    alea_node_id_t comp1 = alea_create_complement(sys, pos);
    alea_node_id_t comp2 = alea_create_complement(sys, comp1);
    alea_node_id_t result = alea_tree_to_nnf(sys, comp2);

    ASSERT_EQ(result, pos);

    alea_destroy(sys);
}

TEST(nnf_demorgan_intersection) {
    alea_system_t* sys = alea_create();
    alea_node_id_t pos1, neg1, pos2, neg2;
    create_sphere(sys, 1, 0, 0, 0, 5.0, &pos1, &neg1);
    create_sphere(sys, 2, 3, 0, 0, 2.0, &pos2, &neg2);

    /* NOT(A AND B) should become NOT(A) OR NOT(B) */
    alea_node_id_t inter = alea_create_intersection(sys, neg1, neg2);
    alea_node_id_t comp = alea_create_complement(sys, inter);
    alea_node_id_t result = alea_tree_to_nnf(sys, comp);

    ASSERT(result < alea_vec_count(&sys->nodes));
    alea_operation_t op = ALEA_GET_OPERATION(&sys->nodes.data[result]);
    ASSERT_EQ(op, ALEA_OP_UNION);
    ASSERT_TRUE(trees_equivalent(sys, comp, result));

    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Flatten and Optimize Tests                                                 */
/* ------------------------------------------------------------------------- */

TEST(flatten_simple_intersection) {
    alea_system_t* sys = alea_create();
    alea_node_id_t pos1, neg1, pos2, neg2;
    create_sphere(sys, 1, 0, 0, 0, 5.0, &pos1, &neg1);
    create_sphere(sys, 2, 3, 0, 0, 2.0, &pos2, &neg2);

    alea_node_id_t original = alea_create_intersection(sys, neg1, neg2);

    alea_simplify_stats_t stats = {0};
    alea_node_id_t result = alea_tree_simplify(sys, original, &stats);

    ASSERT_NE(result, ALEA_NODE_ID_INVALID);
    ASSERT_TRUE(trees_equivalent(sys, original, result));

    alea_destroy(sys);
}

TEST(flatten_nested_intersection) {
    alea_system_t* sys = alea_create();
    alea_node_id_t pos1, neg1, pos2, neg2, pos3, neg3;
    create_sphere(sys, 1, 0, 0, 0, 10.0, &pos1, &neg1);
    create_sphere(sys, 2, 2, 0, 0, 8.0, &pos2, &neg2);
    create_sphere(sys, 3, 1, 1, 0, 6.0, &pos3, &neg3);

    /* (A AND B) AND C */
    alea_node_id_t ab = alea_create_intersection(sys, neg1, neg2);
    alea_node_id_t original = alea_create_intersection(sys, ab, neg3);

    alea_simplify_stats_t stats = {0};
    alea_node_id_t result = alea_tree_simplify(sys, original, &stats);

    ASSERT_NE(result, ALEA_NODE_ID_INVALID);
    ASSERT_TRUE(trees_equivalent(sys, original, result));

    alea_destroy(sys);
}

TEST(flatten_contradiction) {
    alea_system_t* sys = alea_create();
    alea_node_id_t pos, neg;
    create_sphere(sys, 1, 0, 0, 0, 5.0, &pos, &neg);

    /* +S AND -S = empty (inside AND outside) */
    alea_node_id_t original = alea_create_intersection(sys, pos, neg);

    alea_simplify_stats_t stats = {0};
    alea_node_id_t result = alea_tree_simplify(sys, original, &stats);

    ASSERT_EQ(result, ALEA_NODE_ID_INVALID);
    ASSERT(stats.contradictions_found > 0);

    alea_destroy(sys);
}

TEST(flatten_nested_contradiction) {
    alea_system_t* sys = alea_create();
    alea_node_id_t pos1, neg1, pos2, neg2;
    create_sphere(sys, 1, 0, 0, 0, 5.0, &pos1, &neg1);
    create_sphere(sys, 2, 3, 0, 0, 2.0, &pos2, &neg2);

    /* (inside_s1 AND inside_s2) AND outside_s1 = empty */
    alea_node_id_t ab = alea_create_intersection(sys, neg1, neg2);
    alea_node_id_t original = alea_create_intersection(sys, ab, pos1);

    alea_simplify_stats_t stats = {0};
    alea_node_id_t result = alea_tree_simplify(sys, original, &stats);

    ASSERT_EQ(result, ALEA_NODE_ID_INVALID);
    ASSERT(stats.contradictions_found > 0);

    alea_destroy(sys);
}

TEST(flatten_idempotent) {
    alea_system_t* sys = alea_create();
    alea_node_id_t pos, neg;
    create_sphere(sys, 1, 0, 0, 0, 5.0, &pos, &neg);

    /* A AND A should become just A */
    alea_node_id_t original = alea_create_intersection(sys, neg, neg);

    alea_simplify_stats_t stats = {0};
    alea_node_id_t result = alea_tree_simplify(sys, original, &stats);

    ASSERT_EQ(result, neg);
    ASSERT(stats.idempotent_reductions > 0);

    alea_destroy(sys);
}

TEST(flatten_mixed_union_intersection) {
    alea_system_t* sys = alea_create();
    alea_node_id_t pos1, neg1, pos2, neg2, pos3, neg3;
    create_sphere(sys, 1, 0, 0, 0, 10.0, &pos1, &neg1);
    create_sphere(sys, 2, 5, 0, 0, 3.0, &pos2, &neg2);
    create_sphere(sys, 3, -5, 0, 0, 3.0, &pos3, &neg3);

    /* (B OR C) AND A */
    alea_node_id_t bc = alea_create_union(sys, neg2, neg3);
    alea_node_id_t original = alea_create_intersection(sys, bc, neg1);

    alea_simplify_stats_t stats = {0};
    alea_node_id_t result = alea_tree_simplify(sys, original, &stats);

    ASSERT_NE(result, ALEA_NODE_ID_INVALID);
    ASSERT_TRUE(trees_equivalent(sys, original, result));

    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Statistical Simplification Tests                                           */
/* ------------------------------------------------------------------------- */

TEST(statistical_no_false_removal) {
    alea_system_t* sys = alea_create();
    alea_node_id_t pos1, neg1, pos2, neg2;
    create_sphere(sys, 1, 0, 0, 0, 5.0, &pos1, &neg1);
    create_sphere(sys, 2, 3, 0, 0, 4.0, &pos2, &neg2);

    /* Inside sphere 1 AND inside sphere 2 (overlapping region) */
    alea_node_id_t original = alea_create_intersection(sys, neg1, neg2);

    alea_simplify_stats_t stats = {0};
    alea_node_id_t result = alea_tree_simplify(sys, original, &stats);

    ASSERT_NE(result, ALEA_NODE_ID_INVALID);
    ASSERT_TRUE(trees_equivalent(sys, original, result));

    alea_destroy(sys);
}

/* Create a plane surface with both sense nodes */
static void create_plane(alea_system_t* sys, int id, double a, double b, double c, double d,
                          alea_node_id_t* pos, alea_node_id_t* neg) {
    alea_primitive_data_t data;
    data.plane.a = a;
    data.plane.b = b;
    data.plane.c = c;
    data.plane.d = d;

    int8_t inverted = 0;
    alea_primitive_id_t prim_id = alea_get_or_create_primitive(sys, ALEA_PRIMITIVE_PLANE, &data, &inverted);

    *pos = alea_add_primitive_node(sys, prim_id, +1, inverted, id);
    *neg = alea_add_primitive_node(sys, prim_id, -1, inverted, id);

    alea_surface_entry_t* surf = alea_vec_push_uninit(&sys->surfaces, alea_surface_entry_t);
    if (surf) {
        memset(surf, 0, sizeof(*surf));
        surf->mc_surface_id = id;
        surf->primitive_id = prim_id;
        surf->pos_node = *pos;
        surf->neg_node = *neg;
    }
}

/* Create a Z-axis cylinder surface with both sense nodes */
static void create_cylinder_z(alea_system_t* sys, int id, double cx, double cy, double r,
                               alea_node_id_t* pos, alea_node_id_t* neg) {
    alea_primitive_data_t data;
    memset(&data, 0, sizeof(data));
    data.cyl_z.center_x = cx;
    data.cyl_z.center_y = cy;
    data.cyl_z.radius = r;

    int8_t inverted = 0;
    alea_primitive_id_t prim_id = alea_get_or_create_primitive(sys, ALEA_PRIMITIVE_CYLINDER_Z, &data, &inverted);

    *pos = alea_add_primitive_node(sys, prim_id, +1, inverted, id);
    *neg = alea_add_primitive_node(sys, prim_id, -1, inverted, id);

    alea_surface_entry_t* surf = alea_vec_push_uninit(&sys->surfaces, alea_surface_entry_t);
    if (surf) {
        memset(surf, 0, sizeof(*surf));
        surf->mc_surface_id = id;
        surf->primitive_id = prim_id;
        surf->pos_node = *pos;
        surf->neg_node = *neg;
    }
}

/* ------------------------------------------------------------------------- */
/* Semantic Deduplication Tests                                                */
/* ------------------------------------------------------------------------- */

TEST(semantic_dedup_same_primitive_different_nodes) {
    alea_system_t* sys = alea_create();
    alea_node_id_t pos1, neg1;
    create_sphere(sys, 1, 0, 0, 0, 5.0, &pos1, &neg1);

    /* Create a second node referencing the same primitive with same sense */
    const alea_node_t* n = &sys->nodes.data[neg1];
    alea_node_id_t neg1_dup = alea_add_primitive_node(sys, n->primitive.primitive_id,
                                                     n->primitive.sense,
                                                     n->primitive.inverted, 1);
    /* neg1 and neg1_dup have different node IDs but same (primitive_id, sense) */
    ASSERT_NE(neg1, neg1_dup);

    /* Build intersection: neg1 AND neg1_dup */
    alea_node_id_t original = alea_create_intersection(sys, neg1, neg1_dup);

    alea_simplify_stats_t stats = {0};
    alea_node_id_t result = alea_tree_simplify(sys, original, &stats);

    /* Should reduce to a single term */
    ASSERT_NE(result, ALEA_NODE_ID_INVALID);
    ASSERT(stats.idempotent_reductions > 0);

    /* Result should be a primitive, not an intersection */
    ASSERT_EQ(ALEA_GET_OPERATION(&sys->nodes.data[result]), ALEA_OP_PRIMITIVE);

    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Pairwise Containment Tests                                                  */
/* ------------------------------------------------------------------------- */

TEST(pairwise_parallel_planes_neg_sense) {
    alea_system_t* sys = alea_create();

    /* Two parallel planes along Z: z+5=0 and z+10=0 */
    /* sense < 0 means z < -d, so z < -5 and z < -10 */
    /* z < -10 is tighter, z < -5 is redundant */
    alea_node_id_t pos1, neg1, pos2, neg2;
    create_plane(sys, 1, 0, 0, 1, 5.0, &pos1, &neg1);   /* z + 5 = 0 */
    create_plane(sys, 2, 0, 0, 1, 10.0, &pos2, &neg2);  /* z + 10 = 0 */

    /* Also need a bounding sphere so the intersection has a finite bbox */
    alea_node_id_t spos, sneg;
    create_sphere(sys, 3, 0, 0, -20, 50.0, &spos, &sneg);

    /* Intersection: inside sphere AND below plane1 AND below plane2 */
    alea_node_id_t ab = alea_create_intersection(sys, neg1, neg2);
    alea_node_id_t original = alea_create_intersection(sys, sneg, ab);

    alea_simplify_stats_t stats = {0};
    alea_node_id_t result = alea_tree_simplify(sys, original, &stats);

    ASSERT_NE(result, ALEA_NODE_ID_INVALID);
    /* Should have removed the looser plane */
    ASSERT(stats.absorption_reductions > 0);
    ASSERT_TRUE(trees_equivalent(sys, original, result));

    alea_destroy(sys);
}

TEST(pairwise_parallel_planes_pos_sense) {
    alea_system_t* sys = alea_create();

    /* Two parallel planes along Z: z-5=0 and z-10=0 */
    /* sense > 0 means z > d_norm, so z > 5 and z > 10 */
    /* z > 10 is tighter, z > 5 is redundant */
    alea_node_id_t pos1, neg1, pos2, neg2;
    create_plane(sys, 1, 0, 0, 1, -5.0, &pos1, &neg1);   /* z - 5 = 0 */
    create_plane(sys, 2, 0, 0, 1, -10.0, &pos2, &neg2);  /* z - 10 = 0 */

    alea_node_id_t spos, sneg;
    create_sphere(sys, 3, 0, 0, 20, 50.0, &spos, &sneg);

    /* Intersection: inside sphere AND above plane1 AND above plane2 */
    alea_node_id_t ab = alea_create_intersection(sys, pos1, pos2);
    alea_node_id_t original = alea_create_intersection(sys, sneg, ab);

    alea_simplify_stats_t stats = {0};
    alea_node_id_t result = alea_tree_simplify(sys, original, &stats);

    ASSERT_NE(result, ALEA_NODE_ID_INVALID);
    ASSERT(stats.absorption_reductions > 0);
    ASSERT_TRUE(trees_equivalent(sys, original, result));

    alea_destroy(sys);
}

TEST(pairwise_concentric_cylinders_inside) {
    alea_system_t* sys = alea_create();

    /* Two concentric Z-axis cylinders: r=3 and r=5 */
    /* sense < 0 = inside: inside(r=3) ⊂ inside(r=5), so r=5 is redundant */
    alea_node_id_t pos1, neg1, pos2, neg2;
    create_cylinder_z(sys, 1, 0, 0, 3.0, &pos1, &neg1);
    create_cylinder_z(sys, 2, 0, 0, 5.0, &pos2, &neg2);

    alea_node_id_t original = alea_create_intersection(sys, neg1, neg2);

    alea_simplify_stats_t stats = {0};
    alea_node_id_t result = alea_tree_simplify(sys, original, &stats);

    ASSERT_NE(result, ALEA_NODE_ID_INVALID);
    ASSERT(stats.absorption_reductions > 0);

    /* Should reduce to just inside(r=3) */
    ASSERT_EQ(ALEA_GET_OPERATION(&sys->nodes.data[result]), ALEA_OP_PRIMITIVE);
    ASSERT_TRUE(trees_equivalent(sys, original, result));

    alea_destroy(sys);
}

TEST(pairwise_concentric_cylinders_outside) {
    alea_system_t* sys = alea_create();

    /* Two concentric Z-axis cylinders: r=3 and r=5 */
    /* sense > 0 = outside: outside(r=5) ⊂ outside(r=3), so r=3 is redundant */
    alea_node_id_t pos1, neg1, pos2, neg2;
    create_cylinder_z(sys, 1, 0, 0, 3.0, &pos1, &neg1);
    create_cylinder_z(sys, 2, 0, 0, 5.0, &pos2, &neg2);

    /* Need bounding sphere for bbox-based checks */
    alea_node_id_t spos, sneg;
    create_sphere(sys, 3, 0, 0, 0, 20.0, &spos, &sneg);

    alea_node_id_t ab = alea_create_intersection(sys, pos1, pos2);
    alea_node_id_t original = alea_create_intersection(sys, sneg, ab);

    alea_simplify_stats_t stats = {0};
    alea_node_id_t result = alea_tree_simplify(sys, original, &stats);

    ASSERT_NE(result, ALEA_NODE_ID_INVALID);
    ASSERT(stats.absorption_reductions > 0);
    ASSERT_TRUE(trees_equivalent(sys, original, result));

    alea_destroy(sys);
}

TEST(pairwise_concentric_spheres_inside) {
    alea_system_t* sys = alea_create();

    /* Two concentric spheres: r=3 and r=5 */
    /* sense < 0 = inside: inside(r=3) ⊂ inside(r=5), so r=5 is redundant */
    alea_node_id_t pos1, neg1, pos2, neg2;
    create_sphere(sys, 1, 0, 0, 0, 3.0, &pos1, &neg1);
    create_sphere(sys, 2, 0, 0, 0, 5.0, &pos2, &neg2);

    alea_node_id_t original = alea_create_intersection(sys, neg1, neg2);

    alea_simplify_stats_t stats = {0};
    alea_node_id_t result = alea_tree_simplify(sys, original, &stats);

    ASSERT_NE(result, ALEA_NODE_ID_INVALID);
    ASSERT(stats.absorption_reductions > 0);

    /* Should reduce to just inside(r=3) */
    ASSERT_EQ(ALEA_GET_OPERATION(&sys->nodes.data[result]), ALEA_OP_PRIMITIVE);
    ASSERT_TRUE(trees_equivalent(sys, original, result));

    alea_destroy(sys);
}

TEST(pairwise_eccentric_spheres_inside) {
    alea_system_t* sys = alea_create();

    /* Sphere A at origin r=2, sphere B at (1,0,0) r=5 */
    /* dist=1, dist+ra=3 <= rb=5, so A ⊂ B, and inside(B) is redundant */
    alea_node_id_t pos1, neg1, pos2, neg2;
    create_sphere(sys, 1, 0, 0, 0, 2.0, &pos1, &neg1);
    create_sphere(sys, 2, 1, 0, 0, 5.0, &pos2, &neg2);

    alea_node_id_t original = alea_create_intersection(sys, neg1, neg2);

    alea_simplify_stats_t stats = {0};
    alea_node_id_t result = alea_tree_simplify(sys, original, &stats);

    ASSERT_NE(result, ALEA_NODE_ID_INVALID);
    ASSERT(stats.absorption_reductions > 0);

    ASSERT_EQ(ALEA_GET_OPERATION(&sys->nodes.data[result]), ALEA_OP_PRIMITIVE);
    ASSERT_TRUE(trees_equivalent(sys, original, result));

    alea_destroy(sys);
}

TEST(pairwise_no_false_removal) {
    alea_system_t* sys = alea_create();

    /* Non-concentric, overlapping cylinders - no subsumption possible */
    alea_node_id_t pos1, neg1, pos2, neg2;
    create_cylinder_z(sys, 1, 0, 0, 3.0, &pos1, &neg1);
    create_cylinder_z(sys, 2, 5, 0, 3.0, &pos2, &neg2);

    alea_node_id_t original = alea_create_intersection(sys, neg1, neg2);

    alea_simplify_stats_t stats = {0};
    alea_node_id_t result = alea_tree_simplify(sys, original, &stats);

    /* Should NOT reduce - both terms are needed */
    ASSERT_NE(result, ALEA_NODE_ID_INVALID);
    ASSERT_TRUE(trees_equivalent(sys, original, result));

    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Union Optimization Tests                                                   */
/* ------------------------------------------------------------------------- */

TEST(union_algebraic_absorption) {
    /* A ∪ (A∩B) → A  (A has fewer constraints, absorbs A∩B) */
    /* Use planes so the intersection A∩B doesn't simplify away A or B */
    alea_system_t* sys = alea_create();

    alea_node_id_t pA_pos, pA_neg, pB_pos, pB_neg;
    create_plane(sys, 1, 0, 0, 1, 0.0, &pA_pos, &pA_neg);  /* z > 0 */
    create_plane(sys, 2, 1, 0, 0, 0.0, &pB_pos, &pB_neg);  /* x > 0 */

    /* Branch 1: +P1 (z > 0) */
    /* Branch 2: +P1 ∩ +P2 (z > 0 AND x > 0) */
    alea_node_id_t branch2 = alea_create_intersection(sys, pA_pos, pB_pos);
    alea_node_id_t original = alea_create_union(sys, pA_pos, branch2);

    alea_simplify_stats_t stats = {0};
    alea_node_id_t result = alea_tree_simplify(sys, original, &stats);

    ASSERT_NE(result, ALEA_NODE_ID_INVALID);
    /* Result should be equivalent to just z > 0 */
    ASSERT_TRUE(trees_equivalent(sys, pA_pos, result));
    ASSERT_TRUE(stats.union_branches_absorbed > 0);

    alea_destroy(sys);
}

TEST(union_no_false_absorption) {
    /* (A∩B) ∪ (A∩C) — neither absorbed (different non-common literals) */
    alea_system_t* sys = alea_create();

    alea_node_id_t sA_pos, sA_neg, sB_pos, sB_neg, sC_pos, sC_neg;
    create_sphere(sys, 1, 0, 0, 0, 10.0, &sA_pos, &sA_neg);
    create_sphere(sys, 2, 5, 0, 0, 3.0, &sB_pos, &sB_neg);
    create_sphere(sys, 3, -5, 0, 0, 3.0, &sC_pos, &sC_neg);

    alea_node_id_t branch1 = alea_create_intersection(sys, sA_neg, sB_neg);
    alea_node_id_t branch2 = alea_create_intersection(sys, sA_neg, sC_neg);
    alea_node_id_t original = alea_create_union(sys, branch1, branch2);

    alea_simplify_stats_t stats = {0};
    alea_node_id_t result = alea_tree_simplify(sys, original, &stats);

    ASSERT_NE(result, ALEA_NODE_ID_INVALID);
    /* Should preserve both branches — neither absorbs the other */
    ASSERT_TRUE(trees_equivalent(sys, original, result));

    alea_destroy(sys);
}

TEST(union_common_factor_extraction) {
    /* (A∩B∩C) ∪ (A∩B∩D) → A∩B∩(C∪D) */
    /* Use all planes so geometric pruning doesn't reduce intersection terms */
    alea_system_t* sys = alea_create();

    alea_node_id_t pa_pos, pa_neg, pb_pos, pb_neg, pc_pos, pc_neg, pd_pos, pd_neg;
    create_plane(sys, 1, 0, 0, 1, 0.0, &pa_pos, &pa_neg);   /* A: z > 0 */
    create_plane(sys, 2, 1, 0, 0, 0.0, &pb_pos, &pb_neg);   /* B: x > 0 */
    create_plane(sys, 3, 0, 1, 0, -1.0, &pc_pos, &pc_neg);  /* C: y > 1 */
    create_plane(sys, 4, 0, 1, 0, 1.0, &pd_pos, &pd_neg);   /* D: y < -1 (neg sense) */

    /* Branch 1: z>0 ∩ x>0 ∩ y>1 */
    alea_node_id_t ab1 = alea_create_intersection(sys, pa_pos, pb_pos);
    alea_node_id_t branch1 = alea_create_intersection(sys, ab1, pc_pos);

    /* Branch 2: z>0 ∩ x>0 ∩ y<-1 */
    alea_node_id_t ab2 = alea_create_intersection(sys, pa_pos, pb_pos);
    alea_node_id_t branch2 = alea_create_intersection(sys, ab2, pd_neg);

    alea_node_id_t original = alea_create_union(sys, branch1, branch2);

    alea_simplify_stats_t stats = {0};
    alea_node_id_t result = alea_tree_simplify(sys, original, &stats);

    ASSERT_NE(result, ALEA_NODE_ID_INVALID);
    /* Result should be equivalent to original */
    ASSERT_TRUE(trees_equivalent(sys, original, result));
    /* Should have factored */
    ASSERT_TRUE(stats.union_common_factors > 0);

    alea_destroy(sys);
}

TEST(union_geometric_subsumption) {
    /* inside(small_sphere) ∪ inside(big_sphere) → inside(big_sphere)
     * (small is geometrically subsumed by big) */
    alea_system_t* sys = alea_create();

    alea_node_id_t sA_pos, sA_neg, sB_pos, sB_neg;
    create_sphere(sys, 1, 0, 0, 0, 3.0, &sA_pos, &sA_neg);   /* small */
    create_sphere(sys, 2, 0, 0, 0, 10.0, &sB_pos, &sB_neg);  /* big */

    alea_node_id_t original = alea_create_union(sys, sA_neg, sB_neg);

    alea_simplify_stats_t stats = {0};
    alea_node_id_t result = alea_tree_simplify(sys, original, &stats);

    ASSERT_NE(result, ALEA_NODE_ID_INVALID);
    /* Result should be equivalent to just -S2 (big sphere) */
    ASSERT_TRUE(trees_equivalent(sys, sB_neg, result));
    /* Either algebraic absorption or geometric subsumption should fire */
    ASSERT_TRUE(stats.union_branches_absorbed + stats.union_branches_subsumed > 0);

    alea_destroy(sys);
}

TEST(union_all_absorbed_except_one) {
    /* A ∪ (A∩B) ∪ (A∩B∩C) → A */
    /* Use perpendicular planes so intersections don't simplify */
    alea_system_t* sys = alea_create();

    alea_node_id_t pA_pos, pA_neg, pB_pos, pB_neg, pC_pos, pC_neg;
    create_plane(sys, 1, 0, 0, 1, 0.0, &pA_pos, &pA_neg);  /* z > 0 */
    create_plane(sys, 2, 1, 0, 0, 0.0, &pB_pos, &pB_neg);  /* x > 0 */
    create_plane(sys, 3, 0, 1, 0, 0.0, &pC_pos, &pC_neg);  /* y > 0 */

    alea_node_id_t ab = alea_create_intersection(sys, pA_pos, pB_pos);
    alea_node_id_t abc = alea_create_intersection(sys, ab, pC_pos);
    alea_node_id_t u1 = alea_create_union(sys, pA_pos, ab);
    alea_node_id_t original = alea_create_union(sys, u1, abc);

    alea_simplify_stats_t stats = {0};
    alea_node_id_t result = alea_tree_simplify(sys, original, &stats);

    ASSERT_NE(result, ALEA_NODE_ID_INVALID);
    /* Result should be equivalent to just z > 0 */
    ASSERT_TRUE(trees_equivalent(sys, pA_pos, result));
    ASSERT_TRUE(stats.union_branches_absorbed >= 2);

    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Cell Splitting Tests                                                       */
/* ------------------------------------------------------------------------- */

TEST(split_union_cell) {
    /* Cell with A∪B∪C → 3 cells, original removed */
    alea_system_t* sys = alea_create();

    alea_node_id_t sA_pos, sA_neg, sB_pos, sB_neg, sC_pos, sC_neg;
    create_sphere(sys, 1, -5, 0, 0, 3.0, &sA_pos, &sA_neg);
    create_sphere(sys, 2, 0, 0, 0, 3.0, &sB_pos, &sB_neg);
    create_sphere(sys, 3, 5, 0, 0, 3.0, &sC_pos, &sC_neg);

    alea_node_id_t u1 = alea_create_union(sys, sA_neg, sB_neg);
    alea_node_id_t root = alea_create_union(sys, u1, sC_neg);

    int m1 = alea_add_material(sys, 1);

    alea_add_cell(sys, 100, root, m1, -1.0, 0);
    ASSERT_EQ((int)alea_vec_count(&sys->cells), 1);

    int created = alea_split_union_cells(sys);
    ASSERT_EQ(created, 3);
    ASSERT_EQ((int)alea_vec_count(&sys->cells), 3);

    /* All new cells should have material 1 (MCNP ID) */
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        ASSERT_EQ(sys->cells.data[i].material_id, 1);
    }

    alea_destroy(sys);
}

TEST(split_no_union_cell) {
    /* Pure intersection cell unchanged by split */
    alea_system_t* sys = alea_create();

    alea_node_id_t sA_pos, sA_neg, sB_pos, sB_neg;
    create_sphere(sys, 1, 0, 0, 0, 10.0, &sA_pos, &sA_neg);
    create_sphere(sys, 2, 0, 0, 0, 5.0, &sB_pos, &sB_neg);

    alea_node_id_t root = alea_create_intersection(sys, sA_neg, sB_neg);

    int m2 = alea_add_material(sys, 2);

    alea_add_cell(sys, 100, root, m2, -2.0, 0);
    ASSERT_EQ((int)alea_vec_count(&sys->cells), 1);

    int created = alea_split_union_cells(sys);
    ASSERT_EQ(created, 0);
    ASSERT_EQ((int)alea_vec_count(&sys->cells), 1);
    ASSERT_EQ(sys->cells.data[0].material_id, 2);

    alea_destroy(sys);
}

/* ------------------------------------------------------------------------- */
/* Inverted Flag Regression Tests                                              */
/* ------------------------------------------------------------------------- */

TEST(pairwise_inverted_plane_containment) {
    /* Regression: pairwise containment must use effective_sense (accounting for
     * inverted flag).  Create planes with negative leading coefficient so
     * canonicalization sets inverted=1.  The simplifier must still correctly
     * identify containment. */
    alea_system_t* sys = alea_create();

    /* Plane A: -z + 5 = 0  →  canonical z - 5 = 0, inverted=1
     * Plane B: -z + 10 = 0 →  canonical z - 10 = 0, inverted=1
     *
     * With inverted=1 and raw sense=-1:
     *   effective_sense = -(-1) = +1, i.e. z > 5 and z > 10
     *   z > 10 is tighter, z > 5 is redundant.
     */
    alea_primitive_data_t dataA = { .plane = { 0, 0, -1, 5.0 } };
    alea_primitive_data_t dataB = { .plane = { 0, 0, -1, 10.0 } };
    int8_t invA = 0, invB = 0;
    alea_primitive_id_t pidA = alea_get_or_create_primitive(sys, ALEA_PRIMITIVE_PLANE, &dataA, &invA);
    alea_primitive_id_t pidB = alea_get_or_create_primitive(sys, ALEA_PRIMITIVE_PLANE, &dataB, &invB);

    /* Both should be canonicalized with inverted=1 */
    ASSERT_EQ(invA, 1);
    ASSERT_EQ(invB, 1);

    alea_node_id_t nA = alea_add_primitive_node(sys, pidA, -1, invA, 1);
    alea_node_id_t nB = alea_add_primitive_node(sys, pidB, -1, invB, 2);

    /* Add a bounding sphere for finite bbox */
    alea_node_id_t spos, sneg;
    create_sphere(sys, 3, 0, 0, 20, 50.0, &spos, &sneg);

    alea_node_id_t ab = alea_create_intersection(sys, nA, nB);
    alea_node_id_t original = alea_create_intersection(sys, sneg, ab);

    alea_simplify_stats_t stats = {0};
    alea_node_id_t result = alea_tree_simplify(sys, original, &stats);

    ASSERT_NE(result, ALEA_NODE_ID_INVALID);
    /* Should remove the looser plane (z>5 is redundant given z>10) */
    ASSERT(stats.absorption_reductions > 0);
    ASSERT_TRUE(trees_equivalent(sys, original, result));

    alea_destroy(sys);
}

TEST(geometric_redundancy_inverted_flag) {
    /* Regression: check_geometric_redundancy must use effective_sense.
     * Create a cell with an inverted plane that is geometrically redundant
     * given the intersection bbox. Verify it is removed without changing
     * the cell's geometry. */
    alea_system_t* sys = alea_create();

    /* Sphere at origin r=5 (interior) — provides finite bbox */
    alea_node_id_t spos, sneg;
    create_sphere(sys, 1, 0, 0, 0, 5.0, &spos, &sneg);

    /* Plane: -z + 100 = 0  →  canonical z - 100 = 0, inverted=1
     * With inverted=1, sense=-1: effective_sense=+1, i.e. z > 100
     * But the sphere bbox is [-5,5] on Z, so z > 100 contradicts the bbox.
     *
     * With correct effective sense handling, this should be detected as
     * a contradiction (cell is empty). With wrong raw sense, the code would
     * see "z < 100" (which is always true in bbox) and wrongly remove it
     * as redundant, keeping a non-empty cell. */
    alea_primitive_data_t dataP = { .plane = { 0, 0, -1, 100.0 } };
    int8_t invP = 0;
    alea_primitive_id_t pidP = alea_get_or_create_primitive(sys, ALEA_PRIMITIVE_PLANE, &dataP, &invP);
    ASSERT_EQ(invP, 1);
    alea_node_id_t nP = alea_add_primitive_node(sys, pidP, -1, invP, 2);

    /* Intersection: inside sphere AND z > 100 → should be empty */
    alea_node_id_t original = alea_create_intersection(sys, sneg, nP);

    alea_simplify_stats_t stats = {0};
    alea_node_id_t result = alea_tree_simplify(sys, original, &stats);

    /* The cell should be detected as empty (contradiction) */
    ASSERT_EQ(result, ALEA_NODE_ID_INVALID);

    alea_destroy(sys);
}

TEST(proof_simplify_dry_run_and_transactional_apply) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    alea_node_id_t small_pos, small_neg, big_pos, big_neg;
    create_sphere(sys, 101, 0, 0, 0, 2.0, &small_pos, &small_neg);
    create_sphere(sys, 102, 0, 0, 0, 5.0, &big_pos, &big_neg);
    alea_node_id_t root = alea_create_intersection(sys, small_neg, big_neg);
    ASSERT(alea_add_cell(sys, 1, root, ALEA_MATERIAL_VOID, 0.0, 9) >= 0);
    size_t node_count = alea_vec_count(&sys->nodes);

    alea_cell_simplify_proof_options_t options;
    alea_cell_simplify_proof_options_init(&options);
    options.max_depth = 5;
    alea_cell_simplify_proof_result_t dry;
    ASSERT_EQ(alea_cell_simplify_proven(sys, 0, &options, &dry), 0);
    ASSERT(dry.changed);
    ASSERT(!dry.applied);
    ASSERT(dry.complete);
    ASSERT_EQ(dry.nodes_before, 3);
    ASSERT_EQ(dry.nodes_after, 1);
    ASSERT_EQ(alea_vec_count(&sys->nodes), node_count);
    ASSERT_EQ(sys->cells.data[0].root_node_id, root);

    options.apply = true;
    alea_cell_simplify_proof_result_t applied;
    ASSERT_EQ(alea_cell_simplify_proven(sys, 0, &options, &applied), 0);
    ASSERT(applied.changed);
    ASSERT(applied.applied);
    ASSERT_EQ(sys->cells.data[0].root_node_id, small_neg);
    ASSERT(alea_contains_point(sys, small_neg, 1.0, 0.0, 0.0));
    ASSERT(!alea_contains_point(sys, small_neg, 3.0, 0.0, 0.0));
    alea_destroy(sys);
}

TEST(proof_simplify_reports_witness_for_non_equivalent_candidate) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    alea_node_id_t a_pos, a_neg, b_pos, b_neg;
    create_sphere(sys, 111, -1.5, 0, 0, 2.0, &a_pos, &a_neg);
    create_sphere(sys, 112, 1.5, 0, 0, 2.0, &b_pos, &b_neg);
    alea_node_id_t root = alea_create_intersection(sys, a_neg, b_neg);
    ASSERT(alea_add_cell(sys, 2, root, ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    alea_cell_simplify_proof_options_t options;
    alea_cell_simplify_proof_options_init(&options);
    alea_cell_simplify_proof_result_t result;
    ASSERT_EQ(alea_cell_simplify_proven(sys, 0, &options, &result), 0);
    ASSERT(!result.changed);
    ASSERT(result.complete);
    ASSERT_EQ(result.candidates_disproven, 2);
    ASSERT(result.has_witness);
    alea_destroy(sys);
}

TEST(proof_simplify_explicit_bounds_are_a_caller_domain_assertion) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    alea_primitive_data_t p0 = {.plane = {1, 0, 0, 0}};
    alea_primitive_data_t p10 = {.plane = {1, 0, 0, -10}};
    int8_t inv0 = 0, inv10 = 0;
    alea_primitive_id_t id0 = alea_get_or_create_primitive(
        sys, ALEA_PRIMITIVE_PLANE, &p0, &inv0);
    alea_primitive_id_t id10 = alea_get_or_create_primitive(
        sys, ALEA_PRIMITIVE_PLANE, &p10, &inv10);
    alea_node_id_t n0 = alea_add_primitive_node(sys, id0, -1, inv0, 121);
    alea_node_id_t n10 = alea_add_primitive_node(sys, id10, -1, inv10, 122);
    alea_node_id_t root = alea_create_intersection(sys, n0, n10);
    ASSERT(alea_add_cell(sys, 3, root, ALEA_MATERIAL_VOID, 0.0, 0) >= 0);

    alea_cell_simplify_proof_options_t options;
    alea_cell_simplify_proof_options_init(&options);
    options.has_bounds = true;
    options.bounds = (alea_bbox_t){-1, 1, -1, 1, -1, 1};
    alea_cell_simplify_proof_result_t result;
    ASSERT_EQ(alea_cell_simplify_proven(sys, 0, &options, &result), 0);
    ASSERT(result.changed);
    ASSERT_EQ(result.bounds_source, ALEA_PROOF_BOUNDS_EXPLICIT);
    ASSERT(!result.bounds_verified);
    ASSERT_EQ(sys->cells.data[0].root_node_id, root);
    alea_destroy(sys);
}

TEST(proof_simplify_reports_certified_empty_without_deleting_cell) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    alea_node_id_t a_pos, a_neg, b_pos, b_neg;
    create_sphere(sys, 131, -4, 0, 0, 1.0, &a_pos, &a_neg);
    create_sphere(sys, 132, 4, 0, 0, 1.0, &b_pos, &b_neg);
    alea_node_id_t root = alea_create_intersection(sys, a_neg, b_neg);
    ASSERT(alea_add_cell(sys, 4, root, ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    alea_cell_simplify_proof_result_t result;
    ASSERT_EQ(alea_cell_simplify_proven(sys, 0, NULL, &result), 0);
    ASSERT(result.proven_empty);
    ASSERT(!result.changed);
    ASSERT(!result.applied);
    ASSERT_EQ(alea_vec_count(&sys->cells), 1);
    ASSERT_EQ(sys->cells.data[0].root_node_id, root);
    alea_destroy(sys);
}

TEST(proof_simplify_worker_count_preserves_deterministic_receipt) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    alea_node_id_t a_pos, a_neg, b_pos, b_neg;
    create_sphere(sys, 141, 0, 0, 0, 2.0, &a_pos, &a_neg);
    create_sphere(sys, 142, 0, 0, 0, 2.0, &b_pos, &b_neg);
    alea_node_id_t root = alea_create_intersection(sys, a_neg, b_neg);
    ASSERT(alea_add_cell(sys, 5, root, ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    alea_cell_simplify_proof_options_t options;
    alea_cell_simplify_proof_options_init(&options);
    options.max_depth = 3;
    options.requested_workers = 1;
    alea_cell_simplify_proof_result_t serial;
    ASSERT_EQ(alea_cell_simplify_proven(sys, 0, &options, &serial), 0);
    options.requested_workers = 2;
    alea_cell_simplify_proof_result_t parallel;
    ASSERT_EQ(alea_cell_simplify_proven(sys, 0, &options, &parallel), 0);
    ASSERT_EQ(serial.changed, parallel.changed);
    ASSERT_EQ(serial.complete, parallel.complete);
    ASSERT_EQ(serial.last_limit, parallel.last_limit);
    ASSERT_EQ(serial.proof_nodes, parallel.proof_nodes);
    ASSERT_EQ(serial.mixed_leaf_nodes, parallel.mixed_leaf_nodes);
    ASSERT_EQ(serial.has_witness, parallel.has_witness);
    alea_destroy(sys);
}

TEST(proof_simplify_single_item_batch_preserves_inner_parallelism) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    alea_node_id_t small_pos, small_neg, big_pos, big_neg;
    create_sphere(sys, 145, 0, 0, 0, 2.0, &small_pos, &small_neg);
    create_sphere(sys, 146, 0, 0, 0, 5.0, &big_pos, &big_neg);
    alea_node_id_t root = alea_create_intersection(sys, small_neg, big_neg);
    ASSERT(alea_add_cell(sys, 6, root, ALEA_MATERIAL_VOID, 0.0, 0) >= 0);

    alea_cell_simplify_request_t request = {
        .cell_index = 0, .has_bounds = false
    };
    alea_cell_simplify_proof_result_t result;
    alea_cells_simplify_proof_summary_t summary;
    alea_cells_simplify_proof_options_t options;
    alea_cells_simplify_proof_options_init(&options);
    options.max_depth = 5;
    options.requested_workers = 2;

    ASSERT_EQ(alea_cells_simplify_proven(
        sys, &request, 1, &options, &result, &summary), 0);
    const size_t expected_workers =
        alea_parallel_max_threads() >= 2 ? 2u : 1u;
    ASSERT_EQ(result.actual_workers, expected_workers);
    ASSERT_EQ(summary.actual_workers, expected_workers);
    ASSERT_EQ(result.parallel_batch_count, expected_workers > 1 ? 1 : 0);
    ASSERT_EQ(summary.parallel_batch_count, result.parallel_batch_count);
    alea_destroy(sys);
}

TEST(proof_simplify_batch_is_ordered_read_only_and_transactional) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    alea_node_id_t a_small_pos, a_small_neg, a_big_pos, a_big_neg;
    alea_node_id_t b_small_pos, b_small_neg, b_big_pos, b_big_neg;
    create_sphere(sys, 151, 0, 0, 0, 2.0,
                  &a_small_pos, &a_small_neg);
    create_sphere(sys, 152, 0, 0, 0, 5.0,
                  &a_big_pos, &a_big_neg);
    create_sphere(sys, 153, 20, 0, 0, 1.0,
                  &b_small_pos, &b_small_neg);
    create_sphere(sys, 154, 20, 0, 0, 4.0,
                  &b_big_pos, &b_big_neg);
    alea_node_id_t root_a = alea_create_intersection(
        sys, a_small_neg, a_big_neg);
    alea_node_id_t root_b = alea_create_intersection(
        sys, b_small_neg, b_big_neg);
    ASSERT(alea_add_cell(sys, 21, root_a, ALEA_MATERIAL_VOID, 0.0, 7) >= 0);
    ASSERT(alea_add_cell(sys, 22, root_b, ALEA_MATERIAL_VOID, 0.0, 9) >= 0);

    alea_cell_simplify_request_t requests[2] = {
        {.cell_index = 1, .has_bounds = false},
        {.cell_index = 0, .has_bounds = false},
    };
    alea_cell_simplify_proof_result_t results[2];
    alea_cells_simplify_proof_summary_t summary;
    alea_cells_simplify_proof_options_t options;
    alea_cells_simplify_proof_options_init(&options);
    options.max_depth = 5;
    options.requested_workers = 2;
    size_t node_count = alea_vec_count(&sys->nodes);
    uint64_t generation = alea_system_geometry_generation(sys);

    ASSERT_EQ(alea_cells_simplify_proven(
        sys, requests, 2, &options, results, &summary), 0);
    ASSERT_EQ(summary.selected_cells, 2);
    ASSERT_EQ(summary.changed_cells, 2);
    ASSERT_EQ(summary.applied_cells, 0);
    ASSERT(results[0].changed);
    ASSERT(results[1].changed);
    ASSERT_EQ(results[0].root_node_id, root_b);
    ASSERT_EQ(results[1].root_node_id, root_a);
    ASSERT(results[0].bounds.min_x > 10.0);
    ASSERT(results[1].bounds.max_x < 10.0);
    ASSERT_EQ(alea_vec_count(&sys->nodes), node_count);
    ASSERT_EQ(sys->cells.data[0].root_node_id, root_a);
    ASSERT_EQ(sys->cells.data[1].root_node_id, root_b);
    ASSERT_EQ(alea_system_geometry_generation(sys), generation);

    options.apply = true;
    ASSERT_EQ(alea_cells_simplify_proven(
        sys, requests, 2, &options, results, &summary), 0);
    ASSERT_EQ(summary.applied_cells, 2);
    ASSERT(results[0].applied);
    ASSERT(results[1].applied);
    ASSERT_EQ(sys->cells.data[0].root_node_id, a_small_neg);
    ASSERT_EQ(sys->cells.data[1].root_node_id, b_small_neg);
    ASSERT_EQ(alea_system_geometry_generation(sys), generation + 1);
    alea_destroy(sys);
}

TEST(proof_simplify_batch_rejects_duplicate_cells_before_work) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    alea_node_id_t small_pos, small_neg, big_pos, big_neg;
    create_sphere(sys, 161, 0, 0, 0, 2.0, &small_pos, &small_neg);
    create_sphere(sys, 162, 0, 0, 0, 5.0, &big_pos, &big_neg);
    alea_node_id_t root = alea_create_intersection(sys, small_neg, big_neg);
    ASSERT(alea_add_cell(sys, 31, root, ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    alea_cell_simplify_request_t requests[2] = {
        {.cell_index = 0, .has_bounds = false},
        {.cell_index = 0, .has_bounds = false},
    };
    alea_cell_simplify_proof_result_t results[2];
    alea_cells_simplify_proof_summary_t summary;
    size_t node_count = alea_vec_count(&sys->nodes);
    uint64_t generation = alea_system_geometry_generation(sys);
    ASSERT_EQ(alea_cells_simplify_proven(
        sys, requests, 2, NULL, results, &summary), -1);
    ASSERT_EQ(alea_vec_count(&sys->nodes), node_count);
    ASSERT_EQ(sys->cells.data[0].root_node_id, root);
    ASSERT_EQ(alea_system_geometry_generation(sys), generation);
    alea_destroy(sys);
}

TEST_MAIN()
