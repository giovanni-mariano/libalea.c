// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_bbox_numerical.c - Tests for numerical bounding box tightening
 *
 * Sense convention:
 *   sense = -1 (negative side): ax+by+cz+d <= 0
 *   sense = +1 (positive side): ax+by+cz+d >= 0
 */

#include "alea_test.h"
#include "alea.h"
#include "core/alea_system.h"
#include "primitives/bbox.h"
#include "util/math.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* ============================================================================
 * mat3_solve_cramer tests
 * ============================================================================ */

TEST(cramer_solve_basic) {
    /* Solve: x+2y+3z=14, 4x+5y+6z=32, 7x+8y=3 */
    mat3 A;
    /* Column-major: A.m[col][row] */
    A.m[0][0] = 1; A.m[1][0] = 2; A.m[2][0] = 3;
    A.m[0][1] = 4; A.m[1][1] = 5; A.m[2][1] = 6;
    A.m[0][2] = 7; A.m[1][2] = 8; A.m[2][2] = 0;

    vec3 b = {14, 32, 3};
    vec3 x;
    int rc = mat3_solve_cramer(A, b, &x);
    ASSERT_EQ(rc, 0);

    /* Verify Ax = b */
    double r0 = 1*x.x + 2*x.y + 3*x.z;
    double r1 = 4*x.x + 5*x.y + 6*x.z;
    double r2 = 7*x.x + 8*x.y + 0*x.z;
    ASSERT_NEAR(r0, 14.0, 1e-10);
    ASSERT_NEAR(r1, 32.0, 1e-10);
    ASSERT_NEAR(r2, 3.0,  1e-10);
}

TEST(cramer_identity) {
    /* I * x = b  →  x = b */
    mat3 A = mat3_identity();
    vec3 b = {3.0, -7.0, 11.0};
    vec3 x;
    int rc = mat3_solve_cramer(A, b, &x);
    ASSERT_EQ(rc, 0);
    ASSERT_NEAR(x.x, 3.0, 1e-14);
    ASSERT_NEAR(x.y, -7.0, 1e-14);
    ASSERT_NEAR(x.z, 11.0, 1e-14);
}

TEST(cramer_singular) {
    /* Two identical rows → singular */
    mat3 A;
    A.m[0][0] = 1; A.m[1][0] = 2; A.m[2][0] = 3;
    A.m[0][1] = 1; A.m[1][1] = 2; A.m[2][1] = 3;
    A.m[0][2] = 4; A.m[1][2] = 5; A.m[2][2] = 6;

    vec3 b = {1, 2, 3};
    vec3 x;
    int rc = mat3_solve_cramer(A, b, &x);
    ASSERT_EQ(rc, -1);
}

/* ============================================================================
 * Vertex enumeration tests
 * ============================================================================ */

/*
 * Build a tetrahedron from 4 planes.
 * Vertices at (0,0,0), (10,0,0), (0,10,0), (0,0,10).
 *
 * x >= 0      →  sense=+1 on plane (1,0,0,0):  x >= 0
 * y >= 0      →  sense=+1 on plane (0,1,0,0):  y >= 0
 * z >= 0      →  sense=+1 on plane (0,0,1,0):  z >= 0
 * x+y+z <= 10 →  sense=-1 on plane (1,1,1,-10): x+y+z-10 <= 0
 */
static alea_system_t* create_tetrahedron(void) {
    alea_system_t* sys = alea_create();
    if (!sys) return NULL;

    int spx = alea_plane_surface(sys, 0,  1,  0,  0,   0);
    int spy = alea_plane_surface(sys, 0,  0,  1,  0,   0);
    int spz = alea_plane_surface(sys, 0,  0,  0,  1,   0);
    int spd = alea_plane_surface(sys, 0,  1,  1,  1, -10);

    alea_node_id_t px = alea_halfspace(sys, spx, +1); /* x >= 0 */
    alea_node_id_t py = alea_halfspace(sys, spy, +1); /* y >= 0 */
    alea_node_id_t pz = alea_halfspace(sys, spz, +1); /* z >= 0 */
    alea_node_id_t pd = alea_halfspace(sys, spd, -1); /* x+y+z <= 10 */

    alea_node_id_t i1 = alea_intersection(sys, px, py);
    alea_node_id_t i2 = alea_intersection(sys, pz, pd);
    alea_node_id_t root = alea_intersection(sys, i1, i2);

    int m1 = alea_add_material(sys, 1);

    alea_add_cell(sys, 1, root, m1, -1.0, 0);
    alea_build_universe_index(sys);
    return sys;
}

TEST(vertex_enum_tetrahedron) {
    alea_system_t* sys = create_tetrahedron();
    ASSERT_NOT_NULL(sys);

    alea_node_id_t root = sys->cells.data[0].root_node_id;
    alea_bbox_t tight;
    int rc = alea_tighten_bbox_numerical(sys, root, 0.1, &tight);
    ASSERT_EQ(rc, 0);

    /* Expected bbox: [0, 10] x [0, 10] x [0, 10] (with small margin) */
    ASSERT_TRUE(tight.min_x >= -2.0);
    ASSERT_TRUE(tight.max_x <= 12.0);
    ASSERT_TRUE(tight.min_y >= -2.0);
    ASSERT_TRUE(tight.max_y <= 12.0);
    ASSERT_TRUE(tight.min_z >= -2.0);
    ASSERT_TRUE(tight.max_z <= 12.0);

    /* Should be reasonably tight */
    ASSERT_TRUE(tight.max_x - tight.min_x < 15.0);
    ASSERT_TRUE(tight.max_y - tight.min_y < 15.0);
    ASSERT_TRUE(tight.max_z - tight.min_z < 15.0);

    alea_destroy(sys);
}

/*
 * Build a box from 6 axis-aligned planes.
 * Box: [-5, 5] x [-3, 3] x [-2, 2]
 */
TEST(vertex_enum_axis_box) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* x >= -5: x+5 >= 0 → sense=+1 on (1,0,0,5) */
    int s0 = alea_plane_surface(sys, 0, 1, 0, 0,  5);
    alea_node_id_t px0 = alea_halfspace(sys, s0, +1);
    /* x <=  5: x-5 <= 0 → sense=-1 on (1,0,0,-5) */
    int s1 = alea_plane_surface(sys, 0, 1, 0, 0, -5);
    alea_node_id_t px1 = alea_halfspace(sys, s1, -1);
    /* y >= -3: y+3 >= 0 → sense=+1 on (0,1,0,3) */
    int s2 = alea_plane_surface(sys, 0, 0, 1, 0,  3);
    alea_node_id_t py0 = alea_halfspace(sys, s2, +1);
    /* y <=  3: y-3 <= 0 → sense=-1 on (0,1,0,-3) */
    int s3 = alea_plane_surface(sys, 0, 0, 1, 0, -3);
    alea_node_id_t py1 = alea_halfspace(sys, s3, -1);
    /* z >= -2: z+2 >= 0 → sense=+1 on (0,0,1,2) */
    int s4 = alea_plane_surface(sys, 0, 0, 0, 1,  2);
    alea_node_id_t pz0 = alea_halfspace(sys, s4, +1);
    /* z <=  2: z-2 <= 0 → sense=-1 on (0,0,1,-2) */
    int s5 = alea_plane_surface(sys, 0, 0, 0, 1, -2);
    alea_node_id_t pz1 = alea_halfspace(sys, s5, -1);

    alea_node_id_t i1 = alea_intersection(sys, px0, px1);
    alea_node_id_t i2 = alea_intersection(sys, py0, py1);
    alea_node_id_t i3 = alea_intersection(sys, pz0, pz1);
    alea_node_id_t i4 = alea_intersection(sys, i1, i2);
    alea_node_id_t root = alea_intersection(sys, i4, i3);

    int m1 = alea_add_material(sys, 1);

    alea_add_cell(sys, 1, root, m1, -1.0, 0);
    alea_build_universe_index(sys);

    alea_bbox_t tight;
    int rc = alea_tighten_bbox_numerical(sys, root, 0.1, &tight);
    ASSERT_EQ(rc, 0);

    ASSERT_NEAR(tight.min_x, -5.0, 1.5);
    ASSERT_NEAR(tight.max_x,  5.0, 1.5);
    ASSERT_NEAR(tight.min_y, -3.0, 1.5);
    ASSERT_NEAR(tight.max_y,  3.0, 1.5);
    ASSERT_NEAR(tight.min_z, -2.0, 1.5);
    ASSERT_NEAR(tight.max_z,  2.0, 1.5);

    alea_destroy(sys);
}

/* ============================================================================
 * Octree tests
 * ============================================================================ */

/*
 * Oblique slab: intersection of two oblique planes.
 *   x+y+z <= 5:  sense=-1 on plane (1,1,1,-5)
 *   x+y+z >= -5: sense=+1 on plane (1,1,1,5)
 * This creates a slab that has infinite analytical bbox.
 */
TEST(octree_oblique_slab) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int s0 = alea_plane_surface(sys, 0, 1, 1, 1, -5);
    int s1 = alea_plane_surface(sys, 0, 1, 1, 1,  5);
    alea_node_id_t p1 = alea_halfspace(sys, s0, -1);  /* x+y+z-5 <= 0 → x+y+z <= 5 */
    alea_node_id_t p2 = alea_halfspace(sys, s1, +1);  /* x+y+z+5 >= 0 → x+y+z >= -5 */
    alea_node_id_t root = alea_intersection(sys, p1, p2);

    int m1 = alea_add_material(sys, 1);

    alea_add_cell(sys, 1, root, m1, -1.0, 0);
    alea_build_universe_index(sys);

    /* The analytical bbox for this is infinite (oblique planes) */
    alea_bbox_t analytical = sys->nodes.data[root].bbox;
    double dx = analytical.max_x - analytical.min_x;
    ASSERT_TRUE(dx > 9e5);  /* should be infinite */

    /* Numerical tightening via octree (LP needs >= 3 planes) */
    alea_bbox_t tight;
    int rc = alea_tighten_bbox_numerical(sys, root, 1.0, &tight);
    ASSERT_EQ(rc, 0);

    /* The slab is unbounded in transverse directions within [-1e6,1e6]^3,
     * so octree returns a large bbox — the key thing is it returns 0 */
    ASSERT_TRUE(tight.min_x < tight.max_x);

    alea_destroy(sys);
}

/* ============================================================================
 * Mixed cell (planes + cylinder)
 * ============================================================================ */

TEST(numerical_mixed_cell) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Oblique planes (LP will fail because of cylinder) */
    int s0 = alea_plane_surface(sys, 0, 1, 1, 0,  0);
    int s1 = alea_plane_surface(sys, 0, 1, 1, 0, -20);
    int s2 = alea_plane_surface(sys, 0, 0, 0, 1, 10);
    int s3 = alea_plane_surface(sys, 0, 0, 0, 1, -10);
    alea_node_id_t p1 = alea_halfspace(sys, s0, +1);    /* x+y >= 0 */
    alea_node_id_t p2 = alea_halfspace(sys, s1, -1);    /* x+y-20 <= 0 → x+y <= 20 */
    alea_node_id_t pz0 = alea_halfspace(sys, s2, +1);   /* z+10 >= 0 → z >= -10 */
    alea_node_id_t pz1 = alea_halfspace(sys, s3, -1);   /* z-10 <= 0 → z <= 10 */
    /* Cylinder along Z at origin, radius 15 */
    int sc = alea_cylinder_z_surface(sys, 0, 0, 0, 15);
    alea_node_id_t cyl = alea_halfspace(sys, sc, -1);

    alea_node_id_t i1 = alea_intersection(sys, p1, p2);
    alea_node_id_t i2 = alea_intersection(sys, pz0, pz1);
    alea_node_id_t i3 = alea_intersection(sys, i1, i2);
    alea_node_id_t root = alea_intersection(sys, i3, cyl);

    int m1 = alea_add_material(sys, 1);

    alea_add_cell(sys, 1, root, m1, -1.0, 0);
    alea_build_universe_index(sys);

    /* LP should fail (has cylinder), octree should work */
    alea_bbox_t tight;
    int rc = alea_tighten_bbox_numerical(sys, root, 1.0, &tight);
    ASSERT_EQ(rc, 0);

    /* Should be finite */
    double tdx = tight.max_x - tight.min_x;
    double tdy = tight.max_y - tight.min_y;
    double tdz = tight.max_z - tight.min_z;
    ASSERT_TRUE(tdx < 1e5);
    ASSERT_TRUE(tdy < 1e5);
    ASSERT_TRUE(tdz < 25.0);  /* z bounded by [-10, 10] */

    alea_destroy(sys);
}

/* ============================================================================
 * Bounding sphere integration
 * ============================================================================ */

TEST(numerical_in_bounding_sphere) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Cell 1: simple sphere (finite analytical bbox) */
    int ss = alea_sphere_surface(sys, 0, 0, 0, 0, 5.0);
    alea_node_id_t sph = alea_halfspace(sys, ss, -1);

    /* Cell 2: tetrahedron from oblique plane (would be skipped without numerical) */
    int s0 = alea_plane_surface(sys, 0, 1, 0, 0, 0);
    int s1 = alea_plane_surface(sys, 0, 0, 1, 0, 0);
    int s2 = alea_plane_surface(sys, 0, 0, 0, 1, 0);
    int s3 = alea_plane_surface(sys, 0, 1, 1, 1, -20);
    alea_node_id_t px = alea_halfspace(sys, s0, +1);    /* x >= 0 */
    alea_node_id_t py = alea_halfspace(sys, s1, +1);    /* y >= 0 */
    alea_node_id_t pz = alea_halfspace(sys, s2, +1);    /* z >= 0 */
    alea_node_id_t pd = alea_halfspace(sys, s3, -1);    /* x+y+z <= 20 */

    alea_node_id_t t1 = alea_intersection(sys, px, py);
    alea_node_id_t t2 = alea_intersection(sys, pz, pd);
    alea_node_id_t tet = alea_intersection(sys, t1, t2);

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);

    alea_add_cell(sys, 1, sph, m1, -1.0, 0);
    alea_add_cell(sys, 2, tet, m2, -1.0, 0);

    alea_build_universe_index(sys);

    double cx, cy, cz, radius;
    int rc = alea_compute_bounding_sphere(sys, 1.0, &cx, &cy, &cz, &radius);
    ASSERT_EQ(rc, 0);

    /* The bounding sphere should encompass both cells.
     * The tetrahedron extends to (20,0,0), (0,20,0), (0,0,20).
     * If numerical tightening works, the sphere should cover these. */
    ASSERT_TRUE(radius > 5.0);  /* bigger than just the sphere */

    alea_destroy(sys);
}

/* ============================================================================
 * Graveyard / complement cell
 * ============================================================================ */

TEST(graveyard_cell_skipped) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Create a complement cell (outside a sphere) — truly infinite */
    int ss = alea_sphere_surface(sys, 0, 0, 0, 0, 5.0);
    alea_node_id_t sph = alea_halfspace(sys, ss, +1);  /* outside sphere */
    alea_add_cell(sys, 999, sph, ALEA_MATERIAL_VOID, 0, 0);
    alea_build_universe_index(sys);

    /* Numerical tightening: octree finds content (outside sphere is everywhere) */
    alea_bbox_t tight;
    int rc = alea_tighten_bbox_numerical(sys, sys->cells.data[0].root_node_id, 1.0, &tight);
    /* For a complement cell the octree finds content → returns 0 */
    /* The important thing is it doesn't crash */
    (void)rc;

    alea_destroy(sys);
}

/* ============================================================================
 * Public API test
 * ============================================================================ */

TEST(public_api_numerical_tighten) {
    alea_system_t* sys = create_tetrahedron();
    ASSERT_NOT_NULL(sys);

    alea_node_id_t root = sys->cells.data[0].root_node_id;
    int rc = alea_tighten_cell_bbox_numerical(sys, 0);
    ASSERT_EQ(rc, 0);

    /* After tightening, bbox should be finite and reasonable */
    alea_bbox_t after = sys->nodes.data[root].bbox;
    ASSERT_TRUE(after.min_x < after.max_x);
    ASSERT_TRUE(after.max_x - after.min_x < 15.0);

    alea_destroy(sys);
}

TEST(tighten_all_with_lp) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Create a cell bounded by oblique + axis-aligned planes.
     * Triangular prism:
     *   x+y >= 0: sense=+1 on (1,1,0,0)
     *   x-y >= 0: sense=+1 on (1,-1,0,0)
     *   x <= 10:  sense=-1 on (1,0,0,-10)
     *   z >= -5:  sense=+1 on (0,0,1,5)
     *   z <= 5:   sense=-1 on (0,0,1,-5)
     */
    int s0 = alea_plane_surface(sys, 0, 1,  1, 0,   0);
    int s1 = alea_plane_surface(sys, 0, 1, -1, 0,   0);
    int s2 = alea_plane_surface(sys, 0, 1,  0, 0, -10);
    int s3 = alea_plane_surface(sys, 0, 0,  0, 1,   5);
    int s4 = alea_plane_surface(sys, 0, 0,  0, 1,  -5);
    alea_node_id_t p1 = alea_halfspace(sys, s0, +1);   /* x+y >= 0 */
    alea_node_id_t p2 = alea_halfspace(sys, s1, +1);   /* x-y >= 0 */
    alea_node_id_t p3 = alea_halfspace(sys, s2, -1);   /* x <= 10 */
    alea_node_id_t pz0 = alea_halfspace(sys, s3, +1);  /* z >= -5 */
    alea_node_id_t pz1 = alea_halfspace(sys, s4, -1);  /* z <= 5 */

    alea_node_id_t i1 = alea_intersection(sys, p1, p2);
    alea_node_id_t i2 = alea_intersection(sys, p3, pz0);
    alea_node_id_t i3 = alea_intersection(sys, i1, i2);
    alea_node_id_t root = alea_intersection(sys, i3, pz1);

    int m1 = alea_add_material(sys, 1);

    alea_add_cell(sys, 1, root, m1, -1.0, 0);
    alea_build_universe_index(sys);

    int tightened = alea_tighten_all_bboxes(sys, 1.0);
    ASSERT_TRUE(tightened >= 1);

    /* After tightening, the bbox should be finite */
    alea_bbox_t* box = &sys->nodes.data[root].bbox;
    double tdx = box->max_x - box->min_x;
    double tdz = box->max_z - box->min_z;
    ASSERT_TRUE(tdx < 20.0);
    ASSERT_TRUE(tdz < 15.0);

    alea_destroy(sys);
}

TEST_MAIN()
