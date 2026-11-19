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
    alea_bbox_t analytical = alea_node_bbox_get(&sys->nodes.data[root].bbox);
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
    alea_bbox_t after = alea_node_bbox_get(&sys->nodes.data[root].bbox);
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
    alea_bbox_t box_v = alea_node_bbox_get(&sys->nodes.data[root].bbox);
    alea_bbox_t* box = &box_v;
    double tdx = box->max_x - box->min_x;
    double tdz = box->max_z - box->min_z;
    ASSERT_TRUE(tdx < 20.0);
    ASSERT_TRUE(tdz < 15.0);

    alea_destroy(sys);
}

TEST(cell_volume_box_exact_and_read_only) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int surface = alea_box_surface(sys, 0, -1.0, 1.0, -2.0, 2.0, -3.0, 3.0);
    ASSERT(surface >= 0);
    alea_node_id_t root = alea_halfspace(sys, surface, -1);
    ASSERT(alea_add_cell(sys, 10, root, ALEA_MATERIAL_VOID, 0.0, 7) >= 0);
    alea_bbox_t before = alea_node_bbox_get(&sys->nodes.data[root].bbox);

    alea_cell_volume_options_t options;
    alea_cell_volume_options_init(&options);
    options.requested_workers = 1;
    alea_cell_volume_result_t result;
    ASSERT_EQ(alea_cell_estimate_volume(sys, 0, &options, &result), 0);
    ASSERT_NEAR(result.volume, 48.0, 1e-12);
    ASSERT_NEAR(result.lower_bound, 48.0, 1e-12);
    ASSERT_NEAR(result.upper_bound, 48.0, 1e-12);
    ASSERT(result.converged);
    ASSERT(result.complete_cell_domain);
    ASSERT_EQ(result.bounds_source, ALEA_CELL_VOLUME_BOUNDS_STORED);
    alea_bbox_t after = alea_node_bbox_get(&sys->nodes.data[root].bbox);
    ASSERT_EQ(memcmp(&before, &after, sizeof(before)), 0);
    alea_destroy(sys);
}

TEST(cell_volume_sphere_bounds_and_parallel_determinism) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int surface = alea_sphere_surface(sys, 0, 0.0, 0.0, 0.0, 2.0);
    ASSERT(surface >= 0);
    alea_node_id_t root = alea_surface_at(sys, surface)->neg_node;
    ASSERT(alea_add_cell(sys, 1, root, ALEA_MATERIAL_VOID, 0.0, 0) >= 0);

    alea_cell_volume_options_t options;
    alea_cell_volume_options_init(&options);
    options.relative_tolerance = 0.0;
    options.max_depth = 7;
    options.requested_workers = 1;
    alea_cell_volume_result_t serial, parallel;
    ASSERT_EQ(alea_cell_estimate_volume(sys, 0, &options, &serial), 0);
    options.requested_workers = 4;
    ASSERT_EQ(alea_cell_estimate_volume(sys, 0, &options, &parallel), 0);
    double exact = 4.0 * 3.14159265358979323846 * 8.0 / 3.0;
    ASSERT(serial.lower_bound <= exact);
    ASSERT(serial.upper_bound >= exact);
    ASSERT(serial.volume >= serial.lower_bound);
    ASSERT(serial.volume <= serial.upper_bound);
    ASSERT_NEAR(serial.volume, parallel.volume, 0.0);
    ASSERT_NEAR(serial.lower_bound, parallel.lower_bound, 0.0);
    ASSERT_NEAR(serial.upper_bound, parallel.upper_bound, 0.0);
    ASSERT_EQ(serial.total_nodes, parallel.total_nodes);
    alea_destroy(sys);
}

TEST(cell_volume_budgets_keep_complete_partition) {
    alea_system_t* sys = alea_create();
    int surface = alea_sphere_surface(sys, 0, 0.0, 0.0, 0.0, 1.0);
    ASSERT(surface >= 0);
    ASSERT(alea_add_cell(sys, 1, alea_surface_at(sys, surface)->neg_node,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    alea_cell_volume_options_t options;
    alea_cell_volume_options_init(&options);
    options.relative_tolerance = 0.0;
    options.max_evaluations = 1;
    alea_cell_volume_result_t result;
    ASSERT_EQ(alea_cell_estimate_volume(sys, 0, &options, &result), 0);
    ASSERT_EQ(result.interval_evaluations, 1);
    ASSERT(result.evaluation_limit_reached);
    ASSERT(result.lower_bound <= 4.0 * 3.14159265358979323846 / 3.0);
    ASSERT(result.upper_bound >= 4.0 * 3.14159265358979323846 / 3.0);
    ASSERT_EQ(result.unresolved_volume, 8.0);
    options.max_evaluations = 100000;
    options.max_memory_bytes = 1024;
    ASSERT_EQ(alea_cell_estimate_volume(sys, 0, &options, &result), 0);
    ASSERT(result.memory_limit_reached);
    ASSERT(result.peak_memory_bytes <= options.max_memory_bytes);
    ASSERT(result.lower_bound <= 4.0 * 3.14159265358979323846 / 3.0);
    ASSERT(result.upper_bound >= 4.0 * 3.14159265358979323846 / 3.0);
    options.max_evaluations = 0;
    ASSERT(alea_cell_estimate_volume(sys, 0, &options, &result) != 0);
    alea_destroy(sys);
}

TEST(cell_volume_interruption_returns_no_partial_success) {
    alea_system_t* sys = alea_create();
    int surface = alea_sphere_surface(sys, 0, 0.0, 0.0, 0.0, 1.0);
    ASSERT(surface >= 0);
    ASSERT(alea_add_cell(sys, 1, alea_surface_at(sys, surface)->neg_node,
                         ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    alea_cell_volume_options_t options;
    alea_cell_volume_options_init(&options);
    options.max_evaluations = 1;
    alea_cell_volume_result_t result;
    alea_interrupt();
    int rc = alea_cell_estimate_volume(sys, 0, &options, &result);
    alea_clear_interrupt();
    ASSERT_EQ(rc, -1);
    ASSERT_EQ(alea_cell_estimate_volume(sys, 0, &options, &result), 0);
    ASSERT(result.evaluation_limit_reached);
    alea_destroy(sys);
}

TEST(cell_volume_explicit_bounds_clip_unbounded_cell) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int surface = alea_plane_surface(sys, 0, 1.0, 0.0, 0.0, 0.0);
    ASSERT(surface >= 0);
    alea_node_id_t root = alea_halfspace(sys, surface, -1);
    ASSERT(alea_add_cell(sys, 1, root, ALEA_MATERIAL_VOID, 0.0, 0) >= 0);

    alea_cell_volume_options_t options;
    alea_cell_volume_options_init(&options);
    options.has_bounds = true;
    options.bounds = (alea_bbox_t){-1, 1, -1, 1, -1, 1};
    options.relative_tolerance = 0.0;
    options.max_depth = 4;
    options.requested_workers = 1;
    alea_cell_volume_result_t result;
    ASSERT_EQ(alea_cell_estimate_volume(sys, 0, &options, &result), 0);
    ASSERT_NEAR(result.volume, 4.0, 1e-12);
    ASSERT_NEAR(result.lower_bound, 4.0, 1e-12);
    ASSERT_NEAR(result.upper_bound, 4.0, 1e-12);
    ASSERT_EQ(result.bounds_source, ALEA_CELL_VOLUME_BOUNDS_EXPLICIT);
    ASSERT(!result.complete_cell_domain);

    options.has_bounds = false;
    ASSERT_EQ(alea_cell_estimate_volume(sys, 0, &options, &result), -1);
    alea_destroy(sys);
}

TEST(cell_volume_discovers_bounded_plane_constraints) {
    alea_system_t* sys = create_tetrahedron();
    ASSERT_NOT_NULL(sys);
    alea_node_id_t root = sys->cells.data[0].root_node_id;
    alea_bbox_t unbounded = {-1e30, 1e30, -1e30, 1e30, -1e30, 1e30};
    alea_node_bbox_set(&sys->nodes.data[root].bbox, &unbounded);

    alea_cell_volume_options_t options;
    alea_cell_volume_options_init(&options);
    options.relative_tolerance = 0.0;
    options.max_depth = 7;
    options.requested_workers = 1;
    alea_cell_volume_result_t result;
    ASSERT_EQ(alea_cell_estimate_volume(sys, 0, &options, &result), 0);
    ASSERT_EQ(result.bounds_source, ALEA_CELL_VOLUME_BOUNDS_PLANE_CONSTRAINTS);
    ASSERT(result.complete_cell_domain);
    ASSERT(result.lower_bound <= 1000.0 / 6.0);
    ASSERT(result.upper_bound >= 1000.0 / 6.0);
    alea_destroy(sys);
}

TEST(cell_volume_adaptive_search_recovers_missing_stored_bbox) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int surface = alea_sphere_surface(sys, 0, 0.0, 0.0, 0.0, 1.0);
    ASSERT(surface >= 0);
    alea_node_id_t root = alea_halfspace(sys, surface, -1);
    ASSERT(alea_add_cell(sys, 1, root, ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    alea_bbox_t unbounded = {-1e10, 1e10, -1e10, 1e10, -1e10, 1e10};
    alea_node_bbox_set(&sys->nodes.data[root].bbox, &unbounded);

    alea_cell_volume_options_t options;
    alea_cell_volume_options_init(&options);
    options.max_depth = 5;
    alea_cell_volume_result_t result;
    ASSERT_EQ(alea_cell_estimate_volume(sys, 0, &options, &result), 0);
    ASSERT_EQ(result.bounds_source, ALEA_CELL_VOLUME_BOUNDS_ADAPTIVE_SEARCH);
    ASSERT(!result.complete_cell_domain);
    ASSERT(result.bounds.min_x <= -1.0 && result.bounds.max_x >= 1.0);
    ASSERT(result.bounds.min_y <= -1.0 && result.bounds.max_y >= 1.0);
    ASSERT(result.bounds.min_z <= -1.0 && result.bounds.max_z >= 1.0);
    alea_destroy(sys);
}

/* More than the former 64-plane limit, with redundant oblique faces. */
TEST(plane_lp_many_constraints_and_conservative_enclosure) {
    alea_system_t* sys=create_tetrahedron();
    ASSERT_NOT_NULL(sys);
    alea_node_id_t root=sys->cells.data[0].root_node_id;
    for(int i=0;i<128;i++) {
        double a=1+(i%7),b=1+(i%11),c=1+(i%13);
        int sid=alea_plane_surface(sys,0,a,b,c,-10*fmax(a,fmax(b,c))-1);
        root=alea_intersection(sys,root,alea_halfspace(sys,sid,-1));
    }
    alea_bbox_t out;
    ASSERT_EQ(alea_tighten_bbox_plane_constraints(sys,root,1e-8,&out),0);
    ASSERT_TRUE(out.min_x<=0 && out.min_y<=0 && out.min_z<=0);
    ASSERT_TRUE(out.max_x>=10 && out.max_y>=10 && out.max_z>=10);
    ASSERT_TRUE(out.max_x<10.01 && out.min_x>-.01);
    alea_destroy(sys);
}

TEST(plane_lp_rejects_unbounded_polyhedron_with_vertex) {
    alea_system_t* sys=alea_create();
    int x=alea_plane_surface(sys,0,1,0,0,0);
    int y=alea_plane_surface(sys,0,0,1,0,0);
    int z=alea_plane_surface(sys,0,0,0,1,0);
    alea_node_id_t root=alea_intersection(sys,alea_halfspace(sys,x,1),
        alea_intersection(sys,alea_halfspace(sys,y,1),alea_halfspace(sys,z,1)));
    /* A redundant fourth face gives a vertex but no upper bounds. */
    int q=alea_plane_surface(sys,0,1,1,1,1);
    root=alea_intersection(sys,root,alea_halfspace(sys,q,1));
    alea_bbox_t out={1,2,3,4,5,6},before=out;
    ASSERT_EQ(alea_tighten_bbox_plane_constraints(sys,root,1e-8,&out),-1);
    ASSERT_EQ(memcmp(&before,&out,sizeof(out)),0);
    alea_destroy(sys);
}

TEST(plane_lp_large_translated_rotated_box) {
    alea_system_t* sys=alea_create();
    double center[3]={1e9,-2e9,3e9};
    double axes[3][3]={{.6,.8,0},{-.8,.6,0},{0,0,1}};
    alea_node_id_t root=ALEA_NODE_ID_INVALID;
    for(int a=0;a<3;a++)for(int side=-1;side<=1;side+=2) {
        double d=-(axes[a][0]*center[0]+axes[a][1]*center[1]+axes[a][2]*center[2]+side*10);
        int id=alea_plane_surface(sys,0,axes[a][0],axes[a][1],axes[a][2],d);
        alea_node_id_t half=alea_halfspace(sys,id,-side);
        root=root==ALEA_NODE_ID_INVALID?half:alea_intersection(sys,root,half);
    }
    alea_bbox_t out;
    ASSERT_EQ(alea_tighten_bbox_plane_constraints(sys,root,1e-8,&out),0);
    ASSERT_TRUE(out.min_x<=center[0]-14 && out.max_x>=center[0]+14);
    ASSERT_TRUE(out.min_y<=center[1]-14 && out.max_y>=center[1]+14);
    ASSERT_TRUE(out.min_z<=center[2]-10 && out.max_z>=center[2]+10);
    ASSERT_TRUE(out.max_x-out.min_x<29);
    alea_destroy(sys);
}

TEST(plane_lp_empty_and_interrupted_leave_output_untouched) {
    alea_system_t* sys=create_tetrahedron();
    int sid=alea_plane_surface(sys,0,1,1,1,-20);
    alea_node_id_t root=alea_intersection(sys,sys->cells.data[0].root_node_id,alea_halfspace(sys,sid,1));
    alea_bbox_t out={1,2,3,4,5,6},before=out;
    ASSERT_EQ(alea_tighten_bbox_plane_constraints(sys,root,1e-8,&out),-1);
    ASSERT_EQ(memcmp(&before,&out,sizeof(out)),0);
    alea_interrupt();
    int interrupted=alea_tighten_bbox_plane_constraints(sys,sys->cells.data[0].root_node_id,1e-8,&out);
    alea_clear_interrupt();
    ASSERT_EQ(interrupted,-1);
    ASSERT_EQ(memcmp(&before,&out,sizeof(out)),0);
    alea_destroy(sys);
}

TEST(mixed_cylinder_oblique_caps_enclose_analytic_extrema) {
    alea_system_t* sys=alea_create();
    int cylinder=alea_cylinder_z_surface(sys,0,0,0,2);
    int lower=alea_plane_surface(sys,0,-1,0,1,0);
    int upper=alea_plane_surface(sys,0,-1,0,1,-3);
    alea_node_id_t root=alea_intersection(sys,alea_halfspace(sys,cylinder,-1),
        alea_intersection(sys,alea_halfspace(sys,lower,1),alea_halfspace(sys,upper,-1)));
    alea_bbox_t out={0};
    ASSERT_EQ(alea_tighten_bbox_plane_constraints(sys,root,1e-8,&out),-1);
    ASSERT_EQ(alea_tighten_bbox_mixed_constraints(sys,root,&out),0);
    ASSERT(out.min_x<=-2 && out.max_x>=2 && out.min_y<=-2 && out.max_y>=2);
    ASSERT(out.min_z<=-2 && out.max_z>=5);
    ASSERT(out.max_z<5.000001 && out.min_z> -2.000001);
    /* Independently generated points on both oblique caps and the cylinder. */
    for(int i=0;i<100;i++) {
        double angle=i*6.283185307179586/100,x=2*cos(angle),y=2*sin(angle);
        ASSERT(alea_bbox_contains_point(&out,x,y,x));
        ASSERT(alea_bbox_contains_point(&out,x,y,x+3));
    }
    alea_destroy(sys);
}

TEST(mixed_constraints_boolean_relaxation_and_unbounded_fallback) {
    alea_system_t* sys=alea_create();
    int sphere=alea_sphere_surface(sys,0,10,20,30,2);
    int remote=alea_sphere_surface(sys,0,100,200,300,3);
    alea_node_id_t inside=alea_halfspace(sys,sphere,-1);
    alea_node_id_t other=alea_halfspace(sys,remote,-1);
    alea_bbox_t out={1,2,3,4,5,6},before=out;
    /* Taking just one union branch would incorrectly exclude the other. */
    ASSERT_EQ(alea_tighten_bbox_mixed_constraints(sys,alea_union(sys,inside,other),&out),-1);
    ASSERT_EQ(memcmp(&before,&out,sizeof(out)),0);
    ASSERT_EQ(alea_tighten_bbox_mixed_constraints(sys,alea_halfspace(sys,sphere,1),&out),-1);
    ASSERT_EQ(memcmp(&before,&out,sizeof(out)),0);
    alea_node_id_t diff=alea_difference(sys,inside,other);
    ASSERT_EQ(alea_tighten_bbox_mixed_constraints(sys,diff,&out),0);
    ASSERT(out.min_x<=8 && out.max_x>=12 && out.min_z<=28 && out.max_z>=32);
    ASSERT(out.max_x<12.000001);
    alea_node_id_t restricted=alea_intersection(sys,inside,alea_union(sys,inside,other));
    ASSERT_EQ(alea_tighten_bbox_mixed_constraints(sys,restricted,&out),0);
    ASSERT(out.min_x<=8 && out.max_x>=12);
    int cylinder=alea_cylinder_z_surface(sys,0,0,0,2);
    out=before;
    ASSERT_EQ(alea_tighten_bbox_mixed_constraints(sys,alea_halfspace(sys,cylinder,-1),&out),-1);
    ASSERT_EQ(memcmp(&before,&out,sizeof(out)),0);
    alea_interrupt();
    int status=alea_tighten_bbox_mixed_constraints(sys,inside,&out);
    alea_clear_interrupt();
    ASSERT_EQ(status,-1);
    ASSERT_EQ(memcmp(&before,&out,sizeof(out)),0);
    alea_destroy(sys);
}

TEST(tightened_bounds_preserve_surface_references_and_rebuild_spatial_caches) {
    alea_system_t* sys=create_tetrahedron();
    ASSERT_NOT_NULL(sys);
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    const void* surface_indices=sys->cells.data[0].surface_indices;
    ASSERT_NOT_NULL(surface_indices);
    ASSERT(alea_tighten_all_bboxes(sys,1e-6)>0);
    ASSERT(alea_system_query_cache_ready(sys,ALEA_CACHE_CELL_SURFACES));
    ASSERT(!alea_system_query_cache_ready(sys,ALEA_CACHE_HIER_SPATIAL));
    ASSERT(surface_indices==sys->cells.data[0].surface_indices);
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    int cell,material;
    ASSERT_EQ(alea_find_cell_at(sys,1,1,1,&cell,&material),0);
    ASSERT(alea_find_cell_at(sys,11,1,1,&cell,&material)!=0);
    alea_destroy(sys);
}

TEST(mixed_quadric_slab_encloses_rotated_cylinder) {
    alea_system_t* sys=alea_create();ASSERT_NOT_NULL(sys);
    /* (x-y)^2 + 2*z^2 <= 8, 0 <= x+y <= 10. Exact coefficients avoid
     * an assumed rank or eigenvalue threshold in the independent oracle. */
    int q=alea_quadric_surface(sys,0,1,1,2,-2,0,0,0,0,0,-8);
    int low=alea_plane_surface(sys,0,1,1,0,0);
    int high=alea_plane_surface(sys,0,1,1,0,-10);
    alea_node_id_t cylinder=alea_halfspace(sys,q,-1);
    alea_node_id_t root=alea_intersection(sys,cylinder,
        alea_intersection(sys,alea_halfspace(sys,low,1),alea_halfspace(sys,high,-1)));
    alea_bbox_t untouched={1,2,3,4,5,6},box=untouched;
    ASSERT_EQ(alea_tighten_bbox_mixed_constraints(sys,cylinder,&box),-1);
    ASSERT_EQ(memcmp(&box,&untouched,sizeof(box)),0);
    ASSERT_EQ(alea_tighten_bbox_plane_constraints(sys,root,1e-8,&box),-1);
    ASSERT_EQ(alea_tighten_bbox_mixed_constraints(sys,root,&box),0);
    ASSERT(box.min_x> -3 && box.max_x<8 && box.min_z> -5 && box.max_z<5);
    for(int t=0;t<=10;t++)for(int i=0;i<360;i++) {
        double angle=i*6.283185307179586/360;
        double d=sqrt(8.0)*cos(angle),z=2*sin(angle);
        ASSERT(alea_bbox_contains_point(&box,(t+d)/2,(t-d)/2,z));
    }
    /* Reversing one cap gives an unbounded half-cylinder. */
    root=alea_intersection(sys,cylinder,
        alea_intersection(sys,alea_halfspace(sys,low,1),alea_halfspace(sys,high,1)));
    box=untouched;
    ASSERT_EQ(alea_tighten_bbox_mixed_constraints(sys,root,&box),-1);
    ASSERT_EQ(memcmp(&box,&untouched,sizeof(box)),0);
    alea_destroy(sys);
}

TEST(mixed_quadric_global_ellipsoid_and_exterior_fallback) {
    alea_system_t* sys=alea_create();ASSERT_NOT_NULL(sys);
    /* 4*(x-10)^2 + (y+20)^2 + 16*(z-30)^2 <= 16. */
    int q=alea_quadric_surface(sys,0,4,1,16,0,0,0,-80,40,-960,15184);
    alea_node_id_t root=alea_halfspace(sys,q,-1);
    alea_bbox_t box={0};
    ASSERT_EQ(alea_tighten_bbox_mixed_constraints(sys,root,&box),0);
    ASSERT(box.min_x<=8 && box.max_x>=12 && box.min_y<=-24 && box.max_y>=-16);
    ASSERT(box.min_z<=29 && box.max_z>=31);
    ASSERT(box.min_x>7.999999 && box.max_z<31.000001);
    alea_bbox_t before=box;
    ASSERT_EQ(alea_tighten_bbox_mixed_constraints(sys,alea_halfspace(sys,q,1),&box),-1);
    ASSERT_EQ(memcmp(&box,&before,sizeof(box)),0);
    /* Negated coefficients with positive sense describe the same interior. */
    int reverse=alea_quadric_surface(sys,0,-4,-1,-16,0,0,0,80,-40,960,-15184);
    ASSERT_EQ(alea_tighten_bbox_mixed_constraints(sys,alea_halfspace(sys,reverse,1),&box),0);
    ASSERT(box.min_x<=8 && box.max_z>=31);
    alea_destroy(sys);
}

TEST(mixed_quadric_does_not_erase_negative_eigenvalue) {
    alea_system_t* sys=alea_create();ASSERT_NOT_NULL(sys);
    /* x^2+y^2-1e-12*z^2 <= 1 is unbounded, however close it looks to a cylinder. */
    int q=alea_quadric_surface(sys,0,1,1,-1e-12,0,0,0,0,0,0,-1);
    alea_node_id_t interior=alea_halfspace(sys,q,-1);
    alea_bbox_t box={1,2,3,4,5,6},before=box;
    ASSERT_EQ(alea_tighten_bbox_mixed_constraints(sys,interior,&box),-1);
    ASSERT_EQ(memcmp(&box,&before,sizeof(box)),0);
    /* A slab in x still leaves the z recession direction open. */
    int xl=alea_plane_surface(sys,0,1,0,0,2),xh=alea_plane_surface(sys,0,1,0,0,-2);
    alea_node_id_t wrong=alea_intersection(sys,interior,
        alea_intersection(sys,alea_halfspace(sys,xl,1),alea_halfspace(sys,xh,-1)));
    ASSERT_EQ(alea_tighten_bbox_mixed_constraints(sys,wrong,&box),-1);
    int zl=alea_plane_surface(sys,0,0,0,1,100),zh=alea_plane_surface(sys,0,0,0,1,-100);
    alea_node_id_t capped=alea_intersection(sys,interior,
        alea_intersection(sys,alea_halfspace(sys,zl,1),alea_halfspace(sys,zh,-1)));
    ASSERT_EQ(alea_tighten_bbox_mixed_constraints(sys,capped,&box),0);
    for(int i=-100;i<=100;i++) {
        double r=sqrt(1+1e-12*i*i);
        ASSERT(alea_bbox_contains_point(&box,r,0,i));
        ASSERT(alea_bbox_contains_point(&box,-r,0,i));
    }
    /* Constraints inside a union must not be used to restrict another branch. */
    alea_node_id_t unbounded=alea_union(sys,capped,interior);
    box=before;
    ASSERT_EQ(alea_tighten_bbox_mixed_constraints(sys,unbounded,&box),-1);
    ASSERT_EQ(memcmp(&box,&before,sizeof(box)),0);
    alea_destroy(sys);
}

TEST(mixed_quadric_translated_oblique_slab_samples) {
    for(int translated=0;translated<2;translated++) {
        alea_system_t* sys=alea_create();ASSERT_NOT_NULL(sys);
        double center[3]={translated?1e6:100,translated?-2e6:-200,translated?3e6:300};
        /* Mutually orthogonal integer directions of length three. */
        double u[3]={1,2,2},v[3]={2,1,-2},w[3]={-2,2,-1};
        double uc=0,vc=0,wc=0;
        for(int k=0;k<3;k++){uc+=u[k]*center[k];vc+=v[k]*center[k];wc+=w[k]*center[k];}
        int q=alea_quadric_surface(sys,0,5,5,8,8,4,-4,
            -2*(u[0]*uc+v[0]*vc),-2*(u[1]*uc+v[1]*vc),-2*(u[2]*uc+v[2]*vc),
            uc*uc+vc*vc-16);
        int low=alea_plane_surface(sys,0,w[0],w[1],w[2],-wc+6);
        int high=alea_plane_surface(sys,0,w[0],w[1],w[2],-wc-6);
        alea_node_id_t root=alea_intersection(sys,alea_halfspace(sys,q,-1),
            alea_intersection(sys,alea_halfspace(sys,low,1),alea_halfspace(sys,high,-1)));
        alea_bbox_t box;
        ASSERT_EQ(alea_tighten_bbox_mixed_constraints(sys,root,&box),0);
        ASSERT(box.max_x-box.min_x<20);
        for(int axial=-6;axial<=6;axial++)for(int angle=0;angle<360;angle++) {
            double a=angle*6.283185307179586/360,p[3];
            for(int k=0;k<3;k++)p[k]=center[k]+(4*cos(a)*u[k]+4*sin(a)*v[k]+axial*w[k])/9;
            ASSERT(alea_bbox_contains_point(&box,p[0],p[1],p[2]));
        }
        alea_destroy(sys);
    }
}

TEST(mixed_quadric_uncertain_arithmetic_preserves_output) {
    alea_system_t* sys=alea_create();ASSERT_NOT_NULL(sys);
    int q=alea_quadric_surface(sys,0,1,1,2,-2,0,0,0,0,0,-8);
    int low=alea_plane_surface(sys,0,1,1,0,0);
    int high=alea_plane_surface(sys,0,1+1e-9,1,0,-10);
    alea_node_id_t root=alea_intersection(sys,alea_halfspace(sys,q,-1),
        alea_intersection(sys,alea_halfspace(sys,low,1),alea_halfspace(sys,high,-1)));
    alea_bbox_t box={1,2,3,4,5,6},before=box;
    ASSERT_EQ(alea_tighten_bbox_mixed_constraints(sys,root,&box),-1);
    ASSERT_EQ(memcmp(&box,&before,sizeof(box)),0);
    alea_interrupt();
    ASSERT_EQ(alea_tighten_bbox_mixed_constraints(sys,root,&box),-1);
    alea_clear_interrupt();
    ASSERT_EQ(memcmp(&box,&before,sizeof(box)),0);
    alea_destroy(sys);
}

TEST_MAIN()
