// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_raycast.c - Unit tests for ray casting module
 */

#include "alea_test.h"
#include "alea.h"
#include "raycast/raycast.h"
#include "raycast/ray_intersect.h"
#include "core/alea_system.h"

#define EPS 1e-6

/* ------------------------------------------------------------------------- */
/* Ray-Primitive Intersection Tests                                           */
/* ------------------------------------------------------------------------- */

TEST(ray_sphere_intersection) {
    alea_ray_t ray;
    alea_sphere_data_t sphere = {0, 0, 0, 5.0};
    double t[2];
    int count;

    /* Ray from outside, hitting sphere */
    alea_ray_init(&ray, -10, 0, 0, 1, 0, 0);
    count = ray_intersect_sphere(&ray, &sphere, t);
    ASSERT_EQ(count, 2);
    ASSERT_NEAR(t[0], 5.0, EPS);   /* Enter at x=-5 */
    ASSERT_NEAR(t[1], 15.0, EPS);  /* Exit at x=5 */

    /* Ray from inside sphere */
    alea_ray_init(&ray, 0, 0, 0, 1, 0, 0);
    count = ray_intersect_sphere(&ray, &sphere, t);
    ASSERT_EQ(count, 2);
    ASSERT_NEAR(t[0], -5.0, EPS);  /* Behind origin */
    ASSERT_NEAR(t[1], 5.0, EPS);   /* Forward */

    /* Ray missing sphere */
    alea_ray_init(&ray, -10, 10, 0, 1, 0, 0);
    count = ray_intersect_sphere(&ray, &sphere, t);
    ASSERT_EQ(count, 0);

    /* Ray tangent to sphere (grazing) */
    alea_ray_init(&ray, -10, 5, 0, 1, 0, 0);
    count = ray_intersect_sphere(&ray, &sphere, t);
    ASSERT_EQ(count, 2);
    ASSERT_NEAR(t[0], t[1], EPS);  /* Single touch point */
}

TEST(ray_plane_intersection) {
    alea_ray_t ray;
    alea_plane_data_t plane;
    double t[2], nx, ny, nz;
    int count;

    /* Plane at z=5 (normal pointing +z) */
    plane.a = 0; plane.b = 0; plane.c = 1; plane.d = -5;

    /* Ray pointing towards plane */
    alea_ray_init(&ray, 0, 0, 0, 0, 0, 1);
    count = ray_intersect_plane(&ray, &plane, t, &nx, &ny, &nz);
    ASSERT_EQ(count, 1);
    ASSERT_NEAR(t[0], 5.0, EPS);
    ASSERT_NEAR(nz, 1.0, EPS);

    /* Ray parallel to plane */
    alea_ray_init(&ray, 0, 0, 0, 1, 0, 0);
    count = ray_intersect_plane(&ray, &plane, t, NULL, NULL, NULL);
    ASSERT_EQ(count, 0);

    /* Ray pointing away from plane */
    alea_ray_init(&ray, 0, 0, 0, 0, 0, -1);
    count = ray_intersect_plane(&ray, &plane, t, NULL, NULL, NULL);
    ASSERT_EQ(count, 1);
    ASSERT_NEAR(t[0], -5.0, EPS);  /* Intersection behind origin */
}

TEST(ray_cylinder_z_intersection) {
    alea_ray_t ray;
    alea_cylinder_z_data_t cyl = {0, 0, 3.0};
    double t[2];
    int count;

    /* Ray perpendicular to cylinder axis */
    alea_ray_init(&ray, -10, 0, 5, 1, 0, 0);
    count = ray_intersect_cylinder_z(&ray, &cyl, t);
    ASSERT_EQ(count, 2);
    ASSERT_NEAR(t[0], 7.0, EPS);   /* Enter at x=-3 */
    ASSERT_NEAR(t[1], 13.0, EPS);  /* Exit at x=3 */

    /* Ray parallel to cylinder axis (inside) */
    alea_ray_init(&ray, 1, 1, 0, 0, 0, 1);
    count = ray_intersect_cylinder_z(&ray, &cyl, t);
    ASSERT_EQ(count, 0);  /* No intersection with infinite cylinder wall */

    /* Ray missing cylinder */
    alea_ray_init(&ray, -10, 5, 0, 1, 0, 0);
    count = ray_intersect_cylinder_z(&ray, &cyl, t);
    ASSERT_EQ(count, 0);
}

TEST(ray_box_intersection) {
    alea_ray_t ray;
    alea_box_data_t box = {-1, 1, -2, 2, -3, 3};
    double t[2];
    int count;

    /* Ray through center */
    alea_ray_init(&ray, -10, 0, 0, 1, 0, 0);
    count = ray_intersect_box(&ray, &box, t);
    ASSERT_EQ(count, 2);
    ASSERT_NEAR(t[0], 9.0, EPS);   /* Enter at x=-1 */
    ASSERT_NEAR(t[1], 11.0, EPS);  /* Exit at x=1 */

    /* Ray from inside */
    alea_ray_init(&ray, 0, 0, 0, 1, 0, 0);
    count = ray_intersect_box(&ray, &box, t);
    ASSERT_EQ(count, 2);
    ASSERT_NEAR(t[0], -1.0, EPS);  /* Behind */
    ASSERT_NEAR(t[1], 1.0, EPS);   /* Forward */

    /* Ray missing box */
    alea_ray_init(&ray, -10, 10, 0, 1, 0, 0);
    count = ray_intersect_box(&ray, &box, t);
    ASSERT_EQ(count, 0);

    /* Diagonal ray */
    alea_ray_init(&ray, -10, -10, -10, 1, 1, 1);
    count = ray_intersect_box(&ray, &box, t);
    ASSERT_EQ(count, 2);
    ASSERT(t[0] < t[1]);
}

TEST(ray_cone_z_intersection) {
    alea_ray_t ray;
    alea_cone_z_data_t cone = {0, 0, 0, 1.0, 0};  /* 45-degree cone */
    double t[2];
    int count;

    /* Ray perpendicular to axis */
    alea_ray_init(&ray, -10, 0, 5, 1, 0, 0);  /* At z=5, cone radius is 5 */
    count = ray_intersect_cone_z(&ray, &cone, t);
    ASSERT_EQ(count, 2);
    ASSERT_NEAR(t[0], 5.0, EPS);   /* Enter at x=-5 */
    ASSERT_NEAR(t[1], 15.0, EPS);  /* Exit at x=5 */
}

TEST(ray_rcc_intersection) {
    alea_ray_t ray;
    /* RCC: base at (0,0,0), height vector (0,0,10), radius 2 */
    alea_rcc_data_t rcc = {0, 0, 0, 0, 0, 10, 2.0};
    double t[2];
    int count;

    /* Ray through side */
    alea_ray_init(&ray, -5, 0, 5, 1, 0, 0);
    count = ray_intersect_rcc(&ray, &rcc, t);
    ASSERT_EQ(count, 2);
    ASSERT_NEAR(t[0], 3.0, EPS);   /* Enter at x=-2 */
    ASSERT_NEAR(t[1], 7.0, EPS);   /* Exit at x=2 */

    /* Ray through caps */
    alea_ray_init(&ray, 0, 0, -5, 0, 0, 1);
    count = ray_intersect_rcc(&ray, &rcc, t);
    ASSERT_EQ(count, 2);
    ASSERT_NEAR(t[0], 5.0, EPS);   /* Enter bottom cap */
    ASSERT_NEAR(t[1], 15.0, EPS);  /* Exit top cap */

    /* Ray missing */
    alea_ray_init(&ray, -5, 5, 5, 1, 0, 0);
    count = ray_intersect_rcc(&ray, &rcc, t);
    ASSERT_EQ(count, 0);
}

/* ------------------------------------------------------------------------- */
/* Ray Utilities Tests                                                        */
/* ------------------------------------------------------------------------- */

TEST(ray_normalization) {
    alea_ray_t ray;

    /* Non-unit direction should be normalized */
    alea_ray_init(&ray, 0, 0, 0, 3, 4, 0);
    ASSERT_NEAR(ray.dx, 0.6, EPS);
    ASSERT_NEAR(ray.dy, 0.8, EPS);
    ASSERT_NEAR(ray.dz, 0.0, EPS);

    /* Check magnitude is 1 */
    double mag = sqrt(ray.dx*ray.dx + ray.dy*ray.dy + ray.dz*ray.dz);
    ASSERT_NEAR(mag, 1.0, EPS);
}

/* ------------------------------------------------------------------------- */
/* Full Raycast Tests                                                         */
/* ------------------------------------------------------------------------- */

TEST(raycast_simple_geometry) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Create sphere surface and get interior node */
    int surf_idx = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    alea_node_id_t sphere = alea_surface_at(sys, surf_idx)->neg_node;

    /* Add cell */
    int m1 = alea_add_material(sys, 1);

    alea_add_cell(sys, 1, sphere, m1, -2.7, 0);

    /* Cast ray through sphere */
    alea_raycast_result_t result;
    alea_raycast_result_init(&result);

    int rc = alea_raycast(sys, -10, 0, 0, 1, 0, 0, 100, &result);
    ASSERT_EQ(rc, 0);

    /* Should have 2 hits (enter and exit) */
    ASSERT_EQ(result.hit_count, 2);
    ASSERT_NEAR(result.hits[0].t, 5.0, EPS);
    ASSERT_NEAR(result.hits[1].t, 15.0, EPS);

    /* Should have segments */
    ASSERT(result.segment_count >= 1);

    alea_raycast_result_free(&result);
    alea_destroy(sys);
}

TEST(raycast_path_length) {
    alea_raycast_result_t result;
    alea_raycast_result_init(&result);

    /* Manually create some segments */
    alea_ray_segment_t seg1 = {0, 5, 1, 1, -2.7, -1};
    alea_ray_segment_t seg2 = {5, 10, 2, 2, -8.0, -1};
    alea_ray_segment_t seg3 = {10, 15, 3, 1, -2.7, -1};

    result.segments = malloc(3 * sizeof(alea_ray_segment_t));
    ASSERT_NOT_NULL(result.segments);
    result.segments[0] = seg1;
    result.segments[1] = seg2;
    result.segments[2] = seg3;
    result.segment_count = 3;
    result.segment_capacity = 3;

    /* Total path length through material 1 */
    double len1 = alea_raycast_path_length(&result, 1);
    ASSERT_NEAR(len1, 10.0, EPS);  /* 5 + 5 */

    /* Total path length through material 2 */
    double len2 = alea_raycast_path_length(&result, 2);
    ASSERT_NEAR(len2, 5.0, EPS);

    /* Total path length through all materials */
    double total = alea_raycast_path_length(&result, -1);
    ASSERT_NEAR(total, 15.0, EPS);

    alea_raycast_result_free(&result);
}

TEST_MAIN()
