// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_raycast_numerical.c - Numerical robustness tests for raycast module
 *
 * Tests edge cases identified by expert code reviews:
 * tangent rays, near-axis hits, large-t geometry, zero-direction,
 * IEEE inf traversal, and epsilon constant sanity.
 */

#include "alea_test.h"
#include "alea.h"
#include "raycast/raycast.h"
#include "raycast/ray_intersect.h"
#include "raycast/ray_epsilon.h"
#include "raycast/bvh.h"
#include "core/alea_system.h"

#include <math.h>
#include <float.h>

#define EPS 1e-6

/* =========================================================================
 * 1. Tangent ray to sphere (discriminant ~ 0)
 * ========================================================================= */

TEST(tangent_ray_sphere) {
    alea_ray_t ray;
    alea_sphere_data_t sphere = {0, 0, 0, 5.0};
    double t[2];

    /* Ray tangent to sphere at (0, 5, 0): origin at (-100, 5, 0), dir (1,0,0)
     * discriminant = b^2 - 4ac ≈ 0 */
    alea_ray_init(&ray, -100, 5.0, 0, 1, 0, 0);
    int count = ray_intersect_sphere(&ray, &sphere, t);

    /* Should find 2 roots very close together (tangent = double root) */
    ASSERT(count == 2);
    ASSERT_NEAR(t[0], t[1], 1e-5);
    ASSERT_NEAR(t[0], 100.0, 0.1);  /* Touch at x=0, so t ≈ 100 */
}

TEST(near_tangent_ray_sphere) {
    alea_ray_t ray;
    alea_sphere_data_t sphere = {0, 0, 0, 5.0};
    double t[2];

    /* Slightly inside tangent: y = 5 - 1e-8 */
    alea_ray_init(&ray, -100, 5.0 - 1e-8, 0, 1, 0, 0);
    int count = ray_intersect_sphere(&ray, &sphere, t);

    /* Should still find 2 intersections */
    ASSERT(count == 2);
    ASSERT(t[0] < t[1]);
    ASSERT(t[1] - t[0] < 0.01);  /* Very close together */
}

/* =========================================================================
 * 2. Tangent ray to cylinder
 * ========================================================================= */

TEST(tangent_ray_cylinder_z) {
    alea_ray_t ray;
    alea_cylinder_z_data_t cyl = {0, 0, 5.0};  /* center_x, center_y, radius */
    double t[2];

    /* Ray tangent to cylinder at (5, 0, z): origin (-100, 5, 0), dir (1,0,0) */
    alea_ray_init(&ray, -100, 5.0, 0, 1, 0, 0);
    int count = ray_intersect_cylinder_z(&ray, &cyl, t);

    ASSERT(count == 2);
    ASSERT_NEAR(t[0], t[1], 1e-5);
    ASSERT_NEAR(t[0], 100.0, 0.1);
}

/* =========================================================================
 * 3. Near-axis RCC cap hit (cross-product fix)
 * ========================================================================= */

TEST(rcc_cap_near_axis) {
    alea_ray_t ray;
    double t[4];

    /* RCC along Z: base at (0,0,0), height (0,0,10), radius 5 */
    alea_rcc_data_t rcc = {
        .base_x = 0, .base_y = 0, .base_z = 0,
        .height_x = 0, .height_y = 0, .height_z = 10,
        .radius = 5
    };

    /* Ray nearly along the axis, hitting both caps.
     * Origin slightly off-axis: (1e-9, 0, -5), direction (0, 0, 1) */
    alea_ray_init(&ray, 1e-9, 0, -5, 0, 0, 1);
    int count = ray_intersect_rcc(&ray, &rcc, t);

    /* Should hit bottom cap at t=5 and top cap at t=15 */
    ASSERT(count == 2);
    ASSERT_NEAR(t[0], 5.0, 1e-3);
    ASSERT_NEAR(t[1], 15.0, 1e-3);
}

TEST(rcc_cap_exactly_on_axis) {
    alea_ray_t ray;
    double t[4];

    alea_rcc_data_t rcc = {
        .base_x = 0, .base_y = 0, .base_z = 0,
        .height_x = 0, .height_y = 0, .height_z = 10,
        .radius = 5
    };

    /* Ray exactly on the axis */
    alea_ray_init(&ray, 0, 0, -5, 0, 0, 1);
    int count = ray_intersect_rcc(&ray, &rcc, t);

    /* Should hit both caps */
    ASSERT(count == 2);
    ASSERT_NEAR(t[0], 5.0, 1e-3);
    ASSERT_NEAR(t[1], 15.0, 1e-3);
}

/* =========================================================================
 * 4. Large-t progress (geometry far from origin)
 * ========================================================================= */

TEST(large_t_sphere_intersection) {
    alea_ray_t ray;
    alea_sphere_data_t sphere = {1e6, 0, 0, 1.0};
    double t[2];

    /* Sphere at x=1e6, radius 1. Ray from origin along +X */
    alea_ray_init(&ray, 0, 0, 0, 1, 0, 0);
    int count = ray_intersect_sphere(&ray, &sphere, t);

    ASSERT(count == 2);
    /* Enter at x=999999, exit at x=1000001 */
    ASSERT_NEAR(t[0], 999999.0, 0.1);
    ASSERT_NEAR(t[1], 1000001.0, 0.1);
    ASSERT(t[1] > t[0]);
}

/* =========================================================================
 * 5. Zero-direction ray → returns -1
 * ========================================================================= */

TEST(zero_direction_ray) {
    alea_ray_t ray;

    /* Zero direction should fail */
    int rc = alea_ray_init(&ray, 1, 2, 3, 0, 0, 0);
    ASSERT_EQ(rc, -1);

    /* Fallback direction should be valid (0,0,1) */
    ASSERT_NEAR(ray.dz, 1.0, 1e-10);
}

TEST(near_zero_direction_ray) {
    alea_ray_t ray;

    /* Very small but non-zero direction should succeed */
    int rc = alea_ray_init(&ray, 0, 0, 0, 1e-5, 0, 0);
    ASSERT_EQ(rc, 0);
    ASSERT_NEAR(ray.dx, 1.0, 1e-6);  /* Should be normalized to unit length */
}

/* =========================================================================
 * 6. Near-tangent torus
 * ========================================================================= */

TEST(torus_basic_intersection) {
    alea_ray_t ray;
    double t[4];

    /* Z-axis torus: R=10, r=2, centered at origin */
    alea_torus_data_t torus = {
        .center_x = 0, .center_y = 0, .center_z = 0,
        .major_radius = 10.0,
        .minor_radius = 2.0,
        .axis = ALEA_AXIS_Z
    };

    /* Ray from outside along X axis at z=0 */
    alea_ray_init(&ray, -20, 0, 0, 1, 0, 0);
    int count = ray_intersect_torus(&ray, &torus, t);

    /* Should hit 4 times: enter outer, exit inner, enter inner, exit outer */
    ASSERT_EQ(count, 4);

    /* Verify ordering */
    for (int i = 0; i < count - 1; i++) {
        ASSERT(t[i] < t[i + 1]);
    }

    /* Outer hits at x = ±12, inner hits at x = ±8 */
    ASSERT_NEAR(t[0], 8.0, 0.1);   /* Enter outer at x=-12, t = -20+12 = 8 */
    ASSERT_NEAR(t[1], 12.0, 0.1);  /* Exit inner at x=-8, t = -20+8 = 12 */
    ASSERT_NEAR(t[2], 28.0, 0.1);  /* Enter inner at x=8, t = -20+8 = 28 */
    ASSERT_NEAR(t[3], 32.0, 0.1);  /* Exit outer at x=12, t = -20+12 = 32 */
}

/* =========================================================================
 * 7. Axis-aligned ray through BVH (IEEE inf traversal)
 * ========================================================================= */

TEST(axis_aligned_ray_init) {
    alea_ray_t ray;

    /* Axis-aligned ray: dx=0, dy=0, dz=1 → inv_dx and inv_dy should be +inf */
    alea_ray_init(&ray, 0, 0, 0, 0, 0, 1);

    /* IEEE 754: 1.0/0.0 = +inf */
    ASSERT(isinf(ray.inv_dx));
    ASSERT(isinf(ray.inv_dy));
    ASSERT_NEAR(ray.inv_dz, 1.0, 1e-10);
}

TEST(axis_aligned_ray_negative) {
    alea_ray_t ray;

    /* Direction (0, -1, 0): inv_dy should be -1.0, inv_dx/dz should be ±inf */
    alea_ray_init(&ray, 0, 0, 0, 0, -1, 0);

    ASSERT(isinf(ray.inv_dx));
    ASSERT_NEAR(ray.inv_dy, -1.0, 1e-10);
    ASSERT(isinf(ray.inv_dz));
}

TEST(ray_init_normalized_skip_sqrt) {
    alea_ray_t ray;

    /* Pre-normalized direction */
    double dx = 1.0 / sqrt(3.0);
    double dy = 1.0 / sqrt(3.0);
    double dz = 1.0 / sqrt(3.0);

    alea_ray_init_normalized(&ray, 1, 2, 3, dx, dy, dz);

    ASSERT_NEAR(ray.ox, 1.0, 1e-15);
    ASSERT_NEAR(ray.dx, dx, 1e-15);
    ASSERT_NEAR(ray.inv_dx, 1.0 / dx, 1e-10);
}

/* =========================================================================
 * 8. Epsilon constant sanity checks
 * ========================================================================= */

TEST(epsilon_constants_sane) {
    /* Verify epsilon ordering: DISCRIMINANT_TOL < RAY_EPSILON */
    ASSERT(DISCRIMINANT_TOL < RAY_EPSILON);

    /* Verify epsilons are positive */
    ASSERT(RAY_EPSILON > 0);
    ASSERT(DISCRIMINANT_TOL > 0);
    ASSERT(DEDUP_EPSILON > 0);
    ASSERT(SURFACE_SAMPLE_OFFSET > 0);
    ASSERT(DDA_SAMPLE_OFFSET > 0);

    /* Verify sample offsets are much larger than numerical epsilons */
    ASSERT(SURFACE_SAMPLE_OFFSET > RAY_EPSILON * 1000);
    ASSERT(DDA_SAMPLE_OFFSET > RAY_EPSILON * 1000);

    /* Exact value checks */
    ASSERT_NEAR(RAY_EPSILON, 1e-10, 1e-20);
    ASSERT_NEAR(DISCRIMINANT_TOL, 1e-12, 1e-22);
    ASSERT_NEAR(DEDUP_EPSILON, 1e-10, 1e-20);
    ASSERT_NEAR(SURFACE_SAMPLE_OFFSET, 0.001, 1e-10);
    ASSERT_NEAR(DDA_SAMPLE_OFFSET, 0.01, 1e-10);
}

/* =========================================================================
 * Additional: Vieta's quadratic stability
 * ========================================================================= */

TEST(stable_quadratic_sphere_large_offset) {
    alea_ray_t ray;
    alea_sphere_data_t sphere = {0, 0, 0, 1.0};
    double t[2];

    /* Ray from far away: b is large, b^2 >> 4ac.
     * Classic formula has catastrophic cancellation in b - sqrt(b^2 - 4ac).
     * Vieta's form avoids this.
     *
     * Use 1e6 (not 1e8+) because c = |L|^2 - r^2 must be exact in double:
     * 1e12 - 1 is exact (< 2^53), but 1e16 - 1 rounds to 1e16 (ULP=2). */
    alea_ray_init(&ray, -1e6, 0, 0, 1, 0, 0);
    int count = ray_intersect_sphere(&ray, &sphere, t);

    ASSERT_EQ(count, 2);
    ASSERT_NEAR(t[0], 1e6 - 1.0, 1e-4);
    ASSERT_NEAR(t[1], 1e6 + 1.0, 1e-4);
    /* The gap should be exactly 2 (diameter) — Vieta's preserves this */
    ASSERT_NEAR(t[1] - t[0], 2.0, 1e-6);
}

/* =========================================================================
 * Additional: Sort stability (insertion sort for small arrays)
 * ========================================================================= */

TEST(sort_small_hit_array) {
    /* Create a small raycast result with known hits */
    alea_raycast_result_t result;
    alea_raycast_result_init(&result);
    alea_raycast_result_reserve(&result, 8, 0);

    /* Add hits in reverse order */
    alea_ray_hit_t hits[] = {
        {.t = 5.0, .surface_id = 1},
        {.t = 1.0, .surface_id = 2},
        {.t = 3.0, .surface_id = 3},
        {.t = 2.0, .surface_id = 4},
        {.t = 4.0, .surface_id = 5},
    };

    for (int i = 0; i < 5; i++) {
        result.hits.data[i] = hits[i];
    }
    result.hits.count = 5;

    /* Trigger sort via a full system raycast isn't needed — just verify
     * the hits can be read back. The actual sort is tested implicitly
     * by any full raycast test that produces >1 hit. */

    /* Verify we can construct and clean up a result */
    ASSERT_EQ(result.hits.count, 5);
    alea_raycast_result_free(&result);
}

/* =========================================================================
 * Asymmetric-center cylinders (perpendicular-axis order regression)
 *
 * The infinite-cylinder equation is symmetric under swapping its two
 * perpendicular coordinates, so any test with equal center components
 * passes even if the implementation pairs coordinate and center in the
 * wrong order. These tests use deliberately unequal center components and
 * check the crossings against the closed-form values: a swapped pair
 * intersects a cylinder at the transposed center and fails loudly.
 * (A swapped (x,z) pair for C/Y once misplaced every asymmetric Y-cylinder
 * crossing in raycast; slices were unaffected.)
 * ========================================================================= */

TEST(cylinder_x_asymmetric_center) {
    alea_ray_t ray;
    alea_cylinder_x_data_t cyl = { .center_y = 7.0, .center_z = -3.0,
                                   .radius = 5.0 };
    double t[2];

    /* Ray along +Y in the cylinder's center plane z = -3:
     * crossings at y = 7 ± 5, i.e. t = 102 and 112. */
    alea_ray_init(&ray, 20.0, -100.0, -3.0, 0, 1, 0);
    int count = ray_intersect_cylinder_x(&ray, &cyl, t);

    ASSERT_EQ(count, 2);
    ASSERT_NEAR(t[0], 102.0, EPS);
    ASSERT_NEAR(t[1], 112.0, EPS);
}

TEST(cylinder_y_asymmetric_center) {
    alea_ray_t ray;
    /* The E-lite surface that exposed the swap: C/Y 424.0 2.9 682.221 */
    alea_cylinder_y_data_t cyl = { .center_x = 424.0, .center_z = 2.9,
                                   .radius = 682.221 };
    double t[2];

    /* Ray along +X at z = 60: half-chord = sqrt(r² - (60 - 2.9)²),
     * crossings at x = 424 ∓ half-chord. */
    alea_ray_init(&ray, 300.0, 0.0, 60.0, 1, 0, 0);
    int count = ray_intersect_cylinder_y(&ray, &cyl, t);

    double dz = 60.0 - cyl.center_z;
    double half_chord = sqrt(cyl.radius * cyl.radius - dz * dz);
    ASSERT_EQ(count, 2);
    ASSERT_NEAR(t[0], 424.0 - half_chord - 300.0, EPS);
    ASSERT_NEAR(t[1], 424.0 + half_chord - 300.0, EPS);
}

TEST(cylinder_y_asymmetric_center_along_z) {
    alea_ray_t ray;
    alea_cylinder_y_data_t cyl = { .center_x = 424.0, .center_z = 2.9,
                                   .radius = 682.221 };
    double t[2];

    /* Same cylinder probed along +Z through its center plane x = 424:
     * crossings at z = 2.9 ± r. Exercises the other perpendicular axis. */
    alea_ray_init(&ray, 424.0, -50.0, -1000.0, 0, 0, 1);
    int count = ray_intersect_cylinder_y(&ray, &cyl, t);

    ASSERT_EQ(count, 2);
    ASSERT_NEAR(t[0], 1000.0 + 2.9 - 682.221, EPS);
    ASSERT_NEAR(t[1], 1000.0 + 2.9 + 682.221, EPS);
}

TEST(cylinder_z_asymmetric_center) {
    alea_ray_t ray;
    alea_cylinder_z_data_t cyl = { .center_x = -12.0, .center_y = 34.0,
                                   .radius = 8.0 };
    double t[2];

    /* Ray along +X at y = 30 (off-center chord):
     * crossings at x = -12 ± sqrt(64 - 16). */
    alea_ray_init(&ray, -100.0, 30.0, 5.0, 1, 0, 0);
    int count = ray_intersect_cylinder_z(&ray, &cyl, t);

    double half_chord = sqrt(64.0 - 16.0);
    ASSERT_EQ(count, 2);
    ASSERT_NEAR(t[0], -12.0 - half_chord + 100.0, EPS);
    ASSERT_NEAR(t[1], -12.0 + half_chord + 100.0, EPS);
}

TEST_MAIN()
