// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_slice_numerical.c - Numerical robustness tests for slice module
 *
 * Tests edge cases from the Fabio Luporini-style code review:
 * stable quadratic, cylinder/cone consolidation, TRC dispatch,
 * ellipse caching, box polygon ordering, parallel lines, scanline
 * intersection accuracy.
 *
 * Uses the low-level alea_intersect_primitive_plane API to test
 * intersection math directly without needing full CSG systems.
 */

#include "alea_test.h"
#include "alea.h"
#include "alea_slice.h"
#include "slice/curve_intersect.h"
#include "core/alea_system.h"

#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif

/* =========================================================================
 * Helper: create an axis-aligned slice plane
 * ========================================================================= */

static void make_axis_plane(alea_slice_plane_t* p, int axis, double value) {
    alea_slice_plane_init_axis(p, axis, value);
}

/* =========================================================================
 * 1. Stable quadratic: cylinder far from origin
 *
 * A cylinder at x=1e6 sliced perpendicularly should produce a valid circle.
 * Classical formula has catastrophic cancellation; Vieta's form avoids it.
 * ========================================================================= */

TEST(stable_quadratic_far_cylinder) {
    alea_primitive_data_t data;
    memset(&data, 0, sizeof(data));
    data.cyl_z.center_x = 1e6;
    data.cyl_z.center_y = 0;
    data.cyl_z.radius = 1.0;

    alea_slice_plane_t plane;
    make_axis_plane(&plane, 2, 0.0);

    alea_curve_2d_t curve;
    bool hit = alea_intersect_primitive_plane(ALEA_PRIMITIVE_CYLINDER_Z, &data, &plane, &curve);
    ASSERT(hit);
    ASSERT_EQ(curve.type, ALEA_CURVE_CIRCLE);
    ASSERT_NEAR(curve.data.circle.radius, 1.0, 1e-6);
    /* Center should be at (1e6, 0) in plane coords */
    ASSERT_NEAR(curve.data.circle.center[0], 1e6, 0.1);
}

/* =========================================================================
 * 2. Cylinder consolidation: all three axes produce correct circles
 * ========================================================================= */

TEST(cylinder_x_perpendicular_circle) {
    alea_primitive_data_t data;
    memset(&data, 0, sizeof(data));
    data.cyl_x.center_y = 2.0;
    data.cyl_x.center_z = 3.0;
    data.cyl_x.radius = 4.0;

    alea_slice_plane_t plane;
    make_axis_plane(&plane, 0, 0.0);  /* YZ plane */

    alea_curve_2d_t curve;
    bool hit = alea_intersect_primitive_plane(ALEA_PRIMITIVE_CYLINDER_X, &data, &plane, &curve);
    ASSERT(hit);
    ASSERT_EQ(curve.type, ALEA_CURVE_CIRCLE);
    ASSERT_NEAR(curve.data.circle.radius, 4.0, 1e-6);
}

TEST(cylinder_y_perpendicular_circle) {
    alea_primitive_data_t data;
    memset(&data, 0, sizeof(data));
    data.cyl_y.center_x = 1.0;
    data.cyl_y.center_z = -2.0;
    data.cyl_y.radius = 3.0;

    alea_slice_plane_t plane;
    make_axis_plane(&plane, 1, 0.0);  /* XZ plane */

    alea_curve_2d_t curve;
    bool hit = alea_intersect_primitive_plane(ALEA_PRIMITIVE_CYLINDER_Y, &data, &plane, &curve);
    ASSERT(hit);
    ASSERT_EQ(curve.type, ALEA_CURVE_CIRCLE);
    ASSERT_NEAR(curve.data.circle.radius, 3.0, 1e-6);
}

TEST(cylinder_z_perpendicular_circle) {
    alea_primitive_data_t data;
    memset(&data, 0, sizeof(data));
    data.cyl_z.center_x = 0;
    data.cyl_z.center_y = 0;
    data.cyl_z.radius = 5.0;

    alea_slice_plane_t plane;
    make_axis_plane(&plane, 2, 0.0);  /* XY plane */

    alea_curve_2d_t curve;
    bool hit = alea_intersect_primitive_plane(ALEA_PRIMITIVE_CYLINDER_Z, &data, &plane, &curve);
    ASSERT(hit);
    ASSERT_EQ(curve.type, ALEA_CURVE_CIRCLE);
    ASSERT_NEAR(curve.data.circle.radius, 5.0, 1e-6);
}

/* =========================================================================
 * 3. Cone consolidation: all three axes produce valid conics
 * ========================================================================= */

TEST(cone_z_perpendicular_circle) {
    alea_primitive_data_t data;
    memset(&data, 0, sizeof(data));
    data.cone_z.apex_x = 0;
    data.cone_z.apex_y = 0;
    data.cone_z.apex_z = 0;
    data.cone_z.tan_angle_sq = 1.0;  /* 45° half-angle */

    alea_slice_plane_t plane;
    make_axis_plane(&plane, 2, 3.0);  /* z=3 → circle of radius 3 */

    alea_curve_2d_t curve;
    bool hit = alea_intersect_primitive_plane(ALEA_PRIMITIVE_CONE_Z, &data, &plane, &curve);
    ASSERT(hit);

    /* Perpendicular cone cut → circle or degenerate ellipse (semi_a ≈ semi_b) */
    ASSERT(curve.type == ALEA_CURVE_CIRCLE || curve.type == ALEA_CURVE_ELLIPSE);
    if (curve.type == ALEA_CURVE_CIRCLE) {
        ASSERT_NEAR(curve.data.circle.radius, 3.0, 1e-6);
    } else {
        ASSERT_NEAR(curve.data.ellipse.semi_a, 3.0, 1e-4);
        ASSERT_NEAR(curve.data.ellipse.semi_b, 3.0, 1e-4);
    }
}

TEST(cone_x_perpendicular_circle) {
    alea_primitive_data_t data;
    memset(&data, 0, sizeof(data));
    data.cone_x.apex_x = 0;
    data.cone_x.apex_y = 0;
    data.cone_x.apex_z = 0;
    data.cone_x.tan_angle_sq = 0.25;  /* tan(θ) = 0.5 */

    alea_slice_plane_t plane;
    make_axis_plane(&plane, 0, 4.0);  /* x=4 → radius = 4*0.5 = 2 */

    alea_curve_2d_t curve;
    bool hit = alea_intersect_primitive_plane(ALEA_PRIMITIVE_CONE_X, &data, &plane, &curve);
    ASSERT(hit);

    ASSERT(curve.type == ALEA_CURVE_CIRCLE || curve.type == ALEA_CURVE_ELLIPSE);
    if (curve.type == ALEA_CURVE_CIRCLE) {
        ASSERT_NEAR(curve.data.circle.radius, 2.0, 1e-6);
    } else {
        ASSERT_NEAR(curve.data.ellipse.semi_a, 2.0, 1e-4);
        ASSERT_NEAR(curve.data.ellipse.semi_b, 2.0, 1e-4);
    }
}

/* =========================================================================
 * 4. TRC dispatch: perpendicular slice produces circle
 * ========================================================================= */

TEST(trc_perpendicular_circle) {
    alea_primitive_data_t data;
    memset(&data, 0, sizeof(data));
    data.trc.base_x = 0; data.trc.base_y = 0; data.trc.base_z = 0;
    data.trc.height_x = 0; data.trc.height_y = 0; data.trc.height_z = 10;
    data.trc.base_radius = 5.0;
    data.trc.top_radius = 2.0;

    /* Slice at z=5 → radius = 5 + 5/10*(2-5) = 5 - 1.5 = 3.5 */
    alea_slice_plane_t plane;
    make_axis_plane(&plane, 2, 5.0);

    alea_curve_2d_t curve;
    bool hit = alea_intersect_primitive_plane(ALEA_PRIMITIVE_TRC, &data, &plane, &curve);
    ASSERT(hit);
    ASSERT_EQ(curve.type, ALEA_CURVE_CIRCLE);
    ASSERT_NEAR(curve.data.circle.radius, 3.5, 1e-6);
}

TEST(trc_perpendicular_at_base) {
    alea_primitive_data_t data;
    memset(&data, 0, sizeof(data));
    data.trc.base_x = 0; data.trc.base_y = 0; data.trc.base_z = 0;
    data.trc.height_x = 0; data.trc.height_y = 0; data.trc.height_z = 10;
    data.trc.base_radius = 5.0;
    data.trc.top_radius = 2.0;

    /* Slice at z=0 → radius = 5 */
    alea_slice_plane_t plane;
    make_axis_plane(&plane, 2, 0.0);

    alea_curve_2d_t curve;
    bool hit = alea_intersect_primitive_plane(ALEA_PRIMITIVE_TRC, &data, &plane, &curve);
    ASSERT(hit);
    ASSERT_EQ(curve.type, ALEA_CURVE_CIRCLE);
    ASSERT_NEAR(curve.data.circle.radius, 5.0, 1e-6);
}

TEST(trc_perpendicular_at_top) {
    alea_primitive_data_t data;
    memset(&data, 0, sizeof(data));
    data.trc.base_x = 0; data.trc.base_y = 0; data.trc.base_z = 0;
    data.trc.height_x = 0; data.trc.height_y = 0; data.trc.height_z = 10;
    data.trc.base_radius = 5.0;
    data.trc.top_radius = 2.0;

    /* Slice at z=10 → radius = 2 */
    alea_slice_plane_t plane;
    make_axis_plane(&plane, 2, 10.0);

    alea_curve_2d_t curve;
    bool hit = alea_intersect_primitive_plane(ALEA_PRIMITIVE_TRC, &data, &plane, &curve);
    ASSERT(hit);
    ASSERT_EQ(curve.type, ALEA_CURVE_CIRCLE);
    ASSERT_NEAR(curve.data.circle.radius, 2.0, 1e-6);
}

/* =========================================================================
 * 5. Ellipse caching: oblique cylinder slice → valid canonical form
 * ========================================================================= */

TEST(oblique_cylinder_ellipse) {
    alea_primitive_data_t data;
    memset(&data, 0, sizeof(data));
    data.cyl_z.center_x = 0;
    data.cyl_z.center_y = 0;
    data.cyl_z.radius = 3.0;

    /* 45° oblique plane (normal in XZ plane) → ellipse */
    alea_slice_plane_t plane;
    alea_slice_plane_init(&plane,
                          0, 0, 0,       /* origin */
                          1, 0, 1,       /* normal (45°) */
                          0, 1, 0);      /* up */

    alea_curve_2d_t curve;
    bool hit = alea_intersect_primitive_plane(ALEA_PRIMITIVE_CYLINDER_Z, &data, &plane, &curve);
    ASSERT(hit);
    ASSERT_EQ(curve.type, ALEA_CURVE_ELLIPSE);

    /* Canonical form should be pre-computed in data.ellipse */
    double minor = fmin(curve.data.ellipse.semi_a, curve.data.ellipse.semi_b);
    double major = fmax(curve.data.ellipse.semi_a, curve.data.ellipse.semi_b);

    /* Minor axis = cylinder radius = 3 */
    ASSERT_NEAR(minor, 3.0, 0.05);
    /* Major axis = r / cos(45°) = 3√2 ≈ 4.243 */
    ASSERT_NEAR(major, 3.0 * M_SQRT2, 0.05);

    /* Verify the ellipse eval works (read from canonical form) */
    double u, v;
    ASSERT(alea_curve_eval(&curve, 0, &u, &v));
    /* At t=0: point should be at distance semi_a from center along major axis */
    double du = u - curve.data.ellipse.center[0];
    double dv = v - curve.data.ellipse.center[1];
    double dist = sqrt(du*du + dv*dv);
    ASSERT_NEAR(dist, curve.data.ellipse.semi_a, 1e-6);
}

/* =========================================================================
 * 6. Box polygon: vertices should form a convex polygon
 * ========================================================================= */

TEST(box_polygon_ordering) {
    alea_primitive_data_t data;
    memset(&data, 0, sizeof(data));
    data.box.min_x = -5; data.box.max_x = 5;
    data.box.min_y = -5; data.box.max_y = 5;
    data.box.min_z = -5; data.box.max_z = 5;

    /* XY-plane slice at z=0 → square with 4 vertices */
    alea_slice_plane_t plane;
    make_axis_plane(&plane, 2, 0.0);

    alea_curve_2d_t curve;
    bool hit = alea_intersect_primitive_plane(ALEA_PRIMITIVE_RPP, &data, &plane, &curve);
    ASSERT(hit);
    ASSERT_EQ(curve.type, ALEA_CURVE_POLYGON);
    ASSERT_EQ(curve.data.polygon.vertex_count, 4);

    /* Verify consistent winding: cross product of consecutive edges
     * should have same sign for all consecutive pairs */
    int n = curve.data.polygon.vertex_count;
    int positive = 0, negative = 0;
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        int k = (i + 2) % n;
        double e1u = curve.data.polygon.vertices[j][0] - curve.data.polygon.vertices[i][0];
        double e1v = curve.data.polygon.vertices[j][1] - curve.data.polygon.vertices[i][1];
        double e2u = curve.data.polygon.vertices[k][0] - curve.data.polygon.vertices[j][0];
        double e2v = curve.data.polygon.vertices[k][1] - curve.data.polygon.vertices[j][1];
        double cross = e1u * e2v - e1v * e2u;
        if (cross > 0) positive++;
        else if (cross < 0) negative++;
    }
    /* All same sign = convex and consistently wound */
    ASSERT(positive == 0 || negative == 0);
}

TEST(box_oblique_polygon) {
    alea_primitive_data_t data;
    memset(&data, 0, sizeof(data));
    data.box.min_x = -2; data.box.max_x = 2;
    data.box.min_y = -2; data.box.max_y = 2;
    data.box.min_z = -2; data.box.max_z = 2;

    /* 45° oblique slice → should produce 3-6 vertices */
    alea_slice_plane_t plane;
    alea_slice_plane_init(&plane, 0, 0, 0, 1, 1, 1, 0, 0, 1);

    alea_curve_2d_t curve;
    bool hit = alea_intersect_primitive_plane(ALEA_PRIMITIVE_RPP, &data, &plane, &curve);
    ASSERT(hit);
    ASSERT_EQ(curve.type, ALEA_CURVE_POLYGON);
    ASSERT(curve.data.polygon.vertex_count >= 3);
    ASSERT(curve.data.polygon.vertex_count <= 6);
}

/* =========================================================================
 * 7. Parallel lines: cylinder sliced parallel to its axis
 * ========================================================================= */

TEST(cylinder_parallel_lines) {
    alea_primitive_data_t data;
    memset(&data, 0, sizeof(data));
    data.cyl_z.center_x = 0;
    data.cyl_z.center_y = 0;
    data.cyl_z.radius = 5.0;

    /* XZ-plane at y=0 (parallel to Z axis, through center) → 2 parallel lines */
    alea_slice_plane_t plane;
    make_axis_plane(&plane, 1, 0.0);

    alea_curve_2d_t curve;
    bool hit = alea_intersect_primitive_plane(ALEA_PRIMITIVE_CYLINDER_Z, &data, &plane, &curve);
    ASSERT(hit);
    ASSERT_EQ(curve.type, ALEA_CURVE_PARALLEL_LINES);

    /* The two lines should be separated by 2*radius = 10 */
    double u1 = curve.data.parallel_lines.point1[0];
    double u2 = curve.data.parallel_lines.point2[0];
    double gap = fabs(u2 - u1);
    ASSERT_NEAR(gap, 10.0, 0.01);
}

TEST(cylinder_parallel_lines_offset) {
    alea_primitive_data_t data;
    memset(&data, 0, sizeof(data));
    data.cyl_z.center_x = 0;
    data.cyl_z.center_y = 0;
    data.cyl_z.radius = 5.0;

    /* Slice at y=3 (off-center parallel) → lines separated by 2*sqrt(25-9)=8 */
    alea_slice_plane_t plane;
    make_axis_plane(&plane, 1, 3.0);

    alea_curve_2d_t curve;
    bool hit = alea_intersect_primitive_plane(ALEA_PRIMITIVE_CYLINDER_Z, &data, &plane, &curve);
    ASSERT(hit);
    ASSERT_EQ(curve.type, ALEA_CURVE_PARALLEL_LINES);

    double u1 = curve.data.parallel_lines.point1[0];
    double u2 = curve.data.parallel_lines.point2[0];
    double gap = fabs(u2 - u1);
    ASSERT_NEAR(gap, 8.0, 0.01);  /* 2*sqrt(25-9) = 2*4 = 8 */
}

/* =========================================================================
 * 8. Sphere intersection: basic circle
 * ========================================================================= */

TEST(sphere_perpendicular_circle) {
    alea_primitive_data_t data;
    memset(&data, 0, sizeof(data));
    data.sphere.center_x = 0;
    data.sphere.center_y = 0;
    data.sphere.center_z = 0;
    data.sphere.radius = 5.0;

    alea_slice_plane_t plane;
    make_axis_plane(&plane, 2, 3.0);  /* z=3 → r = sqrt(25-9) = 4 */

    alea_curve_2d_t curve;
    bool hit = alea_intersect_primitive_plane(ALEA_PRIMITIVE_SPHERE, &data, &plane, &curve);
    ASSERT(hit);
    ASSERT_EQ(curve.type, ALEA_CURVE_CIRCLE);
    ASSERT_NEAR(curve.data.circle.radius, 4.0, 1e-6);
}

TEST(sphere_tangent_point) {
    alea_primitive_data_t data;
    memset(&data, 0, sizeof(data));
    data.sphere.center_x = 0;
    data.sphere.center_y = 0;
    data.sphere.center_z = 0;
    data.sphere.radius = 5.0;

    /* Slice at z=5 (tangent) → single point or very small circle */
    alea_slice_plane_t plane;
    make_axis_plane(&plane, 2, 5.0);

    alea_curve_2d_t curve;
    bool hit = alea_intersect_primitive_plane(ALEA_PRIMITIVE_SPHERE, &data, &plane, &curve);
    /* Tangent: might produce POINT or tiny circle, or might return false */
    if (hit) {
        ASSERT(curve.type == ALEA_CURVE_POINT || curve.type == ALEA_CURVE_CIRCLE);
        if (curve.type == ALEA_CURVE_CIRCLE) {
            ASSERT(curve.data.circle.radius < 1e-4);
        }
    }
}

TEST(sphere_no_intersection) {
    alea_primitive_data_t data;
    memset(&data, 0, sizeof(data));
    data.sphere.center_x = 0;
    data.sphere.center_y = 0;
    data.sphere.center_z = 0;
    data.sphere.radius = 5.0;

    /* Slice at z=10 → no intersection */
    alea_slice_plane_t plane;
    make_axis_plane(&plane, 2, 10.0);

    alea_curve_2d_t curve;
    bool hit = alea_intersect_primitive_plane(ALEA_PRIMITIVE_SPHERE, &data, &plane, &curve);
    ASSERT(!hit);
}

/* =========================================================================
 * 9. Scanline intersection: circle produces correct u values
 * ========================================================================= */

TEST(circle_scanline_midplane) {
    alea_curve_2d_t curve;
    memset(&curve, 0, sizeof(curve));
    curve.type = ALEA_CURVE_CIRCLE;
    curve.data.circle.center[0] = 0;
    curve.data.circle.center[1] = 0;
    curve.data.circle.radius = 5.0;

    double u_out[4];
    int n = alea_curve_scanline_intersect(&curve, 0.0, u_out, 4);
    ASSERT_EQ(n, 2);

    double umin = fmin(u_out[0], u_out[1]);
    double umax = fmax(u_out[0], u_out[1]);
    ASSERT_NEAR(umin, -5.0, 1e-10);
    ASSERT_NEAR(umax, 5.0, 1e-10);
}

TEST(circle_scanline_offset) {
    alea_curve_2d_t curve;
    memset(&curve, 0, sizeof(curve));
    curve.type = ALEA_CURVE_CIRCLE;
    curve.data.circle.center[0] = 10;
    curve.data.circle.center[1] = 20;
    curve.data.circle.radius = 3.0;

    double u_out[4];

    /* Scanline through center → u = 10 ± 3 */
    int n = alea_curve_scanline_intersect(&curve, 20.0, u_out, 4);
    ASSERT_EQ(n, 2);
    double umin = fmin(u_out[0], u_out[1]);
    double umax = fmax(u_out[0], u_out[1]);
    ASSERT_NEAR(umin, 7.0, 1e-10);
    ASSERT_NEAR(umax, 13.0, 1e-10);

    /* Scanline outside → 0 intersections */
    n = alea_curve_scanline_intersect(&curve, 25.0, u_out, 4);
    ASSERT_EQ(n, 0);
}

TEST(ellipse_scanline) {
    /* Unrotated ellipse: semi_a=4 (horizontal), semi_b=2 (vertical) */
    alea_curve_2d_t curve;
    memset(&curve, 0, sizeof(curve));
    curve.type = ALEA_CURVE_ELLIPSE;
    curve.data.ellipse.center[0] = 0;
    curve.data.ellipse.center[1] = 0;
    curve.data.ellipse.semi_a = 4.0;
    curve.data.ellipse.semi_b = 2.0;
    curve.data.ellipse.angle = 0.0;

    double u_out[4];

    /* v=0 → u = ±4 */
    int n = alea_curve_scanline_intersect(&curve, 0.0, u_out, 4);
    ASSERT_EQ(n, 2);
    double umin = fmin(u_out[0], u_out[1]);
    double umax = fmax(u_out[0], u_out[1]);
    ASSERT_NEAR(umin, -4.0, 1e-8);
    ASSERT_NEAR(umax, 4.0, 1e-8);

    /* v=1 → u²/16 + 1/4 = 1 → u² = 12 → u = ±2√3 */
    n = alea_curve_scanline_intersect(&curve, 1.0, u_out, 4);
    ASSERT_EQ(n, 2);
    umin = fmin(u_out[0], u_out[1]);
    umax = fmax(u_out[0], u_out[1]);
    ASSERT_NEAR(umin, -2.0 * sqrt(3.0), 1e-8);
    ASSERT_NEAR(umax, 2.0 * sqrt(3.0), 1e-8);

    /* v=3 (outside) → 0 intersections */
    n = alea_curve_scanline_intersect(&curve, 3.0, u_out, 4);
    ASSERT_EQ(n, 0);
}

TEST(rotated_ellipse_scanline) {
    /* 90°-rotated ellipse: horizontal extent = semi_b, vertical = semi_a */
    alea_curve_2d_t curve;
    memset(&curve, 0, sizeof(curve));
    curve.type = ALEA_CURVE_ELLIPSE;
    curve.data.ellipse.center[0] = 0;
    curve.data.ellipse.center[1] = 0;
    curve.data.ellipse.semi_a = 4.0;
    curve.data.ellipse.semi_b = 2.0;
    curve.data.ellipse.angle = M_PI / 2.0;

    double u_out[4];

    /* v=0 → u = ±2 */
    int n = alea_curve_scanline_intersect(&curve, 0.0, u_out, 4);
    ASSERT_EQ(n, 2);
    double umin = fmin(u_out[0], u_out[1]);
    double umax = fmax(u_out[0], u_out[1]);
    ASSERT_NEAR(umin, -2.0, 1e-6);
    ASSERT_NEAR(umax, 2.0, 1e-6);

    /* v=3 → should still intersect (vertical extent = 4) */
    n = alea_curve_scanline_intersect(&curve, 3.0, u_out, 4);
    ASSERT_EQ(n, 2);

    /* v=5 → outside */
    n = alea_curve_scanline_intersect(&curve, 5.0, u_out, 4);
    ASSERT_EQ(n, 0);
}

/* =========================================================================
 * 10. Curve eval: circle and ellipse evaluation
 * ========================================================================= */

TEST(circle_eval_quadrants) {
    alea_curve_2d_t curve;
    memset(&curve, 0, sizeof(curve));
    curve.type = ALEA_CURVE_CIRCLE;
    curve.data.circle.center[0] = 1.0;
    curve.data.circle.center[1] = 2.0;
    curve.data.circle.radius = 3.0;
    curve.bounds.t_min = 0;
    curve.bounds.t_max = 2 * M_PI;

    double u, v;

    /* t=0 → (1+3, 2) = (4, 2) */
    ASSERT(alea_curve_eval(&curve, 0, &u, &v));
    ASSERT_NEAR(u, 4.0, 1e-10);
    ASSERT_NEAR(v, 2.0, 1e-10);

    /* t=π/2 → (1, 2+3) = (1, 5) */
    ASSERT(alea_curve_eval(&curve, M_PI / 2, &u, &v));
    ASSERT_NEAR(u, 1.0, 1e-10);
    ASSERT_NEAR(v, 5.0, 1e-10);

    /* t=π → (1-3, 2) = (-2, 2) */
    ASSERT(alea_curve_eval(&curve, M_PI, &u, &v));
    ASSERT_NEAR(u, -2.0, 1e-10);
    ASSERT_NEAR(v, 2.0, 1e-10);
}

TEST(ellipse_eval_quadrants) {
    alea_curve_2d_t curve;
    memset(&curve, 0, sizeof(curve));
    curve.type = ALEA_CURVE_ELLIPSE;
    curve.data.ellipse.center[0] = 0;
    curve.data.ellipse.center[1] = 0;
    curve.data.ellipse.semi_a = 4.0;
    curve.data.ellipse.semi_b = 2.0;
    curve.data.ellipse.angle = 0.0;
    curve.bounds.t_min = 0;
    curve.bounds.t_max = 2 * M_PI;

    double u, v;

    /* t=0 → (4, 0) */
    ASSERT(alea_curve_eval(&curve, 0, &u, &v));
    ASSERT_NEAR(u, 4.0, 1e-10);
    ASSERT_NEAR(v, 0.0, 1e-10);

    /* t=π/2 → (0, 2) */
    ASSERT(alea_curve_eval(&curve, M_PI / 2, &u, &v));
    ASSERT_NEAR(u, 0.0, 1e-10);
    ASSERT_NEAR(v, 2.0, 1e-10);
}

/* =========================================================================
 * 11. Bounding box: circle, ellipse, polygon
 * ========================================================================= */

TEST(circle_bbox) {
    alea_curve_2d_t curve;
    memset(&curve, 0, sizeof(curve));
    curve.type = ALEA_CURVE_CIRCLE;
    curve.data.circle.center[0] = 5.0;
    curve.data.circle.center[1] = -3.0;
    curve.data.circle.radius = 2.0;

    double u_min, u_max, v_min, v_max;
    alea_curve_bbox(&curve, &u_min, &u_max, &v_min, &v_max);

    ASSERT_NEAR(u_min, 3.0, 1e-10);
    ASSERT_NEAR(u_max, 7.0, 1e-10);
    ASSERT_NEAR(v_min, -5.0, 1e-10);
    ASSERT_NEAR(v_max, -1.0, 1e-10);
}

TEST(ellipse_bbox_unrotated) {
    alea_curve_2d_t curve;
    memset(&curve, 0, sizeof(curve));
    curve.type = ALEA_CURVE_ELLIPSE;
    curve.data.ellipse.center[0] = 0;
    curve.data.ellipse.center[1] = 0;
    curve.data.ellipse.semi_a = 4.0;
    curve.data.ellipse.semi_b = 2.0;
    curve.data.ellipse.angle = 0.0;

    double u_min, u_max, v_min, v_max;
    alea_curve_bbox(&curve, &u_min, &u_max, &v_min, &v_max);

    ASSERT_NEAR(u_min, -4.0, 1e-10);
    ASSERT_NEAR(u_max, 4.0, 1e-10);
    ASSERT_NEAR(v_min, -2.0, 1e-10);
    ASSERT_NEAR(v_max, 2.0, 1e-10);
}

TEST(ellipse_bbox_rotated_90) {
    /* 90° rotation swaps axis extents */
    alea_curve_2d_t curve;
    memset(&curve, 0, sizeof(curve));
    curve.type = ALEA_CURVE_ELLIPSE;
    curve.data.ellipse.center[0] = 0;
    curve.data.ellipse.center[1] = 0;
    curve.data.ellipse.semi_a = 4.0;
    curve.data.ellipse.semi_b = 2.0;
    curve.data.ellipse.angle = M_PI / 2.0;

    double u_min, u_max, v_min, v_max;
    alea_curve_bbox(&curve, &u_min, &u_max, &v_min, &v_max);

    ASSERT_NEAR(u_min, -2.0, 1e-6);
    ASSERT_NEAR(u_max, 2.0, 1e-6);
    ASSERT_NEAR(v_min, -4.0, 1e-6);
    ASSERT_NEAR(v_max, 4.0, 1e-6);
}

/* =========================================================================
 * 12. Plane coordinate transforms: roundtrip 3D → 2D → 3D
 * ========================================================================= */

TEST(plane_coord_roundtrip) {
    alea_slice_plane_t plane;
    alea_slice_plane_init(&plane, 1, 2, 3, 0, 0, 1, 0, 1, 0);

    double x = 4.5, y = -1.3, z = 3.0;  /* z=3 matches plane origin z */
    double u, v;
    alea_plane_to_2d(&plane, x, y, z, &u, &v);

    double rx, ry, rz;
    alea_plane_to_3d(&plane, u, v, &rx, &ry, &rz);

    ASSERT_NEAR(rz, 3.0, 1e-10);
    ASSERT_NEAR(rx, x, 1e-10);
    ASSERT_NEAR(ry, y, 1e-10);
}

TEST(plane_coord_axis_aligned) {
    alea_slice_plane_t plane;
    alea_slice_plane_init_axis(&plane, 2, 5.0);

    double u, v;
    alea_plane_to_2d(&plane, 10, 20, 5, &u, &v);

    double rx, ry, rz;
    alea_plane_to_3d(&plane, u, v, &rx, &ry, &rz);

    ASSERT_NEAR(rz, 5.0, 1e-10);
    ASSERT_NEAR(rx, 10.0, 1e-10);
    ASSERT_NEAR(ry, 20.0, 1e-10);
}

/* =========================================================================
 * 13. RCC intersection: perpendicular → circle
 * ========================================================================= */

TEST(rcc_perpendicular_circle) {
    alea_primitive_data_t data;
    memset(&data, 0, sizeof(data));
    data.rcc.base_x = 0; data.rcc.base_y = 0; data.rcc.base_z = 0;
    data.rcc.height_x = 0; data.rcc.height_y = 0; data.rcc.height_z = 10;
    data.rcc.radius = 3.0;

    /* Slice at z=5 → circle of radius 3 */
    alea_slice_plane_t plane;
    make_axis_plane(&plane, 2, 5.0);

    alea_curve_2d_t curve;
    bool hit = alea_intersect_primitive_plane(ALEA_PRIMITIVE_RCC, &data, &plane, &curve);
    ASSERT(hit);
    ASSERT_EQ(curve.type, ALEA_CURVE_CIRCLE);
    ASSERT_NEAR(curve.data.circle.radius, 3.0, 1e-6);
}

TEST_MAIN()
