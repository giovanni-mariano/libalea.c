// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_primitives.c - Primitive evaluation, bbox, and degenerate case tests
 */

#include "alea_test.h"
#include "alea_types.h"
#include "primitives/primitive_eval.h"
#include "primitives/bbox.h"
#include <string.h>
#include <math.h>

/* ========================================================================= */
/* Point evaluation tests                                                    */
/* ========================================================================= */

TEST(eval_plane_pos_neg) {
    alea_primitive_data_t d;
    memset(&d, 0, sizeof(d));
    /* Plane: z = 5 → z - 5 = 0 → a=0,b=0,c=1,d=-5 */
    d.plane.a = 0; d.plane.b = 0; d.plane.c = 1; d.plane.d = -5;

    /* Below plane: inside (-) */
    double v = alea_primitive_eval(ALEA_PRIMITIVE_PLANE, &d, 0, 0, 3);
    ASSERT(v < 0);  /* 3 - 5 = -2 */

    /* Above plane: outside (+) */
    v = alea_primitive_eval(ALEA_PRIMITIVE_PLANE, &d, 0, 0, 7);
    ASSERT(v > 0);  /* 7 - 5 = 2 */

    /* On plane */
    v = alea_primitive_eval(ALEA_PRIMITIVE_PLANE, &d, 0, 0, 5);
    ASSERT_NEAR(v, 0, 1e-10);
}

TEST(eval_sphere_inside_out) {
    alea_primitive_data_t d;
    memset(&d, 0, sizeof(d));
    d.sphere.center_x = 0; d.sphere.center_y = 0; d.sphere.center_z = 0;
    d.sphere.radius = 5.0;

    /* Inside */
    double v = alea_primitive_eval(ALEA_PRIMITIVE_SPHERE, &d, 1, 1, 1);
    ASSERT(v < 0);

    /* Outside */
    v = alea_primitive_eval(ALEA_PRIMITIVE_SPHERE, &d, 10, 0, 0);
    ASSERT(v > 0);

    /* On surface */
    v = alea_primitive_eval(ALEA_PRIMITIVE_SPHERE, &d, 5, 0, 0);
    ASSERT_NEAR(v, 0, 1e-6);
}

TEST(eval_cylinder_x) {
    alea_primitive_data_t d;
    memset(&d, 0, sizeof(d));
    d.cyl_x.center_y = 0; d.cyl_x.center_z = 0;
    d.cyl_x.radius = 3.0;

    /* Inside: (any x, y=1, z=1) → sqrt(1+1) = 1.41 < 3 */
    double v = alea_primitive_eval(ALEA_PRIMITIVE_CYLINDER_X, &d, 100, 1, 1);
    ASSERT(v < 0);

    /* Outside: (any x, y=5, z=0) → 5 > 3 */
    v = alea_primitive_eval(ALEA_PRIMITIVE_CYLINDER_X, &d, 0, 5, 0);
    ASSERT(v > 0);
}

TEST(eval_cylinder_y) {
    alea_primitive_data_t d;
    memset(&d, 0, sizeof(d));
    d.cyl_y.center_x = 0; d.cyl_y.center_z = 0;
    d.cyl_y.radius = 3.0;

    double v = alea_primitive_eval(ALEA_PRIMITIVE_CYLINDER_Y, &d, 1, 100, 1);
    ASSERT(v < 0);

    v = alea_primitive_eval(ALEA_PRIMITIVE_CYLINDER_Y, &d, 5, 0, 0);
    ASSERT(v > 0);
}

TEST(eval_cylinder_z) {
    alea_primitive_data_t d;
    memset(&d, 0, sizeof(d));
    d.cyl_z.center_x = 0; d.cyl_z.center_y = 0;
    d.cyl_z.radius = 3.0;

    double v = alea_primitive_eval(ALEA_PRIMITIVE_CYLINDER_Z, &d, 1, 1, 100);
    ASSERT(v < 0);

    v = alea_primitive_eval(ALEA_PRIMITIVE_CYLINDER_Z, &d, 5, 0, 0);
    ASSERT(v > 0);
}

TEST(eval_cone_x) {
    alea_primitive_data_t d;
    memset(&d, 0, sizeof(d));
    d.cone_x.apex_x = 0; d.cone_x.apex_y = 0; d.cone_x.apex_z = 0;
    d.cone_x.tan_angle_sq = 1.0;  /* 45° half-angle */
    d.cone_x.sheet_selection = 0;  /* both sheets */

    /* At x=5, cone radius = 5 (tan²=1 → r=|x|) */
    /* Point (5, 2, 0): dist_from_axis = 2 < 5 → inside */
    double v = alea_primitive_eval(ALEA_PRIMITIVE_CONE_X, &d, 5, 2, 0);
    ASSERT(v < 0);

    /* Point (5, 6, 0): 6 > 5 → outside */
    v = alea_primitive_eval(ALEA_PRIMITIVE_CONE_X, &d, 5, 6, 0);
    ASSERT(v > 0);
}

TEST(eval_cone_y) {
    alea_primitive_data_t d;
    memset(&d, 0, sizeof(d));
    d.cone_y.apex_x = 0; d.cone_y.apex_y = 0; d.cone_y.apex_z = 0;
    d.cone_y.tan_angle_sq = 1.0;
    d.cone_y.sheet_selection = 0;

    double v = alea_primitive_eval(ALEA_PRIMITIVE_CONE_Y, &d, 2, 5, 0);
    ASSERT(v < 0);
}

TEST(eval_cone_z) {
    alea_primitive_data_t d;
    memset(&d, 0, sizeof(d));
    d.cone_z.apex_x = 0; d.cone_z.apex_y = 0; d.cone_z.apex_z = 0;
    d.cone_z.tan_angle_sq = 1.0;
    d.cone_z.sheet_selection = 0;

    double v = alea_primitive_eval(ALEA_PRIMITIVE_CONE_Z, &d, 2, 0, 5);
    ASSERT(v < 0);
}

TEST(eval_quadric) {
    /* GQ: x² + y² + z² - 25 = 0 (sphere of radius 5) */
    alea_primitive_data_t d;
    memset(&d, 0, sizeof(d));
    d.quadric.coeffs[0] = 1.0;  /* A: x² */
    d.quadric.coeffs[1] = 1.0;  /* B: y² */
    d.quadric.coeffs[2] = 1.0;  /* C: z² */
    d.quadric.coeffs[9] = -25.0; /* K */

    double v = alea_primitive_eval(ALEA_PRIMITIVE_QUADRIC, &d, 0, 0, 0);
    ASSERT(v < 0);

    v = alea_primitive_eval(ALEA_PRIMITIVE_QUADRIC, &d, 10, 0, 0);
    ASSERT(v > 0);
}

TEST(eval_torus_z) {
    alea_primitive_data_t d;
    memset(&d, 0, sizeof(d));
    d.torus.axis = ALEA_AXIS_Z;
    d.torus.center_x = 0; d.torus.center_y = 0; d.torus.center_z = 0;
    d.torus.major_radius = 5.0;
    d.torus.minor_radius = 1.0;
    d.torus.axial_semiwidth_B = 1.0;

    /* On major circle in tube: (5, 0, 0) → inside tube */
    double v = alea_primitive_eval(ALEA_PRIMITIVE_TORUS_Z, &d, 5, 0, 0);
    ASSERT(v < 0);

    /* At origin → outside tube (distance from axis = 0, major = 5) */
    v = alea_primitive_eval(ALEA_PRIMITIVE_TORUS_Z, &d, 0, 0, 0);
    ASSERT(v > 0);
}

TEST(eval_torus_x) {
    alea_primitive_data_t d;
    memset(&d, 0, sizeof(d));
    d.torus.axis = ALEA_AXIS_X;
    d.torus.center_x = 0; d.torus.center_y = 0; d.torus.center_z = 0;
    d.torus.major_radius = 5.0;
    d.torus.minor_radius = 1.0;
    d.torus.axial_semiwidth_B = 1.0;

    /* On major circle */
    double v = alea_primitive_eval(ALEA_PRIMITIVE_TORUS_X, &d, 0, 5, 0);
    ASSERT(v < 0);

    v = alea_primitive_eval(ALEA_PRIMITIVE_TORUS_X, &d, 0, 0, 0);
    ASSERT(v > 0);
}

TEST(eval_torus_y) {
    alea_primitive_data_t d;
    memset(&d, 0, sizeof(d));
    d.torus.axis = ALEA_AXIS_Y;
    d.torus.center_x = 0; d.torus.center_y = 0; d.torus.center_z = 0;
    d.torus.major_radius = 5.0;
    d.torus.minor_radius = 1.0;
    d.torus.axial_semiwidth_B = 1.0;

    double v = alea_primitive_eval(ALEA_PRIMITIVE_TORUS_Y, &d, 5, 0, 0);
    ASSERT(v < 0);

    v = alea_primitive_eval(ALEA_PRIMITIVE_TORUS_Y, &d, 0, 0, 0);
    ASSERT(v > 0);
}

TEST(eval_rcc) {
    alea_primitive_data_t d;
    memset(&d, 0, sizeof(d));
    d.rcc.base_x = 0; d.rcc.base_y = 0; d.rcc.base_z = 0;
    d.rcc.height_x = 0; d.rcc.height_y = 0; d.rcc.height_z = 10;
    d.rcc.radius = 3.0;

    /* Inside: (0, 0, 5) */
    double v = alea_primitive_eval(ALEA_PRIMITIVE_RCC, &d, 0, 0, 5);
    ASSERT(v < 0);

    /* Outside radially: (5, 0, 5) */
    v = alea_primitive_eval(ALEA_PRIMITIVE_RCC, &d, 5, 0, 5);
    ASSERT(v > 0);

    /* Outside axially: (0, 0, 15) */
    v = alea_primitive_eval(ALEA_PRIMITIVE_RCC, &d, 0, 0, 15);
    ASSERT(v > 0);
}

TEST(eval_box_general) {
    alea_primitive_data_t d;
    memset(&d, 0, sizeof(d));
    d.box_general.corner_x = 0; d.box_general.corner_y = 0; d.box_general.corner_z = 0;
    d.box_general.v1_x = 10; d.box_general.v1_y = 0; d.box_general.v1_z = 0;
    d.box_general.v2_x = 0; d.box_general.v2_y = 5; d.box_general.v2_z = 0;
    d.box_general.v3_x = 0; d.box_general.v3_y = 0; d.box_general.v3_z = 3;

    /* Inside: (5, 2.5, 1.5) */
    double v = alea_primitive_eval(ALEA_PRIMITIVE_BOX, &d, 5, 2.5, 1.5);
    ASSERT(v < 0);

    /* Outside: (12, 0, 0) */
    v = alea_primitive_eval(ALEA_PRIMITIVE_BOX, &d, 12, 0, 0);
    ASSERT(v > 0);
}

TEST(eval_sph) {
    alea_primitive_data_t d;
    memset(&d, 0, sizeof(d));
    d.sph.center_x = 3; d.sph.center_y = 4; d.sph.center_z = 5;
    d.sph.radius = 2.0;

    double v = alea_primitive_eval(ALEA_PRIMITIVE_SPH, &d, 3, 4, 5);
    ASSERT(v < 0);

    v = alea_primitive_eval(ALEA_PRIMITIVE_SPH, &d, 3, 4, 8);
    ASSERT(v > 0);
}

TEST(eval_trc) {
    alea_primitive_data_t d;
    memset(&d, 0, sizeof(d));
    d.trc.base_x = 0; d.trc.base_y = 0; d.trc.base_z = 0;
    d.trc.height_x = 0; d.trc.height_y = 0; d.trc.height_z = 10;
    d.trc.base_radius = 3.0;
    d.trc.top_radius = 1.0;

    /* Inside at mid-height: (0, 0, 5) */
    double v = alea_primitive_eval(ALEA_PRIMITIVE_TRC, &d, 0, 0, 5);
    ASSERT(v < 0);

    /* Outside radially */
    v = alea_primitive_eval(ALEA_PRIMITIVE_TRC, &d, 5, 0, 5);
    ASSERT(v > 0);
}

TEST(eval_wed) {
    alea_primitive_data_t d;
    memset(&d, 0, sizeof(d));
    d.wed.vertex_x = 0; d.wed.vertex_y = 0; d.wed.vertex_z = 0;
    d.wed.v1_x = 4; d.wed.v1_y = 0; d.wed.v1_z = 0;
    d.wed.v2_x = 0; d.wed.v2_y = 4; d.wed.v2_z = 0;
    d.wed.v3_x = 0; d.wed.v3_y = 0; d.wed.v3_z = 4;

    /* Inside the triangular prism: (0.5, 0.5, 2) */
    double v = alea_primitive_eval(ALEA_PRIMITIVE_WED, &d, 0.5, 0.5, 2);
    ASSERT(v < 0);

    /* Outside */
    v = alea_primitive_eval(ALEA_PRIMITIVE_WED, &d, 5, 5, 5);
    ASSERT(v > 0);
}

TEST(eval_rhp) {
    alea_primitive_data_t d;
    memset(&d, 0, sizeof(d));
    d.rhp.base_x = 0; d.rhp.base_y = 0; d.rhp.base_z = 0;
    d.rhp.height_x = 0; d.rhp.height_y = 0; d.rhp.height_z = 10;
    d.rhp.r1_x = 2; d.rhp.r1_y = 0; d.rhp.r1_z = 0;
    d.rhp.r2_x = -1; d.rhp.r2_y = 1.732050808; d.rhp.r2_z = 0;
    d.rhp.r3_x = -1; d.rhp.r3_y = -1.732050808; d.rhp.r3_z = 0;

    /* Inside: center (0, 0, 5) */
    double v = alea_primitive_eval(ALEA_PRIMITIVE_RHP, &d, 0, 0, 5);
    ASSERT(v < 0);

    /* Outside: far away */
    v = alea_primitive_eval(ALEA_PRIMITIVE_RHP, &d, 5, 0, 5);
    ASSERT(v > 0);
}

/* ========================================================================= */
/* Bounding box tests                                                        */
/* ========================================================================= */

TEST(bbox_sphere) {
    alea_primitive_data_t d;
    memset(&d, 0, sizeof(d));
    d.sphere.center_x = 3; d.sphere.center_y = 4; d.sphere.center_z = 5;
    d.sphere.radius = 2.0;

    alea_bbox_t bb = alea_primitive_bbox(ALEA_PRIMITIVE_SPHERE, &d);
    ASSERT_NEAR(bb.min_x, 1.0, 1e-6);
    ASSERT_NEAR(bb.max_x, 5.0, 1e-6);
    ASSERT_NEAR(bb.min_y, 2.0, 1e-6);
    ASSERT_NEAR(bb.max_y, 6.0, 1e-6);
    ASSERT_NEAR(bb.min_z, 3.0, 1e-6);
    ASSERT_NEAR(bb.max_z, 7.0, 1e-6);
}

TEST(bbox_cylinder_z) {
    alea_primitive_data_t d;
    memset(&d, 0, sizeof(d));
    d.cyl_z.center_x = 0; d.cyl_z.center_y = 0;
    d.cyl_z.radius = 3.0;

    alea_bbox_t bb = alea_primitive_bbox(ALEA_PRIMITIVE_CYLINDER_Z, &d);
    ASSERT_NEAR(bb.min_x, -3.0, 1e-6);
    ASSERT_NEAR(bb.max_x, 3.0, 1e-6);
    ASSERT_NEAR(bb.min_y, -3.0, 1e-6);
    ASSERT_NEAR(bb.max_y, 3.0, 1e-6);
    /* Z should be unbounded (very large) */
    ASSERT(bb.max_z > 1e10);
}

TEST(bbox_box) {
    alea_primitive_data_t d;
    memset(&d, 0, sizeof(d));
    d.box.min_x = -5; d.box.max_x = 5;
    d.box.min_y = -3; d.box.max_y = 3;
    d.box.min_z = -2; d.box.max_z = 2;

    alea_bbox_t bb = alea_primitive_bbox(ALEA_PRIMITIVE_RPP, &d);
    ASSERT_NEAR(bb.min_x, -5.0, 1e-6);
    ASSERT_NEAR(bb.max_x, 5.0, 1e-6);
    ASSERT_NEAR(bb.min_y, -3.0, 1e-6);
    ASSERT_NEAR(bb.max_y, 3.0, 1e-6);
    ASSERT_NEAR(bb.min_z, -2.0, 1e-6);
    ASSERT_NEAR(bb.max_z, 2.0, 1e-6);
}

TEST(bbox_rcc) {
    alea_primitive_data_t d;
    memset(&d, 0, sizeof(d));
    d.rcc.base_x = 0; d.rcc.base_y = 0; d.rcc.base_z = 0;
    d.rcc.height_x = 0; d.rcc.height_y = 0; d.rcc.height_z = 10;
    d.rcc.radius = 3.0;

    alea_bbox_t bb = alea_primitive_bbox(ALEA_PRIMITIVE_RCC, &d);
    ASSERT_NEAR(bb.min_x, -3.0, 1e-6);
    ASSERT_NEAR(bb.max_x, 3.0, 1e-6);
    ASSERT_NEAR(bb.min_y, -3.0, 1e-6);
    ASSERT_NEAR(bb.max_y, 3.0, 1e-6);
    ASSERT_NEAR(bb.min_z, 0.0, 1e-6);
    ASSERT_NEAR(bb.max_z, 10.0, 1e-6);
}

TEST(bbox_torus) {
    alea_primitive_data_t d;
    memset(&d, 0, sizeof(d));
    d.torus.axis = ALEA_AXIS_Z;
    d.torus.center_x = 0; d.torus.center_y = 0; d.torus.center_z = 0;
    d.torus.major_radius = 5.0;
    d.torus.minor_radius = 1.0;
    d.torus.axial_semiwidth_B = 1.0;

    alea_bbox_t bb = alea_primitive_bbox(ALEA_PRIMITIVE_TORUS_Z, &d);
    ASSERT_NEAR(bb.min_x, -6.0, 1e-6);
    ASSERT_NEAR(bb.max_x, 6.0, 1e-6);
    ASSERT_NEAR(bb.min_y, -6.0, 1e-6);
    ASSERT_NEAR(bb.max_y, 6.0, 1e-6);
    ASSERT_NEAR(bb.min_z, -1.0, 1e-6);
    ASSERT_NEAR(bb.max_z, 1.0, 1e-6);
}

/* ========================================================================= */
/* Degenerate cases                                                          */
/* ========================================================================= */

TEST(eval_zero_radius_sphere) {
    alea_primitive_data_t d;
    memset(&d, 0, sizeof(d));
    d.sphere.center_x = 0; d.sphere.center_y = 0; d.sphere.center_z = 0;
    d.sphere.radius = 0.0;

    /* At center: r²-R² = 0-0 = 0 → on surface */
    double v = alea_primitive_eval(ALEA_PRIMITIVE_SPHERE, &d, 0, 0, 0);
    ASSERT_NEAR(v, 0, 1e-10);

    /* Any other point → outside */
    v = alea_primitive_eval(ALEA_PRIMITIVE_SPHERE, &d, 1, 0, 0);
    ASSERT(v > 0);
}

TEST(eval_zero_radius_cyl) {
    alea_primitive_data_t d;
    memset(&d, 0, sizeof(d));
    d.cyl_z.center_x = 0; d.cyl_z.center_y = 0;
    d.cyl_z.radius = 0.0;

    /* On axis → on surface */
    double v = alea_primitive_eval(ALEA_PRIMITIVE_CYLINDER_Z, &d, 0, 0, 0);
    ASSERT_NEAR(v, 0, 1e-10);

    /* Off axis → outside */
    v = alea_primitive_eval(ALEA_PRIMITIVE_CYLINDER_Z, &d, 1, 0, 0);
    ASSERT(v > 0);
}

/* ========================================================================= */
/* Sense tests                                                               */
/* ========================================================================= */

TEST(eval_with_sense) {
    alea_primitive_data_t d;
    memset(&d, 0, sizeof(d));
    d.sphere.center_x = 0; d.sphere.center_y = 0; d.sphere.center_z = 0;
    d.sphere.radius = 5.0;

    /* Inside sphere, sense -1 (negative/inside) → negative */
    double v = alea_primitive_eval_with_sense(ALEA_PRIMITIVE_SPHERE, &d, -1, 0, 0, 0);
    ASSERT(v < 0);

    /* Inside sphere, sense +1 (positive/outside) → positive */
    v = alea_primitive_eval_with_sense(ALEA_PRIMITIVE_SPHERE, &d, +1, 0, 0, 0);
    ASSERT(v > 0);
}

TEST_MAIN()
