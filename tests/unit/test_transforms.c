// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_transforms.c - Transform creation, composition, and application tests
 */

#include "alea_test.h"
#include "alea.h"
#include "core/alea_system.h"
#include "core/alea_transform.h"
#include <string.h>
#include <math.h>

/* ========================================================================= */
/* Transform basic tests                                                     */
/* ========================================================================= */

TEST(transform_identity) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Identity transform: (0,0,0) translation, identity rotation */
    double data[12] = {0,0,0, 1,0,0, 0,1,0, 0,0,1};
    int rc = alea_add_transform(sys, 1, data, 12, 0);
    ASSERT_EQ(rc, 0);

    const alea_transform_t* tr = alea_get_transform(sys, 1);
    ASSERT_NOT_NULL(tr);

    /* Apply to point (3, 4, 5) → should remain unchanged */
    double x = 3, y = 4, z = 5;
    alea_transform_point(tr, &x, &y, &z);
    ASSERT_NEAR(x, 3.0, 1e-10);
    ASSERT_NEAR(y, 4.0, 1e-10);
    ASSERT_NEAR(z, 5.0, 1e-10);

    alea_destroy(sys);
}

TEST(transform_translation) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Pure translation (10, 20, 30) */
    double data[3] = {10, 20, 30};
    int rc = alea_add_transform(sys, 1, data, 3, 0);
    ASSERT_EQ(rc, 0);

    const alea_transform_t* tr = alea_get_transform(sys, 1);
    ASSERT_NOT_NULL(tr);
    ASSERT_TRUE(alea_transform_is_translation_only(tr));

    double x = 0, y = 0, z = 0;
    alea_transform_point(tr, &x, &y, &z);
    ASSERT_NEAR(x, 10.0, 1e-10);
    ASSERT_NEAR(y, 20.0, 1e-10);
    ASSERT_NEAR(z, 30.0, 1e-10);

    alea_destroy(sys);
}

TEST(transform_rotation_90z) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* 90° rotation about Z: x→y, y→-x, z→z
     * Cosine matrix:
     *   cos(90°) = 0, cos(0°) = 1, cos(90°) = 0
     *   cos(180°)=-1, cos(90°)= 0, cos(90°) = 0 (Wait, let me think...)
     *
     * MCNP *TRn with angles in degrees:
     * Translation: 0 0 0
     * Angles: θ11 θ21 θ31 θ12 θ22 θ32 θ13 θ23 θ33
     * where θij is angle between i-th aux axis and j-th main axis
     * For 90° rotation about Z:
     *   x_aux → y_main: θ11=90, θ21=0, θ31=90
     *   y_aux → -x_main: θ12=180, θ22=90, θ32=90
     *   z_aux → z_main: θ13=90, θ23=90, θ33=0
     */
    double data[12] = {0,0,0, 90,0,90, 180,90,90, 90,90,0};
    int rc = alea_add_transform(sys, 1, data, 12, 1); /* degrees=1 */
    ASSERT_EQ(rc, 0);

    const alea_transform_t* tr = alea_get_transform(sys, 1);
    ASSERT_NOT_NULL(tr);
    ASSERT_FALSE(alea_transform_is_translation_only(tr));

    /* Point (1, 0, 0) should map to (0, 1, 0) */
    double x = 1, y = 0, z = 0;
    alea_transform_point(tr, &x, &y, &z);
    ASSERT_NEAR(x, 0.0, 1e-10);
    ASSERT_NEAR(y, 1.0, 1e-10);
    ASSERT_NEAR(z, 0.0, 1e-10);

    alea_destroy(sys);
}

TEST(transform_full_12) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Translation (10, 0, 0) with identity rotation */
    double data[12] = {10,0,0, 1,0,0, 0,1,0, 0,0,1};
    int rc = alea_add_transform(sys, 1, data, 12, 0);
    ASSERT_EQ(rc, 0);

    const alea_transform_t* tr = alea_get_transform(sys, 1);
    ASSERT_NOT_NULL(tr);

    double x = 0, y = 0, z = 0;
    alea_transform_point(tr, &x, &y, &z);
    ASSERT_NEAR(x, 10.0, 1e-10);
    ASSERT_NEAR(y, 0.0, 1e-10);
    ASSERT_NEAR(z, 0.0, 1e-10);

    alea_destroy(sys);
}

TEST(transform_degrees_vs_cos) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Identity rotation as direction cosines */
    double data_cos[12] = {5,0,0, 1,0,0, 0,1,0, 0,0,1};
    int rc = alea_add_transform(sys, 1, data_cos, 12, 0);
    ASSERT_EQ(rc, 0);

    /* Identity rotation as angles in degrees */
    double data_deg[12] = {5,0,0, 0,90,90, 90,0,90, 90,90,0};
    rc = alea_add_transform(sys, 2, data_deg, 12, 1);
    ASSERT_EQ(rc, 0);

    const alea_transform_t* tr1 = alea_get_transform(sys, 1);
    const alea_transform_t* tr2 = alea_get_transform(sys, 2);
    ASSERT_NOT_NULL(tr1);
    ASSERT_NOT_NULL(tr2);

    /* Both should produce the same result */
    double x1 = 3, y1 = 4, z1 = 7;
    double x2 = 3, y2 = 4, z2 = 7;
    alea_transform_point(tr1, &x1, &y1, &z1);
    alea_transform_point(tr2, &x2, &y2, &z2);
    ASSERT_NEAR(x1, x2, 1e-10);
    ASSERT_NEAR(y1, y2, 1e-10);
    ASSERT_NEAR(z1, z2, 1e-10);

    alea_destroy(sys);
}

TEST(transform_apply_plane) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Translation (10, 0, 0) */
    double tr_data[3] = {10, 0, 0};
    alea_add_transform(sys, 1, tr_data, 3, 0);
    const alea_transform_t* tr = alea_get_transform(sys, 1);
    ASSERT_NOT_NULL(tr);

    /* Plane PZ 5.0: a=0, b=0, c=1, d=-5 */
    alea_primitive_data_t in_data, out_data;
    memset(&in_data, 0, sizeof(in_data));
    in_data.plane.a = 0; in_data.plane.b = 0; in_data.plane.c = 1; in_data.plane.d = -5;

    alea_primitive_type_t out_type;
    bool ok = alea_apply_transform_to_primitive(tr, ALEA_PRIMITIVE_PLANE, &in_data,
                                                &out_type, &out_data);
    ASSERT_TRUE(ok);
    ASSERT_EQ(out_type, ALEA_PRIMITIVE_PLANE);
    /* PZ is along z-axis. Translation (10,0,0) doesn't change z-intercept.
       So plane should remain at z=5: c=1, d=-5 */
    ASSERT_NEAR(out_data.plane.c, 1.0, 1e-10);
    ASSERT_NEAR(out_data.plane.d, -5.0, 1e-10);

    alea_destroy(sys);
}

TEST(transform_apply_sphere) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    double tr_data[3] = {10, 20, 30};
    alea_add_transform(sys, 1, tr_data, 3, 0);
    const alea_transform_t* tr = alea_get_transform(sys, 1);
    ASSERT_NOT_NULL(tr);

    alea_primitive_data_t in_data, out_data;
    memset(&in_data, 0, sizeof(in_data));
    in_data.sphere.center_x = 0; in_data.sphere.center_y = 0;
    in_data.sphere.center_z = 0; in_data.sphere.radius = 5.0;

    alea_primitive_type_t out_type;
    bool ok = alea_apply_transform_to_primitive(tr, ALEA_PRIMITIVE_SPHERE, &in_data,
                                                &out_type, &out_data);
    ASSERT_TRUE(ok);
    ASSERT_EQ(out_type, ALEA_PRIMITIVE_SPHERE);
    ASSERT_NEAR(out_data.sphere.center_x, 10.0, 1e-10);
    ASSERT_NEAR(out_data.sphere.center_y, 20.0, 1e-10);
    ASSERT_NEAR(out_data.sphere.center_z, 30.0, 1e-10);
    ASSERT_NEAR(out_data.sphere.radius, 5.0, 1e-10);

    alea_destroy(sys);
}

TEST(transform_apply_cylinder) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    double tr_data[3] = {5, 10, 0};
    alea_add_transform(sys, 1, tr_data, 3, 0);
    const alea_transform_t* tr = alea_get_transform(sys, 1);
    ASSERT_NOT_NULL(tr);

    alea_primitive_data_t in_data, out_data;
    memset(&in_data, 0, sizeof(in_data));
    in_data.cyl_z.center_x = 0; in_data.cyl_z.center_y = 0;
    in_data.cyl_z.radius = 3.0;

    alea_primitive_type_t out_type;
    bool ok = alea_apply_transform_to_primitive(tr, ALEA_PRIMITIVE_CYLINDER_Z, &in_data,
                                                &out_type, &out_data);
    ASSERT_TRUE(ok);
    /* Translation-only should preserve cylinder type */
    ASSERT_EQ(out_type, ALEA_PRIMITIVE_CYLINDER_Z);
    ASSERT_NEAR(out_data.cyl_z.center_x, 5.0, 1e-10);
    ASSERT_NEAR(out_data.cyl_z.center_y, 10.0, 1e-10);
    ASSERT_NEAR(out_data.cyl_z.radius, 3.0, 1e-10);

    alea_destroy(sys);
}

TEST(transform_vector) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Translation only: vector should be unchanged */
    double tr_data[3] = {10, 20, 30};
    alea_add_transform(sys, 1, tr_data, 3, 0);
    const alea_transform_t* tr = alea_get_transform(sys, 1);
    ASSERT_NOT_NULL(tr);

    double vx = 1, vy = 0, vz = 0;
    alea_transform_vector(tr, &vx, &vy, &vz);
    ASSERT_NEAR(vx, 1.0, 1e-10);
    ASSERT_NEAR(vy, 0.0, 1e-10);
    ASSERT_NEAR(vz, 0.0, 1e-10);

    alea_destroy(sys);
}

TEST(transform_inline) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    double tr_data[3] = {5, 10, 15};
    int id = alea_add_inline_transform(sys, tr_data, 3, 0);
    ASSERT(id > 0);

    const alea_transform_t* tr = alea_get_transform(sys, id);
    ASSERT_NOT_NULL(tr);
    ASSERT_TRUE(alea_transform_is_translation_only(tr));

    alea_destroy(sys);
}

TEST(transform_surface_with_tr) {
    /* Test surface with TRn reference via MCNP parse */
    const char* input =
        "Test surf TR\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 1 SO 5.0\n"
        "\n"
        "TR1 10 0 0\n"
        "M1 92235.80c 1.0\n";
    alea_system_t* sys = alea_load_mcnp_string(input, strlen(input));
    ASSERT_NOT_NULL(sys);
    alea_build_universe_index(sys);

    /* Sphere should be centered at (10,0,0) */
    ASSERT_EQ(alea_material_at(sys, 10, 0, 0), 1);
    ASSERT_EQ(alea_material_at(sys, 0, 0, 0), 0);

    alea_destroy(sys);
}

TEST_MAIN()
