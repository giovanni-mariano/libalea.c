// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_tolerance.c - Numerical tolerance and dedup tests
 */

#include "alea_test.h"
#include "alea.h"
#include "core/alea_system.h"
#include <string.h>
#include <math.h>

/* ========================================================================= */
/* Near-equality tests                                                       */
/* ========================================================================= */

TEST(tol_planes_equal) {
    alea_config_t cfg = ALEA_CONFIG_DEFAULT;
    alea_primitive_data_t d1, d2;
    memset(&d1, 0, sizeof(d1));
    memset(&d2, 0, sizeof(d2));

    d1.plane.a = 0; d1.plane.b = 0; d1.plane.c = 1; d1.plane.d = -5.0;
    d2.plane.a = 0; d2.plane.b = 0; d2.plane.c = 1; d2.plane.d = -5.0 + 1e-14;

    int8_t match_inv = 0;
    bool eq = alea_primitives_equal(ALEA_PRIMITIVE_PLANE, &d1,
                                    ALEA_PRIMITIVE_PLANE, &d2,
                                    &cfg, &match_inv);
    ASSERT_TRUE(eq);
}

TEST(tol_planes_not_equal) {
    alea_config_t cfg = ALEA_CONFIG_DEFAULT;
    alea_primitive_data_t d1, d2;
    memset(&d1, 0, sizeof(d1));
    memset(&d2, 0, sizeof(d2));

    d1.plane.a = 0; d1.plane.b = 0; d1.plane.c = 1; d1.plane.d = -5.0;
    d2.plane.a = 0; d2.plane.b = 0; d2.plane.c = 1; d2.plane.d = -6.0;

    int8_t match_inv = 0;
    bool eq = alea_primitives_equal(ALEA_PRIMITIVE_PLANE, &d1,
                                    ALEA_PRIMITIVE_PLANE, &d2,
                                    &cfg, &match_inv);
    ASSERT_FALSE(eq);
}

TEST(tol_sphere_equal) {
    alea_config_t cfg = ALEA_CONFIG_DEFAULT;
    alea_primitive_data_t d1, d2;
    memset(&d1, 0, sizeof(d1));
    memset(&d2, 0, sizeof(d2));

    d1.sphere.center_x = 0; d1.sphere.center_y = 0; d1.sphere.center_z = 0;
    d1.sphere.radius = 5.0;

    d2.sphere.center_x = 1e-15; d2.sphere.center_y = 0; d2.sphere.center_z = 0;
    d2.sphere.radius = 5.0 + 1e-14;

    int8_t match_inv = 0;
    bool eq = alea_primitives_equal(ALEA_PRIMITIVE_SPHERE, &d1,
                                    ALEA_PRIMITIVE_SPHERE, &d2,
                                    &cfg, &match_inv);
    ASSERT_TRUE(eq);
}

TEST(tol_cylinder_equal) {
    alea_config_t cfg = ALEA_CONFIG_DEFAULT;
    alea_primitive_data_t d1, d2;
    memset(&d1, 0, sizeof(d1));
    memset(&d2, 0, sizeof(d2));

    d1.cyl_z.center_x = 0; d1.cyl_z.center_y = 0; d1.cyl_z.radius = 3.0;
    d2.cyl_z.center_x = 1e-15; d2.cyl_z.center_y = 0; d2.cyl_z.radius = 3.0;

    int8_t match_inv = 0;
    bool eq = alea_primitives_equal(ALEA_PRIMITIVE_CYLINDER_Z, &d1,
                                    ALEA_PRIMITIVE_CYLINDER_Z, &d2,
                                    &cfg, &match_inv);
    ASSERT_TRUE(eq);
}

TEST(tol_canonical_form) {
    /* Canonicalization: first non-zero coeff should be positive */
    alea_primitive_data_t d;
    memset(&d, 0, sizeof(d));
    d.plane.a = -1; d.plane.b = 0; d.plane.c = 0; d.plane.d = 5;

    int8_t inverted = 0;
    alea_canonicalize_primitive(ALEA_PRIMITIVE_PLANE, &d, &inverted);

    /* Should flip to: a=1, d=-5, inverted=1 */
    ASSERT(d.plane.a > 0);
    ASSERT_EQ(inverted, 1);
}

TEST(tol_dedup_creates_once) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    alea_primitive_data_t d;
    memset(&d, 0, sizeof(d));
    d.sphere.center_x = 0; d.sphere.center_y = 0; d.sphere.center_z = 0;
    d.sphere.radius = 5.0;

    int8_t inv1 = 0, inv2 = 0;
    alea_primitive_id_t id1 = alea_get_or_create_primitive(sys, ALEA_PRIMITIVE_SPHERE, &d, &inv1);
    alea_primitive_id_t id2 = alea_get_or_create_primitive(sys, ALEA_PRIMITIVE_SPHERE, &d, &inv2);

    /* Same geometry → same primitive ID */
    ASSERT_EQ(id1, id2);

    alea_destroy(sys);
}

TEST(tol_dedup_inverted) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Plane with positive normal */
    alea_primitive_data_t d1;
    memset(&d1, 0, sizeof(d1));
    d1.plane.a = 0; d1.plane.b = 0; d1.plane.c = 1; d1.plane.d = -5;

    /* Same plane with negative normal */
    alea_primitive_data_t d2;
    memset(&d2, 0, sizeof(d2));
    d2.plane.a = 0; d2.plane.b = 0; d2.plane.c = -1; d2.plane.d = 5;

    int8_t inv1 = 0, inv2 = 0;
    alea_primitive_id_t id1 = alea_get_or_create_primitive(sys, ALEA_PRIMITIVE_PLANE, &d1, &inv1);
    alea_primitive_id_t id2 = alea_get_or_create_primitive(sys, ALEA_PRIMITIVE_PLANE, &d2, &inv2);

    /* Should deduplicate to same primitive with opposite inverted flag */
    ASSERT_EQ(id1, id2);
    ASSERT_NE(inv1, inv2);

    alea_destroy(sys);
}

TEST(tol_near_coincident) {
    /* Two planes that are close but distinctly different */
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    alea_primitive_data_t d1, d2;
    memset(&d1, 0, sizeof(d1));
    memset(&d2, 0, sizeof(d2));

    d1.plane.a = 0; d1.plane.b = 0; d1.plane.c = 1; d1.plane.d = -5.0;
    d2.plane.a = 0; d2.plane.b = 0; d2.plane.c = 1; d2.plane.d = -5.1;

    int8_t inv1 = 0, inv2 = 0;
    alea_primitive_id_t id1 = alea_get_or_create_primitive(sys, ALEA_PRIMITIVE_PLANE, &d1, &inv1);
    alea_primitive_id_t id2 = alea_get_or_create_primitive(sys, ALEA_PRIMITIVE_PLANE, &d2, &inv2);

    /* 0.1 difference → should NOT deduplicate */
    ASSERT_NE(id1, id2);

    alea_destroy(sys);
}

TEST_MAIN()
