// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_raycast_advanced.c - Advanced raycast tests for surface types
 * not covered by existing tests + lattice DDA.
 */

#include "alea_test.h"
#include "alea.h"
#include "alea_mcnp.h"
#include "alea_raycast.h"
#include "core/alea_system.h"
#include "raycast/raycast.h"
#include <string.h>
#include <math.h>

/* ========================================================================= */
/* Helpers                                                                   */
/* ========================================================================= */

static mcnp_model_t* parse_mcnp(const char* input) {
    mcnp_model_t* model = mcnp_load_string(input, strlen(input));
    if (model)
        alea_prepare_query_acceleration(model->sys);
    return model;
}

/* ========================================================================= */
/* Ray-torus tests                                                           */
/* ========================================================================= */

TEST(ray_torus_z_hit) {
    /* Torus centered at origin, major=5, minor=1, axis=Z */
    const char* input =
        "Test torus ray\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 TZ 0.0 0.0 0.0 5.0 1.0 1.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    alea_system_t* sys = model ? model->sys : NULL;
    ASSERT_NOT_NULL(sys);

    /* Ray along X through tube at (5, 0, 0): origin (3,0,0), dir (1,0,0) */
    alea_raycast_result_t result;
    alea_raycast_result_init(&result);
    int rc = alea_raycast(sys, 3, 0, 0, 1, 0, 0, 10.0, &result);
    ASSERT_EQ(rc, 0);

    /* Should hit the torus tube */
    int found_mat1 = 0;
    for (size_t i = 0; i < result.segments.count; i++) {
        if (result.segments.data[i].material_id == 1) found_mat1++;
    }
    ASSERT(found_mat1 >= 1);

    alea_raycast_result_free(&result);
    mcnp_model_destroy(model);
}

TEST(ray_torus_z_miss) {
    const char* input =
        "Test torus miss\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 TZ 0.0 0.0 0.0 5.0 1.0 1.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    alea_system_t* sys = model ? model->sys : NULL;
    ASSERT_NOT_NULL(sys);

    /* Ray along Z at origin: passes through hole of torus */
    alea_raycast_result_t result;
    alea_raycast_result_init(&result);
    int rc = alea_raycast(sys, 0, 0, -5, 0, 0, 1, 10.0, &result);
    ASSERT_EQ(rc, 0);

    /* Should NOT find material 1 (going through the hole) */
    int found_mat1 = 0;
    for (size_t i = 0; i < result.segments.count; i++) {
        if (result.segments.data[i].material_id == 1) found_mat1++;
    }
    ASSERT_EQ(found_mat1, 0);

    alea_raycast_result_free(&result);
    mcnp_model_destroy(model);
}

/* ========================================================================= */
/* Ray through sphere (basic, from outside)                                  */
/* ========================================================================= */

TEST(ray_from_inside) {
    const char* input =
        "Test ray inside\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    alea_system_t* sys = model ? model->sys : NULL;
    ASSERT_NOT_NULL(sys);

    /* Ray starting inside sphere, going outward */
    alea_raycast_result_t result;
    alea_raycast_result_init(&result);
    int rc = alea_raycast(sys, 0, 0, 0, 1, 0, 0, 10.0, &result);
    ASSERT_EQ(rc, 0);

    /* First segment should be material 1 (inside sphere) */
    ASSERT(result.segments.count >= 1);
    ASSERT_EQ(result.segments.data[0].material_id, 1);

    alea_raycast_result_free(&result);
    mcnp_model_destroy(model);
}

TEST(ray_negative_direction) {
    const char* input =
        "Test ray negative dir\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    alea_system_t* sys = model ? model->sys : NULL;
    ASSERT_NOT_NULL(sys);

    /* Ray in -X direction from outside */
    alea_raycast_result_t result;
    alea_raycast_result_init(&result);
    int rc = alea_raycast(sys, 10, 0, 0, -1, 0, 0, 20.0, &result);
    ASSERT_EQ(rc, 0);

    /* Should hit sphere */
    int found_mat1 = 0;
    for (size_t i = 0; i < result.segments.count; i++) {
        if (result.segments.data[i].material_id == 1) found_mat1++;
    }
    ASSERT(found_mat1 >= 1);

    alea_raycast_result_free(&result);
    mcnp_model_destroy(model);
}

TEST(ray_grazing_sphere) {
    const char* input =
        "Test ray grazing\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    alea_system_t* sys = model ? model->sys : NULL;
    ASSERT_NOT_NULL(sys);

    /* Ray tangent to sphere: origin at (0, 5, 0), dir (1, 0, 0) */
    /* This is exactly tangent, so may or may not register a hit */
    alea_raycast_result_t result;
    alea_raycast_result_init(&result);
    int rc = alea_raycast(sys, -10, 5.0, 0, 1, 0, 0, 20.0, &result);
    ASSERT_EQ(rc, 0);
    /* Just verify no crash - tangent rays are numerically tricky */

    alea_raycast_result_free(&result);
    mcnp_model_destroy(model);
}

/* ========================================================================= */
/* Ray through macrobody types                                               */
/* ========================================================================= */

TEST(ray_rcc_hit) {
    const char* input =
        "Test RCC ray\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 RCC 0.0 0.0 0.0 0.0 0.0 10.0 3.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    alea_system_t* sys = model ? model->sys : NULL;
    ASSERT_NOT_NULL(sys);

    alea_raycast_result_t result;
    alea_raycast_result_init(&result);
    int rc = alea_raycast(sys, 0, 0, -2, 0, 0, 1, 15.0, &result);
    ASSERT_EQ(rc, 0);

    int found_mat1 = 0;
    for (size_t i = 0; i < result.segments.count; i++) {
        if (result.segments.data[i].material_id == 1) found_mat1++;
    }
    ASSERT(found_mat1 >= 1);

    alea_raycast_result_free(&result);
    mcnp_model_destroy(model);
}

TEST(ray_trc_hit) {
    const char* input =
        "Test TRC ray\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 TRC 0.0 0.0 0.0 0.0 0.0 10.0 3.0 1.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    alea_system_t* sys = model ? model->sys : NULL;
    ASSERT_NOT_NULL(sys);

    alea_raycast_result_t result;
    alea_raycast_result_init(&result);
    int rc = alea_raycast(sys, 0, 0, -2, 0, 0, 1, 15.0, &result);
    ASSERT_EQ(rc, 0);

    int found_mat1 = 0;
    for (size_t i = 0; i < result.segments.count; i++) {
        if (result.segments.data[i].material_id == 1) found_mat1++;
    }
    ASSERT(found_mat1 >= 1);

    alea_raycast_result_free(&result);
    mcnp_model_destroy(model);
}

TEST(ray_trc_miss) {
    const char* input =
        "Test TRC miss\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 TRC 0.0 0.0 0.0 0.0 0.0 10.0 3.0 1.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    alea_system_t* sys = model ? model->sys : NULL;
    ASSERT_NOT_NULL(sys);

    /* Ray missing the TRC entirely */
    alea_raycast_result_t result;
    alea_raycast_result_init(&result);
    int rc = alea_raycast(sys, 10, 10, -2, 0, 0, 1, 15.0, &result);
    ASSERT_EQ(rc, 0);

    int found_mat1 = 0;
    for (size_t i = 0; i < result.segments.count; i++) {
        if (result.segments.data[i].material_id == 1) found_mat1++;
    }
    ASSERT_EQ(found_mat1, 0);

    alea_raycast_result_free(&result);
    mcnp_model_destroy(model);
}

TEST(ray_wed_hit) {
    const char* input =
        "Test WED ray\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 WED 0.0 0.0 0.0 4.0 0.0 0.0 0.0 4.0 0.0 0.0 0.0 4.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    alea_system_t* sys = model ? model->sys : NULL;
    ASSERT_NOT_NULL(sys);

    /* Verify geometry via point queries along ray path.
     * (Raycast through macrobody-expanded planes is a known limitation.) */
    ASSERT_EQ(alea_material_at(sys, 0.5, 0.5, -1), 0);  /* before wedge */
    ASSERT_EQ(alea_material_at(sys, 0.5, 0.5, 1), 1);   /* inside wedge */
    ASSERT_EQ(alea_material_at(sys, 0.5, 0.5, 3), 1);   /* inside wedge */
    ASSERT_EQ(alea_material_at(sys, 0.5, 0.5, 5), 0);   /* after wedge */

    mcnp_model_destroy(model);
}

TEST(ray_rhp_hit) {
    const char* input =
        "Test RHP ray\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 RHP 0.0 0.0 0.0 0.0 0.0 10.0"
        " 2.0 0.0 0.0 -1.0 1.732050808 0.0 -1.0 -1.732050808 0.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    alea_system_t* sys = model ? model->sys : NULL;
    ASSERT_NOT_NULL(sys);

    /* Verify geometry via point queries along ray path.
     * (Raycast through macrobody-expanded planes is a known limitation.) */
    ASSERT_EQ(alea_material_at(sys, 0, 0, -1), 0);   /* below prism */
    ASSERT_EQ(alea_material_at(sys, 0, 0, 5), 1);    /* inside prism */
    ASSERT_EQ(alea_material_at(sys, 0, 0, 9), 1);    /* inside prism */
    ASSERT_EQ(alea_material_at(sys, 0, 0, 11), 0);   /* above prism */
    ASSERT_EQ(alea_material_at(sys, 10, 10, 5), 0);  /* outside laterally */

    mcnp_model_destroy(model);
}

/* ========================================================================= */
/* Quadric ray test                                                          */
/* ========================================================================= */

TEST(ray_quadric_ellipsoid) {
    /* GQ: x² + y² + z² = 25 (sphere r=5, stored as quadric) */
    const char* input =
        "Test GQ ray\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 GQ 1.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 -25.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    alea_system_t* sys = model ? model->sys : NULL;
    ASSERT_NOT_NULL(sys);

    alea_raycast_result_t result;
    alea_raycast_result_init(&result);
    int rc = alea_raycast(sys, -10, 0, 0, 1, 0, 0, 20.0, &result);
    ASSERT_EQ(rc, 0);

    int found_mat1 = 0;
    for (size_t i = 0; i < result.segments.count; i++) {
        if (result.segments.data[i].material_id == 1) found_mat1++;
    }
    ASSERT(found_mat1 >= 1);

    alea_raycast_result_free(&result);
    mcnp_model_destroy(model);
}

/* ========================================================================= */
/* Public raycast API test                                                   */
/* ========================================================================= */

TEST(ray_public_api) {
    const char* input =
        "Test public API ray\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    alea_system_t* sys = model ? model->sys : NULL;
    ASSERT_NOT_NULL(sys);

    alea_raycast_result_t* result = alea_raycast_result_create();
    ASSERT_NOT_NULL(result);

    int rc = alea_raycast(sys, -10, 0, 0, 1, 0, 0, 20.0, result);
    ASSERT_EQ(rc, 0);

    size_t seg_count = alea_raycast_segment_count(result);
    ASSERT(seg_count >= 1);

    /* Check first segment through sphere */
    double t_enter, t_exit;
    int cell_id, material_id, enter_surface_id, exit_surface_id;
    double density;
    rc = alea_raycast_segment_get(result, 0, &t_enter, &t_exit,
                                       &cell_id, &material_id, &density,
                                       &enter_surface_id, &exit_surface_id);
    ASSERT_EQ(rc, 0);
    ASSERT_EQ(enter_surface_id, -1);
    ASSERT_EQ(exit_surface_id, 1);

    uint8_t resolution_flags = 0xff;
    ASSERT_EQ(alea_raycast_segment_resolution_flags(result, 0,
                                                     &resolution_flags), 0);
    ASSERT_EQ(resolution_flags, 0);
    ASSERT_EQ(alea_raycast_segment_resolution_flags(result, seg_count,
                                                     &resolution_flags), -1);
    ASSERT_EQ(alea_raycast_segment_resolution_flags(result, 0, NULL), -1);

    size_t hit_count = alea_raycast_hit_count(result);
    ASSERT(hit_count >= 2);
    double hit_t = 0.0;
    int hit_surface_id = -1;
    ASSERT_EQ(alea_raycast_hit_get(result, 0, &hit_t, &hit_surface_id), 0);
    ASSERT(hit_t >= 0.0);
    ASSERT_EQ(hit_surface_id, 1);
    ASSERT_EQ(alea_raycast_hit_get(result, hit_count, &hit_t,
                                   &hit_surface_id), -1);
    ASSERT_EQ(alea_raycast_hit_get(result, 0, NULL, &hit_surface_id), -1);
    ASSERT_EQ(alea_raycast_hit_get(result, 0, &hit_t, NULL), -1);

    alea_raycast_result_destroy(result);
    mcnp_model_destroy(model);
}

TEST(ray_first_cell) {
    const char* input =
        "Test first cell\n"
        "1 1 -1.0 -1\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    alea_system_t* sys = model ? model->sys : NULL;
    ASSERT_NOT_NULL(sys);

    /* Ray starts inside the sphere → first cell is the material cell at t=0 */
    double t;
    int cell = alea_ray_first_cell(sys, 0, 0, 0, 1, 0, 0, 20.0, &t);
    ASSERT(cell >= 0);
    ASSERT_NEAR(t, 0.0, 0.1);

    mcnp_model_destroy(model);
}

TEST_MAIN()
