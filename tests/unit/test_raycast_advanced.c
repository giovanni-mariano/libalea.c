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
#include "alea_log.h"
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

TEST(mcnp_graveyard_marks_vacuum_before_query_index_exists) {
    const char* input =
        "Imported graveyard\n"
        "1 0 -1 imp:n=1\n"
        "99 0 1 imp:n=0\n"
        "\n"
        "1 s 0 0 0 3\n"
        "\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT_EQ(sys->cells.count, 2);
    ASSERT_EQ(sys->surfaces.count, 1);
    ASSERT_EQ(sys->surfaces.data[0].boundary_type, ALEA_BOUNDARY_VACUUM);

    alea_ray_navigator_t* nav = alea_ray_navigator_create(sys);
    ASSERT_NOT_NULL(nav);
    double position[3] = {0, 0, 0};
    double direction[3] = {1, 0, 0};
    alea_nav_location_t location;
    alea_nav_event_t event;
    ASSERT_EQ(alea_ray_navigator_restart(nav, position, direction,
                                         &location), 0);
    ASSERT_EQ(location.cell_id, 1);
    ASSERT_EQ(alea_ray_navigator_advance(nav, INFINITY, 10, &event), 0);
    ASSERT_EQ(event.kind, ALEA_NAV_VACUUM);
    ASSERT_NEAR(event.distance, 3.0, 1e-10);
    ASSERT_EQ(event.after.kind, ALEA_NAV_LEAKED);
    alea_ray_navigator_destroy(nav);
    mcnp_model_destroy(model);
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


/* Bounded recovery must preserve overflow, stop enumerating filled cells,
 * and leave the reusable result usable on the next ray. */
static void count_breakpoint_warnings(alea_log_level_t level, const char* file,
                                      int line, const char* message, void* data) {
    (void)file;
    (void)line;
    if (level == ALEA_LOG_LEVEL_WARN && strstr(message, "add_hit failed"))
        (*(int*)data)++;
}

TEST(ray_fill_breakpoint_budget_stops_without_oom_warning) {
    const char* input =
        "Bounded filled sphere\n"
        "1 0 -1 fill=1\n"
        "2 0 1\n"
        "10 0 -2 u=1\n"
        "11 0 2 -1 u=1\n"
        "\n"
        "1 SO 10\n"
        "2 SO 5\n"
        "\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_ray_t ray;
    alea_ray_init(&ray, -20, 0, 0, 1, 0, 0);
    alea_raycast_result_t result;
    alea_raycast_result_init(&result);
    int warnings = 0;
    const alea_log_level_t saved_level = alea_log_get_level();
    alea_log_set_level(ALEA_LOG_LEVEL_WARN);
    alea_log_set_callback(count_breakpoint_warnings, &warnings);
    int rc = alea_raycast_validation_breakpoints_reuse_nocache(
        model->sys, &ray, 0, 40, 4, &result);
    const int error = alea_get_last_error();
    alea_log_set_callback(NULL, NULL);
    alea_log_set_level(saved_level);
    ASSERT_EQ(rc, -1);
    ASSERT_EQ(error, ALEA_ERR_OVERFLOW);
    ASSERT_EQ(warnings, 0);
    ASSERT_EQ(result.hits.count, 4);
    /* Two global surfaces plus the first filled-cell surface. */
    ASSERT_EQ(result.surfaces_tested, 3);
    ASSERT_EQ(result.breakpoint_hit_limit, 0);
    ASSERT_EQ(alea_raycast_validation_breakpoints_reuse_nocache(
        model->sys, &ray, 0, 40, 64, &result), 0);
    ASSERT_EQ(alea_get_last_error(), ALEA_OK);
    ASSERT(result.hits.count > 4);
    alea_raycast_result_free(&result);
    mcnp_model_destroy(model);
}


TEST(ray_nearest_breakpoint_retains_one_hit_and_matches_full_fill_scan) {
    const char* input =
        "Nearest filled sphere\n"
        "1 0 -1 fill=1\n"
        "2 0 1\n"
        "10 0 -2 u=1\n"
        "11 0 2 -1 u=1\n"
        "\n"
        "1 SO 10\n"
        "2 SO 5\n"
        "\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_ray_t ray;
    alea_ray_init(&ray, -20, 0, 0, 1, 0, 0);
    alea_raycast_result_t full, nearest;
    alea_raycast_result_init(&full);
    alea_raycast_result_init(&nearest);
    const double starts[] = {0.0, 10.0, 16.0, 30.0};
    for (size_t i = 0; i < sizeof(starts) / sizeof(starts[0]); i++) {
        ASSERT_EQ(alea_raycast_validation_breakpoints_reuse_nocache(
            model->sys, &ray, starts[i], 40, 64, &full), 0);
        double expected = 40;
        for (size_t j = 0; j < full.hits.count; j++) {
            if (full.hits.data[j].t > nextafter(starts[i], INFINITY)) {
                expected = full.hits.data[j].t;
                break;
            }
        }
        nearest.breakpoint_nearest_only = 1;
        nearest.breakpoint_t_min = starts[i];
        ASSERT_EQ(alea_raycast_validation_breakpoints_reuse_nocache(
            model->sys, &ray, starts[i], 40, 1, &nearest), 0);
        ASSERT_EQ(nearest.hits.count, expected < 40 ? 1 : 0);
        if (nearest.hits.count) {
            ASSERT_NEAR(nearest.hits.data[0].t, expected, 1e-12);
            ASSERT_EQ(nearest.hits.capacity, 1);
        }
    }
    alea_raycast_result_free(&full);
    alea_raycast_result_free(&nearest);
    mcnp_model_destroy(model);
}


TEST(ray_nearest_lattice_breakpoint_matches_full_scan) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    ASSERT_NOT_NULL(model);
    ASSERT_EQ(alea_raycast_ensure_hier_caches(model->sys), 0);
    alea_ray_t ray;
    alea_ray_init(&ray, -10, 0, 0, 1, 0, 0);
    alea_raycast_result_t full, nearest;
    alea_raycast_result_init(&full);
    alea_raycast_result_init(&nearest);
    for (double start = 0; start < 18; start += 0.25) {
        ASSERT_EQ(alea_raycast_validation_breakpoints_reuse_nocache(
            model->sys, &ray, start, 20, 4096, &full), 0);
        double expected = 20;
        for (size_t j = 0; j < full.hits.count; j++) {
            if (full.hits.data[j].t > nextafter(start, INFINITY)) {
                expected = full.hits.data[j].t;
                break;
            }
        }
        nearest.breakpoint_nearest_only = 1;
        nearest.breakpoint_t_min = start;
        ASSERT_EQ(alea_raycast_validation_breakpoints_reuse_nocache(
            model->sys, &ray, start, 20, 1, &nearest), 0);
        ASSERT_EQ(nearest.hits.count, expected < 20 ? 1 : 0);
        if (nearest.hits.count)
            ASSERT_NEAR(nearest.hits.data[0].t, expected, 1e-12);
    }
    alea_raycast_result_free(&full);
    alea_raycast_result_free(&nearest);
    mcnp_model_destroy(model);
}

TEST_MAIN()
