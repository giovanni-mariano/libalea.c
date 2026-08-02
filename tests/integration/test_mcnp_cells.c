// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_mcnp_cells.c - MCNP cell card tests
 *
 * Tests cell parameters, LIKE BUT, and geometry expressions through parse
 * and roundtrip verification.
 */

#include "alea_test.h"
#include "alea.h"
#include "alea_mcnp.h"
#include "alea_slice.h"
#include "alea_raycast.h"
#include "raycast/raycast.h"
#include "core/alea_system.h"
#include "core/alea_export.h"
#include <string.h>
#include <stdio.h>

/* ========================================================================= */
/* Helper: parse and verify                                                  */
/* ========================================================================= */

static mcnp_model_t* parse_mcnp(const char* input) {
    mcnp_model_t* model = mcnp_load_string(input, strlen(input));
    if (model) alea_build_universe_index(model->sys);
    return model;
}

/* ========================================================================= */
/* Cell parameter tests                                                      */
/* ========================================================================= */

TEST(cell_void) {
    const char* input =
        "Test void\n"
        "1 0 -1\n"
        "\n"
        "1 SO 5.0\n"
        "\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_cell_count(sys), 1);
    /* Void cell: material = 0 */
    alea_cell_info_t info;
    ASSERT_EQ(alea_cell_get_info(sys, 0, &info), 0);
    ASSERT_EQ(info.material_id, 0);
    mcnp_model_destroy(model);
}

TEST(cell_material_atom_dens) {
    const char* input =
        "Test atom density\n"
        "1 1 0.05 -1\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    alea_cell_info_t info;
    ASSERT_EQ(alea_cell_get_info(sys, 0, &info), 0);
    ASSERT_EQ(info.material_id, 1);
    ASSERT_NEAR(info.density, 0.05, 1e-6);
    ASSERT_FALSE(info.is_mass_density);
    mcnp_model_destroy(model);
}

TEST(cell_material_mass_dens) {
    const char* input =
        "Test mass density\n"
        "1 1 -7.8 -1\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    alea_cell_info_t info;
    ASSERT_EQ(alea_cell_get_info(sys, 0, &info), 0);
    ASSERT_EQ(info.material_id, 1);
    ASSERT_NEAR(info.density, 7.8, 1e-6);
    ASSERT_TRUE(info.is_mass_density);
    mcnp_model_destroy(model);
}

TEST(cell_imp_n) {
    const char* input =
        "Test IMP:N\n"
        "1 1 -1.0 -1 IMP:N=0.5\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    ASSERT_EQ(alea_cell_count(model->sys), 2);
    mcnp_model_destroy(model);
}

TEST(cell_imp_combined) {
    const char* input =
        "Test IMP combined\n"
        "1 1 -1.0 -1 IMP:N,P=1 IMP:E=0\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    ASSERT_EQ(alea_cell_count(model->sys), 2);
    mcnp_model_destroy(model);
}

TEST(cell_vol) {
    const char* input =
        "Test VOL\n"
        "1 1 -1.0 -1 VOL=100.0\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    ASSERT_EQ(alea_cell_count(model->sys), 2);
    mcnp_model_destroy(model);
}

TEST(cell_tmp) {
    const char* input =
        "Test TMP\n"
        "1 1 -1.0 -1 TMP=2.53e-8\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    ASSERT_EQ(alea_cell_count(model->sys), 2);
    mcnp_model_destroy(model);
}

TEST(cell_universe) {
    const char* input =
        "Test universe\n"
        "1 1 -1.0 -1 U=5\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    alea_cell_info_t info;
    ASSERT_EQ(alea_cell_get_info(sys, 0, &info), 0);
    ASSERT_EQ(info.universe_id, 5);
    mcnp_model_destroy(model);
}

TEST(cell_fill_simple) {
    const char* input =
        "Test FILL\n"
        "1 0 -1 FILL=3\n"
        "10 1 -1.0 -2 U=3\n"
        "11 0 2 U=3\n"
        "\n"
        "1 SO 10.0\n"
        "2 SO 3.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    alea_cell_info_t info;
    ASSERT_EQ(alea_cell_find_info(sys, 1, &info), 0);
    ASSERT_EQ(info.fill_universe, 3);
    mcnp_model_destroy(model);
}

TEST(cell_fill_transform) {
    /* FILL with inline translation */
    const char* input =
        "Test FILL transform\n"
        "1 0 -1 FILL=3 (5 0 0)\n"
        "10 1 -1.0 -2 U=3\n"
        "11 0 2 U=3\n"
        "\n"
        "1 SO 20.0\n"
        "2 SO 3.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    /* Verify fill_transform is set on cell 1 */
    alea_cell_info_t info;
    ASSERT_EQ(alea_cell_find_info(sys, 1, &info), 0);
    ASSERT_EQ(info.fill_universe, 3);
    ASSERT(info.fill_transform > 0);
    /* Material 1 at (5,0,0) via lazy universe resolution */
    int cell_id, material;
    ASSERT_EQ(alea_find_cell_lazy(sys, 5.0, 0, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 1);
    mcnp_model_destroy(model);
}

TEST(cell_lat_rect) {
    const char* input =
        "Test LAT=1\n"
        "1 0 -1 2 -3 4 -5 6 LAT=1 FILL=0:0 0:0 0:0 2 U=10\n"
        "10 1 -1.0 -7 U=2\n"
        "11 0 7 U=2\n"
        "100 0 -100 FILL=10\n"
        "200 0 100\n"
        "\n"
        "1 PX 1.0\n"
        "2 PX -1.0\n"
        "3 PY 1.0\n"
        "4 PY -1.0\n"
        "5 PZ 1.0\n"
        "6 PZ -1.0\n"
        "7 SO 0.3\n"
        "100 SO 50.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    alea_cell_info_t info;
    ASSERT_EQ(alea_cell_find_info(sys, 1, &info), 0);
    ASSERT_EQ(info.lat_type, 1);
    ASSERT_FALSE(info.lat_fill_repeating);
    ASSERT_TRUE(info.lat_fill_zero_element_coords);

    const char* tmpfile = "test_finite_lat_rt_tmp.mcnp";
    ASSERT_EQ(mcnp_export(model, tmpfile), 0);
    mcnp_model_t* roundtrip = mcnp_load(tmpfile);
    ASSERT_NOT_NULL(roundtrip);
    ASSERT_EQ(alea_cell_find_info(roundtrip->sys, 1, &info), 0);
    ASSERT_FALSE(info.lat_fill_repeating);
    mcnp_model_destroy(roundtrip);
    remove(tmpfile);

    mcnp_model_destroy(model);
}

TEST(cell_lat_rect_simple_fill_repeats) {
    const char* input =
        "Test repeating LAT=1 simple fill\n"
        "1 0 -1 2 -3 4 -5 6 LAT=1 FILL=2 U=10\n"
        "10 1 -1.0 -7 U=2\n"
        "11 0 7 U=2\n"
        "100 0 -100 FILL=10\n"
        "200 0 100\n"
        "\n"
        "1 PX 1.0\n"
        "2 PX -1.0\n"
        "3 PY 1.0\n"
        "4 PY -1.0\n"
        "5 PZ 1.0\n"
        "6 PZ -1.0\n"
        "7 SO 0.3\n"
        "100 SO 50.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;

    alea_cell_info_t info;
    ASSERT_EQ(alea_cell_find_info(sys, 1, &info), 0);
    ASSERT_TRUE(info.lat_fill_repeating);
    ASSERT_TRUE(info.lat_fill_zero_element_coords);

    /* The simple fill repeats the U=2 sphere in adjacent lattice elements. */
    ASSERT_EQ(alea_material_at(sys, 0.0, 0.0, 0.0), 1);
    ASSERT_EQ(alea_material_at(sys, 2.0, 0.0, 0.0), 1);
    ASSERT_EQ(alea_material_at(sys, -2.0, 0.0, 0.0), 1);
    ASSERT_EQ(alea_material_at(sys, 0.0, 2.0, 0.0), 1);
    ASSERT_EQ(alea_material_at(sys, 0.0, 0.0, 2.0), 1);

    /* The hierarchy-backed resolver used by slice grids must follow the
     * same repeated element, rather than returning the LAT container as an
     * undefined fill outside the fundamental element. */
    alea_cell_hit_t hits[8];
    int hit_count = alea_find_all_cells(sys, 2.0, 0.0, 0.0, hits, 8);
    ASSERT(hit_count >= 3);
    ASSERT_EQ(hits[hit_count - 1].material_id, 1);
    ASSERT_EQ(hits[hit_count - 1].resolution_flags & ALEA_RESOLVE_UNDEFINED_FILL, 0);

    alea_slice_view_t view = {
        .plane = {
            .origin = {0.0, 0.0, 0.0},
            .normal = {0.0, 0.0, 1.0},
            .u_axis = {1.0, 0.0, 0.0},
            .v_axis = {0.0, 1.0, 0.0}
        },
        .u_min = 1.9, .u_max = 2.1,
        .v_min = -0.1, .v_max = 0.1
    };
    int grid_cell = -1, grid_material = -1;
    uint8_t grid_error = 0xff;
    ASSERT_EQ(alea_find_cells_grid(sys, &view, 1, 1, -1,
                                   &grid_cell, &grid_material, &grid_error), 0);
    ASSERT_EQ(grid_material, 1);
    ASSERT_EQ(grid_error, 0);

    alea_raycast_result_t ray_result;
    alea_raycast_result_init(&ray_result);
    ASSERT_EQ(alea_raycast(sys, -4.5, 0.0, 0.0, 1.0, 0.0, 0.0,
                           9.0, &ray_result), 0);
    int repeated_sphere_segments = 0;
    for (size_t i = 0; i < ray_result.segments.count; i++) {
        if (ray_result.segments.data[i].material_id == 1)
            repeated_sphere_segments++;
    }
    ASSERT(repeated_sphere_segments >= 5);
    alea_raycast_result_free(&ray_result);

    const char* tmpfile = "test_repeating_lat_rt_tmp.mcnp";
    ASSERT_EQ(mcnp_export(model, tmpfile), 0);
    mcnp_model_t* roundtrip = mcnp_load(tmpfile);
    ASSERT_NOT_NULL(roundtrip);
    ASSERT_EQ(alea_cell_find_info(roundtrip->sys, 1, &info), 0);
    ASSERT_TRUE(info.lat_fill_repeating);
    mcnp_model_destroy(roundtrip);
    remove(tmpfile);

    mcnp_model_destroy(model);
}

TEST(cell_lat_rect_preserves_zero_element_coordinates) {
    /* MCNP filling-universe surfaces use the coordinate system of lattice
     * element (0,0,0); they are not implicitly recentered at the origin. */
    const char* input =
        "Test non-origin repeating LAT=1\n"
        "1 0 -1 2 -3 4 -5 6 LAT=1 FILL=2 U=10\n"
        "10 1 -1.0 -7 U=2\n"
        "11 0 7 U=2\n"
        "100 0 -100 FILL=10\n"
        "200 0 100\n"
        "\n"
        "1 PX 12.0\n"
        "2 PX 10.0\n"
        "3 PY 1.0\n"
        "4 PY -1.0\n"
        "5 PZ 1.0\n"
        "6 PZ -1.0\n"
        "7 S 11.0 0.0 0.0 0.3\n"
        "100 SO 50.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;

    alea_cell_info_t info;
    ASSERT_EQ(alea_cell_find_info(sys, 1, &info), 0);
    ASSERT_TRUE(info.lat_fill_repeating);
    ASSERT_TRUE(info.lat_fill_zero_element_coords);

    /* Element zero remains at the authored x=11 coordinates.  Neighbouring
     * elements remove only their +/-2 index displacement on descent. */
    ASSERT_EQ(alea_material_at(sys, 11.0, 0.0, 0.0), 1);
    ASSERT_EQ(alea_material_at(sys, 13.0, 0.0, 0.0), 1);
    ASSERT_EQ(alea_material_at(sys, 9.0, 0.0, 0.0), 1);
    ASSERT_EQ(alea_material_at(sys, 12.0, 0.0, 0.0), 0);

    alea_cell_hit_t hits[8];
    int n = alea_find_all_cells(sys, 13.0, 0.0, 0.0, hits, 8);
    ASSERT(n >= 3);
    ASSERT_EQ(hits[n - 1].cell_id, 10);
    ASSERT_NEAR(hits[n - 1].local_x, 11.0, 1e-4);
    ASSERT_EQ(hits[n - 1].resolution_flags, 0);

    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    alea_raycast_result_t ray_result;
    alea_raycast_result_init(&ray_result);
    ASSERT_EQ(alea_raycast(sys, 8.5, 0.0, 0.0, 1.0, 0.0, 0.0,
                           5.0, &ray_result), 0);
    int material_segments = 0;
    for (size_t i = 0; i < ray_result.segments.count; i++) {
        if (ray_result.segments.data[i].material_id == 1)
            material_segments++;
    }
    ASSERT(material_segments >= 3);
    alea_raycast_result_free(&ray_result);

    mcnp_model_destroy(model);
}

TEST(cell_trcl_id) {
    const char* input =
        "Test TRCL ID\n"
        "1 1 -1.0 -1 TRCL=1\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "TR1 10 0 0\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    /* Sphere at (10,0,0) */
    ASSERT_EQ(alea_material_at(sys, 10, 0, 0), 1);
    ASSERT_EQ(alea_material_at(sys, 0, 0, 0), 0);
    mcnp_model_destroy(model);
}

TEST(cell_trcl_inline_trans) {
    const char* input =
        "Test TRCL inline\n"
        "1 1 -1.0 -1 TRCL=(10 20 30)\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_material_at(sys, 10, 20, 30), 1);
    ASSERT_EQ(alea_material_at(sys, 0, 0, 0), 0);
    mcnp_model_destroy(model);
}

TEST(cell_multiple_params) {
    const char* input =
        "Test multiple params\n"
        "1 1 -1.0 -1 IMP:N=1 IMP:P=1 U=2 VOL=50 TMP=2.53e-8\n"
        "2 0 1\n"
        "\n"
        "1 SO 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    alea_cell_info_t info;
    ASSERT_EQ(alea_cell_get_info(sys, 0, &info), 0);
    ASSERT_EQ(info.universe_id, 2);
    ASSERT_EQ(info.material_id, 1);
    mcnp_model_destroy(model);
}

/* ========================================================================= */
/* LIKE BUT tests                                                            */
/* ========================================================================= */

TEST(like_but_basic) {
    const char* input =
        "Test LIKE BUT\n"
        "1 1 -1.0 -1\n"
        "2 LIKE 1 BUT TRCL=(10 0 0)\n"
        "3 0 1 2\n"
        "\n"
        "1 SO 3.0\n"
        "2 S 10.0 0.0 0.0 3.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    /* Verify LIKE BUT created cell 2 with same material */
    alea_cell_info_t info;
    ASSERT_EQ(alea_cell_find_info(sys, 2, &info), 0);
    ASSERT_EQ(info.material_id, 1);
    /* Cell 1 at origin */
    ASSERT_EQ(alea_material_at(sys, 0, 0, 0), 1);
    mcnp_model_destroy(model);
}

TEST(like_but_mat_override) {
    const char* input =
        "Test LIKE BUT MAT\n"
        "1 1 -1.0 -1\n"
        "2 LIKE 1 BUT MAT=2 RHO=-8.0 TRCL=(10 0 0)\n"
        "3 0 1 2\n"
        "\n"
        "1 SO 3.0\n"
        "2 S 10.0 0.0 0.0 3.0\n"
        "\n"
        "M1 92235.80c 1.0\n"
        "M2 26056.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    /* Cell 1 at origin has material 1 */
    ASSERT_EQ(alea_material_at(sys, 0, 0, 0), 1);
    /* Verify LIKE BUT overrode material to 2 on cell 2 */
    alea_cell_info_t info;
    ASSERT_EQ(alea_cell_find_info(sys, 2, &info), 0);
    ASSERT_EQ(info.material_id, 2);
    mcnp_model_destroy(model);
}

TEST(like_but_imp_override) {
    const char* input =
        "Test LIKE BUT IMP\n"
        "1 1 -1.0 -1 IMP:N=1\n"
        "2 LIKE 1 BUT IMP:N=0 TRCL=(10 0 0)\n"
        "3 0 1 2\n"
        "\n"
        "1 SO 3.0\n"
        "2 S 10.0 0.0 0.0 3.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    ASSERT_EQ(alea_cell_count(model->sys), 3);
    mcnp_model_destroy(model);
}

/* ========================================================================= */
/* Geometry expression tests                                                 */
/* ========================================================================= */

TEST(geom_simple_intersection) {
    /* Box region: -1 2 -3 means inside surf 1, outside surf 2, inside surf 3 */
    const char* input =
        "Test intersection\n"
        "1 1 -7.8 -1 2 -3\n"
        "2 0 1 : -2 : 3\n"
        "\n"
        "1 PZ 10.0\n"
        "2 PZ 0.0\n"
        "3 CZ 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    /* Inside: z=5, r=0 → inside cz, below pz 10, above pz 0 */
    ASSERT_EQ(alea_material_at(sys, 0, 0, 5), 1);
    /* Outside: z=15 → above pz 10 */
    ASSERT_EQ(alea_material_at(sys, 0, 0, 15), 0);
    /* Outside: r=7 → outside cz */
    ASSERT_EQ(alea_material_at(sys, 7, 0, 5), 0);
    mcnp_model_destroy(model);
}

TEST(geom_union) {
    /* Union of two spheres */
    const char* input =
        "Test union\n"
        "1 1 -7.8 (-1 : -2)\n"
        "2 0 1 2\n"
        "\n"
        "1 S 0.0 0.0 0.0 3.0\n"
        "2 S 5.0 0.0 0.0 3.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    /* Center of first sphere */
    ASSERT_EQ(alea_material_at(sys, 0, 0, 0), 1);
    /* Center of second sphere */
    ASSERT_EQ(alea_material_at(sys, 5, 0, 0), 1);
    /* Far away */
    ASSERT_EQ(alea_material_at(sys, 20, 0, 0), 0);
    mcnp_model_destroy(model);
}

TEST(geom_complement_cell) {
    /* Cell complement: #2 */
    const char* input =
        "Test complement\n"
        "1 1 -1.0 -1\n"
        "2 2 -1.0 -2\n"
        "3 3 -1.0 #1 #2\n"
        "\n"
        "1 SO 3.0\n"
        "2 S 10.0 0.0 0.0 3.0\n"
        "\n"
        "M1 92235.80c 1.0\n"
        "M2 26056.80c 1.0\n"
        "M3 1001.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    /* Inside cell 1 → mat 1 */
    ASSERT_EQ(alea_material_at(sys, 0, 0, 0), 1);
    /* Inside cell 2 → mat 2 */
    ASSERT_EQ(alea_material_at(sys, 10, 0, 0), 2);
    /* Between them → mat 3 (complement of both) */
    ASSERT_EQ(alea_material_at(sys, 5, 0, 0), 3);
    mcnp_model_destroy(model);
}

TEST(geom_complement_surface) {
    /* Complement of surface expression: #(-1 2) */
    const char* input =
        "Test complement expr\n"
        "1 1 -1.0 -1 2\n"
        "2 0 #(-1 2)\n"
        "\n"
        "1 PZ 10.0\n"
        "2 PZ 0.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    /* Between planes: z=5 → mat 1 */
    ASSERT_EQ(alea_material_at(sys, 0, 0, 5), 1);
    /* Below plane 2: z=-5 → void (complement) */
    ASSERT_EQ(alea_material_at(sys, 0, 0, -5), 0);
    mcnp_model_destroy(model);
}

TEST(geom_nested_parens) {
    /* (-1 2 : -3 4) — union of two intersection terms */
    const char* input =
        "Test nested parens\n"
        "1 1 -1.0 (-1 2 : -3 4)\n"
        "2 0 #1\n"
        "\n"
        "1 PZ 10.0\n"
        "2 PZ 0.0\n"
        "3 PZ 20.0\n"
        "4 PZ 15.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    /* z=5 → in first term (0 < z < 10) */
    ASSERT_EQ(alea_material_at(sys, 0, 0, 5), 1);
    /* z=17 → in second term (15 < z < 20) */
    ASSERT_EQ(alea_material_at(sys, 0, 0, 17), 1);
    /* z=12 → outside both */
    ASSERT_EQ(alea_material_at(sys, 0, 0, 12), 0);
    mcnp_model_destroy(model);
}

TEST(geom_deep_nesting) {
    /* Multiple parenthesis levels */
    const char* input =
        "Test deep nesting\n"
        "1 1 -1.0 ((-1 2) (-3 4))\n"
        "2 0 #1\n"
        "\n"
        "1 PZ 10.0\n"
        "2 PZ 0.0\n"
        "3 CZ 5.0\n"
        "4 PX 0.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    /* Inside all: z in (0,10), r < 5, x > 0 → at (2, 0, 5) */
    ASSERT_EQ(alea_material_at(sys, 2, 0, 5), 1);
    /* Outside r: at (7, 0, 5) */
    ASSERT_EQ(alea_material_at(sys, 7, 0, 5), 0);
    mcnp_model_destroy(model);
}

/* ========================================================================= */
/* Roundtrip test: cells with parameters survive export/re-parse             */
/* ========================================================================= */

TEST(cell_roundtrip_params) {
    const char* input =
        "Test cell params roundtrip\n"
        "1 1 -7.8 -1 2 -3 IMP:N=1\n"
        "2 0 1 : -2 : 3 IMP:N=0\n"
        "\n"
        "1 PZ 10.0\n"
        "2 PZ 0.0\n"
        "3 CZ 5.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);

    const char* tmpfile = "test_cell_rt_tmp.mcnp";
    int rc = mcnp_export(model, tmpfile);
    ASSERT_EQ(rc, 0);
    mcnp_model_destroy(model);

    mcnp_model_t* model2 = mcnp_load(tmpfile);
    ASSERT_NOT_NULL(model2);
    alea_system_t* sys2 = model2->sys;
    alea_build_universe_index(sys2);

    /* Geometry preserved */
    ASSERT_EQ(alea_material_at(sys2, 0, 0, 5), 1);
    ASSERT_EQ(alea_material_at(sys2, 0, 0, 15), 0);

    mcnp_model_destroy(model2);
    remove(tmpfile);
}

TEST(cell_difference_roundtrip) {
    /* A difference operation: sphere with hole */
    const char* input =
        "Test difference\n"
        "1 1 -1.0 -1 2\n"
        "2 0 1 : -2\n"
        "\n"
        "1 SO 5.0\n"
        "2 SO 2.0\n"
        "\n"
        "M1 92235.80c 1.0\n";
    mcnp_model_t* model = parse_mcnp(input);
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;

    /* Between radii: mat 1 */
    ASSERT_EQ(alea_material_at(sys, 3.0, 0, 0), 1);
    /* Inside hole: void */
    ASSERT_EQ(alea_material_at(sys, 0, 0, 0), 0);
    /* Outside: void */
    ASSERT_EQ(alea_material_at(sys, 7.0, 0, 0), 0);

    /* Roundtrip */
    const char* tmpfile = "test_diff_rt_tmp.mcnp";
    int rc = mcnp_export_system(sys, tmpfile);
    ASSERT_EQ(rc, 0);
    mcnp_model_destroy(model);

    mcnp_model_t* model2 = mcnp_load(tmpfile);
    ASSERT_NOT_NULL(model2);
    alea_system_t* sys2 = model2->sys;
    alea_build_universe_index(sys2);
    ASSERT_EQ(alea_material_at(sys2, 3.0, 0, 0), 1);
    ASSERT_EQ(alea_material_at(sys2, 0, 0, 0), 0);
    mcnp_model_destroy(model2);
    remove(tmpfile);
}

TEST_MAIN()
