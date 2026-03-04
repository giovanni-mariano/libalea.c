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
