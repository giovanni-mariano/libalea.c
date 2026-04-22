// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_mcnp.c - MCNP parsing, conversion, and export tests
 */

#include "alea_test.h"
#include "alea.h"
#include "core/alea_system.h"
#include "core/alea_export.h"
#include "core/alea_universe.h"
#include "alea_mcnp.h"
#include "alea_openmc.h"
#include "raycast/raycast.h"

/* ------------------------------------------------------------------------- */
/* TRCL (Cell Transformation) Tests                                          */
/* ------------------------------------------------------------------------- */

static const char* TRCL_INLINE_INPUT =
"Test TRCL\n"
"c Cell at origin with TRCL to translate to (10, 0, 0)\n"
"1 1 -1.0 -1 TRCL=(10 0 0)\n"
"2 0 1\n"
"\n"
"1 SO 5\n"
"\n"
"M1 92235.80c 1.0\n";

static const char* TRCL_TR_INPUT =
"Test TRCL with TR card\n"
"c Cell at origin with TRCL=1 referencing TR1\n"
"1 1 -1.0 -1 TRCL=1\n"
"2 0 1\n"
"\n"
"1 SO 5\n"
"\n"
"TR1 20 0 0\n"
"M1 92235.80c 1.0\n";

TEST(trcl_inline_translation) {
    FILE* f = fopen("test_trcl_inline_tmp.mcnp", "w");
    ASSERT_NOT_NULL(f);
    fprintf(f, "%s", TRCL_INLINE_INPUT);
    fclose(f);

    mcnp_model_t* model = mcnp_load("test_trcl_inline_tmp.mcnp");
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;

    /* Origin should be void (sphere translated away) */
    int origin_mat = alea_material_at(sys, 0, 0, 0);
    ASSERT_EQ(origin_mat, 0);

    /* New center (10,0,0) should have material 1 */
    int center_mat = alea_material_at(sys, 10, 0, 0);
    ASSERT_EQ(center_mat, 1);

    /* Point inside translated sphere should have material 1 */
    int near_mat = alea_material_at(sys, 10, 4, 0);
    ASSERT_EQ(near_mat, 1);

    /* Point outside translated sphere should be void */
    int far_mat = alea_material_at(sys, 10, 6, 0);
    ASSERT_EQ(far_mat, 0);

    mcnp_model_destroy(model);
    remove("test_trcl_inline_tmp.mcnp");
}

TEST(trcl_with_tr_card) {
    FILE* f = fopen("test_trcl_tr_tmp.mcnp", "w");
    ASSERT_NOT_NULL(f);
    fprintf(f, "%s", TRCL_TR_INPUT);
    fclose(f);

    mcnp_model_t* model = mcnp_load("test_trcl_tr_tmp.mcnp");
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;

    /* Origin should be void (sphere translated away) */
    int origin_mat = alea_material_at(sys, 0, 0, 0);
    ASSERT_EQ(origin_mat, 0);

    /* (20,0,0) should have material 1 */
    int center_mat = alea_material_at(sys, 20, 0, 0);
    ASSERT_EQ(center_mat, 1);

    mcnp_model_destroy(model);
    remove("test_trcl_tr_tmp.mcnp");
}

/* ------------------------------------------------------------------------- */
/* MCNP File Parsing Tests                                                    */
/* ------------------------------------------------------------------------- */

TEST(parse_simple_sphere) {
    mcnp_model_t* model = mcnp_load("tests/data/simple_sphere.mcnp");
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT(alea_cell_count(sys) > 0);
    mcnp_model_destroy(model);
}

TEST(parse_simple_box) {
    mcnp_model_t* model = mcnp_load("tests/data/simple_box.mcnp");
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT(alea_cell_count(sys) > 0);
    mcnp_model_destroy(model);
}

/* ------------------------------------------------------------------------- */
/* Deduplication Export Tests                                                 */
/* ------------------------------------------------------------------------- */

TEST(dedup_identical_spheres) {
    mcnp_model_t* model = mcnp_load("tests/data/dedup_test.mcnp");
    if (!model) {
        SKIP("Test data file not found");
    }
    alea_system_t* sys = model->sys;

    alea_build_universe_index(sys);

    /* Verify geometry before export */
    int cell_lower = alea_identify_cell_at_point(sys, 0, 0, -5);
    int cell_upper = alea_identify_cell_at_point(sys, 0, 0, 5);
    ASSERT(cell_lower >= 0);
    ASSERT(cell_upper >= 0);

    /* Export with dedup enabled */
    alea_config_t cfg = alea_get_config(sys);
    cfg.dedup = true;
    alea_set_config(sys, &cfg);
    int rc = mcnp_export_system(sys, "test_dedup_tmp.mcnp");
    ASSERT_EQ(rc, 0);

    mcnp_model_destroy(model);

    /* Roundtrip: load exported file and verify */
    mcnp_model_t* model2 = mcnp_load("test_dedup_tmp.mcnp");
    ASSERT_NOT_NULL(model2);
    alea_system_t* sys2 = model2->sys;

    alea_build_universe_index(sys2);
    cell_lower = alea_identify_cell_at_point(sys2, 0, 0, -5);
    cell_upper = alea_identify_cell_at_point(sys2, 0, 0, 5);
    ASSERT(cell_lower >= 0);
    ASSERT(cell_upper >= 0);

    mcnp_model_destroy(model2);
    remove("test_dedup_tmp.mcnp");
}

TEST(dedup_opposite_planes) {
    mcnp_model_t* model = mcnp_load("tests/data/dedup_opposite_signs.mcnp");
    if (!model) {
        SKIP("Test data file not found");
    }
    alea_system_t* sys = model->sys;

    alea_build_universe_index(sys);

    int cell_left = alea_identify_cell_at_point(sys, 0, 0, 0);
    int cell_right = alea_identify_cell_at_point(sys, 10, 0, 0);
    ASSERT(cell_left >= 0);
    ASSERT(cell_right >= 0);

    alea_config_t cfg = alea_get_config(sys);
    cfg.dedup = true;
    alea_set_config(sys, &cfg);
    int rc = mcnp_export_system(sys, "test_dedup_opposite_tmp.mcnp");
    ASSERT_EQ(rc, 0);

    mcnp_model_destroy(model);

    /* Roundtrip verification */
    mcnp_model_t* model2 = mcnp_load("test_dedup_opposite_tmp.mcnp");
    ASSERT_NOT_NULL(model2);
    alea_system_t* sys2 = model2->sys;

    alea_build_universe_index(sys2);
    cell_left = alea_identify_cell_at_point(sys2, 0, 0, 0);
    cell_right = alea_identify_cell_at_point(sys2, 10, 0, 0);
    ASSERT(cell_left >= 0);
    ASSERT(cell_right >= 0);

    mcnp_model_destroy(model2);
    remove("test_dedup_opposite_tmp.mcnp");
}

/* ------------------------------------------------------------------------- */
/* Negated plane roundtrip test                                               */
/*                                                                            */
/* Surface 1 is "P -1 0 0 -5" (x = 5 with negated normal).                  */
/* Canonicalization flips to (1,0,0,-5) with inverted=1.                     */
/* The MCNP exporter un-canonicalizes surface coefficients, so the cell      */
/* sense must be computed relative to the un-canonicalized surface.           */
/* This test verifies materials are correct after roundtrip (non-dedup).     */
/* ------------------------------------------------------------------------- */

TEST(negated_plane_roundtrip) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_negated_plane.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    alea_build_universe_index(sys);

    /* Verify original parse: x=0 should be in cell 1 (mat 1, x < 5) */
    int mat_left = alea_material_at(sys, 0, 0, 0);
    ASSERT_EQ(mat_left, 1);

    /* x=10 should be in cell 2 (mat 2, x > 5) */
    int mat_right = alea_material_at(sys, 10, 0, 0);
    ASSERT_EQ(mat_right, 2);

    /* Export without dedup (exercises un-canonicalization of negated plane) */
    int rc = mcnp_export_system(sys, "test_negated_plane_tmp.mcnp");
    ASSERT_EQ(rc, 0);
    mcnp_model_destroy(model);

    /* Roundtrip: re-parse the exported file */
    mcnp_model_t* model2 = mcnp_load("test_negated_plane_tmp.mcnp");
    ASSERT_NOT_NULL(model2);
    alea_system_t* sys2 = model2->sys;

    alea_build_universe_index(sys2);

    /* Materials must be preserved after roundtrip */
    int mat_left2 = alea_material_at(sys2, 0, 0, 0);
    ASSERT_MSG(mat_left2 == 1,
        "x=0 should have mat 1 (x < 5) after roundtrip");

    int mat_right2 = alea_material_at(sys2, 10, 0, 0);
    ASSERT_MSG(mat_right2 == 2,
        "x=10 should have mat 2 (x > 5) after roundtrip");

    mcnp_model_destroy(model2);
    remove("test_negated_plane_tmp.mcnp");
}

/* Also test with dedup enabled */
TEST(negated_plane_roundtrip_dedup) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_negated_plane.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    alea_build_universe_index(sys);

    /* Export with dedup */
    alea_config_t cfg = alea_get_config(sys);
    cfg.dedup = true;
    alea_set_config(sys, &cfg);
    int rc = mcnp_export_system(sys, "test_negated_plane_dedup_tmp.mcnp");
    ASSERT_EQ(rc, 0);
    mcnp_model_destroy(model);

    /* Roundtrip verification */
    mcnp_model_t* model2 = mcnp_load("test_negated_plane_dedup_tmp.mcnp");
    ASSERT_NOT_NULL(model2);
    alea_system_t* sys2 = model2->sys;

    alea_build_universe_index(sys2);

    int mat_left = alea_material_at(sys2, 0, 0, 0);
    ASSERT_MSG(mat_left == 1,
        "x=0 should have mat 1 (x < 5) after dedup roundtrip");

    int mat_right = alea_material_at(sys2, 10, 0, 0);
    ASSERT_MSG(mat_right == 2,
        "x=10 should have mat 2 (x > 5) after dedup roundtrip");

    mcnp_model_destroy(model2);
    remove("test_negated_plane_dedup_tmp.mcnp");
}

/* ------------------------------------------------------------------------- */
/* Lattice Tests                                                              */
/* ------------------------------------------------------------------------- */

TEST(lattice_mcnp_parse) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice.mcnp");
    if (!model) {
        SKIP("Test data file not found");
    }
    alea_system_t* sys = model->sys;

    ASSERT(alea_cell_count(sys) > 0);
    mcnp_model_destroy(model);
}

/*
 * Lattice point query test.
 *
 * 3x3 rectangular lattice, pitch 2x2x2 (cell bounded by ±1 planes).
 * Universe 1: mat 1 inside cz 0.3, mat 2 outside.
 * Universe 3: mat 3 inside cz 0.3, mat 4 outside.
 * Checkerboard fill: 1 3 1 / 3 1 3 / 1 3 1
 */
TEST(lattice_eval_point_query) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    alea_build_universe_index(sys);

    int cell_id, material;

    /* Element (0,0,0) center = origin, universe 1 */
    ASSERT_EQ(alea_find_cell_lazy(sys, 0, 0, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 1);

    /* Element (0,0,0), outside cylinder → material 2 */
    ASSERT_EQ(alea_find_cell_lazy(sys, 0.5, 0.5, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 2);

    /* Element (1,0,0) center = (2,0,0), universe 3 */
    ASSERT_EQ(alea_find_cell_lazy(sys, 2, 0, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 3);

    /* Element (1,0,0), outside cylinder → material 4 */
    ASSERT_EQ(alea_find_cell_lazy(sys, 2.5, 0.5, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 4);

    /* Element (2,0,0) center = (4,0,0), universe 1 again */
    ASSERT_EQ(alea_find_cell_lazy(sys, 4, 0, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 1);

    /* Element (0,1,0) center = (0,2,0), universe 3 */
    ASSERT_EQ(alea_find_cell_lazy(sys, 0, 2, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 3);

    /* Outside lattice bounds → void */
    ASSERT_EQ(alea_find_cell_lazy(sys, 10, 0, 0, &cell_id, &material, NULL), -1);

    mcnp_model_destroy(model);
}

TEST(lattice_public_point_query_matches_deepest_hit) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    alea_build_universe_index(sys);

    int cell_id = -1;
    int material = -1;

    ASSERT_EQ(alea_material_at(sys, 2, 0, 0), 3);
    ASSERT_EQ(alea_find_cell_at(sys, 2, 0, 0, &cell_id, &material), 0);
    ASSERT_EQ(cell_id, 3);
    ASSERT_EQ(material, 3);

    ASSERT_EQ(alea_material_at(sys, 2.5, 0.5, 0), 4);
    ASSERT_EQ(alea_find_cell_at(sys, 2.5, 0.5, 0, &cell_id, &material), 0);
    ASSERT_EQ(cell_id, 4);
    ASSERT_EQ(material, 4);

    mcnp_model_destroy(model);
}

/*
 * Hex lattice point query via OpenMC.
 *
 * 3x3 hex lattice, pitch 2.0, flat-top orientation.
 * Basis: e1=(2,0), e2=(1, sqrt(3)).
 * Element (0,0) at origin → univ 1, (1,0) at (2,0) → univ 3, etc.
 */
TEST(lattice_hex_eval) {
    openmc_model_t* omc = openmc_load("tests/data/openmc_hex_lattice.xml");
    if (!omc) SKIP("Test data file not found");
    alea_system_t* sys = omc->sys;

    alea_build_universe_index(sys);

    int cell_id, material;

    /* Element (0,0): center at origin, universe 1 → mat 1 inside cyl */
    ASSERT_EQ(alea_find_cell_lazy(sys, 0, 0, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 1);

    /* Element (1,0): center at (2, 0), universe 3 → mat 3 inside cyl */
    ASSERT_EQ(alea_find_cell_lazy(sys, 2, 0, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 3);

    /* Element (0,1): center at (1, sqrt(3)≈1.732), universe 3 → mat 3 */
    ASSERT_EQ(alea_find_cell_lazy(sys, 1.0, 1.732, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 3);

    /* Element (1,1): center at (3, sqrt(3)), universe 1 → mat 1 */
    ASSERT_EQ(alea_find_cell_lazy(sys, 3.0, 1.732, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 1);

    /* Outside cylinder in element (0,0) → mat 2 */
    ASSERT_EQ(alea_find_cell_lazy(sys, 0.5, 0.5, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 2);

    openmc_model_destroy(omc);
}

/*
 * Raycast through a rectangular lattice.
 *
 * Ray along x-axis at y=0, z=0 through a 3x3 lattice (pitch=2).
 * Elements (0,0), (1,0), (2,0) have universes 1, 3, 1.
 * Each has a cylinder r=0.3 at the center.
 * Expected: mat2, mat1, mat2, mat4, mat3, mat4, mat2, mat1, mat2
 */
TEST(lattice_raycast) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    alea_build_universe_index(sys);

    alea_raycast_result_t result;
    alea_raycast_result_init(&result);

    /* Ray along +x through all 3 lattice elements, y=z=0 */
    int rc = alea_raycast(sys, -1.5, 0, 0, 1, 0, 0, 7.0, &result);
    ASSERT_EQ(rc, 0);

    /* Debug: print what we got */
    printf("  hits=%zu segments=%zu\n", result.hits.count, result.segments.count);
    for (size_t i = 0; i < result.segments.count; i++) {
        printf("  seg[%zu]: t=[%.4f, %.4f] cell=%d mat=%d\n",
               i, result.segments.data[i].t_enter, result.segments.data[i].t_exit,
               result.segments.data[i].cell_id, result.segments.data[i].material_id);
    }
    for (size_t i = 0; i < result.hits.count; i++) {
        printf("  hit[%zu]: t=%.6f surf=%d\n",
               i, result.hits.data[i].t, result.hits.data[i].surface_id);
    }

    /* We should have segments with correct materials */
    ASSERT(result.segments.count >= 5);

    /* Verify path through cylinder intersections:
     * Cylinder r=0.3 at each element center (0, 2, 4).
     * Entry at x=-0.3 (t=1.2), exit at x=0.3 (t=1.8) for element 0
     * Entry at x=1.7 (t=3.2), exit at x=2.3 (t=3.8) for element 1
     * Entry at x=3.7 (t=5.2), exit at x=4.3 (t=5.8) for element 2
     */
    int found_mat1 = 0, found_mat3 = 0;
    for (size_t i = 0; i < result.segments.count; i++) {
        if (result.segments.data[i].material_id == 1) found_mat1++;
        if (result.segments.data[i].material_id == 3) found_mat3++;
    }

    /* Should find material 1 (univ 1 cylinder) and material 3 (univ 3 cylinder) */
    ASSERT(found_mat1 >= 2);  /* Elements (0,0) and (2,0) */
    ASSERT(found_mat3 >= 1);  /* Element (1,0) */

    alea_raycast_result_free(&result);
    mcnp_model_destroy(model);
}

/*
 * Raycast through a hex lattice (via OpenMC input).
 *
 * Ray along x-axis at y=0, z=0 through a 3x3 hex lattice (pitch=2).
 * Element (0,0) at origin → univ 1, element (1,0) at (2,0) → univ 3.
 * Each has a cylinder r=0.3.
 * Expected: outside cyl (univ 1), inside cyl (univ 1), outside,
 *           boundary, outside cyl (univ 3), inside cyl (univ 3), ...
 */
TEST(lattice_hex_raycast) {
    openmc_model_t* omc = openmc_load("tests/data/openmc_hex_lattice.xml");
    if (!omc) SKIP("Test data file not found");
    alea_system_t* sys = omc->sys;

    alea_build_universe_index(sys);

    alea_raycast_result_t result;
    alea_raycast_result_init(&result);

    /* Ray along +x at y=0, z=0 */
    int rc = alea_raycast(sys, -0.5, 0, 0, 1, 0, 0, 5.0, &result);
    ASSERT_EQ(rc, 0);

    printf("  hex hits=%zu segments=%zu\n",
           result.hits.count, result.segments.count);
    for (size_t i = 0; i < result.segments.count; i++) {
        printf("  seg[%zu]: t=[%.4f, %.4f] cell=%d mat=%d\n",
               i, result.segments.data[i].t_enter, result.segments.data[i].t_exit,
               result.segments.data[i].cell_id, result.segments.data[i].material_id);
    }

    /* Should find both material 1 (univ 1 cyl) and material 3 (univ 3 cyl) */
    int found_mat1 = 0, found_mat3 = 0;
    for (size_t i = 0; i < result.segments.count; i++) {
        if (result.segments.data[i].material_id == 1) found_mat1++;
        if (result.segments.data[i].material_id == 3) found_mat3++;
    }
    ASSERT(found_mat1 >= 1);
    ASSERT(found_mat3 >= 1);

    alea_raycast_result_free(&result);
    openmc_model_destroy(omc);
}

/*
 * Nested lattice (lattice of lattices) point query.
 *
 * Outer: 2x1 rectangular lattice, pitch 2 (cell 100, fill=0:1 0:0 0:0).
 *   Element (0,0,0) center at (0,0,0) → universe 10
 *   Element (1,0,0) center at (2,0,0) → universe 10
 *
 * Inner (universe 10): 2x2 rectangular lattice, pitch 1 (cell 200,
 *   fill=-1:0 -1:0 0:0).  Covers [-1,1]x[-1,1] in element-local coords.
 *   fill[0]=1 (-1,-1) center (-0.5,-0.5)
 *   fill[1]=3 (-1, 0) center (-0.5, 0.5)
 *   fill[2]=3 ( 0,-1) center ( 0.5,-0.5)
 *   fill[3]=1 ( 0, 0) center ( 0.5, 0.5)
 *
 * Universe 1: mat 1 inside cz 0.2, mat 2 outside.
 * Universe 3: mat 3 inside cz 0.2, mat 4 outside.
 */
TEST(lattice_nested) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_nested_lattice.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    alea_build_universe_index(sys);

    int cell_id, material;

    /* ---- Outer element (0,0,0), inner element (-1,-1) → univ 1 ---- */
    /* Center of inner element at local (-0.5,-0.5) = global (-0.5,-0.5) */
    ASSERT_EQ(alea_find_cell_lazy(sys, -0.5, -0.5, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 1);   /* inside cylinder */

    /* Same element, outside cylinder */
    ASSERT_EQ(alea_find_cell_lazy(sys, -0.8, -0.8, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 2);

    /* ---- Outer element (0,0,0), inner element (0,-1) → univ 3 ---- */
    /* Center at local (0.5,-0.5) = global (0.5,-0.5) */
    ASSERT_EQ(alea_find_cell_lazy(sys, 0.5, -0.5, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 3);

    /* ---- Outer element (0,0,0), inner element (-1,0) → univ 3 ---- */
    /* Center at local (-0.5,0.5) = global (-0.5,0.5) */
    ASSERT_EQ(alea_find_cell_lazy(sys, -0.5, 0.5, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 3);

    /* ---- Outer element (0,0,0), inner element (0,0) → univ 1 ---- */
    /* Center at local (0.5,0.5) = global (0.5,0.5) */
    ASSERT_EQ(alea_find_cell_lazy(sys, 0.5, 0.5, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 1);

    /* ---- Outer element (1,0,0), inner element (-1,-1) → univ 1 ---- */
    /* Outer center at (2,0,0). Inner center at local (-0.5,-0.5) = global (1.5,-0.5) */
    ASSERT_EQ(alea_find_cell_lazy(sys, 1.5, -0.5, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 1);

    /* ---- Outer element (1,0,0), inner element (0,-1) → univ 3 ---- */
    /* Global (2.5,-0.5) */
    ASSERT_EQ(alea_find_cell_lazy(sys, 2.5, -0.5, 0, &cell_id, &material, NULL), 0);
    ASSERT_EQ(material, 3);

    /* ---- Outside lattice bounds → void ---- */
    ASSERT_EQ(alea_find_cell_lazy(sys, 10, 0, 0, &cell_id, &material, NULL), -1);

    mcnp_model_destroy(model);
}

TEST_MAIN()
