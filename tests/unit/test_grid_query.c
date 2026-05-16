// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_grid_query.c - Tests for alea_find_cells_grid
 *
 * Covers:
 *   - Uniform region: every pixel lands in the expected cell
 *   - Two-cell X split: boundary pixel gets the correct cell even when
 *     the previous pixel's hint is stale (tests the verify+re-scan path)
 *   - Three-cell row: two transitions, both corrected
 *   - Overlap detection: ALEA_GRID_OVERLAP flagged for geometry errors
 *   - No false overlaps: clean model produces no error pixels
 *   - Flat/hier parity: simple_fill.mcnp gives identical cell arrays in
 *     both spatial modes
 */

#define ALEA_TEST_IMPLEMENTATION
#include "alea_test.h"
#include "alea.h"
#include "alea_slice.h"
#include "alea_mcnp.h"

/* =========================================================================
 * Helpers
 * ========================================================================= */

/* Build a system with a plane at x=split_x dividing a sphere of radius R.
 * Cell 1 (mat 1): x < split_x, inside sphere
 * Cell 2 (mat 2): x >= split_x, inside sphere
 * Cell 3 (void):  outside sphere
 */
static alea_system_t* make_x_split(double split_x, double R) {
    alea_system_t* sys = alea_create();
    if (!sys) return NULL;

    /* d = -split_x: plane equation x + d = 0 → plane at x = split_x */
    int p_idx = alea_plane_surface(sys, 1, 1.0, 0.0, 0.0, -split_x);
    int s_idx = alea_sphere_surface(sys, 2, 0.0, 0.0, 0.0, R);
    if (p_idx < 0 || s_idx < 0) { alea_destroy(sys); return NULL; }

    /* half-spaces: plane neg_node = x < split_x, pos = x >= split_x */
    alea_node_id_t plane_neg = alea_halfspace(sys, p_idx, -1); /* x < split_x */
    alea_node_id_t plane_pos = alea_halfspace(sys, p_idx,  1); /* x >= split_x */
    alea_node_id_t sph_in    = alea_halfspace(sys, s_idx, -1); /* inside sphere */

    alea_node_id_t cell1_node = alea_intersection(sys, sph_in,    plane_neg);
    alea_node_id_t cell2_node = alea_intersection(sys, sph_in,    plane_pos);
    alea_node_id_t void_node  = alea_halfspace(sys, s_idx, 1);   /* outside sphere */

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);

    alea_add_cell(sys, 1, cell1_node, m1, -1.0, 0);
    alea_add_cell(sys, 2, cell2_node, m2, -2.0, 0);
    alea_add_cell(sys, 3, void_node,  ALEA_MATERIAL_VOID, 0.0, 0);

    alea_build_universe_index(sys);
    return sys;
}

/* MCNP input for three X-aligned cells split by planes at x=4 and x=8.
 * Cell 1 (mat 1): x < 4, inside sphere SO 12
 * Cell 2 (mat 2): 4 <= x < 8
 * Cell 3 (mat 3): x >= 8
 * Cell 4 (void):  outside sphere
 */
static const char* THREE_CELLS_MCNP =
    "three cells\n"
    "1 1 -1.0 -3 -1\n"
    "2 2 -2.0 -3 1 -2\n"
    "3 3 -3.0 -3 2\n"
    "4 0 3\n"
    "\n"
    "1 PX 4\n"
    "2 PX 8\n"
    "3 SO 12\n"
    "\n"
    "M1 1001.80c 1.0\n"
    "M2 8016.80c 1.0\n"
    "M3 2004.80c 1.0\n";

static mcnp_model_t* load_three_cells(void) {
    FILE* f = fopen("_three_cells_tmp.mcnp", "w");
    if (!f) return NULL;
    fputs(THREE_CELLS_MCNP, f);
    fclose(f);
    mcnp_model_t* m = mcnp_load("_three_cells_tmp.mcnp");
    remove("_three_cells_tmp.mcnp");
    return m;
}

/* Count pixels with error flag set. */
static int count_errors(const uint8_t* errors, int n) {
    int c = 0;
    for (int i = 0; i < n; i++) if (errors[i]) c++;
    return c;
}

static int count_coverage(const uint8_t* coverage, int n, uint8_t value) {
    int c = 0;
    for (int i = 0; i < n; i++) if (coverage[i] == value) c++;
    return c;
}

static int count_components_of_kind(
    const alea_plot_error_component_result_t* result,
    alea_plot_error_kind_t kind) {
    int c = 0;
    if (!result) return 0;
    for (size_t i = 0; i < result->component_count; i++)
        if (result->components[i].kind == kind) c++;
    return c;
}

/* =========================================================================
 * Test 1: uniform region — single cell, all pixels land in it
 * ========================================================================= */
TEST(grid_uniform) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int s = alea_sphere_surface(sys, 1, 0, 0, 0, 20.0);
    ASSERT(s >= 0);
    alea_node_id_t inner = alea_halfspace(sys, s, -1);
    alea_node_id_t outer = alea_halfspace(sys, s,  1);
    int mat = alea_add_material(sys, 1);
    alea_add_cell(sys, 1, inner, mat,  -1.0, 0);
    alea_add_cell(sys, 2, outer, ALEA_MATERIAL_VOID, 0.0, 0);
    alea_build_universe_index(sys);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    const int nu = 8, nv = 4;
    int cell_ids[32]; int mat_ids[32];
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -10.0, 10.0, -5.0, 5.0);

    int rc = alea_find_cells_grid(sys, &view, nu, nv, -1, cell_ids, mat_ids, NULL);
    ASSERT_EQ(rc, 0);

    /* All 32 pixels should be inside the sphere → material 1 */
    for (int i = 0; i < nu * nv; i++) {
        ASSERT_EQ(mat_ids[i], 1);
    }

    alea_destroy(sys);
}

/* =========================================================================
 * Test 2: two-cell X split — boundary pixel must get the correct cell
 *
 * The plane sits at x=0; grid spans [-10, 10].  With 20 pixels, each is
 * 1 unit wide.  Pixel centres: -9.5, -8.5, …, -0.5 | 0.5, 1.5, …, 9.5
 * The first right-side pixel (x=0.5) gets a stale hint from x=-0.5 in
 * the first pass; the verify/re-scan step must correct it to mat=2.
 * ========================================================================= */
TEST(grid_two_cells_boundary) {
    alea_system_t* sys = make_x_split(0.0, 15.0);
    ASSERT_NOT_NULL(sys);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    const int nu = 20, nv = 4;
    int cell_ids[80]; int mat_ids[80];
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -10.0, 10.0, -2.0, 2.0);

    int rc = alea_find_cells_grid(sys, &view, nu, nv, -1, cell_ids, mat_ids, NULL);
    ASSERT_EQ(rc, 0);

    /* Row 0: pixels 0-9 (x<0) → mat 1; pixels 10-19 (x>=0) → mat 2 */
    for (int j = 0; j < nv; j++) {
        for (int i = 0; i < nu; i++) {
            int expected_mat = (i < nu / 2) ? 1 : 2;
            ASSERT_MSG(mat_ids[j * nu + i] == expected_mat,
                       "boundary pixel has wrong material");
        }
    }

    alea_destroy(sys);
}

/* =========================================================================
 * Test 3: three cells in a row — two transitions, both corrected
 * ========================================================================= */
TEST(grid_three_cells_transitions) {
    /* planes at x=4 and x=8 via MCNP inline model */
    mcnp_model_t* model = load_three_cells();
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    /* Sanity-check: x<4 → mat1, 4<=x<8 → mat2, x>=8 → mat3 */
    ASSERT_EQ(alea_material_at(sys, -4.0, 0.0, 0.0), 1);
    ASSERT_EQ(alea_material_at(sys,  6.0, 0.0, 0.0), 2);
    ASSERT_EQ(alea_material_at(sys, 10.0, 0.0, 0.0), 3);

    /* 24 pixels over [-12, 12]: 1 pixel per unit.
     * Centres: -11.5..3.5 → mat1 (pixels 0-15)
     *            4.5..7.5 → mat2 (pixels 16-19)
     *            8.5..11.5 → mat3 (pixels 20-23) */
    const int nu = 24, nv = 2;
    int cell_ids[48]; int mat_ids[48];
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -12.0, 12.0, -1.0, 1.0);

    ASSERT_EQ(alea_find_cells_grid(sys, &view, nu, nv, -1, cell_ids, mat_ids, NULL), 0);

    for (int j = 0; j < nv; j++) {
        for (int i = 0; i < nu; i++) {
            int idx = j * nu + i;
            double x_centre = -12.0 + (i + 0.5) * (24.0 / nu);
            int expected = (x_centre < 4.0) ? 1 : (x_centre < 8.0) ? 2 : 3;
            ASSERT_EQ(mat_ids[idx], expected);
        }
    }

    mcnp_model_destroy(model);
}

/* =========================================================================
 * Test 4: overlap detection — two partially overlapping spheres
 * ========================================================================= */
TEST(grid_overlap_detected) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Two spheres that overlap in x in [0, 5] */
    int s1 = alea_sphere_surface(sys, 1, -2.5, 0, 0, 7.5); /* centre -2.5, r=7.5 */
    int s2 = alea_sphere_surface(sys, 2,  2.5, 0, 0, 7.5); /* centre +2.5, r=7.5 */
    ASSERT(s1 >= 0 && s2 >= 0);

    alea_node_id_t in1 = alea_halfspace(sys, s1, -1);
    alea_node_id_t in2 = alea_halfspace(sys, s2, -1);

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);
    alea_add_cell(sys, 1, in1, m1, -1.0, 0);
    alea_add_cell(sys, 2, in2, m2, -2.0, 0);
    /* No void cell — this is an intentional overlap model */

    alea_build_universe_index(sys);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    const int nu = 20, nv = 4;
    int cell_ids[80]; int mat_ids[80];
    uint8_t errors[80];
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -10.0, 10.0, -2.0, 2.0);

    ASSERT_EQ(alea_find_cells_grid(sys, &view, nu, nv, -1, cell_ids, mat_ids, errors), 0);

    /* At least some overlap pixels must be flagged */
    ASSERT_MSG(count_errors(errors, nu * nv) > 0,
               "expected overlap pixels in centre zone");

    /* Pixels well to the left (x ~ -9) should be clean (only sphere 1) */
    for (int j = 0; j < nv; j++)
        ASSERT_EQ(errors[j * nu + 0], 0); /* leftmost pixel */

    /* Pixels well to the right (x ~ 9) should be clean (only sphere 2) */
    for (int j = 0; j < nv; j++)
        ASSERT_EQ(errors[j * nu + nu - 1], 0); /* rightmost pixel */

    alea_destroy(sys);
}

/* =========================================================================
 * Test 4b: coverage grid classifies undefined/one/multi pixels
 * ========================================================================= */
TEST(grid_coverage_classes) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int s1 = alea_sphere_surface(sys, 1, -2.5, 0, 0, 7.5);
    int s2 = alea_sphere_surface(sys, 2,  2.5, 0, 0, 7.5);
    ASSERT(s1 >= 0 && s2 >= 0);

    alea_node_id_t in1 = alea_halfspace(sys, s1, -1);
    alea_node_id_t in2 = alea_halfspace(sys, s2, -1);

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);
    alea_add_cell(sys, 1, in1, m1, -1.0, 0);
    alea_add_cell(sys, 2, in2, m2, -2.0, 0);

    alea_build_universe_index(sys);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    const int nu = 40, nv = 8;
    int cell_ids[320]; int mat_ids[320]; int secondary[320];
    uint8_t errors[320]; uint8_t coverage[320];
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -14.0, 14.0, -4.0, 4.0);

    ASSERT_EQ(alea_find_cells_grid_coverage(
                  sys, &view, nu, nv, -1, ALEA_GRID_COVERAGE_FAST,
                  cell_ids, mat_ids, secondary, coverage, errors), 0);

    ASSERT_MSG(count_coverage(coverage, nu * nv, ALEA_COVERAGE_NONE) > 0,
               "expected undefined pixels outside both spheres");
    ASSERT_MSG(count_coverage(coverage, nu * nv, ALEA_COVERAGE_ONE) > 0,
               "expected single-coverage pixels");
    ASSERT_MSG(count_coverage(coverage, nu * nv, ALEA_COVERAGE_MULTI) > 0,
               "expected overlap pixels");

    for (int i = 0; i < nu * nv; i++) {
        ASSERT_EQ(secondary[i], -1);
        if (coverage[i] == ALEA_COVERAGE_NONE) {
            ASSERT(cell_ids[i] < 0 || errors[i] == ALEA_GRID_UNDEFINED);
        } else if (coverage[i] == ALEA_COVERAGE_MULTI) {
            ASSERT_EQ(errors[i], ALEA_GRID_OVERLAP);
        }
    }

    alea_destroy(sys);
}

/* =========================================================================
 * Test 4c: exact coverage catches total/nested overlap with no winner boundary
 * ========================================================================= */
TEST(grid_exact_coverage_nested_overlap) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int outer = alea_sphere_surface(sys, 1, 0, 0, 0, 10.0);
    int inner = alea_sphere_surface(sys, 2, 0, 0, 0, 3.0);
    ASSERT(outer >= 0 && inner >= 0);

    alea_node_id_t outer_in = alea_halfspace(sys, outer, -1);
    alea_node_id_t inner_in = alea_halfspace(sys, inner, -1);

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);
    alea_add_cell(sys, 1, outer_in, m1, -1.0, 0);
    alea_add_cell(sys, 2, inner_in, m2, -2.0, 0);

    alea_build_universe_index(sys);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    const int nu = 80, nv = 80;
    int* cell_ids = calloc((size_t)nu * nv, sizeof(int));
    int* secondary = calloc((size_t)nu * nv, sizeof(int));
    uint8_t* errors = calloc((size_t)nu * nv, sizeof(uint8_t));
    uint8_t* coverage = calloc((size_t)nu * nv, sizeof(uint8_t));
    ASSERT_NOT_NULL(cell_ids);
    ASSERT_NOT_NULL(secondary);
    ASSERT_NOT_NULL(errors);
    ASSERT_NOT_NULL(coverage);

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -12.0, 12.0, -12.0, 12.0);

    ASSERT_EQ(alea_find_cells_grid_coverage(
                  sys, &view, nu, nv, -1, ALEA_GRID_COVERAGE_EXACT,
                  cell_ids, NULL, secondary, coverage, errors), 0);

    int multi = count_coverage(coverage, nu * nv, ALEA_COVERAGE_MULTI);
    ASSERT_MSG(multi > 0, "expected exact coverage to catch nested overlap");
    ASSERT_MSG(count_errors(errors, nu * nv) >= multi,
               "overlap errors should cover multi-coverage pixels");

    int center = (nv / 2) * nu + (nu / 2);
    ASSERT_EQ(coverage[center], ALEA_COVERAGE_MULTI);
    ASSERT(secondary[center] > 0);
    ASSERT_EQ(errors[center], ALEA_GRID_OVERLAP);

    free(coverage);
    free(errors);
    free(secondary);
    free(cell_ids);
    alea_destroy(sys);
}

/* =========================================================================
 * Test 4d: tile exact coverage catches nested overlap from fast grid baseline
 * ========================================================================= */
TEST(grid_tile_coverage_nested_overlap) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int outer = alea_sphere_surface(sys, 1, 0, 0, 0, 10.0);
    int inner = alea_sphere_surface(sys, 2, 0, 0, 0, 3.0);
    ASSERT(outer >= 0 && inner >= 0);

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);
    alea_add_cell(sys, 1, alea_halfspace(sys, outer, -1), m1, -1.0, 0);
    alea_add_cell(sys, 2, alea_halfspace(sys, inner, -1), m2, -2.0, 0);

    alea_build_universe_index(sys);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    const int nu = 64, nv = 64;
    int* cell_ids = calloc((size_t)nu * nv, sizeof(int));
    int* secondary = calloc((size_t)nu * nv, sizeof(int));
    uint8_t* errors = calloc((size_t)nu * nv, sizeof(uint8_t));
    uint8_t* coverage = calloc((size_t)nu * nv, sizeof(uint8_t));
    ASSERT_NOT_NULL(cell_ids);
    ASSERT_NOT_NULL(secondary);
    ASSERT_NOT_NULL(errors);
    ASSERT_NOT_NULL(coverage);

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -12.0, 12.0, -12.0, 12.0);

    ASSERT_EQ(alea_find_cells_grid_coverage(
                  sys, &view, nu, nv, -1, ALEA_GRID_COVERAGE_FAST,
                  cell_ids, NULL, secondary, coverage, errors), 0);
    int before = count_coverage(coverage, nu * nv, ALEA_COVERAGE_MULTI);

    int refined = alea_refine_grid_coverage_tiles_exact(
        sys, &view, nu, nv, -1, 16, 16, secondary, coverage, errors);
    ASSERT(refined > 0);
    alea_tile_coverage_stats_t stats = alea_tile_coverage_stats_get();
    ASSERT(stats.tiles > 0);
    ASSERT(stats.pixels > 0);
    ASSERT(stats.dedup_candidate_max > 0);
    ASSERT_EQ((int)stats.refined_pixels, refined);

    int after = count_coverage(coverage, nu * nv, ALEA_COVERAGE_MULTI);
    ASSERT(after > before);
    int center = (nv / 2) * nu + (nu / 2);
    ASSERT_EQ(coverage[center], ALEA_COVERAGE_MULTI);
    ASSERT(secondary[center] > 0);
    alea_plot_error_component_result_t* comps =
        alea_classify_plot_error_components(
            cell_ids, secondary, coverage, nu, nv);
    ASSERT_NOT_NULL(comps);
    ASSERT(count_components_of_kind(comps, ALEA_PLOT_ERR_TOTAL_OVERLAP) > 0);
    alea_plot_error_components_free(comps);

    free(coverage);
    free(errors);
    free(secondary);
    free(cell_ids);
    alea_destroy(sys);
}

/* =========================================================================
 * Test 4e: component classification distinguishes partial overlap
 * ========================================================================= */
TEST(grid_component_partial_overlap) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int s1 = alea_sphere_surface(sys, 1, -2, 0, 0, 5.0);
    int s2 = alea_sphere_surface(sys, 2,  2, 0, 0, 5.0);
    ASSERT(s1 >= 0 && s2 >= 0);

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);
    alea_add_cell(sys, 1, alea_halfspace(sys, s1, -1), m1, -1.0, 0);
    alea_add_cell(sys, 2, alea_halfspace(sys, s2, -1), m2, -2.0, 0);

    alea_build_universe_index(sys);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    const int nu = 80, nv = 48;
    int* cell_ids = calloc((size_t)nu * nv, sizeof(int));
    int* secondary = calloc((size_t)nu * nv, sizeof(int));
    uint8_t* errors = calloc((size_t)nu * nv, sizeof(uint8_t));
    uint8_t* coverage = calloc((size_t)nu * nv, sizeof(uint8_t));
    ASSERT_NOT_NULL(cell_ids);
    ASSERT_NOT_NULL(secondary);
    ASSERT_NOT_NULL(errors);
    ASSERT_NOT_NULL(coverage);

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -8.0, 8.0, -6.0, 6.0);

    ASSERT_EQ(alea_find_cells_grid_coverage(
                  sys, &view, nu, nv, -1,
                  ALEA_GRID_COVERAGE_EXACT | ALEA_GRID_SECONDARY_CELL_IDS,
                  cell_ids, NULL, secondary, coverage, errors), 0);
    ASSERT(count_coverage(coverage, nu * nv, ALEA_COVERAGE_MULTI) > 0);

    alea_plot_error_component_result_t* comps =
        alea_classify_plot_error_components(
            cell_ids, secondary, coverage, nu, nv);
    ASSERT_NOT_NULL(comps);
    ASSERT(count_components_of_kind(comps, ALEA_PLOT_ERR_PARTIAL_OVERLAP) > 0);
    alea_plot_error_components_free(comps);

    free(coverage);
    free(errors);
    free(secondary);
    free(cell_ids);
    alea_destroy(sys);
}

/* =========================================================================
 * Test 4f: component classification reports undefined regions
 * ========================================================================= */
TEST(grid_component_undefined_region) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int s1 = alea_sphere_surface(sys, 1, -4, 0, 0, 2.0);
    int s2 = alea_sphere_surface(sys, 2,  4, 0, 0, 2.0);
    ASSERT(s1 >= 0 && s2 >= 0);

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);
    alea_add_cell(sys, 1, alea_halfspace(sys, s1, -1), m1, -1.0, 0);
    alea_add_cell(sys, 2, alea_halfspace(sys, s2, -1), m2, -2.0, 0);

    alea_build_universe_index(sys);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    const int nu = 48, nv = 32;
    int* cell_ids = calloc((size_t)nu * nv, sizeof(int));
    uint8_t* errors = calloc((size_t)nu * nv, sizeof(uint8_t));
    uint8_t* coverage = calloc((size_t)nu * nv, sizeof(uint8_t));
    ASSERT_NOT_NULL(cell_ids);
    ASSERT_NOT_NULL(errors);
    ASSERT_NOT_NULL(coverage);

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -8.0, 8.0, -5.0, 5.0);

    ASSERT_EQ(alea_find_cells_grid_coverage(
                  sys, &view, nu, nv, -1, ALEA_GRID_COVERAGE_FAST,
                  cell_ids, NULL, NULL, coverage, errors), 0);
    ASSERT(count_coverage(coverage, nu * nv, ALEA_COVERAGE_NONE) > 0);

    alea_plot_error_component_result_t* comps =
        alea_classify_plot_error_components(
            cell_ids, NULL, coverage, nu, nv);
    ASSERT_NOT_NULL(comps);
    ASSERT(count_components_of_kind(comps, ALEA_PLOT_ERR_UNDEFINED_REGION) > 0);
    alea_plot_error_components_free(comps);

    free(coverage);
    free(errors);
    free(cell_ids);
    alea_destroy(sys);
}

/* =========================================================================
 * Test 4g: exact coverage honors explicit universe depth
 * ========================================================================= */
TEST(grid_exact_coverage_explicit_universe_depth) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int root_s = alea_sphere_surface(sys, 1, 0, 0, 0, 10.0);
    int c1_s = alea_sphere_surface(sys, 2, -1, 0, 0, 4.0);
    int c2_s = alea_sphere_surface(sys, 3,  1, 0, 0, 4.0);
    ASSERT(root_s >= 0 && c1_s >= 0 && c2_s >= 0);

    int root_mat = alea_add_material(sys, 10);
    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);

    int root_cell = alea_add_cell(
        sys, 10, alea_halfspace(sys, root_s, -1), root_mat, -1.0, 0);
    ASSERT(root_cell >= 0);
    ASSERT_EQ(alea_set_fill(sys, root_cell, 1, 0), 0);

    alea_add_cell(sys, 11, alea_halfspace(sys, c1_s, -1), m1, -1.0, 1);
    alea_add_cell(sys, 12, alea_halfspace(sys, c2_s, -1), m2, -2.0, 1);

    alea_build_universe_index(sys);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    const int nu = 48, nv = 48;
    int* cell_ids = calloc((size_t)nu * nv, sizeof(int));
    int* secondary = calloc((size_t)nu * nv, sizeof(int));
    uint8_t* errors = calloc((size_t)nu * nv, sizeof(uint8_t));
    uint8_t* coverage = calloc((size_t)nu * nv, sizeof(uint8_t));
    ASSERT_NOT_NULL(cell_ids);
    ASSERT_NOT_NULL(secondary);
    ASSERT_NOT_NULL(errors);
    ASSERT_NOT_NULL(coverage);

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -6.0, 6.0, -6.0, 6.0);

    ASSERT_EQ(alea_find_cells_grid_coverage(
                  sys, &view, nu, nv, 1,
                  ALEA_GRID_COVERAGE_EXACT | ALEA_GRID_SECONDARY_CELL_IDS,
                  cell_ids, NULL, secondary, coverage, errors), 0);

    int center = (nv / 2) * nu + (nu / 2);
    ASSERT_EQ(coverage[center], ALEA_COVERAGE_MULTI);
    ASSERT_EQ(errors[center], ALEA_GRID_OVERLAP);
    ASSERT(secondary[center] > 0);
    alea_point_coverage_stats_t pc_stats = alea_point_coverage_stats_get();
    ASSERT_EQ((int)pc_stats.queries, nu * nv);
    ASSERT(pc_stats.spatial_queries > 0);
    ASSERT(pc_stats.spatial_multi_early_exit > 0);

    ASSERT_EQ(alea_find_cells_grid_coverage(
                  sys, &view, nu, nv, 0,
                  ALEA_GRID_COVERAGE_EXACT | ALEA_GRID_SECONDARY_CELL_IDS,
                  cell_ids, NULL, secondary, coverage, errors), 0);
    ASSERT_EQ(coverage[center], ALEA_COVERAGE_ONE);
    pc_stats = alea_point_coverage_stats_get();
    ASSERT_EQ((int)pc_stats.queries, nu * nv);
    ASSERT(pc_stats.spatial_queries > 0);

    free(coverage);
    free(errors);
    free(secondary);
    free(cell_ids);
    alea_destroy(sys);
}

/* =========================================================================
 * Test 5: no false overlaps — two spheres that do NOT overlap
 * ========================================================================= */
TEST(grid_no_false_overlaps) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int s1 = alea_sphere_surface(sys, 1, -6, 0, 0, 4.0);
    int s2 = alea_sphere_surface(sys, 2,  6, 0, 0, 4.0);
    int so = alea_sphere_surface(sys, 3,  0, 0, 0, 20.0);
    ASSERT(s1 >= 0 && s2 >= 0 && so >= 0);

    alea_node_id_t in1  = alea_halfspace(sys, s1, -1);
    alea_node_id_t in2  = alea_halfspace(sys, s2, -1);
    alea_node_id_t out1 = alea_halfspace(sys, s1,  1);
    alea_node_id_t out2 = alea_halfspace(sys, s2,  1);
    alea_node_id_t ino  = alea_halfspace(sys, so, -1);
    alea_node_id_t outo = alea_halfspace(sys, so,  1);

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);
    alea_node_id_t gap = alea_intersection(sys, alea_intersection(sys, ino, out1), out2);

    alea_add_cell(sys, 1, in1,  m1, -1.0, 0);
    alea_add_cell(sys, 2, in2,  m2, -2.0, 0);
    alea_add_cell(sys, 3, gap,  ALEA_MATERIAL_VOID, 0.0, 0);
    alea_add_cell(sys, 4, outo, ALEA_MATERIAL_VOID, 0.0, 0);

    alea_build_universe_index(sys);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    const int nu = 20, nv = 4;
    int cell_ids[80]; int mat_ids[80];
    uint8_t errors[80];
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -10.0, 10.0, -2.0, 2.0);

    ASSERT_EQ(alea_find_cells_grid(sys, &view, nu, nv, -1, cell_ids, mat_ids, errors), 0);
    ASSERT_EQ(count_errors(errors, nu * nv), 0);

    alea_destroy(sys);
}

/* =========================================================================
 * Test 6: flat/hier parity — simple_fill.mcnp produces identical cell
 * arrays in both spatial modes (regression test for the optimisation)
 * ========================================================================= */
static const char* SIMPLE_FILL_PATH = "tests/data/simple_fill.mcnp";

TEST(grid_flat_hier_parity) {
    /* flat mode */
    mcnp_model_t* mf = mcnp_load(SIMPLE_FILL_PATH);
    if (!mf) {
        /* skip gracefully if fixture not found from this working dir */
        return;
    }
    alea_system_t* sf = mf->sys;
    alea_config_t cfg = alea_get_config(sf);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_FLAT;
    alea_set_config(sf, &cfg);
    ASSERT_EQ(alea_prepare_query_acceleration(sf), 0);

    /* hier mode */
    mcnp_model_t* mh = mcnp_load(SIMPLE_FILL_PATH);
    ASSERT_NOT_NULL(mh);
    alea_system_t* sh = mh->sys;
    cfg = alea_get_config(sh);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_HIER;
    alea_set_config(sh, &cfg);
    ASSERT_EQ(alea_prepare_query_acceleration(sh), 0);

    const int nu = 32, nv = 32;
    int cells_flat[1024]; int mats_flat[1024];
    int cells_hier[1024]; int mats_hier[1024];

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -60.0, 60.0, -60.0, 60.0);

    ASSERT_EQ(alea_find_cells_grid(sf, &view, nu, nv, -1, cells_flat, mats_flat, NULL), 0);
    ASSERT_EQ(alea_find_cells_grid(sh, &view, nu, nv, -1, cells_hier, mats_hier, NULL), 0);

    int mismatches = 0;
    for (int i = 0; i < nu * nv; i++) {
        if (mats_flat[i] != mats_hier[i]) mismatches++;
    }
    ASSERT_MSG(mismatches == 0, "flat and hier modes gave different material maps");

    mcnp_model_destroy(mf);
    mcnp_model_destroy(mh);
}

TEST_MAIN()
