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
#include "core/alea_system.h"
#include "core/alea_universe.h"

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

static alea_system_t* make_z_split_sphere(double R) {
    alea_system_t* sys = alea_create();
    if (!sys) return NULL;

    int sphere = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, R);
    int plane = alea_plane_surface(sys, 2, 0.0, 0.0, 1.0, 0.0);
    if (sphere < 0 || plane < 0) { alea_destroy(sys); return NULL; }

    alea_node_id_t sphere_in = alea_halfspace(sys, sphere, -1);
    alea_node_id_t sphere_out = alea_halfspace(sys, sphere, 1);
    alea_node_id_t lower_half = alea_halfspace(sys, plane, -1);
    alea_node_id_t upper_half = alea_halfspace(sys, plane, 1);
    alea_node_id_t lower = alea_intersection(sys, sphere_in, lower_half);
    alea_node_id_t upper = alea_intersection(sys, sphere_in, upper_half);
    if (sphere_in == ALEA_NODE_ID_INVALID ||
        sphere_out == ALEA_NODE_ID_INVALID ||
        lower_half == ALEA_NODE_ID_INVALID ||
        upper_half == ALEA_NODE_ID_INVALID ||
        lower == ALEA_NODE_ID_INVALID ||
        upper == ALEA_NODE_ID_INVALID) {
        alea_destroy(sys);
        return NULL;
    }

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);
    alea_add_cell(sys, 1, lower, m1, -1.0, 0);
    alea_add_cell(sys, 2, upper, m2, -2.0, 0);
    alea_add_cell(sys, 3, sphere_out, ALEA_MATERIAL_VOID, 0.0, 0);
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

TEST(grid_coherence_reports_swept_owner_and_flags_overlap) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int first_sphere = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 2.0);
    int later_sphere = alea_sphere_surface(sys, 2, -3.0, 0.0, 0.0, 2.0);
    ASSERT(first_sphere >= 0 && later_sphere >= 0);
    int first_mat = alea_add_material(sys, 1);
    int later_mat = alea_add_material(sys, 2);
    ASSERT(first_mat >= 0 && later_mat >= 0);
    /* Cell 1 is first in deck order, but a left-to-right row starts in cell 2
     * alone and then enters their overlap, so the sweep carries cell 2 in.
     * The grid resolves ownership coherently rather than re-deriving deck order
     * per pixel; what it owes the caller in an overlap is the error flag, not a
     * particular one of the two claimants. */
    ASSERT(alea_add_cell(sys, 1, alea_halfspace(sys, first_sphere, -1),
                         first_mat, -1.0, 0) >= 0);
    ASSERT(alea_add_cell(sys, 2, alea_halfspace(sys, later_sphere, -1),
                         later_mat, -1.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    int cell_ids[6];
    int material_ids[6];
    uint8_t errors[6];
    alea_slice_view_t view;
    /* Centres: -4.5, -3.5, -2.5, -1.5, -0.5, 0.5. */
    alea_slice_view_axis(&view, 2, 0.0, -5.0, 1.0, -1.0, 1.0);
    ASSERT_EQ(alea_find_cells_grid(sys, &view, 6, 1, -1,
                                   cell_ids, material_ids, errors), 0);
    /* x = -4.5: cell 2 alone, before the overlap. */
    ASSERT_EQ(material_ids[0], 2);
    ASSERT_EQ(errors[0], ALEA_GRID_OK);
    /* x = -1.5: claimed by both. The swept cell is kept... */
    ASSERT_EQ(cell_ids[3], 2);
    ASSERT_EQ(material_ids[3], 2);
    /* ...and the pixel is reported as an overlap regardless of which won. */
    ASSERT_EQ(errors[3], ALEA_GRID_OVERLAP);
    /* x = -0.5: cell 1 alone, past the overlap. */
    ASSERT_EQ(cell_ids[4], 1);
    ASSERT_EQ(errors[4], ALEA_GRID_OK);

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

    alea_tile_refinement_options_t options;
    alea_tile_refinement_options_init(&options);
    ASSERT_EQ(options.max_candidates, 4096);
    ASSERT_EQ(options.max_evaluated_candidates, 4096);
    ASSERT_EQ(options.cap_storm_threshold, 16);
    /* Exercise the explicit contract rather than hidden process environment. */
    options.use_hier_chain_candidates = true;
    options.use_chain_bitset = true;
    int refined = alea_refine_grid_coverage_tiles_exact_ex(
        sys, &view, nu, nv, -1, 16, 16, &options,
        secondary, coverage, errors);
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

    /* A bounded fallback must retain its uncertainty instead of completing an
     * unplanned exact raster pass.  Two sphere candidates exceed this tiny
     * chain-evaluation cap, while one 16x16 tile exceeds the fallback budget. */
    memset(secondary, 0, (size_t)nu * nv * sizeof(*secondary));
    memset(errors, 0, (size_t)nu * nv * sizeof(*errors));
    memset(coverage, 0, (size_t)nu * nv * sizeof(*coverage));
    ASSERT_EQ(alea_find_cells_grid_coverage(
                  sys, &view, nu, nv, -1, ALEA_GRID_COVERAGE_FAST,
                  cell_ids, NULL, secondary, coverage, errors), 0);
    alea_tile_refinement_options_init(&options);
    options.use_hier_chain_candidates = true;
    options.max_evaluated_candidates = 1;
    options.max_exact_fallback_pixels = 10;
    ASSERT_EQ(alea_refine_grid_coverage_tiles_exact_ex(
                  sys, &view, nu, nv, -1, 16, 16, &options,
                  secondary, coverage, errors), 0);
    stats = alea_tile_coverage_stats_get();
    ASSERT(stats.incomplete);
    ASSERT(stats.skipped_tiles > 0);
    ASSERT(stats.skipped_pixels > 0);
    ASSERT_EQ(stats.exact_fallback_pixels, 0);

    free(coverage);
    free(errors);
    free(secondary);
    free(cell_ids);
    alea_destroy(sys);
}

TEST(grid_tile_coverage_honors_cooperative_interrupt) {
    alea_system_t* sys = make_x_split(0.0, 10.0);
    ASSERT_NOT_NULL(sys);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    int secondary[1] = {-1};
    uint8_t coverage[1] = {ALEA_COVERAGE_ONE};
    uint8_t errors[1] = {0};
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -1.0, 1.0, -1.0, 1.0);
    alea_tile_refinement_options_t options;
    alea_tile_refinement_options_init(&options);

    alea_interrupt();
    ASSERT_EQ(alea_refine_grid_coverage_tiles_exact_ex(
                  sys, &view, 1, 1, -1, 1, 1, &options,
                  secondary, coverage, errors), -1);
    alea_clear_interrupt();
    alea_destroy(sys);
}

TEST(grid_path_coverage_nested_overlap) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    alea_config_t cfg = alea_get_config(sys);
    alea_set_config(sys, &cfg);

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
    uint32_t* path_ids = calloc((size_t)nu * nv, sizeof(uint32_t));
    alea_slice_path_table_t paths = {0};
    ASSERT_NOT_NULL(cell_ids);
    ASSERT_NOT_NULL(secondary);
    ASSERT_NOT_NULL(errors);
    ASSERT_NOT_NULL(coverage);
    ASSERT_NOT_NULL(path_ids);

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -12.0, 12.0, -12.0, 12.0);

    ASSERT_EQ(alea_find_cells_grid_coverage_paths(
                  sys, &view, nu, nv, -1,
                  ALEA_GRID_COVERAGE_FAST | ALEA_GRID_PATH_IDS,
                  cell_ids, NULL, secondary, coverage, errors,
                  path_ids, &paths), 0);
    ASSERT(paths.count > 0);
    int before = count_coverage(coverage, nu * nv, ALEA_COVERAGE_MULTI);

    int refined = alea_refine_grid_coverage_paths_exact(
        sys, &view, nu, nv, -1, 16, 16, cell_ids, path_ids, &paths,
        secondary, coverage, errors);
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

    /* The path engine must honor the same selected-tile mask as the global
     * tile engine; interactive diagnostics must not scan the full raster. */
    uint8_t selected_tiles[16] = {0};
    selected_tiles[0] = 1;
    alea_tile_refinement_options_t selected_options;
    alea_tile_refinement_options_init(&selected_options);
    selected_options.tile_mask = selected_tiles;
    selected_options.tile_mask_count = 16;
    ASSERT(alea_refine_grid_coverage_paths_exact_ex(
               sys, &view, nu, nv, -1, 16, 16, cell_ids, path_ids, &paths,
               &selected_options, secondary, coverage, errors) >= 0);
    stats = alea_tile_coverage_stats_get();
    ASSERT_EQ(stats.tiles, 1);

    /* The interactive producer must resolve paths only inside the requested
     * tiles.  Outside entries remain explicit missing paths and are never
     * mistaken for clean ownership evidence. */
    uint32_t* selected_path_ids =
        calloc((size_t)nu * nv, sizeof(*selected_path_ids));
    alea_slice_path_table_t selected_paths = {0};
    ASSERT_NOT_NULL(selected_path_ids);
    ASSERT_EQ(alea_find_cells_grid_paths_selected(
                  sys, &view, nu, nv, -1, 16, 16,
                  selected_tiles, 16, cell_ids,
                  selected_path_ids, &selected_paths), 0);
    ASSERT(selected_paths.count > 0);
    int selected_valid_paths = 0;
    for (int j = 0; j < nv; j++) {
        for (int i = 0; i < nu; i++) {
            size_t pidx = (size_t)j * (size_t)nu + (size_t)i;
            if (i < 16 && j < 16) {
                if (selected_path_ids[pidx] != UINT32_MAX)
                    selected_valid_paths++;
            } else {
                ASSERT_EQ(selected_path_ids[pidx], UINT32_MAX);
            }
        }
    }
    ASSERT(selected_valid_paths > 0);

    memset(secondary, 0, (size_t)nu * nv * sizeof(*secondary));
    memset(errors, 0, (size_t)nu * nv * sizeof(*errors));
    memset(coverage, 0, (size_t)nu * nv * sizeof(*coverage));
    ASSERT(alea_refine_grid_coverage_paths_exact_ex(
               sys, &view, nu, nv, -1, 16, 16,
               cell_ids, selected_path_ids, &selected_paths,
               &selected_options, secondary, coverage, errors) >= 0);
    stats = alea_tile_coverage_stats_get();
    ASSERT_EQ(stats.tiles, 1);
    ASSERT(stats.fallback_tiles <= stats.tiles);
    ASSERT(stats.exact_fallback_pixels > 0);
    alea_slice_path_table_free(&selected_paths);
    free(selected_path_ids);

    /* Missing concrete paths normally use recursive exact fallback.  Under a
     * finite budget, retain the second missing-path pixel as incomplete rather
     * than silently growing the diagnostic pass. */
    memset(secondary, 0, (size_t)nu * nv * sizeof(*secondary));
    memset(errors, 0, (size_t)nu * nv * sizeof(*errors));
    memset(coverage, 0, (size_t)nu * nv * sizeof(*coverage));
    path_ids[0] = UINT32_MAX;
    path_ids[1] = UINT32_MAX;
    alea_tile_refinement_options_t options;
    alea_tile_refinement_options_init(&options);
    options.max_exact_fallback_pixels = 1;
    ASSERT(alea_refine_grid_coverage_paths_exact_ex(
               sys, &view, nu, nv, -1, 16, 16, cell_ids, path_ids, &paths,
               &options, secondary, coverage, errors) >= 0);
    stats = alea_tile_coverage_stats_get();
    ASSERT(stats.incomplete);
    ASSERT(stats.skipped_pixels > 0);
    ASSERT_EQ(stats.exact_fallback_pixels, 1);

    alea_slice_path_table_free(&paths);
    free(path_ids);
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

TEST(grid_path_ids_filled_universe) {
    mcnp_model_t* model = mcnp_load(SIMPLE_FILL_PATH);
    if (!model) return;
    alea_system_t* sys = model->sys;
    alea_config_t cfg = alea_get_config(sys);
    alea_set_config(sys, &cfg);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    const int nu = 32, nv = 32;
    int cell_ids[1024];
    uint8_t coverage[1024];
    uint8_t errors[1024];
    uint32_t path_ids[1024];
    alea_slice_path_table_t paths = {0};

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -60.0, 60.0, -60.0, 60.0);

    ASSERT_EQ(alea_find_cells_grid_coverage_paths(
        sys, &view, nu, nv, -1,
        ALEA_GRID_COVERAGE_FAST | ALEA_GRID_PATH_IDS,
        cell_ids, NULL, NULL, coverage, errors, path_ids, &paths), 0);

    ASSERT(paths.count > 0);
    uint32_t fill_path = UINT32_MAX;
    int fill_pixels = 0;
    for (int i = 0; i < nu * nv; i++) {
        if (cell_ids[i] == 2 || cell_ids[i] == 3) {
            ASSERT(path_ids[i] != UINT32_MAX);
            if (fill_path == UINT32_MAX) fill_path = path_ids[i];
            ASSERT_EQ(path_ids[i], fill_path);
            fill_pixels++;
        }
    }
    ASSERT(fill_pixels > 0);
    ASSERT(fill_path < paths.count);
    ASSERT_EQ(paths.records[fill_path].universe_id, 1);
    ASSERT_EQ(paths.records[fill_path].depth, 1);

    alea_slice_path_table_free(&paths);
    mcnp_model_destroy(model);
}

TEST(grid_path_ids_distinguish_lattice_element_placements) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) return;
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    int cell_ids[3];
    uint8_t coverage[3];
    uint8_t errors[3];
    uint32_t path_ids[3];
    alea_slice_path_table_t paths = {0};
    alea_slice_view_t view;
    /* Pixel centers are x = 0, 2, 4. The first and third select the same
     * child universe but different lattice elements. */
    alea_slice_view_axis(&view, 2, 0.0, -1.0, 5.0, -1.0, 1.0);

    ASSERT_EQ(alea_find_cells_grid_coverage_paths(
                  sys, &view, 3, 1, -1,
                  ALEA_GRID_COVERAGE_FAST | ALEA_GRID_PATH_IDS,
                  cell_ids, NULL, NULL, coverage, errors, path_ids, &paths), 0);
    ASSERT_EQ(cell_ids[0], 1);
    ASSERT_EQ(cell_ids[2], 1);
    ASSERT(path_ids[0] != UINT32_MAX);
    ASSERT(path_ids[2] != UINT32_MAX);
    ASSERT(path_ids[0] != path_ids[2]);

    const alea_slice_path_record_t* left = &paths.records[path_ids[0]];
    const alea_slice_path_record_t* right = &paths.records[path_ids[2]];
    ASSERT_EQ(left->universe_id, 1);
    ASSERT_EQ(right->universe_id, 1);
    ASSERT_EQ(left->depth, 1);
    ASSERT_EQ(right->depth, 1);
    ASSERT(left->world_to_local[3] != right->world_to_local[3]);

    alea_slice_path_table_free(&paths);
    mcnp_model_destroy(model);
}

TEST(grid_path_coverage_matches_recursive_lattice) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) return;
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    const int nu = 18, nv = 18;
    const size_t count = (size_t)nu * (size_t)nv;
    int* exact_cells = calloc(count, sizeof(*exact_cells));
    int* exact_secondary = calloc(count, sizeof(*exact_secondary));
    uint8_t* exact_coverage = calloc(count, sizeof(*exact_coverage));
    uint8_t* exact_errors = calloc(count, sizeof(*exact_errors));
    int* tile_cells = calloc(count, sizeof(*tile_cells));
    int* tile_secondary = calloc(count, sizeof(*tile_secondary));
    uint8_t* tile_coverage = calloc(count, sizeof(*tile_coverage));
    uint8_t* tile_errors = calloc(count, sizeof(*tile_errors));
    uint32_t* path_ids = calloc(count, sizeof(*path_ids));
    alea_slice_path_table_t paths = {0};
    ASSERT_NOT_NULL(exact_cells);
    ASSERT_NOT_NULL(exact_secondary);
    ASSERT_NOT_NULL(exact_coverage);
    ASSERT_NOT_NULL(exact_errors);
    ASSERT_NOT_NULL(tile_cells);
    ASSERT_NOT_NULL(tile_secondary);
    ASSERT_NOT_NULL(tile_coverage);
    ASSERT_NOT_NULL(tile_errors);
    ASSERT_NOT_NULL(path_ids);

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -1.0, 5.0, -1.0, 5.0);
    ASSERT_EQ(alea_find_cells_grid_coverage(
                  sys, &view, nu, nv, -1, ALEA_GRID_COVERAGE_EXACT,
                  exact_cells, NULL, exact_secondary,
                  exact_coverage, exact_errors), 0);
    ASSERT_EQ(alea_find_cells_grid_coverage_paths(
                  sys, &view, nu, nv, -1,
                  ALEA_GRID_COVERAGE_FAST | ALEA_GRID_PATH_IDS,
                  tile_cells, NULL, tile_secondary,
                  tile_coverage, tile_errors, path_ids, &paths), 0);

    alea_tile_refinement_options_t options;
    alea_tile_refinement_options_init(&options);
    ASSERT(alea_refine_grid_coverage_paths_exact_ex(
               sys, &view, nu, nv, -1, 6, 6, tile_cells, path_ids, &paths,
               &options,
               tile_secondary, tile_coverage, tile_errors) >= 0);

    for (size_t i = 0; i < count; i++) {
        if (tile_coverage[i] != exact_coverage[i] ||
            tile_errors[i] != exact_errors[i] ||
            tile_secondary[i] != exact_secondary[i]) {
            fprintf(stderr,
                    "lattice tile mismatch i=%zu px=%zu py=%zu exact=(%u,%u,%d) tile=(%u,%u,%d)\n",
                    i, i % (size_t)nu, i / (size_t)nu,
                    (unsigned)exact_coverage[i], (unsigned)exact_errors[i],
                    exact_secondary[i], (unsigned)tile_coverage[i],
                    (unsigned)tile_errors[i], tile_secondary[i]);
        }
        ASSERT_EQ(tile_coverage[i], exact_coverage[i]);
        ASSERT_EQ(tile_errors[i], exact_errors[i]);
        ASSERT_EQ(tile_secondary[i], exact_secondary[i]);
    }

    alea_slice_path_table_free(&paths);
    free(path_ids);
    free(tile_errors);
    free(tile_coverage);
    free(tile_secondary);
    free(tile_cells);
    free(exact_errors);
    free(exact_coverage);
    free(exact_secondary);
    free(exact_cells);
    mcnp_model_destroy(model);
}

TEST(grid_local_coverage_matches_exact_lattice_oracle) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) return;
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    const int nu = 18, nv = 18;
    const size_t count = (size_t)nu * (size_t)nv;
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -1.0, 5.0, -1.0, 5.0);
    int* primary = calloc(count, sizeof(*primary));
    int* secondary = calloc(count, sizeof(*secondary));
    uint8_t* coverage = calloc(count, sizeof(*coverage));
    uint8_t* errors = calloc(count, sizeof(*errors));
    ASSERT_NOT_NULL(primary); ASSERT_NOT_NULL(secondary);
    ASSERT_NOT_NULL(coverage); ASSERT_NOT_NULL(errors);
    ASSERT_EQ(alea_find_cells_grid_coverage(
        sys, &view, nu, nv, -1,
        ALEA_GRID_COVERAGE_EXACT | ALEA_GRID_SECONDARY_CELL_IDS,
        primary, NULL, secondary, coverage, errors), 0);
    ASSERT_EQ(alea_filter_grid_boundary_ambiguities(
        sys, &view, nu, nv, -1, primary, secondary, coverage, errors, NULL), 0);
    alea_plot_error_component_result_t* oracle =
        alea_classify_plot_error_components(primary, secondary, coverage, nu, nv);
    ASSERT_NOT_NULL(oracle);

    alea_plot_error_component_result_t* compact = NULL;
    alea_local_coverage_stats_t stats;
    ASSERT_EQ(alea_find_local_coverage_components(
        sys, &view, nu, nv, -1, count, 16 * 1024 * 1024, 1,
        &compact, &stats), 0);
    ASSERT_NOT_NULL(compact);
    ASSERT_EQ((int)compact->component_count, (int)oracle->component_count);
    ASSERT_EQ(count_components_of_kind(compact, ALEA_PLOT_ERR_UNDEFINED_REGION),
              count_components_of_kind(oracle, ALEA_PLOT_ERR_UNDEFINED_REGION));
    ASSERT_EQ(stats.point_coverage.lattice_fallbacks, 0);
    ASSERT_EQ(stats.point_coverage.recursive_fallbacks, 0);
    ASSERT_EQ(stats.point_coverage.spatial_queries, count);

    alea_plot_error_components_free(compact);
    alea_plot_error_components_free(oracle);
    free(errors); free(coverage); free(secondary); free(primary);
    mcnp_model_destroy(model);
}

/* The compact local query must retain multiple terminal claimants that occur
 * within one concrete lattice element.  This is specifically different from
 * merely seeing the same universe definition elsewhere in a lattice. */
TEST(grid_local_coverage_detects_lattice_terminal_overlap) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_overlap.mcnp");
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    /* The recursive reference path establishes the intended ownership at the
     * centre of the concrete lattice element: both child terminals claim it. */
    alea_cell_hit_t recursive_hits[16];
    int recursive_count = alea_find_all_cells_at_point_recursive(
        sys, 0.0, 0.0, 0.0, recursive_hits, 16);
    ASSERT(recursive_count >= 2);
    bool has_cell_1 = false, has_cell_2 = false;
    for (int i = 0; i < recursive_count; i++) {
        has_cell_1 |= recursive_hits[i].cell_id == 1;
        has_cell_2 |= recursive_hits[i].cell_id == 2;
    }
    ASSERT(has_cell_1);
    ASSERT(has_cell_2);

    const int nu = 32, nv = 32;
    const size_t count = (size_t)nu * (size_t)nv;
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -1.0, 1.0, -1.0, 1.0);

    int* primary = calloc(count, sizeof(*primary));
    int* secondary = calloc(count, sizeof(*secondary));
    uint8_t* coverage = calloc(count, sizeof(*coverage));
    uint8_t* errors = calloc(count, sizeof(*errors));
    ASSERT_NOT_NULL(primary); ASSERT_NOT_NULL(secondary);
    ASSERT_NOT_NULL(coverage); ASSERT_NOT_NULL(errors);
    ASSERT_EQ(alea_find_cells_grid_coverage(
        sys, &view, nu, nv, -1,
        ALEA_GRID_COVERAGE_EXACT | ALEA_GRID_SECONDARY_CELL_IDS,
        primary, NULL, secondary, coverage, errors), 0);
    ASSERT_EQ(alea_filter_grid_boundary_ambiguities(
        sys, &view, nu, nv, -1, primary, secondary, coverage, errors, NULL), 0);
    alea_plot_error_component_result_t* global =
        alea_classify_plot_error_components(primary, secondary, coverage, nu, nv);
    ASSERT_NOT_NULL(global);
    ASSERT(count_components_of_kind(global, ALEA_PLOT_ERR_TOTAL_OVERLAP) > 0);

    alea_plot_error_component_result_t* compact = NULL;
    alea_local_coverage_stats_t stats;
    ASSERT_EQ(alea_find_local_coverage_components(
        sys, &view, nu, nv, -1, count, 16 * 1024 * 1024, 1,
        &compact, &stats), 0);
    ASSERT_NOT_NULL(compact);
    ASSERT_EQ((int)compact->component_count, (int)global->component_count);
    ASSERT_EQ(count_components_of_kind(compact, ALEA_PLOT_ERR_TOTAL_OVERLAP),
              count_components_of_kind(global, ALEA_PLOT_ERR_TOTAL_OVERLAP));
    ASSERT_EQ(stats.point_coverage.lattice_fallbacks, 0);
    ASSERT_EQ(stats.point_coverage.recursive_fallbacks, 0);
    /* Boundary classification probes selected pixels around the component, so
     * this can exceed the raster size; all probes must remain spatial. */
    ASSERT(stats.point_coverage.spatial_queries >= count);

    alea_plot_error_components_free(compact);
    alea_plot_error_components_free(global);
    free(errors); free(coverage); free(secondary); free(primary);
    mcnp_model_destroy(model);
}

/* Repeated transforms of the same lattice universe must remain distinct
 * occurrences.  Each placement contains the same two conflicting child
 * definitions, but the local diagnostic must retain two separate overlaps. */
TEST(grid_local_coverage_distinguishes_transformed_lattice_overlaps) {
    mcnp_model_t* model = mcnp_load(
        "tests/data/mcnp_transformed_lattice_overlap.mcnp");
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    alea_cell_hit_t left_hits[16], right_hits[16];
    uint64_t left_keys[16], right_keys[16];
    uint64_t left_parents[16], right_parents[16];
    int left_count = alea_find_all_cells_at_point_coverage_chain_recursive(
        sys, -5.0, 0.0, 0.0, left_hits, left_keys, left_parents, 16);
    int right_count = alea_find_all_cells_at_point_coverage_chain_recursive(
        sys, 5.0, 0.0, 0.0, right_hits, right_keys, right_parents, 16);
    ASSERT(left_count >= 2);
    ASSERT(right_count >= 2);
    uint64_t left_cell_1_key = 0, right_cell_1_key = 0;
    bool left_has_cell_2 = false, right_has_cell_2 = false;
    for (int i = 0; i < left_count; i++) {
        if (left_hits[i].cell_id == 1) left_cell_1_key = left_keys[i];
        if (left_hits[i].cell_id == 2) left_has_cell_2 = true;
    }
    for (int i = 0; i < right_count; i++) {
        if (right_hits[i].cell_id == 1) right_cell_1_key = right_keys[i];
        if (right_hits[i].cell_id == 2) right_has_cell_2 = true;
    }
    ASSERT(left_cell_1_key != 0);
    ASSERT(right_cell_1_key != 0);
    ASSERT(left_cell_1_key != right_cell_1_key);
    ASSERT(left_has_cell_2);
    ASSERT(right_has_cell_2);

    const int nu = 112, nv = 32;
    const size_t count = (size_t)nu * (size_t)nv;
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -7.0, 7.0, -1.0, 1.0);
    int* primary = calloc(count, sizeof(*primary));
    int* secondary = calloc(count, sizeof(*secondary));
    uint8_t* coverage = calloc(count, sizeof(*coverage));
    uint8_t* errors = calloc(count, sizeof(*errors));
    ASSERT_NOT_NULL(primary); ASSERT_NOT_NULL(secondary);
    ASSERT_NOT_NULL(coverage); ASSERT_NOT_NULL(errors);
    ASSERT_EQ(alea_find_cells_grid_coverage(
        sys, &view, nu, nv, -1,
        ALEA_GRID_COVERAGE_EXACT | ALEA_GRID_SECONDARY_CELL_IDS,
        primary, NULL, secondary, coverage, errors), 0);
    ASSERT_EQ(alea_filter_grid_boundary_ambiguities(
        sys, &view, nu, nv, -1, primary, secondary, coverage, errors, NULL), 0);
    alea_plot_error_component_result_t* global =
        alea_classify_plot_error_components(primary, secondary, coverage, nu, nv);
    ASSERT_NOT_NULL(global);
    ASSERT_EQ(count_components_of_kind(global, ALEA_PLOT_ERR_TOTAL_OVERLAP), 2);

    alea_plot_error_component_result_t* compact = NULL;
    alea_local_coverage_stats_t stats;
    ASSERT_EQ(alea_find_local_coverage_components(
        sys, &view, nu, nv, -1, count, 16 * 1024 * 1024, 1,
        &compact, &stats), 0);
    ASSERT_NOT_NULL(compact);
    ASSERT_EQ((int)compact->component_count, (int)global->component_count);
    ASSERT_EQ(count_components_of_kind(compact, ALEA_PLOT_ERR_TOTAL_OVERLAP), 2);
    ASSERT_EQ(stats.point_coverage.lattice_fallbacks, 0);
    ASSERT_EQ(stats.point_coverage.recursive_fallbacks, 0);

    alea_plot_error_components_free(compact);
    alea_plot_error_components_free(global);
    free(errors); free(coverage); free(secondary); free(primary);
    mcnp_model_destroy(model);
}

TEST(grid_local_coverage_detects_nested_lattice_terminal_overlap) {
    mcnp_model_t* model = mcnp_load(
        "tests/data/mcnp_nested_lattice_overlap.mcnp");
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    alea_cell_hit_t recursive_hits[16];
    int recursive_count = alea_find_all_cells_at_point_recursive(
        sys, 0.0, 0.0, 0.0, recursive_hits, 16);
    ASSERT(recursive_count >= 2);
    bool has_cell_1 = false, has_cell_2 = false;
    for (int i = 0; i < recursive_count; i++) {
        has_cell_1 |= recursive_hits[i].cell_id == 1;
        has_cell_2 |= recursive_hits[i].cell_id == 2;
    }
    ASSERT(has_cell_1);
    ASSERT(has_cell_2);

    const int nu = 32, nv = 32;
    const size_t count = (size_t)nu * (size_t)nv;
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -1.0, 1.0, -1.0, 1.0);
    int* primary = calloc(count, sizeof(*primary));
    int* secondary = calloc(count, sizeof(*secondary));
    uint8_t* coverage = calloc(count, sizeof(*coverage));
    uint8_t* errors = calloc(count, sizeof(*errors));
    ASSERT_NOT_NULL(primary); ASSERT_NOT_NULL(secondary);
    ASSERT_NOT_NULL(coverage); ASSERT_NOT_NULL(errors);
    ASSERT_EQ(alea_find_cells_grid_coverage(
        sys, &view, nu, nv, -1,
        ALEA_GRID_COVERAGE_EXACT | ALEA_GRID_SECONDARY_CELL_IDS,
        primary, NULL, secondary, coverage, errors), 0);
    ASSERT_EQ(alea_filter_grid_boundary_ambiguities(
        sys, &view, nu, nv, -1, primary, secondary, coverage, errors, NULL), 0);
    alea_plot_error_component_result_t* global =
        alea_classify_plot_error_components(primary, secondary, coverage, nu, nv);
    ASSERT_NOT_NULL(global);
    ASSERT(count_components_of_kind(global, ALEA_PLOT_ERR_TOTAL_OVERLAP) > 0);

    alea_plot_error_component_result_t* compact = NULL;
    alea_local_coverage_stats_t stats;
    ASSERT_EQ(alea_find_local_coverage_components(
        sys, &view, nu, nv, -1, count, 16 * 1024 * 1024, 1,
        &compact, &stats), 0);
    ASSERT_NOT_NULL(compact);
    ASSERT_EQ((int)compact->component_count, (int)global->component_count);
    ASSERT_EQ(count_components_of_kind(compact, ALEA_PLOT_ERR_TOTAL_OVERLAP),
              count_components_of_kind(global, ALEA_PLOT_ERR_TOTAL_OVERLAP));
    ASSERT_EQ(stats.point_coverage.lattice_fallbacks, 0);
    ASSERT_EQ(stats.point_coverage.recursive_fallbacks, 0);

    alea_plot_error_components_free(compact);
    alea_plot_error_components_free(global);
    free(errors); free(coverage); free(secondary); free(primary);
    mcnp_model_destroy(model);
}

/* A terminal root claimant and a terminal claimant below a competing fill are
 * simultaneous concrete owners even though their leaf depths differ. The
 * compact local scanner must use occurrence-tree leaves, not one chosen depth. */
TEST(grid_local_coverage_detects_unequal_depth_terminal_overlap) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int sphere = alea_sphere_surface(sys, 1401, 0.0, 0.0, 0.0, 2.0);
    ASSERT(sphere >= 0);
    alea_node_id_t inside = alea_halfspace(sys, sphere, -1);
    int root_terminal = alea_add_cell(
        sys, 141, inside, ALEA_MATERIAL_VOID, 0.0, 0);
    int root_fill = alea_add_cell(
        sys, 142, inside, ALEA_MATERIAL_VOID, 0.0, 0);
    int child_terminal = alea_add_cell(
        sys, 143, inside, ALEA_MATERIAL_VOID, 0.0, 8);
    ASSERT(root_terminal >= 0);
    ASSERT(root_fill >= 0);
    ASSERT(child_terminal >= 0);
    ASSERT_EQ(alea_set_cell_fill(sys, root_fill, 8, 0), 0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -1.0, 1.0, -1.0, 1.0);
    alea_plot_error_component_result_t* components = NULL;
    alea_local_coverage_stats_t stats;
    ASSERT_EQ(alea_find_local_coverage_components(
        sys, &view, 8, 8, -1, 64, 1024 * 1024, 1,
        &components, &stats), 0);
    ASSERT_NOT_NULL(components);
    ASSERT_EQ(stats.incomplete_points, (size_t)0);
    ASSERT_EQ(components->component_count, (size_t)1);
    ASSERT_EQ(components->components[0].kind, ALEA_PLOT_ERR_TOTAL_OVERLAP);
    ASSERT_EQ(components->components[0].primary_cell_id, 141);
    ASSERT_EQ(components->components[0].secondary_cell_id, 143);

    alea_plot_error_components_free(components);
    alea_destroy(sys);
}

TEST(grid_local_coverage_reports_saturated_owner_queries_incomplete) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int sphere = alea_sphere_surface(sys, 1501, 0.0, 0.0, 0.0, 2.0);
    ASSERT(sphere >= 0);
    alea_node_id_t inside = alea_halfspace(sys, sphere, -1);
    for (int cell = 0; cell < 33; cell++)
        ASSERT(alea_add_cell(
            sys, 150 + cell, inside, ALEA_MATERIAL_VOID, 0.0, 0) >= 0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -1.0, 1.0, -1.0, 1.0);
    alea_plot_error_component_result_t* components = NULL;
    alea_local_coverage_stats_t stats;
    ASSERT_EQ(alea_find_local_coverage_components(
        sys, &view, 4, 4, -1, 16, 1024 * 1024, 1,
        &components, &stats), 0);
    ASSERT_NOT_NULL(components);
    ASSERT_EQ(stats.incomplete_points, (size_t)16);
    ASSERT_EQ(stats.point_coverage.truncated_fallbacks, (size_t)16);
    ASSERT_EQ(components->component_count, (size_t)0);

    alea_plot_error_components_free(components);
    alea_destroy(sys);
}

TEST(grid_hex_lattice_coherent_walk_resolves_child_cells) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_hex_lattice.mcnp");
    if (!model) return;
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    int cell_ids[3];
    uint8_t errors[3];
    alea_slice_view_t view;
    /* Pixel centers are x = -2, 0, 2 at y = z = 0. The latter two cover
     * the canonical hex-lattice regression points in one coherent row. */
    alea_slice_view_axis(&view, 2, 0.0, -3.0, 3.0, -1.0, 1.0);
    ASSERT_EQ(alea_find_cells_grid(sys, &view, 3, 1, -1,
                                   cell_ids, NULL, errors), 0);
    ASSERT_EQ(cell_ids[1], 1);
    ASSERT_EQ(cell_ids[2], 3);
    ASSERT_EQ(errors[0], 0);
    ASSERT_EQ(errors[1], 0);
    ASSERT_EQ(errors[2], 0);

    mcnp_model_destroy(model);
}

TEST(grid_boundary_filter_suppresses_shared_surface) {
    alea_system_t* sys = make_x_split(0.0, 10.0);
    ASSERT_NOT_NULL(sys);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    int cell_ids[1];
    int secondary_ids[1] = {-1};
    uint8_t coverage[1];
    uint8_t errors[1];
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -1.0, 1.0, -1.0, 1.0);

    ASSERT_EQ(alea_find_cells_grid_coverage(
        sys, &view, 1, 1, -1,
        ALEA_GRID_COVERAGE_EXACT | ALEA_GRID_SECONDARY_CELL_IDS,
        cell_ids, NULL, secondary_ids, coverage, errors), 0);
    ASSERT_EQ(coverage[0], ALEA_COVERAGE_MULTI);
    ASSERT_EQ(errors[0], ALEA_GRID_OVERLAP);

    alea_boundary_filter_stats_t stats;
    ASSERT_EQ(alea_filter_grid_boundary_ambiguities(
        sys, &view, 1, 1, -1, cell_ids, secondary_ids,
        coverage, errors, &stats), 1);
    ASSERT_EQ(stats.checked, 1);
    ASSERT_EQ(stats.suppressed, 1);
    ASSERT_EQ(coverage[0], ALEA_COVERAGE_ONE);
    ASSERT_EQ(errors[0], ALEA_GRID_OK);
    ASSERT_EQ(secondary_ids[0], -1);

    alea_destroy(sys);
}

TEST(grid_fast_suppresses_coplanar_split_surface_overlap) {
    alea_system_t* sys = make_z_split_sphere(5.0);
    ASSERT_NOT_NULL(sys);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    const int nu = 199, nv = 199;
    const int n = nu * nv;
    int* cell_ids = calloc((size_t)n, sizeof(int));
    uint8_t* coverage = calloc((size_t)n, sizeof(uint8_t));
    uint8_t* errors = calloc((size_t)n, sizeof(uint8_t));
    ASSERT_NOT_NULL(cell_ids);
    ASSERT_NOT_NULL(coverage);
    ASSERT_NOT_NULL(errors);

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -12.0, 12.0, -12.0, 12.0);

    ASSERT_EQ(alea_find_cells_grid_coverage(
        sys, &view, nu, nv, -1, ALEA_GRID_COVERAGE_FAST,
        cell_ids, NULL, NULL, coverage, errors), 0);
    ASSERT_EQ(count_coverage(coverage, n, ALEA_COVERAGE_MULTI), 0);
    ASSERT_EQ(count_errors(errors, n), 0);

    alea_slice_view_axis(&view, 2, 2.0, -12.0, 12.0, -12.0, 12.0);
    ASSERT_EQ(alea_find_cells_grid_coverage(
        sys, &view, nu, nv, -1, ALEA_GRID_COVERAGE_FAST,
        cell_ids, NULL, NULL, coverage, errors), 0);
    ASSERT_EQ(count_coverage(coverage, n, ALEA_COVERAGE_MULTI), 0);
    ASSERT_EQ(count_errors(errors, n), 0);

    free(errors);
    free(coverage);
    free(cell_ids);
    alea_destroy(sys);
}

TEST(grid_boundary_filter_keeps_finite_overlap) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int s1 = alea_sphere_surface(sys, 1, -0.25, 0.0, 0.0, 1.0);
    int s2 = alea_sphere_surface(sys, 2,  0.25, 0.0, 0.0, 1.0);
    ASSERT(s1 >= 0);
    ASSERT(s2 >= 0);
    alea_node_id_t n1 = alea_halfspace(sys, s1, -1);
    alea_node_id_t n2 = alea_halfspace(sys, s2, -1);
    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);
    alea_add_cell(sys, 1, n1, m1, -1.0, 0);
    alea_add_cell(sys, 2, n2, m2, -1.0, 0);
    alea_build_universe_index(sys);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    int cell_ids[1];
    int secondary_ids[1] = {-1};
    uint8_t coverage[1];
    uint8_t errors[1];
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -1.0, 1.0, -1.0, 1.0);

    ASSERT_EQ(alea_find_cells_grid_coverage(
        sys, &view, 1, 1, -1,
        ALEA_GRID_COVERAGE_EXACT | ALEA_GRID_SECONDARY_CELL_IDS,
        cell_ids, NULL, secondary_ids, coverage, errors), 0);
    ASSERT_EQ(coverage[0], ALEA_COVERAGE_MULTI);
    ASSERT_EQ(errors[0], ALEA_GRID_OVERLAP);

    alea_boundary_filter_stats_t stats;
    ASSERT_EQ(alea_filter_grid_boundary_ambiguities(
        sys, &view, 1, 1, -1, cell_ids, secondary_ids,
        coverage, errors, &stats), 0);
    ASSERT_EQ(stats.checked, 1);
    ASSERT_EQ(stats.retained, 1);
    ASSERT_EQ(coverage[0], ALEA_COVERAGE_MULTI);
    ASSERT_EQ(errors[0], ALEA_GRID_OVERLAP);

    alea_destroy(sys);
}

TEST(grid_boundary_filter_keeps_nested_overlap) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int outer = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 10.0);
    int inner = alea_sphere_surface(sys, 2, 0.0, 0.0, 0.0, 3.0);
    ASSERT(outer >= 0);
    ASSERT(inner >= 0);

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);
    alea_add_cell(sys, 1, alea_halfspace(sys, outer, -1), m1, -1.0, 0);
    alea_add_cell(sys, 2, alea_halfspace(sys, inner, -1), m2, -2.0, 0);
    alea_build_universe_index(sys);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    int cell_ids[1];
    int secondary_ids[1] = {-1};
    uint8_t coverage[1];
    uint8_t errors[1];
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -1.0, 1.0, -1.0, 1.0);

    ASSERT_EQ(alea_find_cells_grid_coverage(
        sys, &view, 1, 1, -1,
        ALEA_GRID_COVERAGE_EXACT | ALEA_GRID_SECONDARY_CELL_IDS,
        cell_ids, NULL, secondary_ids, coverage, errors), 0);
    ASSERT_EQ(coverage[0], ALEA_COVERAGE_MULTI);
    ASSERT_EQ(errors[0], ALEA_GRID_OVERLAP);

    alea_boundary_filter_stats_t stats;
    ASSERT_EQ(alea_filter_grid_boundary_ambiguities(
        sys, &view, 1, 1, -1, cell_ids, secondary_ids,
        coverage, errors, &stats), 0);
    ASSERT_EQ(stats.checked, 1);
    ASSERT_EQ(stats.retained, 1);
    ASSERT_EQ(coverage[0], ALEA_COVERAGE_MULTI);
    ASSERT_EQ(errors[0], ALEA_GRID_OVERLAP);
    ASSERT(secondary_ids[0] > 0);

    alea_destroy(sys);
}

TEST(grid_local_coverage_returns_compact_components_with_budget) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int outer = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 10.0);
    int inner = alea_sphere_surface(sys, 2, 0.0, 0.0, 0.0, 3.0);
    ASSERT(outer >= 0 && inner >= 0);
    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);
    alea_add_cell(sys, 1, alea_halfspace(sys, outer, -1), m1, -1.0, 0);
    alea_add_cell(sys, 2, alea_halfspace(sys, inner, -1), m2, -2.0, 0);
    alea_build_universe_index(sys);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0.0, -12.0, 12.0, -12.0, 12.0);
    alea_plot_error_component_result_t* components = NULL;
    alea_local_coverage_stats_t stats;
    ASSERT_EQ(alea_find_local_coverage_components(
        sys, &view, 32, 32, -1, 1024, 16 * 1024 * 1024, 1,
        &components, &stats), 0);
    ASSERT_NOT_NULL(components);
    ASSERT(count_components_of_kind(components, ALEA_PLOT_ERR_TOTAL_OVERLAP) > 0);
    ASSERT_EQ((int)stats.pixels, 1024);
    ASSERT_EQ((int)stats.scratch_bytes, 15 * 1024);
    ASSERT_EQ(stats.worker_limit, 1);
    ASSERT(stats.point_coverage.queries > 0);
    ASSERT(components->components[0].representative_i >= 0);
    ASSERT(components->components[0].representative_j >= 0);
    /* The compact query must preserve the exact global-grid verdict inside
     * the same local domain; only its raster publication is omitted. */
    int cell_ids[1024], secondary_ids[1024];
    uint8_t coverage[1024], errors[1024];
    ASSERT_EQ(alea_find_cells_grid_coverage(
        sys, &view, 32, 32, -1,
        ALEA_GRID_COVERAGE_EXACT | ALEA_GRID_SECONDARY_CELL_IDS,
        cell_ids, NULL, secondary_ids, coverage, errors), 0);
    ASSERT_EQ(alea_filter_grid_boundary_ambiguities(
        sys, &view, 32, 32, -1, cell_ids, secondary_ids, coverage, errors,
        NULL), 0);
    alea_plot_error_component_result_t* global_components =
        alea_classify_plot_error_components(cell_ids, secondary_ids, coverage, 32, 32);
    ASSERT_NOT_NULL(global_components);
    ASSERT_EQ((int)components->component_count, (int)global_components->component_count);
    ASSERT_EQ(count_components_of_kind(components, ALEA_PLOT_ERR_TOTAL_OVERLAP),
              count_components_of_kind(global_components, ALEA_PLOT_ERR_TOTAL_OVERLAP));
    alea_plot_error_components_free(global_components);
    alea_plot_error_components_free(components);

    components = NULL;
    ASSERT_EQ(alea_find_local_coverage_components(
        sys, &view, 32, 32, -1, 1023, 16 * 1024 * 1024, 1,
        &components, &stats), ALEA_LOCAL_COVERAGE_BUDGET_EXCEEDED);
    ASSERT_NULL(components);

    alea_interrupt();
    ASSERT_EQ(alea_find_local_coverage_components(
        sys, &view, 32, 32, -1, 1024, 16 * 1024 * 1024, 1,
        &components, &stats), -1);
    ASSERT_NULL(components);
    alea_clear_interrupt();
    alea_destroy(sys);
}

TEST_MAIN()
