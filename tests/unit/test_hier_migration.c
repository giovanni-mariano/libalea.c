// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_hier_migration.c - Tests that all migrated APIs work in hierarchical
 * spatial mode and do not build the flat spatial index as a side effect.
 *
 * Coverage:
 *   - alea_get_slice_curves():  flat/hier parity, no flat index in hier mode
 *   - alea_find_cells_grid():   hier mode on fill model
 *   - alea_mesh_sample():       hier mode, no flat index built
 *   - render3d batch render:    hier mode, no flat index built
 */

#define _POSIX_C_SOURCE 200112L
#define ALEA_TEST_IMPLEMENTATION
#include "alea_test.h"
#include "alea.h"
#include "alea_mcnp.h"
#include "alea_slice.h"
#include "alea_mesh.h"
#include "render/render3d.h"
#include "core/alea_system.h"
#include "core/alea_universe.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ============================================================================
 * Helpers
 * ============================================================================ */

/* Build a simple sphere scene programmatically. */
static alea_system_t* make_sphere_sys(void) {
    alea_system_t* sys = alea_create();
    if (!sys) return NULL;
    int si = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    alea_node_id_t sph = alea_surface_at(sys, si)->neg_node;
    int m1 = alea_add_material(sys, 1);
    alea_add_cell(sys, 1, sph, m1, -2.7, 0);
    return sys;
}

/* ============================================================================
 * Slice curve tests
 * ============================================================================ */

/* Flat and hier modes produce the same number of curves and the same types. */
TEST(slice_curves_flat_hier_parity) {
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0, -10, 10, -10, 10);

    /* --- flat --- */
    alea_system_t* flat = make_sphere_sys();
    ASSERT_NOT_NULL(flat);
    alea_config_t cfg = alea_get_config(flat);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_FLAT;
    alea_set_config(flat, &cfg);
    alea_build_universe_index(flat);

    alea_slice_curves_t* flat_c = alea_get_slice_curves(flat, &view);
    ASSERT_NOT_NULL(flat_c);
    size_t flat_count = alea_slice_curves_count(flat_c);
    ASSERT(flat_count > 0);

    /* --- hier --- */
    alea_system_t* hier = make_sphere_sys();
    ASSERT_NOT_NULL(hier);
    cfg = alea_get_config(hier);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_HIER;
    alea_set_config(hier, &cfg);
    alea_build_universe_index(hier);

    alea_slice_curves_t* hier_c = alea_get_slice_curves(hier, &view);
    ASSERT_NOT_NULL(hier_c);
    size_t hier_count = alea_slice_curves_count(hier_c);

    ASSERT_EQ((int)flat_count, (int)hier_count);

    /* Same curve type for each curve */
    for (size_t i = 0; i < flat_count; i++) {
        alea_curve_t fc, hc;
        ASSERT_EQ(alea_slice_curves_get(flat_c, i, &fc), 0);
        ASSERT_EQ(alea_slice_curves_get(hier_c, i, &hc), 0);
        ASSERT_EQ(fc.type, hc.type);
        if (fc.type == ALEA_CURVE_CIRCLE) {
            ASSERT_NEAR(fc.data.circle.radius, hc.data.circle.radius, 0.01);
        }
    }

    alea_slice_curves_free(flat_c);
    alea_slice_curves_free(hier_c);
    alea_destroy(flat);
    alea_destroy(hier);
}

/* alea_get_slice_curves() in hier mode must NOT build the flat spatial index. */
TEST(slice_curves_hier_no_flat_spatial_index) {
    alea_system_t* sys = make_sphere_sys();
    ASSERT_NOT_NULL(sys);

    alea_config_t cfg = alea_get_config(sys);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_HIER;
    alea_set_config(sys, &cfg);
    alea_build_universe_index(sys);

    ASSERT_NULL(sys->spatial_index);

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0, -10, 10, -10, 10);
    alea_slice_curves_t* curves = alea_get_slice_curves(sys, &view);
    ASSERT_NOT_NULL(curves);
    ASSERT(alea_slice_curves_count(curves) > 0);

    ASSERT_NULL(sys->spatial_index);
    ASSERT_NOT_NULL(sys->hier_spatial_index);

    alea_slice_curves_free(curves);
    alea_destroy(sys);
}

/* alea_get_slice_curves() in auto mode must NOT build the flat spatial index
 * (auto selects hier for all models with this threshold override). */
TEST(slice_curves_auto_no_flat_spatial_index) {
    alea_setenv("ALEA_SPATIAL_AUTO_CELL_THRESHOLD", "1", 1);

    alea_system_t* sys = make_sphere_sys();
    ASSERT_NOT_NULL(sys);

    alea_config_t cfg = alea_get_config(sys);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_AUTO;
    alea_set_config(sys, &cfg);
    alea_build_universe_index(sys);

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0, -10, 10, -10, 10);
    alea_slice_curves_t* curves = alea_get_slice_curves(sys, &view);
    ASSERT_NOT_NULL(curves);
    ASSERT(alea_slice_curves_count(curves) > 0);

    ASSERT_NULL(sys->spatial_index);

    alea_slice_curves_free(curves);
    alea_destroy(sys);
    alea_unsetenv("ALEA_SPATIAL_AUTO_CELL_THRESHOLD");
}

/* alea_check_slice_errors() (the --errors overlay path) works in hier mode. */
TEST(slice_errors_hier_mode) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Two non-overlapping spheres */
    int s1 = alea_sphere_surface(sys, 1, 0, 0, 0, 3.0);
    int s2 = alea_sphere_surface(sys, 2, 8, 0, 0, 3.0);
    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);
    alea_add_cell(sys, 1, alea_surface_at(sys, s1)->neg_node, m1, -2.7, 0);
    alea_add_cell(sys, 2, alea_surface_at(sys, s2)->neg_node, m2, -8.0, 0);

    alea_config_t cfg = alea_get_config(sys);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_HIER;
    alea_set_config(sys, &cfg);
    alea_build_universe_index(sys);

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0, -6, 14, -6, 6);

    alea_slice_curves_t* curves = alea_get_slice_curves(sys, &view);
    ASSERT_NOT_NULL(curves);
    ASSERT(alea_slice_curves_count(curves) >= 2);

    /* Check errors — should complete without crash, no overlaps */
    alea_slice_error_result_t* errs = alea_check_slice_errors(sys, &view, curves, -1);
    ASSERT_NOT_NULL(errs);

    size_t overlaps = 0;
    for (size_t i = 0; i < errs->error_count; i++) {
        if (errs->errors[i].type == ALEA_SLICE_ERR_OVERLAP) overlaps++;
    }
    ASSERT_EQ((int)overlaps, 0);

    ASSERT_NULL(sys->spatial_index);

    alea_slice_errors_free(errs);
    alea_slice_curves_free(curves);
    alea_destroy(sys);
}

/* alea_find_surface_label_positions() (the --labels=surfaces path) works in
 * hier mode. */
TEST(surface_labels_hier_mode) {
    alea_system_t* sys = make_sphere_sys();
    ASSERT_NOT_NULL(sys);

    alea_config_t cfg = alea_get_config(sys);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_HIER;
    alea_set_config(sys, &cfg);
    alea_build_universe_index(sys);

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0, -10, 10, -10, 10);

    alea_slice_curves_t* curves = alea_get_slice_curves(sys, &view);
    ASSERT_NOT_NULL(curves);

    alea_label_position_t* labels = NULL;
    int label_count = 0;
    int rc = alea_find_surface_label_positions(curves, -10, 10, -10, 10,
                                               64, 64, 2,
                                               &labels, &label_count);
    ASSERT_EQ(rc, 0);
    ASSERT(label_count > 0);

    ASSERT_NULL(sys->spatial_index);

    free(labels);
    alea_slice_curves_free(curves);
    alea_destroy(sys);
}

/* ============================================================================
 * Grid query tests
 * ============================================================================ */

/* alea_find_cells_grid() in hier mode on a fill model returns correct results
 * and does NOT build the flat spatial index. */
TEST(find_cells_grid_hier_no_flat_index) {
    mcnp_model_t* model = mcnp_load("tests/data/simple_fill.mcnp");
    if (!model) SKIP("test data not found");
    alea_system_t* sys = model->sys;

    alea_config_t cfg = alea_get_config(sys);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_HIER;
    alea_set_config(sys, &cfg);

    ASSERT_NULL(sys->spatial_index);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    ASSERT_NULL(sys->spatial_index);
    ASSERT_NOT_NULL(sys->hier_spatial_index);

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0, -110, 110, -110, 110);

    int W = 64, H = 64;
    int* cell_ids = malloc(W * H * sizeof(int));
    int* mat_ids  = malloc(W * H * sizeof(int));
    ASSERT_NOT_NULL(cell_ids);
    ASSERT_NOT_NULL(mat_ids);

    int rc = alea_find_cells_grid(sys, &view, W, H, -1, cell_ids, mat_ids, NULL);
    ASSERT_EQ(rc, 0);

    /* simple_fill: inner sphere r=5 with mat 1; outer r=100 is void fill.
     * Center pixel should hit mat 1. */
    int cx = W/2, cy = H/2;
    ASSERT_EQ(mat_ids[cy * W + cx], 1);

    /* No flat index was built as a side effect */
    ASSERT_NULL(sys->spatial_index);

    free(cell_ids);
    free(mat_ids);
    mcnp_model_destroy(model);
}

/* Flat and hier grid modes produce the same material layout. */
TEST(find_cells_grid_flat_hier_parity) {
    /* load twice to get two independent systems */
    mcnp_model_t* mf = mcnp_load("tests/data/simple_fill.mcnp");
    mcnp_model_t* mh = mcnp_load("tests/data/simple_fill.mcnp");
    if (!mf || !mh) {
        if (mf) mcnp_model_destroy(mf);
        if (mh) mcnp_model_destroy(mh);
        SKIP("test data not found");
    }

    alea_config_t cfg;

    cfg = alea_get_config(mf->sys);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_FLAT;
    alea_set_config(mf->sys, &cfg);
    alea_prepare_query_acceleration(mf->sys);

    cfg = alea_get_config(mh->sys);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_HIER;
    alea_set_config(mh->sys, &cfg);
    alea_prepare_query_acceleration(mh->sys);

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0, -110, 110, -110, 110);
    int W = 32, H = 32;
    int* flat_cell = malloc(W * H * sizeof(int));
    int* hier_cell = malloc(W * H * sizeof(int));
    int* flat_mat  = malloc(W * H * sizeof(int));
    int* hier_mat  = malloc(W * H * sizeof(int));
    ASSERT_NOT_NULL(flat_cell); ASSERT_NOT_NULL(hier_cell);
    ASSERT_NOT_NULL(flat_mat);  ASSERT_NOT_NULL(hier_mat);

    ASSERT_EQ(alea_find_cells_grid(mf->sys, &view, W, H, -1, flat_cell, flat_mat, NULL), 0);
    ASSERT_EQ(alea_find_cells_grid(mh->sys, &view, W, H, -1, hier_cell, hier_mat, NULL), 0);

    for (int i = 0; i < W * H; i++) {
        ASSERT_EQ(flat_mat[i], hier_mat[i]);
    }

    free(flat_cell); free(hier_cell);
    free(flat_mat);  free(hier_mat);
    mcnp_model_destroy(mf);
    mcnp_model_destroy(mh);
}

/* ============================================================================
 * Mesh export tests
 * ============================================================================ */

/* alea_mesh_sample() in hier mode produces correct results and does NOT build
 * the flat spatial index. */
TEST(mesh_hier_no_flat_spatial_index) {
    alea_system_t* sys = make_sphere_sys();
    ASSERT_NOT_NULL(sys);

    alea_config_t cfg = alea_get_config(sys);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_HIER;
    alea_set_config(sys, &cfg);
    alea_build_universe_index(sys);

    ASSERT_NULL(sys->spatial_index);

    alea_mesh_config_t mcfg;
    alea_mesh_config_init(&mcfg);
    mcfg.nx = mcfg.ny = mcfg.nz = 8;
    mcfg.x_min = mcfg.y_min = mcfg.z_min = -10.0;
    mcfg.x_max = mcfg.y_max = mcfg.z_max =  10.0;

    alea_mesh_result_t* mesh = alea_mesh_sample(sys, &mcfg);
    ASSERT_NOT_NULL(mesh);

    /* Center voxel should be inside sphere (mat 1) */
    int ci = 4 * 8 * 8 + 4 * 8 + 4;
    ASSERT_EQ(mesh->material_ids[ci], 1);

    ASSERT_NULL(sys->spatial_index);

    alea_mesh_result_free(mesh);
    alea_destroy(sys);
}

/* Flat and hier mesh sampling produce the same material layout. */
TEST(mesh_flat_hier_parity) {
    alea_mesh_config_t mcfg;
    alea_mesh_config_init(&mcfg);
    mcfg.nx = mcfg.ny = mcfg.nz = 6;
    mcfg.x_min = mcfg.y_min = mcfg.z_min = -8.0;
    mcfg.x_max = mcfg.y_max = mcfg.z_max =  8.0;

    alea_system_t* flat = make_sphere_sys();
    alea_system_t* hier = make_sphere_sys();
    ASSERT_NOT_NULL(flat);
    ASSERT_NOT_NULL(hier);

    alea_config_t cfg;
    cfg = alea_get_config(flat);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_FLAT;
    alea_set_config(flat, &cfg);
    alea_build_universe_index(flat);

    cfg = alea_get_config(hier);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_HIER;
    alea_set_config(hier, &cfg);
    alea_build_universe_index(hier);

    alea_mesh_result_t* rm_flat = alea_mesh_sample(flat, &mcfg);
    alea_mesh_result_t* rm_hier = alea_mesh_sample(hier, &mcfg);
    ASSERT_NOT_NULL(rm_flat);
    ASSERT_NOT_NULL(rm_hier);

    int n = 6 * 6 * 6;
    for (int i = 0; i < n; i++) {
        ASSERT_EQ(rm_flat->material_ids[i], rm_hier->material_ids[i]);
    }

    alea_mesh_result_free(rm_flat);
    alea_mesh_result_free(rm_hier);
    alea_destroy(flat);
    alea_destroy(hier);
}

/* ============================================================================
 * Render3D tests
 * ============================================================================ */

static render_framebuffer_t* do_render(alea_system_t* sys, int w, int h) {
    render_config_t rcfg;
    render_config_init(&rcfg);
    rcfg.width  = w;
    rcfg.height = h;
    rcfg.render_mode = RENDER_MODE_SOLID;

    render_camera_t cam;
    if (render_camera_setup(&cam, &rcfg, sys) != 0) {
        render_config_free(&rcfg);
        return NULL;
    }

    render_framebuffer_t* fb = render_framebuffer_create(w, h, 0);
    if (!fb) { render_config_free(&rcfg); return NULL; }

    if (render_scene(sys, &rcfg, &cam, fb) != 0) {
        render_framebuffer_free(fb);
        render_config_free(&rcfg);
        return NULL;
    }
    render_config_free(&rcfg);
    return fb;
}

/* render3d batch render in hier mode does NOT build the flat spatial index. */
TEST(render3d_hier_no_flat_spatial_index) {
    alea_system_t* sys = make_sphere_sys();
    ASSERT_NOT_NULL(sys);

    alea_config_t cfg = alea_get_config(sys);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_HIER;
    alea_set_config(sys, &cfg);
    alea_build_universe_index(sys);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    ASSERT_NULL(sys->spatial_index);
    ASSERT_NOT_NULL(sys->hier_spatial_index);

    render_framebuffer_t* fb = do_render(sys, 16, 16);
    ASSERT_NOT_NULL(fb);

    /* At least some non-background pixels (sphere should be visible) */
    int hit = 0;
    for (int i = 0; i < 16 * 16; i++) {
        if (fb->cell_id[i] > 0) hit++;
    }
    ASSERT(hit > 0);

    ASSERT_NULL(sys->spatial_index);

    render_framebuffer_free(fb);
    alea_destroy(sys);
}

/* Flat and hier renders produce identical pixel cell IDs on a simple scene. */
TEST(render3d_flat_hier_parity) {
    alea_system_t* flat = make_sphere_sys();
    alea_system_t* hier = make_sphere_sys();
    ASSERT_NOT_NULL(flat);
    ASSERT_NOT_NULL(hier);

    alea_config_t cfg;
    cfg = alea_get_config(flat);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_FLAT;
    alea_set_config(flat, &cfg);
    alea_build_universe_index(flat);
    alea_prepare_query_acceleration(flat);

    cfg = alea_get_config(hier);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_HIER;
    alea_set_config(hier, &cfg);
    alea_build_universe_index(hier);
    alea_prepare_query_acceleration(hier);

    render_framebuffer_t* fb_flat = do_render(flat, 16, 16);
    render_framebuffer_t* fb_hier = do_render(hier, 16, 16);
    ASSERT_NOT_NULL(fb_flat);
    ASSERT_NOT_NULL(fb_hier);

    for (int i = 0; i < 16 * 16; i++) {
        ASSERT_EQ(fb_flat->cell_id[i], fb_hier->cell_id[i]);
    }

    render_framebuffer_free(fb_flat);
    render_framebuffer_free(fb_hier);
    alea_destroy(flat);
    alea_destroy(hier);
}

TEST_MAIN()
