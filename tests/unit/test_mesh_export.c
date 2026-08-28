// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_mesh_export.c - Unit tests for mesh export module
 */

#include "alea_test.h"
#include "alea.h"
#include "alea_mesh.h"
#include "core/alea_system.h"
#include "core/alea_universe.h"
#include "util/compat.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

/* ============================================================================
 * Helper: create a sphere-in-box scene
 *   Cell 1: sphere at origin, radius 5, material 1
 *   Void outside
 * ============================================================================ */

static alea_system_t *create_sphere_scene(void) {
    alea_system_t *sys = alea_create();
    if (!sys) return NULL;

    int s1 = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    alea_node_id_t sphere = alea_surface_at(sys, s1)->neg_node;

    int m1 = alea_add_material(sys, 1);

    alea_add_cell(sys, 1, sphere, m1, -2.7, 0);

    alea_build_universe_index(sys);
    return sys;
}

static alea_system_t *create_small_inclusion_scene(void) {
    alea_system_t *sys = alea_create();
    if (!sys) return NULL;

    int s1 = alea_sphere_surface(sys, 1, 0.25, 0.25, 0.25, 0.20);
    alea_node_id_t sphere = alea_surface_at(sys, s1)->neg_node;
    int m1 = alea_add_material(sys, 1);

    alea_add_cell(sys, 1, sphere, m1, -1.0, 0);
    alea_build_universe_index(sys);
    return sys;
}

static alea_system_t *create_tied_split_scene(void) {
    alea_system_t *sys = alea_create();
    if (!sys) return NULL;

    int plane = alea_plane_surface(sys, 1, 1.0, 0.0, 0.0, 0.0);
    if (plane < 0) { alea_destroy(sys); return NULL; }

    int mat2 = alea_add_material(sys, 2);
    int mat1 = alea_add_material(sys, 1);
    alea_add_cell(sys, 10, alea_halfspace(sys, plane, -1), mat2, -1.0, 0);
    alea_add_cell(sys, 20, alea_halfspace(sys, plane,  1), mat1, -1.0, 0);
    alea_build_universe_index(sys);
    return sys;
}

static alea_system_t *create_tied_cell_scene(void) {
    alea_system_t *sys = alea_create();
    if (!sys) return NULL;
    int plane = alea_plane_surface(sys, 1, 1.0, 0.0, 0.0, 0.0);
    if (plane < 0) { alea_destroy(sys); return NULL; }
    int mat = alea_add_material(sys, 7);
    alea_add_cell(sys, 30, alea_halfspace(sys, plane, -1), mat, -1.0, 0);
    alea_add_cell(sys, 10, alea_halfspace(sys, plane,  1), mat, -1.0, 0);
    alea_build_universe_index(sys);
    return sys;
}

static double mesh_fraction_for_material(const alea_mesh_result_t *mesh,
                                         size_t voxel_index,
                                         int material_id) {
    alea_mesh_fraction_span_t span = mesh->fraction_spans[voxel_index];
    for (uint32_t i = 0; i < span.count; i++) {
        const alea_mesh_material_fraction_t *f =
            &mesh->fractions[(size_t)span.offset + (size_t)i];
        if (f->material_id == material_id)
            return f->fraction;
    }
    return 0.0;
}

/* ============================================================================
 * Tests
 * ============================================================================ */

TEST(mesh_config_defaults) {
    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);

    ASSERT_EQ(cfg.nx, 10);
    ASSERT_EQ(cfg.ny, 10);
    ASSERT_EQ(cfg.nz, 10);
    ASSERT_EQ(cfg.format, ALEA_MESH_GMSH);
    ASSERT_EQ(cfg.void_material_id, 0);
    ASSERT_NEAR(cfg.auto_pad, 0.01, 1e-15);
    ASSERT_EQ(cfg.sampling_mode, ALEA_MESH_SAMPLE_SUBCELL);
    ASSERT_EQ(cfg.subsamples_per_axis, 2);
    ASSERT_NEAR(cfg.mixed_threshold, 0.0, 1e-15);
    ASSERT_NEAR(cfg.x_min, 0.0, 1e-15);
    ASSERT_NEAR(cfg.x_max, 0.0, 1e-15);
    ASSERT_NULL(cfg.x_nodes);
    ASSERT_NULL(cfg.y_nodes);
    ASSERT_NULL(cfg.z_nodes);
}

TEST(mesh_sample_sphere) {
    alea_system_t *sys = create_sphere_scene();
    ASSERT_NOT_NULL(sys);

    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    cfg.nx = cfg.ny = cfg.nz = 10;
    cfg.x_min = cfg.y_min = cfg.z_min = -10.0;
    cfg.x_max = cfg.y_max = cfg.z_max =  10.0;

    alea_mesh_result_t *mesh = alea_mesh_sample(sys, &cfg);
    ASSERT_NOT_NULL(mesh);
    ASSERT_EQ(mesh->nx, 10);
    ASSERT_EQ(mesh->ny, 10);
    ASSERT_EQ(mesh->nz, 10);

    /* Center voxel (5,5,5) should have material 1 */
    int center_idx = 5 * 10 * 10 + 5 * 10 + 5;
    ASSERT_EQ(mesh->material_ids[center_idx], 1);
    ASSERT_EQ(mesh->cell_ids[center_idx], 1);

    /* Corner voxel (0,0,0) center is at (-9,-9,-9) — outside sphere → void */
    int corner_idx = 0;
    ASSERT_EQ(mesh->material_ids[corner_idx], 0);
    ASSERT_EQ(mesh->cell_ids[corner_idx], -1);

    /* unique materials should include 0 and 1 */
    ASSERT(mesh->num_materials >= 2);

    alea_mesh_result_free(mesh);
    alea_destroy(sys);
}

TEST(mesh_fractions_populated_by_default) {
    alea_system_t *sys = create_sphere_scene();
    ASSERT_NOT_NULL(sys);

    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    cfg.nx = cfg.ny = cfg.nz = 2;
    cfg.x_min = cfg.y_min = cfg.z_min = -6.0;
    cfg.x_max = cfg.y_max = cfg.z_max =  6.0;

    alea_mesh_result_t *mesh = alea_mesh_sample(sys, &cfg);
    ASSERT_NOT_NULL(mesh);
    ASSERT_NOT_NULL(mesh->mixed_flags);
    ASSERT_NOT_NULL(mesh->dominant_fractions);
    ASSERT_NOT_NULL(mesh->sample_counts);
    ASSERT_NOT_NULL(mesh->tie_flags);
    ASSERT_NOT_NULL(mesh->fraction_spans);
    ASSERT_NOT_NULL(mesh->fractions);
    ASSERT(mesh->fraction_count > 0);
    ASSERT(mesh->dominant_fractions[0] > 0.0);
    ASSERT(mesh->dominant_fractions[0] <= 1.0);
    ASSERT_EQ(mesh->sample_counts[0], 8);

    alea_mesh_result_free(mesh);
    alea_destroy(sys);
}

TEST(mesh_subsample_detects_inclusion_missed_by_center) {
    alea_system_t *sys = create_small_inclusion_scene();
    ASSERT_NOT_NULL(sys);

    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    cfg.nx = cfg.ny = cfg.nz = 1;
    cfg.x_min = cfg.y_min = cfg.z_min = 0.0;
    cfg.x_max = cfg.y_max = cfg.z_max = 1.0;
    cfg.sampling_mode = ALEA_MESH_SAMPLE_SUBCELL;
    cfg.subsamples_per_axis = 2;

    alea_mesh_result_t *mesh = alea_mesh_sample(sys, &cfg);
    ASSERT_NOT_NULL(mesh);

    /* The voxel center is void, but a 2x2x2 sample records the inclusion. */
    ASSERT_EQ(mesh->material_ids[0], 0);
    ASSERT_NOT_NULL(mesh->mixed_flags);
    ASSERT_NOT_NULL(mesh->dominant_fractions);
    ASSERT_NOT_NULL(mesh->fraction_spans);
    ASSERT_NOT_NULL(mesh->fractions);
    ASSERT_EQ(mesh->mixed_flags[0], 1);
    ASSERT_EQ(mesh->mixed_count, 1);
    ASSERT_EQ(mesh->fraction_spans[0].count, 2);
    ASSERT_EQ(mesh->fraction_count, 2);
    ASSERT_EQ(mesh->num_materials, 2);
    ASSERT_EQ(mesh->unique_materials[0], 0);
    ASSERT_EQ(mesh->unique_materials[1], 1);
    ASSERT_NEAR(mesh_fraction_for_material(mesh, 0, 1), 0.125, 1e-15);
    ASSERT_NEAR(mesh_fraction_for_material(mesh, 0, 0), 0.875, 1e-15);
    ASSERT_NEAR(mesh->dominant_fractions[0], 0.875, 1e-15);

    FILE *f = tmpfile();
    ASSERT_NOT_NULL(f);
    int rc = alea_mesh_export_stream(mesh, ALEA_MESH_VTK, f);
    ASSERT_EQ(rc, 0);
    rewind(f);
    char buf[4096];
    size_t n = fread(buf, 1, sizeof(buf) - 1, f);
    buf[n] = '\0';
    fclose(f);
    ASSERT(strstr(buf, "SCALARS mixed_flag int") != NULL);
    ASSERT(strstr(buf, "SCALARS dominant_sampled_fraction double") != NULL);
    ASSERT(strstr(buf, "SCALARS sample_count int") != NULL);
    ASSERT(strstr(buf, "SCALARS tie_flag int") != NULL);

    alea_mesh_result_free(mesh);
    alea_destroy(sys);
}

TEST(mesh_tie_is_explicit_and_stable) {
    alea_system_t *sys = create_tied_split_scene();
    ASSERT_NOT_NULL(sys);

    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    cfg.nx = cfg.ny = cfg.nz = 1;
    cfg.x_min = cfg.y_min = cfg.z_min = -1.0;
    cfg.x_max = cfg.y_max = cfg.z_max = 1.0;
    cfg.sampling_mode = ALEA_MESH_SAMPLE_SUBCELL;
    cfg.subsamples_per_axis = 2;

    alea_mesh_result_t *mesh = alea_mesh_sample(sys, &cfg);
    ASSERT_NOT_NULL(mesh);
    ASSERT_EQ(mesh->sample_counts[0], 8);
    ASSERT(mesh->tie_flags[0] & ALEA_MESH_TIE_MATERIAL);
    ASSERT_EQ(mesh->material_ids[0], 1);
    ASSERT_EQ(mesh->cell_ids[0], 20);
    ASSERT_NEAR(mesh->dominant_fractions[0], 0.5, 1e-15);

    alea_mesh_result_free(mesh);
    alea_destroy(sys);
}

TEST(mesh_cell_tie_is_explicit_and_stable) {
    alea_system_t *sys = create_tied_cell_scene();
    ASSERT_NOT_NULL(sys);

    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    cfg.nx = cfg.ny = cfg.nz = 1;
    cfg.x_min = cfg.y_min = cfg.z_min = -1.0;
    cfg.x_max = cfg.y_max = cfg.z_max = 1.0;

    alea_mesh_result_t *mesh = alea_mesh_sample(sys, &cfg);
    ASSERT_NOT_NULL(mesh);
    ASSERT(!(mesh->tie_flags[0] & ALEA_MESH_TIE_MATERIAL));
    ASSERT(mesh->tie_flags[0] & ALEA_MESH_TIE_CELL);
    ASSERT_EQ(mesh->material_ids[0], 7);
    ASSERT_EQ(mesh->cell_ids[0], 10);

    alea_mesh_result_free(mesh);
    alea_destroy(sys);
}

TEST(mesh_rejects_invalid_dimensions) {
    alea_system_t *sys = create_sphere_scene();
    ASSERT_NOT_NULL(sys);

    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    cfg.nx = 0;
    cfg.x_min = cfg.y_min = cfg.z_min = -1.0;
    cfg.x_max = cfg.y_max = cfg.z_max =  1.0;

    alea_mesh_result_t *mesh = alea_mesh_sample(sys, &cfg);
    ASSERT_NULL(mesh);

    cfg.nx = 2;
    cfg.sampling_mode = ALEA_MESH_SAMPLE_SUBCELL;
    cfg.subsamples_per_axis = 9;
    mesh = alea_mesh_sample(sys, &cfg);
    ASSERT_NULL(mesh);

    cfg.subsamples_per_axis = 2;
    cfg.sampling_mode = (alea_mesh_sampling_mode_t)99;
    mesh = alea_mesh_sample(sys, &cfg);
    ASSERT_NULL(mesh);

    alea_destroy(sys);
}

TEST(mesh_rejects_invalid_bounds) {
    alea_system_t *sys = create_sphere_scene();
    ASSERT_NOT_NULL(sys);

    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    cfg.x_min = 1.0;  cfg.x_max = -1.0;
    cfg.y_min = -1.0; cfg.y_max =  1.0;
    cfg.z_min = -1.0; cfg.z_max =  1.0;

    alea_mesh_result_t *mesh = alea_mesh_sample(sys, &cfg);
    ASSERT_NULL(mesh);

    cfg.x_min = -1.0; cfg.x_max = 1.0;
    cfg.y_min = NAN;  cfg.y_max = 1.0;
    mesh = alea_mesh_sample(sys, &cfg);
    ASSERT_NULL(mesh);

    alea_destroy(sys);
}

TEST(mesh_rejects_invalid_custom_nodes) {
    alea_system_t *sys = create_sphere_scene();
    ASSERT_NOT_NULL(sys);

    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    cfg.nx = cfg.ny = cfg.nz = 2;
    cfg.x_min = cfg.y_min = cfg.z_min = -1.0;
    cfg.x_max = cfg.y_max = cfg.z_max =  1.0;

    double bad_x[] = { -1.0, 0.0, 0.0 };
    cfg.x_nodes = bad_x;
    alea_mesh_result_t *mesh = alea_mesh_sample(sys, &cfg);
    ASSERT_NULL(mesh);

    double nan_x[] = { -1.0, NAN, 1.0 };
    cfg.x_nodes = nan_x;
    mesh = alea_mesh_sample(sys, &cfg);
    ASSERT_NULL(mesh);

    alea_destroy(sys);
}

TEST(mesh_rejects_huge_dimensions_before_allocation) {
    alea_system_t *sys = create_sphere_scene();
    ASSERT_NOT_NULL(sys);

    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    cfg.nx = INT_MAX;
    cfg.ny = 2;
    cfg.nz = 2;
    cfg.x_min = cfg.y_min = cfg.z_min = -1.0;
    cfg.x_max = cfg.y_max = cfg.z_max =  1.0;

    alea_mesh_result_t *mesh = alea_mesh_sample(sys, &cfg);
    ASSERT_NULL(mesh);

    alea_destroy(sys);
}

TEST(mesh_auto_bounds_empty_system_fails) {
    alea_system_t *sys = alea_create();
    ASSERT_NOT_NULL(sys);

    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);

    alea_mesh_result_t *mesh = alea_mesh_sample(sys, &cfg);
    ASSERT_NULL(mesh);

    alea_destroy(sys);
}

TEST(mesh_export_gmsh) {
    alea_system_t *sys = create_sphere_scene();
    ASSERT_NOT_NULL(sys);

    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    cfg.nx = cfg.ny = cfg.nz = 3;
    cfg.x_min = cfg.y_min = cfg.z_min = -6.0;
    cfg.x_max = cfg.y_max = cfg.z_max =  6.0;

    alea_mesh_result_t *mesh = alea_mesh_sample(sys, &cfg);
    ASSERT_NOT_NULL(mesh);

    /* Write to tmpfile and verify header */
    FILE *f = tmpfile();
    ASSERT_NOT_NULL(f);

    int rc = alea_mesh_export_stream(mesh, ALEA_MESH_GMSH, f);
    ASSERT_EQ(rc, 0);

    /* Read back and check key markers */
    rewind(f);
    char buf[4096];
    size_t n = fread(buf, 1, sizeof(buf) - 1, f);
    buf[n] = '\0';
    fclose(f);

    ASSERT(strstr(buf, "$MeshFormat") != NULL);
    ASSERT(strstr(buf, "2.2 0 8") != NULL);
    ASSERT(strstr(buf, "$Nodes") != NULL);
    ASSERT(strstr(buf, "$Elements") != NULL);
    ASSERT(strstr(buf, "$PhysicalNames") != NULL);
    ASSERT(strstr(buf, "3 0 \"material_") == NULL);
    ASSERT(strstr(buf, "3 1 \"material_0\"") != NULL);

    /* Node count: (3+1)^3 = 64 */
    ASSERT(strstr(buf, "\n64\n") != NULL);
    /* Element count: 3^3 = 27 */
    ASSERT(strstr(buf, "\n27\n") != NULL);

    alea_mesh_result_free(mesh);
    alea_destroy(sys);
}

TEST(mesh_export_vtk) {
    alea_system_t *sys = create_sphere_scene();
    ASSERT_NOT_NULL(sys);

    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    cfg.nx = cfg.ny = cfg.nz = 4;
    cfg.x_min = cfg.y_min = cfg.z_min = -6.0;
    cfg.x_max = cfg.y_max = cfg.z_max =  6.0;

    alea_mesh_result_t *mesh = alea_mesh_sample(sys, &cfg);
    ASSERT_NOT_NULL(mesh);

    FILE *f = tmpfile();
    ASSERT_NOT_NULL(f);

    int rc = alea_mesh_export_stream(mesh, ALEA_MESH_VTK, f);
    ASSERT_EQ(rc, 0);

    rewind(f);
    char buf[8192];
    size_t n = fread(buf, 1, sizeof(buf) - 1, f);
    buf[n] = '\0';
    fclose(f);

    ASSERT(strstr(buf, "# vtk DataFile Version 3.0") != NULL);
    ASSERT(strstr(buf, "ASCII") != NULL);
    /* Uniform grid → STRUCTURED_POINTS */
    ASSERT(strstr(buf, "STRUCTURED_POINTS") != NULL);
    ASSERT(strstr(buf, "DIMENSIONS 5 5 5") != NULL);
    ASSERT(strstr(buf, "CELL_DATA 64") != NULL);
    ASSERT(strstr(buf, "SCALARS material_id int") != NULL);
    ASSERT(strstr(buf, "SCALARS cell_id int") != NULL);

    alea_mesh_result_free(mesh);
    alea_destroy(sys);
}

TEST(mesh_auto_bounds) {
    alea_system_t *sys = create_sphere_scene();
    ASSERT_NOT_NULL(sys);

    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    cfg.nx = cfg.ny = cfg.nz = 5;
    /* Leave bounds at 0 → auto-detect */

    alea_mesh_result_t *mesh = alea_mesh_sample(sys, &cfg);
    ASSERT_NOT_NULL(mesh);

    /* Auto-bounds should encompass the sphere (radius 5) */
    ASSERT(mesh->x_nodes[0] < -4.0);
    ASSERT(mesh->x_nodes[5] > 4.0);
    ASSERT(mesh->y_nodes[0] < -4.0);
    ASSERT(mesh->y_nodes[5] > 4.0);
    ASSERT(mesh->z_nodes[0] < -4.0);
    ASSERT(mesh->z_nodes[5] > 4.0);

    alea_mesh_result_free(mesh);
    alea_destroy(sys);
}

TEST(mesh_non_uniform_grid) {
    alea_system_t *sys = create_sphere_scene();
    ASSERT_NOT_NULL(sys);

    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    cfg.nx = 3;  cfg.ny = 3;  cfg.nz = 3;

    /* Non-uniform X nodes: finer near center */
    double xn[] = { -6.0, -1.0, 1.0, 6.0 };
    double yn[] = { -6.0, -1.0, 1.0, 6.0 };
    double zn[] = { -6.0, -1.0, 1.0, 6.0 };
    cfg.x_nodes = xn;
    cfg.y_nodes = yn;
    cfg.z_nodes = zn;

    alea_mesh_result_t *mesh = alea_mesh_sample(sys, &cfg);
    ASSERT_NOT_NULL(mesh);

    /* Verify node arrays were copied correctly */
    ASSERT_NEAR(mesh->x_nodes[0], -6.0, 1e-15);
    ASSERT_NEAR(mesh->x_nodes[1], -1.0, 1e-15);
    ASSERT_NEAR(mesh->x_nodes[2],  1.0, 1e-15);
    ASSERT_NEAR(mesh->x_nodes[3],  6.0, 1e-15);

    /* Center voxel (1,1,1) should be inside sphere */
    int center = 1 * 3 * 3 + 1 * 3 + 1;
    ASSERT_EQ(mesh->material_ids[center], 1);

    /* Write VTK — should use RECTILINEAR_GRID since non-uniform */
    FILE *f = tmpfile();
    ASSERT_NOT_NULL(f);
    int rc = alea_mesh_export_stream(mesh, ALEA_MESH_VTK, f);
    ASSERT_EQ(rc, 0);

    rewind(f);
    char buf[4096];
    size_t n = fread(buf, 1, sizeof(buf) - 1, f);
    buf[n] = '\0';
    fclose(f);

    ASSERT(strstr(buf, "RECTILINEAR_GRID") != NULL);
    ASSERT(strstr(buf, "X_COORDINATES") != NULL);

    alea_mesh_result_free(mesh);
    alea_destroy(sys);
}

TEST(mesh_export_system_oneshot) {
    alea_system_t *sys = create_sphere_scene();
    ASSERT_NOT_NULL(sys);

    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    cfg.nx = cfg.ny = cfg.nz = 3;
    cfg.x_min = cfg.y_min = cfg.z_min = -6.0;
    cfg.x_max = cfg.y_max = cfg.z_max =  6.0;
    cfg.format = ALEA_MESH_GMSH;

    /* Test the oneshot function by sampling and exporting to stream */
    alea_mesh_result_t *mesh = alea_mesh_sample(sys, &cfg);
    ASSERT_NOT_NULL(mesh);

    FILE *f = tmpfile();
    ASSERT_NOT_NULL(f);

    int rc = alea_mesh_export_stream(mesh, cfg.format, f);
    ASSERT_EQ(rc, 0);

    /* Verify content */
    rewind(f);
    char buf[256];
    ASSERT_NOT_NULL(fgets(buf, sizeof(buf), f));
    ASSERT(strstr(buf, "$MeshFormat") != NULL);
    fclose(f);

    alea_mesh_result_free(mesh);
    alea_destroy(sys);
}

TEST(mesh_sampling_does_not_validate_unused_export_format) {
    alea_system_t *sys = create_sphere_scene();
    ASSERT_NOT_NULL(sys);
    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    cfg.nx = cfg.ny = cfg.nz = 1;
    cfg.x_min = cfg.y_min = cfg.z_min = -1.0;
    cfg.x_max = cfg.y_max = cfg.z_max = 1.0;
    cfg.format = (alea_mesh_format_t)99;
    alea_mesh_result_t *mesh = alea_mesh_sample(sys, &cfg);
    ASSERT_NOT_NULL(mesh);
    alea_mesh_result_free(mesh);
    alea_destroy(sys);
}

TEST(mesh_export_material_fractions_opt_in) {
    alea_system_t *sys = create_small_inclusion_scene();
    ASSERT_NOT_NULL(sys);
    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    cfg.nx = cfg.ny = cfg.nz = 1;
    cfg.x_min = cfg.y_min = cfg.z_min = 0.0;
    cfg.x_max = cfg.y_max = cfg.z_max = 1.0;
    alea_mesh_result_t *mesh = alea_mesh_sample(sys, &cfg);
    ASSERT_NOT_NULL(mesh);

    alea_mesh_export_options_t options;
    alea_mesh_export_options_init(&options);
    options.fields |= ALEA_MESH_EXPORT_MATERIAL_FRACTIONS;

    FILE *f = tmpfile();
    ASSERT_NOT_NULL(f);
    ASSERT_EQ(alea_mesh_export_stream_ex(mesh, ALEA_MESH_VTK, f, &options), 0);
    rewind(f);
    char buf[16384];
    size_t n = fread(buf, 1, sizeof(buf) - 1, f);
    buf[n] = '\0';
    fclose(f);
    ASSERT(strstr(buf, "SCALARS sampled_fraction_material_0 double") != NULL);
    ASSERT(strstr(buf, "SCALARS sampled_fraction_material_1 double") != NULL);

    f = tmpfile();
    ASSERT_NOT_NULL(f);
    ASSERT_EQ(alea_mesh_export_stream_ex(mesh, ALEA_MESH_GMSH, f, &options), 0);
    rewind(f);
    n = fread(buf, 1, sizeof(buf) - 1, f);
    buf[n] = '\0';
    fclose(f);
    ASSERT(strstr(buf, "$ElementData") != NULL);
    ASSERT(strstr(buf, "\"sampled_fraction_material_0\"") != NULL);
    ASSERT(strstr(buf, "\"sampled_fraction_material_1\"") != NULL);

    options.max_fraction_materials = 1;
    f = tmpfile();
    ASSERT_NOT_NULL(f);
    ASSERT_EQ(alea_mesh_export_stream_ex(mesh, ALEA_MESH_VTK, f, &options), -1);
    fclose(f);

    alea_mesh_result_free(mesh);
    alea_destroy(sys);
}

TEST(mesh_export_rejects_malformed_result) {
    alea_system_t *sys = create_sphere_scene();
    ASSERT_NOT_NULL(sys);
    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    cfg.nx = cfg.ny = cfg.nz = 2;
    cfg.x_min = cfg.y_min = cfg.z_min = -6.0;
    cfg.x_max = cfg.y_max = cfg.z_max = 6.0;
    alea_mesh_result_t *mesh = alea_mesh_sample(sys, &cfg);
    ASSERT_NOT_NULL(mesh);

    alea_mesh_result_t malformed = *mesh;
    malformed.material_ids = NULL;
    FILE *f = tmpfile();
    ASSERT_NOT_NULL(f);
    ASSERT_EQ(alea_mesh_export_stream(&malformed, ALEA_MESH_VTK, f), -1);
    fclose(f);

    malformed = *mesh;
    uint32_t original_count = mesh->fraction_spans[0].count;
    malformed.fraction_spans[0].count = UINT32_MAX;
    f = tmpfile();
    ASSERT_NOT_NULL(f);
    ASSERT_EQ(alea_mesh_export_stream(&malformed, ALEA_MESH_GMSH, f), -1);
    fclose(f);
    mesh->fraction_spans[0].count = original_count;

    alea_mesh_result_free(mesh);
    alea_destroy(sys);
}

#ifndef _WIN32
TEST(mesh_export_reports_stream_write_failure) {
    alea_system_t *sys = create_sphere_scene();
    ASSERT_NOT_NULL(sys);
    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    cfg.nx = cfg.ny = cfg.nz = 2;
    cfg.x_min = cfg.y_min = cfg.z_min = -6.0;
    cfg.x_max = cfg.y_max = cfg.z_max = 6.0;
    alea_mesh_result_t *mesh = alea_mesh_sample(sys, &cfg);
    ASSERT_NOT_NULL(mesh);

    FILE *f = fopen("/dev/full", "w");
    ASSERT_NOT_NULL(f);
    ASSERT_EQ(alea_mesh_export_stream(mesh, ALEA_MESH_VTK, f), -1);
    fclose(f);

    alea_mesh_result_free(mesh);
    alea_destroy(sys);
}
#endif

TEST(mesh_filename_export_replaces_only_on_success) {
    alea_system_t *sys = create_sphere_scene();
    ASSERT_NOT_NULL(sys);
    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    cfg.nx = cfg.ny = cfg.nz = 2;
    cfg.x_min = cfg.y_min = cfg.z_min = -6.0;
    cfg.x_max = cfg.y_max = cfg.z_max = 6.0;
    alea_mesh_result_t *mesh = alea_mesh_sample(sys, &cfg);
    ASSERT_NOT_NULL(mesh);

    char path[256];
    FILE *f = alea_tmpfile(path);
    ASSERT_NOT_NULL(f);
    ASSERT(fputs("original", f) >= 0);
    ASSERT_EQ(fclose(f), 0);

    ASSERT_EQ(alea_mesh_export(mesh, (alea_mesh_format_t)99, path), -1);
    f = fopen(path, "r");
    ASSERT_NOT_NULL(f);
    char buf[32] = {0};
    ASSERT(fread(buf, 1, 8, f) == 8);
    fclose(f);
    ASSERT(strncmp(buf, "original", 8) == 0);

    ASSERT_EQ(alea_mesh_export(mesh, ALEA_MESH_VTK, path), 0);
    f = fopen(path, "r");
    ASSERT_NOT_NULL(f);
    ASSERT_NOT_NULL(fgets(buf, sizeof(buf), f));
    fclose(f);
    ASSERT(strstr(buf, "# vtk DataFile") != NULL);
    remove(path);

    alea_mesh_result_free(mesh);
    alea_destroy(sys);
}

TEST_MAIN()
