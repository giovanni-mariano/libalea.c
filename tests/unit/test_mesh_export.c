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

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

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

TEST_MAIN()
