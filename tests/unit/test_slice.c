// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file test_slice.c
 * @brief Unit tests for 2D slice API
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "alea.h"
#include "alea_slice.h"
#include "core/alea_system.h"
#include "core/alea_universe.h"
#include "primitives/primitive_create.h"

#define ASSERT_NEAR(a, b, eps) assert(fabs((a) - (b)) < (eps))

/* ============================================================================
 * HELPER FUNCTIONS
 * ============================================================================ */

static alea_system_t* create_test_geometry(void) {
    alea_system_t* sys = alea_create();
    if (!sys) return NULL;

    /* Create a simple geometry: sphere + box */
    /* Using public API which properly registers surfaces */

    /* Cell 1: Sphere at origin, radius 5 */
    int s1_idx = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    alea_node_id_t sphere = alea_surface_at(sys, s1_idx)->neg_node;  /* interior */

    /* Cell 2: Box from (8, -3, -3) to (14, 3, 3) */
    int s2_idx = alea_box_surface(sys, 2, 8, 14, -3, 3, -3, 3);
    alea_node_id_t box = alea_surface_at(sys, s2_idx)->neg_node;  /* interior */

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);

    alea_add_cell(sys, 1, sphere, m1, -2.7, 0);
    alea_add_cell(sys, 2, box, m2, -8.0, 0);

    return sys;
}

/* ============================================================================
 * GRID-BASED TESTS (migrated from old render API)
 * ============================================================================ */

static void test_slice_basic(void) {
    printf("  test_slice_basic... ");

    alea_system_t* sys = create_test_geometry();
    assert(sys != NULL);

    /* Slice at z=0 */
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0, -10, 20, -10, 10);

    int width = 300, height = 200;
    int* cell_ids = malloc(width * height * sizeof(int));
    assert(cell_ids != NULL);

    int rc = alea_find_cells_grid(sys, &view, width, height, -1,
                                       cell_ids, NULL, NULL);
    assert(rc == 0);

    /* Check that we found some cells */
    int cells_found = 0;
    for (int i = 0; i < width * height; i++) {
        if (cell_ids[i] > 0) cells_found++;
    }
    assert(cells_found > 0);

    /* Check pixel at center (should be inside sphere, cell 1)
     * x=0 in world -> u = (0 - (-10)) / 30 * 300 = 100
     * y=0 in world -> v = (0 - (-10)) / 20 * 200 = 100 */
    int center_cell = cell_ids[100 * width + 100];
    assert(center_cell == 1);

    /* Check pixel at right side (should be inside box, cell 2)
     * x=11 in world -> u = (11 - (-10)) / 30 * 300 = 210 */
    int box_x = (int)((11.0 - (-10)) / 30.0 * 300);
    int right_cell = cell_ids[100 * width + box_x];
    assert(right_cell == 2);

    free(cell_ids);
    alea_destroy(sys);

    printf("OK\n");
}

static void test_slice_different_planes(void) {
    printf("  test_slice_different_planes... ");

    alea_system_t* sys = create_test_geometry();
    assert(sys != NULL);

    int width = 200, height = 150;
    int* cell_ids = malloc(width * height * sizeof(int));
    assert(cell_ids != NULL);

    alea_slice_view_t view;
    int rc;

    /* XY plane (z=0) */
    alea_slice_view_axis(&view, 2, 0, -10, 20, -10, 10);
    rc = alea_find_cells_grid(sys, &view, width, height, -1,
                                   cell_ids, NULL, NULL);
    assert(rc == 0);
    int cells_found = 0;
    for (int i = 0; i < width * height; i++) {
        if (cell_ids[i] > 0) cells_found++;
    }
    assert(cells_found > 0);

    /* XZ plane (y=0) */
    alea_slice_view_axis(&view, 1, 0, -10, 20, -10, 10);
    rc = alea_find_cells_grid(sys, &view, width, height, -1,
                                   cell_ids, NULL, NULL);
    assert(rc == 0);
    cells_found = 0;
    for (int i = 0; i < width * height; i++) {
        if (cell_ids[i] > 0) cells_found++;
    }
    assert(cells_found > 0);

    /* YZ plane (x=0) */
    alea_slice_view_axis(&view, 0, 0, -10, 10, -10, 10);
    rc = alea_find_cells_grid(sys, &view, width, height, -1,
                                   cell_ids, NULL, NULL);
    assert(rc == 0);
    cells_found = 0;
    for (int i = 0; i < width * height; i++) {
        if (cell_ids[i] > 0) cells_found++;
    }
    assert(cells_found > 0);

    free(cell_ids);
    alea_destroy(sys);

    printf("OK\n");
}

static void test_slice_arbitrary_plane(void) {
    printf("  test_slice_arbitrary_plane... ");

    alea_system_t* sys = create_test_geometry();
    assert(sys != NULL);

    /* Slice at 45 degree angle through origin */
    alea_slice_view_t view;
    alea_slice_view_init(&view,
                              0, 0, 0,        /* origin */
                              1, 0, 1,        /* normal (will be normalized) */
                              0, 1, 0,        /* up */
                              -15, 15,        /* u bounds */
                              -10, 10);       /* v bounds */

    int width = 200, height = 150;
    int* cell_ids = malloc(width * height * sizeof(int));
    assert(cell_ids != NULL);

    int rc = alea_find_cells_grid(sys, &view, width, height, -1,
                                       cell_ids, NULL, NULL);
    assert(rc == 0);

    /* Should find the sphere (at least) */
    int cells_found = 0;
    for (int i = 0; i < width * height; i++) {
        if (cell_ids[i] > 0) cells_found++;
    }
    assert(cells_found > 0);

    free(cell_ids);
    alea_destroy(sys);

    printf("OK\n");
}

/* ============================================================================
 * ERROR DETECTION TESTS
 * ============================================================================ */

static alea_system_t* create_overlapping_geometry(void) {
    alea_system_t* sys = alea_create();
    if (!sys) return NULL;

    /* Cell 1: Sphere at origin, radius 5 */
    int s1_idx = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    alea_node_id_t sphere1 = alea_surface_at(sys, s1_idx)->neg_node;

    /* Cell 2: Sphere at (3, 0, 0), radius 5 - overlaps with sphere1 */
    int s2_idx = alea_sphere_surface(sys, 2, 3, 0, 0, 5.0);
    alea_node_id_t sphere2 = alea_surface_at(sys, s2_idx)->neg_node;

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);

    alea_add_cell(sys, 1, sphere1, m1, -2.7, 0);
    alea_add_cell(sys, 2, sphere2, m2, -8.0, 0);

    return sys;
}

static void test_slice_overlap_detection(void) {
    printf("  test_slice_overlap_detection... ");

    alea_system_t* sys = create_overlapping_geometry();
    assert(sys != NULL);

    /* Point (1.5, 0, 0) is inside both spheres:
     *   sphere1 at origin r=5: dist = 1.5 < 5  ✓
     *   sphere2 at (3,0,0) r=5: dist = 1.5 < 5 ✓ */
    alea_cell_hit_t hits[8];
    int num_hits = alea_find_all_cells(sys, 1.5, 0, 0, hits, 8);
    assert(num_hits == 2);  /* Both cells claim this point */

    /* Verify both cells are present */
    int found_cell1 = 0, found_cell2 = 0;
    for (int i = 0; i < num_hits; i++) {
        if (hits[i].cell_id == 1) found_cell1 = 1;
        if (hits[i].cell_id == 2) found_cell2 = 1;
    }
    assert(found_cell1 && found_cell2);

    /* Also verify grid detects some overlaps (best-effort with adjacency opt) */
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0, -10, 10, -10, 10);

    int width = 200, height = 200;
    int* cell_ids = malloc(width * height * sizeof(int));
    uint8_t* errors = malloc(width * height * sizeof(uint8_t));
    assert(cell_ids != NULL && errors != NULL);

    int rc = alea_find_cells_grid(sys, &view, width, height, -1,
                                       cell_ids, NULL, errors);
    assert(rc == 0);

    /* Boundary recheck should detect overlaps where the two spheres meet */
    int overlap_count = 0;
    for (int i = 0; i < width * height; i++) {
        if (errors[i] == 1) overlap_count++;
    }
    assert(overlap_count > 0);  /* Boundary recheck must catch offset-sphere overlaps */

    free(cell_ids);
    free(errors);
    alea_destroy(sys);

    printf("OK (point query: 2 cells, grid: %d overlap pixels)\n", overlap_count);
}

static void test_slice_nested_overlap(void) {
    printf("  test_slice_nested_overlap... ");

    /* Two concentric spheres at origin — fully nested overlap */
    alea_system_t* sys = alea_create();
    assert(sys != NULL);

    /* Cell 1: Sphere r=5 at origin */
    int s1_idx = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    alea_node_id_t sphere1 = alea_surface_at(sys, s1_idx)->neg_node;

    /* Cell 2: Sphere r=3 at origin — fully inside cell 1 */
    int s2_idx = alea_sphere_surface(sys, 2, 0, 0, 0, 3.0);
    alea_node_id_t sphere2 = alea_surface_at(sys, s2_idx)->neg_node;

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);

    alea_add_cell(sys, 1, sphere1, m1, -2.7, 0);
    alea_add_cell(sys, 2, sphere2, m2, -8.0, 0);

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0, -10, 10, -10, 10);

    int width = 200, height = 200;
    int* cell_ids = malloc(width * height * sizeof(int));
    uint8_t* errors = calloc(width * height, sizeof(uint8_t));
    assert(cell_ids != NULL && errors != NULL);

    int rc = alea_find_cells_grid(sys, &view, width, height, -1,
                                       cell_ids, NULL, errors);
    assert(rc == 0);

    /* Boundary recheck may miss fully-nested overlaps (no cell-boundary transition
     * between cell 1 and cell 2 if the adjacency opt picks the same cell everywhere) */
    int boundary_overlaps = 0;
    for (int i = 0; i < width * height; i++) {
        if (errors[i] == 1) boundary_overlaps++;
    }

    /* Full overlap check should detect them */
    rc = alea_check_grid_overlaps(sys, &view, width, height, -1,
                                       cell_ids, errors);
    assert(rc == 0);

    int full_overlaps = 0;
    for (int i = 0; i < width * height; i++) {
        if (errors[i] == 1) full_overlaps++;
    }
    assert(full_overlaps > 0);  /* Must detect concentric overlap */
    assert(full_overlaps > boundary_overlaps);  /* Full check finds more */

    free(cell_ids);
    free(errors);
    alea_destroy(sys);

    printf("OK (boundary: %d, full: %d overlap pixels)\n",
           boundary_overlaps, full_overlaps);
}

static void test_slice_undefined_detection(void) {
    printf("  test_slice_undefined_detection... ");

    /* Create a simple geometry with gaps (undefined regions) */
    alea_system_t* sys = alea_create();
    assert(sys != NULL);

    /* Single small sphere - leaves most of the slice undefined */
    int s_idx = alea_sphere_surface(sys, 1, 0, 0, 0, 2.0);
    alea_node_id_t sphere = alea_surface_at(sys, s_idx)->neg_node;

    int m1 = alea_add_material(sys, 1);

    alea_add_cell(sys, 1, sphere, m1, -2.7, 0);

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0, -10, 10, -10, 10);

    int width = 100, height = 100;
    int* cell_ids = malloc(width * height * sizeof(int));
    uint8_t* errors = malloc(width * height * sizeof(uint8_t));
    assert(cell_ids != NULL && errors != NULL);

    int rc = alea_find_cells_grid(sys, &view, width, height, -1,
                                       cell_ids, NULL, errors);
    assert(rc == 0);

    int undefined_count = 0;
    int cells_found = 0;
    for (int i = 0; i < width * height; i++) {
        if (errors[i] == 2) undefined_count++;
        if (cell_ids[i] > 0) cells_found++;
    }

    /* Most pixels should be undefined (outside the sphere) */
    assert(undefined_count > 0);
    assert(undefined_count > cells_found);

    /* Check a point outside the sphere: u=5, v=5
     * pixel_u = (5 - (-10)) / 20 * 100 = 75
     * pixel_v = (5 - (-10)) / 20 * 100 = 75 */
    int px = 75;
    int py = 75;
    assert(errors[py * width + px] == 2);  /* 2 = undefined */

    free(cell_ids);
    free(errors);
    alea_destroy(sys);

    printf("OK (detected %d undefined pixels)\n", undefined_count);
}

static void test_slice_real_boundaries(void) {
    printf("  test_slice_real_boundaries... ");

    alea_system_t* sys = create_test_geometry();
    assert(sys != NULL);

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0, -10, 20, -10, 10);

    int width = 300, height = 200;
    int* cell_ids = malloc(width * height * sizeof(int));
    assert(cell_ids != NULL);

    int rc = alea_find_cells_grid(sys, &view, width, height, -1,
                                       cell_ids, NULL, NULL);
    assert(rc == 0);

    /* Center of slice (pixel 100) = x = -10 + (100.5/300)*30 = 0.05 -> inside sphere */
    int center_cell = cell_ids[100 * width + 100];
    assert(center_cell == 1);

    /* Box center at x=11: pixel = (11 - (-10)) / 30 * 300 = 210 */
    int right_cell = cell_ids[100 * width + 210];
    assert(right_cell == 2);

    /* Far left (pixel 10) = x = -10 + (10.5/300)*30 = -8.95 -> void */
    int left_cell = cell_ids[100 * width + 10];
    assert(left_cell == -1);  /* void */

    free(cell_ids);
    alea_destroy(sys);

    printf("OK\n");
}

/* ============================================================================
 * CURVE GENERATION TESTS
 * ============================================================================ */

static void test_curve_sphere(void) {
    printf("  test_curve_sphere... ");

    alea_system_t* sys = alea_create();
    assert(sys != NULL);

    /* Create a sphere cell at origin */
    int s_idx = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    alea_node_id_t sphere = alea_surface_at(sys, s_idx)->neg_node;

    int m1 = alea_add_material(sys, 1);

    alea_add_cell(sys, 1, sphere, m1, -2.7, 0);

    /* Get curves from Z=0 slice */
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0, -10, 10, -10, 10);
    alea_slice_curves_t* curves = alea_get_slice_curves(sys, &view);
    assert(curves != NULL);

    size_t count = alea_slice_curves_count(curves);
    assert(count > 0);  /* Should have at least the circle from the sphere */

    /* Get the first curve (should be the sphere's circle) */
    alea_curve_t curve;
    int rc = alea_slice_curves_get(curves, 0, &curve);
    assert(rc == 0);

    /* The curve should be a circle */
    assert(curve.type == ALEA_CURVE_CIRCLE);
    ASSERT_NEAR(curve.data.circle.radius, 5.0, 0.01);

    printf("OK (circle r=%.2f)\n", curve.data.circle.radius);

    alea_slice_curves_free(curves);
    alea_destroy(sys);
}

static void test_curve_box(void) {
    printf("  test_curve_box... ");

    alea_system_t* sys = alea_create();
    assert(sys != NULL);

    /* Create a box cell */
    int s_idx = alea_box_surface(sys, 1, -5, 5, -5, 5, -5, 5);
    alea_node_id_t box = alea_surface_at(sys, s_idx)->neg_node;

    int m1 = alea_add_material(sys, 1);

    alea_add_cell(sys, 10, box, m1, -2.7, 0);

    /* Get curves from Z=0 slice */
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0, -10, 10, -10, 10);
    alea_slice_curves_t* curves = alea_get_slice_curves(sys, &view);
    assert(curves != NULL);

    size_t count = alea_slice_curves_count(curves);
    assert(count > 0);  /* Should have the box polygon */

    /* Get the first curve (should be the box's polygon) */
    alea_curve_t curve;
    int rc = alea_slice_curves_get(curves, 0, &curve);
    assert(rc == 0);

    /* The curve should be a polygon (rectangle) */
    assert(curve.type == ALEA_CURVE_POLYGON);

    printf("OK (polygon with %d vertices)\n", curve.data.polygon.count);

    alea_slice_curves_free(curves);
    alea_destroy(sys);
}

static void test_curve_cylinder(void) {
    printf("  test_curve_cylinder... ");

    alea_system_t* sys = alea_create();
    assert(sys != NULL);

    /* Create a Z-axis cylinder */
    int s_idx = alea_cylinder_z_surface(sys, 1, 0, 0, 3.0);
    alea_node_id_t cyl = alea_surface_at(sys, s_idx)->neg_node;

    /* Cap the cylinder with planes */
    int plane_top_idx = alea_plane_surface(sys, 2, 0, 0, 1, -5.0);
    int plane_bot_idx = alea_plane_surface(sys, 3, 0, 0, 1, 5.0);
    alea_node_id_t top = alea_surface_at(sys, plane_top_idx)->neg_node;
    alea_node_id_t bot = alea_surface_at(sys, plane_bot_idx)->pos_node;

    alea_node_id_t capped = alea_intersection(sys, cyl, top);
    capped = alea_intersection(sys, capped, bot);

    int m1 = alea_add_material(sys, 1);

    alea_add_cell(sys, 5, capped, m1, -2.7, 0);

    /* Get curves from Z=0 slice (should produce a circle) */
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0, -10, 10, -10, 10);
    alea_slice_curves_t* curves = alea_get_slice_curves(sys, &view);
    assert(curves != NULL);

    size_t count = alea_slice_curves_count(curves);
    assert(count > 0);

    /* Find the circle curve (cylinder cross-section) */
    int found_circle = 0;
    for (size_t i = 0; i < count; i++) {
        alea_curve_t curve;
        alea_slice_curves_get(curves, i, &curve);
        if (curve.type == ALEA_CURVE_CIRCLE) {
            found_circle = 1;
            ASSERT_NEAR(curve.data.circle.radius, 3.0, 0.01);
            printf("OK (circle r=%.2f)\n", curve.data.circle.radius);
            break;
        }
    }
    assert(found_circle);

    alea_slice_curves_free(curves);
    alea_destroy(sys);
}

static void test_curve_two_cells(void) {
    printf("  test_curve_two_cells... ");

    alea_system_t* sys = alea_create();
    assert(sys != NULL);

    /* Create two adjacent cells sharing a boundary */
    int plane_idx = alea_plane_surface(sys, 1, 1, 0, 0, 0);
    alea_node_id_t left = alea_surface_at(sys, plane_idx)->neg_node;
    alea_node_id_t right = alea_surface_at(sys, plane_idx)->pos_node;

    int box_idx = alea_box_surface(sys, 2, -10, 10, -10, 10, -10, 10);
    alea_node_id_t box = alea_surface_at(sys, box_idx)->neg_node;

    alea_node_id_t cell1_geom = alea_intersection(sys, left, box);
    alea_node_id_t cell2_geom = alea_intersection(sys, right, box);

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);

    alea_add_cell(sys, 1, cell1_geom, m1, -1.0, 0);
    alea_add_cell(sys, 2, cell2_geom, m2, -2.0, 0);

    /* Get curves from Z=0 slice */
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0, -15, 15, -15, 15);
    alea_slice_curves_t* curves = alea_get_slice_curves(sys, &view);
    assert(curves != NULL);

    size_t count = alea_slice_curves_count(curves);
    assert(count > 0);

    /* Find the central dividing line (surface 1, the x=0 plane) */
    int found_interface = 0;
    for (size_t i = 0; i < count; i++) {
        alea_curve_t curve;
        alea_slice_curves_get(curves, i, &curve);
        if (curve.surface_id == 1) {
            found_interface = 1;
            break;
        }
    }

    printf("OK (found %zu curves, interface=%s)\n", count,
           found_interface ? "yes" : "no");

    alea_slice_curves_free(curves);
    alea_destroy(sys);
}

/* ============================================================================
 * MAIN
 * ============================================================================ */

int main(void) {
    printf("\n=== Slice Unit Tests ===\n\n");

    printf("Grid Query Tests:\n");
    test_slice_basic();
    test_slice_different_planes();
    test_slice_arbitrary_plane();

    printf("\nError Detection Tests:\n");
    test_slice_overlap_detection();
    test_slice_nested_overlap();
    test_slice_undefined_detection();
    test_slice_real_boundaries();

    printf("\nCurve Generation Tests:\n");
    test_curve_sphere();
    test_curve_box();
    test_curve_cylinder();
    test_curve_two_cells();

    printf("\n=== All slice tests passed! ===\n\n");
    return 0;
}
