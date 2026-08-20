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
    assert(alea_prepare_query_acceleration(sys) == 0);
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

static int curves_have_surface_type(const alea_slice_curves_t* curves,
                                    int surface_id, alea_curve_type_t type) {
    size_t count = alea_slice_curves_count(curves);
    for (size_t i = 0; i < count; i++) {
        alea_curve_t curve;
        assert(alea_slice_curves_get(curves, i, &curve) == 0);
        if (curve.surface_id == surface_id && curve.type == type) return 1;
    }
    return 0;
}

static int labels_have_surface(const alea_label_position_t* labels,
                               int count, int surface_id) {
    for (int i = 0; i < count; i++) {
        if (labels[i].id == surface_id) return 1;
    }
    return 0;
}

static void assert_surface_can_be_labelled(const alea_slice_curves_t* curves,
                                           int surface_id) {
    alea_label_position_t* labels = NULL;
    int count = 0;
    assert(alea_find_surface_label_positions(
        curves, -10, 10, -10, 10, 160, 160, 2,
        &labels, &count) == 0);
    assert(labels_have_surface(labels, count, surface_id));
    free(labels);
}

static void test_curve_noncanonical_type_identity(void) {
    printf("  test_curve_noncanonical_type_identity... ");

    alea_system_t* sys = alea_create();
    assert(sys != NULL);
    int cone_idx = alea_cone_z_surface(sys, 10, 0, 0, 0, 1.0);
    int torus_idx = alea_torus_z_surface(sys, 20, 0, 0, 0, 4.0, 1.0);
    int material = alea_add_material(sys, 1);
    alea_add_cell(sys, 10, alea_surface_at(sys, cone_idx)->neg_node,
                  material, -1.0, 0);
    alea_add_cell(sys, 20, alea_surface_at(sys, torus_idx)->neg_node,
                  material, -1.0, 0);

    alea_slice_view_t view;
    alea_slice_curves_t* curves;

    /* y=1 cuts x^2 + y^2 = z^2 as a hyperbola. */
    alea_slice_view_axis(&view, 1, 1, -10, 10, -10, 10);
    curves = alea_get_slice_curves(sys, &view);
    assert(curves != NULL);
    assert(curves_have_surface_type(curves, 10, ALEA_CURVE_HYPERBOLA));
    assert_surface_can_be_labelled(curves, 10);
    alea_slice_curves_free(curves);

    /* z=x+1 cuts the same cone as y^2=2x+1, a parabola. */
    alea_slice_view_init(&view, 0, 0, 1, -1, 0, 1, 0, 1, 0,
                         -10, 10, -10, 10);
    curves = alea_get_slice_curves(sys, &view);
    assert(curves != NULL);
    assert(curves_have_surface_type(curves, 10, ALEA_CURVE_PARABOLA));
    assert_surface_can_be_labelled(curves, 10);
    alea_slice_curves_free(curves);

    /* A central Z-torus slice is represented by one quartic curve object. */
    alea_slice_view_axis(&view, 2, 0, -10, 10, -10, 10);
    curves = alea_get_slice_curves(sys, &view);
    assert(curves != NULL);
    assert(curves_have_surface_type(curves, 20, ALEA_CURVE_QUARTIC));
    assert_surface_can_be_labelled(curves, 20);
    alea_slice_curves_free(curves);

    alea_destroy(sys);
    printf("OK\n");
}

static void test_surface_labels_filter_hidden_csg_boundary(void) {
    printf("  test_surface_labels_filter_hidden_csg_boundary... ");

    alea_system_t* sys = alea_create();
    assert(sys != NULL);
    int inner_idx = alea_sphere_surface(sys, 31, 0, 0, 0, 2.0);
    int outer_idx = alea_sphere_surface(sys, 32, 0, 0, 0, 5.0);
    alea_node_id_t inner = alea_surface_at(sys, inner_idx)->neg_node;
    alea_node_id_t outer = alea_surface_at(sys, outer_idx)->neg_node;
    alea_node_id_t redundant_union = alea_union(sys, outer, inner);
    int material = alea_add_material(sys, 1);
    alea_add_cell(sys, 1, redundant_union, material, -1.0, 0);

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0, -10, 10, -10, 10);
    alea_slice_curves_t* curves = alea_get_slice_curves(sys, &view);
    assert(curves != NULL);
    assert(curves_have_surface_type(curves, 31, ALEA_CURVE_CIRCLE));
    assert(curves_have_surface_type(curves, 32, ALEA_CURVE_CIRCLE));

    const int width = 160, height = 160;
    int* cell_ids = malloc((size_t)width * height * sizeof(int));
    assert(cell_ids != NULL);
    assert(alea_find_cells_grid(sys, &view, width, height, -1,
                                cell_ids, NULL, NULL) == 0);

    alea_label_position_t* labels = NULL;
    int count = 0;
    assert(alea_find_surface_label_positions_with_provenance(
        sys, &view, curves, cell_ids, -10, 10, -10, 10, width, height, 2,
        0, &labels, &count) == 0);
    assert(labels_have_surface(labels, count, 32));
    assert(!labels_have_surface(labels, count, 31));

    free(labels);
    free(cell_ids);
    alea_slice_curves_free(curves);
    alea_destroy(sys);
    printf("OK\n");
}

/* A plane split exercises a sparse label candidate that crosses the exact
 * finite RIGHT edge at its midpoint.  It must survive canonical provenance,
 * not only the geometry-only boundary-grid placement path. */
static void test_surface_labels_verify_plane_interface(void) {
    printf("  test_surface_labels_verify_plane_interface... ");
    alea_system_t* sys = alea_create();
    assert(sys != NULL);
    int plane = alea_plane_surface(sys, 71, 1, 0, 0, 0);
    int box = alea_box_surface(sys, 72, -8, 8, -8, 8, -2, 2);
    assert(plane >= 0 && box >= 0);
    alea_node_id_t left = alea_intersection(
        sys, alea_surface_at(sys, plane)->neg_node,
        alea_surface_at(sys, box)->neg_node);
    alea_node_id_t right = alea_intersection(
        sys, alea_surface_at(sys, plane)->pos_node,
        alea_surface_at(sys, box)->neg_node);
    int material = alea_add_material(sys, 1);
    assert(material >= 0);
    assert(alea_add_cell(sys, 1, left, material, -1.0, 0) >= 0);
    assert(alea_add_cell(sys, 2, right, material, -2.0, 0) >= 0);
    assert(alea_prepare_query_acceleration(sys) == 0);

    const int width = 160, height = 160;
    int* ids = malloc((size_t)width * height * sizeof(*ids));
    assert(ids != NULL);
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0, -10, 10, -10, 10);
    assert(alea_find_cells_grid(sys, &view, width, height, -1,
                                ids, NULL, NULL) == 0);
    alea_slice_curves_t* curves = alea_get_slice_curves(sys, &view);
    assert(curves != NULL);
    alea_label_position_t* labels = NULL;
    int count = 0;
    assert(alea_find_surface_label_positions_on_boundaries(
        curves, ids, -10, 10, -10, 10, width, height, 2,
        &labels, &count) == 0);
    assert(labels_have_surface(labels, count, 71));
    free(labels);
    labels = NULL;
    count = 0;
    assert(alea_find_surface_label_positions_with_provenance(
        sys, &view, curves, ids, -10, 10, -10, 10, width, height, 2,
        0, &labels, &count) == 0);
    assert(labels_have_surface(labels, count, 71));
    free(labels);
    alea_slice_curves_free(curves);
    free(ids);
    alea_destroy(sys);
    printf("OK\n");
}

/* Full circles use a fixed geometric anchor.  With a rendered boundary grid,
 * the analytical candidate search must move only a colliding label so two
 * concentric boundaries remain readable. */
static void test_concentric_surface_labels_do_not_stack(void) {
    printf("  test_concentric_surface_labels_do_not_stack... ");

    alea_system_t* sys = alea_create();
    assert(sys != NULL);
    int inner_idx = alea_sphere_surface(sys, 61, 0, 0, 0, 2.0);
    int outer_idx = alea_sphere_surface(sys, 62, 0, 0, 0, 5.0);
    int material = alea_add_material(sys, 1);
    alea_node_id_t inside_inner = alea_surface_at(sys, inner_idx)->neg_node;
    alea_node_id_t shell = alea_intersection(
        sys, alea_surface_at(sys, inner_idx)->pos_node,
        alea_surface_at(sys, outer_idx)->neg_node);
    alea_add_cell(sys, 1, inside_inner, material, -1.0, 0);
    alea_add_cell(sys, 2, shell, material, -1.0, 0);
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0, -10, 10, -10, 10);
    alea_slice_curves_t* curves = alea_get_slice_curves(sys, &view);
    assert(curves != NULL);

    const int width = 200, height = 200;
    int* cell_ids = malloc((size_t)width * height * sizeof(*cell_ids));
    assert(cell_ids != NULL);
    assert(alea_find_cells_grid(sys, &view, width, height, -1,
                                cell_ids, NULL, NULL) == 0);
    alea_label_position_t* labels = NULL;
    int count = 0;
    assert(alea_find_surface_label_positions_on_boundaries(
        curves, cell_ids, -10, 10, -10, 10, width, height, 2,
        &labels, &count) == 0);

    const alea_label_position_t *inner = NULL, *outer = NULL;
    for (int i = 0; i < count; i++) {
        if (labels[i].id == 61) inner = &labels[i];
        if (labels[i].id == 62) outer = &labels[i];
    }
    assert(inner != NULL && outer != NULL);
    double dx = inner->px - outer->px;
    double dy = inner->py - outer->py;
    assert(dx * dx + dy * dy >= 50.0 * 50.0);

    free(labels);
    free(cell_ids);
    alea_slice_curves_free(curves);
    alea_destroy(sys);
    printf("OK\n");
}

static void assert_surface_labelled_on_grid(alea_system_t* sys,
                                            const alea_slice_view_t* view,
                                            int surface_id) {
    const int width = 200, height = 200;
    int* cell_ids = malloc((size_t)width * height * sizeof(int));
    assert(cell_ids != NULL);
    assert(alea_find_cells_grid(sys, view, width, height, -1,
                                cell_ids, NULL, NULL) == 0);
    alea_slice_curves_t* curves = alea_get_slice_curves(sys, view);
    assert(curves != NULL);
    alea_label_position_t* labels = NULL;
    int count = 0;
    assert(alea_find_surface_label_positions_with_provenance(
        sys, view, curves, cell_ids, -10, 10, -10, 10, width, height, 2,
        0, &labels, &count) == 0);
    assert(labels_have_surface(labels, count, surface_id));
    free(labels);
    alea_slice_curves_free(curves);
    free(cell_ids);
}

static void test_noncanonical_surface_labels_on_boundaries(void) {
    printf("  test_noncanonical_surface_labels_on_boundaries... ");

    alea_system_t* cone_sys = alea_create();
    assert(cone_sys != NULL);
    int cone_idx = alea_cone_z_surface(cone_sys, 41, 0, 0, 0, 1.0);
    int material = alea_add_material(cone_sys, 1);
    alea_add_cell(cone_sys, 1, alea_surface_at(cone_sys, cone_idx)->neg_node,
                  material, -1.0, 0);
    alea_slice_view_t view;
    alea_slice_view_init(&view, 0, 0, 1, -1, 0, 1, 0, 1, 0,
                         -10, 10, -10, 10);
    assert_surface_labelled_on_grid(cone_sys, &view, 41);
    alea_slice_view_axis(&view, 1, 1, -10, 10, -10, 10);
    assert_surface_labelled_on_grid(cone_sys, &view, 41);
    alea_destroy(cone_sys);

    alea_system_t* torus_sys = alea_create();
    assert(torus_sys != NULL);
    int torus_idx = alea_torus_z_surface(torus_sys, 42, 0, 0, 0, 4.0, 1.0);
    material = alea_add_material(torus_sys, 1);
    alea_add_cell(torus_sys, 1,
                  alea_surface_at(torus_sys, torus_idx)->neg_node,
                  material, -1.0, 0);
    alea_slice_view_axis(&view, 2, 0, -10, 10, -10, 10);
    assert_surface_labelled_on_grid(torus_sys, &view, 42);
    alea_destroy(torus_sys);

    printf("OK\n");
}

/* A surface that CSG splits into several separated pieces must be labelled on
 * each piece: one label averaged over all of them would land between them, on
 * no boundary at all. */
static void test_boundary_map_labels_each_disconnected_arc(void) {
    printf("  test_boundary_map_labels_each_disconnected_arc... ");

    alea_system_t* sys = alea_create();
    assert(sys != NULL);
    int world_idx = alea_sphere_surface(sys, 50, 0, 0, 0, 9.5);
    int sphere_idx = alea_sphere_surface(sys, 51, 0, 0, 0, 8.0);
    int upper_idx = alea_plane_surface(sys, 52, 0, 1, 0, -2.0);
    int lower_idx = alea_plane_surface(sys, 53, 0, 1, 0, 2.0);
    alea_node_id_t inside = alea_surface_at(sys, sphere_idx)->neg_node;
    int material = alea_add_material(sys, 1);
    /* Two polar caps of one sphere: surface 51 bounds both, as two arcs. */
    alea_node_id_t upper_cap =
        alea_intersection(sys, inside, alea_surface_at(sys, upper_idx)->pos_node);
    alea_node_id_t lower_cap =
        alea_intersection(sys, inside, alea_surface_at(sys, lower_idx)->neg_node);
    alea_add_cell(sys, 1, upper_cap, material, -1.0, 0);
    alea_add_cell(sys, 2, lower_cap, material, -1.0, 0);
    /* The caps sit in a filled world, so their arcs are real transitions
     * between two defined cells rather than coverage gaps. */
    alea_add_cell(sys, 3,
                  alea_intersection(
                      sys, alea_surface_at(sys, world_idx)->neg_node,
                      alea_complement(sys, alea_union(sys, upper_cap, lower_cap))),
                  material, -1.0, 0);
    assert(alea_prepare_query_acceleration(sys) == 0);

    const int width = 200, height = 200;
    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0, -10, 10, -10, 10);
    int* cell_ids = malloc((size_t)width * height * sizeof(int));
    assert(cell_ids != NULL);
    assert(alea_find_cells_grid(sys, &view, width, height, -1,
                                cell_ids, NULL, NULL) == 0);

    alea_slice_surface_boundary_map_t* map = NULL;
    assert(alea_slice_surface_boundary_map_create(
        sys, &view, width, height, cell_ids,
        alea_slice_classify_cell, NULL, &map) == 0);
    assert(map != NULL);

    alea_label_position_t* labels = NULL;
    int count = 0;
    assert(alea_find_surface_labels_on_boundary_map(map, 2, &labels, &count) == 0);

    int arcs = 0, above = 0, below = 0;
    for (int i = 0; i < count; i++) {
        if (labels[i].id != 51) continue;
        arcs++;
        if (labels[i].py > height / 2) above++;
        if (labels[i].py < height / 2) below++;
    }
    assert(arcs == 2);
    assert(above == 1 && below == 1);

    free(labels);
    alea_slice_surface_boundary_map_free(map);
    free(cell_ids);
    alea_destroy(sys);
    printf("OK\n");
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
    test_curve_noncanonical_type_identity();
    test_surface_labels_filter_hidden_csg_boundary();
    test_surface_labels_verify_plane_interface();
    test_concentric_surface_labels_do_not_stack();
    test_noncanonical_surface_labels_on_boundaries();
    test_boundary_map_labels_each_disconnected_arc();

    printf("\n=== All slice tests passed! ===\n\n");
    return 0;
}
