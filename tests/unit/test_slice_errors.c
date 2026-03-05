// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file test_slice_errors.c
 * @brief Unit tests for analytical error line checking (alea_check_slice_errors)
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "alea.h"
#include "alea_slice.h"
#include "core/alea_system.h"
#include "primitives/primitive_create.h"

/* ============================================================================
 * TEST: Two non-overlapping cells — no errors expected
 * ============================================================================ */

static void test_no_errors(void) {
    printf("  test_no_errors... ");

    alea_system_t* sys = alea_create();
    assert(sys != NULL);

    /* Cell 1: Sphere at origin, radius 5 */
    int s1_idx = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    alea_node_id_t sphere1 = alea_surface_at(sys, s1_idx)->neg_node;

    /* Cell 2: Sphere at (15, 0, 0), radius 5 — no overlap */
    int s2_idx = alea_sphere_surface(sys, 2, 15, 0, 0, 5.0);
    alea_node_id_t sphere2 = alea_surface_at(sys, s2_idx)->neg_node;

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);

    alea_add_cell(sys, 1, sphere1, m1, -2.7, 0);
    alea_add_cell(sys, 2, sphere2, m2, -8.0, 0);

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0, -10, 25, -10, 10);

    alea_slice_curves_t* curves = alea_get_slice_curves(sys, &view);
    assert(curves != NULL);

    alea_slice_error_result_t* result = alea_check_slice_errors(
        sys, &view, curves, -1);

    /* Boundary of each sphere: one side has a cell, other side is void.
     * Void is still "no cell" → could be classified as GAP.
     * This is expected: a sphere boundary in open space has gap outside. */
    /* We just verify the function runs without crashing and returns a result */
    assert(result != NULL);

    /* No overlaps should be detected */
    size_t overlap_count = 0;
    for (size_t i = 0; i < result->error_count; i++) {
        if (result->errors[i].type == ALEA_SLICE_ERR_OVERLAP)
            overlap_count++;
    }
    assert(overlap_count == 0);

    printf("OK (errors: %zu, overlaps: %zu)\n",
           result->error_count, overlap_count);

    alea_slice_errors_free(result);
    alea_slice_curves_free(curves);
    alea_destroy(sys);
}

/* ============================================================================
 * TEST: Overlapping cells — OVERLAP errors expected
 * ============================================================================ */

static void test_overlap_detection(void) {
    printf("  test_overlap_detection... ");

    alea_system_t* sys = alea_create();
    assert(sys != NULL);

    /* Cell 1: Sphere at origin, radius 5 */
    int s1_idx = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    alea_node_id_t sphere1 = alea_surface_at(sys, s1_idx)->neg_node;

    /* Cell 2: Sphere at (3, 0, 0), radius 5 — overlaps with sphere1 */
    int s2_idx = alea_sphere_surface(sys, 2, 3, 0, 0, 5.0);
    alea_node_id_t sphere2 = alea_surface_at(sys, s2_idx)->neg_node;

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);

    alea_add_cell(sys, 1, sphere1, m1, -2.7, 0);
    alea_add_cell(sys, 2, sphere2, m2, -8.0, 0);

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0, -10, 13, -10, 10);

    alea_slice_curves_t* curves = alea_get_slice_curves(sys, &view);
    assert(curves != NULL);

    alea_slice_error_result_t* result = alea_check_slice_errors(
        sys, &view, curves, -1);
    assert(result != NULL);

    /* Should detect overlap errors along the curves in the overlap region */
    size_t overlap_count = 0;
    for (size_t i = 0; i < result->error_count; i++) {
        if (result->errors[i].type == ALEA_SLICE_ERR_OVERLAP)
            overlap_count++;
    }
    assert(overlap_count > 0);

    printf("OK (total errors: %zu, overlaps: %zu)\n",
           result->error_count, overlap_count);

    alea_slice_errors_free(result);
    alea_slice_curves_free(curves);
    alea_destroy(sys);
}

/* ============================================================================
 * TEST: Gap in geometry — GAP errors expected
 * ============================================================================ */

static void test_gap_detection(void) {
    printf("  test_gap_detection... ");

    alea_system_t* sys = alea_create();
    assert(sys != NULL);

    /* Single sphere — boundary has void (no cell) on outside */
    int s_idx = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    alea_node_id_t sphere = alea_surface_at(sys, s_idx)->neg_node;

    int m1 = alea_add_material(sys, 1);
    alea_add_cell(sys, 1, sphere, m1, -2.7, 0);

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0, -10, 10, -10, 10);

    alea_slice_curves_t* curves = alea_get_slice_curves(sys, &view);
    assert(curves != NULL);

    alea_slice_error_result_t* result = alea_check_slice_errors(
        sys, &view, curves, -1);
    assert(result != NULL);

    /* The outside of the sphere has no cell → GAP errors */
    size_t gap_count = 0;
    for (size_t i = 0; i < result->error_count; i++) {
        if (result->errors[i].type == ALEA_SLICE_ERR_GAP)
            gap_count++;
    }
    assert(gap_count > 0);

    printf("OK (total errors: %zu, gaps: %zu)\n",
           result->error_count, gap_count);

    alea_slice_errors_free(result);
    alea_slice_curves_free(curves);
    alea_destroy(sys);
}

/* ============================================================================
 * TEST: Two adjacent cells sharing a boundary — no errors
 * ============================================================================ */

static void test_adjacent_cells_no_errors(void) {
    printf("  test_adjacent_cells_no_errors... ");

    alea_system_t* sys = alea_create();
    assert(sys != NULL);

    /* Cell 1: inside sphere, cell 2: outside sphere (complement) */
    int s_idx = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    alea_node_id_t inside  = alea_surface_at(sys, s_idx)->neg_node;
    alea_node_id_t outside = alea_surface_at(sys, s_idx)->pos_node;

    /* Bounding box to contain the "outside" cell */
    int s2_idx = alea_sphere_surface(sys, 2, 0, 0, 0, 15.0);
    alea_node_id_t bound = alea_surface_at(sys, s2_idx)->neg_node;

    alea_node_id_t shell = alea_intersection(sys, outside, bound);

    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);

    alea_add_cell(sys, 1, inside, m1, -2.7, 0);
    alea_add_cell(sys, 2, shell, m2, -8.0, 0);

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, 0, -10, 10, -10, 10);

    alea_slice_curves_t* curves = alea_get_slice_curves(sys, &view);
    assert(curves != NULL);

    alea_slice_error_result_t* result = alea_check_slice_errors(
        sys, &view, curves, -1);
    assert(result != NULL);

    /* The inner sphere boundary should have exactly 1 cell on each side → no errors */
    /* The outer sphere boundary will have GAP outside, so only check inner surface */
    size_t overlap_count = 0;
    for (size_t i = 0; i < result->error_count; i++) {
        if (result->errors[i].type == ALEA_SLICE_ERR_OVERLAP)
            overlap_count++;
    }
    assert(overlap_count == 0);

    /* Check that the inner sphere surface (ID=1) has no errors at all */
    size_t surf1_errors = 0;
    for (size_t i = 0; i < result->error_count; i++) {
        if (result->errors[i].surface_id == 1)
            surf1_errors++;
    }
    assert(surf1_errors == 0);

    printf("OK (total errors: %zu, surf1 errors: %zu, overlaps: %zu)\n",
           result->error_count, surf1_errors, overlap_count);

    alea_slice_errors_free(result);
    alea_slice_curves_free(curves);
    alea_destroy(sys);
}

/* ============================================================================
 * TEST: NULL/empty input handling
 * ============================================================================ */

static void test_null_inputs(void) {
    printf("  test_null_inputs... ");

    /* NULL system */
    alea_slice_error_result_t* r = alea_check_slice_errors(NULL, NULL, NULL, -1);
    assert(r == NULL);

    /* Free NULL is safe */
    alea_slice_errors_free(NULL);

    printf("OK\n");
}

/* ============================================================================
 * MAIN
 * ============================================================================ */

int main(void) {
    printf("\n=== Slice Error Line Tests ===\n\n");

    test_null_inputs();
    test_no_errors();
    test_gap_detection();
    test_overlap_detection();
    test_adjacent_cells_no_errors();

    printf("\n=== All slice error tests passed! ===\n\n");
    return 0;
}
