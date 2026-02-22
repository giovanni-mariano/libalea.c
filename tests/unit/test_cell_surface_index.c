// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file test_cell_surface_index.c
 * @brief Unit tests for cell surface index functionality
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "alea.h"
#include "core/alea_system.h"
#include "core/alea_universe.h"
#include "core/alea_ops.h"
#include "primitives/primitive_create.h"
#include "raycast/raycast.h"

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { \
    printf("  Testing %s... ", #name); \
    if (test_##name()) { \
        printf("PASSED\n"); \
        tests_passed++; \
    } else { \
        printf("FAILED\n"); \
    } \
    tests_run++; \
} while(0)

/* ============================================================================
 * TEST: Empty system
 * ============================================================================ */

static int test_empty_system(void) {
    alea_system_t* sys = alea_create();
    if (!sys) return 0;

    int rc = alea_build_cell_surface_index(sys);
    int ok = (rc == 0);

    alea_destroy(sys);
    return ok;
}

/* ============================================================================
 * TEST: Single cell single surface
 * ============================================================================ */

static int test_single_cell_single_surface(void) {
    alea_system_t* sys = alea_create();
    if (!sys) return 0;

    int s1_idx = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    alea_node_id_t s1 = alea_surface_at(sys, s1_idx)->pos_node;
    alea_add_cell(sys, 1, s1, 0, 0.0, 0);

    int rc = alea_build_cell_surface_index(sys);
    int ok = (rc == 0 &&
              alea_vec_count(&sys->cells) == 1 &&
              sys->cells.data[0].surface_index_count == 1);

    alea_destroy(sys);
    return ok;
}

/* ============================================================================
 * TEST: Single cell multiple surfaces
 * ============================================================================ */

static int test_single_cell_multiple_surfaces(void) {
    alea_system_t* sys = alea_create();
    if (!sys) return 0;

    int s1_idx = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    int s2_idx = alea_sphere_surface(sys, 2, 2, 0, 0, 3.0);
    alea_node_id_t s1 = alea_surface_at(sys, s1_idx)->pos_node;
    alea_node_id_t s2 = alea_surface_at(sys, s2_idx)->pos_node;

    alea_node_id_t region = alea_create_intersection(sys, s1, s2);
    alea_add_cell(sys, 1, region, 0, 0.0, 0);

    int rc = alea_build_cell_surface_index(sys);
    int ok = (rc == 0 &&
              alea_vec_count(&sys->cells) == 1 &&
              sys->cells.data[0].surface_index_count == 2);

    alea_destroy(sys);
    return ok;
}

/* ============================================================================
 * TEST: Multiple cells shared surfaces
 * ============================================================================ */

static int test_multiple_cells_shared_surfaces(void) {
    alea_system_t* sys = alea_create();
    if (!sys) return 0;

    // Three spheres
    int s1_idx = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    int s2_idx = alea_sphere_surface(sys, 2, 10, 0, 0, 3.0);
    int s3_idx = alea_sphere_surface(sys, 3, 0, 10, 0, 3.0);
    alea_node_id_t s1 = alea_surface_at(sys, s1_idx)->pos_node;
    alea_node_id_t s2 = alea_surface_at(sys, s2_idx)->pos_node;
    alea_node_id_t s3 = alea_surface_at(sys, s3_idx)->pos_node;

    // Cell 1: s1 AND s2 -> 2 surfaces
    alea_node_id_t region1 = alea_create_intersection(sys, s1, s2);
    alea_add_cell(sys, 1, region1, 0, 0.0, 0);

    // Cell 2: s1 AND s3 -> 2 surfaces
    alea_node_id_t region2 = alea_create_intersection(sys, s1, s3);
    alea_add_cell(sys, 2, region2, 0, 0.0, 0);

    // Cell 3: just s2 -> 1 surface
    alea_add_cell(sys, 3, s2, 0, 0.0, 0);

    int rc = alea_build_cell_surface_index(sys);

    int ok = (rc == 0 &&
              alea_vec_count(&sys->cells) == 3 &&
              sys->cells.data[0].surface_index_count == 2 &&
              sys->cells.data[1].surface_index_count == 2 &&
              sys->cells.data[2].surface_index_count == 1);

    alea_destroy(sys);
    return ok;
}

/* ============================================================================
 * TEST: Complex CSG tree (union + intersection + complement)
 * ============================================================================ */

static int test_complex_csg_tree(void) {
    alea_system_t* sys = alea_create();
    if (!sys) return 0;

    // Four surfaces
    int s1_idx = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    int s2_idx = alea_sphere_surface(sys, 2, 2, 0, 0, 3.0);
    int s3_idx = alea_sphere_surface(sys, 3, 0, 2, 0, 3.0);
    int s4_idx = alea_sphere_surface(sys, 4, 0, 0, 2, 3.0);
    alea_node_id_t s1 = alea_surface_at(sys, s1_idx)->pos_node;
    alea_node_id_t s2 = alea_surface_at(sys, s2_idx)->pos_node;
    alea_node_id_t s3 = alea_surface_at(sys, s3_idx)->pos_node;
    alea_node_id_t s4 = alea_surface_at(sys, s4_idx)->pos_node;

    // Complex region: ((s1 AND s2) OR s3) AND NOT(s4)
    alea_node_id_t r1 = alea_create_intersection(sys, s1, s2);
    alea_node_id_t r2 = alea_create_union(sys, r1, s3);
    alea_node_id_t r3 = alea_create_complement(sys, s4);
    alea_node_id_t region = alea_create_intersection(sys, r2, r3);

    alea_add_cell(sys, 1, region, 0, 0.0, 0);

    int rc = alea_build_cell_surface_index(sys);

    // Should have all 4 surfaces
    int ok = (rc == 0 &&
              alea_vec_count(&sys->cells) == 1 &&
              sys->cells.data[0].surface_index_count == 4);

    alea_destroy(sys);
    return ok;
}

/* ============================================================================
 * TEST: Rebuild index (call twice)
 * ============================================================================ */

static int test_rebuild_index(void) {
    alea_system_t* sys = alea_create();
    if (!sys) return 0;

    int s1_idx = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    alea_node_id_t s1 = alea_surface_at(sys, s1_idx)->pos_node;
    alea_add_cell(sys, 1, s1, 0, 0.0, 0);

    // Build first time
    int rc1 = alea_build_cell_surface_index(sys);
    size_t count1 = sys->cells.data[0].surface_index_count;

    // Build second time (should free old and rebuild)
    int rc2 = alea_build_cell_surface_index(sys);
    size_t count2 = sys->cells.data[0].surface_index_count;

    int ok = (rc1 == 0 && rc2 == 0 &&
              count1 == 1 && count2 == 1);

    alea_destroy(sys);
    return ok;
}

/* ============================================================================
 * TEST: No duplicate surfaces
 * ============================================================================ */

static int test_no_duplicate_surfaces(void) {
    alea_system_t* sys = alea_create();
    if (!sys) return 0;

    // Create one surface
    int s1_idx = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    alea_node_id_t s1 = alea_surface_at(sys, s1_idx)->pos_node;

    // Use same surface twice in tree: s1 AND s1 (weird but possible)
    alea_node_id_t region = alea_create_intersection(sys, s1, s1);
    alea_add_cell(sys, 1, region, 0, 0.0, 0);

    int rc = alea_build_cell_surface_index(sys);

    // Should have 1 surface, not 2
    int ok = (rc == 0 &&
              sys->cells.data[0].surface_index_count == 1);

    alea_destroy(sys);
    return ok;
}

/* ============================================================================
 * TEST: Cell-aware raycast basic
 * ============================================================================ */

static int test_raycast_cell_aware_basic(void) {
    alea_system_t* sys = alea_create();
    if (!sys) return 0;

    // Create a sphere at origin with radius 5
    // neg_node gives interior (sense=-1)
    int s1_idx = alea_sphere_surface(sys, 1, 0, 0, 0, 5.0);
    alea_node_id_t inside_s1 = alea_surface_at(sys, s1_idx)->neg_node;
    alea_add_cell(sys, 1, inside_s1, 1, 1.0, 0);  // interior of sphere

    // Build universe index first (needed for find_all_cells)
    alea_build_universe_index(sys);

    // Build surface index
    alea_build_cell_surface_index(sys);

    // Cast ray through sphere center
    alea_raycast_result_t result;
    alea_raycast_result_init(&result);

    int rc = alea_raycast_cell_aware(sys, -10, 0, 0, 1, 0, 0, 100, &result);

    // Should have segments: void -> cell -> void
    // Ray enters sphere at t=5, exits at t=15
    int ok = (rc == 0 && result.segment_count >= 1);

    // Find the cell segment
    int found_cell = 0;
    for (size_t i = 0; i < result.segment_count; i++) {
        if (result.segments[i].cell_id == 1) {
            found_cell = 1;
            // Check approximate entry/exit (sphere at origin, r=5)
            // Entry should be around t=5, exit around t=15
            double t_enter = result.segments[i].t_enter;
            double t_exit = result.segments[i].t_exit;
            if (fabs(t_enter - 5.0) > 0.1 || fabs(t_exit - 15.0) > 0.1) {
                printf("(t_enter=%.2f, t_exit=%.2f) ", t_enter, t_exit);
                ok = 0;
            }
        }
    }

    if (!found_cell) {
        printf("(cell segment not found) ");
        ok = 0;
    }

    alea_raycast_result_free(&result);
    alea_destroy(sys);
    return ok;
}

/* ============================================================================
 * TEST: Cell-aware raycast multiple cells
 * ============================================================================ */

static int test_raycast_cell_aware_multiple(void) {
    alea_system_t* sys = alea_create();
    if (!sys) return 0;

    // Two non-overlapping spheres - neg_node gives interior
    int s1_idx = alea_sphere_surface(sys, 1, -10, 0, 0, 3.0);
    int s2_idx = alea_sphere_surface(sys, 2, 10, 0, 0, 3.0);
    alea_node_id_t inside_s1 = alea_surface_at(sys, s1_idx)->neg_node;
    alea_node_id_t inside_s2 = alea_surface_at(sys, s2_idx)->neg_node;

    alea_add_cell(sys, 1, inside_s1, 1, 1.0, 0);
    alea_add_cell(sys, 2, inside_s2, 2, 2.0, 0);

    alea_build_universe_index(sys);
    alea_build_cell_surface_index(sys);

    // Cast ray through both spheres
    alea_raycast_result_t result;
    alea_raycast_result_init(&result);

    int rc = alea_raycast_cell_aware(sys, -20, 0, 0, 1, 0, 0, 40, &result);

    // Should hit cell 1 then cell 2
    int found_cell1 = 0, found_cell2 = 0;
    for (size_t i = 0; i < result.segment_count; i++) {
        if (result.segments[i].cell_id == 1) found_cell1 = 1;
        if (result.segments[i].cell_id == 2) found_cell2 = 1;
    }

    int ok = (rc == 0 && found_cell1 && found_cell2);

    alea_raycast_result_free(&result);
    alea_destroy(sys);
    return ok;
}

/* ============================================================================
 * MAIN
 * ============================================================================ */

int main(void) {
    printf("\n=== Cell Surface Index Unit Tests ===\n\n");

    TEST(empty_system);
    TEST(single_cell_single_surface);
    TEST(single_cell_multiple_surfaces);
    TEST(multiple_cells_shared_surfaces);
    TEST(complex_csg_tree);
    TEST(rebuild_index);
    TEST(no_duplicate_surfaces);
    TEST(raycast_cell_aware_basic);
    TEST(raycast_cell_aware_multiple);

    printf("\n=== Results: %d/%d tests passed ===\n\n",
           tests_passed, tests_run);

    return (tests_passed == tests_run) ? 0 : 1;
}
