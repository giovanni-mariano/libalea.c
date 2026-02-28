// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * test_overlap.c - Verify overlap detection works
 *
 * Creates two overlapping spheres and verifies that alea_find_all_cells
 * returns BOTH cells when querying a point in the overlap region.
 */

#include <stdio.h>
#include <stdlib.h>
#include "alea.h"
#include "core/alea_system.h"
#include "core/alea_universe.h"

#define PASS "\033[32mPASS\033[0m"
#define FAIL "\033[31mFAIL\033[0m"

int main(void) {
    printf("=== Overlap Detection Test ===\n\n");

    /* Create system */
    alea_system_t* sys = alea_create();
    if (!sys) {
        printf("[%s] Failed to create system\n", FAIL);
        return 1;
    }

    /* Create two overlapping spheres using surface API:
     * Sphere 1: center (0,0,0), radius 2
     * Sphere 2: center (1,0,0), radius 2
     *
     * They overlap in the region roughly x in [-1, 2]
     * Point (0.5, 0, 0) is inside BOTH spheres
     */

    /* Sphere 1 - returns stable index */
    int s1 = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 2.0);
    if (s1 < 0) {
        printf("[%s] Failed to create sphere 1 surface\n", FAIL);
        alea_destroy(sys);
        return 1;
    }

    /* Sphere 2 - returns stable index */
    int s2 = alea_sphere_surface(sys, 2, 1.0, 0.0, 0.0, 2.0);
    if (s2 < 0) {
        printf("[%s] Failed to create sphere 2 surface\n", FAIL);
        alea_destroy(sys);
        return 1;
    }

    /* Access via index - always valid */
    alea_node_id_t sphere1 = alea_surface_at(sys, s1)->neg_node;
    alea_node_id_t sphere2 = alea_surface_at(sys, s2)->neg_node;

    /* Register materials */
    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);

    /* Add cells using copied node IDs - both in universe 0 */
    int cell1_idx = alea_add_cell(sys, 1, sphere1, m1, 1.0, 0);
    int cell2_idx = alea_add_cell(sys, 2, sphere2, m2, 2.0, 0);

    if (cell1_idx < 0 || cell2_idx < 0) {
        printf("[%s] Failed to add cells\n", FAIL);
        alea_destroy(sys);
        return 1;
    }

    printf("Created two overlapping spheres:\n");
    printf("  Cell 1: sphere at (0,0,0) r=2, material 1\n");
    printf("  Cell 2: sphere at (1,0,0) r=2, material 2\n\n");

    /* Build universe index */
    if (alea_build_universe_index(sys) < 0) {
        printf("[%s] Failed to build universe index\n", FAIL);
        alea_destroy(sys);
        return 1;
    }

    /* Test 1: Point in overlap region - should find BOTH cells */
    printf("Test 1: Query point (0.5, 0, 0) in overlap region\n");

    alea_cell_hit_t hits[16];
    int hit_count = alea_find_all_cells(sys, 0.5, 0.0, 0.0, hits, 16);

    printf("  Found %d cells\n", hit_count);

    if (hit_count != 2) {
        printf("[%s] Expected 2 cells, got %d\n", FAIL, hit_count);
        printf("  This is the overlap detection bug!\n");
        alea_destroy(sys);
        return 1;
    }

    /* Verify we got both cells */
    int found_cell1 = 0, found_cell2 = 0;
    for (int i = 0; i < hit_count; i++) {
        printf("  Hit %d: cell_id=%d, material=%d\n",
               i, hits[i].cell_id, hits[i].material_id);
        if (hits[i].cell_id == 1) found_cell1 = 1;
        if (hits[i].cell_id == 2) found_cell2 = 1;
    }

    if (!found_cell1 || !found_cell2) {
        printf("[%s] Missing cell in results\n", FAIL);
        alea_destroy(sys);
        return 1;
    }

    printf("[%s] Found both overlapping cells\n\n", PASS);

    /* Test 2: Point only in sphere 1 - should find only cell 1 */
    printf("Test 2: Query point (-1.5, 0, 0) only in sphere 1\n");

    hit_count = alea_find_all_cells(sys, -1.5, 0.0, 0.0, hits, 16);
    printf("  Found %d cells\n", hit_count);

    if (hit_count != 1 || hits[0].cell_id != 1) {
        printf("[%s] Expected only cell 1, got %d cells\n", FAIL, hit_count);
        alea_destroy(sys);
        return 1;
    }
    printf("[%s] Found only cell 1\n\n", PASS);

    /* Test 3: Point only in sphere 2 - should find only cell 2 */
    printf("Test 3: Query point (2.5, 0, 0) only in sphere 2\n");

    hit_count = alea_find_all_cells(sys, 2.5, 0.0, 0.0, hits, 16);
    printf("  Found %d cells\n", hit_count);

    if (hit_count != 1 || hits[0].cell_id != 2) {
        printf("[%s] Expected only cell 2, got %d cells\n", FAIL, hit_count);
        alea_destroy(sys);
        return 1;
    }
    printf("[%s] Found only cell 2\n\n", PASS);

    /* Test 4: Point outside both - should find 0 cells */
    printf("Test 4: Query point (10, 0, 0) outside both\n");

    hit_count = alea_find_all_cells(sys, 10.0, 0.0, 0.0, hits, 16);
    printf("  Found %d cells\n", hit_count);

    if (hit_count != 0) {
        printf("[%s] Expected 0 cells, got %d\n", FAIL, hit_count);
        alea_destroy(sys);
        return 1;
    }
    printf("[%s] Found no cells (void)\n\n", PASS);

    alea_destroy(sys);

    printf("========================================\n");
    printf("All overlap detection tests PASSED\n");
    printf("========================================\n");

    return 0;
}
