// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * test_flatten.c - Universe flattening tests
 */

#include "alea_test.h"
#include "alea.h"
#include "core/alea_system.h"
#include "core/alea_universe.h"
#include "core/alea_simplify.h"

/* Helper: get terminal material at point using hierarchical traversal */
static int get_terminal_material(alea_system_t* sys, double x, double y, double z) {
    alea_cell_hit_t hits[16];
    int count = alea_find_all_cells(sys, x, y, z, hits, 16);
    if (count <= 0) return 0;
    /* Last hit is the innermost (terminal) cell */
    return hits[count - 1].material_id;
}

/*
 * Test geometry: Two-level universe hierarchy
 *
 * Universe 0 (base):
 *   Cell 1: Box [-10,10]^3, FILL=1
 *
 * Universe 1:
 *   Cell 10: Sphere at (-3,0,0), radius 2, material 1
 *   Cell 11: Sphere at (+3,0,0), radius 2, material 2
 *
 * After flattening universe 0:
 *   - Cell 1 (FILL) disappears
 *   - Cells 10 and 11 become terminal cells in universe 0
 */

TEST(flatten_nested_universes) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Universe 1: Two spheres at different positions */
    int si1 = alea_sphere_surface(sys, 0, -3, 0, 0, 2.0);
    alea_node_id_t s1 = alea_halfspace(sys, si1, -1);
    ASSERT(alea_add_cell(sys, 10, s1, 1, -1.0, 1) >= 0);

    int si2 = alea_sphere_surface(sys, 0, 3, 0, 0, 2.0);
    alea_node_id_t s2 = alea_halfspace(sys, si2, -1);
    ASSERT(alea_add_cell(sys, 11, s2, 2, -1.0, 1) >= 0);

    /* Universe 0: Container with FILL */
    int sib = alea_box_surface(sys, 0, -10, 10, -10, 10, -10, 10);
    alea_node_id_t box = alea_halfspace(sys, sib, -1);
    int cell1_idx = alea_add_cell(sys, 1, box, 0, 0.0, 0);
    ASSERT(cell1_idx >= 0);
    ASSERT_EQ(alea_set_cell_fill(sys, cell1_idx, 1, 0), 0);

    /* Build universe index */
    ASSERT_EQ(alea_build_universe_index(sys), 0);

    /* Test point queries before flattening (hierarchical traversal) */
    int mat = get_terminal_material(sys, -3, 0, 0);
    ASSERT_EQ(mat, 1);

    mat = get_terminal_material(sys, 3, 0, 0);
    ASSERT_EQ(mat, 2);

    /* Flatten universe 0 */
    int cell_count = alea_flatten_in_place(sys, 0);
    ASSERT(cell_count > 0);

    /* Rebuild universe index */
    ASSERT_EQ(alea_build_universe_index(sys), 0);

    /* Verify all cells are now in universe 0 with no FILL */
    for (size_t i = 0; i < alea_cell_count(sys); i++) {
        alea_cell_info_t info;
        if (alea_cell_get_info(sys, i, &info) == 0) {
            ASSERT_EQ(info.universe_id, 0);
            ASSERT_EQ(info.fill_universe, 0);
        }
    }

    /* Test point queries after flattening - same results */
    mat = get_terminal_material(sys, -3, 0, 0);
    ASSERT_EQ(mat, 1);

    mat = get_terminal_material(sys, 3, 0, 0);
    ASSERT_EQ(mat, 2);

    alea_destroy(sys);
}

TEST(flatten_simple_fill) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Universe 1: single sphere */
    int si = alea_sphere_surface(sys, 0, 0, 0, 0, 5.0);
    alea_node_id_t s = alea_halfspace(sys, si, -1);
    ASSERT(alea_add_cell(sys, 10, s, 1, -1.0, 1) >= 0);

    /* Universe 0: container with fill */
    int sib = alea_box_surface(sys, 0, -10, 10, -10, 10, -10, 10);
    alea_node_id_t box = alea_halfspace(sys, sib, -1);
    int cell_idx = alea_add_cell(sys, 1, box, 0, 0.0, 0);
    ASSERT(cell_idx >= 0);
    ASSERT_EQ(alea_set_cell_fill(sys, cell_idx, 1, 0), 0);

    ASSERT_EQ(alea_build_universe_index(sys), 0);

    /* Before flatten: material at origin is 1 (hierarchical traversal) */
    ASSERT_EQ(get_terminal_material(sys, 0, 0, 0), 1);

    /* Flatten */
    int count = alea_flatten_in_place(sys, 0);
    ASSERT(count > 0);

    ASSERT_EQ(alea_build_universe_index(sys), 0);

    /* After flatten: same result */
    ASSERT_EQ(get_terminal_material(sys, 0, 0, 0), 1);

    alea_destroy(sys);
}

/*
 * Test: Empty cells from fill universes are removed after flatten + simplify.
 *
 * Universe 0: Sphere at origin, radius 3, FILL=1
 * Universe 1:
 *   Cell A: Sphere at (0,0,0), radius 1, material 1  (inside parent)
 *   Cell B: Sphere at (100,0,0), radius 1, material 2 (outside parent → empty)
 *
 * After flatten + alea_flatten_all_cells, Cell B should be removed.
 */
TEST(flatten_removes_empty_cells) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Universe 1: two spheres, one far away */
    int siA = alea_sphere_surface(sys, 0, 0, 0, 0, 1.0);
    alea_node_id_t sA = alea_halfspace(sys, siA, -1);
    ASSERT(alea_add_cell(sys, 10, sA, 1, -1.0, 1) >= 0);

    int siB = alea_sphere_surface(sys, 0, 100, 0, 0, 1.0);
    alea_node_id_t sB = alea_halfspace(sys, siB, -1);
    ASSERT(alea_add_cell(sys, 11, sB, 2, -2.0, 1) >= 0);

    /* Universe 0: parent sphere filled with universe 1 */
    int sip = alea_sphere_surface(sys, 0, 0, 0, 0, 3.0);
    alea_node_id_t parent = alea_halfspace(sys, sip, -1);
    int cell_idx = alea_add_cell(sys, 1, parent, 0, 0.0, 0);
    ASSERT(cell_idx >= 0);
    ASSERT_EQ(alea_set_cell_fill(sys, cell_idx, 1, 0), 0);

    ASSERT_EQ(alea_build_universe_index(sys), 0);

    /* Flatten */
    int flat_count = alea_flatten_in_place(sys, 0);
    ASSERT(flat_count > 0);

    /* Before simplification: both fill-cells exist */
    size_t before = alea_cell_count(sys);
    ASSERT(before >= 2);

    /* Simplify — should detect and remove the empty cell */
    alea_simplify_stats_t stats = {0};
    alea_flatten_all_cells(sys, &stats);

    size_t after = alea_cell_count(sys);
    ASSERT(stats.empty_cells_removed > 0);
    ASSERT(after < before);

    /* Rebuild universe index and verify surviving cell works */
    ASSERT_EQ(alea_build_universe_index(sys), 0);
    ASSERT_EQ(get_terminal_material(sys, 0, 0, 0), 1);

    /* Point far away should find nothing */
    ASSERT_EQ(get_terminal_material(sys, 100, 0, 0), 0);

    alea_destroy(sys);
}

TEST_MAIN()
