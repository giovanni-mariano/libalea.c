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

    int si2 = alea_sphere_surface(sys, 0, 3, 0, 0, 2.0);
    alea_node_id_t s2 = alea_halfspace(sys, si2, -1);

    /* Universe 0: Container with FILL */
    int sib = alea_box_surface(sys, 0, -10, 10, -10, 10, -10, 10);
    alea_node_id_t box = alea_halfspace(sys, sib, -1);

    /* Register materials */
    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);

    ASSERT(alea_add_cell(sys, 10, s1, m1, -1.0, 1) >= 0);
    ASSERT(alea_add_cell(sys, 11, s2, m2, -1.0, 1) >= 0);

    int cell1_idx = alea_add_cell(sys, 1, box, ALEA_MATERIAL_VOID, 0.0, 0);
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

    /* Universe 0: container with fill */
    int sib = alea_box_surface(sys, 0, -10, 10, -10, 10, -10, 10);
    alea_node_id_t box = alea_halfspace(sys, sib, -1);

    /* Register materials */
    int m1 = alea_add_material(sys, 1);

    ASSERT(alea_add_cell(sys, 10, s, m1, -1.0, 1) >= 0);

    int cell_idx = alea_add_cell(sys, 1, box, ALEA_MATERIAL_VOID, 0.0, 0);
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

TEST(point_query_public_api_uses_deepest_fill_hit) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int child_surface = alea_sphere_surface(sys, 0, 0, 0, 0, 5.0);
    alea_node_id_t child_node = alea_halfspace(sys, child_surface, -1);

    int box_surface = alea_box_surface(sys, 0, -10, 10, -10, 10, -10, 10);
    alea_node_id_t box_node = alea_halfspace(sys, box_surface, -1);

    int mat1 = alea_add_material(sys, 1);
    ASSERT(mat1 >= 0);

    int child_idx = alea_add_cell(sys, 10, child_node, mat1, -1.0, 1);
    ASSERT(child_idx >= 0);

    int parent_idx = alea_add_cell(sys, 1, box_node, ALEA_MATERIAL_VOID, 0.0, 0);
    ASSERT(parent_idx >= 0);
    ASSERT_EQ(alea_set_cell_fill(sys, parent_idx, 1, 0), 0);

    ASSERT_EQ(alea_build_universe_index(sys), 0);

    int cell_id = -1;
    int material = -1;

    ASSERT_EQ(alea_find_cell(sys, 0, 0, 0), child_idx);
    ASSERT_EQ(alea_material_at(sys, 0, 0, 0), 1);
    ASSERT_EQ(alea_find_cell_at(sys, 0, 0, 0, &cell_id, &material), 0);
    ASSERT_EQ(cell_id, 10);
    ASSERT_EQ(material, 1);

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
 * After hierarchy flattening and simplification, Cell B should be removed.
 */
TEST(flatten_removes_empty_cells) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    /* Universe 1: two spheres, one far away */
    int siA = alea_sphere_surface(sys, 0, 0, 0, 0, 1.0);
    alea_node_id_t sA = alea_halfspace(sys, siA, -1);

    int siB = alea_sphere_surface(sys, 0, 100, 0, 0, 1.0);
    alea_node_id_t sB = alea_halfspace(sys, siB, -1);

    /* Universe 0: parent sphere filled with universe 1 */
    int sip = alea_sphere_surface(sys, 0, 0, 0, 0, 3.0);
    alea_node_id_t parent = alea_halfspace(sys, sip, -1);

    /* Register materials */
    int m1 = alea_add_material(sys, 1);
    int m2 = alea_add_material(sys, 2);

    ASSERT(alea_add_cell(sys, 10, sA, m1, -1.0, 1) >= 0);
    ASSERT(alea_add_cell(sys, 11, sB, m2, -2.0, 1) >= 0);

    int cell_idx = alea_add_cell(sys, 1, parent, ALEA_MATERIAL_VOID, 0.0, 0);
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
    alea_simplify_and_prune_cells(sys, &stats);

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

TEST(flatten_preserves_material_indices_and_cell_bookkeeping) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int hydrogen = alea_add_material(sys, 10);
    int uranium = alea_add_material(sys, 20);
    ASSERT_EQ(alea_material_add_nuclide(sys, hydrogen, 1001, "80c", 1), 0);
    ASSERT_EQ(alea_material_add_nuclide(sys, uranium, 92235, "80c", 1), 0);
    for (int i = 0; i < 3; i++) {
        int surface = alea_sphere_surface(sys, i + 1, i * 10, 0, 0, 1);
        ASSERT_EQ(alea_add_cell(sys, 10 + i, alea_halfspace(sys, surface, -1),
                  i == 2 ? ALEA_MATERIAL_VOID : i, -1, 0), i);
    }
    ASSERT_EQ(alea_cell_set_comment(sys, 1, "uranium"), 0);
    for (int copy = 0; copy < 2; copy++) {
        alea_flatten_config_t cfg = ALEA_FLATTEN_DEFAULT;
        cfg.copy_materials = copy;
        alea_system_t* flat = alea_flatten_to_new_system(sys, &cfg);
        ASSERT_NOT_NULL(flat);
        ASSERT_EQ(flat->cells.count, 3);
        for (int i = 0; i < 3; i++) {
            const alea_cell_entry_t* cell = &flat->cells.data[i];
            ASSERT_EQ(alea_find_cell_by_id(flat, cell->mc_cell_id), i);
            ASSERT_EQ(cell->material_id, i == 2 ? 0 : (i + 1) * 10);
            ASSERT_EQ(cell->material_index, copy && i < 2 ? i : -1);
            ASSERT_EQ(cell->original_root_node_id, ALEA_NODE_ID_INVALID);
            ASSERT_EQ(cell->lat_outer_universe, -1);
        }
        ASSERT_STR_EQ(flat->cells.data[1].comments, "uranium");
        if (copy) {
            ASSERT_EQ(alea_material_get_id(flat, flat->cells.data[1].material_index), 20);
            int added = alea_add_material(flat, 0);
            ASSERT(alea_material_get_id(flat, added) > 20);
        }
        ASSERT_EQ(alea_add_cell(flat, 0, flat->cells.data[0].root_node_id,
                               ALEA_MATERIAL_VOID, 0, 0), 3);
        ASSERT(flat->cells.data[3].mc_cell_id > 12);
        ASSERT_EQ(alea_find_cell_by_id(flat, flat->cells.data[3].mc_cell_id), 3);
        alea_destroy(flat);
    }
    ASSERT_EQ(alea_flatten_in_place(sys, NULL), 3);
    ASSERT_EQ(sys->cells.data[1].material_index, uranium);
    ASSERT_EQ(alea_add_cell(sys, 0, sys->cells.data[0].root_node_id,
                           ALEA_MATERIAL_VOID, 0, 0), 3);
    ASSERT(sys->cells.data[3].mc_cell_id > 12);
    alea_destroy(sys);
}

TEST_MAIN()
