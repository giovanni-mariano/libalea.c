// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_test.h"
#include "alea.h"
#include "core/alea_system.h"

static int terminal_material(alea_system_t* sys, double x, double y, double z) {
    alea_cell_hit_t hits[16];
    int n = alea_find_all_cells(sys, x, y, z, hits, 16);
    return n > 0 ? hits[n - 1].material_id : 0;
}

TEST(extract_universe_rebases_requested_universe_to_root) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int surface = alea_sphere_surface(sys, 0, 0, 0, 0, 1.0);
    alea_node_id_t sphere = alea_halfspace(sys, surface, -1);
    int material = alea_add_material(sys, 7);
    ASSERT(alea_add_cell(sys, 40, sphere, material, -1.0, 4) >= 0);

    alea_system_t* extracted = alea_extract_universe(sys, 4);
    ASSERT_NOT_NULL(extracted);
    ASSERT_EQ(terminal_material(extracted, 0, 0, 0), 7);

    alea_cell_info_t info;
    ASSERT_EQ(alea_cell_get_info(extracted, 0, &info), 0);
    ASSERT_EQ(info.universe_id, 0);

    alea_destroy(extracted);
    alea_destroy(sys);
}

TEST(extract_universe_keeps_nested_fills_and_drops_unrelated_cells) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int box_surface = alea_box_surface(sys, 0, -10, 10, -10, 10, -10, 10);
    int sphere_surface = alea_sphere_surface(sys, 0, 0, 0, 0, 1.0);
    alea_node_id_t box = alea_halfspace(sys, box_surface, -1);
    alea_node_id_t sphere = alea_halfspace(sys, sphere_surface, -1);
    int child_material = alea_add_material(sys, 5);
    int unrelated_material = alea_add_material(sys, 9);

    ASSERT(alea_add_cell(sys, 50, sphere, child_material, -1.0, 5) >= 0);
    int container = alea_add_cell(sys, 40, box, ALEA_MATERIAL_VOID, 0.0, 4);
    ASSERT(container >= 0);
    ASSERT_EQ(alea_set_cell_fill(sys, container, 5, 0), 0);
    ASSERT(alea_add_cell(sys, 1, sphere, unrelated_material, -1.0, 0) >= 0);

    alea_system_t* extracted = alea_extract_universe(sys, 4);
    ASSERT_NOT_NULL(extracted);
    ASSERT_EQ(alea_cell_count(extracted), 2);
    ASSERT_EQ(terminal_material(extracted, 0, 0, 0), 5);

    alea_destroy(extracted);
    alea_destroy(sys);
}

TEST_MAIN()
