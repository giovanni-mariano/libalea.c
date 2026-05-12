// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#define ALEA_TEST_IMPLEMENTATION
#include "alea_test.h"
#include "alea.h"
#include "core/alea_system.h"
#include "core/alea_universe.h"

TEST(universe_point_bvh_preserves_all_cells_order) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    ASSERT(mat >= 0);

    for (int i = 0; i < 20; i++) {
        int s = alea_sphere_surface(sys, i + 1, 0.0, 0.0, 0.0, 10.0);
        ASSERT(s >= 0);
        const alea_surface_entry_t* surf = alea_surface_at(sys, s);
        ASSERT_NOT_NULL(surf);
        int c = alea_add_cell(sys, i + 1, surf->neg_node, mat, 1.0, 0);
        ASSERT(c >= 0);
    }

    ASSERT_EQ(alea_build_universe_index(sys), 0);

    alea_cell_hit_t hits[32];
    int n = alea_find_all_cells_at_point_recursive(sys, 0.0, 0.0, 0.0,
                                                   hits, 32);
    ASSERT_EQ(n, 20);
    for (int i = 0; i < n; i++) {
        ASSERT_EQ(hits[i].cell_id, i + 1);
        ASSERT_EQ(hits[i].depth, 0);
    }

    alea_destroy(sys);
}

TEST(universe_point_bvh_filters_by_bbox) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    ASSERT(mat >= 0);

    for (int i = 0; i < 24; i++) {
        double x = (double)i * 10.0;
        int s = alea_sphere_surface(sys, i + 1, x, 0.0, 0.0, 2.0);
        ASSERT(s >= 0);
        const alea_surface_entry_t* surf = alea_surface_at(sys, s);
        ASSERT_NOT_NULL(surf);
        int c = alea_add_cell(sys, i + 1, surf->neg_node, mat, 1.0, 0);
        ASSERT(c >= 0);
    }

    ASSERT_EQ(alea_build_universe_index(sys), 0);

    alea_cell_hit_t hits[8];
    int n = alea_find_all_cells_at_point_recursive(sys, 50.0, 0.0, 0.0,
                                                   hits, 8);
    ASSERT_EQ(n, 1);
    ASSERT_EQ(hits[0].cell_id, 6);
    ASSERT_EQ(hits[0].depth, 0);

    alea_destroy(sys);
}

TEST_MAIN()
