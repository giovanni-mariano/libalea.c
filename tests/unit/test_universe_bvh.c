// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#define _POSIX_C_SOURCE 200112L
#define ALEA_TEST_IMPLEMENTATION
#include "alea_test.h"
#include "alea.h"
#include "core/alea_system.h"
#include "core/alea_universe.h"

TEST(universe_point_bvh_preserves_all_cells_order) {
    setenv("ALEA_UNIVERSE_POINT_BVH_THRESHOLD", "16", 1);
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
    unsetenv("ALEA_UNIVERSE_POINT_BVH_THRESHOLD");
}

TEST(universe_point_bvh_filters_by_bbox) {
    setenv("ALEA_UNIVERSE_POINT_BVH_THRESHOLD", "16", 1);
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
    unsetenv("ALEA_UNIVERSE_POINT_BVH_THRESHOLD");
}

/* In hier mode, alea_prepare_query_acceleration must NOT build per-universe
 * point BVHs — they are unused (BLAS handles all find_all_cells queries). */
TEST(hier_mode_skips_per_universe_bvh_build) {
    setenv("ALEA_UNIVERSE_POINT_BVH_THRESHOLD", "1", 1);

    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    alea_config_t cfg = alea_get_config(sys);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_HIER;
    alea_set_config(sys, &cfg);

    int mat = alea_add_material(sys, 1);
    for (int i = 0; i < 5; i++) {
        int s = alea_sphere_surface(sys, i + 1, (double)i * 10.0, 0.0, 0.0, 2.0);
        alea_add_cell(sys, i + 1, alea_surface_at(sys, s)->neg_node, mat, 1.0, 0);
    }

    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    /* No per-universe BVH should have been built */
    for (size_t i = 0; i < alea_vec_count(&sys->universes); i++) {
        ASSERT_FALSE(sys->universes.data[i].point_bvh_built);
    }

    /* But queries must still work (via BLAS) */
    alea_cell_hit_t hits[8];
    int n = alea_find_all_cells_at_point(sys, 0.0, 0.0, 0.0, hits, 8);
    ASSERT_EQ(n, 1);
    ASSERT_EQ(hits[0].cell_id, 1);

    alea_destroy(sys);
    unsetenv("ALEA_UNIVERSE_POINT_BVH_THRESHOLD");
}

TEST_MAIN()
