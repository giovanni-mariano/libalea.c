// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#define _POSIX_C_SOURCE 200112L
#define ALEA_TEST_IMPLEMENTATION
#include "alea_test.h"
#include "alea.h"
#include "core/alea_spatial_hier.h"
#include "core/alea_system.h"

TEST(hier_spatial_builds_universe_blas_stats) {
    setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    ASSERT(mat >= 0);

    for (int i = 0; i < 6; i++) {
        int s = alea_sphere_surface(sys, i + 1, (double)i * 4.0, 0.0, 0.0, 1.0);
        ASSERT(s >= 0);
        const alea_surface_entry_t* surf = alea_surface_at(sys, s);
        ASSERT_NOT_NULL(surf);
        int c = alea_add_cell(sys, i + 1, surf->neg_node, mat, 1.0, 0);
        ASSERT(c >= 0);
    }

    ASSERT_EQ(alea_hier_spatial_index_build(sys), 0);
    ASSERT_NOT_NULL(sys->hier_spatial_index);

    const alea_hier_spatial_stats_t* stats =
        alea_hier_spatial_index_stats(sys->hier_spatial_index);
    ASSERT_NOT_NULL(stats);
    ASSERT_EQ(stats->universe_count, 1);
    ASSERT_EQ(stats->blas_count, 1);
    ASSERT_EQ(stats->linear_universe_count, 0);
    ASSERT_EQ(stats->blas_cell_count, 6);
    ASSERT(stats->blas_node_count > 1);
    ASSERT_EQ(stats->placement_count, 1);
    ASSERT_EQ(stats->root_placement_count, 1);
    ASSERT_EQ(stats->fill_placement_count, 0);
    ASSERT_EQ(stats->lattice_placement_count, 0);
    ASSERT_EQ(stats->max_placement_depth, 0);
    ASSERT_EQ(stats->fill_cell_count, 0);
    ASSERT_EQ(stats->lattice_cell_count, 0);
    ASSERT(stats->memory_bytes > 0);

    alea_destroy(sys);
    unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_spatial_respects_blas_threshold) {
    setenv("ALEA_HIER_BLAS_THRESHOLD", "99", 1);

    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    ASSERT(mat >= 0);

    int s = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 1.0);
    ASSERT(s >= 0);
    const alea_surface_entry_t* surf = alea_surface_at(sys, s);
    ASSERT_NOT_NULL(surf);
    int c = alea_add_cell(sys, 1, surf->neg_node, mat, 1.0, 0);
    ASSERT(c >= 0);

    ASSERT_EQ(alea_hier_spatial_index_build(sys), 0);
    ASSERT_NOT_NULL(sys->hier_spatial_index);

    const alea_hier_spatial_stats_t* stats =
        alea_hier_spatial_index_stats(sys->hier_spatial_index);
    ASSERT_NOT_NULL(stats);
    ASSERT_EQ(stats->universe_count, 1);
    ASSERT_EQ(stats->blas_count, 0);
    ASSERT_EQ(stats->linear_universe_count, 1);
    ASSERT_EQ(stats->blas_cell_count, 0);
    ASSERT_EQ(stats->placement_count, 1);
    ASSERT_EQ(stats->root_placement_count, 1);
    ASSERT_EQ(stats->fill_placement_count, 0);
    ASSERT_EQ(stats->lattice_placement_count, 0);

    alea_destroy(sys);
    unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_spatial_collects_fill_placements_without_terminals) {
    setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    ASSERT(mat >= 0);

    int outer_s = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 10.0);
    int inner_s = alea_sphere_surface(sys, 2, 0.0, 0.0, 0.0, 1.0);
    ASSERT(outer_s >= 0);
    ASSERT(inner_s >= 0);

    const alea_surface_entry_t* outer = alea_surface_at(sys, outer_s);
    const alea_surface_entry_t* inner = alea_surface_at(sys, inner_s);
    ASSERT_NOT_NULL(outer);
    ASSERT_NOT_NULL(inner);

    int container = alea_add_cell(sys, 1, outer->neg_node, mat, 1.0, 0);
    int terminal = alea_add_cell(sys, 2, inner->neg_node, mat, 1.0, 10);
    ASSERT(container >= 0);
    ASSERT(terminal >= 0);
    ASSERT_EQ(alea_set_cell_fill(sys, container, 10, 0), 0);

    ASSERT_EQ(alea_hier_spatial_index_build(sys), 0);

    const alea_hier_spatial_stats_t* stats =
        alea_hier_spatial_index_stats(sys->hier_spatial_index);
    ASSERT_NOT_NULL(stats);
    ASSERT_EQ(stats->placement_count, 2);
    ASSERT_EQ(stats->root_placement_count, 1);
    ASSERT_EQ(stats->fill_placement_count, 1);
    ASSERT_EQ(stats->lattice_placement_count, 0);
    ASSERT_EQ(stats->fill_cell_count, 1);
    ASSERT_EQ(stats->max_placement_depth, 1);

    alea_destroy(sys);
    unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST_MAIN()
