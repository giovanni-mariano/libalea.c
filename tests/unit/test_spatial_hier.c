// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#define _POSIX_C_SOURCE 200112L
#define ALEA_TEST_IMPLEMENTATION
#include "alea_test.h"
#include "alea.h"
#include "alea_mcnp.h"
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

static void assert_hits_match(const alea_cell_hit_t* a,
                              const alea_cell_hit_t* b,
                              int count) {
    for (int i = 0; i < count; i++) {
        ASSERT_EQ(a[i].cell_id, b[i].cell_id);
        ASSERT_EQ(a[i].cell_index, b[i].cell_index);
        ASSERT_EQ(a[i].material_id, b[i].material_id);
        ASSERT_EQ(a[i].universe_id, b[i].universe_id);
        ASSERT_EQ(a[i].fill_universe, b[i].fill_universe);
        ASSERT_EQ(a[i].depth, b[i].depth);
    }
}

static int spatial_hits_contain_cell(const alea_spatial_hit_t* hits,
                                     int count,
                                     int cell_id,
                                     int depth) {
    for (int i = 0; i < count; i++) {
        if (hits[i].cell_id == cell_id && (depth < 0 || hits[i].depth == depth)) return 1;
    }
    return 0;
}

TEST(hier_spatial_point_query_matches_recursive_simple) {
    setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    ASSERT(mat >= 0);

    for (int i = 0; i < 8; i++) {
        int s = alea_sphere_surface(sys, i + 1, (double)i * 5.0, 0.0, 0.0, 1.5);
        ASSERT(s >= 0);
        const alea_surface_entry_t* surf = alea_surface_at(sys, s);
        ASSERT_NOT_NULL(surf);
        int c = alea_add_cell(sys, i + 1, surf->neg_node, mat, 1.0, 0);
        ASSERT(c >= 0);
    }

    ASSERT_EQ(alea_hier_spatial_index_build(sys), 0);

    alea_cell_hit_t ref[8];
    alea_cell_hit_t got[8];
    int nr = alea_find_all_cells_at_point_recursive(sys, 10.0, 0.0, 0.0, ref, 8);
    int ng = alea_hier_spatial_find_cells_at_point(sys, 10.0, 0.0, 0.0, got, 8);
    ASSERT_EQ(ng, nr);
    assert_hits_match(ref, got, nr);

    alea_destroy(sys);
    unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_spatial_point_query_matches_recursive_fill) {
    setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    ASSERT(mat >= 0);

    int outer_s = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 10.0);
    int inner_s = alea_sphere_surface(sys, 2, 0.0, 0.0, 0.0, 2.0);
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

    alea_cell_hit_t ref[8];
    alea_cell_hit_t got[8];
    int nr = alea_find_all_cells_at_point_recursive(sys, 0.0, 0.0, 0.0, ref, 8);
    int ng = alea_hier_spatial_find_cells_at_point(sys, 0.0, 0.0, 0.0, got, 8);
    ASSERT_EQ(ng, nr);
    assert_hits_match(ref, got, nr);

    alea_destroy(sys);
    unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_spatial_point_query_matches_recursive_lattice) {
    setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    ASSERT_EQ(alea_hier_spatial_index_build(sys), 0);

    const double points[][3] = {
        {0.0, 0.0, 0.0},
        {0.5, 0.5, 0.0},
        {2.0, 0.0, 0.0},
        {2.5, 0.5, 0.0},
        {4.0, 0.0, 0.0},
        {0.0, 2.0, 0.0}
    };

    for (size_t i = 0; i < sizeof(points) / sizeof(points[0]); i++) {
        alea_cell_hit_t ref[8];
        alea_cell_hit_t got[8];
        int nr = alea_find_all_cells_at_point_recursive(sys,
                                                        points[i][0],
                                                        points[i][1],
                                                        points[i][2],
                                                        ref, 8);
        int ng = alea_hier_spatial_find_cells_at_point(sys,
                                                       points[i][0],
                                                       points[i][1],
                                                       points[i][2],
                                                       got, 8);
        ASSERT_EQ(ng, nr);
        assert_hits_match(ref, got, nr);
    }

    mcnp_model_destroy(model);
    unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_spatial_region_query_matches_flat_simple) {
    setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    ASSERT(mat >= 0);

    for (int i = 0; i < 4; i++) {
        int s = alea_sphere_surface(sys, i + 1, (double)i * 5.0, 0.0, 0.0, 1.0);
        ASSERT(s >= 0);
        const alea_surface_entry_t* surf = alea_surface_at(sys, s);
        ASSERT_NOT_NULL(surf);
        int c = alea_add_cell(sys, i + 1, surf->neg_node, mat, 1.0, 0);
        ASSERT(c >= 0);
    }

    ASSERT_EQ(alea_spatial_index_build(sys), 0);
    ASSERT_EQ(alea_hier_spatial_index_build(sys), 0);

    alea_bbox_t query = {-1.5, 6.5, -2.0, 2.0, -2.0, 2.0};
    alea_spatial_hit_t flat[16];
    alea_spatial_hit_t hier[16];
    int nf = alea_spatial_query_region(sys, &query, flat, 16);
    int nh = alea_hier_spatial_query_region(sys, &query, hier, 16);
    ASSERT_EQ(nh, nf);
    ASSERT(spatial_hits_contain_cell(hier, nh, 1, 0));
    ASSERT(spatial_hits_contain_cell(hier, nh, 2, 0));

    alea_destroy(sys);
    unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_spatial_slice_query_matches_flat_fill) {
    setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    ASSERT(mat >= 0);

    int outer_s = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 10.0);
    int inner_s = alea_sphere_surface(sys, 2, 0.0, 0.0, 0.0, 2.0);
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

    ASSERT_EQ(alea_spatial_index_build(sys), 0);
    ASSERT_EQ(alea_hier_spatial_index_build(sys), 0);

    alea_spatial_hit_t flat[16];
    alea_spatial_hit_t hier[16];
    int nf = alea_spatial_query_slice_z(sys, 0.0, -3.0, 3.0, -3.0, 3.0,
                                        flat, 16);
    int nh = alea_hier_spatial_query_slice_z(sys, 0.0, -3.0, 3.0, -3.0, 3.0,
                                             hier, 16);
    ASSERT_EQ(nh, nf);
    ASSERT(spatial_hits_contain_cell(hier, nh, 2, 1));

    alea_destroy(sys);
    unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_spatial_slice_query_resolves_lattice_terminals) {
    setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    ASSERT_EQ(alea_hier_spatial_index_build(sys), 0);

    alea_spatial_hit_t hits[32];
    int n = alea_hier_spatial_query_slice_z(sys, 0.0, -0.75, 0.75, -0.75, 0.75,
                                            hits, 32);
    ASSERT(n > 0);
    ASSERT(spatial_hits_contain_cell(hits, n, 2, -1));

    n = alea_hier_spatial_query_slice_z(sys, 0.0, 1.25, 2.75, -0.75, 0.75,
                                        hits, 32);
    ASSERT(n > 0);
    ASSERT(spatial_hits_contain_cell(hits, n, 4, -1));

    mcnp_model_destroy(model);
    unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST_MAIN()
