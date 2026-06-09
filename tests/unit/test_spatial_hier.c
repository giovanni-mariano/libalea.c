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
#include "core/alea_transform.h"
#include "raycast/raycast.h"

/* This file intentionally exercises deprecated flat-spatial compatibility APIs
 * as part of the migration gate. */
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

TEST(hier_spatial_builds_universe_blas_stats) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

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
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_spatial_respects_blas_threshold) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "99", 1);

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
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_spatial_collects_fill_placements_without_terminals) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

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
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
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

static int cell_hits_contain_cell(const alea_cell_hit_t* hits,
                                  int count,
                                  int cell_id,
                                  int depth) {
    for (int i = 0; i < count; i++) {
        if (hits[i].cell_id == cell_id && (depth < 0 || hits[i].depth == depth)) return 1;
    }
    return 0;
}

static void assert_raycast_segments_match(const alea_raycast_result_t* flat,
                                          const alea_raycast_result_t* hier) {
    ASSERT_EQ(hier->segments.count, flat->segments.count);
    for (size_t i = 0; i < flat->segments.count; i++) {
        const alea_ray_segment_t* a = &flat->segments.data[i];
        const alea_ray_segment_t* b = &hier->segments.data[i];
        ASSERT_NEAR(b->t_enter, a->t_enter, 1e-9);
        ASSERT_NEAR(b->t_exit, a->t_exit, 1e-9);
        ASSERT_EQ(b->cell_id, a->cell_id);
        ASSERT_EQ(b->material_id, a->material_id);
        ASSERT_NEAR(b->density, a->density, 1e-12);
        ASSERT_EQ(b->enter_surface_id, a->enter_surface_id);
        ASSERT_EQ(b->exit_surface_id, a->exit_surface_id);
    }
}

static void assert_raycast_material_segments_match(const alea_raycast_result_t* flat,
                                                   const alea_raycast_result_t* hier) {
    ASSERT_EQ(hier->segments.count, flat->segments.count);
    for (size_t i = 0; i < flat->segments.count; i++) {
        const alea_ray_segment_t* a = &flat->segments.data[i];
        const alea_ray_segment_t* b = &hier->segments.data[i];
        ASSERT_NEAR(b->t_enter, a->t_enter, 1e-9);
        ASSERT_NEAR(b->t_exit, a->t_exit, 1e-9);
        ASSERT_EQ(b->cell_id, a->cell_id);
        ASSERT_EQ(b->material_id, a->material_id);
        ASSERT_NEAR(b->density, a->density, 1e-12);
    }
}

TEST(hier_spatial_point_query_matches_recursive_simple) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

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
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_spatial_point_query_matches_recursive_fill) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

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
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_spatial_point_query_matches_recursive_lattice) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

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
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_spatial_hex_lattice_points_return_expected_terminals) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    mcnp_model_t* model = mcnp_load("tests/data/mcnp_hex_lattice.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    ASSERT_EQ(alea_hier_spatial_index_build(sys), 0);

    struct {
        double x, y, z;
        int cell_id;
    } cases[] = {
        {0.0, 0.0, 0.0, 1},
        {2.0, 0.0, 0.0, 3},
        {0.5, 0.5, 0.0, 2}
    };

    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        alea_cell_hit_t hits[8];
        int n = alea_hier_spatial_find_cells_at_point(sys,
                                                       cases[i].x,
                                                       cases[i].y,
                                                       cases[i].z,
                                                       hits, 8);
        ASSERT(n > 0);
        ASSERT(cell_hits_contain_cell(hits, n, cases[i].cell_id, -1));
    }
    ASSERT_NULL(sys->spatial_index);

    mcnp_model_destroy(model);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_spatial_region_query_returns_expected_root_cells) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

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

    ASSERT_EQ(alea_hier_spatial_index_build(sys), 0);

    alea_bbox_t query = {-1.5, 6.5, -2.0, 2.0, -2.0, 2.0};
    alea_spatial_hit_t hier[16];
    int nh = alea_hier_spatial_query_region(sys, &query, hier, 16);
    ASSERT_EQ(nh, 2);
    ASSERT(spatial_hits_contain_cell(hier, nh, 1, 0));
    ASSERT(spatial_hits_contain_cell(hier, nh, 2, 0));
    ASSERT(!spatial_hits_contain_cell(hier, nh, 3, 0));
    ASSERT(!spatial_hits_contain_cell(hier, nh, 4, 0));
    ASSERT_NULL(sys->spatial_index);

    alea_destroy(sys);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

static void build_three_sphere_universe(alea_system_t* sys) {
    int mat = alea_add_material(sys, 1);
    ASSERT(mat >= 0);

    for (int i = 0; i < 3; i++) {
        int s = alea_sphere_surface(sys, i + 1, (double)i * 5.0, 0.0, 0.0, 1.0);
        ASSERT(s >= 0);
        const alea_surface_entry_t* surf = alea_surface_at(sys, s);
        ASSERT_NOT_NULL(surf);
        int c = alea_add_cell(sys, i + 1, surf->neg_node, mat, 1.0, 0);
        ASSERT(c >= 0);
    }
}

static void assert_universe_ray_query_three_spheres(const char* threshold) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", threshold, 1);

    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    build_three_sphere_universe(sys);
    ASSERT_EQ(alea_hier_spatial_index_build(sys), 0);

    alea_ray_t ray;
    ASSERT_EQ(alea_ray_init(&ray, -3.0, 0.0, 0.0, 1.0, 0.0, 0.0), 0);

    alea_hier_ray_candidate_t hits[8];
    int n = alea_hier_spatial_query_universe_ray(
        sys, 0,
        ray.ox, ray.oy, ray.oz,
        ray.dx, ray.dy, ray.dz,
        ray.inv_dx, ray.inv_dy, ray.inv_dz,
        0.0, 20.0,
        hits, 8);
    ASSERT_EQ(n, 3);
    ASSERT_EQ(hits[0].cell_id, 1);
    ASSERT_EQ(hits[1].cell_id, 2);
    ASSERT_EQ(hits[2].cell_id, 3);
    ASSERT(hits[0].t_enter <= hits[1].t_enter);
    ASSERT(hits[1].t_enter <= hits[2].t_enter);
    ASSERT_NEAR(hits[0].t_enter, 2.0, 1e-5);
    ASSERT_NEAR(hits[1].t_enter, 7.0, 1e-5);
    ASSERT_NEAR(hits[2].t_enter, 12.0, 1e-5);

    alea_destroy(sys);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_spatial_universe_ray_query_returns_ordered_candidates) {
    assert_universe_ray_query_three_spheres("1");
    assert_universe_ray_query_three_spheres("99");
}

TEST(hier_spatial_ordered_universe_lookup_preserves_cell_order) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    ASSERT(mat >= 0);

    int small_s = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 1.0);
    int large_s = alea_sphere_surface(sys, 2, 0.0, 0.0, 0.0, 10.0);
    ASSERT(small_s >= 0);
    ASSERT(large_s >= 0);

    const alea_surface_entry_t* small = alea_surface_at(sys, small_s);
    const alea_surface_entry_t* large = alea_surface_at(sys, large_s);
    ASSERT_NOT_NULL(small);
    ASSERT_NOT_NULL(large);

    int small_cell = alea_add_cell(sys, 1, small->neg_node, mat, 1.0, 0);
    int large_cell = alea_add_cell(sys, 2, large->neg_node, mat, 1.0, 0);
    ASSERT(small_cell >= 0);
    ASSERT(large_cell >= 0);
    ASSERT_EQ(alea_hier_spatial_index_build(sys), 0);

    int ordered = alea_hier_spatial_find_ordered_cell_in_universe(
        sys, 0, 0.0, 0.0, 0.0);
    ASSERT_EQ(ordered, large_cell);

    alea_destroy(sys);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_spatial_placement_ray_query_returns_root_and_fill) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

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

    alea_ray_t ray;
    ASSERT_EQ(alea_ray_init(&ray, -12.0, 0.0, 0.0, 1.0, 0.0, 0.0), 0);

    alea_hier_placement_ray_candidate_t hits[8];
    int n = alea_hier_spatial_query_placements_ray(
        sys,
        ray.ox, ray.oy, ray.oz,
        ray.dx, ray.dy, ray.dz,
        ray.inv_dx, ray.inv_dy, ray.inv_dz,
        0.0, 24.0,
        hits, 8);
    ASSERT(n >= 2);

    int saw_root = 0;
    int saw_fill = 0;
    for (int i = 0; i < n; i++) {
        ASSERT(hits[i].t_enter <= hits[i].t_exit);
        if (hits[i].universe_id == 0 && hits[i].depth == 0) saw_root = 1;
        if (hits[i].universe_id == 10 && hits[i].depth == 1) saw_fill = 1;
    }
    ASSERT(saw_root);
    ASSERT(saw_fill);

    alea_destroy(sys);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_spatial_slice_query_returns_expected_fill_cell) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

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

    alea_spatial_hit_t hier[16];
    int nh = alea_hier_spatial_query_slice_z(sys, 0.0, -3.0, 3.0, -3.0, 3.0,
                                             hier, 16);
    ASSERT_EQ(nh, 1);
    ASSERT(spatial_hits_contain_cell(hier, nh, 2, 1));
    ASSERT(!spatial_hits_contain_cell(hier, nh, 1, 0));
    ASSERT_NULL(sys->spatial_index);

    alea_destroy(sys);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_spatial_path_query_returns_explicit_fill_chain) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

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

    alea_hier_cell_hit_t hit;
    alea_hier_ray_path_t path;
    ASSERT_EQ(alea_hier_spatial_find_path_at_point(sys, 0.0, 0.0, 0.0,
                                                   &hit, &path), 1);
    ASSERT_EQ(hit.hit.cell_id, 2);
    ASSERT_EQ(path.count, 2);
    ASSERT_EQ(path.entries[0].cell_id, 1);
    ASSERT_EQ(path.entries[1].cell_id, 2);

    alea_matrix_t transform;
    int lattice_cell_index = -2;
    alea_matrix_t lattice_transform;
    ASSERT_EQ(alea_hier_spatial_check_path_containment(
                  sys, &path, (uint32_t)terminal, 0.1, 0.0, 0.0,
                  &transform, &lattice_cell_index, &lattice_transform), 1);
    ASSERT_EQ(lattice_cell_index, -1);
    ASSERT_EQ(alea_hier_spatial_check_path_containment(
                  sys, &path, (uint32_t)terminal, 20.0, 0.0, 0.0,
                  NULL, NULL, NULL), 0);

    alea_hier_cell_hit_t refreshed_hit;
    alea_hier_ray_path_t refreshed_path;
    ASSERT_EQ(alea_hier_spatial_find_path_from_parent(
                  sys, &path, 0, 0.1, 0.0, 0.0,
                  &refreshed_hit, &refreshed_path), 1);
    ASSERT_EQ(refreshed_hit.hit.cell_id, 2);
    ASSERT_EQ(refreshed_path.count, 2);
    ASSERT_EQ(refreshed_path.entries[0].cell_id, 1);
    ASSERT_EQ(refreshed_path.entries[1].cell_id, 2);

    alea_destroy(sys);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_spatial_direct_region_query_resolves_fill) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

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

    alea_bbox_t query = {-3.0, 3.0, -3.0, 3.0, -0.1, 0.1};
    alea_spatial_hit_t hits[16];
    int n = alea_hier_spatial_query_region_direct(sys, &query, hits, 16);
    ASSERT(n > 0);
    ASSERT(spatial_hits_contain_cell(hits, n, 2, 1));

    alea_destroy(sys);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_spatial_translated_fill_returns_expected_terminal) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    ASSERT(mat >= 0);

    int outer_s = alea_sphere_surface(sys, 1, 10.0, 0.0, 0.0, 10.0);
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

    const double translation[3] = {10.0, 0.0, 0.0};
    ASSERT_EQ(alea_add_transform(sys, 99, translation, 3, 0), 0);
    ASSERT_EQ(alea_set_cell_fill(sys, container, 10, 99), 0);

    ASSERT_EQ(alea_hier_spatial_index_build(sys), 0);

    alea_cell_hit_t hits[8];
    int n = alea_hier_spatial_find_cells_at_point(sys, 10.0, 0.0, 0.0,
                                                  hits, 8);
    ASSERT_EQ(n, 2);
    ASSERT_EQ(hits[0].cell_id, 1);
    ASSERT_EQ(hits[0].depth, 0);
    ASSERT_EQ(hits[1].cell_id, 2);
    ASSERT_EQ(hits[1].depth, 1);

    n = alea_hier_spatial_find_cells_at_point(sys, -5.0, 0.0, 0.0, hits, 8);
    ASSERT(n <= 0 || !cell_hits_contain_cell(hits, n, 2, 1));
    ASSERT_NULL(sys->spatial_index);

    alea_destroy(sys);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_spatial_chain_region_query_carries_fill_ancestor) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

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

    alea_bbox_t query = {-3.0, 3.0, -3.0, 3.0, -0.1, 0.1};
    alea_hier_spatial_chain_hit_t hits[16];
    int n = alea_hier_spatial_query_region_chain(sys, &query, hits, 16);
    ASSERT(n > 0);

    int found = 0;
    for (int i = 0; i < n; i++) {
        if (hits[i].hit.cell_id == 2 && hits[i].hit.depth == 1) {
            ASSERT_EQ(hits[i].ancestor_count, 1);
            ASSERT_EQ(hits[i].ancestor_cell_indices[0], (uint32_t)container);
            ASSERT_EQ(hits[i].chain_truncated, 0);
            found = 1;
        }
    }
    ASSERT(found);

    alea_destroy(sys);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_spatial_universe_region_query_filters_to_one_universe) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

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

    alea_bbox_t query = {-3.0, 3.0, -3.0, 3.0, -0.1, 0.1};
    alea_spatial_hit_t hits[16];
    int n = alea_hier_spatial_query_universe_region(sys, 10, &query, hits, 16);
    ASSERT(n > 0);
    ASSERT(spatial_hits_contain_cell(hits, n, 2, 0));
    ASSERT(!spatial_hits_contain_cell(hits, n, 1, 0));
    for (int i = 0; i < n; i++) {
        ASSERT_EQ(hits[i].universe_id, 10);
    }

    alea_destroy(sys);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_spatial_slice_query_resolves_lattice_terminals) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

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
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_raycast_matches_flat_root_universe) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

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

    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    alea_raycast_result_t flat;
    alea_raycast_result_t hier;
    alea_raycast_result_init(&flat);
    alea_raycast_result_init(&hier);

    ASSERT_EQ(alea_raycast(sys, -5.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                           10.0, &flat), 0);
    ASSERT_EQ(alea_raycast_hier(sys, -5.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                                10.0, &hier), 0);

    assert_raycast_segments_match(&flat, &hier);

    alea_raycast_result_free(&flat);
    alea_raycast_result_free(&hier);
    alea_destroy(sys);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_blas_raycast_matches_hier_root_universe) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    ASSERT(mat >= 0);

    for (int i = 0; i < 3; i++) {
        int s = alea_sphere_surface(sys, i + 1, (double)i * 4.0, 0.0, 0.0, 1.0);
        ASSERT(s >= 0);
        const alea_surface_entry_t* surf = alea_surface_at(sys, s);
        ASSERT_NOT_NULL(surf);
        int c = alea_add_cell(sys, i + 1, surf->neg_node, mat, 1.0, 0);
        ASSERT(c >= 0);
    }

    alea_raycast_result_t hier;
    alea_raycast_result_t blas;
    alea_raycast_result_init(&hier);
    alea_raycast_result_init(&blas);

    ASSERT_EQ(alea_raycast_hier_blas_experimental(
                  sys, -3.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                  14.0, &blas), 0);
    ASSERT_NULL(sys->spatial_index);
    ASSERT_NOT_NULL(sys->hier_spatial_index);

    ASSERT_EQ(alea_raycast_hier(sys, -3.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                                14.0, &hier), 0);

    assert_raycast_material_segments_match(&hier, &blas);

    alea_raycast_result_free(&hier);
    alea_raycast_result_free(&blas);
    alea_destroy(sys);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_blas_raycast_matches_hier_fill_universe) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    mcnp_model_t* model = mcnp_load("tests/data/simple_fill.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    alea_raycast_result_t hier;
    alea_raycast_result_t blas;
    alea_raycast_result_init(&hier);
    alea_raycast_result_init(&blas);

    ASSERT_EQ(alea_raycast_hier_blas_experimental(
                  sys, -10.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                  25.0, &blas), 0);
    ASSERT_NULL(sys->spatial_index);
    ASSERT_NOT_NULL(sys->hier_spatial_index);

    ASSERT_EQ(alea_raycast_hier(sys, -10.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                                25.0, &hier), 0);

    assert_raycast_material_segments_match(&hier, &blas);

    alea_raycast_result_free(&hier);
    alea_raycast_result_free(&blas);
    mcnp_model_destroy(model);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_blas_raycast_matches_hier_lattice_universe) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    alea_raycast_result_t hier;
    alea_raycast_result_t blas;
    alea_raycast_result_init(&hier);
    alea_raycast_result_init(&blas);

    ASSERT_EQ(alea_raycast_hier_blas_experimental(
                  sys, -1.5, 0.0, 0.0, 1.0, 0.0, 0.0,
                  7.0, &blas), 0);
    ASSERT_NULL(sys->spatial_index);
    ASSERT_NOT_NULL(sys->hier_spatial_index);

    ASSERT_EQ(alea_raycast_hier(sys, -1.5, 0.0, 0.0, 1.0, 0.0, 0.0,
                                7.0, &hier), 0);

    assert_raycast_material_segments_match(&hier, &blas);

    alea_raycast_result_free(&hier);
    alea_raycast_result_free(&blas);
    mcnp_model_destroy(model);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_raycast_does_not_build_flat_spatial_index) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    mcnp_model_t* model = mcnp_load("tests/data/simple_fill.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    ASSERT_NULL(sys->spatial_index);

    alea_raycast_result_t result;
    alea_raycast_result_init(&result);
    ASSERT_EQ(alea_raycast_hier(sys, -5.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                                10.0, &result), 0);

    ASSERT_NULL(sys->spatial_index);
    ASSERT_NOT_NULL(sys->hier_spatial_index);
    ASSERT(result.segments.count > 0);

    alea_raycast_result_free(&result);
    mcnp_model_destroy(model);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_fast_segments_raycast_does_not_build_flat_spatial_index) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    ASSERT_NULL(sys->spatial_index);

    alea_raycast_result_t flat;
    alea_raycast_result_t hier;
    alea_raycast_result_init(&flat);
    alea_raycast_result_init(&hier);

    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    ASSERT_EQ(alea_raycast(sys, -1.5, 0.0, 0.0, 1.0, 0.0, 0.0,
                           7.0, &flat), 0);

    if (sys->spatial_index) {
        alea_spatial_index_free(sys->spatial_index);
        sys->spatial_index = NULL;
    }

    ASSERT_EQ(alea_raycast_hier_fast_segments(sys, -1.5, 0.0, 0.0,
                                              1.0, 0.0, 0.0,
                                              7.0, &hier), 0);

    ASSERT_NULL(sys->spatial_index);
    ASSERT_NOT_NULL(sys->hier_spatial_index);
    assert_raycast_material_segments_match(&flat, &hier);

    alea_raycast_result_free(&flat);
    alea_raycast_result_free(&hier);
    mcnp_model_destroy(model);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(query_acceleration_respects_hier_spatial_mode) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    mcnp_model_t* model = mcnp_load("tests/data/simple_fill.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    alea_config_t cfg = alea_get_config(sys);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_HIER;
    alea_set_config(sys, &cfg);

    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    ASSERT_NULL(sys->spatial_index);
    ASSERT_NOT_NULL(sys->hier_spatial_index);
    ASSERT(alea_system_query_cache_ready(sys, ALEA_CACHE_RAYCAST_HIER));

    mcnp_model_destroy(model);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(public_point_query_hier_mode_builds_no_flat_spatial_index) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    ASSERT(mat >= 0);
    int s = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 5.0);
    ASSERT(s >= 0);
    const alea_surface_entry_t* surf = alea_surface_at(sys, s);
    ASSERT_NOT_NULL(surf);
    int c = alea_add_cell(sys, 1, surf->neg_node, mat, 1.0, 0);
    ASSERT(c >= 0);

    alea_config_t cfg = alea_get_config(sys);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_HIER;
    alea_set_config(sys, &cfg);

    ASSERT_NULL(sys->spatial_index);
    ASSERT_NULL(sys->hier_spatial_index);

    alea_cell_hit_t hits[8];
    int n = alea_find_all_cells(sys, 0.0, 0.0, 0.0, hits, 8);
    ASSERT(n > 0);
    ASSERT_NULL(sys->spatial_index);
    ASSERT_NOT_NULL(sys->hier_spatial_index);

    alea_destroy(sys);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(query_acceleration_stats_reports_hier_index_shape) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    ASSERT(mat >= 0);
    for (int i = 0; i < 3; i++) {
        int s = alea_sphere_surface(sys, i + 1, (double)i * 3.0, 0.0, 0.0, 1.0);
        ASSERT(s >= 0);
        const alea_surface_entry_t* surf = alea_surface_at(sys, s);
        ASSERT_NOT_NULL(surf);
        int c = alea_add_cell(sys, i + 1, surf->neg_node, mat, 1.0, 0);
        ASSERT(c >= 0);
    }

    alea_config_t cfg = alea_get_config(sys);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_HIER;
    alea_set_config(sys, &cfg);

    alea_query_acceleration_stats_t stats;
    ASSERT_EQ(alea_query_acceleration_stats(sys, &stats), 0);
    ASSERT_EQ(stats.configured_mode, ALEA_SPATIAL_MODE_HIER);
    ASSERT_EQ(stats.resolved_mode, ALEA_SPATIAL_MODE_HIER);
    ASSERT(!stats.built);

    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    ASSERT_EQ(alea_query_acceleration_stats(sys, &stats), 0);
    ASSERT(stats.built);
    ASSERT_EQ(stats.flat_instance_count, 0);
    ASSERT_EQ(stats.hier_universe_count, 1);
    ASSERT_EQ(stats.hier_blas_count, 1);
    ASSERT_EQ(stats.hier_blas_cell_count, 3);
    ASSERT_EQ(stats.hier_placement_count, 1);
    ASSERT(stats.memory_bytes > 0);
    ASSERT_NULL(sys->spatial_index);
    ASSERT_NOT_NULL(sys->hier_spatial_index);

    alea_destroy(sys);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(query_acceleration_stats_reports_flat_index_shape) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    ASSERT(mat >= 0);
    int s = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 2.0);
    ASSERT(s >= 0);
    const alea_surface_entry_t* surf = alea_surface_at(sys, s);
    ASSERT_NOT_NULL(surf);
    int c = alea_add_cell(sys, 1, surf->neg_node, mat, 1.0, 0);
    ASSERT(c >= 0);

    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);

    alea_query_acceleration_stats_t stats;
    ASSERT_EQ(alea_query_acceleration_stats(sys, &stats), 0);
    ASSERT_EQ(stats.configured_mode, ALEA_SPATIAL_MODE_FLAT);
    ASSERT_EQ(stats.resolved_mode, ALEA_SPATIAL_MODE_FLAT);
    ASSERT(stats.built);
    ASSERT_EQ(stats.flat_instance_count, 1);
    ASSERT_EQ(stats.hier_blas_count, 0);
    ASSERT(stats.memory_bytes > 0);
    ASSERT_NOT_NULL(sys->spatial_index);
    ASSERT_NULL(sys->hier_spatial_index);

    alea_destroy(sys);
}

TEST(flat_path_api_wraps_flat_instance_volumes) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    ASSERT(mat >= 0);
    int s = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 2.0);
    ASSERT(s >= 0);
    const alea_surface_entry_t* surf = alea_surface_at(sys, s);
    ASSERT_NOT_NULL(surf);
    int c = alea_add_cell(sys, 1, surf->neg_node, mat, 1.0, 0);
    ASSERT(c >= 0);

    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    ASSERT_NOT_NULL(sys->spatial_index);

    size_t flat_instances = alea_spatial_index_instance_count(sys);
    ASSERT_EQ(flat_instances, 1);
    ASSERT_EQ(alea_volume_path_count(sys), flat_instances);

    alea_volume_path_t path;
    ASSERT_EQ(alea_volume_paths_get(sys, &path, 1), 1);
    ASSERT_EQ(path.path_id, 0);
    ASSERT_EQ(path.terminal_cell_index, c);

    double volumes[1] = {0.0};
    ASSERT_EQ(alea_estimate_path_volumes(sys, 8, volumes, NULL), 0);
    ASSERT(volumes[0] >= 0.0);
    ASSERT_NOT_NULL(sys->spatial_index);
    ASSERT_NULL(sys->hier_spatial_index);

    alea_destroy(sys);
}

TEST(volume_path_index_reuses_and_invalidates) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    ASSERT(mat >= 0);
    int s = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 2.0);
    ASSERT(s >= 0);
    const alea_surface_entry_t* surf = alea_surface_at(sys, s);
    ASSERT_NOT_NULL(surf);

    int c0 = alea_add_cell(sys, 1, surf->neg_node, mat, 1.0, 0);
    ASSERT(c0 >= 0);

    ASSERT_EQ(alea_volume_path_count(sys), 1);
    void* first_index = sys->volume_path_index;
    ASSERT_NOT_NULL(first_index);
    ASSERT_EQ(alea_volume_path_count(sys), 1);
    ASSERT_EQ((void*)sys->volume_path_index, first_index);

    int c1 = alea_add_cell(sys, 2, surf->neg_node, mat, 1.0, 0);
    ASSERT(c1 >= 0);
    ASSERT_NULL(sys->volume_path_index);
    ASSERT_EQ(alea_volume_path_count(sys), 2);
    ASSERT_NOT_NULL(sys->volume_path_index);

    alea_volume_path_t paths[2];
    ASSERT_EQ(alea_volume_paths_get(sys, paths, 2), 2);
    ASSERT_EQ(paths[0].path_id, 0);
    ASSERT_EQ(paths[1].path_id, 1);

    alea_destroy(sys);
}

TEST(volume_path_index_respects_max_count_guard) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);
    alea_setenv("ALEA_VOLUME_PATH_MAX_COUNT", "1", 1);

    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    ASSERT(mat >= 0);
    for (int i = 0; i < 2; i++) {
        int s = alea_sphere_surface(sys, i + 1, (double)i * 5.0, 0.0, 0.0, 1.0);
        ASSERT(s >= 0);
        const alea_surface_entry_t* surf = alea_surface_at(sys, s);
        ASSERT_NOT_NULL(surf);
        int c = alea_add_cell(sys, i + 1, surf->neg_node, mat, 1.0, 0);
        ASSERT(c >= 0);
    }

    alea_config_t cfg = alea_get_config(sys);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_HIER;
    alea_set_config(sys, &cfg);

    alea_error_clear();
    ASSERT_EQ(alea_volume_path_count(sys), 0);
    ASSERT_EQ(alea_error_code(), (int)ALEA_ERR_INVALID_STATE);
    ASSERT_NULL(sys->volume_path_index);

    double volumes[2] = {0.0, 0.0};
    ASSERT_EQ(alea_estimate_path_volumes(sys, 8, volumes, NULL), -1);
    ASSERT_EQ(alea_error_code(), (int)ALEA_ERR_INVALID_STATE);
    ASSERT_NULL(sys->spatial_index);

    alea_destroy(sys);
    alea_unsetenv("ALEA_VOLUME_PATH_MAX_COUNT");
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(volume_path_index_rejects_duplicate_structural_keys) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

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

    alea_config_t cfg = alea_get_config(sys);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_HIER;
    alea_set_config(sys, &cfg);

    ASSERT_EQ(alea_build_universe_index(sys), 0);
    alea_universe_t* root = NULL;
    for (size_t i = 0; i < sys->universes.count; i++) {
        if (sys->universes.data[i].universe_id == 0) {
            root = &sys->universes.data[i];
            break;
        }
    }
    ASSERT_NOT_NULL(root);
    ASSERT_EQ(alea_vec_push(&root->cell_indices, (size_t)c, size_t), 0);

    alea_error_clear();
    ASSERT_EQ(alea_volume_path_count(sys), 0);
    ASSERT_EQ(alea_error_code(), (int)ALEA_ERR_INVALID_STATE);
    ASSERT_NULL(sys->volume_path_index);

    alea_destroy(sys);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_volume_paths_distinguish_repeated_transformed_fills) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    ASSERT(mat >= 0);

    int left_s = alea_sphere_surface(sys, 1, -5.0, 0.0, 0.0, 3.0);
    int right_s = alea_sphere_surface(sys, 2, 5.0, 0.0, 0.0, 3.0);
    int child_s = alea_sphere_surface(sys, 3, 0.0, 0.0, 0.0, 1.0);
    ASSERT(left_s >= 0);
    ASSERT(right_s >= 0);
    ASSERT(child_s >= 0);

    int left = alea_add_cell(sys, 10, alea_surface_at(sys, left_s)->neg_node,
                             mat, 1.0, 0);
    int right = alea_add_cell(sys, 20, alea_surface_at(sys, right_s)->neg_node,
                              mat, 1.0, 0);
    int child = alea_add_cell(sys, 30, alea_surface_at(sys, child_s)->neg_node,
                              mat, 1.0, 10);
    ASSERT(left >= 0);
    ASSERT(right >= 0);
    ASSERT(child >= 0);

    const double left_tr[3] = {-5.0, 0.0, 0.0};
    const double right_tr[3] = {5.0, 0.0, 0.0};
    ASSERT_EQ(alea_add_transform(sys, 101, left_tr, 3, 0), 0);
    ASSERT_EQ(alea_add_transform(sys, 102, right_tr, 3, 0), 0);
    ASSERT_EQ(alea_set_fill(sys, left, 10, 101), 0);
    ASSERT_EQ(alea_set_fill(sys, right, 10, 102), 0);

    alea_config_t cfg = alea_get_config(sys);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_HIER;
    alea_set_config(sys, &cfg);

    size_t count = alea_volume_path_count(sys);
    ASSERT_EQ(count, 2);
    ASSERT_NULL(sys->spatial_index);
    ASSERT_NOT_NULL(sys->hier_spatial_index);

    alea_volume_path_t paths[2];
    ASSERT_EQ(alea_volume_paths_get(sys, paths, 2), 2);
    ASSERT_EQ(paths[0].terminal_cell_index, child);
    ASSERT_EQ(paths[1].terminal_cell_index, child);
    ASSERT_EQ(paths[0].ancestor_count, 1);
    ASSERT_EQ(paths[1].ancestor_count, 1);
    ASSERT(paths[0].ancestor_cell_indices[0] != paths[1].ancestor_cell_indices[0]);

    alea_volume_path_t hit;
    ASSERT_EQ(alea_volume_path_at_point(sys, -5.0, 0.0, 0.0, &hit), 1);
    ASSERT_EQ(hit.terminal_cell_index, child);
    ASSERT_EQ(hit.ancestor_cell_indices[0], left);
    ASSERT_EQ(alea_volume_path_at_point(sys, 5.0, 0.0, 0.0, &hit), 1);
    ASSERT_EQ(hit.terminal_cell_index, child);
    ASSERT_EQ(hit.ancestor_cell_indices[0], right);
    ASSERT_NULL(sys->spatial_index);

    alea_destroy(sys);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_volume_paths_distinguish_rect_lattice_elements) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    alea_config_t cfg = alea_get_config(sys);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_HIER;
    alea_set_config(sys, &cfg);

    size_t count = alea_volume_path_count(sys);
    ASSERT_EQ(count, 18);
    ASSERT_NULL(sys->spatial_index);
    ASSERT_NOT_NULL(sys->hier_spatial_index);

    alea_volume_path_t paths[18];
    ASSERT_EQ(alea_volume_paths_get(sys, paths, 18), 18);

    int fuel_u1_paths = 0;
    int fuel_u3_paths = 0;
    int fuel_u1_linear_seen[9] = {0};
    int fuel_u3_linear_seen[9] = {0};
    for (size_t i = 0; i < count; i++) {
        ASSERT_EQ(paths[i].lattice_step_count, 1);
        int linear = paths[i].lattice_steps[0].linear_index;
        ASSERT(linear >= 0 && linear < 9);
        if (paths[i].terminal_cell_id == 1) {
            fuel_u1_paths++;
            fuel_u1_linear_seen[linear]++;
        } else if (paths[i].terminal_cell_id == 3) {
            fuel_u3_paths++;
            fuel_u3_linear_seen[linear]++;
        }
    }
    ASSERT_EQ(fuel_u1_paths, 5);
    ASSERT_EQ(fuel_u3_paths, 4);
    ASSERT_EQ(fuel_u1_linear_seen[0], 1);
    ASSERT_EQ(fuel_u1_linear_seen[2], 1);
    ASSERT_EQ(fuel_u1_linear_seen[4], 1);
    ASSERT_EQ(fuel_u1_linear_seen[6], 1);
    ASSERT_EQ(fuel_u1_linear_seen[8], 1);
    ASSERT_EQ(fuel_u3_linear_seen[1], 1);
    ASSERT_EQ(fuel_u3_linear_seen[3], 1);
    ASSERT_EQ(fuel_u3_linear_seen[5], 1);
    ASSERT_EQ(fuel_u3_linear_seen[7], 1);

    double volumes[18] = {0.0};
    ASSERT_EQ(alea_estimate_path_volumes(sys, 4096, volumes, NULL), 0);

    int nonzero_u1_moderator = 0;
    int nonzero_u3_moderator = 0;
    for (size_t i = 0; i < count; i++) {
        if (paths[i].terminal_cell_id == 2 && volumes[i] > 0.0) nonzero_u1_moderator++;
        if (paths[i].terminal_cell_id == 4 && volumes[i] > 0.0) nonzero_u3_moderator++;
    }
    ASSERT(nonzero_u1_moderator >= 2);
    ASSERT(nonzero_u3_moderator >= 2);
    ASSERT_NULL(sys->spatial_index);

    mcnp_model_destroy(model);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_path_volume_estimation_no_flat_spatial_index) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    ASSERT(mat >= 0);

    int outer_s = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 3.0);
    int inner_s = alea_sphere_surface(sys, 2, 0.0, 0.0, 0.0, 1.0);
    ASSERT(outer_s >= 0);
    ASSERT(inner_s >= 0);

    int outer = alea_add_cell(sys, 10, alea_surface_at(sys, outer_s)->neg_node,
                              mat, 1.0, 0);
    int inner = alea_add_cell(sys, 20, alea_surface_at(sys, inner_s)->neg_node,
                              mat, 1.0, 10);
    ASSERT(outer >= 0);
    ASSERT(inner >= 0);
    ASSERT_EQ(alea_set_fill(sys, outer, 10, 0), 0);

    alea_config_t cfg = alea_get_config(sys);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_HIER;
    alea_set_config(sys, &cfg);

    ASSERT_EQ(alea_volume_path_count(sys), 1);
    double volumes[1] = {0.0};
    double errors[1] = {0.0};
    ASSERT_EQ(alea_estimate_path_volumes(sys, 16, volumes, errors), 0);
    ASSERT(volumes[0] >= 0.0);
    ASSERT_NULL(sys->spatial_index);
    ASSERT_NOT_NULL(sys->hier_spatial_index);

    alea_destroy(sys);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_path_volume_estimation_sphere_smoke) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

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

    alea_config_t cfg = alea_get_config(sys);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_HIER;
    alea_set_config(sys, &cfg);

    ASSERT_EQ(alea_volume_path_count(sys), 1);

    double volumes[1] = {0.0};
    double errors[1] = {0.0};
    ASSERT_EQ(alea_estimate_path_volumes(sys, 4096, volumes, errors), 0);
    double expected = 4.0 * 3.14159265358979323846 / 3.0;
    ASSERT(volumes[0] > expected * 0.65);
    ASSERT(volumes[0] < expected * 1.35);
    ASSERT(errors[0] >= 0.0);
    ASSERT_NULL(sys->spatial_index);
    ASSERT_NOT_NULL(sys->hier_spatial_index);

    alea_destroy(sys);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(flat_spatial_index_fails_in_hier_mode) {
    mcnp_model_t* model = mcnp_load("tests/data/simple_fill.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    alea_config_t cfg = alea_get_config(sys);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_HIER;
    alea_set_config(sys, &cfg);

    ASSERT_EQ(alea_build_spatial_index(sys), -1);
    ASSERT_NULL(sys->spatial_index);

    mcnp_model_destroy(model);
}

TEST(auto_spatial_mode_avoids_flat_index_above_threshold) {
    alea_setenv("ALEA_SPATIAL_AUTO_CELL_THRESHOLD", "1", 1);
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    mcnp_model_t* model = mcnp_load("tests/data/simple_fill.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    alea_config_t cfg = alea_get_config(sys);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_AUTO;
    alea_set_config(sys, &cfg);

    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    ASSERT_NULL(sys->spatial_index);
    ASSERT_NOT_NULL(sys->hier_spatial_index);
    ASSERT_EQ(alea_build_spatial_index(sys), -1);
    ASSERT_EQ(alea_spatial_index_instance_count(sys), 0);

    double volumes[1] = {0.0};
    ASSERT_EQ(alea_estimate_instance_volumes(sys, 8, volumes, NULL), -1);

    mcnp_model_destroy(model);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
    alea_unsetenv("ALEA_SPATIAL_AUTO_CELL_THRESHOLD");
}

TEST(public_raycast_uses_hier_mode_without_flat_spatial_index) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    mcnp_model_t* model = mcnp_load("tests/data/simple_fill.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    alea_raycast_result_t flat;
    alea_raycast_result_t routed;
    alea_raycast_result_init(&flat);
    alea_raycast_result_init(&routed);

    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    ASSERT_EQ(alea_raycast(sys, -5.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                           10.0, &flat), 0);

    if (sys->spatial_index) {
        alea_spatial_index_free(sys->spatial_index);
        sys->spatial_index = NULL;
    }

    alea_config_t cfg = alea_get_config(sys);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_HIER;
    alea_set_config(sys, &cfg);

    ASSERT_EQ(alea_raycast(sys, -5.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                           10.0, &routed), 0);
    ASSERT_NULL(sys->spatial_index);
    ASSERT_NOT_NULL(sys->hier_spatial_index);
    assert_raycast_material_segments_match(&flat, &routed);
    ASSERT_EQ(alea_ray_is_occluded(sys, -5.0, 0.0, 0.0,
                                   1.0, 0.0, 0.0, 10.0), 1);
    ASSERT_NULL(sys->spatial_index);

    alea_raycast_result_free(&flat);
    alea_raycast_result_free(&routed);
    mcnp_model_destroy(model);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(spatial_mode_change_invalidates_flat_spatial_index) {
    mcnp_model_t* model = mcnp_load("tests/data/simple_fill.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    ASSERT_NOT_NULL(sys->spatial_index);

    alea_config_t cfg = alea_get_config(sys);
    cfg.spatial_mode = ALEA_SPATIAL_MODE_HIER;
    alea_set_config(sys, &cfg);

    ASSERT_NULL(sys->spatial_index);
    ASSERT_EQ(alea_spatial_index_instance_count(sys), 0);
    ASSERT_EQ(alea_build_spatial_index(sys), -1);

    double volumes[1] = {0.0};
    ASSERT_EQ(alea_estimate_instance_volumes(sys, 8, volumes, NULL), -1);

    mcnp_model_destroy(model);
}

static void assert_raycast_hits_match(const alea_raycast_result_t* flat,
                                      const alea_raycast_result_t* hier) {
    ASSERT_EQ(hier->hits.count, flat->hits.count);
    for (size_t i = 0; i < flat->hits.count; i++) {
        const alea_ray_hit_t* a = &flat->hits.data[i];
        const alea_ray_hit_t* b = &hier->hits.data[i];
        ASSERT_NEAR(b->t, a->t, 1e-9);
        ASSERT_EQ(b->surface_id, a->surface_id);
        ASSERT_EQ(b->primitive_id, a->primitive_id);
        ASSERT_NEAR(b->nx, a->nx, 1e-9);
        ASSERT_NEAR(b->ny, a->ny, 1e-9);
        ASSERT_NEAR(b->nz, a->nz, 1e-9);
    }
}

/* Phase 1 (PLAN_HIER_RAY_NEXT_STEPS): alea_raycast_hier_with_hits() returns
 * the same path segments as the segment-only fast path, plus a boundary hit
 * list. For a simple non-lattice model (root universe, identity transforms)
 * its hits match the flat alea_raycast() path exactly. */
TEST(hier_with_hits_matches_flat_for_disjoint_spheres) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    ASSERT(mat >= 0);

    /* Disjoint spheres strung along x so the ray crosses each exactly twice
     * with void gaps between them — unambiguous, all in the root universe. */
    for (int i = 0; i < 5; i++) {
        int s = alea_sphere_surface(sys, i + 1, (double)i * 5.0, 0.0, 0.0, 1.5);
        ASSERT(s >= 0);
        const alea_surface_entry_t* surf = alea_surface_at(sys, s);
        ASSERT_NOT_NULL(surf);
        int c = alea_add_cell(sys, i + 1, surf->neg_node, mat, 1.0, 0);
        ASSERT(c >= 0);
    }

    alea_raycast_result_t flat;
    alea_raycast_result_t hier;
    alea_raycast_result_init(&flat);
    alea_raycast_result_init(&hier);

    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    ASSERT_EQ(alea_raycast(sys, -5.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                           30.0, &flat), 0);

    if (sys->spatial_index) {
        alea_spatial_index_free(sys->spatial_index);
        sys->spatial_index = NULL;
    }

    ASSERT_EQ(alea_raycast_hier_with_hits(sys, -5.0, 0.0, 0.0,
                                          1.0, 0.0, 0.0, 30.0, &hier), 0);

    ASSERT_NULL(sys->spatial_index);
    ASSERT_NOT_NULL(sys->hier_spatial_index);

    /* Every sphere is crossed twice -> 10 boundary hits. */
    ASSERT_EQ(hier.hits.count, 10u);
    assert_raycast_material_segments_match(&flat, &hier);
    assert_raycast_hits_match(&flat, &hier);

    alea_raycast_result_free(&flat);
    alea_raycast_result_free(&hier);
    alea_destroy(sys);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

/* with_hits must not change the stepped segments relative to fast_segments,
 * even for a fill/lattice model with non-identity transforms. */
TEST(hier_with_hits_segments_match_fast_segments) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    alea_raycast_result_t segs;
    alea_raycast_result_t hits;
    alea_raycast_result_init(&segs);
    alea_raycast_result_init(&hits);

    ASSERT_EQ(alea_raycast_hier_fast_segments(sys, -1.5, 0.0, 0.0,
                                              1.0, 0.0, 0.0, 7.0, &segs), 0);
    ASSERT_EQ(alea_raycast_hier_with_hits(sys, -1.5, 0.0, 0.0,
                                          1.0, 0.0, 0.0, 7.0, &hits), 0);

    ASSERT_NULL(sys->spatial_index);
    assert_raycast_segments_match(&segs, &hits);

    /* fast_segments must not populate hits; with_hits should. */
    ASSERT_EQ(segs.hits.count, 0u);

    alea_raycast_result_free(&segs);
    alea_raycast_result_free(&hits);
    mcnp_model_destroy(model);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

/* Synthetic lattice DDA boundaries (surface_id == 0) are reported on segment
 * boundaries but never emitted as physical surface hits. */
TEST(hier_with_hits_skips_synthetic_lattice_boundaries) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    alea_raycast_result_t hits;
    alea_raycast_result_init(&hits);

    ASSERT_EQ(alea_raycast_hier_with_hits(sys, -1.5, 0.0, 0.0,
                                          1.0, 0.0, 0.0, 7.0, &hits), 0);

    /* No emitted hit may carry a synthetic (0) or "none" (-1) surface id. */
    for (size_t i = 0; i < hits.hits.count; i++) {
        ASSERT(hits.hits.data[i].surface_id > 0);
        ASSERT(hits.hits.data[i].primitive_id != ALEA_PRIMITIVE_ID_INVALID);
    }

    /* Segments may still record synthetic boundaries; confirm the model
     * actually exercises at least one so the assertion above is meaningful. */
    int saw_synthetic = 0;
    for (size_t i = 0; i < hits.segments.count; i++) {
        if (hits.segments.data[i].exit_surface_id == 0 ||
            hits.segments.data[i].enter_surface_id == 0) {
            saw_synthetic = 1;
            break;
        }
    }
    ASSERT(saw_synthetic);

    alea_raycast_result_free(&hits);
    mcnp_model_destroy(model);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST_MAIN()
