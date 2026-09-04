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
#include "core/alea_universe.h"
#include "raycast/raycast.h"

/* This file intentionally exercises deprecated flat-spatial compatibility APIs
 * as part of the migration gate. */
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

typedef struct {
    alea_system_t* sys;
    double* lengths;
    size_t path_count;
} volume_interval_length_probe_t;

typedef struct {
    size_t calls;
    size_t cancel_at;
    int invalid_receipt;
} volume_progress_probe_t;

static int cancel_volume_progress(size_t completed, size_t maximum,
                                  double maximum_error, void* context) {
    volume_progress_probe_t* probe = context;
    if (completed > maximum || maximum_error < 0.0)
        probe->invalid_receipt = 1;
    probe->calls++;
    return completed >= probe->cancel_at;
}

static int accumulate_indexed_interval_length(
        void* context,
        const alea_raycast_selected_interval_view_t* interval) {
    volume_interval_length_probe_t* probe = context;
    if (interval->cell_id < 0) return 0;
    uint64_t path_id = UINT64_MAX;
    int found = alea_volume_path_id_from_hier_path(
        probe->sys, interval->path, &path_id);
    if (found < 0 || (found > 0 && path_id >= probe->path_count)) return -1;
    if (found > 0)
        probe->lengths[path_id] += interval->t_exit - interval->t_enter;
    return 0;
}

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

TEST(hier_paths_record_logical_lattice_occurrences) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    const alea_cell_entry_t* lattice = NULL;
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        if (sys->cells.data[i].lat_type != 0 && sys->cells.data[i].lat_fill) {
            lattice = &sys->cells.data[i];
            break;
        }
    }
    ASSERT_NOT_NULL(lattice);

    const double points[2][3] = {{0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}};
    alea_lattice_location_t locations[2];
    alea_hier_ray_path_t paths[2];
    for (int n = 0; n < 2; n++) {
        alea_hier_cell_hit_t hit;
        ASSERT_EQ(alea_lattice_locate_point(sys, lattice,
                                             points[n][0], points[n][1], points[n][2],
                                             &locations[n]), 1);
        ASSERT_EQ(alea_hier_spatial_find_path_at_point(
                      sys, points[n][0], points[n][1], points[n][2],
                      &hit, &paths[n]), 1);
        ASSERT(paths[n].count >= 2);
        ASSERT_EQ(paths[n].entries[0].is_lattice, 1);
        ASSERT_EQ(paths[n].entries[0].lat_i, locations[n].i);
        ASSERT_EQ(paths[n].entries[0].lat_j, locations[n].j);
        ASSERT_EQ(paths[n].entries[0].lat_k, locations[n].k);

        uint64_t path_id = UINT64_MAX;
        ASSERT_EQ(alea_volume_path_id_from_hier_path(
                      sys, &paths[n], &path_id), 1);
        alea_volume_path_t canonical;
        ASSERT_EQ(alea_volume_path_at_point(
                      sys, points[n][0], points[n][1], points[n][2],
                      &canonical), 1);
        ASSERT_EQ(path_id, canonical.path_id);
    }
    ASSERT(locations[0].i != locations[1].i ||
           locations[0].j != locations[1].j ||
           locations[0].k != locations[1].k);

    mcnp_model_destroy(model);
}

TEST(hier_spatial_candidate_lookup_resolves_lattice_container) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    int lattice_cell = -1;
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        if (sys->cells.data[i].universe_id == 0 &&
            sys->cells.data[i].lat_type != 0 && sys->cells.data[i].lat_fill) {
            lattice_cell = (int)i;
            break;
        }
    }
    ASSERT(lattice_cell >= 0);
    ASSERT_EQ(alea_hier_spatial_index_build(sys), 0);
    ASSERT_EQ(alea_hier_spatial_find_cell_in_universe(sys, 0,
                                                       0.0, 0.0, 0.0),
              lattice_cell);

    mcnp_model_destroy(model);
}

/* Sweep a straight run of points through the coherent resolver, the way a grid
 * row does, and record the cell reported at each one. */
static void coherent_scan_cell_ids(alea_system_t* sys,
                                   alea_hier_coherence_ownership_t ownership,
                                   double x0, double dx, double y, double z,
                                   int n, int* out_cell_ids) {
    alea_hier_coherence_state_t state_a, state_b;
    alea_hier_coherence_state_t* previous = &state_a;
    alea_hier_coherence_state_t* current = &state_b;
    alea_hier_coherence_state_clear(previous);
    alea_hier_coherence_state_clear(current);

    for (int i = 0; i < n; i++) {
        double x = x0 + i * dx;
        alea_hier_cell_hit_t hit;
        alea_hier_coherence_kind_t kind;
        int rc = alea_hier_spatial_resolve_coherent(
            sys, x, y, z, i == 0 ? NULL : previous, ownership,
            current, &hit, &kind);
        out_cell_ids[i] = rc > 0 ? hit.hit.cell_id : -1;
        if (rc > 0) {
            /* Whatever the mode reports must actually contain the point. */
            ASSERT_EQ(alea_hier_spatial_check_path_containment(
                          sys, &current->path, current->path.count - 1,
                          x, y, z, NULL, NULL, NULL), 1);
        }
        alea_hier_coherence_state_t* tmp = previous;
        previous = current;
        current = tmp;
    }
}

/* Where every point has a unique owner, the cheap mode must be indistinguishable
 * from deck-order re-derivation and from a cold query. */
TEST(hier_coherent_ownership_matches_canonical_without_overlaps) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_hier_spatial_index_build(sys), 0);

    enum { N = 64 };
    const double x0 = -0.93, dx = 0.09, y = 0.17, z = 0.0;
    int coherent[N], canonical[N];
    coherent_scan_cell_ids(sys, ALEA_HIER_COH_OWNERSHIP_COHERENT,
                           x0, dx, y, z, N, coherent);
    coherent_scan_cell_ids(sys, ALEA_HIER_COH_OWNERSHIP_CANONICAL,
                           x0, dx, y, z, N, canonical);

    int saw_transition = 0;
    for (int i = 0; i < N; i++) {
        alea_hier_cell_hit_t cold;
        ASSERT_EQ(alea_hier_spatial_find_path_at_point(sys, x0 + i * dx, y, z,
                                                       &cold, NULL), 1);
        ASSERT_EQ(coherent[i], canonical[i]);
        ASSERT_EQ(coherent[i], cold.hit.cell_id);
        if (i > 0 && coherent[i] != coherent[i - 1]) saw_transition = 1;
    }
    /* A scan that never left one cell would prove nothing about path reuse. */
    ASSERT(saw_transition);

    mcnp_model_destroy(model);
}

/* In an overlap the modes are meant to disagree: coherent keeps the cell the
 * sweep entered, canonical re-derives the deck-first owner. */
TEST(hier_coherent_ownership_keeps_entered_cell_in_overlap) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_overlap.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_hier_spatial_index_build(sys), 0);

    /* Cells 1 and 2 are spheres at x=0 and x=3, both r=5, so on the axis they
     * share -2 < x < 5.  Sweeping leftwards enters the shared region from the
     * cell-2 side, against deck order. */
    enum { N = 5 };
    const double x0 = 6.0, dx = -1.0, y = 0.0, z = 0.0;
    int coherent[N], canonical[N];
    coherent_scan_cell_ids(sys, ALEA_HIER_COH_OWNERSHIP_COHERENT,
                           x0, dx, y, z, N, coherent);
    coherent_scan_cell_ids(sys, ALEA_HIER_COH_OWNERSHIP_CANONICAL,
                           x0, dx, y, z, N, canonical);

    /* x = 6: outside sphere 1, both modes must say cell 2. */
    ASSERT_EQ(coherent[0], 2);
    ASSERT_EQ(canonical[0], 2);
    /* x = 4, 3, 2: inside both.  Coherent carries cell 2 in, canonical resets
     * to the deck-first owner. */
    for (int i = 1; i < N; i++) {
        ASSERT_EQ(coherent[i], 2);
        ASSERT_EQ(canonical[i], 1);
    }

    mcnp_model_destroy(model);
}

TEST(hier_path_restart_rebuilds_changed_lattice_element) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    alea_hier_cell_hit_t start_hit;
    alea_hier_ray_path_t start_path;
    ASSERT_EQ(alea_hier_spatial_find_path_at_point(sys, 0.0, 0.0, 0.0,
                                                   &start_hit, &start_path), 1);
    ASSERT(start_path.count >= 2);
    ASSERT_EQ(start_path.entries[0].is_lattice, 1);

    alea_hier_cell_hit_t restarted_hit;
    alea_hier_ray_path_t restarted_path;
    ASSERT_EQ(alea_hier_spatial_find_path_from_parent(
                  sys, &start_path, 0, 2.0, 0.0, 0.0,
                  &restarted_hit, &restarted_path), 1);

    const alea_cell_entry_t* lattice =
        &sys->cells.data[start_path.entries[0].cell_index];
    alea_lattice_location_t location;
    ASSERT_EQ(alea_lattice_locate_point(sys, lattice, 2.0, 0.0, 0.0,
                                         &location), 1);
    alea_hier_cell_hit_t direct_hit;
    alea_hier_ray_path_t direct_path;
    ASSERT_EQ(alea_hier_path_enter_lattice_location(
                  sys, &start_path, 0, 2.0, 0.0, 0.0, &location,
                  &direct_hit, &direct_path), 1);

    alea_hier_cell_hit_t reference_hit;
    alea_hier_ray_path_t reference_path;
    ASSERT_EQ(alea_hier_spatial_find_path_at_point(sys, 2.0, 0.0, 0.0,
                                                   &reference_hit, &reference_path), 1);

    /* The worker-local resolver is the grid/DDA-facing entry point. It must
     * retain the lattice container and rebuild just the changed placement. */
    alea_hier_coherence_state_t state_a, state_b;
    alea_hier_coherence_state_clear(&state_a);
    alea_hier_coherence_state_clear(&state_b);
    alea_hier_cell_hit_t coherent_hit;
    alea_hier_coherence_kind_t coherence_kind;
    ASSERT_EQ(alea_hier_spatial_resolve_coherent(
                  sys, 0.0, 0.0, 0.0, NULL, ALEA_HIER_COH_OWNERSHIP_CANONICAL,
                  &state_a, &coherent_hit, &coherence_kind), 1);
    ASSERT_EQ(coherence_kind, ALEA_HIER_COH_ROOT_QUERY);
    ASSERT_EQ(alea_hier_spatial_resolve_coherent(
                  sys, 2.0, 0.0, 0.0, &state_a,
                  ALEA_HIER_COH_OWNERSHIP_CANONICAL, &state_b, &coherent_hit,
                  &coherence_kind), 1);
    ASSERT_EQ(coherence_kind, ALEA_HIER_COH_LATTICE_TRANSITION);
    ASSERT_EQ(coherent_hit.hit.cell_id, reference_hit.hit.cell_id);
    ASSERT_EQ(state_b.path.count, reference_path.count);

    ASSERT_EQ(restarted_hit.hit.cell_id, reference_hit.hit.cell_id);
    ASSERT_EQ(direct_hit.hit.cell_id, reference_hit.hit.cell_id);
    ASSERT_EQ(restarted_path.count, reference_path.count);
    ASSERT_EQ(direct_path.count, reference_path.count);
    for (int i = 0; i < reference_path.count; i++) {
        ASSERT_EQ(restarted_path.entries[i].cell_index,
                  reference_path.entries[i].cell_index);
        ASSERT_EQ(restarted_path.entries[i].lat_i,
                  reference_path.entries[i].lat_i);
        ASSERT_EQ(restarted_path.entries[i].lat_j,
                  reference_path.entries[i].lat_j);
        ASSERT_EQ(restarted_path.entries[i].lat_k,
                  reference_path.entries[i].lat_k);
        ASSERT_EQ(direct_path.entries[i].cell_index,
                  reference_path.entries[i].cell_index);
        ASSERT_EQ(direct_path.entries[i].lat_i,
                  reference_path.entries[i].lat_i);
        ASSERT_EQ(direct_path.entries[i].lat_j,
                  reference_path.entries[i].lat_j);
        ASSERT_EQ(direct_path.entries[i].lat_k,
                  reference_path.entries[i].lat_k);
        ASSERT_EQ(state_b.path.entries[i].cell_index,
                  reference_path.entries[i].cell_index);
        ASSERT_EQ(state_b.path.entries[i].lat_i,
                  reference_path.entries[i].lat_i);
        ASSERT_EQ(state_b.path.entries[i].lat_j,
                  reference_path.entries[i].lat_j);
        ASSERT_EQ(state_b.path.entries[i].lat_k,
                  reference_path.entries[i].lat_k);
        for (int m = 0; m < 12; m++) {
            ASSERT_EQ(restarted_path.entries[i].transform.m[m],
                      reference_path.entries[i].transform.m[m]);
            ASSERT_EQ(restarted_path.entries[i].transform.inv[m],
                      reference_path.entries[i].transform.inv[m]);
            ASSERT_EQ(direct_path.entries[i].transform.m[m],
                      reference_path.entries[i].transform.m[m]);
            ASSERT_EQ(direct_path.entries[i].transform.inv[m],
                      reference_path.entries[i].transform.inv[m]);
            ASSERT_EQ(state_b.path.entries[i].transform.m[m],
                      reference_path.entries[i].transform.m[m]);
            ASSERT_EQ(state_b.path.entries[i].transform.inv[m],
                      reference_path.entries[i].transform.inv[m]);
        }
    }

    mcnp_model_destroy(model);
}

TEST(lattice_location_round_trips_indices_and_enforces_window) {
    const struct {
        const char* path;
        int lat_type;
    } cases[] = {
        {"tests/data/mcnp_lattice_eval.mcnp", 1},
        {"tests/data/mcnp_hex_lattice.mcnp", 2},
    };

    for (size_t n = 0; n < sizeof(cases) / sizeof(cases[0]); n++) {
        mcnp_model_t* model = mcnp_load(cases[n].path);
        if (!model) SKIP("Test data file not found");
        alea_system_t* sys = model->sys;
        const alea_cell_entry_t* lattice = NULL;
        for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
            if (sys->cells.data[i].lat_type == cases[n].lat_type) {
                lattice = &sys->cells.data[i];
                break;
            }
        }
        ASSERT_NOT_NULL(lattice);

        alea_lattice_location_t point_location;
        ASSERT_EQ(alea_lattice_locate_point(sys, lattice, 0.0, 0.0, 0.0,
                                             &point_location), 1);

        alea_lattice_location_t indexed_location;
        ASSERT_EQ(alea_lattice_location_from_indices(
                      lattice, point_location.i, point_location.j,
                      point_location.k, &indexed_location), 1);
        ASSERT_EQ(indexed_location.fill_universe, point_location.fill_universe);
        ASSERT_EQ(indexed_location.linear_index, point_location.linear_index);
        ASSERT_EQ(indexed_location.ox, point_location.ox);
        ASSERT_EQ(indexed_location.oy, point_location.oy);
        ASSERT_EQ(indexed_location.oz, point_location.oz);

        ASSERT_EQ(alea_lattice_locate_point(sys, lattice, 0.0, 0.0, 2.0,
                                             &point_location), 0);
        mcnp_model_destroy(model);
    }
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

    /* Overlapping cells: the first containing cell in definition order wins,
     * matching the canonical point resolver (find_cell_recursive). */
    int ordered = alea_hier_spatial_find_ordered_cell_in_universe(
        sys, 0, 0.0, 0.0, 0.0, -1);
    ASSERT_EQ(ordered, small_cell);
    (void)large_cell;

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
                  sys, &path, 1, 0.1, 0.0, 0.0,
                  &transform, &lattice_cell_index, &lattice_transform), 1);
    ASSERT_EQ(lattice_cell_index, -1);
    ASSERT_EQ(alea_hier_spatial_check_path_containment(
                  sys, &path, 1, 20.0, 0.0, 0.0,
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

    /* Ray restarts use the coherent form: a valid ordinary fill prefix must
     * rebuild the same child path without a deck-order owner query. */
    ASSERT_EQ(alea_hier_spatial_find_path_from_parent_coherent(
                  sys, &path, 0, 0.1, 0.0, 0.0,
                  &refreshed_hit, &refreshed_path), 1);
    ASSERT_EQ(refreshed_hit.hit.cell_id, 2);
    ASSERT_EQ(refreshed_path.count, 2);
    ASSERT_EQ(refreshed_path.entries[0].cell_id, 1);
    ASSERT_EQ(refreshed_path.entries[1].cell_id, 2);

    alea_destroy(sys);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_path_containment_targets_exact_repeated_entry) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    ASSERT(mat >= 0);
    int sphere = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 100.0);
    ASSERT(sphere >= 0);
    int cell = alea_add_cell(sys, 1, alea_surface_at(sys, sphere)->neg_node,
                             mat, 1.0, 0);
    ASSERT(cell >= 0);

    /* Two occurrences of the same definition have different placements.
     * A cell index alone cannot identify which transform the caller means. */
    alea_hier_ray_path_t path;
    memset(&path, 0, sizeof(path));
    path.count = 2;
    for (int i = 0; i < path.count; i++) {
        path.entries[i].cell_index = (uint32_t)cell;
        path.entries[i].cell_id = 1;
        path.entries[i].universe_id = 0;
        path.entries[i].depth = i;
        alea_matrix_identity(&path.entries[i].transform);
    }
    path.entries[1].transform.m[3] = 10.0;
    path.entries[1].transform.has_inverse = false;
    ASSERT(alea_matrix_invert(&path.entries[1].transform));

    alea_matrix_t first_transform;
    alea_matrix_t second_transform;
    ASSERT_EQ(alea_hier_spatial_check_path_containment(
                  sys, &path, 0, 0.0, 0.0, 0.0,
                  &first_transform, NULL, NULL), 1);
    ASSERT_EQ(alea_hier_spatial_check_path_containment(
                  sys, &path, 1, 0.0, 0.0, 0.0,
                  &second_transform, NULL, NULL), 1);
    ASSERT_EQ(first_transform.inv[3], 0.0);
    ASSERT_EQ(second_transform.inv[3], -10.0);
    ASSERT_EQ(alea_hier_spatial_check_path_containment(
                  sys, &path, 2, 0.0, 0.0, 0.0,
                  NULL, NULL, NULL), -1);

    alea_destroy(sys);
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

TEST(hier_spatial_chain_region_query_carries_lattice_occurrence) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_hier_spatial_index_build(sys), 0);

    alea_bbox_t query = {-5.0 / 3.0, 5.0 / 3.0,
                         -5.0 / 3.0, 5.0 / 3.0, -0.1, 0.1};
    alea_hier_spatial_chain_hit_t hits[64];
    int n = alea_hier_spatial_query_region_chain(sys, &query, hits, 64);
    ASSERT(n > 0);

    int first_i = 0;
    int distinct_i = 0;
    int terminals = 0;
    int found_bounded_cylinder = 0;
    for (int h = 0; h < n; h++) {
        if (hits[h].hit.depth != 1) continue;
        ASSERT_EQ(hits[h].ancestor_count, 1);
        ASSERT_EQ(hits[h].ancestor_is_lattice[0], 1);
        ASSERT(hits[h].ancestor_lattice_fill_universes[0] > 0);
        if (hits[h].hit.cell_id == 1)
            found_bounded_cylinder = 1;
        if (terminals == 0) {
            first_i = hits[h].ancestor_lattice_i[0];
        } else if (hits[h].ancestor_lattice_i[0] != first_i) {
            distinct_i = 1;
        }
        terminals++;
    }
    ASSERT(terminals > 0);
    ASSERT_MSG(found_bounded_cylinder,
               "an infinite-axis cylinder must retain its finite XY bounds");
    ASSERT_MSG(distinct_i,
               "region chain should distinguish repeated lattice elements");

    mcnp_model_destroy(model);
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


    alea_raycast_result_t result;
    alea_raycast_result_init(&result);
    ASSERT_EQ(alea_raycast_hier(sys, -5.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                                10.0, &result), 0);

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


    alea_raycast_result_t flat;
    alea_raycast_result_t hier;
    alea_raycast_result_init(&flat);
    alea_raycast_result_init(&hier);

    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    ASSERT_EQ(alea_raycast(sys, -1.5, 0.0, 0.0, 1.0, 0.0, 0.0,
                           7.0, &flat), 0);


    ASSERT_EQ(alea_raycast_hier_fast_segments(sys, -1.5, 0.0, 0.0,
                                              1.0, 0.0, 0.0,
                                              7.0, &hier), 0);

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
    alea_set_config(sys, &cfg);

    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
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
    alea_set_config(sys, &cfg);

    ASSERT_NULL(sys->hier_spatial_index);

    alea_cell_hit_t hits[8];
    int n = alea_find_all_cells(sys, 0.0, 0.0, 0.0, hits, 8);
    ASSERT(n > 0);
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
    alea_set_config(sys, &cfg);

    alea_query_acceleration_stats_t stats;
    ASSERT_EQ(alea_query_acceleration_stats(sys, &stats), 0);
    ASSERT(!stats.built);

    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    ASSERT_EQ(alea_query_acceleration_stats(sys, &stats), 0);
    ASSERT(stats.built);
    ASSERT_EQ(stats.hier_universe_count, 1);
    ASSERT_EQ(stats.hier_blas_count, 1);
    ASSERT_EQ(stats.hier_blas_cell_count, 3);
    ASSERT_EQ(stats.hier_placement_count, 1);
    ASSERT(stats.memory_bytes > 0);
    ASSERT_NOT_NULL(sys->hier_spatial_index);

    alea_destroy(sys);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
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

TEST(volume_path_structural_resolver_avoids_global_index) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    ASSERT(mat >= 0);
    int surface = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 2.0);
    ASSERT(surface >= 0);
    alea_node_id_t interior = alea_halfspace(sys, surface, -1);
    ASSERT(interior != ALEA_NODE_ID_INVALID);
    int cell = alea_add_cell(sys, 1, interior, mat, 1.0, 0);
    ASSERT(cell >= 0);
    ASSERT_NULL(sys->volume_path_index);

    alea_volume_path_t path;
    ASSERT_EQ(alea_volume_path_resolve_at_point(
        sys, 0.0, 0.0, 0.0, &path), 1);
    ASSERT_EQ(path.terminal_cell_index, cell);
    ASSERT_EQ(path.path_id, UINT64_MAX);
    ASSERT_NULL(sys->volume_path_index);

    alea_destroy(sys);
}

TEST(volume_path_target_resolver_finds_ancestor_and_terminal) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    ASSERT(mat >= 0);
    int outer_surface = alea_sphere_surface(sys, 1, 10.0, 0.0, 0.0, 4.0);
    int child_surface = alea_sphere_surface(sys, 2, 0.0, 0.0, 0.0, 1.0);
    ASSERT(outer_surface >= 0);
    ASSERT(child_surface >= 0);
    alea_node_id_t outer = alea_halfspace(sys, outer_surface, -1);
    alea_node_id_t child_region = alea_halfspace(sys, child_surface, -1);
    ASSERT(outer != ALEA_NODE_ID_INVALID);
    ASSERT(child_region != ALEA_NODE_ID_INVALID);

    int container = alea_add_cell(sys, 10, outer, mat, 1.0, 0);
    int child = alea_add_cell(sys, 20, child_region, mat, 1.0, 7);
    ASSERT(container >= 0);
    ASSERT(child >= 0);
    const double translation[3] = {10.0, 0.0, 0.0};
    ASSERT_EQ(alea_add_transform(sys, 70, translation, 3, 0), 0);
    ASSERT_EQ(alea_set_cell_fill(sys, container, 7, 70), 0);
    ASSERT_NULL(sys->volume_path_index);

    alea_volume_path_t path;
    ASSERT_EQ(alea_volume_path_resolve_cell_at_point(
        sys, 10.0, 0.0, 0.0, 10, 0, &path), 1);
    ASSERT_EQ(path.terminal_cell_index, container);
    ASSERT_EQ(path.depth, 0);
    ASSERT_EQ(path.ancestor_count, 0);
    ASSERT_EQ(path.world_to_local[3], 0.0);

    ASSERT_EQ(alea_volume_path_resolve_cell_at_point(
        sys, 10.0, 0.0, 0.0, 20, 7, &path), 1);
    ASSERT_EQ(path.terminal_cell_index, child);
    ASSERT_EQ(path.depth, 1);
    ASSERT_EQ(path.ancestor_count, 1);
    ASSERT_EQ(path.ancestor_cell_indices[0], container);
    ASSERT_EQ(path.path_id, UINT64_MAX);
    ASSERT_EQ(path.world_to_local[3], -10.0);
    ASSERT_NULL(sys->volume_path_index);

    ASSERT_EQ(alea_volume_path_resolve_cell_at_point(
        sys, 30.0, 0.0, 0.0, 20, 7, &path), 0);

    alea_destroy(sys);
}

TEST(volume_path_target_resolver_reports_repeated_overlap) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    ASSERT(mat >= 0);
    int outer_surface = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 4.0);
    int child_surface = alea_sphere_surface(sys, 2, 0.0, 0.0, 0.0, 1.0);
    ASSERT(outer_surface >= 0);
    ASSERT(child_surface >= 0);
    alea_node_id_t outer = alea_halfspace(sys, outer_surface, -1);
    alea_node_id_t child_region = alea_halfspace(sys, child_surface, -1);
    ASSERT(outer != ALEA_NODE_ID_INVALID);
    ASSERT(child_region != ALEA_NODE_ID_INVALID);

    int left = alea_add_cell(sys, 10, outer, mat, 1.0, 0);
    int right = alea_add_cell(sys, 11, outer, mat, 1.0, 0);
    int child = alea_add_cell(sys, 20, child_region, mat, 1.0, 7);
    ASSERT(left >= 0);
    ASSERT(right >= 0);
    ASSERT(child >= 0);
    ASSERT_EQ(alea_set_cell_fill(sys, left, 7, 0), 0);
    ASSERT_EQ(alea_set_cell_fill(sys, right, 7, 0), 0);

    alea_volume_path_t path;
    ASSERT_EQ(alea_volume_path_resolve_cell_at_point(
        sys, 0.0, 0.0, 0.0, 20, 7, &path), 2);
    ASSERT_NULL(sys->volume_path_index);

    alea_destroy(sys);
}

TEST(volume_path_point_sets_aggregate_structural_occurrences) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    int left_surface = alea_sphere_surface(sys, 1, -10.0, 0.0, 0.0, 4.0);
    int right_surface = alea_sphere_surface(sys, 2, 10.0, 0.0, 0.0, 4.0);
    int child_surface = alea_sphere_surface(sys, 3, 0.0, 0.0, 0.0, 1.0);
    ASSERT(mat >= 0);
    ASSERT(left_surface >= 0);
    ASSERT(right_surface >= 0);
    ASSERT(child_surface >= 0);
    int left = alea_add_cell(
        sys, 10, alea_halfspace(sys, left_surface, -1), mat, 1.0, 0);
    int right = alea_add_cell(
        sys, 11, alea_halfspace(sys, right_surface, -1), mat, 1.0, 0);
    int child = alea_add_cell(
        sys, 20, alea_halfspace(sys, child_surface, -1), mat, 1.0, 7);
    ASSERT(left >= 0);
    ASSERT(right >= 0);
    ASSERT(child >= 0);
    const double left_translation[3] = {-10.0, 0.0, 0.0};
    const double right_translation[3] = {10.0, 0.0, 0.0};
    ASSERT_EQ(alea_add_transform(sys, 70, left_translation, 3, 0), 0);
    ASSERT_EQ(alea_add_transform(sys, 71, right_translation, 3, 0), 0);
    ASSERT_EQ(alea_set_cell_fill(sys, left, 7, 70), 0);
    ASSERT_EQ(alea_set_cell_fill(sys, right, 7, 71), 0);

    const double points[] = {
        -10.0, 0.0, 0.0,
        -10.1, 0.0, 0.0,
         10.0, 0.0, 0.0,
         30.0, 0.0, 0.0,
    };
    const alea_volume_path_point_set_t sets[] = {
        {.point_offset = 0, .point_count = 2,
         .target_cell_id = 20, .target_universe_id = 7},
        {.point_offset = 0, .point_count = 3,
         .target_cell_id = 20, .target_universe_id = 7},
        {.point_offset = 3, .point_count = 1,
         .target_cell_id = 20, .target_universe_id = 7},
        {.point_offset = 0, .point_count = 0,
         .target_cell_id = 20, .target_universe_id = 7},
    };
    alea_volume_path_point_set_result_t serial[4];
    alea_volume_path_point_set_batch_stats_t serial_stats;
    ASSERT_EQ(alea_volume_path_resolve_cell_point_sets(
        sys, points, 4, sets, 4, 4, 0, serial, &serial_stats), 0);
    ASSERT_EQ(serial_stats.point_set_count, 4);
    ASSERT_EQ(serial_stats.completed_point_set_count, 4);
    ASSERT_EQ(serial_stats.requested_workers, 4);
    ASSERT_EQ(serial_stats.actual_workers, 1);
    ASSERT(serial_stats.reserved_scratch_bytes_per_worker >=
           3 * sizeof(alea_volume_path_t));
    ASSERT_EQ(serial_stats.reserved_parallel_scratch_bytes,
              serial_stats.reserved_scratch_bytes_per_worker);

    ASSERT_EQ(serial[0].status, ALEA_VOLUME_PATH_UNIQUE);
    ASSERT_EQ(serial[0].tested_point_count, 2);
    ASSERT_EQ(serial[0].matching_occurrence_count, 1);
    ASSERT_EQ(serial[0].unique_path.terminal_cell_index, child);
    ASSERT_EQ(serial[0].unique_path.ancestor_cell_indices[0], left);

    ASSERT_EQ(serial[1].status, ALEA_VOLUME_PATH_MULTIPLE);
    ASSERT_EQ(serial[1].tested_point_count, 3);
    ASSERT_EQ(serial[1].matching_occurrence_count, 2);
    ASSERT_EQ(serial[2].status, ALEA_VOLUME_PATH_NONE);
    ASSERT_EQ(serial[2].tested_point_count, 1);
    ASSERT_EQ(serial[2].matching_occurrence_count, 0);
    ASSERT_EQ(serial[3].status, ALEA_VOLUME_PATH_NONE);
    ASSERT_EQ(serial[3].tested_point_count, 0);

    alea_volume_path_point_set_result_t parallel[4];
    alea_volume_path_point_set_batch_stats_t parallel_stats;
    ASSERT_EQ(alea_volume_path_resolve_cell_point_sets(
        sys, points, 4, sets, 4, 4, UINT64_MAX,
        parallel, &parallel_stats), 0);
    for (size_t set = 0; set < 4; set++) {
        ASSERT_EQ(parallel[set].status, serial[set].status);
        ASSERT_EQ(parallel[set].tested_point_count,
                  serial[set].tested_point_count);
        ASSERT_EQ(parallel[set].matching_occurrence_count,
                  serial[set].matching_occurrence_count);
    }
    ASSERT(parallel_stats.actual_workers >= 1);
    ASSERT(parallel_stats.actual_workers <= 4);
    ASSERT_NULL(sys->volume_path_index);

    alea_volume_path_point_set_t invalid = {
        .point_offset = 4, .point_count = 1,
        .target_cell_id = 20, .target_universe_id = 7};
    alea_volume_path_point_set_result_t invalid_result;
    alea_volume_path_point_set_batch_stats_t invalid_stats;
    ASSERT_EQ(alea_volume_path_resolve_cell_point_sets(
        sys, points, 4, &invalid, 1, 1, 0,
        &invalid_result, &invalid_stats), -1);

    alea_destroy(sys);
}

TEST(volume_path_point_sets_preserve_individual_ambiguity) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    int outer_surface = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 4.0);
    int child_surface = alea_sphere_surface(sys, 2, 0.0, 0.0, 0.0, 1.0);
    ASSERT(mat >= 0);
    int left = alea_add_cell(
        sys, 10, alea_halfspace(sys, outer_surface, -1), mat, 1.0, 0);
    int right = alea_add_cell(
        sys, 11, alea_halfspace(sys, outer_surface, -1), mat, 1.0, 0);
    ASSERT(left >= 0);
    ASSERT(right >= 0);
    ASSERT(alea_add_cell(
        sys, 20, alea_halfspace(sys, child_surface, -1), mat, 1.0, 7) >= 0);
    ASSERT_EQ(alea_set_cell_fill(sys, left, 7, 0), 0);
    ASSERT_EQ(alea_set_cell_fill(sys, right, 7, 0), 0);

    const double point[] = {0.0, 0.0, 0.0};
    const alea_volume_path_point_set_t set = {
        .point_offset = 0, .point_count = 1,
        .target_cell_id = 20, .target_universe_id = 7};
    alea_volume_path_point_set_result_t result;
    alea_volume_path_point_set_batch_stats_t stats;
    ASSERT_EQ(alea_volume_path_resolve_cell_point_sets(
        sys, point, 1, &set, 1, 1, 0, &result, &stats), 0);
    ASSERT_EQ(result.status, ALEA_VOLUME_PATH_MULTIPLE);
    ASSERT_EQ(result.tested_point_count, 1);
    ASSERT_EQ(result.matching_occurrence_count, 2);

    alea_destroy(sys);
}

TEST(volume_path_target_resolver_enforces_lattice_window) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    alea_volume_path_t path;
    ASSERT_EQ(alea_volume_path_resolve_cell_at_point(
        sys, 0.0, 0.0, 0.0, 1, 1, &path), 1);
    ASSERT_EQ(path.terminal_cell_id, 1);
    ASSERT_EQ(path.universe_id, 1);
    ASSERT_EQ(path.lattice_step_count, 1);

    ASSERT_EQ(alea_volume_path_resolve_cell_at_point(
        sys, 0.0, 0.0, 2.0, 1, 1, &path), 0);
    ASSERT_EQ(alea_volume_path_resolve_cell_at_point(
        sys, 0.0, 0.0, 2.0, 100, 0, &path), 0);
    ASSERT_NULL(sys->volume_path_index);

    mcnp_model_destroy(model);
}

TEST(volume_path_transform_evidence_resolves_outside_reported_cell) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    int outer_surface = alea_sphere_surface(sys, 1, 10.0, 0.0, 0.0, 4.0);
    int child_surface = alea_sphere_surface(sys, 2, 0.0, 0.0, 0.0, 1.0);
    ASSERT(mat >= 0);
    ASSERT(outer_surface >= 0);
    ASSERT(child_surface >= 0);
    int container = alea_add_cell(
        sys, 10, alea_halfspace(sys, outer_surface, -1), mat, 1.0, 0);
    int child = alea_add_cell(
        sys, 20, alea_halfspace(sys, child_surface, -1), mat, 1.0, 7);
    ASSERT(container >= 0);
    ASSERT(child >= 0);
    const double translation[3] = {10.0, 0.0, 0.0};
    ASSERT_EQ(alea_add_transform(sys, 70, translation, 3, 0), 0);
    ASSERT_EQ(alea_set_cell_fill(sys, container, 7, 70), 0);

    const double world_point[3] = {12.0, 0.0, 0.0};
    const double local_point[3] = {2.0, 0.0, 0.0};
    const double direction[3] = {1.0, 0.0, 0.0};
    const double precision[3] = {1e-3, 1e-3, 1e-3};
    alea_volume_path_t path;
    ASSERT_EQ(alea_volume_path_resolve_cell_at_point(
        sys, 12.0, 0.0, 0.0, 20, 7, &path), 0);
    ASSERT_EQ(alea_volume_path_resolve_cell_from_transform_evidence(
        sys, 20, 7, world_point, direction, local_point, direction,
        precision, precision, precision, precision, &path), 1);
    ASSERT_EQ(path.terminal_cell_index, child);
    ASSERT_EQ(path.depth, 1);
    ASSERT_EQ(path.ancestor_count, 1);
    ASSERT_EQ(path.ancestor_cell_indices[0], container);
    ASSERT_EQ(path.world_to_local[3], -10.0);
    ASSERT_NULL(sys->volume_path_index);

    alea_destroy(sys);
}

TEST(volume_path_transform_evidence_reports_ambiguous_placements) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    int outer_surface = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 4.0);
    int child_surface = alea_sphere_surface(sys, 2, 0.0, 0.0, 0.0, 1.0);
    ASSERT(mat >= 0);
    int left = alea_add_cell(
        sys, 10, alea_halfspace(sys, outer_surface, -1), mat, 1.0, 0);
    int right = alea_add_cell(
        sys, 11, alea_halfspace(sys, outer_surface, -1), mat, 1.0, 0);
    ASSERT(left >= 0);
    ASSERT(right >= 0);
    ASSERT(alea_add_cell(
        sys, 20, alea_halfspace(sys, child_surface, -1), mat, 1.0, 7) >= 0);
    ASSERT_EQ(alea_set_cell_fill(sys, left, 7, 0), 0);
    ASSERT_EQ(alea_set_cell_fill(sys, right, 7, 0), 0);

    const double point[3] = {2.0, 0.0, 0.0};
    const double direction[3] = {1.0, 0.0, 0.0};
    const double precision[3] = {1e-3, 1e-3, 1e-3};
    alea_volume_path_t path;
    ASSERT_EQ(alea_volume_path_resolve_cell_from_transform_evidence(
        sys, 20, 7, point, direction, point, direction,
        precision, precision, precision, precision, &path), 2);
    ASSERT_NULL(sys->volume_path_index);

    alea_destroy(sys);
}

TEST(volume_path_transform_evidence_uses_ancestor_to_disambiguate) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);

    int mat = alea_add_material(sys, 1);
    int left_surface = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 4.0);
    int right_surface = alea_sphere_surface(sys, 2, 10.0, 0.0, 0.0, 4.0);
    int child_surface = alea_sphere_surface(sys, 3, 0.0, 0.0, 0.0, 1.0);
    ASSERT(mat >= 0);
    int left = alea_add_cell(
        sys, 10, alea_halfspace(sys, left_surface, -1), mat, 1.0, 0);
    int right = alea_add_cell(
        sys, 11, alea_halfspace(sys, right_surface, -1), mat, 1.0, 0);
    ASSERT(left >= 0);
    ASSERT(right >= 0);
    ASSERT(alea_add_cell(
        sys, 20, alea_halfspace(sys, child_surface, -1), mat, 1.0, 7) >= 0);
    ASSERT_EQ(alea_set_cell_fill(sys, left, 7, 0), 0);
    ASSERT_EQ(alea_set_cell_fill(sys, right, 7, 0), 0);

    const double point[3] = {2.0, 0.0, 0.0};
    const double direction[3] = {1.0, 0.0, 0.0};
    const double precision[3] = {1e-3, 1e-3, 1e-3};
    alea_volume_path_t path;
    ASSERT_EQ(alea_volume_path_resolve_cell_from_transform_evidence(
        sys, 20, 7, point, direction, point, direction,
        precision, precision, precision, precision, &path), 1);
    ASSERT_EQ(path.ancestor_count, 1);
    ASSERT_EQ(path.ancestor_cell_indices[0], left);

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
    alea_set_config(sys, &cfg);

    alea_error_clear();
    ASSERT_EQ(alea_volume_path_count(sys), 0);
    ASSERT_EQ(alea_error_code(), (int)ALEA_ERR_INVALID_STATE);
    ASSERT_NULL(sys->volume_path_index);

    double volumes[2] = {0.0, 0.0};
    ASSERT_EQ(alea_estimate_volumes(sys, 8, volumes, NULL), -1);
    ASSERT_EQ(alea_error_code(), (int)ALEA_ERR_INVALID_STATE);

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
    alea_set_config(sys, &cfg);

    size_t count = alea_volume_path_count(sys);
    ASSERT_EQ(count, 2);
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

    alea_destroy(sys);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_volume_interval_paths_resolve_repeated_fills) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int mat = alea_add_material(sys, 1);
    ASSERT(mat >= 0);

    int left_s = alea_sphere_surface(sys, 1, -5.0, 0.0, 0.0, 3.0);
    int right_s = alea_sphere_surface(sys, 2, 5.0, 0.0, 0.0, 3.0);
    int child_s = alea_sphere_surface(sys, 3, 0.0, 0.0, 0.0, 1.0);
    ASSERT(left_s >= 0 && right_s >= 0 && child_s >= 0);

    int left = alea_add_cell(sys, 10, alea_surface_at(sys, left_s)->neg_node,
                             ALEA_MATERIAL_VOID, 0.0, 0);
    int right = alea_add_cell(sys, 20, alea_surface_at(sys, right_s)->neg_node,
                              ALEA_MATERIAL_VOID, 0.0, 0);
    ASSERT(left >= 0 && right >= 0);
    const alea_node_id_t child_region =
        alea_surface_at(sys, child_s)->neg_node;
    int child = alea_add_cell(sys, 30, child_region, mat, 1.0, 10);
    ASSERT(child >= 0);
    ASSERT(alea_add_cell(sys, 31, alea_complement(sys, child_region),
                         ALEA_MATERIAL_VOID, 0.0, 10) >= 0);

    const double left_tr[3] = {-5.0, 0.0, 0.0};
    const double right_tr[3] = {5.0, 0.0, 0.0};
    ASSERT_EQ(alea_add_transform(sys, 101, left_tr, 3, 0), 0);
    ASSERT_EQ(alea_add_transform(sys, 102, right_tr, 3, 0), 0);
    ASSERT_EQ(alea_set_fill(sys, left, 10, 101), 0);
    ASSERT_EQ(alea_set_fill(sys, right, 10, 102), 0);

    const size_t path_count = alea_volume_path_count(sys);
    ASSERT_EQ(path_count, 4);
    double* streamed = calloc(path_count, sizeof(*streamed));
    alea_volume_path_t* paths = calloc(path_count, sizeof(*paths));
    ASSERT_NOT_NULL(streamed);
    ASSERT_NOT_NULL(paths);
    ASSERT_EQ(alea_volume_paths_get(sys, paths, path_count), path_count);

    alea_ray_t ray;
    ASSERT_EQ(alea_ray_init(&ray, -10.0, 0.0, 0.0, 1.0, 0.0, 0.0), 0);
    alea_raycast_result_t scratch;
    alea_raycast_result_init(&scratch);
    volume_interval_length_probe_t probe = {
        .sys = sys, .lengths = streamed, .path_count = path_count
    };
    ASSERT_EQ(alea_raycast_hier_visit_intervals_nocache(
                  sys, &ray, 20.0, &scratch,
                  accumulate_indexed_interval_length, &probe), 0);

    size_t child_occurrences = 0;
    for (size_t i = 0; i < path_count; i++) {
        if (paths[i].terminal_cell_index != child) continue;
        ASSERT_NEAR(streamed[paths[i].path_id], 2.0, 1e-12);
        child_occurrences++;
    }
    ASSERT_EQ(child_occurrences, 2);
    ASSERT_EQ(scratch.segments.count, 0);
    ASSERT_EQ(scratch.hits.count, 0);

    alea_raycast_result_free(&scratch);
    free(paths);
    free(streamed);
    alea_destroy(sys);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_volume_paths_distinguish_rect_lattice_elements) {
    alea_setenv("ALEA_HIER_BLAS_THRESHOLD", "1", 1);

    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_eval.mcnp");
    if (!model) SKIP("Test data file not found");
    alea_system_t* sys = model->sys;

    alea_config_t cfg = alea_get_config(sys);
    alea_set_config(sys, &cfg);

    size_t count = alea_volume_path_count(sys);
    ASSERT_EQ(count, 18);
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
    ASSERT_EQ(alea_estimate_volumes(sys, 4096, volumes, NULL), 0);

    int nonzero_u1_moderator = 0;
    int nonzero_u3_moderator = 0;
    for (size_t i = 0; i < count; i++) {
        if (paths[i].terminal_cell_id == 2 && volumes[i] > 0.0) nonzero_u1_moderator++;
        if (paths[i].terminal_cell_id == 4 && volumes[i] > 0.0) nonzero_u3_moderator++;
    }
    ASSERT(nonzero_u1_moderator >= 2);
    ASSERT(nonzero_u3_moderator >= 2);

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
    alea_set_config(sys, &cfg);

    ASSERT_EQ(alea_volume_path_count(sys), 1);
    double volumes[1] = {0.0};
    double errors[1] = {0.0};
    ASSERT_EQ(alea_estimate_volumes(sys, 16, volumes, errors), 0);
    ASSERT(volumes[0] >= 0.0);
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
    alea_set_config(sys, &cfg);

    ASSERT_EQ(alea_volume_path_count(sys), 1);

    double volumes[1] = {0.0};
    double errors[1] = {0.0};
    ASSERT_EQ(alea_estimate_volumes(sys, 4096, volumes, errors), 0);
    double expected = 4.0 * 3.14159265358979323846 / 3.0;
    ASSERT(volumes[0] > expected * 0.65);
    ASSERT(volumes[0] < expected * 1.35);
    ASSERT(errors[0] >= 0.0);
    ASSERT_NOT_NULL(sys->hier_spatial_index);

    alea_destroy(sys);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
}

TEST(hier_volume_estimation_options_are_reproducible_and_cancellable) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int mat = alea_add_material(sys, 1);
    int surface = alea_sphere_surface(sys, 1, 0.0, 0.0, 0.0, 1.0);
    ASSERT(mat >= 0 && surface >= 0);
    ASSERT(alea_add_cell(sys, 1, alea_surface_at(sys, surface)->neg_node,
                         mat, 1.0, 0) >= 0);

    alea_volume_estimate_options_t options;
    alea_volume_estimate_options_init(&options);
    options.max_rays = 400;
    options.batch_size = 100;
    options.seed = UINT64_C(123456789);
    options.requested_workers = 1;
    volume_progress_probe_t progress = {.cancel_at = 200};
    options.progress = cancel_volume_progress;
    options.progress_user_data = &progress;

    double first_volume[1], first_error[1];
    alea_volume_estimate_stats_t first_stats;
    ASSERT_EQ(alea_estimate_volumes_ex(
        sys, &options, first_volume, first_error, &first_stats), 0);
    ASSERT_EQ(progress.calls, 2);
    ASSERT_EQ(progress.invalid_receipt, 0);
    ASSERT_EQ(first_stats.rays_completed, 200);
    ASSERT(first_stats.cancelled);
    ASSERT(!first_stats.converged);

    progress.calls = 0;
    double second_volume[1], second_error[1];
    alea_volume_estimate_stats_t second_stats;
    ASSERT_EQ(alea_estimate_volumes_ex(
        sys, &options, second_volume, second_error, &second_stats), 0);
    ASSERT_EQ(first_volume[0], second_volume[0]);
    ASSERT_EQ(first_error[0], second_error[0]);
    ASSERT_EQ(second_stats.rays_completed, first_stats.rays_completed);

    options.progress = NULL;
    options.progress_user_data = NULL;
    options.target_rel_error = 1.0;
    alea_volume_estimate_stats_t converged_stats;
    ASSERT_EQ(alea_estimate_volumes_ex(
        sys, &options, second_volume, second_error, &converged_stats), 0);
    ASSERT(converged_stats.converged);
    ASSERT(!converged_stats.cancelled);
    ASSERT_EQ(converged_stats.rays_completed, options.batch_size);

    alea_destroy(sys);
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


    alea_config_t cfg = alea_get_config(sys);
    alea_set_config(sys, &cfg);

    ASSERT_EQ(alea_raycast(sys, -5.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                           10.0, &routed), 0);
    ASSERT_NOT_NULL(sys->hier_spatial_index);
    assert_raycast_material_segments_match(&flat, &routed);
    ASSERT_EQ(alea_ray_is_occluded(sys, -5.0, 0.0, 0.0,
                                   1.0, 0.0, 0.0, 10.0), 1);

    alea_raycast_result_free(&flat);
    alea_raycast_result_free(&routed);
    mcnp_model_destroy(model);
    alea_unsetenv("ALEA_HIER_BLAS_THRESHOLD");
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


    ASSERT_EQ(alea_raycast_hier_with_hits(sys, -5.0, 0.0, 0.0,
                                          1.0, 0.0, 0.0, 30.0, &hier), 0);

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
