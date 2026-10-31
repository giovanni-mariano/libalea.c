// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#include "alea_test.h"
#include "alea.h"
#include "alea_mcnp.h"
#include "alea_raycast.h"
#include "raycast/raycast.h"
#include "raycast/ray_intersect.h"
#include "raycast/recovery_context.h"
#include <float.h>
#include "core/alea_system.h"

/* Explicitly disable both prepared levels when comparing tree traversal. */
static uint32_t** detach_recovery_lists(alea_system_t* sys) {
    uint32_t** saved = malloc(sys->universes.count * sizeof(*saved));
    if (!saved) abort();
    for (size_t u = 0; u < sys->universes.count; u++) {
        saved[u] = sys->universes.data[u].recovery_surfaces;
        sys->universes.data[u].recovery_surfaces = NULL;
    }
    return saved;
}
static void restore_recovery_lists(alea_system_t* sys, uint32_t** saved) {
    for (size_t u = 0; u < sys->universes.count; u++)
        sys->universes.data[u].recovery_surfaces = saved[u];
    free(saved);
}

TEST(ownership_tracks_preserve_thin_same_material_regions_and_occurrences) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int material = alea_add_material(sys, 1);
    const double edges[] = {0, 1, 1.00001, 2};
    for (int i = 0; i < 3; ++i) {
        int surface = alea_box_surface(sys, i+1, edges[i], edges[i+1], -1, 1, -1, 1);
        ASSERT(surface >= 0);
        ASSERT(alea_add_cell(sys, i+1, alea_halfspace(sys, surface, -1),
                             material, -1, 0) >= 0);
    }
    const double origins[] = {0.25, 0, 0, 1.75, 0, 0};
    const double directions[] = {1, 0, 0, -1, 0, 0};
    alea_raycast_batch_options_t options = {
        .struct_size = sizeof(options), .fields = ALEA_RAY_BATCH_FULL_PATHS,
        .projected_depth = -1, .max_segments = 6, .max_path_entries = 6
    };
    alea_raycast_batch_result_t* tracks = alea_raycast_batch_result_create();
    ASSERT_NOT_NULL(tracks);
    ASSERT_EQ(alea_raycast_hier_batch(sys, origins, directions, 2, 1.5,
                                      &options, tracks), 0);
    const uint64_t* rows = alea_raycast_batch_ray_offsets(tracks);
    const uint64_t* paths = alea_raycast_batch_segment_path_offsets(tracks);
    const uint64_t* keys = alea_raycast_batch_path_occurrence_keys(tracks);
    const int32_t* cells = alea_raycast_batch_cell_ids(tracks);
    const double* enter = alea_raycast_batch_t_enter(tracks);
    const double* exit = alea_raycast_batch_t_exit(tracks);
    ASSERT_EQ(rows[1], 3);
    ASSERT_EQ(rows[2], 6);
    for (size_t i = 0; i < 3; ++i) {
        ASSERT_EQ(cells[i], (int)i+1);
        ASSERT_EQ(cells[5-i], (int)i+1);
        ASSERT_EQ(paths[i+1] - paths[i], 1);
        ASSERT_EQ(paths[6-i] - paths[5-i], 1);
        ASSERT(keys[paths[i]] != 0);
        ASSERT_EQ(keys[paths[i]], keys[paths[5-i]]);
    }
    ASSERT_NEAR(exit[1]-enter[1], 0.00001, 1e-12);
    ASSERT_NEAR(exit[4]-enter[4], 0.00001, 1e-12);
    ASSERT_NEAR(exit[2], 1.5, 1e-12);
    ASSERT_NEAR(exit[5], 1.5, 1e-12);
    alea_raycast_batch_result_destroy(tracks);
    alea_destroy(sys);
}

TEST(ownership_tracks_match_surface_tracks_for_lattice_paths) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_lattice_repeating.mcnp");
    ASSERT_NOT_NULL(model);
    const double origins[] = {-4, 0.1, 0.1, 4, 0.1, 0.1, 0.1, -4, 0.1};
    const double directions[] = {1, 0, 0, -1, 0, 0, 0, 1, 0};
    alea_raycast_batch_result_t* lean = alea_raycast_batch_result_create();
    alea_raycast_batch_result_t* rich = alea_raycast_batch_result_create();
    ASSERT_NOT_NULL(lean);
    ASSERT_NOT_NULL(rich);
    for (int paths = 0; paths < 2; ++paths) {
        alea_raycast_batch_options_t options = {
            .struct_size = sizeof(options),
            .fields = ALEA_RAY_BATCH_MATERIAL | ALEA_RAY_BATCH_DENSITY |
                      ALEA_RAY_BATCH_RESOLUTION_FLAGS,
            .projected_depth = -1
        };
        if (paths) options.fields |= ALEA_RAY_BATCH_FULL_PATHS |
                                    ALEA_RAY_BATCH_PROJECTED_OWNER;
        ASSERT_EQ(alea_raycast_hier_batch(model->sys, origins, directions,
                                          3, 8, &options, lean), 0);
        options.fields |= ALEA_RAY_BATCH_SURFACES;
        ASSERT_EQ(alea_raycast_hier_batch(model->sys, origins, directions,
                                          3, 8, &options, rich), 0);
        size_t count = alea_raycast_batch_segment_count(lean);
        ASSERT_EQ(count, alea_raycast_batch_segment_count(rich));
#define SAME_ARRAY(accessor, length) \
        ASSERT(memcmp(accessor(lean), accessor(rich), \
                      (length) * sizeof(*accessor(lean))) == 0)
        SAME_ARRAY(alea_raycast_batch_ray_offsets, 4);
        SAME_ARRAY(alea_raycast_batch_t_enter, count);
        SAME_ARRAY(alea_raycast_batch_t_exit, count);
        SAME_ARRAY(alea_raycast_batch_cell_ids, count);
        SAME_ARRAY(alea_raycast_batch_material_ids, count);
        SAME_ARRAY(alea_raycast_batch_densities, count);
        SAME_ARRAY(alea_raycast_batch_resolution_flags, count);
        if (paths) {
            SAME_ARRAY(alea_raycast_batch_segment_path_offsets, count+1);
            const uint64_t entries = alea_raycast_batch_segment_path_offsets(lean)[count];
            SAME_ARRAY(alea_raycast_batch_path_occurrence_keys, entries);
            SAME_ARRAY(alea_raycast_batch_projected_occurrence_keys, count);
        }
#undef SAME_ARRAY
        alea_raycast_batch_work_stats_t a, b;
        ASSERT_EQ(alea_raycast_batch_result_get_work_stats_internal(lean, &a), 0);
        ASSERT_EQ(alea_raycast_batch_result_get_work_stats_internal(rich, &b), 0);
        ASSERT_EQ(a.max_boundary_event_enrichments, 0);
        ASSERT(b.max_boundary_event_enrichments > 0);
    }
    alea_raycast_batch_result_destroy(lean);
    alea_raycast_batch_result_destroy(rich);
    mcnp_model_destroy(model);
}

TEST(deep_translated_tracks_and_navigation_match_analytic_world_interface) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_deep_translated_fills.mcnp");
    ASSERT_NOT_NULL(model);
    const double origins[] = {-4, 0, 0, 6, 0, 0};
    const double directions[] = {1, 0, 0, -1, 0, 0};
    alea_raycast_batch_options_t options = {
        .struct_size = sizeof(options), .fields = ALEA_RAY_BATCH_FULL_PATHS |
                                                  ALEA_RAY_BATCH_SURFACES,
        .projected_depth = -1
    };
    alea_raycast_batch_result_t* tracks = alea_raycast_batch_result_create();
    ASSERT_NOT_NULL(tracks);
    ASSERT_EQ(alea_raycast_hier_batch(model->sys, origins, directions, 2, 10,
                                      &options, tracks), 0);
    ASSERT_EQ(alea_raycast_batch_segment_count(tracks), 4);
    const int32_t* cells = alea_raycast_batch_cell_ids(tracks);
    const double* enter = alea_raycast_batch_t_enter(tracks);
    const double* exit = alea_raycast_batch_t_exit(tracks);
    const uint64_t* paths = alea_raycast_batch_segment_path_offsets(tracks);
    const uint64_t* keys = alea_raycast_batch_path_occurrence_keys(tracks);
    const int expected_cells[] = {1, 2, 2, 1};
    const double lengths[] = {6, 4, 4, 6};
    for (size_t i = 0; i < 4; ++i) {
        ASSERT_EQ(cells[i], expected_cells[i]);
        ASSERT_NEAR(exit[i]-enter[i], lengths[i], 1e-12);
        ASSERT_EQ(paths[i+1]-paths[i], 9);
    }
    ASSERT_EQ(keys[paths[1]-1], keys[paths[4]-1]);
    ASSERT_EQ(keys[paths[2]-1], keys[paths[3]-1]);
    alea_ray_navigator_t* nav = alea_ray_navigator_create(model->sys);
    ASSERT_NOT_NULL(nav);
    const double start[] = {1, 0, 0};
    alea_nav_location_t location;
    ASSERT_EQ(alea_ray_navigator_restart(nav, start, directions, &location), 0);
    ASSERT_EQ(location.kind, ALEA_NAV_MATERIAL);
    ASSERT_EQ(location.cell_id, 1);
    ASSERT_EQ(location.occurrence_key, keys[paths[1]-1]);
    alea_nav_event_t event;
    ASSERT_EQ(alea_ray_navigator_advance(nav, INFINITY, 10, &event), 0);
    ASSERT_EQ(event.kind, ALEA_NAV_BOUNDARY);
    ASSERT_NEAR(event.distance, 1, 1e-12);
    ASSERT_NEAR(event.position[0], 2, 1e-12);
    ASSERT_EQ(event.after.cell_id, 2);
    ASSERT_EQ(event.surface_id, 1);
    ASSERT_EQ(alea_ray_navigator_set_direction(nav, directions+3), 0);
    ASSERT_EQ(alea_ray_navigator_advance(nav, 0.125, 1, &event), 0);
    ASSERT_EQ(event.kind, ALEA_NAV_COLLISION);
    ASSERT_EQ(event.before.cell_id, 1);
    ASSERT_NEAR(event.position[0], 1.875, 1e-12);
    alea_ray_navigator_destroy(nav);
    alea_raycast_batch_result_destroy(tracks);
    mcnp_model_destroy(model);
}

TEST(navigation_preserves_short_flights_before_translated_interface) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_deep_translated_fills.mcnp");
    ASSERT_NOT_NULL(model);
    alea_ray_navigator_t* nav = alea_ray_navigator_create(model->sys);
    ASSERT_NOT_NULL(nav);
    /* The analytic world interface is x=2. Both starts and the midpoint of
     * each short flight are representable, but the forward location sample
     * crosses the plane. Preserve that short flight and its incident owner. */
    const double delta = 0x1p-49;
    for (int mode = ALEA_NAV_VALIDATE_STRICT; mode <= ALEA_NAV_VALIDATE_INTERVAL; ++mode) {
        ASSERT_EQ(alea_ray_navigator_set_validation_mode(nav, mode), 0);
        for (int sign = -1; sign <= 1; sign += 2) {
            double position[] = {2 - sign*delta, 0, 0};
            double direction[] = {sign, 0, 0};
            const int incident = sign > 0 ? 1 : 2;
            const int outgoing = sign > 0 ? 2 : 1;
            alea_nav_location_t location;
            ASSERT_EQ(alea_ray_navigator_restart(nav, position, direction, &location), 0);
            ASSERT_EQ(location.kind, ALEA_NAV_MATERIAL);
            ASSERT_EQ(location.cell_id, incident);
            /* Re-locating after a direction change must preserve the same
             * pre-crossing owner, just as a fresh restart does. */
            direction[0] = -sign;
            ASSERT_EQ(alea_ray_navigator_set_direction(nav, direction), 0);
            direction[0] = sign;
            ASSERT_EQ(alea_ray_navigator_set_direction(nav, direction), 0);
            alea_nav_event_t event;
            ASSERT_EQ(alea_ray_navigator_advance(nav, INFINITY, 1, &event), 0);
            ASSERT_EQ(event.kind, ALEA_NAV_BOUNDARY);
            ASSERT_EQ(event.before.cell_id, incident);
            ASSERT_EQ(event.after.cell_id, outgoing);
            ASSERT_EQ(event.surface_id, 1);
            ASSERT(event.distance > 0);
            ASSERT_NEAR(event.distance / delta, 1, 1e-12);
            ASSERT_EQ(event.position[0], 2);
            ASSERT_EQ(alea_ray_navigator_advance(nav, 0.125, 1, &event), 0);
            ASSERT_EQ(event.kind, ALEA_NAV_COLLISION);
            ASSERT_EQ(event.before.cell_id, outgoing);
            ASSERT_NEAR(event.position[0], 2 + sign*0.125, 1e-12);
        }
    }
    alea_ray_navigator_destroy(nav);
    mcnp_model_destroy(model);
}

TEST(navigation_restarts_near_transformed_lattice_cylinder_exit) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_nested_lattice_transformed.mcnp");
    ASSERT_NOT_NULL(model);
    alea_ray_navigator_t* nav = alea_ray_navigator_create(model->sys);
    ASSERT_NOT_NULL(nav);
    /* Public workload reproducer: two placements of the same cylinder.
     * Require successful boundary handling and subsequent outward flight;
     * the translated-plane test above supplies the analytic length oracle. */
    const double starts[] = {6.5318492462042892, 4.5318478156928146};
    const double direction[] = {-1, -0.25, 0.0625};
    for (size_t i = 0; i < 2; ++i) {
        double position[] = {starts[i], -0.69744837142309479, 0.030286056323686663};
        alea_nav_location_t location;
        ASSERT_EQ(alea_ray_navigator_restart(nav, position, direction, &location), 0);
        ASSERT_EQ(location.kind, ALEA_NAV_MATERIAL);
        alea_nav_event_t event;
        ASSERT_EQ(alea_ray_navigator_advance(nav, 0.01, 0.1, &event), 0);
        if (event.kind == ALEA_NAV_BOUNDARY) {
            ASSERT_EQ(event.surface_id, 10);
            ASSERT(event.distance < 1e-12);
            ASSERT_EQ(event.after.cell_id, 2);
            ASSERT_EQ(alea_ray_navigator_advance(nav, 0.01, 0.1, &event), 0);
        }
        ASSERT_EQ(event.kind, ALEA_NAV_COLLISION);
        ASSERT_EQ(event.before.cell_id, 2);
    }
    alea_ray_navigator_destroy(nav);
    mcnp_model_destroy(model);
}

TEST(nearest_recovery_preserves_translated_aliases_and_close_surfaces) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_deep_translated_fills.mcnp");
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    /* Independently specified stress geometry: shared Boolean subexpressions
     * repeat two coincident cards and a distinct nearby plane. Eight existing
     * translations place the planes at world x=2 and x=2.00001. */
    int alias = alea_plane_surface(sys, 99, 1, 0, 0, 0);
    int close = alea_plane_surface(sys, 98, 1, 0, 0, -0.00001);
    ASSERT(alias >= 0 && close >= 0);
    int index = alea_cell_find(sys, 1);
    ASSERT(index >= 0);
    alea_node_id_t region = alea_intersection(sys,
        alea_halfspace(sys, alias, -1),
        alea_intersection(sys, sys->cells.data[index].root_node_id,
                          alea_halfspace(sys, close, -1)));
    for (int i = 0; i < 12; ++i) region = alea_intersection(sys, region, region);
    ASSERT_EQ(alea_cell_set_region(sys, index, region), 0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    alea_cell_entry_t* cell = &sys->cells.data[index];
    ASSERT_EQ(cell->surface_index_count, 3);
    const double origins[] = {1, 3, 2.000005};
    const double directions[] = {1, -1, 1};
    const double distances[] = {1, 0.99999, 0.000005};
    const int surfaces[] = {99, 98, 98};
    for (size_t i = 0; i < 3; ++i) {
        alea_ray_t ray;
        ASSERT_EQ(alea_ray_init(&ray, origins[i], 0, 0, directions[i], 0, 0), 0);
        alea_raycast_result_t indexed, tree;
        alea_raycast_result_init(&indexed);
        alea_raycast_result_init(&tree);
        indexed.breakpoint_nearest_only = tree.breakpoint_nearest_only = 1;
        ASSERT_EQ(alea_raycast_validation_breakpoints_reuse_nocache(
            sys, &ray, 0, 3, 1, &indexed), 0);
        uint32_t* saved = cell->surface_indices;
        cell->surface_indices = NULL;
        uint32_t** lists = detach_recovery_lists(sys);
        const int rc = alea_raycast_validation_breakpoints_reuse_nocache(
            sys, &ray, 0, 3, 1, &tree);
        restore_recovery_lists(sys, lists);
        cell->surface_indices = saved;
        ASSERT_EQ(rc, 0);
        ASSERT_EQ(indexed.hits.count, 1);
        ASSERT_EQ(tree.hits.count, 1);
        ASSERT_NEAR(indexed.hits.data[0].t, distances[i], 1e-12);
        ASSERT_EQ(indexed.hits.data[0].surface_id, surfaces[i]);
        ASSERT_NEAR(indexed.hits.data[0].nx, 1, 1e-12);
        ASSERT_EQ(indexed.hits.data[0].ny, 0);
        ASSERT_EQ(indexed.hits.data[0].nz, 0);
        ASSERT_EQ(indexed.hits.data[0].t, tree.hits.data[0].t);
        ASSERT_EQ(indexed.hits.data[0].surface_id, tree.hits.data[0].surface_id);
        ASSERT_EQ(indexed.hits.data[0].primitive_id, tree.hits.data[0].primitive_id);
        ASSERT_EQ(indexed.hits.data[0].nx, tree.hits.data[0].nx);
        ASSERT_EQ(indexed.hits.data[0].ny, tree.hits.data[0].ny);
        ASSERT_EQ(indexed.hits.data[0].nz, tree.hits.data[0].nz);
        ASSERT(indexed.surfaces_tested * 10 < tree.surfaces_tested);
        alea_raycast_result_free(&indexed);
        alea_raycast_result_free(&tree);
    }
    mcnp_model_destroy(model);
}

TEST(nearest_recovery_matches_tree_fallback_in_repeated_geometry) {
    const char* models[] = {
        "tests/data/mcnp_lattice_repeating.mcnp",
        "tests/data/mcnp_nested_lattice_transformed.mcnp",
        "tests/data/mcnp_deep_translated_fills.mcnp"
    };
    const double starts[] = {-4, 3, -4};
    const double lengths[] = {8, 6, 10};
    for (size_t m = 0; m < 3; ++m) {
        mcnp_model_t* model = mcnp_load(models[m]);
        ASSERT_NOT_NULL(model);
        alea_system_t* sys = model->sys;
        ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
        uint32_t** saved = malloc(sys->cells.count * sizeof(*saved));
        ASSERT_NOT_NULL(saved);
        for (size_t c = 0; c < sys->cells.count; ++c)
            saved[c] = sys->cells.data[c].surface_indices;
        for (int sign = -1; sign <= 1; sign += 2) {
            alea_ray_t ray;
            ASSERT_EQ(alea_ray_init(&ray, starts[m] + (sign < 0 ? lengths[m] : 0),
                                    0.125, 0.0625, sign, 0, 0), 0);
            alea_raycast_result_t indexed, tree;
            alea_raycast_result_init(&indexed);
            alea_raycast_result_init(&tree);
            indexed.breakpoint_nearest_only = tree.breakpoint_nearest_only = 1;
            const int indexed_rc = alea_raycast_validation_breakpoints_reuse_nocache(
                sys, &ray, 0, lengths[m], 1, &indexed);
            for (size_t c = 0; c < sys->cells.count; ++c)
                sys->cells.data[c].surface_indices = NULL;
            uint32_t** lists = detach_recovery_lists(sys);
            const int tree_rc = alea_raycast_validation_breakpoints_reuse_nocache(
                sys, &ray, 0, lengths[m], 1, &tree);
            restore_recovery_lists(sys, lists);
            for (size_t c = 0; c < sys->cells.count; ++c)
                sys->cells.data[c].surface_indices = saved[c];
            ASSERT_EQ(indexed_rc, 0);
            ASSERT_EQ(tree_rc, 0);
            ASSERT_EQ(indexed.hits.count, 1);
            ASSERT_EQ(tree.hits.count, 1);
            ASSERT_EQ(indexed.hits.data[0].t, tree.hits.data[0].t);
            ASSERT_EQ(indexed.hits.data[0].surface_id, tree.hits.data[0].surface_id);
            ASSERT_EQ(indexed.hits.data[0].primitive_id, tree.hits.data[0].primitive_id);
            ASSERT_EQ(indexed.hits.data[0].nx, tree.hits.data[0].nx);
            ASSERT_EQ(indexed.hits.data[0].ny, tree.hits.data[0].ny);
            ASSERT_EQ(indexed.hits.data[0].nz, tree.hits.data[0].nz);
            alea_raycast_result_free(&indexed);
            alea_raycast_result_free(&tree);
        }
        free(saved);
        mcnp_model_destroy(model);
    }
}

TEST(nearest_recovery_shares_solves_within_each_occurrence_only) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_deep_translated_fills.mcnp");
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    int terminal = alea_cell_find(sys, 1);
    ASSERT(terminal >= 0);
    alea_node_id_t region = sys->cells.data[terminal].root_node_id;
    /* Deliberate duplicate diagnostic cells stress mathematical crossing
     * enumeration, rather than define a valid physical ownership model. */
    for (int i = 0; i < 256; ++i)
        ASSERT(alea_add_cell(sys, 1000+i, region, ALEA_MATERIAL_VOID, 0, 8) >= 0);
    int box = alea_box_surface(sys, 999, -1, 5, -1, 1, -1, 1);
    ASSERT(box >= 0);
    int parent = alea_add_cell(sys, 999, alea_halfspace(sys, box, -1),
                               ALEA_MATERIAL_VOID, 0, 0);
    const double translation[] = {1, 0, 0};
    ASSERT_EQ(alea_add_transform(sys, 999, translation, 3, 0), 0);
    ASSERT_EQ(alea_set_fill(sys, parent, 8, 999), 0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    /* The original occurrence has its plane at x=2. The later occurrence
     * translates the same cards to x=1 and must replace the farther hit. */
    alea_ray_t ray;
    ASSERT_EQ(alea_ray_init(&ray, 0.5, 0, 0, 1, 0, 0), 0);
    alea_raycast_result_t indexed, fallback;
    alea_raycast_result_init(&indexed);
    alea_raycast_result_init(&fallback);
    indexed.breakpoint_nearest_only = fallback.breakpoint_nearest_only = 1;
    ASSERT_EQ(alea_raycast_validation_breakpoints_reuse_nocache(
        sys, &ray, 0, 3, 1, &indexed), 0);
    uint32_t** saved = malloc(sys->cells.count * sizeof(*saved));
    ASSERT_NOT_NULL(saved);
    for (size_t c = 0; c < sys->cells.count; ++c) {
        saved[c] = sys->cells.data[c].surface_indices;
        sys->cells.data[c].surface_indices = NULL;
    }
    uint32_t** lists = detach_recovery_lists(sys);
    int rc = alea_raycast_validation_breakpoints_reuse_nocache(
        sys, &ray, 0, 3, 1, &fallback);
    restore_recovery_lists(sys, lists);
    for (size_t c = 0; c < sys->cells.count; ++c)
        sys->cells.data[c].surface_indices = saved[c];
    free(saved);
    ASSERT_EQ(rc, 0);
    ASSERT_EQ(indexed.hits.count, 1);
    ASSERT_EQ(fallback.hits.count, 1);
    ASSERT_EQ(indexed.hits.data[0].t, 0.5);
    ASSERT_EQ(indexed.hits.data[0].surface_id, 1);
    ASSERT_EQ(indexed.hits.data[0].nx, 1);
    ASSERT_EQ(indexed.hits.data[0].t, fallback.hits.data[0].t);
    ASSERT_EQ(indexed.hits.data[0].surface_id, fallback.hits.data[0].surface_id);
    ASSERT(indexed.surfaces_tested * 3 < fallback.surfaces_tested);
    alea_raycast_result_free(&indexed);
    alea_raycast_result_free(&fallback);
    mcnp_model_destroy(model);
}

TEST(recovery_surface_cache_survives_clone_rebuild_and_geometry_edits) {
    mcnp_model_t* model = mcnp_load("tests/data/mcnp_deep_translated_fills.mcnp");
    ASSERT_NOT_NULL(model);
    alea_system_t* sys = model->sys;
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    const alea_universe_t* univ = alea_get_universe(sys, 8);
    ASSERT_NOT_NULL(univ);
    ASSERT_NOT_NULL(univ->recovery_surfaces);
    ASSERT_EQ(univ->recovery_surface_count, 1);
    alea_system_t* clone = alea_clone(sys);
    ASSERT_NOT_NULL(clone);
    ASSERT_EQ(alea_prepare_query_acceleration(clone), 0);
    const alea_universe_t* copied = alea_get_universe(clone, 8);
    ASSERT_NOT_NULL(copied->recovery_surfaces);
    ASSERT(copied->recovery_surfaces != univ->recovery_surfaces);
    ASSERT_EQ(alea_build_cell_surface_index(sys), 0);
    ASSERT(!alea_system_query_cache_ready(sys, ALEA_CACHE_RECOVERY_SURFACES));
    ASSERT(univ->recovery_surfaces == NULL);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    ASSERT_NOT_NULL(alea_get_universe(sys, 8)->recovery_surfaces);
    int surface = alea_plane_surface(sys, 999, 1, 0, 0, -0.25);
    ASSERT(surface >= 0);
    ASSERT_EQ(alea_cell_set_region(sys, alea_cell_find(sys, 1),
                                  alea_halfspace(sys, surface, -1)), 0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys), 0);
    ASSERT_EQ(alea_get_universe(sys, 8)->recovery_surface_count, 2);
    ASSERT_EQ(alea_get_universe(clone, 8)->recovery_surface_count, 1);
    /* Destroying either cache must leave the other usable. */
    mcnp_model_destroy(model);
    alea_ray_t ray;
    ASSERT_EQ(alea_ray_init(&ray, 1, 0, 0, 1, 0, 0), 0);
    alea_raycast_result_t result;
    alea_raycast_result_init(&result);
    result.breakpoint_nearest_only = 1;
    ASSERT_EQ(alea_raycast_validation_breakpoints_reuse_nocache(
        clone, &ray, 0, 3, 1, &result), 0);
    ASSERT_EQ(result.hits.count, 1);
    ASSERT_EQ(result.hits.data[0].t, 1);
    alea_raycast_result_free(&result);
    alea_destroy(clone);
}

/* Independent analytic torus cases; the rejector encloses the complete
 * surface, including roots outside any Boolean cell region. */
TEST(torus_window_certificate_preserves_crossings_and_boundary_cases) {
    for (int axis = ALEA_AXIS_X; axis <= ALEA_AXIS_Z; axis++) {
        int radial = (axis + 1) % 3, other = (axis + 2) % 3;
        alea_torus_data_t torus = {.axis=axis, .major_radius=5,
                                   .minor_radius=1, .axial_semiwidth_B=1};
        double origin[3]={0}, direction[3]={0};
        origin[radial]=-10; direction[radial]=1;
        alea_ray_t ray;
        ASSERT_EQ(alea_ray_init(&ray,origin[0],origin[1],origin[2],
                               direction[0],direction[1],direction[2]),0);
        double roots[4];
        ASSERT_EQ(ray_intersect_torus(&ray,&torus,roots),4);
        const double expected[]={4,6,14,16};
        for (int j=0;j<4;j++) {
            ASSERT_NEAR(roots[j],expected[j],1e-8);
            ASSERT(ray_torus_window_may_hit(&ray,&torus,
                                            expected[j]-1e-6,expected[j]+1e-6));
        }
        ASSERT(!ray_torus_window_may_hit(&ray,&torus,0,2));
        ASSERT(ray_torus_window_may_hit(&ray,&torus,nextafter(4,0),nextafter(4,INFINITY)));
        origin[other]=6; /* Outer tangent at radial coordinate zero. */
        ASSERT_EQ(alea_ray_init(&ray,origin[0],origin[1],origin[2],
                               direction[0],direction[1],direction[2]),0);
        ASSERT(ray_torus_window_may_hit(&ray,&torus,0,20));
        origin[other]=nextafter(6,INFINITY);
        ASSERT_EQ(alea_ray_init(&ray,origin[0],origin[1],origin[2],
                               direction[0],direction[1],direction[2]),0);
        ASSERT(ray_torus_window_may_hit(&ray,&torus,0,20));
        origin[other]=6.001;
        ASSERT_EQ(alea_ray_init(&ray,origin[0],origin[1],origin[2],
                               direction[0],direction[1],direction[2]),0);
        ASSERT(!ray_torus_window_may_hit(&ray,&torus,0,20));
        origin[other]=0; origin[axis]=2;
        ASSERT_EQ(alea_ray_init(&ray,origin[0],origin[1],origin[2],
                               direction[0],direction[1],direction[2]),0);
        ASSERT(!ray_torus_window_may_hit(&ray,&torus,0,20));
        /* Keep the elliptical surface when the axial width differs. */
        torus.axial_semiwidth_B=3;
        ASSERT(ray_torus_window_may_hit(&ray,&torus,0,20));
        torus.axial_semiwidth_B=.25; origin[axis]=.75;
        ASSERT_EQ(alea_ray_init(&ray,origin[0],origin[1],origin[2],
                               direction[0],direction[1],direction[2]),0);
        ASSERT(!ray_torus_window_may_hit(&ray,&torus,0,20));
    }
}

TEST(torus_window_certificate_handles_tiny_directions_large_coordinates_and_fallback) {
    alea_torus_data_t torus={.axis=ALEA_AXIS_Z,.major_radius=5,
                             .minor_radius=1,.axial_semiwidth_B=1};
    alea_ray_t ray;
    ASSERT_EQ(alea_ray_init(&ray,5,0,2,0,0,-1),0);
    ray.dz=-1e-20; /* Algebraic non-unit ray: never treat tiny as zero. */
    ASSERT(ray_torus_window_may_hit(&ray,&torus,1e20,3e20));
    torus.center_x=1e12; torus.center_y=-1e12; torus.center_z=1e12;
    ASSERT_EQ(alea_ray_init(&ray,1e12-10,-1e12,1e12,1,0,0),0);
    ASSERT(ray_torus_window_may_hit(&ray,&torus,3.999,4.001));
    ASSERT(!ray_torus_window_may_hit(&ray,&torus,0,1));
    /* Degenerate and non-finite inputs retain solver fallback. */
    ASSERT_EQ(alea_ray_init(&ray,0,0,2,1,0,0),0);
    torus.center_x=torus.center_y=torus.center_z=0;
    torus.axial_semiwidth_B=0;
    ASSERT(ray_torus_window_may_hit(&ray,&torus,0,20));
    torus.axial_semiwidth_B=1; torus.major_radius=0;
    ASSERT(ray_torus_window_may_hit(&ray,&torus,0,20));
    torus.major_radius=5; torus.minor_radius=NAN;
    ASSERT(ray_torus_window_may_hit(&ray,&torus,0,20));
    torus.minor_radius=1;
    ASSERT(ray_torus_window_may_hit(&ray,&torus,0,INFINITY));
    torus.minor_radius=torus.axial_semiwidth_B=1e-200;
    ASSERT_EQ(alea_ray_init(&ray,-10,0,1e-10,1,0,0),0);
    /* Existing quartic can publish roots after r^2 underflows: retain fallback. */
    double thin_roots[4];
    ASSERT(ray_intersect_torus(&ray,&torus,thin_roots)>0);
    ASSERT(ray_torus_window_may_hit(&ray,&torus,0,20));
    torus.minor_radius=torus.axial_semiwidth_B=1;
    torus.major_radius=DBL_MAX; torus.minor_radius=DBL_MAX/2;
    ASSERT(ray_torus_window_may_hit(&ray,&torus,0,20));
    torus.major_radius=5; torus.minor_radius=1;
    ray.ox=INFINITY;
    ASSERT(ray_torus_window_may_hit(&ray,&torus,0,20));
}

TEST(nearest_torus_recovery_matches_unpruned_tree_across_transformed_occurrences) {
    alea_system_t* sys=alea_create();
    ASSERT_NOT_NULL(sys);
    int torus=alea_torus_z_surface(sys,1,0,0,0,5,1);
    ASSERT(torus>=0);
    ASSERT(alea_add_cell(sys,1,alea_halfspace(sys,torus,-1),ALEA_MATERIAL_VOID,0,1)>=0);
    ASSERT(alea_add_cell(sys,2,alea_halfspace(sys,torus,1),ALEA_MATERIAL_VOID,0,1)>=0);
    int box=alea_box_surface(sys,2,-100,100,-100,100,-100,100);
    ASSERT(box>=0);
    /* Overlapping diagnostic parents test geometric enumeration only. */
    for (int i=0;i<2;i++) {
        int parent=alea_add_cell(sys,10+i,alea_halfspace(sys,box,-1),ALEA_MATERIAL_VOID,0,0);
        const double transform[]={i*20.0,0,0, 1,0,0, 0,0,-1, 0,1,0};
        ASSERT_EQ(alea_add_transform(sys,10+i,transform,12,0),0);
        ASSERT_EQ(alea_set_fill(sys,parent,1,10+i),0);
    }
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    for (int reverse=0;reverse<2;reverse++) for (int miss=0;miss<2;miss++) {
        alea_ray_t ray;
        ASSERT_EQ(alea_ray_init(&ray,reverse?30:-10,miss?2:0,miss?2:0,reverse?-1:1,0,0),0);
        alea_raycast_result_t indexed,tree;
        alea_raycast_result_init(&indexed); alea_raycast_result_init(&tree);
        indexed.breakpoint_nearest_only=tree.breakpoint_nearest_only=1;
        ASSERT_EQ(alea_raycast_validation_breakpoints_reuse_nocache(sys,&ray,0,40,1,&indexed),0);
        uint32_t** lists=detach_recovery_lists(sys);
        uint32_t* saved[2]={sys->cells.data[0].surface_indices,sys->cells.data[1].surface_indices};
        sys->cells.data[0].surface_indices=sys->cells.data[1].surface_indices=NULL;
        int rc=alea_raycast_validation_breakpoints_reuse_nocache(sys,&ray,0,40,1,&tree);
        sys->cells.data[0].surface_indices=saved[0];sys->cells.data[1].surface_indices=saved[1];
        restore_recovery_lists(sys,lists);
        ASSERT_EQ(rc,0);ASSERT_EQ(indexed.hits.count,tree.hits.count);
        ASSERT_EQ(indexed.hits.count,miss?0:1);
        if (!miss) {
            ASSERT_NEAR(indexed.hits.data[0].t,4,1e-8);
            ASSERT_EQ(indexed.hits.data[0].t,tree.hits.data[0].t);
            ASSERT_EQ(indexed.hits.data[0].surface_id,tree.hits.data[0].surface_id);
            ASSERT_EQ(indexed.hits.data[0].nx,tree.hits.data[0].nx);
            ASSERT_EQ(indexed.hits.data[0].ny,tree.hits.data[0].ny);
            ASSERT_EQ(indexed.hits.data[0].nz,tree.hits.data[0].nz);
        } else ASSERT(indexed.surfaces_tested<tree.surfaces_tested);
        alea_raycast_result_free(&indexed);alea_raycast_result_free(&tree);
    }
    alea_destroy(sys);
}

TEST(torus_window_rejection_preserves_solver_roots_in_deterministic_sweep) {
    uint64_t seed=123456789;
    const double scales[]={1e-9,1e-6,1e-3,1,1e3,1e6};
    for (int axis=ALEA_AXIS_X;axis<=ALEA_AXIS_Z;axis++)
        for (size_t scale=0;scale<6;scale++) for (int i=0;i<512;i++) {
            double values[6];
            for (int j=0;j<6;j++) {
                seed=seed*UINT64_C(6364136223846793005)+1;
                values[j]=2.0*(double)(seed>>11)/9007199254740992.0-1;
            }
            double size=scales[scale];
            alea_torus_data_t torus={.axis=axis,.major_radius=5*size,
                .minor_radius=size,.axial_semiwidth_B=(i&1?0.5:2)*size};
            alea_ray_t ray;
            ASSERT_EQ(alea_ray_init(&ray,values[0]*12*size,values[1]*12*size,
                values[2]*12*size,values[3],values[4],values[5]),0);
            double roots[4];int count=ray_intersect_torus(&ray,&torus,roots);
            double end=30*size;
            if (!ray_torus_window_may_hit(&ray,&torus,0,end))
                for (int j=0;j<count;j++) ASSERT(!(roots[j]>=0 && roots[j]<=end));
            for (int j=0;j<count;j++) if (roots[j]>=0 && roots[j]<=end)
                ASSERT(ray_torus_window_may_hit(&ray,&torus,
                    nextafter(roots[j],0),nextafter(roots[j],INFINITY)));
        }
}

TEST(recovery_context_requires_exact_keys_and_contained_windows) {
    recovery_context_t* table = recovery_context_create(1);
    ASSERT_NOT_NULL(table);
    ASSERT(recovery_context_create(0)==NULL);
    ASSERT(recovery_context_create(8192)==NULL);
    uint32_t cards[]={1,2};
    alea_ray_t ray;
    ASSERT_EQ(alea_ray_init(&ray,0,0,0,1,0,0),0);
    alea_matrix_t matrix;
    alea_matrix_identity(&matrix);
    recovery_context_key_t key, changed;
    ASSERT(recovery_context_key(&key,cards,2,&ray,&matrix));
    recovery_context_entry_t* slot=recovery_context_find(table,&key);
    ASSERT_NOT_NULL(slot);
    ASSERT(!recovery_context_contains(slot,1,2)); /* Not yet completed. */
    slot->key=key; slot->lo=1; slot->hi=3; slot->complete=true;
    ASSERT(recovery_context_contains(slot,1,3));
    ASSERT(recovery_context_contains(slot,2,2));
    ASSERT(!recovery_context_contains(slot,nextafter(1,0),2));
    ASSERT(!recovery_context_contains(slot,2,nextafter(3,INFINITY)));
    ASSERT(!recovery_context_contains(slot,4,5));
    ASSERT(!recovery_context_contains(slot,NAN,2));
    ASSERT(!recovery_context_contains(slot,2,INFINITY));
    changed=key; changed.count=1;
    ASSERT(recovery_context_find(table,&changed)==NULL); /* Collision + full. */
    changed=key; changed.values[0]=nextafter(0,1);
    ASSERT(recovery_context_find(table,&changed)==NULL);
    matrix.m[3]=1;
    ASSERT(recovery_context_key(&changed,cards,2,&ray,&matrix));
    ASSERT(recovery_context_find(table,&changed)==NULL);
    ray.ox=NAN;
    ASSERT(!recovery_context_key(&changed,cards,2,&ray,&matrix));
    ASSERT(recovery_context_find(NULL,&key)==NULL);
    free(table);
}

TEST(recovery_context_reuses_equal_parents_preserving_aliases_and_nearer_transforms) {
    alea_system_t* sys=alea_create();
    ASSERT_NOT_NULL(sys);
    int plane=alea_plane_surface(sys,1,1,0,0,-2);
    int alias=alea_plane_surface(sys,99,1,0,0,-2);
    ASSERT(plane>=0 && alias>=0);
    ASSERT(alea_add_cell(sys,1,alea_halfspace(sys,plane,-1),ALEA_MATERIAL_VOID,0,8)>=0);
    ASSERT(alea_add_cell(sys,2,alea_halfspace(sys,alias,1),ALEA_MATERIAL_VOID,0,8)>=0);
    int box=alea_box_surface(sys,3,-10,10,-2,2,-2,2);
    for(int i=0;i<32;i++) {
        int parent=alea_add_cell(sys,100+i,alea_halfspace(sys,box,-1),ALEA_MATERIAL_VOID,0,0);
        ASSERT(parent>=0);
        ASSERT_EQ(alea_set_fill(sys,parent,8,0),0);
    }
    const double translation[]={-1,0,0};
    ASSERT_EQ(alea_add_transform(sys,9,translation,3,0),0);
    int parent=alea_add_cell(sys,200,alea_halfspace(sys,box,-1),ALEA_MATERIAL_VOID,0,0);
    ASSERT_EQ(alea_set_fill(sys,parent,8,9),0);
    ASSERT_EQ(alea_prepare_query_acceleration(sys),0);
    for(int pass=0;pass<4;pass++) {
        alea_ray_t ray;
        ASSERT_EQ(alea_ray_init(&ray,pass==1 ? 4 : 0.5,0,0,pass==1 ? -1 : 1,0,0),0);
        alea_raycast_result_t cached, plain;
        alea_raycast_result_init(&cached); alea_raycast_result_init(&plain);
        cached.breakpoint_nearest_only=plain.breakpoint_nearest_only=1;
        if(pass==3) cached.recovery_context_warmup=16;
        else {
            cached.recovery_context=recovery_context_create(pass==2 ? 1 : 2048);
            ASSERT_NOT_NULL(cached.recovery_context);
        }
        if(pass==0) {
            alea_interrupt();
            ASSERT(alea_raycast_validation_breakpoints_reuse_nocache(
                sys,&ray,0,5,1,&cached)!=0);
            alea_clear_interrupt();
        }
        ASSERT_EQ(alea_raycast_validation_breakpoints_reuse_nocache(sys,&ray,0,5,1,&cached),0);
        ASSERT_EQ(alea_raycast_validation_breakpoints_reuse_nocache(sys,&ray,0,5,1,&plain),0);
        ASSERT_NOT_NULL(cached.recovery_context);
        ASSERT_EQ(cached.hits.count,1); ASSERT_EQ(plain.hits.count,1);
        ASSERT_EQ(cached.hits.data[0].t,pass==1 ? 2 : 0.5);
        ASSERT_EQ(cached.hits.data[0].surface_id,1);
        ASSERT_EQ(cached.hits.data[0].nx,1);
        ASSERT_EQ(cached.hits.data[0].t,plain.hits.data[0].t);
        ASSERT_EQ(cached.hits.data[0].surface_id,plain.hits.data[0].surface_id);
        ASSERT(cached.surfaces_tested<plain.surfaces_tested);
        free(cached.recovery_context); cached.recovery_context=NULL;
        alea_raycast_result_free(&cached); alea_raycast_result_free(&plain);
    }
    alea_destroy(sys);
}

TEST(recovery_context_preserves_elliptical_hits_in_independent_queries) {
    mcnp_model_t* model=mcnp_load("tests/data/mcnp_elliptical_torus.mcnp");
    ASSERT_NOT_NULL(model);
    ASSERT_EQ(alea_prepare_query_acceleration(model->sys),0);
    for(int pass=0;pass<3;pass++) {
        alea_ray_t ray;
        ASSERT_EQ(alea_ray_init(&ray,5,0,pass==1 ? 3 : -3,0,0,pass==1 ? -1 : 1),0);
        alea_raycast_result_t result;
        alea_raycast_result_init(&result);
        result.breakpoint_nearest_only=1;
        result.recovery_context=recovery_context_create(2048);
        ASSERT_NOT_NULL(result.recovery_context);
        ASSERT_EQ(alea_raycast_validation_breakpoints_reuse_nocache(
            model->sys,&ray,0,6,1,&result),0);
        ASSERT_EQ(result.hits.count,1);
        ASSERT_NEAR(result.hits.data[0].t,1,1e-12);
        ASSERT_NEAR(result.hits.data[0].nz,pass==1 ? 1 : -1,1e-12);
        free(result.recovery_context); result.recovery_context=NULL;
        alea_raycast_result_free(&result);
    }
    mcnp_model_destroy(model);
}

TEST_MAIN()
