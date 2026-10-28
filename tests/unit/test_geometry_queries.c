// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#include "alea_test.h"
#include "alea.h"
#include "alea_mcnp.h"
#include "alea_raycast.h"

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

TEST_MAIN()
