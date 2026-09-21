// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_test.h"
#include "alea_tally.h"
#include "transport/tally_internal.h"

TEST(tallies_group_terminal_universes_and_commit_source_history_moments) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    int sphere = alea_sphere_surface(sys, 1, 0, 0, 0, 1);
    ASSERT_TRUE(sphere >= 0);
    alea_node_id_t region = alea_halfspace(sys, sphere, -1);
    ASSERT_EQ(alea_add_cell(sys, 11, region, ALEA_MATERIAL_VOID, 0, 7), 0);
    ASSERT_EQ(alea_add_cell(sys, 12, region, ALEA_MATERIAL_VOID, 0, 7), 1);
    ASSERT_EQ(alea_add_cell(sys, 13, region, ALEA_MATERIAL_VOID, 0, 9), 2);
    alea_tally_plan_t* plan = alea_tally_plan_create(sys);
    ASSERT_NOT_NULL(plan);
    alea_tally_spec_t spec = {
        .score = ALEA_TALLY_TRACK_LENGTH,
        .domain = ALEA_TALLY_UNIVERSE,
        .particle_mask = ALEA_TALLY_NEUTRON,
        .energy_min = 1, .energy_max = 3
    };
    ASSERT_EQ(alea_tally_plan_add(plan, &spec, NULL), ALEA_OK);
    spec.domain = ALEA_TALLY_CARTESIAN_MESH;
    spec.lower[0] = 0; spec.lower[1] = -1; spec.lower[2] = -1;
    spec.upper[0] = 2; spec.upper[1] = 1; spec.upper[2] = 1;
    spec.dimensions[0] = 2;
    spec.dimensions[1] = spec.dimensions[2] = 1;
    ASSERT_EQ(alea_tally_plan_add(plan, &spec, NULL), ALEA_OK);
    spec.score = ALEA_TALLY_REACTION_EVENT;
    spec.domain = ALEA_TALLY_UNIVERSE;
    spec.reaction_mt = 16;
    spec.energy_min = 13; spec.energy_max = 15;
    spec.nuclide_zaid = 1001;
    ASSERT_EQ(alea_tally_plan_add(plan, &spec, NULL), ALEA_OK);
    spec.nuclide_zaid = 3006;
    ASSERT_EQ(alea_tally_plan_add(plan, &spec, NULL), ALEA_OK);
    alea_tally_results_t* results = alea_tally_results_create(plan);
    ASSERT_NOT_NULL(results);
    alea_nav_location_t first = {.kind = ALEA_NAV_VOID,
        .cell_id = 11, .cell_index = 0};
    alea_nav_location_t second = {.kind = ALEA_NAV_VOID,
        .cell_id = 12, .cell_index = 1};
    const double start0[3] = {0, 0, 0};
    const double start1[3] = {0, 0.5, 0};
    const double direction[3] = {1, 0, 0};
    ASSERT_EQ(alea_tally_record_track(results, plan, &first,
        start0, direction, 2, 1, 2, 0, 1, ALEA_TALLY_NEUTRON, NULL), ALEA_OK);
    ASSERT_EQ(alea_tally_record_track(results, plan, &second,
        start1, direction, 1, 1, 2, 0, 1, ALEA_TALLY_NEUTRON, NULL), ALEA_OK);
    ASSERT_EQ(alea_tally_record_collision(results, plan, &first,
        start0, 1, 14, 0, ALEA_TALLY_NEUTRON, 16, 1001), ALEA_OK);
    ASSERT_EQ(alea_tally_record_collision(results, plan, &second,
        start1, 1, 14, 0, ALEA_TALLY_NEUTRON, 16, 1001), ALEA_OK);
    ASSERT_EQ(alea_tally_commit_history(results), ALEA_OK);
    ASSERT_EQ(alea_tally_record_track(results, plan, &first,
        start0, direction, 1, 1, 2, 0, 1, ALEA_TALLY_NEUTRON, NULL), ALEA_OK);
    ASSERT_EQ(alea_tally_record_track(results, plan, &first,
        start0, direction, 1, 1, 4, 0, 1, ALEA_TALLY_NEUTRON, NULL), ALEA_OK);
    ASSERT_EQ(alea_tally_record_collision(results, plan, &first,
        start0, 1, 14, 0, ALEA_TALLY_NEUTRON, 102, 1001), ALEA_OK);
    ASSERT_EQ(alea_tally_commit_history(results), ALEA_OK);
    alea_tally_results_finalize(results);
    alea_tally_view_t view;
    ASSERT_EQ(alea_tally_results_view(results, 0, &view), ALEA_OK);
    ASSERT_EQ(view.histories, 2);
    ASSERT_EQ(view.bin_count, 2);
    ASSERT_EQ(view.bin_ids[0], 7);
    ASSERT_EQ(view.bin_ids[1], 9);
    ASSERT_NEAR(view.sum[0], 4, 1e-12);
    ASSERT_NEAR(view.sum_squared[0], 10, 1e-12);
    ASSERT_NEAR(view.sum[1], 0, 1e-12);
    ASSERT_EQ(alea_tally_results_view(results, 1, &view), ALEA_OK);
    ASSERT_NEAR(view.sum[0], 3, 1e-12);
    ASSERT_NEAR(view.sum_squared[0], 5, 1e-12);
    ASSERT_NEAR(view.sum[1], 1, 1e-12);
    ASSERT_EQ(alea_tally_results_view(results, 2, &view), ALEA_OK);
    ASSERT_NEAR(view.sum[0], 2, 1e-12);
    ASSERT_NEAR(view.sum_squared[0], 4, 1e-12);
    ASSERT_EQ(alea_tally_results_view(results, 3, &view), ALEA_OK);
    ASSERT_NEAR(view.sum[0], 0, 1e-12);
    alea_tally_results_free(results);
    alea_tally_plan_free(plan);
    alea_destroy(sys);
}

TEST(tally_plan_rejects_invalid_bins_and_mutated_geometry) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    alea_tally_plan_t* plan = alea_tally_plan_create(sys);
    ASSERT_NOT_NULL(plan);
    alea_tally_spec_t spec = {
        .score = ALEA_TALLY_TRACK_LENGTH,
        .domain = ALEA_TALLY_CARTESIAN_MESH,
        .lower = {0, 0, 0}, .upper = {1, 1, 1},
        .dimensions = {0, 1, 1}
    };
    ASSERT_EQ(alea_tally_plan_add(plan, &spec, NULL), ALEA_ERR_INVALID_ARG);
    spec.dimensions[0] = 1;
    spec.energy_min = 4; spec.energy_max = 2;
    ASSERT_EQ(alea_tally_plan_add(plan, &spec, NULL), ALEA_ERR_INVALID_ARG);
    spec.energy_min = spec.energy_max = 0;
    ASSERT_EQ(alea_tally_plan_add(plan, &spec, NULL), ALEA_OK);
    spec.time_min = 2; spec.time_max = 1;
    ASSERT_EQ(alea_tally_plan_add(plan, &spec, NULL), ALEA_ERR_INVALID_ARG);
    spec.time_min = spec.time_max = 0;
    spec.score = ALEA_TALLY_REACTION_RATE;
    spec.reaction_mt = -1;
    ASSERT_EQ(alea_tally_plan_add(plan, &spec, NULL), ALEA_ERR_INVALID_ARG);
    spec.reaction_mt = 16;
    ASSERT_EQ(alea_tally_plan_add(plan, &spec, NULL), ALEA_OK);
    double bad_edges[3] = {0, 2, 2};
    spec.score = ALEA_TALLY_TRACK_LENGTH;
    spec.reaction_mt = 0;
    spec.energy_edges = bad_edges;
    spec.energy_group_count = 2;
    ASSERT_EQ(alea_tally_plan_add(plan, &spec, NULL), ALEA_ERR_INVALID_ARG);
    ASSERT_TRUE(alea_sphere_surface(sys, 1, 0, 0, 0, 1) >= 0);
    ASSERT_TRUE(!alea_tally_plan_matches_geometry(plan, sys));
    ASSERT_EQ(alea_tally_plan_add(plan, &spec, NULL), ALEA_ERR_INVALID_STATE);
    alea_tally_plan_free(plan);
    alea_destroy(sys);
}

TEST(mesh_tally_clips_reverse_and_corner_crossing_segments) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    alea_tally_plan_t* plan = alea_tally_plan_create(sys);
    ASSERT_NOT_NULL(plan);
    alea_tally_spec_t spec = {
        .score = ALEA_TALLY_TRACK_LENGTH,
        .domain = ALEA_TALLY_CARTESIAN_MESH,
        .lower = {0, 0, -1}, .upper = {2, 2, 1},
        .dimensions = {2, 2, 1}
    };
    ASSERT_EQ(alea_tally_plan_add(plan, &spec, NULL), ALEA_OK);
    alea_tally_results_t* results = alea_tally_results_create(plan);
    ASSERT_NOT_NULL(results);
    alea_nav_location_t void_location = {
        .kind = ALEA_NAV_VOID, .cell_id = -1, .cell_index = -1
    };
    const double start[3] = {-1, 0.5, 0};
    const double positive[3] = {1, 0, 0};
    const double reverse_start[3] = {2, 0.5, 0};
    const double negative[3] = {-1, 0, 0};
    const double corner_start[3] = {0, 0, 0};
    const double diagonal[3] = {0.7071067811865476, 0.7071067811865476, 0};
    ASSERT_EQ(alea_tally_record_track(results, plan, &void_location,
        start, positive, 4, 1, 2, 0, 1, ALEA_TALLY_NEUTRON, NULL), ALEA_OK);
    ASSERT_EQ(alea_tally_record_track(results, plan, &void_location,
        reverse_start, negative, 2, 1, 2, 0, 1, ALEA_TALLY_NEUTRON, NULL), ALEA_OK);
    ASSERT_EQ(alea_tally_record_track(results, plan, &void_location,
        corner_start, diagonal, 2.8284271247461903,
        1, 2, 0, 1, ALEA_TALLY_NEUTRON, NULL), ALEA_OK);
    ASSERT_EQ(alea_tally_commit_history(results), ALEA_OK);
    alea_tally_view_t view;
    ASSERT_EQ(alea_tally_results_view(results, 0, &view), ALEA_OK);
    ASSERT_NEAR(view.sum[0], 2 + 1.4142135623730951, 1e-10);
    ASSERT_NEAR(view.sum[1], 2, 1e-10);
    ASSERT_NEAR(view.sum[2], 0, 1e-10);
    ASSERT_NEAR(view.sum[3], 1.4142135623730951, 1e-10);
    alea_tally_results_free(results);
    alea_tally_plan_free(plan);
    alea_destroy(sys);
}

TEST(time_filter_clips_flights_and_gates_collision_events) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    alea_tally_plan_t* plan = alea_tally_plan_create(sys);
    ASSERT_NOT_NULL(plan);
    alea_tally_spec_t spec = {
        .score = ALEA_TALLY_TRACK_LENGTH,
        .domain = ALEA_TALLY_CARTESIAN_MESH,
        .lower = {0, -1, -1}, .upper = {4, 1, 1},
        .dimensions = {4, 1, 1},
        .time_min = 0.5, .time_max = 1.5
    };
    ASSERT_EQ(alea_tally_plan_add(plan, &spec, NULL), ALEA_OK);
    spec.score = ALEA_TALLY_COLLISION;
    ASSERT_EQ(alea_tally_plan_add(plan, &spec, NULL), ALEA_OK);
    alea_tally_results_t* results = alea_tally_results_create(plan);
    ASSERT_NOT_NULL(results);
    alea_nav_location_t location = {.kind = ALEA_NAV_VOID,
        .cell_id = -1, .cell_index = -1};
    const double start[3] = {0, 0, 0}, direction[3] = {1, 0, 0};
    const double inside[3] = {2, 0, 0};
    ASSERT_EQ(alea_tally_record_track(results, plan, &location,
        start, direction, 4, 2, 1, 0, 2, ALEA_TALLY_NEUTRON, NULL), ALEA_OK);
    ASSERT_EQ(alea_tally_record_collision(results, plan, &location,
        inside, 1, 1, 0.5, ALEA_TALLY_NEUTRON, 2, 0), ALEA_OK);
    ASSERT_EQ(alea_tally_record_collision(results, plan, &location,
        inside, 1, 1, 1.5, ALEA_TALLY_NEUTRON, 2, 0), ALEA_OK);
    ASSERT_EQ(alea_tally_commit_history(results), ALEA_OK);
    alea_tally_view_t view;
    ASSERT_EQ(alea_tally_results_view(results, 0, &view), ALEA_OK);
    ASSERT_NEAR(view.sum[0], 0, 1e-12);
    ASSERT_NEAR(view.sum[1], 2, 1e-12);
    ASSERT_NEAR(view.sum[2], 2, 1e-12);
    ASSERT_NEAR(view.sum[3], 0, 1e-12);
    ASSERT_EQ(alea_tally_results_view(results, 1, &view), ALEA_OK);
    ASSERT_NEAR(view.sum[2], 1, 1e-12);
    ASSERT_EQ(view.time_min, 0.5);
    ASSERT_EQ(view.time_max, 1.5);
    alea_tally_results_free(results);
    alea_tally_plan_free(plan);
    alea_destroy(sys);
}

TEST(energy_axis_groups_mesh_track_and_collision_scores) {
    alea_system_t* sys = alea_create();
    ASSERT_NOT_NULL(sys);
    alea_tally_plan_t* plan = alea_tally_plan_create(sys);
    ASSERT_NOT_NULL(plan);
    double edges[3] = {0, 1, 3};
    alea_tally_spec_t spec = {
        .score = ALEA_TALLY_TRACK_LENGTH,
        .domain = ALEA_TALLY_CARTESIAN_MESH,
        .lower = {0, -1, -1}, .upper = {2, 1, 1},
        .dimensions = {2, 1, 1},
        .energy_edges = edges, .energy_group_count = 2
    };
    ASSERT_EQ(alea_tally_plan_add(plan, &spec, NULL), ALEA_OK);
    spec.score = ALEA_TALLY_COLLISION;
    ASSERT_EQ(alea_tally_plan_add(plan, &spec, NULL), ALEA_OK);
    edges[1] = 9; /* The plan owns its own edge array. */
    alea_tally_results_t* results = alea_tally_results_create(plan);
    ASSERT_NOT_NULL(results);
    alea_nav_location_t location = {.kind = ALEA_NAV_VOID,
        .cell_id = -1, .cell_index = -1};
    const double start[3] = {0, 0, 0}, direction[3] = {1, 0, 0};
    const double collision[3] = {1.5, 0, 0};
    ASSERT_EQ(alea_tally_record_track(results, plan, &location,
        start, direction, 2, 1, 0.5, 0, 1, ALEA_TALLY_NEUTRON, NULL), ALEA_OK);
    ASSERT_EQ(alea_tally_record_track(results, plan, &location,
        start, direction, 1, 2, 2, 0, 1, ALEA_TALLY_NEUTRON, NULL), ALEA_OK);
    ASSERT_EQ(alea_tally_record_track(results, plan, &location,
        start, direction, 2, 1, 3, 0, 1, ALEA_TALLY_NEUTRON, NULL), ALEA_OK);
    ASSERT_EQ(alea_tally_record_collision(results, plan, &location,
        collision, 1, 1, 0, ALEA_TALLY_NEUTRON, 2, 0), ALEA_OK);
    ASSERT_EQ(alea_tally_record_collision(results, plan, &location,
        collision, 1, 3, 0, ALEA_TALLY_NEUTRON, 2, 0), ALEA_OK);
    ASSERT_EQ(alea_tally_commit_history(results), ALEA_OK);
    alea_tally_view_t view;
    ASSERT_EQ(alea_tally_results_view(results, 0, &view), ALEA_OK);
    ASSERT_EQ(view.spatial_bin_count, 2);
    ASSERT_EQ(view.energy_group_count, 2);
    ASSERT_EQ(view.bin_count, 4);
    ASSERT_NEAR(view.energy_edges[1], 1, 1e-12);
    ASSERT_NEAR(view.sum[0], 1, 1e-12);
    ASSERT_NEAR(view.sum[1], 1, 1e-12);
    ASSERT_NEAR(view.sum[2], 2, 1e-12);
    ASSERT_NEAR(view.sum[3], 0, 1e-12);
    ASSERT_NEAR(view.sum_squared[2], 4, 1e-12);
    ASSERT_EQ(alea_tally_results_view(results, 1, &view), ALEA_OK);
    ASSERT_NEAR(view.sum[3], 1, 1e-12);
    alea_tally_results_free(results);
    alea_tally_plan_free(plan);
    alea_destroy(sys);
}

TEST_MAIN()
