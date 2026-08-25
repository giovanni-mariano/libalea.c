// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/*
 * Compile-only check: every installed header must resolve with -Iinclude and
 * without src/ on the include path. Keep this list in sync with include/.
 */
#include "alea.h"
#include "alea_geo_validator.h"
#include "alea_log.h"
#include "alea_mcnp.h"
#include "alea_mesh.h"
#include "alea_nucdata.h"
#include "alea_nucdata_types.h"
#include "alea_openmc.h"
#include "alea_raycast.h"
#include "alea_render.h"
#include "alea_serpent.h"
#include "alea_slice.h"
#include "alea_types.h"

#ifdef alea_surface_at
#error "alea_surface_at must remain an internal-only helper"
#endif

void public_api_smoke(void) {
    alea_system_t* sys = alea_create();
    int surface = alea_sphere_surface(sys, 0, 0.0, 0.0, 0.0, 1.0);
    alea_node_id_t inside = alea_halfspace(sys, surface, -1);
    alea_node_id_t outside = alea_halfspace(sys, surface, +1);
    alea_raycast_result_t* result = alea_raycast_result_create();
    alea_ray_coverage_slice_result_t* coverage =
        alea_ray_coverage_slice_result_create();
    alea_ray_coverage_slice_options_t coverage_options;
    alea_ray_boundary_event_options_t event_options;
    alea_ray_boundary_event_provenance_t event_provenance;
    alea_ray_boundary_event_query_result_t* event_result =
        alea_ray_boundary_event_query_result_create();
    alea_transition_slice_options_t transition_slice_options;
    alea_transition_slice_result_t* transition_slice =
        alea_transition_slice_result_create();
    alea_transition_slice_finding_t transition_finding;
    alea_transition_slice_component_t transition_component;
    alea_transition_slice_coverage_finding_t coverage_finding;
    alea_transition_slice_coverage_component_t coverage_component;
    alea_transition_slice_component_link_t component_link;
    alea_transition_slice_refinement_frontier_t refinement_frontier;
    alea_transition_slice_critical_tile_t critical_tile;
    alea_transition_slice_critical_tile_source_t critical_tile_source;
    alea_transition_slice_stats_t transition_stats;
    uint8_t flags = 0;
    double hit_t = 0.0;
    int hit_surface_id = -1;

    (void)inside;
    (void)outside;
    (void)alea_raycast_segment_resolution_flags(result, 0, &flags);
    (void)alea_raycast_hit_count(result);
    (void)alea_raycast_hit_get(result, 0, &hit_t, &hit_surface_id);
    alea_ray_coverage_slice_options_init(&coverage_options);
    alea_ray_boundary_event_options_init(&event_options);
    alea_transition_slice_options_init(&transition_slice_options);
    (void)alea_transition_slice_finding_get(
        transition_slice, 0, &transition_finding);
    (void)alea_transition_slice_stats(
        transition_slice, &transition_stats);
    (void)alea_transition_slice_finding_count(transition_slice);
    (void)alea_transition_slice_component_count(transition_slice);
    (void)alea_transition_slice_component_get(
        transition_slice, 0, &transition_component);
    (void)alea_transition_slice_coverage_finding_count(transition_slice);
    (void)alea_transition_slice_coverage_finding_get(
        transition_slice, 0, &coverage_finding);
    (void)alea_transition_slice_coverage_component_count(transition_slice);
    (void)alea_transition_slice_coverage_component_get(
        transition_slice, 0, &coverage_component);
    (void)alea_transition_slice_component_link_count(transition_slice);
    (void)alea_transition_slice_component_link_get(
        transition_slice, 0, &component_link);
    (void)alea_transition_slice_refinement_frontier_count(transition_slice);
    (void)alea_transition_slice_refinement_frontier_get(
        transition_slice, 0, &refinement_frontier);
    (void)alea_transition_slice_critical_tile_count(transition_slice);
    (void)alea_transition_slice_critical_tile_get(
        transition_slice, 0, &critical_tile);
    (void)alea_transition_slice_critical_tile_source_count(transition_slice);
    (void)alea_transition_slice_critical_tile_source_get(
        transition_slice, 0, &critical_tile_source);
    (void)alea_transition_slice_critical_stop_reason_name(
        ALEA_TRANSITION_SLICE_CRITICAL_DISABLED);
    (void)alea_transition_slice_stop_reason_name(
        ALEA_TRANSITION_SLICE_STOP_NONE);
    (void)alea_transition_slice_refinement_status_name(
        ALEA_TRANSITION_SLICE_REFINEMENT_NOT_REQUESTED);
    (void)alea_ray_boundary_event_provenance_get(
        event_result, 0, &event_provenance);
    (void)alea_ray_coverage_slice_row_count(coverage);
    (void)alea_ray_coverage_slice_interval_count(coverage);
    (void)alea_ray_coverage_slice_owner_count(coverage);
    (void)alea_ray_coverage_slice_refinement_status(coverage);
    (void)alea_ray_coverage_slice_row_offsets(coverage);
    (void)alea_ray_coverage_slice_row_direction_tags(coverage);
    (void)alea_ray_coverage_slice_row_transverse_coordinates(coverage);
    (void)alea_ray_coverage_slice_t_enter(coverage);
    (void)alea_ray_coverage_slice_t_exit(coverage);
    (void)alea_ray_coverage_slice_kinds(coverage);
    (void)alea_ray_coverage_slice_owner_offsets(coverage);
    (void)alea_ray_coverage_slice_owner_count_lower_bounds(coverage);
    (void)alea_ray_coverage_slice_owner_cell_ids(coverage);
    (void)alea_ray_coverage_slice_owner_material_ids(coverage);
    (void)alea_ray_coverage_slice_owner_universe_ids(coverage);
    (void)alea_ray_coverage_slice_owner_fill_universes(coverage);
    (void)alea_ray_coverage_slice_owner_depths(coverage);
    (void)alea_ray_coverage_slice_owner_occurrence_keys(coverage);
    (void)alea_ray_coverage_slice_owner_parent_occurrence_keys(coverage);
    (void)alea_ray_coverage_slice_owner_resolution_flags(coverage);
    alea_log_set_callback(NULL, NULL);
    alea_raycast_result_destroy(result);
    alea_ray_boundary_event_query_result_destroy(event_result);
    alea_transition_slice_result_destroy(transition_slice);
    alea_ray_coverage_slice_result_destroy(coverage);
    alea_destroy(sys);
}
