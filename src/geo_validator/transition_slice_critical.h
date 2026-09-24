// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_TRANSITION_SLICE_CRITICAL_H
#define ALEA_TRANSITION_SLICE_CRITICAL_H

#include "alea_geo_validator.h"

typedef int (*alea_transition_slice_critical_finding_sink_t)(
    const alea_transition_slice_critical_finding_t* finding, void* userdata);

typedef struct {
    size_t tile_index;
    double uv_min[2];
    double uv_max[2];
    double uncertainty;
    alea_slice_error_numerical_cause_t cause;
    int source_cell_id;
    int source_surface_id;
    uint32_t source_primitive_id;
    uint64_t source_occurrence_key;
    uint64_t source_universe_occurrence_key;
} alea_transition_slice_numerical_region_t;

typedef int (*alea_transition_slice_numerical_region_sink_t)(
    const alea_transition_slice_numerical_region_t* region, void* userdata);

int alea_transition_slice_enumerate_critical_tiles(
    alea_system_t* sys,
    const alea_slice_view_t* view,
    const alea_transition_slice_options_t* options,
    const alea_transition_slice_critical_tile_t* tiles,
    size_t tile_count,
    alea_transition_slice_critical_finding_sink_t finding_sink,
    void* finding_sink_userdata,
    alea_transition_slice_stats_t* stats);

/* As above, with one stop reason per input tile. The caller supplies at
 * least tile_count entries. A NONE entry means that tile's requested scan
 * completed, even when another tile stopped the call-wide receipt. */
int alea_transition_slice_enumerate_critical_tiles_with_status(
    alea_system_t* sys,
    const alea_slice_view_t* view,
    const alea_transition_slice_options_t* options,
    const alea_transition_slice_critical_tile_t* tiles,
    size_t tile_count,
    alea_transition_slice_critical_finding_sink_t finding_sink,
    void* finding_sink_userdata,
    alea_transition_slice_stats_t* stats,
    alea_transition_slice_critical_stop_reason_t* tile_stop_reasons,
    size_t tile_stop_reason_capacity);

int alea_transition_slice_enumerate_critical_tiles_with_regions(
    alea_system_t* sys,
    const alea_slice_view_t* view,
    const alea_transition_slice_options_t* options,
    const alea_transition_slice_critical_tile_t* tiles,
    size_t tile_count,
    alea_transition_slice_critical_finding_sink_t finding_sink,
    void* finding_sink_userdata,
    alea_transition_slice_stats_t* stats,
    alea_transition_slice_critical_stop_reason_t* tile_stop_reasons,
    size_t tile_stop_reason_capacity,
    alea_transition_slice_numerical_region_sink_t numerical_region_sink,
    void* numerical_region_sink_userdata);

/* Internal tile-local key conversion used by critical point deduplication. */
int alea_transition_slice_quantize_point_for_tile(
    const alea_transition_slice_critical_tile_t* tile,
    double tolerance, double u, double v,
    int64_t* qu, int64_t* qv);

#endif
