// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_geo_validator.h"
#include "transition_slice_critical.h"

#include "raycast/raycast.h"
#include "raycast/ray_epsilon.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

struct alea_transition_slice_result {
    alea_transition_slice_finding_t* findings;
    size_t finding_count;
    size_t finding_capacity;
    alea_transition_slice_component_t* components;
    size_t component_count;
    size_t component_capacity;
    alea_transition_slice_coverage_finding_t* coverage_findings;
    size_t coverage_finding_count;
    size_t coverage_finding_capacity;
    alea_transition_slice_coverage_component_t* coverage_components;
    size_t coverage_component_count;
    size_t coverage_component_capacity;
    alea_transition_slice_component_link_t* component_links;
    size_t component_link_count;
    size_t component_link_capacity;
    alea_transition_slice_refinement_frontier_t* refinement_frontiers;
    size_t refinement_frontier_count;
    size_t refinement_frontier_capacity;
    alea_transition_slice_critical_tile_t* critical_tiles;
    size_t critical_tile_count;
    alea_transition_slice_critical_tile_source_t* critical_tile_sources;
    size_t critical_tile_source_count;
    alea_transition_slice_critical_finding_t* critical_findings;
    size_t critical_finding_count;
    size_t critical_finding_capacity;
    alea_transition_slice_stats_t stats;
};

typedef struct {
    double transverse_coordinate;
    size_t base_ray_index;
    uint32_t refinement_depth;
    uint64_t signature_a;
    uint64_t signature_b;
    size_t event_count;
    size_t finding_count;
} transition_slice_row_t;

typedef struct {
    alea_cell_hit_t* hits;
    uint64_t* occurrence_keys;
    uint64_t* parent_occurrence_keys;
    uint8_t* owner_mask;
    size_t capacity;
} transition_slice_coverage_scratch_t;

void alea_transition_slice_options_init(
    alea_transition_slice_options_t* options) {
    if (!options) return;
    memset(options, 0, sizeof(*options));
    options->struct_size = sizeof(*options);
    options->horizontal_rays = 64;
    options->vertical_rays = 64;
    options->max_rays = 4096;
    options->max_events = 1000000;
    options->max_events_per_ray = 8192;
    options->max_findings = 10000;
    options->max_components = 10000;
    options->max_output_bytes = 64u * 1024u * 1024u;
    options->max_scratch_bytes = 16u * 1024u * 1024u;
    options->max_coverage_fallbacks = 10000;
    options->max_coverage_hits = 256;
    options->refine_signals = ALEA_TRANSITION_SLICE_REFINE_SIGNATURE;
    options->max_row_scratch_bytes = 1024u * 1024u;
    options->max_coverage_probes = 10000;
    options->max_coverage_findings = 10000;
    options->max_coverage_components = 10000;
    options->max_component_links = 10000;
    options->max_refinement_frontiers = 1024;
    options->max_critical_tiles = 256;
    options->max_critical_tile_sources = 1024;
    options->max_critical_scratch_bytes = 4u * 1024u * 1024u;
    options->max_curves_per_tile = 1024;
    options->max_critical_points = 2048;
    options->max_active_boundary_tests = 100000;
    options->max_critical_probes = 4096;
    options->max_critical_findings = 1024;
    options->max_curve_pairs = 100000;
    options->max_critical_sector_witnesses = 8192;
}

const char* alea_transition_slice_critical_stop_reason_name(
    alea_transition_slice_critical_stop_reason_t reason) {
    switch (reason) {
    case ALEA_TRANSITION_SLICE_CRITICAL_DISABLED: return "disabled";
    case ALEA_TRANSITION_SLICE_CRITICAL_NONE: return "none";
    case ALEA_TRANSITION_SLICE_CRITICAL_MAX_FRONTIERS:
        return "max_refinement_frontiers";
    case ALEA_TRANSITION_SLICE_CRITICAL_MAX_TILE_SOURCES:
        return "max_critical_tile_sources";
    case ALEA_TRANSITION_SLICE_CRITICAL_MAX_TILES:
        return "max_critical_tiles";
    case ALEA_TRANSITION_SLICE_CRITICAL_MAX_SCRATCH_BYTES:
        return "max_critical_scratch_bytes";
    case ALEA_TRANSITION_SLICE_CRITICAL_MAX_OUTPUT_BYTES:
        return "max_output_bytes";
    case ALEA_TRANSITION_SLICE_CRITICAL_MAX_CURVES:
        return "max_curves_per_tile";
    case ALEA_TRANSITION_SLICE_CRITICAL_CHAIN_TRUNCATED:
        return "occurrence_chain_truncated";
    case ALEA_TRANSITION_SLICE_CRITICAL_MAX_POINTS:
        return "max_critical_points";
    case ALEA_TRANSITION_SLICE_CRITICAL_MAX_PROBES:
        return "max_critical_probes";
    case ALEA_TRANSITION_SLICE_CRITICAL_MAX_FINDINGS:
        return "max_critical_findings";
    case ALEA_TRANSITION_SLICE_CRITICAL_UNSUPPORTED_CURVE:
        return "unsupported_curve";
    case ALEA_TRANSITION_SLICE_CRITICAL_MAX_CURVE_PAIRS:
        return "max_curve_pairs";
    case ALEA_TRANSITION_SLICE_CRITICAL_MAX_SECTOR_WITNESSES:
        return "max_critical_sector_witnesses";
    }
    return "unknown";
}

const char* alea_transition_slice_refinement_status_name(
    alea_transition_slice_refinement_status_t status) {
    switch (status) {
    case ALEA_TRANSITION_SLICE_REFINEMENT_NOT_REQUESTED:
        return "not_requested";
    case ALEA_TRANSITION_SLICE_REFINEMENT_CONVERGED: return "converged";
    case ALEA_TRANSITION_SLICE_REFINEMENT_MAX_DEPTH: return "max_depth";
    case ALEA_TRANSITION_SLICE_REFINEMENT_MIN_SPACING: return "min_spacing";
    case ALEA_TRANSITION_SLICE_REFINEMENT_STOPPED: return "stopped";
    }
    return "unknown";
}

const char* alea_transition_slice_stop_reason_name(
    alea_transition_slice_stop_reason_t reason) {
    switch (reason) {
    case ALEA_TRANSITION_SLICE_STOP_NONE: return "none";
    case ALEA_TRANSITION_SLICE_STOP_MAX_RAYS: return "max_rays";
    case ALEA_TRANSITION_SLICE_STOP_MAX_EVENTS: return "max_events";
    case ALEA_TRANSITION_SLICE_STOP_MAX_FINDINGS: return "max_findings";
    case ALEA_TRANSITION_SLICE_STOP_MAX_COMPONENTS: return "max_components";
    case ALEA_TRANSITION_SLICE_STOP_MAX_OUTPUT_BYTES: return "max_output_bytes";
    case ALEA_TRANSITION_SLICE_STOP_MAX_COVERAGE_FALLBACKS:
        return "max_coverage_fallbacks";
    case ALEA_TRANSITION_SLICE_STOP_MAX_COVERAGE_PROBES:
        return "max_coverage_probes";
    case ALEA_TRANSITION_SLICE_STOP_MAX_SCRATCH_BYTES:
        return "max_scratch_bytes";
    case ALEA_TRANSITION_SLICE_STOP_INTERRUPTED: return "interrupted";
    case ALEA_TRANSITION_SLICE_STOP_MAX_COMPONENT_LINKS:
        return "max_component_links";
    }
    return "unknown";
}

alea_transition_slice_result_t* alea_transition_slice_result_create(void) {
    return calloc(1, sizeof(alea_transition_slice_result_t));
}

void alea_transition_slice_result_destroy(
    alea_transition_slice_result_t* result) {
    if (!result) return;
    free(result->findings);
    free(result->components);
    free(result->coverage_findings);
    free(result->coverage_components);
    free(result->component_links);
    free(result->refinement_frontiers);
    free(result->critical_tiles);
    free(result->critical_tile_sources);
    free(result->critical_findings);
    free(result);
}

size_t alea_transition_slice_component_count(
    const alea_transition_slice_result_t* result) {
    return result ? result->component_count : 0;
}

int alea_transition_slice_component_get(
    const alea_transition_slice_result_t* result, size_t index,
    alea_transition_slice_component_t* out_component) {
    if (!result || !out_component || index >= result->component_count)
        return -1;
    *out_component = result->components[index];
    return 0;
}

size_t alea_transition_slice_coverage_finding_count(
    const alea_transition_slice_result_t* result) {
    return result ? result->coverage_finding_count : 0;
}

int alea_transition_slice_coverage_finding_get(
    const alea_transition_slice_result_t* result, size_t index,
    alea_transition_slice_coverage_finding_t* out_finding) {
    if (!result || !out_finding || index >= result->coverage_finding_count)
        return -1;
    *out_finding = result->coverage_findings[index];
    return 0;
}

size_t alea_transition_slice_coverage_component_count(
    const alea_transition_slice_result_t* result) {
    return result ? result->coverage_component_count : 0;
}

int alea_transition_slice_coverage_component_get(
    const alea_transition_slice_result_t* result, size_t index,
    alea_transition_slice_coverage_component_t* out_component) {
    if (!result || !out_component || index >= result->coverage_component_count)
        return -1;
    *out_component = result->coverage_components[index];
    return 0;
}

size_t alea_transition_slice_component_link_count(
    const alea_transition_slice_result_t* result) {
    return result ? result->component_link_count : 0;
}

int alea_transition_slice_component_link_get(
    const alea_transition_slice_result_t* result, size_t index,
    alea_transition_slice_component_link_t* out_link) {
    if (!result || !out_link || index >= result->component_link_count)
        return -1;
    *out_link = result->component_links[index];
    return 0;
}

size_t alea_transition_slice_refinement_frontier_count(
    const alea_transition_slice_result_t* result) {
    return result ? result->refinement_frontier_count : 0;
}

int alea_transition_slice_refinement_frontier_get(
    const alea_transition_slice_result_t* result, size_t index,
    alea_transition_slice_refinement_frontier_t* out_frontier) {
    if (!result || !out_frontier ||
        index >= result->refinement_frontier_count) return -1;
    *out_frontier = result->refinement_frontiers[index];
    return 0;
}

size_t alea_transition_slice_critical_tile_count(
    const alea_transition_slice_result_t* result) {
    return result ? result->critical_tile_count : 0;
}

int alea_transition_slice_critical_tile_get(
    const alea_transition_slice_result_t* result, size_t index,
    alea_transition_slice_critical_tile_t* out_tile) {
    if (!result || !out_tile || index >= result->critical_tile_count) return -1;
    *out_tile = result->critical_tiles[index];
    return 0;
}

size_t alea_transition_slice_critical_tile_source_count(
    const alea_transition_slice_result_t* result) {
    return result ? result->critical_tile_source_count : 0;
}

int alea_transition_slice_critical_tile_source_get(
    const alea_transition_slice_result_t* result, size_t index,
    alea_transition_slice_critical_tile_source_t* out_source) {
    if (!result || !out_source ||
        index >= result->critical_tile_source_count) return -1;
    *out_source = result->critical_tile_sources[index];
    return 0;
}

size_t alea_transition_slice_critical_finding_count(
    const alea_transition_slice_result_t* result) {
    return result ? result->critical_finding_count : 0;
}

int alea_transition_slice_critical_finding_get(
    const alea_transition_slice_result_t* result, size_t index,
    alea_transition_slice_critical_finding_t* out_finding) {
    if (!result || !out_finding || index >= result->critical_finding_count)
        return -1;
    *out_finding = result->critical_findings[index];
    return 0;
}

size_t alea_transition_slice_finding_count(
    const alea_transition_slice_result_t* result) {
    return result ? result->finding_count : 0;
}

int alea_transition_slice_finding_get(
    const alea_transition_slice_result_t* result, size_t index,
    alea_transition_slice_finding_t* out_finding) {
    if (!result || !out_finding || index >= result->finding_count) return -1;
    *out_finding = result->findings[index];
    return 0;
}

int alea_transition_slice_stats(
    const alea_transition_slice_result_t* result,
    alea_transition_slice_stats_t* out_stats) {
    if (!result || !out_stats) return -1;
    *out_stats = result->stats;
    return 0;
}

static int transition_slice_add_bytes(size_t* total, size_t count,
                                      size_t item_size) {
    if (count > SIZE_MAX / item_size ||
        *total > SIZE_MAX - count * item_size)
        return -1;
    *total += count * item_size;
    return 0;
}

static size_t transition_slice_retained_bytes(
    const alea_transition_slice_result_t* result) {
    size_t bytes = 0;
    if (transition_slice_add_bytes(&bytes, result->finding_count,
                                   sizeof(*result->findings)) ||
        transition_slice_add_bytes(&bytes, result->component_count,
                                   sizeof(*result->components)) ||
        transition_slice_add_bytes(&bytes, result->coverage_finding_count,
                                   sizeof(*result->coverage_findings)) ||
        transition_slice_add_bytes(&bytes, result->coverage_component_count,
                                   sizeof(*result->coverage_components)) ||
        transition_slice_add_bytes(&bytes, result->component_link_count,
                                   sizeof(*result->component_links)) ||
        transition_slice_add_bytes(&bytes, result->refinement_frontier_count,
                                   sizeof(*result->refinement_frontiers)) ||
        transition_slice_add_bytes(&bytes, result->critical_tile_count,
                                   sizeof(*result->critical_tiles)) ||
        transition_slice_add_bytes(&bytes, result->critical_tile_source_count,
                                   sizeof(*result->critical_tile_sources)) ||
        transition_slice_add_bytes(&bytes, result->critical_finding_count,
                                   sizeof(*result->critical_findings)))
        return SIZE_MAX;
    return bytes;
}

typedef struct {
    alea_transition_slice_result_t* result;
    const alea_transition_slice_options_t* options;
} transition_slice_critical_sink_context_t;

static int transition_slice_append_critical_finding(
    const alea_transition_slice_critical_finding_t* finding,
    void* userdata) {
    transition_slice_critical_sink_context_t* context = userdata;
    alea_transition_slice_result_t* result = context->result;
    const alea_transition_slice_options_t* options = context->options;
    if (options->max_critical_findings &&
        result->critical_finding_count >= options->max_critical_findings) {
        result->stats.omitted_critical_findings++;
        return 1;
    }
    const size_t retained = transition_slice_retained_bytes(result);
    if (options->max_output_bytes &&
        (retained > options->max_output_bytes ||
         sizeof(*result->critical_findings) >
            options->max_output_bytes - retained)) {
        result->stats.omitted_critical_findings++;
        return 1;
    }
    if (result->critical_finding_count == result->critical_finding_capacity) {
        size_t capacity = result->critical_finding_capacity
            ? result->critical_finding_capacity * 2u : 16u;
        if (options->max_critical_findings &&
            capacity > options->max_critical_findings)
            capacity = (size_t)options->max_critical_findings;
        void* memory = realloc(
            result->critical_findings,
            capacity * sizeof(*result->critical_findings));
        if (!memory) return -1;
        result->critical_findings = memory;
        result->critical_finding_capacity = capacity;
    }
    result->critical_findings[result->critical_finding_count++] = *finding;
    result->stats.critical_findings = result->critical_finding_count;
    result->stats.retained_output_bytes = transition_slice_retained_bytes(result);
    return 0;
}

static size_t transition_slice_coverage_scratch_bytes(size_t capacity) {
    size_t bytes = 0;
    if (transition_slice_add_bytes(&bytes, capacity, sizeof(alea_cell_hit_t)) ||
        transition_slice_add_bytes(&bytes, capacity, sizeof(uint64_t)) ||
        transition_slice_add_bytes(&bytes, capacity, sizeof(uint64_t)) ||
        transition_slice_add_bytes(&bytes, capacity, sizeof(uint8_t)))
        return SIZE_MAX;
    return bytes;
}

static int transition_slice_coverage_scratch_init(
    transition_slice_coverage_scratch_t* scratch, size_t capacity) {
    memset(scratch, 0, sizeof(*scratch));
    if (capacity == 0) return 0;
    scratch->hits = calloc(capacity, sizeof(*scratch->hits));
    scratch->occurrence_keys = calloc(capacity, sizeof(*scratch->occurrence_keys));
    scratch->parent_occurrence_keys =
        calloc(capacity, sizeof(*scratch->parent_occurrence_keys));
    scratch->owner_mask = calloc(capacity, sizeof(*scratch->owner_mask));
    if (!scratch->hits || !scratch->occurrence_keys ||
        !scratch->parent_occurrence_keys || !scratch->owner_mask) {
        free(scratch->hits); free(scratch->occurrence_keys);
        free(scratch->parent_occurrence_keys); free(scratch->owner_mask);
        memset(scratch, 0, sizeof(*scratch));
        return -1;
    }
    scratch->capacity = capacity;
    return 0;
}

static void transition_slice_coverage_scratch_free(
    transition_slice_coverage_scratch_t* scratch) {
    free(scratch->hits); free(scratch->occurrence_keys);
    free(scratch->parent_occurrence_keys); free(scratch->owner_mask);
    memset(scratch, 0, sizeof(*scratch));
}

static int transition_slice_append(
    alea_transition_slice_result_t* result,
    const alea_transition_slice_options_t* options,
    const alea_transition_slice_finding_t* finding) {
    if (options->max_findings &&
        result->finding_count >= options->max_findings) {
        result->stats.stop_reason = ALEA_TRANSITION_SLICE_STOP_MAX_FINDINGS;
        return 1;
    }
    size_t retained = transition_slice_retained_bytes(result);
    if (options->max_output_bytes &&
        (retained > options->max_output_bytes ||
         sizeof(*result->findings) > options->max_output_bytes - retained)) {
        result->stats.stop_reason =
            ALEA_TRANSITION_SLICE_STOP_MAX_OUTPUT_BYTES;
        return 1;
    }
    if (result->finding_count == result->finding_capacity) {
        size_t capacity = result->finding_capacity
            ? result->finding_capacity * 2 : 16;
        size_t hard_capacity = SIZE_MAX / sizeof(*result->findings);
        if (options->max_findings && hard_capacity > options->max_findings)
            hard_capacity = (size_t)options->max_findings;
        if (options->max_output_bytes) {
            size_t byte_capacity = (size_t)(options->max_output_bytes /
                                             sizeof(*result->findings));
            if (hard_capacity > byte_capacity) hard_capacity = byte_capacity;
        }
        if (capacity > hard_capacity) capacity = hard_capacity;
        if (capacity <= result->finding_capacity) {
            result->stats.stop_reason = options->max_findings &&
                result->finding_capacity >= options->max_findings
                ? ALEA_TRANSITION_SLICE_STOP_MAX_FINDINGS
                : ALEA_TRANSITION_SLICE_STOP_MAX_OUTPUT_BYTES;
            return 1;
        }
        void* memory = realloc(
            result->findings, capacity * sizeof(*result->findings));
        if (!memory) return -1;
        result->findings = memory;
        result->finding_capacity = capacity;
    }
    result->findings[result->finding_count++] = *finding;
    result->stats.findings = result->finding_count;
    result->stats.retained_output_bytes = transition_slice_retained_bytes(result);
    return 0;
}

static int transition_slice_coverage_append(
    alea_transition_slice_result_t* result,
    const alea_transition_slice_options_t* options,
    const alea_transition_slice_coverage_finding_t* finding) {
    if (options->max_coverage_findings &&
        result->coverage_finding_count >= options->max_coverage_findings) {
        result->stats.stop_reason = ALEA_TRANSITION_SLICE_STOP_MAX_FINDINGS;
        return 1;
    }
    size_t retained = transition_slice_retained_bytes(result);
    if (options->max_output_bytes &&
        (retained > options->max_output_bytes ||
         sizeof(*result->coverage_findings) >
            options->max_output_bytes - retained)) {
        result->stats.stop_reason = ALEA_TRANSITION_SLICE_STOP_MAX_OUTPUT_BYTES;
        return 1;
    }
    if (result->coverage_finding_count == result->coverage_finding_capacity) {
        size_t capacity = result->coverage_finding_capacity
            ? result->coverage_finding_capacity * 2 : 16;
        size_t hard_capacity = SIZE_MAX / sizeof(*result->coverage_findings);
        if (options->max_coverage_findings &&
            hard_capacity > options->max_coverage_findings)
            hard_capacity = (size_t)options->max_coverage_findings;
        if (capacity > hard_capacity) capacity = hard_capacity;
        if (capacity <= result->coverage_finding_capacity) {
            result->stats.stop_reason = ALEA_TRANSITION_SLICE_STOP_MAX_FINDINGS;
            return 1;
        }
        void* memory = realloc(
            result->coverage_findings,
            capacity * sizeof(*result->coverage_findings));
        if (!memory) return -1;
        result->coverage_findings = memory;
        result->coverage_finding_capacity = capacity;
    }
    result->coverage_findings[result->coverage_finding_count++] = *finding;
    result->stats.coverage_findings = result->coverage_finding_count;
    result->stats.retained_output_bytes = transition_slice_retained_bytes(result);
    return 0;
}

static void transition_slice_world_point(
    const alea_slice_view_t* view, double u, double v, double point[3]) {
    for (int axis = 0; axis < 3; axis++)
        point[axis] = view->plane.origin[axis] +
            u * view->plane.u_axis[axis] +
            v * view->plane.v_axis[axis];
}

static int transition_slice_probe_coverage(
    alea_system_t* sys, const alea_slice_view_t* view,
    const alea_transition_slice_options_t* options,
    alea_transition_slice_orientation_t orientation,
    size_t ray_index, size_t base_ray_index, uint32_t refinement_depth,
    double transverse_coordinate, const alea_ray_t* ray,
    double ray_t, double bracket_t_enter, double bracket_t_exit,
    transition_slice_coverage_scratch_t* scratch,
    alea_transition_slice_result_t* result) {
    if (options->max_coverage_probes &&
        result->stats.coverage_probes >= options->max_coverage_probes) {
        result->stats.stop_reason =
            ALEA_TRANSITION_SLICE_STOP_MAX_COVERAGE_PROBES;
        return 1;
    }
    double point[3];
    alea_ray_point_at(ray, ray_t, &point[0], &point[1], &point[2]);
    int hit_count = alea_find_all_cells_coverage_chain(
        sys, point[0], point[1], point[2], scratch->hits,
        scratch->occurrence_keys, scratch->parent_occurrence_keys,
        scratch->capacity);
    if (hit_count < 0) return -1;
    result->stats.coverage_probes++;

    alea_transition_slice_coverage_finding_t finding;
    memset(&finding, 0, sizeof(finding));
    finding.orientation = orientation;
    finding.ray_index = ray_index;
    finding.base_ray_index = base_ray_index;
    finding.refinement_depth = refinement_depth;
    finding.transverse_coordinate = transverse_coordinate;
    finding.ray_t = ray_t;
    finding.bracket_t_enter = bracket_t_enter;
    finding.bracket_t_exit = bracket_t_exit;
    memcpy(finding.world_point, point, sizeof(point));
    if (orientation == ALEA_TRANSITION_SLICE_HORIZONTAL) {
        finding.uv[0] = view->u_min + ray_t;
        finding.uv[1] = transverse_coordinate;
    } else {
        finding.uv[0] = transverse_coordinate;
        finding.uv[1] = view->v_min + ray_t;
    }

    if ((size_t)hit_count >= scratch->capacity) {
        finding.kind = ALEA_POINT_COVERAGE_UNRESOLVED;
        finding.truncated = 1;
        finding.target_depth = -1;
        finding.owner_count_lower_bound = scratch->capacity;
        result->stats.truncated_coverage_probes++;
    } else {
        alea_point_coverage_classification_t classification;
        memset(scratch->owner_mask, 0,
               scratch->capacity * sizeof(*scratch->owner_mask));
        if (alea_classify_point_coverage_chain(
                scratch->hits, scratch->occurrence_keys,
                scratch->parent_occurrence_keys, (size_t)hit_count, -1,
                scratch->owner_mask, &classification) != 0)
            return -1;
        finding.kind = classification.kind;
        finding.target_depth = classification.target_depth;
        finding.owner_count_lower_bound = classification.owner_count;
        if (classification.kind == ALEA_POINT_COVERAGE_UNIQUE) {
            result->stats.unique_coverage_probes++;
            return 0;
        }
        if (classification.kind == ALEA_POINT_COVERAGE_GAP &&
            !options->report_unowned_coverage) {
            result->stats.skipped_unowned_coverage_probes++;
            return 0;
        }
        for (int hit = 0; hit < hit_count &&
             finding.owner_count <
                ALEA_TRANSITION_SLICE_COVERAGE_OWNER_CAPACITY; hit++) {
            if (!scratch->owner_mask[hit]) continue;
            size_t owner = finding.owner_count++;
            finding.owner_cell_ids[owner] = scratch->hits[hit].cell_id;
            finding.owner_universe_ids[owner] = scratch->hits[hit].universe_id;
            finding.owner_depths[owner] = scratch->hits[hit].depth;
            finding.owner_occurrence_keys[owner] =
                scratch->occurrence_keys[hit];
            finding.owner_parent_occurrence_keys[owner] =
                scratch->parent_occurrence_keys[hit];
        }
    }
    return transition_slice_coverage_append(result, options, &finding);
}

static size_t transition_slice_event_cap(
    const alea_transition_slice_options_t* options,
    const alea_transition_slice_stats_t* stats,
    alea_transition_slice_stop_reason_t* limiting_reason) {
    uint64_t cap = options->max_events_per_ray;
    *limiting_reason = ALEA_TRANSITION_SLICE_STOP_MAX_EVENTS;
    if (options->max_events) {
        uint64_t remaining = stats->events_checked < options->max_events
            ? options->max_events - stats->events_checked : 0;
        if (!cap || remaining < cap) cap = remaining;
    }
    if (options->max_scratch_bytes) {
        /* Vector growth is geometric. Half the byte-derived count guarantees
         * its capacity remains under the caller's byte ceiling. */
        uint64_t scratch_cap = options->max_scratch_bytes /
            (2u * sizeof(alea_ray_boundary_event_t));
        if (scratch_cap == 0 && options->max_scratch_bytes >=
                sizeof(alea_ray_boundary_event_t)) scratch_cap = 1;
        if (!cap || scratch_cap < cap) {
            cap = scratch_cap;
            *limiting_reason = ALEA_TRANSITION_SLICE_STOP_MAX_SCRATCH_BYTES;
        }
    }
    return cap > SIZE_MAX ? SIZE_MAX : (size_t)cap;
}

static uint64_t transition_slice_hash(uint64_t hash, uint64_t value) {
    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    value ^= value >> 31;
    return hash ^ (value + (hash << 6) + (hash >> 2));
}

static void transition_slice_row_signature(
    const alea_ray_boundary_event_result_t* events,
    transition_slice_row_t* row) {
    uint64_t a = UINT64_C(0xcbf29ce484222325);
    uint64_t b = UINT64_C(0x84222325cbf29ce4);
    for (size_t i = 0; i < events->events.count; i++) {
        const alea_ray_boundary_event_t* event = &events->events.data[i];
#define HASH_EVENT(VALUE) \
        do { \
            uint64_t value_ = (uint64_t)(VALUE); \
            a = transition_slice_hash(a, value_); \
            b = transition_slice_hash(b, value_ ^ (uint64_t)i); \
        } while (0)
        HASH_EVENT(event->kind);
        HASH_EVENT((uint32_t)event->surface_id);
        HASH_EVENT((uint32_t)event->cell_before);
        HASH_EVENT((uint32_t)event->cell_after);
        HASH_EVENT((uint32_t)event->active_cell_id);
        HASH_EVENT((uint32_t)event->active_universe_id);
        HASH_EVENT((uint32_t)event->active_depth);
        HASH_EVENT(event->active_occurrence_key);
        HASH_EVENT(event->active_parent_occurrence_key);
        HASH_EVENT(event->before_occurrence_key);
        HASH_EVENT(event->after_occurrence_key);
        HASH_EVENT(event->local_surface_complete);
        HASH_EVENT(event->local_surface_count);
        for (size_t surface = 0; surface < event->local_surface_count;
             surface++)
            HASH_EVENT((uint32_t)event->local_surface_ids[surface]);
#undef HASH_EVENT
    }
    row->signature_a = transition_slice_hash(a, events->events.count);
    row->signature_b = transition_slice_hash(b, events->events.count);
    row->event_count = events->events.count;
}

static int transition_slice_scan_ray(
    alea_system_t* sys, const alea_slice_view_t* view,
    const alea_transition_slice_options_t* options,
    alea_transition_slice_orientation_t orientation, size_t ray_index,
    size_t base_ray_index, uint32_t refinement_depth,
    double transverse_coordinate, alea_raycast_result_t* scratch,
    alea_ray_boundary_event_result_t* events,
    transition_slice_coverage_scratch_t* coverage_scratch,
    alea_transition_slice_result_t* result, transition_slice_row_t* row) {
    const double u_span = view->u_max - view->u_min;
    const double v_span = view->v_max - view->v_min;
    double u, v, origin[3], direction[3], length;
    if (orientation == ALEA_TRANSITION_SLICE_HORIZONTAL) {
        u = view->u_min;
        v = transverse_coordinate;
        direction[0] = view->plane.u_axis[0];
        direction[1] = view->plane.u_axis[1];
        direction[2] = view->plane.u_axis[2];
        length = u_span;
    } else {
        u = transverse_coordinate;
        v = view->v_min;
        direction[0] = view->plane.v_axis[0];
        direction[1] = view->plane.v_axis[1];
        direction[2] = view->plane.v_axis[2];
        length = v_span;
    }
    transition_slice_world_point(view, u, v, origin);
    alea_ray_t ray;
    if (alea_ray_init(&ray, origin[0], origin[1], origin[2],
                      direction[0], direction[1], direction[2]) != 0)
        return -1;

    alea_transition_slice_stop_reason_t limiting_reason;
    size_t event_cap = transition_slice_event_cap(
        options, &result->stats, &limiting_reason);
    if (event_cap == 0) {
        result->stats.stop_reason = limiting_reason;
        return 1;
    }
    const alea_ray_boundary_event_options_internal_t event_options = {
        .max_events = event_cap,
        .max_output_bytes = options->max_scratch_bytes,
        .skip_open_side_coverage = true
    };
    if (alea_raycast_selected_boundary_events_with_options_nocache(
            sys, &ray, length, &event_options, scratch, events) != 0) {
        if (alea_error_code() != ALEA_ERR_OVERFLOW) return -1;
        result->stats.stop_reason = limiting_reason;
        return 1;
    }
    result->stats.executed_rays++;
    if (orientation == ALEA_TRANSITION_SLICE_HORIZONTAL)
        result->stats.horizontal_rays_executed++;
    else
        result->stats.vertical_rays_executed++;
    if (events->events.count > result->stats.peak_live_events)
        result->stats.peak_live_events = events->events.count;
    const size_t live_bytes = events->events.capacity * sizeof(*events->events.data);
    if (live_bytes > result->stats.peak_live_event_bytes)
        result->stats.peak_live_event_bytes = live_bytes;
    row->transverse_coordinate = transverse_coordinate;
    row->base_ray_index = base_ray_index;
    row->refinement_depth = refinement_depth;
    row->finding_count = 0;
    transition_slice_row_signature(events, row);
    const size_t findings_before = result->finding_count;
    const size_t coverage_findings_before = result->coverage_finding_count;

    if (options->coverage_uniform_probes_per_ray) {
        const size_t count = options->coverage_uniform_probes_per_ray;
        size_t selected_boundary = 0;
        for (size_t probe = 0; probe < count; probe++) {
            double bin_enter = length * (double)probe / (double)count;
            double bin_exit = length * (double)(probe + 1) / (double)count;
            double ray_t = bin_enter + 0.381966011250105 *
                (bin_exit - bin_enter);
            while (selected_boundary < events->events.count &&
                   events->events.data[selected_boundary].t < ray_t)
                selected_boundary++;
            double t_enter = selected_boundary
                ? events->events.data[selected_boundary - 1].t : 0.0;
            double t_exit = selected_boundary < events->events.count
                ? events->events.data[selected_boundary].t : length;
            int rc = transition_slice_probe_coverage(
                sys, view, options, orientation, ray_index, base_ray_index,
                refinement_depth, transverse_coordinate, &ray,
                ray_t, t_enter, t_exit, coverage_scratch, result);
            if (rc != 0) return rc;
        }
    }
    if (options->coverage_probe_selected_intervals) {
        double t_enter = 0.0;
        for (size_t boundary = 0; boundary <= events->events.count; boundary++) {
            double t_exit = boundary < events->events.count
                ? events->events.data[boundary].t : length;
            if (t_exit - t_enter > 64.0 * RAY_EPSILON) {
                double ray_t = t_enter + 0.381966011250105 *
                    (t_exit - t_enter);
                int rc = transition_slice_probe_coverage(
                    sys, view, options, orientation, ray_index, base_ray_index,
                    refinement_depth, transverse_coordinate, &ray,
                    ray_t, t_enter, t_exit, coverage_scratch, result);
                if (rc != 0) return rc;
            }
            if (t_exit > t_enter) t_enter = t_exit;
        }
    }

    for (size_t event_index = 0; event_index < events->events.count;
         event_index++) {
        const alea_ray_boundary_event_t* event = &events->events.data[event_index];
        result->stats.events_checked++;
        if (event->kind != ALEA_RAY_BOUNDARY_EVENT_PHYSICAL) continue;
        result->stats.physical_events_seen++;
        if (!options->include_void_transitions &&
            (event->cell_before < 0 || event->cell_after < 0)) {
            result->stats.skipped_void_transitions++;
            continue;
        }
        if (options->max_coverage_fallbacks &&
            result->stats.coverage_fallbacks >=
                options->max_coverage_fallbacks) {
            result->stats.stop_reason =
                ALEA_TRANSITION_SLICE_STOP_MAX_COVERAGE_FALLBACKS;
            return 1;
        }

        double previous_t = event_index ? events->events.data[event_index - 1].t : 0.0;
        double next_t = event_index + 1 < events->events.count
            ? events->events.data[event_index + 1].t : length;
        double safe_max = 0.25 * fmin(event->t - previous_t, next_t - event->t);
        alea_transition_result_t transition;
        memset(&transition, 0, sizeof(transition));
        if (!(safe_max > 64.0 * RAY_EPSILON)) {
            transition.kind = ALEA_TRANSITION_AMBIGUOUS_BOUNDARY;
            transition.after_coverage_kind = ALEA_POINT_COVERAGE_UNRESOLVED;
            transition.universe_id = event->active_universe_id;
            transition.current_cell_id = event->active_cell_id;
            transition.primary_surface_id = event->surface_id;
            transition.connecting_surface_id = -1;
            transition.after_cell_id = event->cell_after;
            transition.occurrence_depth = event->active_depth;
            memcpy(transition.crossing_point, event->local_point,
                   sizeof(transition.crossing_point));
            memcpy(transition.direction, event->local_direction,
                   sizeof(transition.direction));
        } else {
            alea_transition_options_t transition_options;
            alea_transition_options_init(&transition_options);
            if (options->probe_distance > 0.0)
                transition_options.probe_distance = options->probe_distance;
            if (transition_options.probe_distance > safe_max)
                transition_options.probe_distance = safe_max;
            transition_options.max_probe_distance = safe_max;
            if (options->max_probe_distance > 0.0 &&
                transition_options.max_probe_distance >
                    options->max_probe_distance)
                transition_options.max_probe_distance =
                    options->max_probe_distance;
            transition_options.max_coverage_hits = options->max_coverage_hits;
            if (options->max_coverage_fallbacks)
                transition_options.max_coverage_fallbacks =
                    options->max_coverage_fallbacks -
                    result->stats.coverage_fallbacks;
            if (alea_check_selected_boundary_event_transition_nocache(
                    sys, event, &transition_options, &transition) != 0)
                return -1;
        }
        result->stats.coverage_fallbacks += transition.coverage_fallbacks;
        if (transition.kind == ALEA_TRANSITION_VALID) {
            result->stats.valid_transitions++;
            continue;
        }
        alea_transition_slice_finding_t finding;
        memset(&finding, 0, sizeof(finding));
        finding.transition = transition;
        finding.orientation = orientation;
        finding.ray_index = ray_index;
        finding.event_index = event_index;
        finding.base_ray_index = base_ray_index;
        finding.refinement_depth = refinement_depth;
        finding.transverse_coordinate = transverse_coordinate;
        finding.ray_t = event->t;
        if (orientation == ALEA_TRANSITION_SLICE_HORIZONTAL) {
            finding.uv[0] = view->u_min + event->t;
            finding.uv[1] = v;
        } else {
            finding.uv[0] = u;
            finding.uv[1] = view->v_min + event->t;
        }
        alea_ray_point_at(&ray, event->t,
                          &finding.world_point[0], &finding.world_point[1],
                          &finding.world_point[2]);
        int append_rc = transition_slice_append(result, options, &finding);
        if (append_rc != 0) return append_rc;
        if (transition.kind == ALEA_TRANSITION_TRUNCATED &&
            options->max_coverage_fallbacks &&
            result->stats.coverage_fallbacks >=
                options->max_coverage_fallbacks) {
            result->stats.stop_reason =
                ALEA_TRANSITION_SLICE_STOP_MAX_COVERAGE_FALLBACKS;
            return 1;
        }
    }
    row->finding_count = result->finding_count - findings_before +
        result->coverage_finding_count - coverage_findings_before;
    return 0;
}

static int transition_slice_rows_differ(
    const transition_slice_row_t* first,
    const transition_slice_row_t* second, uint32_t signals) {
    if ((signals & ALEA_TRANSITION_SLICE_REFINE_SIGNATURE) &&
        (first->event_count != second->event_count ||
         first->signature_a != second->signature_a ||
         first->signature_b != second->signature_b))
        return 1;
    return (signals & ALEA_TRANSITION_SLICE_REFINE_FINDING) &&
        (first->finding_count != 0 || second->finding_count != 0);
}

static int transition_slice_append_frontier(
    const alea_slice_view_t* view,
    const alea_transition_slice_options_t* options,
    alea_transition_slice_orientation_t orientation,
    const transition_slice_row_t* first,
    const transition_slice_row_t* second,
    alea_transition_slice_result_t* result) {
    if (!options->enable_critical_refinement) return 0;
    if (options->max_refinement_frontiers &&
        result->refinement_frontier_count >=
            options->max_refinement_frontiers) {
        result->stats.omitted_refinement_frontiers++;
        if (result->stats.critical_stop_reason ==
            ALEA_TRANSITION_SLICE_CRITICAL_NONE)
            result->stats.critical_stop_reason =
                ALEA_TRANSITION_SLICE_CRITICAL_MAX_FRONTIERS;
        return 0;
    }
    size_t retained = transition_slice_retained_bytes(result);
    if (options->max_output_bytes &&
        (retained > options->max_output_bytes ||
         sizeof(*result->refinement_frontiers) >
             options->max_output_bytes - retained)) {
        result->stats.omitted_refinement_frontiers++;
        if (result->stats.critical_stop_reason ==
            ALEA_TRANSITION_SLICE_CRITICAL_NONE)
            result->stats.critical_stop_reason =
                ALEA_TRANSITION_SLICE_CRITICAL_MAX_OUTPUT_BYTES;
        return 0;
    }
    if (result->refinement_frontier_count ==
        result->refinement_frontier_capacity) {
        size_t capacity = result->refinement_frontier_capacity
            ? result->refinement_frontier_capacity * 2 : 16;
        if (options->max_refinement_frontiers &&
            capacity > options->max_refinement_frontiers)
            capacity = options->max_refinement_frontiers;
        void* memory = realloc(
            result->refinement_frontiers,
            capacity * sizeof(*result->refinement_frontiers));
        if (!memory) return -1;
        result->refinement_frontiers = memory;
        result->refinement_frontier_capacity = capacity;
    }
    alea_transition_slice_refinement_frontier_t* frontier =
        &result->refinement_frontiers[result->refinement_frontier_count++];
    memset(frontier, 0, sizeof(*frontier));
    frontier->orientation = orientation;
    frontier->refinement_depth = first->refinement_depth >
            second->refinement_depth
        ? first->refinement_depth : second->refinement_depth;
    frontier->transverse_min = first->transverse_coordinate;
    frontier->transverse_max = second->transverse_coordinate;
    frontier->signature_a[0] = first->signature_a;
    frontier->signature_a[1] = second->signature_a;
    frontier->signature_b[0] = first->signature_b;
    frontier->signature_b[1] = second->signature_b;
    frontier->max_event_count = first->event_count > second->event_count
        ? first->event_count : second->event_count;
    if (orientation == ALEA_TRANSITION_SLICE_HORIZONTAL) {
        frontier->uv_min[0] = view->u_min;
        frontier->uv_max[0] = view->u_max;
        frontier->uv_min[1] = frontier->transverse_min;
        frontier->uv_max[1] = frontier->transverse_max;
    } else {
        frontier->uv_min[0] = frontier->transverse_min;
        frontier->uv_max[0] = frontier->transverse_max;
        frontier->uv_min[1] = view->v_min;
        frontier->uv_max[1] = view->v_max;
    }
    result->stats.refinement_frontiers = result->refinement_frontier_count;
    result->stats.retained_output_bytes = transition_slice_retained_bytes(result);
    return 0;
}

static int transition_slice_retain_differing_frontiers(
    const alea_slice_view_t* view,
    const alea_transition_slice_options_t* options,
    alea_transition_slice_orientation_t orientation,
    const transition_slice_row_t* rows, size_t row_count,
    alea_transition_slice_result_t* result) {
    for (size_t row = 0; row + 1 < row_count; row++) {
        if (!transition_slice_rows_differ(
                &rows[row], &rows[row + 1], options->refine_signals))
            continue;
        if (transition_slice_append_frontier(
                view, options, orientation, &rows[row], &rows[row + 1],
                result) != 0)
            return -1;
    }
    return 0;
}

static void transition_slice_record_row_scratch(
    alea_transition_slice_result_t* result, size_t row_bytes) {
    if (row_bytes > result->stats.peak_row_scratch_bytes)
        result->stats.peak_row_scratch_bytes = row_bytes;
    size_t total = row_bytes;
    if (SIZE_MAX - total < result->stats.peak_live_event_bytes)
        total = SIZE_MAX;
    else
        total += result->stats.peak_live_event_bytes;
    if (total > result->stats.peak_scratch_bytes)
        result->stats.peak_scratch_bytes = total;
}

static int transition_slice_status_rank(
    alea_transition_slice_refinement_status_t status) {
    switch (status) {
    case ALEA_TRANSITION_SLICE_REFINEMENT_STOPPED: return 4;
    case ALEA_TRANSITION_SLICE_REFINEMENT_MAX_DEPTH: return 3;
    case ALEA_TRANSITION_SLICE_REFINEMENT_MIN_SPACING: return 2;
    case ALEA_TRANSITION_SLICE_REFINEMENT_CONVERGED: return 1;
    case ALEA_TRANSITION_SLICE_REFINEMENT_NOT_REQUESTED: return 0;
    }
    return 4;
}

static void transition_slice_merge_status(
    alea_transition_slice_result_t* result,
    alea_transition_slice_refinement_status_t status) {
    if (transition_slice_status_rank(status) >
        transition_slice_status_rank(result->stats.refinement_status))
        result->stats.refinement_status = status;
}

static int transition_slice_compare_int(int first, int second) {
    return first < second ? -1 : first > second ? 1 : 0;
}

static int transition_slice_compare_size(size_t first, size_t second) {
    return first < second ? -1 : first > second ? 1 : 0;
}

static int transition_slice_finding_compare(const void* first_pointer,
                                            const void* second_pointer) {
    const alea_transition_slice_finding_t* first = first_pointer;
    const alea_transition_slice_finding_t* second = second_pointer;
#define COMPARE_INT(FIELD) \
    do { \
        int comparison = transition_slice_compare_int( \
            (int)first->FIELD, (int)second->FIELD); \
        if (comparison) return comparison; \
    } while (0)
    COMPARE_INT(transition.kind);
    COMPARE_INT(orientation);
    COMPARE_INT(transition.universe_id);
    COMPARE_INT(transition.current_cell_id);
    COMPARE_INT(transition.after_cell_id);
    COMPARE_INT(transition.primary_surface_id);
    COMPARE_INT(transition.connecting_surface_id);
#undef COMPARE_INT
    if (first->transition.current_occurrence_key !=
        second->transition.current_occurrence_key)
        return first->transition.current_occurrence_key <
                second->transition.current_occurrence_key ? -1 : 1;
    for (int axis = 0; axis < 2; axis++) {
        if (first->uv[axis] != second->uv[axis])
            return first->uv[axis] < second->uv[axis] ? -1 : 1;
    }
    int comparison = transition_slice_compare_size(
        first->ray_index, second->ray_index);
    if (comparison) return comparison;
    return transition_slice_compare_size(
        first->event_index, second->event_index);
}

static int transition_slice_findings_same_component(
    const alea_transition_slice_finding_t* first,
    const alea_transition_slice_finding_t* second) {
    return first->transition.kind == second->transition.kind &&
        first->orientation == second->orientation &&
        first->transition.universe_id == second->transition.universe_id &&
        first->transition.current_cell_id == second->transition.current_cell_id &&
        first->transition.after_cell_id == second->transition.after_cell_id &&
        first->transition.primary_surface_id ==
            second->transition.primary_surface_id &&
        first->transition.connecting_surface_id ==
            second->transition.connecting_surface_id &&
        first->transition.current_occurrence_key ==
            second->transition.current_occurrence_key;
}

static int transition_slice_build_components(
    const alea_transition_slice_options_t* options,
    alea_transition_slice_result_t* result) {
    if (result->finding_count == 0) return 0;
    qsort(result->findings, result->finding_count,
          sizeof(*result->findings), transition_slice_finding_compare);
    size_t required = 1;
    for (size_t i = 1; i < result->finding_count; i++)
        if (!transition_slice_findings_same_component(
                &result->findings[i - 1], &result->findings[i]))
            required++;
    size_t capacity = required;
    if (options->max_components && capacity > options->max_components)
        capacity = (size_t)options->max_components;
    if (options->max_output_bytes) {
        size_t retained = transition_slice_retained_bytes(result);
        size_t available = retained < options->max_output_bytes
            ? (size_t)(options->max_output_bytes - retained) : 0;
        size_t byte_capacity = available / sizeof(*result->components);
        if (capacity > byte_capacity) capacity = byte_capacity;
    }
    if (capacity > 0) {
        result->components = calloc(capacity, sizeof(*result->components));
        if (!result->components) return -1;
        result->component_capacity = capacity;
    }
    size_t finding = 0;
    while (finding < result->finding_count &&
           result->component_count < capacity) {
        size_t end = finding + 1;
        while (end < result->finding_count &&
               transition_slice_findings_same_component(
                   &result->findings[finding], &result->findings[end]))
            end++;
        const alea_transition_slice_finding_t* representative =
            &result->findings[finding];
        alea_transition_slice_component_t* component =
            &result->components[result->component_count++];
        component->kind = representative->transition.kind;
        component->orientation = representative->orientation;
        component->universe_id = representative->transition.universe_id;
        component->current_cell_id = representative->transition.current_cell_id;
        component->after_cell_id = representative->transition.after_cell_id;
        component->primary_surface_id =
            representative->transition.primary_surface_id;
        component->connecting_surface_id =
            representative->transition.connecting_surface_id;
        component->current_occurrence_key =
            representative->transition.current_occurrence_key;
        component->first_finding_index = finding;
        component->finding_count = end - finding;
        for (int axis = 0; axis < 2; axis++)
            component->uv_min[axis] = component->uv_max[axis] =
                representative->uv[axis];
        for (int axis = 0; axis < 3; axis++)
            component->world_min[axis] = component->world_max[axis] =
                representative->world_point[axis];
        for (size_t i = finding; i < end; i++) {
            const alea_transition_slice_finding_t* item = &result->findings[i];
            if (component->max_refinement_depth < item->refinement_depth)
                component->max_refinement_depth = item->refinement_depth;
            for (int axis = 0; axis < 2; axis++) {
                if (component->uv_min[axis] > item->uv[axis])
                    component->uv_min[axis] = item->uv[axis];
                if (component->uv_max[axis] < item->uv[axis])
                    component->uv_max[axis] = item->uv[axis];
            }
            for (int axis = 0; axis < 3; axis++) {
                if (component->world_min[axis] > item->world_point[axis])
                    component->world_min[axis] = item->world_point[axis];
                if (component->world_max[axis] < item->world_point[axis])
                    component->world_max[axis] = item->world_point[axis];
            }
        }
        finding = end;
    }
    result->stats.components = result->component_count;
    result->stats.retained_output_bytes = transition_slice_retained_bytes(result);
    if (result->component_count < required &&
        result->stats.stop_reason == ALEA_TRANSITION_SLICE_STOP_NONE) {
        result->stats.stop_reason = options->max_components &&
                result->component_count >= options->max_components
            ? ALEA_TRANSITION_SLICE_STOP_MAX_COMPONENTS
            : ALEA_TRANSITION_SLICE_STOP_MAX_OUTPUT_BYTES;
        return 1;
    }
    return 0;
}

static int transition_slice_coverage_key_compare(
    const alea_transition_slice_coverage_finding_t* first,
    const alea_transition_slice_coverage_finding_t* second) {
    int comparison = transition_slice_compare_int(
        (int)first->kind, (int)second->kind);
    if (comparison) return comparison;
    comparison = transition_slice_compare_int(first->truncated, second->truncated);
    if (comparison) return comparison;
    comparison = transition_slice_compare_int(
        (int)first->orientation, (int)second->orientation);
    if (comparison) return comparison;
    comparison = transition_slice_compare_size(
        first->owner_count_lower_bound, second->owner_count_lower_bound);
    if (comparison) return comparison;
    comparison = transition_slice_compare_size(
        first->owner_count, second->owner_count);
    if (comparison) return comparison;
    for (size_t owner = 0; owner < first->owner_count; owner++) {
        comparison = transition_slice_compare_int(
            first->owner_cell_ids[owner], second->owner_cell_ids[owner]);
        if (comparison) return comparison;
        comparison = transition_slice_compare_int(
            first->owner_universe_ids[owner],
            second->owner_universe_ids[owner]);
        if (comparison) return comparison;
        if (first->owner_occurrence_keys[owner] !=
            second->owner_occurrence_keys[owner])
            return first->owner_occurrence_keys[owner] <
                second->owner_occurrence_keys[owner] ? -1 : 1;
    }
    return 0;
}

static int transition_slice_coverage_finding_compare(
    const void* first_pointer, const void* second_pointer) {
    const alea_transition_slice_coverage_finding_t* first = first_pointer;
    const alea_transition_slice_coverage_finding_t* second = second_pointer;
    int comparison = transition_slice_coverage_key_compare(first, second);
    if (comparison) return comparison;
    for (int axis = 0; axis < 2; axis++) {
        if (first->uv[axis] != second->uv[axis])
            return first->uv[axis] < second->uv[axis] ? -1 : 1;
    }
    comparison = transition_slice_compare_size(
        first->ray_index, second->ray_index);
    if (comparison) return comparison;
    return first->ray_t < second->ray_t ? -1 : first->ray_t > second->ray_t;
}

static int transition_slice_build_coverage_components(
    const alea_transition_slice_options_t* options,
    alea_transition_slice_result_t* result) {
    if (result->coverage_finding_count == 0) return 0;
    qsort(result->coverage_findings, result->coverage_finding_count,
          sizeof(*result->coverage_findings),
          transition_slice_coverage_finding_compare);
    size_t required = 1;
    for (size_t i = 1; i < result->coverage_finding_count; i++)
        if (transition_slice_coverage_key_compare(
                &result->coverage_findings[i - 1],
                &result->coverage_findings[i]) != 0)
            required++;
    size_t capacity = required;
    if (options->max_coverage_components &&
        capacity > options->max_coverage_components)
        capacity = (size_t)options->max_coverage_components;
    if (options->max_output_bytes) {
        size_t retained = transition_slice_retained_bytes(result);
        size_t available = retained < options->max_output_bytes
            ? (size_t)(options->max_output_bytes - retained) : 0;
        size_t byte_capacity = available / sizeof(*result->coverage_components);
        if (capacity > byte_capacity) capacity = byte_capacity;
    }
    if (capacity) {
        result->coverage_components =
            calloc(capacity, sizeof(*result->coverage_components));
        if (!result->coverage_components) return -1;
        result->coverage_component_capacity = capacity;
    }
    size_t finding = 0;
    while (finding < result->coverage_finding_count &&
           result->coverage_component_count < capacity) {
        size_t end = finding + 1;
        while (end < result->coverage_finding_count &&
               transition_slice_coverage_key_compare(
                   &result->coverage_findings[finding],
                   &result->coverage_findings[end]) == 0)
            end++;
        const alea_transition_slice_coverage_finding_t* representative =
            &result->coverage_findings[finding];
        alea_transition_slice_coverage_component_t* component =
            &result->coverage_components[result->coverage_component_count++];
        component->kind = representative->kind;
        component->truncated = representative->truncated;
        component->orientation = representative->orientation;
        component->first_finding_index = finding;
        component->finding_count = end - finding;
        component->owner_count_lower_bound =
            representative->owner_count_lower_bound;
        for (int axis = 0; axis < 2; axis++)
            component->uv_min[axis] = component->uv_max[axis] =
                representative->uv[axis];
        for (int axis = 0; axis < 3; axis++)
            component->world_min[axis] = component->world_max[axis] =
                representative->world_point[axis];
        for (size_t i = finding; i < end; i++) {
            const alea_transition_slice_coverage_finding_t* item =
                &result->coverage_findings[i];
            if (component->max_refinement_depth < item->refinement_depth)
                component->max_refinement_depth = item->refinement_depth;
            for (int axis = 0; axis < 2; axis++) {
                if (component->uv_min[axis] > item->uv[axis])
                    component->uv_min[axis] = item->uv[axis];
                if (component->uv_max[axis] < item->uv[axis])
                    component->uv_max[axis] = item->uv[axis];
            }
            for (int axis = 0; axis < 3; axis++) {
                if (component->world_min[axis] > item->world_point[axis])
                    component->world_min[axis] = item->world_point[axis];
                if (component->world_max[axis] < item->world_point[axis])
                    component->world_max[axis] = item->world_point[axis];
            }
        }
        finding = end;
    }
    result->stats.coverage_components = result->coverage_component_count;
    result->stats.retained_output_bytes = transition_slice_retained_bytes(result);
    if (result->coverage_component_count < required &&
        result->stats.stop_reason == ALEA_TRANSITION_SLICE_STOP_NONE) {
        result->stats.stop_reason = options->max_coverage_components &&
                result->coverage_component_count >=
                    options->max_coverage_components
            ? ALEA_TRANSITION_SLICE_STOP_MAX_COMPONENTS
            : ALEA_TRANSITION_SLICE_STOP_MAX_OUTPUT_BYTES;
        return 1;
    }
    return 0;
}

typedef struct {
    alea_transition_slice_orientation_t orientation;
    size_t ray_index;
    double ray_t;
    size_t component_index;
} transition_slice_component_boundary_ref_t;

typedef struct {
    alea_transition_slice_orientation_t orientation;
    size_t ray_index;
    double ray_t;
    size_t component_index;
    uint32_t side;
} transition_slice_coverage_endpoint_ref_t;

static int transition_slice_boundary_ref_compare(const void* first_pointer,
                                                 const void* second_pointer) {
    const transition_slice_component_boundary_ref_t* first = first_pointer;
    const transition_slice_component_boundary_ref_t* second = second_pointer;
    int comparison = transition_slice_compare_int(
        (int)first->orientation, (int)second->orientation);
    if (comparison) return comparison;
    comparison = transition_slice_compare_size(first->ray_index,
                                                second->ray_index);
    if (comparison) return comparison;
    if (first->ray_t != second->ray_t)
        return first->ray_t < second->ray_t ? -1 : 1;
    return transition_slice_compare_size(first->component_index,
                                         second->component_index);
}

static int transition_slice_endpoint_ref_compare(const void* first_pointer,
                                                 const void* second_pointer) {
    const transition_slice_coverage_endpoint_ref_t* first = first_pointer;
    const transition_slice_coverage_endpoint_ref_t* second = second_pointer;
    int comparison = transition_slice_compare_int(
        (int)first->orientation, (int)second->orientation);
    if (comparison) return comparison;
    comparison = transition_slice_compare_size(first->ray_index,
                                                second->ray_index);
    if (comparison) return comparison;
    if (first->ray_t != second->ray_t)
        return first->ray_t < second->ray_t ? -1 : 1;
    comparison = transition_slice_compare_size(first->component_index,
                                                second->component_index);
    if (comparison) return comparison;
    return transition_slice_compare_int((int)first->side, (int)second->side);
}

static int transition_slice_link_compare(const void* first_pointer,
                                         const void* second_pointer) {
    const alea_transition_slice_component_link_t* first = first_pointer;
    const alea_transition_slice_component_link_t* second = second_pointer;
    int comparison = transition_slice_compare_size(
        first->transition_component_index,
        second->transition_component_index);
    if (comparison) return comparison;
    comparison = transition_slice_compare_size(
        first->coverage_component_index, second->coverage_component_index);
    if (comparison) return comparison;
    return transition_slice_compare_int((int)first->boundary_sides,
                                        (int)second->boundary_sides);
}

static int transition_slice_same_ray(
    alea_transition_slice_orientation_t first_orientation,
    size_t first_ray, alea_transition_slice_orientation_t second_orientation,
    size_t second_ray) {
    return first_orientation == second_orientation && first_ray == second_ray;
}

static int transition_slice_ray_key_compare(
    alea_transition_slice_orientation_t first_orientation,
    size_t first_ray, alea_transition_slice_orientation_t second_orientation,
    size_t second_ray) {
    int comparison = transition_slice_compare_int(
        (int)first_orientation, (int)second_orientation);
    return comparison ? comparison :
        transition_slice_compare_size(first_ray, second_ray);
}

static int transition_slice_build_component_links(
    const alea_transition_slice_options_t* options,
    alea_transition_slice_result_t* result) {
    if (result->component_count == 0 ||
        result->coverage_component_count == 0)
        return 0;
    if (result->coverage_finding_count > SIZE_MAX / 2) return -1;
    const size_t boundary_count = result->finding_count;
    const size_t endpoint_count = result->coverage_finding_count * 2;
    size_t scratch_bytes = 0;
    if (transition_slice_add_bytes(
            &scratch_bytes, boundary_count,
            sizeof(transition_slice_component_boundary_ref_t)) ||
        transition_slice_add_bytes(
            &scratch_bytes, endpoint_count,
            sizeof(transition_slice_coverage_endpoint_ref_t)) ||
        transition_slice_add_bytes(
            &scratch_bytes, endpoint_count,
            sizeof(alea_transition_slice_component_link_t)))
        return -1;
    if (options->max_scratch_bytes &&
        scratch_bytes > options->max_scratch_bytes) {
        if (result->stats.stop_reason == ALEA_TRANSITION_SLICE_STOP_NONE)
            result->stats.stop_reason =
                ALEA_TRANSITION_SLICE_STOP_MAX_SCRATCH_BYTES;
        return 1;
    }
    if (scratch_bytes > result->stats.peak_scratch_bytes)
        result->stats.peak_scratch_bytes = scratch_bytes;

    transition_slice_component_boundary_ref_t* boundaries =
        malloc(boundary_count * sizeof(*boundaries));
    transition_slice_coverage_endpoint_ref_t* endpoints =
        malloc(endpoint_count * sizeof(*endpoints));
    alea_transition_slice_component_link_t* raw_links =
        malloc(endpoint_count * sizeof(*raw_links));
    if (!boundaries || !endpoints || !raw_links) {
        free(boundaries); free(endpoints); free(raw_links);
        return -1;
    }
    size_t boundary = 0;
    for (size_t component_index = 0;
         component_index < result->component_count; component_index++) {
        const alea_transition_slice_component_t* component =
            &result->components[component_index];
        const size_t end = component->first_finding_index +
            component->finding_count;
        for (size_t finding_index = component->first_finding_index;
             finding_index < end; finding_index++) {
            const alea_transition_slice_finding_t* finding =
                &result->findings[finding_index];
            boundaries[boundary++] =
                (transition_slice_component_boundary_ref_t){
                    .orientation = finding->orientation,
                    .ray_index = finding->ray_index,
                    .ray_t = finding->ray_t,
                    .component_index = component_index
                };
        }
    }
    size_t endpoint = 0;
    for (size_t component_index = 0;
         component_index < result->coverage_component_count;
         component_index++) {
        const alea_transition_slice_coverage_component_t* component =
            &result->coverage_components[component_index];
        const size_t end = component->first_finding_index +
            component->finding_count;
        for (size_t finding_index = component->first_finding_index;
             finding_index < end; finding_index++) {
            const alea_transition_slice_coverage_finding_t* finding =
                &result->coverage_findings[finding_index];
            endpoints[endpoint++] =
                (transition_slice_coverage_endpoint_ref_t){
                    .orientation = finding->orientation,
                    .ray_index = finding->ray_index,
                    .ray_t = finding->bracket_t_enter,
                    .component_index = component_index,
                    .side = ALEA_TRANSITION_SLICE_LINK_ENTER
                };
            endpoints[endpoint++] =
                (transition_slice_coverage_endpoint_ref_t){
                    .orientation = finding->orientation,
                    .ray_index = finding->ray_index,
                    .ray_t = finding->bracket_t_exit,
                    .component_index = component_index,
                    .side = ALEA_TRANSITION_SLICE_LINK_EXIT
                };
        }
    }
    qsort(boundaries, boundary, sizeof(*boundaries),
          transition_slice_boundary_ref_compare);
    qsort(endpoints, endpoint, sizeof(*endpoints),
          transition_slice_endpoint_ref_compare);

    size_t boundary_cursor = 0, raw_count = 0;
    for (size_t endpoint_index = 0; endpoint_index < endpoint;
         endpoint_index++) {
        const transition_slice_coverage_endpoint_ref_t* current =
            &endpoints[endpoint_index];
        while (boundary_cursor < boundary &&
               transition_slice_ray_key_compare(
                   boundaries[boundary_cursor].orientation,
                   boundaries[boundary_cursor].ray_index,
                   current->orientation, current->ray_index) < 0)
            boundary_cursor++;
        size_t candidate = boundary_cursor;
        while (candidate < boundary &&
               transition_slice_same_ray(
                   boundaries[candidate].orientation,
                   boundaries[candidate].ray_index,
                   current->orientation, current->ray_index)) {
            const double scale = 1.0 + fmax(fabs(current->ray_t),
                                            fabs(boundaries[candidate].ray_t));
            const double tolerance = 64.0 * RAY_EPSILON * scale;
            if (boundaries[candidate].ray_t < current->ray_t - tolerance) {
                candidate++;
                continue;
            }
            if (boundaries[candidate].ray_t > current->ray_t + tolerance)
                break;
            raw_links[raw_count++] =
                (alea_transition_slice_component_link_t){
                    .transition_component_index =
                        boundaries[candidate].component_index,
                    .coverage_component_index = current->component_index,
                    .boundary_sides = current->side,
                    .witness_pair_count = 1
                };
            break;
        }
        boundary_cursor = candidate;
    }
    free(boundaries);
    free(endpoints);
    if (raw_count == 0) {
        free(raw_links);
        return 0;
    }
    qsort(raw_links, raw_count, sizeof(*raw_links),
          transition_slice_link_compare);
    size_t required = 1;
    for (size_t i = 1; i < raw_count; i++)
        if (raw_links[i - 1].transition_component_index !=
                raw_links[i].transition_component_index ||
            raw_links[i - 1].coverage_component_index !=
                raw_links[i].coverage_component_index)
            required++;
    size_t capacity = required;
    if (options->max_component_links &&
        capacity > options->max_component_links)
        capacity = (size_t)options->max_component_links;
    if (options->max_output_bytes) {
        const size_t retained = transition_slice_retained_bytes(result);
        const size_t available = retained < options->max_output_bytes
            ? (size_t)(options->max_output_bytes - retained) : 0;
        const size_t byte_capacity =
            available / sizeof(*result->component_links);
        if (capacity > byte_capacity) capacity = byte_capacity;
    }
    if (capacity) {
        result->component_links =
            calloc(capacity, sizeof(*result->component_links));
        if (!result->component_links) {
            free(raw_links);
            return -1;
        }
        result->component_link_capacity = capacity;
    }
    size_t raw = 0;
    while (raw < raw_count && result->component_link_count < capacity) {
        size_t end = raw + 1;
        alea_transition_slice_component_link_t link = raw_links[raw];
        while (end < raw_count &&
               raw_links[end].transition_component_index ==
                   link.transition_component_index &&
               raw_links[end].coverage_component_index ==
                   link.coverage_component_index) {
            link.boundary_sides |= raw_links[end].boundary_sides;
            link.witness_pair_count += raw_links[end].witness_pair_count;
            end++;
        }
        result->component_links[result->component_link_count++] = link;
        raw = end;
    }
    free(raw_links);
    result->stats.component_links = result->component_link_count;
    result->stats.retained_output_bytes = transition_slice_retained_bytes(result);
    if (result->component_link_count < required &&
        result->stats.stop_reason == ALEA_TRANSITION_SLICE_STOP_NONE) {
        result->stats.stop_reason = options->max_component_links &&
                result->component_link_count >= options->max_component_links
            ? ALEA_TRANSITION_SLICE_STOP_MAX_COMPONENT_LINKS
            : ALEA_TRANSITION_SLICE_STOP_MAX_OUTPUT_BYTES;
        return 1;
    }
    return 0;
}

typedef struct {
    double uv_min[2];
    double uv_max[2];
    alea_transition_slice_tile_source_kind_t kind;
    size_t source_index;
    int priority;
} transition_slice_tile_seed_t;

static int transition_slice_tile_seed_compare(const void* first_pointer,
                                              const void* second_pointer) {
    const transition_slice_tile_seed_t* first = first_pointer;
    const transition_slice_tile_seed_t* second = second_pointer;
    int comparison = transition_slice_compare_int(first->priority,
                                                  second->priority);
    if (comparison) return comparison;
    for (int axis = 0; axis < 2; axis++) {
        if (first->uv_min[axis] != second->uv_min[axis])
            return first->uv_min[axis] < second->uv_min[axis] ? -1 : 1;
    }
    for (int axis = 0; axis < 2; axis++) {
        if (first->uv_max[axis] != second->uv_max[axis])
            return first->uv_max[axis] < second->uv_max[axis] ? -1 : 1;
    }
    comparison = transition_slice_compare_int((int)first->kind,
                                               (int)second->kind);
    return comparison ? comparison :
        transition_slice_compare_size(first->source_index,
                                      second->source_index);
}

static int transition_slice_rectangles_touch(
    const double first_min[2], const double first_max[2],
    const double second_min[2], const double second_max[2]) {
    return first_min[0] <= second_max[0] && second_min[0] <= first_max[0] &&
        first_min[1] <= second_max[1] && second_min[1] <= first_max[1];
}

static void transition_slice_critical_stop(
    alea_transition_slice_result_t* result,
    alea_transition_slice_critical_stop_reason_t reason) {
    if (result->stats.critical_stop_reason ==
        ALEA_TRANSITION_SLICE_CRITICAL_NONE)
        result->stats.critical_stop_reason = reason;
}

static int transition_slice_build_critical_tiles(
    const alea_slice_view_t* view,
    const alea_transition_slice_options_t* options,
    alea_transition_slice_result_t* result) {
    if (!options->enable_critical_refinement) return 0;
    size_t seed_count = result->component_count;
    if (seed_count > SIZE_MAX - result->coverage_component_count ||
        seed_count + result->coverage_component_count >
            SIZE_MAX - result->refinement_frontier_count)
        return -1;
    seed_count += result->coverage_component_count +
        result->refinement_frontier_count;
    result->stats.critical_tile_seeds = seed_count;
    if (seed_count == 0) return 0;
    size_t scratch_bytes = 0;
    if (transition_slice_add_bytes(&scratch_bytes, seed_count,
                                   sizeof(transition_slice_tile_seed_t)) ||
        transition_slice_add_bytes(&scratch_bytes, seed_count,
                                   sizeof(alea_transition_slice_critical_tile_t)))
        return -1;
    if (options->max_critical_scratch_bytes &&
        scratch_bytes > options->max_critical_scratch_bytes) {
        transition_slice_critical_stop(
            result, ALEA_TRANSITION_SLICE_CRITICAL_MAX_SCRATCH_BYTES);
        result->stats.omitted_critical_tile_sources = seed_count;
        return 0;
    }
    result->stats.peak_critical_scratch_bytes = scratch_bytes;
    if (result->stats.peak_scratch_bytes < scratch_bytes)
        result->stats.peak_scratch_bytes = scratch_bytes;
    transition_slice_tile_seed_t* seeds =
        calloc(seed_count, sizeof(*seeds));
    alea_transition_slice_critical_tile_t* merged =
        calloc(seed_count, sizeof(*merged));
    if (!seeds || !merged) { free(seeds); free(merged); return -1; }
    size_t seed = 0;
    for (size_t i = 0; i < result->component_count; i++) {
        const alea_transition_slice_component_t* source =
            &result->components[i];
        transition_slice_tile_seed_t* item = &seeds[seed++];
        memcpy(item->uv_min, source->uv_min, sizeof(item->uv_min));
        memcpy(item->uv_max, source->uv_max, sizeof(item->uv_max));
        item->kind = ALEA_TRANSITION_SLICE_TILE_SOURCE_TRANSITION_COMPONENT;
        item->source_index = i;
        item->priority = source->kind == ALEA_TRANSITION_TRUNCATED ||
                source->kind == ALEA_TRANSITION_UNRESOLVED ||
                source->kind == ALEA_TRANSITION_AMBIGUOUS_BOUNDARY
            ? 0 : 1;
    }
    for (size_t i = 0; i < result->coverage_component_count; i++) {
        const alea_transition_slice_coverage_component_t* source =
            &result->coverage_components[i];
        transition_slice_tile_seed_t* item = &seeds[seed++];
        memcpy(item->uv_min, source->uv_min, sizeof(item->uv_min));
        memcpy(item->uv_max, source->uv_max, sizeof(item->uv_max));
        const size_t end = source->first_finding_index + source->finding_count;
        for (size_t j = source->first_finding_index; j < end; j++) {
            const alea_transition_slice_coverage_finding_t* finding =
                &result->coverage_findings[j];
            const int axis = finding->orientation ==
                ALEA_TRANSITION_SLICE_HORIZONTAL ? 0 : 1;
            const double offset = axis == 0 ? view->u_min : view->v_min;
            const double low = offset + finding->bracket_t_enter;
            const double high = offset + finding->bracket_t_exit;
            if (item->uv_min[axis] > low) item->uv_min[axis] = low;
            if (item->uv_max[axis] < high) item->uv_max[axis] = high;
        }
        item->kind = ALEA_TRANSITION_SLICE_TILE_SOURCE_COVERAGE_COMPONENT;
        item->source_index = i;
        item->priority = source->truncated ||
            source->kind == ALEA_POINT_COVERAGE_UNRESOLVED ? 0 : 2;
    }
    for (size_t i = 0; i < result->refinement_frontier_count; i++) {
        const alea_transition_slice_refinement_frontier_t* source =
            &result->refinement_frontiers[i];
        transition_slice_tile_seed_t* item = &seeds[seed++];
        memcpy(item->uv_min, source->uv_min, sizeof(item->uv_min));
        memcpy(item->uv_max, source->uv_max, sizeof(item->uv_max));
        item->kind = ALEA_TRANSITION_SLICE_TILE_SOURCE_REFINEMENT_FRONTIER;
        item->source_index = i;
        item->priority = 0;
    }
    const double padding = options->critical_tile_padding > 0.0
        ? options->critical_tile_padding : 0.0;
    for (size_t i = 0; i < seed_count; i++) {
        seeds[i].uv_min[0] = fmax(view->u_min, seeds[i].uv_min[0] - padding);
        seeds[i].uv_max[0] = fmin(view->u_max, seeds[i].uv_max[0] + padding);
        seeds[i].uv_min[1] = fmax(view->v_min, seeds[i].uv_min[1] - padding);
        seeds[i].uv_max[1] = fmin(view->v_max, seeds[i].uv_max[1] + padding);
    }
    qsort(seeds, seed_count, sizeof(*seeds),
          transition_slice_tile_seed_compare);
    size_t selected = seed_count;
    if (options->max_critical_tile_sources &&
        selected > options->max_critical_tile_sources) {
        selected = options->max_critical_tile_sources;
        result->stats.omitted_critical_tile_sources += seed_count - selected;
        transition_slice_critical_stop(
            result, ALEA_TRANSITION_SLICE_CRITICAL_MAX_TILE_SOURCES);
    }
    size_t merged_count = 0;
    for (size_t i = 0; i < selected; i++) {
        size_t target = SIZE_MAX;
        for (size_t tile = 0; tile < merged_count; tile++) {
            if (transition_slice_rectangles_touch(
                    merged[tile].uv_min, merged[tile].uv_max,
                    seeds[i].uv_min, seeds[i].uv_max)) {
                target = tile;
                break;
            }
        }
        if (target == SIZE_MAX) {
            if (options->max_critical_tiles &&
                merged_count >= options->max_critical_tiles) {
                transition_slice_critical_stop(
                    result, ALEA_TRANSITION_SLICE_CRITICAL_MAX_TILES);
                continue;
            }
            target = merged_count++;
            memcpy(merged[target].uv_min, seeds[i].uv_min,
                   sizeof(merged[target].uv_min));
            memcpy(merged[target].uv_max, seeds[i].uv_max,
                   sizeof(merged[target].uv_max));
        } else {
            for (int axis = 0; axis < 2; axis++) {
                if (merged[target].uv_min[axis] > seeds[i].uv_min[axis])
                    merged[target].uv_min[axis] = seeds[i].uv_min[axis];
                if (merged[target].uv_max[axis] < seeds[i].uv_max[axis])
                    merged[target].uv_max[axis] = seeds[i].uv_max[axis];
            }
        }
        merged[target].source_flags |= 1u << (unsigned)seeds[i].kind;
        for (size_t tile = 0; tile < merged_count;) {
            if (tile == target || !transition_slice_rectangles_touch(
                    merged[target].uv_min, merged[target].uv_max,
                    merged[tile].uv_min, merged[tile].uv_max)) {
                tile++;
                continue;
            }
            for (int axis = 0; axis < 2; axis++) {
                if (merged[target].uv_min[axis] > merged[tile].uv_min[axis])
                    merged[target].uv_min[axis] = merged[tile].uv_min[axis];
                if (merged[target].uv_max[axis] < merged[tile].uv_max[axis])
                    merged[target].uv_max[axis] = merged[tile].uv_max[axis];
            }
            merged[target].source_flags |= merged[tile].source_flags;
            if (tile < target) target--;
            memmove(&merged[tile], &merged[tile + 1],
                    (merged_count - tile - 1) * sizeof(*merged));
            merged_count--;
            tile = 0;
        }
    }
    size_t retained = transition_slice_retained_bytes(result);
    size_t tile_capacity = merged_count;
    if (options->max_output_bytes) {
        size_t available = retained < options->max_output_bytes
            ? options->max_output_bytes - retained : 0;
        size_t cap = available / sizeof(*result->critical_tiles);
        if (tile_capacity > cap) tile_capacity = cap;
    }
    if (tile_capacity < merged_count)
        transition_slice_critical_stop(
            result, ALEA_TRANSITION_SLICE_CRITICAL_MAX_OUTPUT_BYTES);
    if (tile_capacity) {
        result->critical_tiles = calloc(
            tile_capacity, sizeof(*result->critical_tiles));
        if (!result->critical_tiles) { free(seeds); free(merged); return -1; }
        memcpy(result->critical_tiles, merged,
               tile_capacity * sizeof(*result->critical_tiles));
        result->critical_tile_count = tile_capacity;
    }
    retained = transition_slice_retained_bytes(result);
    size_t source_capacity = selected;
    if (options->max_output_bytes) {
        size_t available = retained < options->max_output_bytes
            ? options->max_output_bytes - retained : 0;
        size_t cap = available / sizeof(*result->critical_tile_sources);
        if (source_capacity > cap) source_capacity = cap;
    }
    if (source_capacity) {
        result->critical_tile_sources = calloc(
            source_capacity, sizeof(*result->critical_tile_sources));
        if (!result->critical_tile_sources) {
            free(seeds); free(merged); return -1;
        }
    }
    for (size_t tile = 0; tile < result->critical_tile_count; tile++) {
        result->critical_tiles[tile].first_source_index =
            result->critical_tile_source_count;
        for (size_t i = 0; i < selected; i++) {
            if (result->critical_tile_source_count >= source_capacity) break;
            if (seeds[i].uv_min[0] < merged[tile].uv_min[0] ||
                seeds[i].uv_max[0] > merged[tile].uv_max[0] ||
                seeds[i].uv_min[1] < merged[tile].uv_min[1] ||
                seeds[i].uv_max[1] > merged[tile].uv_max[1]) continue;
            result->critical_tile_sources[result->critical_tile_source_count++] =
                (alea_transition_slice_critical_tile_source_t){
                    .kind = seeds[i].kind,
                    .source_index = seeds[i].source_index
                };
            result->critical_tiles[tile].source_count++;
        }
    }
    if (result->critical_tile_source_count < selected) {
        result->stats.omitted_critical_tile_sources +=
            selected - result->critical_tile_source_count;
        transition_slice_critical_stop(
            result, ALEA_TRANSITION_SLICE_CRITICAL_MAX_OUTPUT_BYTES);
    }
    result->stats.critical_tiles = result->critical_tile_count;
    result->stats.critical_tile_sources = result->critical_tile_source_count;
    result->stats.retained_output_bytes = transition_slice_retained_bytes(result);
    free(seeds);
    free(merged);
    return 0;
}

static int transition_slice_scan_orientation(
    alea_system_t* sys, const alea_slice_view_t* view,
    const alea_transition_slice_options_t* options,
    alea_transition_slice_orientation_t orientation,
    size_t reserved_base_rays,
    alea_raycast_result_t* scratch,
    alea_ray_boundary_event_result_t* events,
    transition_slice_coverage_scratch_t* coverage_scratch,
    size_t coverage_scratch_bytes,
    alea_transition_slice_result_t* result) {
    const size_t base_count = orientation == ALEA_TRANSITION_SLICE_HORIZONTAL
        ? options->horizontal_rays : options->vertical_rays;
    if (base_count == 0) return 0;
    if (base_count > SIZE_MAX / sizeof(transition_slice_row_t)) return -1;
    size_t row_bytes = base_count * sizeof(transition_slice_row_t);
    size_t base_live_bytes = row_bytes > SIZE_MAX - coverage_scratch_bytes
        ? SIZE_MAX : row_bytes + coverage_scratch_bytes;
    if (options->max_row_scratch_bytes &&
        base_live_bytes > options->max_row_scratch_bytes) {
        result->stats.stop_reason =
            ALEA_TRANSITION_SLICE_STOP_MAX_SCRATCH_BYTES;
        return 1;
    }
    transition_slice_row_t* rows = calloc(base_count, sizeof(*rows));
    if (!rows) return -1;
    transition_slice_record_row_scratch(result, base_live_bytes);
    const double low = orientation == ALEA_TRANSITION_SLICE_HORIZONTAL
        ? view->v_min : view->u_min;
    const double span = orientation == ALEA_TRANSITION_SLICE_HORIZONTAL
        ? view->v_max - view->v_min : view->u_max - view->u_min;
    for (size_t base = 0; base < base_count; base++) {
        if (alea_interrupted()) {
            result->stats.stop_reason = ALEA_TRANSITION_SLICE_STOP_INTERRUPTED;
            free(rows);
            return 1;
        }
        const double coordinate = low + ((double)base + 0.5) * span /
            (double)base_count;
        const size_t ray_index = orientation == ALEA_TRANSITION_SLICE_HORIZONTAL
            ? result->stats.horizontal_rays_executed
            : result->stats.vertical_rays_executed;
        int rc = transition_slice_scan_ray(
            sys, view, options, orientation, ray_index, base, 0, coordinate,
            scratch, events, coverage_scratch, result, &rows[base]);
        transition_slice_record_row_scratch(result, base_live_bytes);
        if (rc != 0) { free(rows); return rc; }
    }
    if (options->max_refinement_depth == 0) {
        if (transition_slice_retain_differing_frontiers(
                view, options, orientation, rows, base_count, result) != 0) {
            free(rows);
            return -1;
        }
        free(rows);
        return 0;
    }

    size_t row_count = base_count;
    for (uint32_t depth = 0;; depth++) {
        size_t marked = 0, spacing_limited = 0;
        for (size_t row = 0; row + 1 < row_count; row++) {
            if (!transition_slice_rows_differ(
                    &rows[row], &rows[row + 1], options->refine_signals))
                continue;
            const double spacing = rows[row + 1].transverse_coordinate -
                rows[row].transverse_coordinate;
            if (options->min_transverse_spacing > 0.0 &&
                spacing < 2.0 * options->min_transverse_spacing) {
                spacing_limited++;
                continue;
            }
            marked++;
        }
        if (marked == 0) {
            if (spacing_limited &&
                transition_slice_retain_differing_frontiers(
                    view, options, orientation, rows, row_count,
                    result) != 0) {
                free(rows);
                return -1;
            }
            transition_slice_merge_status(
                result, spacing_limited
                    ? ALEA_TRANSITION_SLICE_REFINEMENT_MIN_SPACING
                    : ALEA_TRANSITION_SLICE_REFINEMENT_CONVERGED);
            free(rows);
            return 0;
        }
        if (depth >= options->max_refinement_depth) {
            if (transition_slice_retain_differing_frontiers(
                    view, options, orientation, rows, row_count,
                    result) != 0) {
                free(rows);
                return -1;
            }
            transition_slice_merge_status(
                result, ALEA_TRANSITION_SLICE_REFINEMENT_MAX_DEPTH);
            free(rows);
            return 0;
        }
        uint64_t available_rays = options->max_rays &&
                result->stats.executed_rays + reserved_base_rays <
                    options->max_rays
            ? options->max_rays - result->stats.executed_rays -
                reserved_base_rays
            : 0;
        if (options->max_rays && marked > available_rays) {
            if (transition_slice_retain_differing_frontiers(
                    view, options, orientation, rows, row_count,
                    result) != 0) {
                free(rows);
                return -1;
            }
            result->stats.stop_reason = ALEA_TRANSITION_SLICE_STOP_MAX_RAYS;
            free(rows);
            return 1;
        }
        if (marked > SIZE_MAX - row_count) { free(rows); return -1; }
        size_t next_count = row_count + marked;
        if (next_count > SIZE_MAX / sizeof(*rows)) { free(rows); return -1; }
        size_t next_bytes = next_count * sizeof(*rows);
        size_t live_bytes = row_bytes > SIZE_MAX - next_bytes
            ? SIZE_MAX : row_bytes + next_bytes;
        live_bytes = live_bytes > SIZE_MAX - coverage_scratch_bytes
            ? SIZE_MAX : live_bytes + coverage_scratch_bytes;
        if (options->max_row_scratch_bytes &&
            live_bytes > options->max_row_scratch_bytes) {
            if (transition_slice_retain_differing_frontiers(
                    view, options, orientation, rows, row_count,
                    result) != 0) {
                free(rows);
                return -1;
            }
            result->stats.stop_reason =
                ALEA_TRANSITION_SLICE_STOP_MAX_SCRATCH_BYTES;
            free(rows);
            return 1;
        }
        transition_slice_row_t* next = calloc(next_count, sizeof(*next));
        if (!next) { free(rows); return -1; }
        transition_slice_record_row_scratch(result, live_bytes);
        size_t output = 0;
        for (size_t row = 0; row < row_count; row++) {
            next[output++] = rows[row];
            if (row + 1 == row_count ||
                !transition_slice_rows_differ(
                    &rows[row], &rows[row + 1], options->refine_signals))
                continue;
            const double spacing = rows[row + 1].transverse_coordinate -
                rows[row].transverse_coordinate;
            if (options->min_transverse_spacing > 0.0 &&
                spacing < 2.0 * options->min_transverse_spacing)
                continue;
            const double coordinate = 0.5 *
                (rows[row].transverse_coordinate +
                 rows[row + 1].transverse_coordinate);
            const size_t ray_index =
                orientation == ALEA_TRANSITION_SLICE_HORIZONTAL
                    ? result->stats.horizontal_rays_executed
                    : result->stats.vertical_rays_executed;
            int rc = transition_slice_scan_ray(
                sys, view, options, orientation, ray_index, SIZE_MAX,
                depth + 1, coordinate, scratch, events, coverage_scratch,
                result, &next[output]);
            transition_slice_record_row_scratch(result, live_bytes);
            if (rc != 0) {
                free(next); free(rows); return rc;
            }
            result->stats.refined_rays_executed++;
            if (result->stats.max_refinement_depth_reached < depth + 1)
                result->stats.max_refinement_depth_reached = depth + 1;
            output++;
        }
        free(rows);
        rows = next;
        row_count = next_count;
        row_bytes = next_bytes;
    }
}

int alea_transition_slice_screen(
    alea_system_t* sys, const alea_slice_view_t* view,
    const alea_transition_slice_options_t* input,
    alea_transition_slice_result_t* result) {
    if (!sys || !view || !result) return -1;
    alea_transition_slice_options_t defaults, options;
    alea_transition_slice_options_init(&defaults);
    options = defaults;
    if (input) {
        if (input->struct_size < sizeof(input->struct_size)) return -1;
        size_t bytes = input->struct_size < sizeof(options)
            ? input->struct_size : sizeof(options);
        memcpy(&options, input, bytes);
    }
    if (!(view->u_max > view->u_min) || !(view->v_max > view->v_min) ||
        options.max_coverage_hits == 0 ||
        options.max_coverage_hits > 16384 ||
        options.horizontal_rays > SIZE_MAX - options.vertical_rays ||
        !isfinite(options.min_transverse_spacing) ||
        options.min_transverse_spacing < 0.0 ||
        !isfinite(options.critical_tile_padding) ||
        options.critical_tile_padding < 0.0 ||
        (options.refine_signals &
         ~(ALEA_TRANSITION_SLICE_REFINE_SIGNATURE |
           ALEA_TRANSITION_SLICE_REFINE_FINDING)))
        return -1;

    alea_transition_slice_result_t candidate = {0};
    candidate.stats.critical_enabled = options.enable_critical_refinement != 0;
    candidate.stats.critical_stop_reason = options.enable_critical_refinement
        ? ALEA_TRANSITION_SLICE_CRITICAL_NONE
        : ALEA_TRANSITION_SLICE_CRITICAL_DISABLED;
    candidate.stats.requested_rays =
        options.horizontal_rays + options.vertical_rays;
    candidate.stats.refinement_status = options.max_refinement_depth
        ? ALEA_TRANSITION_SLICE_REFINEMENT_CONVERGED
        : ALEA_TRANSITION_SLICE_REFINEMENT_NOT_REQUESTED;
    int scan_incomplete = 0;
    if (options.max_rays && candidate.stats.requested_rays > options.max_rays) {
        candidate.stats.stop_reason = ALEA_TRANSITION_SLICE_STOP_MAX_RAYS;
        scan_incomplete = 1;
    } else if (alea_raycast_ensure_hier_caches(sys) != 0) {
        return -1;
    } else {
        alea_raycast_result_t scratch;
        alea_ray_boundary_event_result_t events;
        transition_slice_coverage_scratch_t coverage_scratch;
        memset(&coverage_scratch, 0, sizeof(coverage_scratch));
        const int coverage_enabled =
            options.coverage_uniform_probes_per_ray != 0 ||
            options.coverage_probe_selected_intervals;
        size_t coverage_scratch_bytes = coverage_enabled
            ? transition_slice_coverage_scratch_bytes(
                options.max_coverage_hits) : 0;
        alea_raycast_result_init(&scratch);
        alea_ray_boundary_event_result_init(&events);
        int stop = 0, failed = 0;
        if (coverage_scratch_bytes == SIZE_MAX ||
            (options.max_row_scratch_bytes && coverage_enabled &&
             coverage_scratch_bytes > options.max_row_scratch_bytes)) {
            candidate.stats.stop_reason =
                ALEA_TRANSITION_SLICE_STOP_MAX_SCRATCH_BYTES;
            stop = 1;
            scan_incomplete = 1;
        } else if (coverage_enabled &&
                   transition_slice_coverage_scratch_init(
                       &coverage_scratch, options.max_coverage_hits) != 0) {
            failed = 1;
        }
        for (int orientation = ALEA_TRANSITION_SLICE_HORIZONTAL;
             orientation <= ALEA_TRANSITION_SLICE_VERTICAL && !stop && !failed;
             orientation++) {
            int rc = transition_slice_scan_orientation(
                sys, view, &options,
                (alea_transition_slice_orientation_t)orientation,
                orientation == ALEA_TRANSITION_SLICE_HORIZONTAL
                    ? options.vertical_rays : 0,
                &scratch, &events, &coverage_scratch,
                coverage_scratch_bytes, &candidate);
            if (rc < 0) failed = 1;
            if (rc > 0) { stop = 1; scan_incomplete = 1; }
            if (failed) break;
        }
        alea_ray_boundary_event_result_free(&events);
        alea_raycast_result_free(&scratch);
        transition_slice_coverage_scratch_free(&coverage_scratch);
        if (!failed) {
            int component_rc = transition_slice_build_components(
                &options, &candidate);
            if (component_rc < 0) failed = 1;
            if (component_rc > 0) stop = 1;
        }
        if (!failed) {
            int component_rc = transition_slice_build_coverage_components(
                &options, &candidate);
            if (component_rc < 0) failed = 1;
            if (component_rc > 0) stop = 1;
        }
        if (!failed) {
            int component_rc = transition_slice_build_component_links(
                &options, &candidate);
            if (component_rc < 0) failed = 1;
            if (component_rc > 0) stop = 1;
        }
        if (!failed && transition_slice_build_critical_tiles(
                view, &options, &candidate) != 0)
            failed = 1;
        if (!failed && options.enable_critical_refinement) {
            transition_slice_critical_sink_context_t critical_sink = {
                .result = &candidate, .options = &options
            };
            if (alea_transition_slice_enumerate_critical_tiles(
                    sys, view, &options, candidate.critical_tiles,
                    candidate.critical_tile_count,
                    transition_slice_append_critical_finding,
                    &critical_sink, &candidate.stats) != 0)
                failed = 1;
        }
        if (failed) {
            free(candidate.findings);
            free(candidate.components);
            free(candidate.coverage_findings);
            free(candidate.coverage_components);
            free(candidate.component_links);
            free(candidate.refinement_frontiers);
            free(candidate.critical_tiles);
            free(candidate.critical_tile_sources);
            free(candidate.critical_findings);
            return -1;
        }
    }
    candidate.stats.critical_complete = options.enable_critical_refinement &&
        candidate.stats.critical_stop_reason ==
            ALEA_TRANSITION_SLICE_CRITICAL_NONE;
    candidate.stats.complete =
        candidate.stats.stop_reason == ALEA_TRANSITION_SLICE_STOP_NONE &&
        (!options.enable_critical_refinement ||
         candidate.stats.critical_complete);
    if (scan_incomplete && options.max_refinement_depth) {
        candidate.stats.refinement_status =
            ALEA_TRANSITION_SLICE_REFINEMENT_STOPPED;
    }
    candidate.stats.converged = candidate.stats.refinement_status ==
            ALEA_TRANSITION_SLICE_REFINEMENT_CONVERGED;
    free(result->findings);
    free(result->components);
    free(result->coverage_findings);
    free(result->coverage_components);
    free(result->component_links);
    free(result->refinement_frontiers);
    free(result->critical_tiles);
    free(result->critical_tile_sources);
    free(result->critical_findings);
    *result = candidate;
    return 0;
}
