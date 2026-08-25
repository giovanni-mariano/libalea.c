// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_geo_validator.h"

#include "core/alea_eval.h"
#include "core/alea_system.h"
#include "raycast/ray_epsilon.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    alea_transition_result_t result;
} transition_probe_t;

void alea_transition_options_init(alea_transition_options_t* options) {
    if (!options) return;
    memset(options, 0, sizeof(*options));
    options->struct_size = sizeof(*options);
    options->probe_distance = SURFACE_SAMPLE_OFFSET;
    options->max_coverage_hits = 256;
}

const char* alea_transition_kind_name(alea_transition_kind_t kind) {
    switch (kind) {
        case ALEA_TRANSITION_VALID: return "valid";
        case ALEA_TRANSITION_GAP: return "gap";
        case ALEA_TRANSITION_OVERLAP: return "overlap";
        case ALEA_TRANSITION_UNDEFINED_FILL: return "undefined_fill";
        case ALEA_TRANSITION_MISSING_NEIGHBOR: return "missing_neighbor";
        case ALEA_TRANSITION_NON_ADJACENT: return "non_adjacent_transition";
        case ALEA_TRANSITION_AMBIGUOUS_NEIGHBOR: return "ambiguous_neighbor";
        case ALEA_TRANSITION_SURFACE_CHAIN_CORNER: return "surface_chain_corner";
        case ALEA_TRANSITION_AMBIGUOUS_BOUNDARY: return "ambiguous_boundary";
        case ALEA_TRANSITION_UNRESOLVED: return "unresolved";
        case ALEA_TRANSITION_TRUNCATED: return "truncated";
    }
    return "unresolved";
}

static uint32_t transition_surface_index(const alea_system_t* sys,
                                         int surface_id) {
    if (!sys || surface_id <= 0 || !sys->mc_id_to_surface ||
        (size_t)surface_id >= sys->mc_id_to_surface_size)
        return UINT32_MAX;
    return sys->mc_id_to_surface[surface_id];
}

static int transition_cell_contains(const alea_system_t* sys,
                                    uint32_t cell_index,
                                    const double point[3]) {
    if (!sys || cell_index >= alea_vec_count(&sys->cells)) return 0;
    const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
    if (cell->root_node_id == ALEA_NODE_ID_INVALID) return 0;
    return alea_contains_point(sys, cell->root_node_id,
                               point[0], point[1], point[2]) ? 1 : 0;
}

/* Returns -1/1 for one oriented reference, 0 when absent, and 2 when the
 * current cell uses both orientations of the exact card. */
static int transition_current_sense(const alea_system_t* sys,
                                    uint32_t surface_index,
                                    uint32_t current_cell_index) {
    if (!sys || !sys->surface_cell_offsets ||
        surface_index >= alea_vec_count(&sys->surfaces)) return 0;
    size_t begin = sys->surface_cell_offsets[surface_index];
    size_t end = sys->surface_cell_offsets[surface_index + 1];
    int sense = 0;
    for (size_t i = begin; i < end; i++) {
        const alea_surface_cell_ref_t* ref = &sys->surface_cell_refs[i];
        if (ref->cell_index != current_cell_index) continue;
        if (sense != 0 && sense != ref->sense) return 2;
        sense = ref->sense;
    }
    return sense;
}

static void transition_append_id(int id, int* values, size_t* count,
                                 uint32_t* flags, uint32_t truncated_flag) {
    for (size_t i = 0; i < *count; i++)
        if (values[i] == id) return;
    if (*count < ALEA_TRANSITION_EVIDENCE_CAPACITY) {
        values[(*count)++] = id;
    } else {
        *flags |= truncated_flag;
    }
}

static void transition_surface_candidates(
    const alea_system_t* sys,
    uint32_t surface_index,
    uint32_t current_cell_index,
    int universe_id,
    int current_sense,
    const double after_point[3],
    size_t* out_candidate_count,
    size_t* out_containing_count,
    int* out_first_containing,
    alea_transition_result_t* result) {
    *out_candidate_count = 0;
    *out_containing_count = 0;
    *out_first_containing = -1;
    if (!sys || current_sense == 0 || current_sense == 2 ||
        !sys->surface_cell_offsets ||
        surface_index >= alea_vec_count(&sys->surfaces)) return;

    size_t begin = sys->surface_cell_offsets[surface_index];
    size_t end = sys->surface_cell_offsets[surface_index + 1];
    for (size_t i = begin; i < end; i++) {
        const alea_surface_cell_ref_t* ref = &sys->surface_cell_refs[i];
        if (ref->cell_index == current_cell_index ||
            ref->cell_index >= alea_vec_count(&sys->cells) ||
            ref->sense == current_sense) continue;
        const alea_cell_entry_t* cell = &sys->cells.data[ref->cell_index];
        if (cell->universe_id != universe_id) continue;
        (*out_candidate_count)++;
        transition_append_id(cell->mc_cell_id,
                             result->candidate_cell_ids,
                             &result->candidate_cell_count,
                             &result->flags,
                             ALEA_TRANSITION_FLAG_CANDIDATES_TRUNCATED);
        if (transition_cell_contains(sys, ref->cell_index, after_point)) {
            if (*out_first_containing < 0)
                *out_first_containing = cell->mc_cell_id;
            (*out_containing_count)++;
        }
    }
}

static int transition_after_coverage(
    alea_system_t* sys, int universe_id, const double point[3],
    size_t capacity, alea_transition_result_t* result) {
    alea_cell_hit_t* hits = calloc(capacity, sizeof(*hits));
    uint64_t* keys = calloc(capacity, sizeof(*keys));
    uint64_t* parents = calloc(capacity, sizeof(*parents));
    uint8_t* owners = calloc(capacity, sizeof(*owners));
    if (!hits || !keys || !parents || !owners) {
        free(hits); free(keys); free(parents); free(owners);
        return -1;
    }
    int count = alea_find_all_cells_in_universe_coverage_chain(
        sys, universe_id, point[0], point[1], point[2],
        hits, keys, parents, capacity);
    if (count < 0) {
        free(hits); free(keys); free(parents); free(owners);
        return -1;
    }
    result->flags |= ALEA_TRANSITION_FLAG_COVERAGE_FALLBACK;
    result->coverage_fallbacks++;
    if ((size_t)count >= capacity) {
        result->flags |= ALEA_TRANSITION_FLAG_OWNERS_TRUNCATED;
        result->after_coverage_kind = ALEA_POINT_COVERAGE_UNRESOLVED;
        free(hits); free(keys); free(parents); free(owners);
        return 1;
    }

    alea_point_coverage_classification_t classification;
    if (alea_classify_point_coverage_chain(
            hits, keys, parents, (size_t)count, -1,
            owners, &classification) != 0) {
        free(hits); free(keys); free(parents); free(owners);
        return -1;
    }
    result->after_coverage_kind = classification.kind;
    result->after_owner_count = classification.owner_count;
    for (int i = 0; i < count; i++) {
        if (!owners[i]) continue;
        transition_append_id(hits[i].cell_id,
                             result->owner_cell_ids,
                             &result->owner_cell_count,
                             &result->flags,
                             ALEA_TRANSITION_FLAG_OWNERS_TRUNCATED);
        if (classification.owner_count == 1)
            result->after_cell_id = hits[i].cell_id;
    }
    free(hits); free(keys); free(parents); free(owners);
    return 0;
}

static int transition_one_surface_connects(
    const alea_system_t* sys,
    int surface_id,
    uint32_t current_cell_index, int universe_id,
    const double after_point[3], int after_cell_id,
    alea_transition_result_t* result) {
    if (surface_id <= 0 || surface_id == result->primary_surface_id) return 0;
    uint32_t surface_index = transition_surface_index(sys, surface_id);
    if (surface_index == UINT32_MAX) return 0;
    int sense = transition_current_sense(sys, surface_index, current_cell_index);
    size_t candidate_count = 0, containing_count = 0;
    int first = -1;
    transition_surface_candidates(
        sys, surface_index, current_cell_index, universe_id, sense,
        after_point, &candidate_count, &containing_count, &first, result);
    return containing_count == 1 && first == after_cell_id;
}

static int transition_tied_connects(
    const alea_system_t* sys,
    const int* tied_surface_ids, size_t tied_surface_count,
    uint32_t current_cell_index, int universe_id,
    const double before_point[3], const double after_point[3],
    int after_cell_id, alea_transition_result_t* result) {
    const alea_cell_entry_t* current = &sys->cells.data[current_cell_index];
    int connecting = -1;

    /* Discover simultaneous local cards directly from the current cell: an
     * oriented halfspace whose side changes across this probe is part of the
     * physical crossing group at the selected offset scale. */
    for (size_t i = 0; i < current->surface_index_count; i++) {
        uint32_t surface_index = current->surface_indices[i];
        if (surface_index >= alea_vec_count(&sys->surfaces)) continue;
        const alea_surface_entry_t* surface = &sys->surfaces.data[surface_index];
        if (surface->mc_surface_id == result->primary_surface_id) continue;
        int before_negative = alea_contains_point(
            sys, surface->neg_node,
            before_point[0], before_point[1], before_point[2]);
        int after_negative = alea_contains_point(
            sys, surface->neg_node,
            after_point[0], after_point[1], after_point[2]);
        if (before_negative == after_negative) continue;
        if (transition_one_surface_connects(
                sys, surface->mc_surface_id, current_cell_index, universe_id,
                after_point, after_cell_id, result) && connecting < 0)
            connecting = surface->mc_surface_id;
    }

    /* Event-provided ties supplement automatic discovery, particularly when
     * numerical coincidence is known by the ray walker at a tighter scale. */
    for (size_t i = 0; i < tied_surface_count; i++) {
        int surface_id = tied_surface_ids[i];
        if (transition_one_surface_connects(
                sys, surface_id, current_cell_index, universe_id,
                after_point, after_cell_id, result) && connecting < 0)
            connecting = surface_id;
    }
    return connecting;
}

static int transition_probe(
    alea_system_t* sys, int universe_id, uint32_t current_cell_index,
    int current_cell_id, uint32_t primary_surface_index,
    int primary_surface_id, const int* tied_surface_ids,
    size_t tied_surface_count, const double point[3], const double direction[3],
    double distance, size_t coverage_capacity, transition_probe_t* probe) {
    alea_transition_result_t* result = &probe->result;
    memset(result, 0, sizeof(*result));
    result->kind = ALEA_TRANSITION_UNRESOLVED;
    result->after_coverage_kind = ALEA_POINT_COVERAGE_UNRESOLVED;
    result->universe_id = universe_id;
    result->current_cell_id = current_cell_id;
    result->primary_surface_id = primary_surface_id;
    result->connecting_surface_id = -1;
    result->after_cell_id = -1;
    result->occurrence_depth = -1;
    result->probe_distance = distance;
    memcpy(result->crossing_point, point, sizeof(result->crossing_point));
    memcpy(result->direction, direction, sizeof(result->direction));
    for (int axis = 0; axis < 3; axis++) {
        result->before_point[axis] = point[axis] - distance * direction[axis];
        result->after_point[axis] = point[axis] + distance * direction[axis];
    }

    if (transition_cell_contains(sys, current_cell_index, result->before_point))
        result->flags |= ALEA_TRANSITION_FLAG_CURRENT_BEFORE_CONTAINS;
    if (transition_cell_contains(sys, current_cell_index, result->after_point))
        result->flags |= ALEA_TRANSITION_FLAG_CURRENT_AFTER_CONTAINS;

    int current_sense = transition_current_sense(
        sys, primary_surface_index, current_cell_index);
    result->current_sense = current_sense;
    if (!(result->flags & ALEA_TRANSITION_FLAG_CURRENT_BEFORE_CONTAINS) ||
        (result->flags & ALEA_TRANSITION_FLAG_CURRENT_AFTER_CONTAINS)) {
        result->kind = ALEA_TRANSITION_AMBIGUOUS_BOUNDARY;
        return 0;
    }
    /* A Boolean cell may reference the same exact surface card with both
     * orientations in different branches.  Card-level sense is then
     * ambiguous, but containment on the two probe sides still proves that
     * this concrete witness exits the cell.  Skip the exact-card fast path
     * (which requires one orientation) and let complete after-side coverage
     * classify the transition. */
    if (current_sense == 0)
        result->flags |= ALEA_TRANSITION_FLAG_PRIMARY_MISSING;

    int first_primary = -1;
    transition_surface_candidates(
        sys, primary_surface_index, current_cell_index, universe_id,
        current_sense, result->after_point,
        &result->primary_candidate_count,
        &result->primary_containing_count,
        &first_primary, result);
    if (result->primary_candidate_count == 0)
        result->flags |= ALEA_TRANSITION_FLAG_PRIMARY_MISSING;

    if (result->primary_containing_count == 1) {
        result->kind = ALEA_TRANSITION_VALID;
        result->after_cell_id = first_primary;
        return 0;
    }

    int coverage_rc = transition_after_coverage(
        sys, universe_id, result->after_point, coverage_capacity, result);
    if (coverage_rc < 0) return -1;
    if (coverage_rc > 0) {
        result->kind = ALEA_TRANSITION_TRUNCATED;
        return 0;
    }
    switch (result->after_coverage_kind) {
        case ALEA_POINT_COVERAGE_GAP:
            result->kind = ALEA_TRANSITION_GAP;
            return 0;
        case ALEA_POINT_COVERAGE_OVERLAP:
            result->kind = ALEA_TRANSITION_OVERLAP;
            return 0;
        case ALEA_POINT_COVERAGE_UNDEFINED_FILL:
            result->kind = ALEA_TRANSITION_UNDEFINED_FILL;
            return 0;
        case ALEA_POINT_COVERAGE_UNRESOLVED:
            result->kind = ALEA_TRANSITION_UNRESOLVED;
            return 0;
        case ALEA_POINT_COVERAGE_UNIQUE:
            break;
    }

    int tied = transition_tied_connects(
        sys, tied_surface_ids, tied_surface_count,
        current_cell_index, universe_id, result->before_point,
        result->after_point,
        result->after_cell_id, result);
    if (tied >= 0) {
        result->kind = ALEA_TRANSITION_SURFACE_CHAIN_CORNER;
        result->connecting_surface_id = tied;
        result->flags |= ALEA_TRANSITION_FLAG_TIED_SURFACE_CONNECTS;
    } else if (result->primary_containing_count > 1) {
        result->kind = ALEA_TRANSITION_AMBIGUOUS_NEIGHBOR;
    } else {
        result->kind = ALEA_TRANSITION_NON_ADJACENT;
    }
    return 0;
}

static int transition_probe_same(const transition_probe_t* a,
                                 const transition_probe_t* b) {
    return a->result.kind == b->result.kind &&
           a->result.after_cell_id == b->result.after_cell_id &&
           a->result.connecting_surface_id == b->result.connecting_surface_id &&
           a->result.primary_containing_count ==
               b->result.primary_containing_count &&
           a->result.after_owner_count == b->result.after_owner_count;
}

int alea_check_transition_local(
    alea_system_t* sys, int universe_id, int current_cell_id,
    int primary_surface_id, const int* tied_surface_ids,
    size_t tied_surface_count, const double point[3],
    const double direction_input[3], const alea_transition_options_t* input,
    alea_transition_result_t* result) {
    if (!sys || !point || !direction_input || !result ||
        (tied_surface_count > 0 && !tied_surface_ids)) return -1;
    if (!sys->cell_adjacency_built || !sys->surface_cell_offsets) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              "transition check: query acceleration is not prepared");
        return -1;
    }
    int current_cell_index = alea_find_cell_by_id(sys, current_cell_id);
    if (current_cell_index < 0 ||
        sys->cells.data[current_cell_index].universe_id != universe_id) {
        alea_set_error_detail(ALEA_ERR_NOT_FOUND,
                              "transition check: cell %d is not in universe %d",
                              current_cell_id, universe_id);
        return -1;
    }
    uint32_t surface_index = transition_surface_index(sys, primary_surface_id);
    if (surface_index == UINT32_MAX) {
        alea_set_error_detail(ALEA_ERR_NOT_FOUND,
                              "transition check: surface %d not found",
                              primary_surface_id);
        return -1;
    }
    double norm = sqrt(direction_input[0] * direction_input[0] +
                       direction_input[1] * direction_input[1] +
                       direction_input[2] * direction_input[2]);
    if (!(norm > 0.0) || !isfinite(norm)) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "transition check: direction must be finite and non-zero");
        return -1;
    }
    double direction[3] = {
        direction_input[0] / norm,
        direction_input[1] / norm,
        direction_input[2] / norm
    };

    alea_transition_options_t options;
    alea_transition_options_init(&options);
    if (input) {
        if (input->probe_distance > 0.0)
            options.probe_distance = input->probe_distance;
        if (input->max_probe_distance > 0.0)
            options.max_probe_distance = input->max_probe_distance;
        if (input->max_coverage_hits > 0)
            options.max_coverage_hits = input->max_coverage_hits;
        if (input->max_coverage_fallbacks > 0)
            options.max_coverage_fallbacks = input->max_coverage_fallbacks;
    }
    if (options.max_probe_distance <= 0.0)
        options.max_probe_distance = 8.0 * options.probe_distance;
    if (!isfinite(options.probe_distance) || options.probe_distance <= 0.0 ||
        !isfinite(options.max_probe_distance) ||
        options.max_probe_distance <= 0.0 ||
        options.max_coverage_hits == 0 || options.max_coverage_hits > 16384) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "transition check: invalid resource or offset options");
        return -1;
    }

    double distances[3] = {
        options.probe_distance / 8.0,
        options.probe_distance,
        options.probe_distance * 8.0
    };
    size_t distance_count = 0;
    for (size_t i = 0; i < 3; i++) {
        double distance = distances[i];
        if (distance > options.max_probe_distance) continue;
        if (distance_count > 0 &&
            fabs(distance - distances[distance_count - 1]) <=
                1e-15 * fmax(1.0, distance)) continue;
        distances[distance_count++] = distance;
    }
    if (distance_count == 0) distances[distance_count++] =
        options.max_probe_distance;

    transition_probe_t previous, current;
    int have_previous = 0;
    size_t total_fallbacks = 0;
    for (size_t i = 0; i < distance_count; i++) {
        if (options.max_coverage_fallbacks != 0 &&
            total_fallbacks >= options.max_coverage_fallbacks) {
            *result = previous.result;
            result->kind = ALEA_TRANSITION_TRUNCATED;
            result->coverage_fallbacks = total_fallbacks;
            result->offset_attempts = i;
            return 0;
        }
        if (transition_probe(
                sys, universe_id, (uint32_t)current_cell_index,
                current_cell_id, surface_index, primary_surface_id,
                tied_surface_ids, tied_surface_count, point, direction,
                distances[i], options.max_coverage_hits, &current) != 0)
            return -1;
        current.result.offset_attempts = i + 1;
        total_fallbacks += current.result.coverage_fallbacks;
        if (have_previous && transition_probe_same(&previous, &current)) {
            *result = current.result;
            result->coverage_fallbacks = total_fallbacks;
            result->flags |= ALEA_TRANSITION_FLAG_OFFSET_STABLE;
            return 0;
        }
        previous = current;
        have_previous = 1;
    }

    *result = previous.result;
    result->coverage_fallbacks = total_fallbacks;
    result->offset_attempts = distance_count;
    if (distance_count == 1) {
        result->flags |= ALEA_TRANSITION_FLAG_OFFSET_STABLE;
    } else {
        result->kind = ALEA_TRANSITION_AMBIGUOUS_BOUNDARY;
    }
    return 0;
}
