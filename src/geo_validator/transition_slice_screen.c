// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_geo_validator.h"
#include "transition_slice_critical.h"
#include "transition_validation.h"

#include "raycast/raycast.h"
#include "raycast/ray_epsilon.h"
#include "core/alea_system.h"
#include "core/alea_spatial_hier.h"
#include "core/alea_universe.h"
#include "primitives/bbox.h"
#include "primitives/primitive_desc.h"
#include "util/alea_atomic.h"
#include "util/alea_parallel.h"
#include "util/compat.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    size_t cell_index;
    alea_bbox_t box;
} slice_error_index_cell_t;

struct alea_slice_error_query {
    alea_system_t* sys;
    alea_slice_error_query_options_t options;
    uint64_t query_id;
    uint64_t geometry_generation;
    size_t page_count;
    slice_error_index_cell_t* indexed_cells;
    size_t indexed_count;
    size_t* uncertain_cells;
    size_t uncertain_count;
    size_t uncached_cell_begin;
    size_t index_bytes;
    int has_hierarchy;
};

struct alea_slice_error_page {
    alea_slice_error_page_receipt_t receipt;
    alea_transition_slice_critical_finding_t* findings;
    size_t finding_count;
    size_t finding_capacity;
    alea_slice_error_witness_t* witnesses;
    size_t witness_count;
    size_t witness_capacity;
    alea_slice_error_interval_t* intervals;
    size_t interval_count;
    alea_slice_error_circle_t* circles;
    size_t circle_count;
    alea_slice_error_region_t* regions;
    size_t region_count;
    int populated;
};

static size_t slice_error_page_output_bytes(
    const alea_slice_error_page_t* page);

static atomic_uint_fast64_t slice_error_next_query_id = 1;
static int slice_error_build_index(alea_slice_error_query_t* query);
static int slice_error_slice_axis_for_world(const alea_slice_view_t* view,
                                           int world_axis);
static double slice_error_slice_axis_sign(const alea_slice_view_t* view,
                                          int world_axis);

void alea_slice_error_query_options_init(
    alea_slice_error_query_options_t* options) {
    if (!options) return;
    memset(options, 0, sizeof(*options));
    options->struct_size = sizeof(*options);
    options->tile_columns = 1;
    options->tile_rows = 1;
    options->max_index_bytes = 32u * 1024u * 1024u;
    alea_transition_slice_options_init(&options->scan_options);
    /* A small deterministic interior grid finds defects that do not intersect
     * a selected boundary transition. The existing option remains zero for
     * the general transition screen unless its caller explicitly requests it. */
    options->scan_options.coverage_uniform_probes_per_ray = 3;
    options->scan_options.occurrence_discovery =
        ALEA_TRANSITION_SLICE_OCCURRENCE_EXHAUSTIVE;
}

static int slice_error_finite_view(const alea_slice_view_t* view) {
    for (size_t i = 0; i < 3; ++i) {
        if (!isfinite(view->plane.origin[i]) ||
            !isfinite(view->plane.normal[i]) ||
            !isfinite(view->plane.u_axis[i]) ||
            !isfinite(view->plane.v_axis[i])) return 0;
    }
    return isfinite(view->u_min) && isfinite(view->u_max) &&
        isfinite(view->v_min) && isfinite(view->v_max) &&
        view->u_min < view->u_max && view->v_min < view->v_max;
}

alea_slice_error_query_t* alea_slice_error_query_create(
    alea_system_t* sys, const alea_slice_error_query_options_t* input) {
    if (!sys || !input || input->struct_size <
            offsetof(alea_slice_error_query_options_t, scan_options))
        return NULL;
    alea_slice_error_query_options_t options;
    alea_slice_error_query_options_init(&options);
    const size_t copy_size = input->struct_size < sizeof(options)
        ? input->struct_size : sizeof(options);
    memcpy(&options, input, copy_size);
    if (!options.max_index_bytes)
        options.max_index_bytes = 32u * 1024u * 1024u;
    const alea_slice_view_t* view = &options.view;
    if (!slice_error_finite_view(view) ||
        !isfinite(options.required_uv_min[0]) ||
        !isfinite(options.required_uv_min[1]) ||
        !isfinite(options.required_uv_max[0]) ||
        !isfinite(options.required_uv_max[1]) ||
        !(options.required_uv_max[0] > options.required_uv_min[0]) ||
        !(options.required_uv_max[1] > options.required_uv_min[1]) ||
        !isfinite(options.required_uv_max[0] -
                  options.required_uv_min[0]) ||
        !isfinite(options.required_uv_max[1] -
                  options.required_uv_min[1]) ||
        !isfinite(options.scan_options.critical_relative_distance_tolerance) ||
        options.scan_options.critical_relative_distance_tolerance < 0.0 ||
        options.scan_options.coverage_uniform_probes_per_ray > 64 ||
        !options.tile_columns || !options.tile_rows ||
        options.tile_columns > SIZE_MAX / options.tile_rows)
        return NULL;
    /* A candidate scan may use the existing critical machinery, but its
     * sampled discovery mode cannot establish exhaustive occurrence coverage. */
    options.scan_options.occurrence_discovery =
        ALEA_TRANSITION_SLICE_OCCURRENCE_EXHAUSTIVE;
    if (!options.scan_options.max_curves_per_tile ||
        !options.scan_options.max_critical_points ||
        !options.scan_options.max_coverage_hits ||
        !options.scan_options.max_exhaustive_occurrence_hits ||
        !options.scan_options.max_critical_scratch_bytes ||
        !options.scan_options.max_critical_findings ||
        !options.scan_options.max_output_bytes)
        return NULL;
    alea_slice_error_query_t* query = calloc(1, sizeof(*query));
    if (!query) return NULL;
    query->sys = sys;
    query->options = options;
    query->query_id = atomic_fetch_add(&slice_error_next_query_id, 1);
    query->geometry_generation = alea_system_geometry_generation(sys);
    query->page_count = options.tile_columns * options.tile_rows;
    if (slice_error_build_index(query) != 0) {
        alea_slice_error_query_destroy(query);
        return NULL;
    }
    return query;
}

void alea_slice_error_query_destroy(alea_slice_error_query_t* query) {
    if (!query) return;
    free(query->indexed_cells);
    free(query->uncertain_cells);
    free(query);
}

size_t alea_slice_error_query_page_count(const alea_slice_error_query_t* query) {
    return query ? query->page_count : 0;
}

uint64_t alea_slice_error_query_id(const alea_slice_error_query_t* query) {
    return query ? query->query_id : 0;
}

alea_slice_error_page_t* alea_slice_error_page_create(void) {
    return calloc(1, sizeof(alea_slice_error_page_t));
}

void alea_slice_error_page_destroy(alea_slice_error_page_t* page) {
    if (!page) return;
    free(page->findings);
    free(page->witnesses);
    free(page->intervals);
    free(page->circles);
    free(page->regions);
    free(page);
}

int alea_slice_error_page_receipt(const alea_slice_error_page_t* page,
                                  alea_slice_error_page_receipt_t* out_receipt) {
    if (!page || !out_receipt || !page->populated) return -1;
    *out_receipt = page->receipt;
    return 0;
}

size_t alea_slice_error_page_context_finding_count(
    const alea_slice_error_page_t* page) {
    return page && page->populated ? page->finding_count : 0;
}

int alea_slice_error_page_context_finding_get(
    const alea_slice_error_page_t* page, size_t index,
    alea_transition_slice_critical_finding_t* out_finding) {
    if (!page || !page->populated || !out_finding ||
        index >= page->finding_count) return -1;
    *out_finding = page->findings[index];
    return 0;
}

size_t alea_slice_error_page_witness_count(
    const alea_slice_error_page_t* page) {
    return page && page->populated ? page->witness_count : 0;
}

int alea_slice_error_page_witness_get(
    const alea_slice_error_page_t* page, size_t index,
    alea_slice_error_witness_t* out_witness) {
    if (!page || !page->populated || !out_witness ||
        index >= page->witness_count) return -1;
    *out_witness = page->witnesses[index];
    return 0;
}

size_t alea_slice_error_page_interval_count(const alea_slice_error_page_t* page) {
    return page && page->populated ? page->interval_count : 0;
}

int alea_slice_error_page_interval_get(const alea_slice_error_page_t* page,
                                       size_t index,
                                       alea_slice_error_interval_t* out_interval) {
    if (!page || !page->populated || !out_interval ||
        index >= page->interval_count) return -1;
    *out_interval = page->intervals[index];
    return 0;
}

size_t alea_slice_error_page_circle_count(const alea_slice_error_page_t* page) {
    return page && page->populated ? page->circle_count : 0;
}

int alea_slice_error_page_circle_get(const alea_slice_error_page_t* page,
                                     size_t index,
                                     alea_slice_error_circle_t* out_circle) {
    if (!page || !page->populated || !out_circle ||
        index >= page->circle_count) return -1;
    *out_circle = page->circles[index];
    return 0;
}

size_t alea_slice_error_page_region_count(const alea_slice_error_page_t* page) {
    return page && page->populated ? page->region_count : 0;
}

int alea_slice_error_page_region_get(const alea_slice_error_page_t* page,
                                     size_t index,
                                     alea_slice_error_region_t* out_region) {
    if (!page || !page->populated || !out_region ||
        index >= page->region_count) return -1;
    *out_region = page->regions[index];
    return 0;
}

typedef struct {
    alea_slice_error_page_t* page;
    size_t max_findings;
    size_t max_output_bytes;
    size_t omitted_findings;
} slice_error_finding_sink_t;

static int slice_error_retain_finding(
    const alea_transition_slice_critical_finding_t* finding, void* userdata) {
    slice_error_finding_sink_t* sink = userdata;
    alea_slice_error_page_t* page = sink->page;
    const size_t retained = slice_error_page_output_bytes(page);
    if (page->finding_count >= sink->max_findings ||
        retained > sink->max_output_bytes ||
        sizeof(*page->findings) > sink->max_output_bytes - retained) {
        if (sink->omitted_findings < SIZE_MAX) sink->omitted_findings++;
        return 1;
    }
    if (page->finding_count == page->finding_capacity) {
        size_t capacity = page->finding_capacity
            ? (page->finding_capacity > SIZE_MAX / 2u
                ? SIZE_MAX : page->finding_capacity * 2u) : 4u;
        if (capacity > sink->max_findings) capacity = sink->max_findings;
        const size_t byte_capacity =
            sink->max_output_bytes / sizeof(*page->findings);
        if (capacity > byte_capacity) capacity = byte_capacity;
        void* next = realloc(page->findings,
                             capacity * sizeof(*page->findings));
        if (!next) return -1;
        page->findings = next;
        page->finding_capacity = capacity;
    }
    page->findings[page->finding_count++] = *finding;
    return 0;
}

typedef struct {
    alea_slice_error_page_t* page;
    uint8_t* localized_tiles;
    size_t tile_count;
    size_t max_output_bytes;
    size_t omitted_regions;
} slice_error_numerical_region_sink_t;

static int slice_error_retain_numerical_region(
    const alea_transition_slice_numerical_region_t* input, void* userdata) {
    slice_error_numerical_region_sink_t* sink = userdata;
    if (!input || !sink || input->tile_index >= sink->tile_count) return -1;
    alea_slice_error_page_t* page = sink->page;
    const size_t search_begin = page->region_count > 256u
        ? page->region_count - 256u : 0u;
    for (size_t i = page->region_count; i-- > search_begin;) {
        alea_slice_error_region_t* existing = &page->regions[i];
        if (existing->kind != ALEA_POINT_COVERAGE_UNRESOLVED ||
            existing->owner_count != 0 ||
            existing->numerical_cause != input->cause ||
            existing->source_cell_id != input->source_cell_id ||
            existing->source_surface_id != input->source_surface_id ||
            existing->source_primitive_id != input->source_primitive_id ||
            existing->source_occurrence_key !=
                input->source_occurrence_key ||
            existing->source_universe_occurrence_key !=
                input->source_universe_occurrence_key) continue;
        int existing_contains = 1;
        int input_contains = 1;
        for (int axis = 0; axis < 2; ++axis) {
            const double tolerance = fmax(
                input->uncertainty,
                fmax(existing->uv_min_uncertainty[axis],
                     existing->uv_max_uncertainty[axis]));
            existing_contains &=
                existing->uv_min[axis] <= input->uv_min[axis] + tolerance &&
                existing->uv_max[axis] >= input->uv_max[axis] - tolerance;
            input_contains &=
                input->uv_min[axis] <= existing->uv_min[axis] + tolerance &&
                input->uv_max[axis] >= existing->uv_max[axis] - tolerance;
        }
        if (!existing_contains && !input_contains) continue;
        if (input_contains && !existing_contains) {
            memcpy(existing->uv_min, input->uv_min,
                   sizeof(existing->uv_min));
            memcpy(existing->uv_max, input->uv_max,
                   sizeof(existing->uv_max));
            for (int axis = 0; axis < 2; ++axis) {
                existing->uv_min_uncertainty[axis] = input->uncertainty;
                existing->uv_max_uncertainty[axis] = input->uncertainty;
            }
        }
        sink->localized_tiles[input->tile_index] = 1;
        return 0;
    }
    const size_t retained = slice_error_page_output_bytes(page);
    if (retained > sink->max_output_bytes ||
        sizeof(*page->regions) > sink->max_output_bytes - retained ||
        page->region_count >= SIZE_MAX / sizeof(*page->regions)) {
        if (sink->omitted_regions < SIZE_MAX) sink->omitted_regions++;
        return 1;
    }
    void* next = realloc(
        page->regions, (page->region_count + 1u) * sizeof(*page->regions));
    if (!next) return -1;
    page->regions = next;
    alea_slice_error_region_t* region =
        &page->regions[page->region_count++];
    memset(region, 0, sizeof(*region));
    memcpy(region->uv_min, input->uv_min, sizeof(region->uv_min));
    memcpy(region->uv_max, input->uv_max, sizeof(region->uv_max));
    for (int axis = 0; axis < 2; ++axis) {
        region->uv_min_uncertainty[axis] = input->uncertainty;
        region->uv_max_uncertainty[axis] = input->uncertainty;
    }
    region->kind = ALEA_POINT_COVERAGE_UNRESOLVED;
    region->numerical_cause = input->cause;
    region->source_cell_id = input->source_cell_id;
    region->source_surface_id = input->source_surface_id;
    region->source_primitive_id = input->source_primitive_id;
    region->source_occurrence_key = input->source_occurrence_key;
    region->source_universe_occurrence_key =
        input->source_universe_occurrence_key;
    sink->localized_tiles[input->tile_index] = 1;
    return 0;
}

typedef struct {
    alea_cell_hit_t* hits;
    uint64_t* occurrence_keys;
    uint64_t* parent_occurrence_keys;
    uint8_t* owner_mask;
    size_t capacity;
    size_t attempts;
    size_t failures;
    size_t confirmed_gaps;
    size_t confirmed_overlaps;
    size_t omitted;
} slice_error_confirmation_t;

static void slice_error_confirmation_free(slice_error_confirmation_t* scratch) {
    if (!scratch) return;
    free(scratch->hits);
    free(scratch->occurrence_keys);
    free(scratch->parent_occurrence_keys);
    free(scratch->owner_mask);
    memset(scratch, 0, sizeof(*scratch));
}

static int slice_error_confirmation_init(
    slice_error_confirmation_t* scratch, size_t capacity) {
    if (!scratch || !capacity ||
        capacity > SIZE_MAX / sizeof(*scratch->hits) ||
        capacity > SIZE_MAX / sizeof(*scratch->occurrence_keys) ||
        capacity > SIZE_MAX / sizeof(*scratch->parent_occurrence_keys) ||
        capacity > SIZE_MAX / sizeof(*scratch->owner_mask)) return -1;
    memset(scratch, 0, sizeof(*scratch));
    scratch->hits = calloc(capacity, sizeof(*scratch->hits));
    scratch->occurrence_keys = calloc(
        capacity, sizeof(*scratch->occurrence_keys));
    scratch->parent_occurrence_keys = calloc(
        capacity, sizeof(*scratch->parent_occurrence_keys));
    scratch->owner_mask = calloc(capacity, sizeof(*scratch->owner_mask));
    if (!scratch->hits || !scratch->occurrence_keys ||
        !scratch->parent_occurrence_keys || !scratch->owner_mask) {
        slice_error_confirmation_free(scratch);
        return -1;
    }
    scratch->capacity = capacity;
    return 0;
}

static int slice_error_append_witness(
    alea_slice_error_page_t* page, const alea_slice_error_witness_t* witness,
    size_t output_limit, slice_error_confirmation_t* scratch) {
    const size_t retained = slice_error_page_output_bytes(page);
    if (retained > output_limit ||
        sizeof(*page->witnesses) > output_limit - retained ||
        page->witness_count >= SIZE_MAX / sizeof(*page->witnesses)) {
        if (scratch->omitted < SIZE_MAX) scratch->omitted++;
        return 0;
    }
    if (page->witness_count == page->witness_capacity) {
        size_t capacity = page->witness_capacity
            ? (page->witness_capacity > SIZE_MAX / 2u
                ? SIZE_MAX : page->witness_capacity * 2u) : 8u;
        if (capacity > SIZE_MAX / sizeof(*page->witnesses)) {
            if (scratch->failures < SIZE_MAX) scratch->failures++;
            if (scratch->omitted < SIZE_MAX) scratch->omitted++;
            return 0;
        }
        const size_t byte_capacity =
            output_limit / sizeof(*page->witnesses);
        if (capacity > byte_capacity) capacity = byte_capacity;
        if (capacity <= page->witness_capacity) {
            if (scratch->omitted < SIZE_MAX) scratch->omitted++;
            return 0;
        }
        void* next = realloc(
            page->witnesses, capacity * sizeof(*page->witnesses));
        if (!next) {
            if (scratch->failures < SIZE_MAX) scratch->failures++;
            if (scratch->omitted < SIZE_MAX) scratch->omitted++;
            return 0;
        }
        page->witnesses = next;
        page->witness_capacity = capacity;
    }
    page->witnesses[page->witness_count++] = *witness;
    if (witness->kind == ALEA_POINT_COVERAGE_GAP)
        scratch->confirmed_gaps++;
    else if (witness->kind == ALEA_POINT_COVERAGE_OVERLAP)
        scratch->confirmed_overlaps++;
    return 1;
}

static int slice_error_confirm_point(
    const alea_slice_error_query_t* query,
    alea_slice_error_page_t* page,
    slice_error_confirmation_t* scratch,
    const double uv[2], const double world[3],
    alea_slice_error_witness_source_t source,
    const alea_transition_slice_critical_finding_t* finding) {
    scratch->attempts++;
    if (!scratch->capacity) {
        scratch->failures++;
        return 0;
    }
    int hit_count =
        alea_find_all_cells_at_point_coverage_chain_recursive(
            query->sys, world[0], world[1], world[2], scratch->hits,
            scratch->occurrence_keys, scratch->parent_occurrence_keys,
            scratch->capacity);
    if (hit_count < 0 && alea_vec_count(&query->sys->cells) == 0)
        hit_count = 0;
    if (hit_count < 0 || (size_t)hit_count >= scratch->capacity) {
        scratch->failures++;
        return 0;
    }
    memset(scratch->owner_mask, 0,
           scratch->capacity * sizeof(*scratch->owner_mask));
    alea_point_coverage_classification_t classification;
    if (alea_classify_point_coverage_chain(
            scratch->hits, scratch->occurrence_keys,
            scratch->parent_occurrence_keys, (size_t)hit_count, -1,
            scratch->owner_mask, &classification) != 0) {
        scratch->failures++;
        return 0;
    }

    /* A selected fill container followed by no child is an actual gap when
     * the referenced universe exists.  Keep a genuinely missing universe as
     * UNDEFINED_FILL: that is malformed hierarchy metadata rather than proof
     * of empty geometric coverage. */
    if (classification.kind == ALEA_POINT_COVERAGE_UNDEFINED_FILL) {
        int container = -1;
        for (int i = 0; i < hit_count; ++i) {
            if (scratch->owner_mask[i]) {
                container = i;
                break;
            }
        }
        if (container >= 0 && scratch->hits[container].fill_universe > 0 &&
            alea_get_universe(query->sys,
                              scratch->hits[container].fill_universe)) {
            classification.kind = ALEA_POINT_COVERAGE_GAP;
            classification.target_depth = scratch->hits[container].depth + 1;
            classification.owner_count = 0;
            memset(scratch->owner_mask, 0,
                   scratch->capacity * sizeof(*scratch->owner_mask));
        }
    }
    if (classification.kind != ALEA_POINT_COVERAGE_GAP &&
        classification.kind != ALEA_POINT_COVERAGE_OVERLAP) return 0;
    if (classification.kind == ALEA_POINT_COVERAGE_GAP &&
        classification.target_depth < 0)
        classification.target_depth = 0;

    alea_slice_error_witness_t witness;
    memset(&witness, 0, sizeof(witness));
    witness.evidence_scope = ALEA_SLICE_BOUNDARY_EVIDENCE_VERIFIED_POINT;
    witness.kind = classification.kind;
    witness.source = source;
    memcpy(witness.uv, uv, sizeof(witness.uv));
    memcpy(witness.world_point, world, sizeof(witness.world_point));
    witness.target_depth = classification.target_depth;
    witness.owner_count_lower_bound = classification.owner_count;
    for (int i = 0; i < hit_count; ++i) {
        if (!scratch->owner_mask[i]) continue;
        if (witness.owner_count == ALEA_SLICE_ERROR_OWNER_CAPACITY) continue;
        const size_t owner = witness.owner_count++;
        witness.owner_cell_ids[owner] = scratch->hits[i].cell_id;
        witness.owner_universe_ids[owner] = scratch->hits[i].universe_id;
        witness.owner_depths[owner] = scratch->hits[i].depth;
        witness.owner_occurrence_keys[owner] = scratch->occurrence_keys[i];
        witness.owner_parent_occurrence_keys[owner] =
            scratch->parent_occurrence_keys[i];
    }
    witness.owners_complete =
        witness.owner_count == classification.owner_count;
    if (finding) {
        witness.source_cell_id = finding->source_cell_id;
        witness.source_surface_id = finding->source_surface_id;
        witness.source_occurrence_key = finding->source_occurrence_key;
        witness.source_universe_occurrence_key =
            finding->source_universe_occurrence_key;
    }
    return slice_error_append_witness(
        page, &witness, query->options.scan_options.max_output_bytes, scratch);
}

typedef struct {
    int axis;
    int coefficient_sign;
    double coordinate;
    double uncertainty;
    int surface_id;
    uint32_t primitive_id;
} slice_error_axis_plane_t;

/* The analytic CSG bbox is conservative for these operations. Reject
 * complement and malformed trees before using it to exclude a cell. */
static int slice_error_bbox_tree_safe(const alea_system_t* sys,
                                     alea_node_id_t id, size_t depth,
                                     size_t* work) {
    if (!*work || depth > 128 || id >= sys->nodes.count) return 0;
    --*work;
    const alea_node_t* node = &sys->nodes.data[id];
    const alea_operation_t op = ALEA_GET_OPERATION(node);
    if (op == ALEA_OP_PRIMITIVE) {
        const uint32_t id = node->primitive.primitive_id;
        if (id >= sys->primitives.count) return 0;
        const alea_primitive_entry_t* primitive = &sys->primitives.data[id];
        if (node->primitive.prim_type != primitive->type) return 0;
        if (primitive->type == ALEA_PRIMITIVE_SPHERE ||
            primitive->type == ALEA_PRIMITIVE_SPH) return 1;
        if (primitive->type != ALEA_PRIMITIVE_PLANE ||
            primitive->payload_index >= sys->primitive_planes.count)
            return 0;
        const alea_plane_data_t* p =
            &sys->primitive_planes.data[primitive->payload_index];
        if (!isfinite(p->a) || !isfinite(p->b) ||
            !isfinite(p->c) || !isfinite(p->d)) return 0;
        return (p->a != 0.0 && p->b == 0.0 && p->c == 0.0) ||
               (p->b != 0.0 && p->a == 0.0 && p->c == 0.0) ||
               (p->c != 0.0 && p->a == 0.0 && p->b == 0.0);
    }
    if (op != ALEA_OP_UNION && op != ALEA_OP_INTERSECTION &&
        op != ALEA_OP_DIFFERENCE) return 0;
    return slice_error_bbox_tree_safe(sys, node->operation.left,
                                     depth + 1, work) &&
           slice_error_bbox_tree_safe(sys, node->operation.right,
                                     depth + 1, work);
}

static int slice_error_index_cell_compare(const void* lhs, const void* rhs) {
    const slice_error_index_cell_t* a = lhs;
    const slice_error_index_cell_t* b = rhs;
    if (a->box.min_x < b->box.min_x) return -1;
    if (a->box.min_x > b->box.min_x) return 1;
    return (a->cell_index > b->cell_index) -
           (a->cell_index < b->cell_index);
}

static int slice_error_size_compare(const void* lhs, const void* rhs) {
    const size_t a = *(const size_t*)lhs;
    const size_t b = *(const size_t*)rhs;
    return (a > b) - (a < b);
}

static double slice_error_trusted_lower(double value) {
    return isfinite(value) && value > -1e10 && value < 1e10
        ? value : -INFINITY;
}

static double slice_error_trusted_upper(double value) {
    return isfinite(value) && value > -1e10 && value < 1e10
        ? value : INFINITY;
}

/* Build once per model revision. A cell with uncertain analytic support is
 * always visited; it is never excluded by the sorted bounded-cell index. */
static int slice_error_build_index(alea_slice_error_query_t* query) {
    const alea_system_t* sys = query->sys;
    size_t roots = 0;
    for (size_t i = 0; i < sys->cells.count; ++i)
        roots += sys->cells.data[i].universe_id == 0;
    for (size_t i = 0; i < sys->cells.count; ++i)
        query->has_hierarchy |= sys->cells.data[i].fill_universe > 0 ||
            sys->cells.data[i].lat_type != 0;
    const size_t per_root = sizeof(*query->indexed_cells) +
        sizeof(*query->uncertain_cells);
    size_t capacity = query->options.max_index_bytes / per_root;
    if (capacity > roots) capacity = roots;
    query->indexed_cells = capacity
        ? malloc(capacity * sizeof(*query->indexed_cells)) : NULL;
    query->uncertain_cells = capacity
        ? malloc(capacity * sizeof(*query->uncertain_cells)) : NULL;
    if (capacity && (!query->indexed_cells || !query->uncertain_cells))
        return -1;
    query->index_bytes = capacity * per_root;
    query->uncached_cell_begin = sys->cells.count;
    for (size_t i = 0; i < sys->cells.count; ++i) {
        const alea_cell_entry_t* cell = &sys->cells.data[i];
        if (cell->universe_id != 0) continue;
        if (query->indexed_count + query->uncertain_count == capacity) {
            query->uncached_cell_begin = i;
            break;
        }
        size_t work = query->options.scan_options.max_active_boundary_tests;
        if (work > 4096) work = 4096;
        if (cell->original_root_node_id != ALEA_NODE_ID_INVALID ||
            !slice_error_bbox_tree_safe(sys, cell->root_node_id, 0,
                                        &work)) {
            query->uncertain_cells[query->uncertain_count++] = i;
            continue;
        }
        alea_bbox_t box = alea_get_bbox(sys, cell->root_node_id);
        if (!(box.min_x <= box.max_x && box.min_y <= box.max_y &&
              box.min_z <= box.max_z)) {
            query->uncertain_cells[query->uncertain_count++] = i;
            continue;
        }
        box.min_x = slice_error_trusted_lower(box.min_x);
        box.max_x = slice_error_trusted_upper(box.max_x);
        box.min_y = slice_error_trusted_lower(box.min_y);
        box.max_y = slice_error_trusted_upper(box.max_y);
        box.min_z = slice_error_trusted_lower(box.min_z);
        box.max_z = slice_error_trusted_upper(box.max_z);
        query->indexed_cells[query->indexed_count++] =
            (slice_error_index_cell_t){i, box};
    }
    if (query->indexed_count > 1)
        qsort(query->indexed_cells, query->indexed_count,
              sizeof(*query->indexed_cells), slice_error_index_cell_compare);
    return 0;
}

/* Finite analytic bounds exclude a root cell only when separated from the
 * whole tile slab by more than coordinate rounding. */
static int slice_error_box_outside_tile(
    const alea_bbox_t* box, const alea_slice_error_query_options_t* o,
    const alea_transition_slice_critical_tile_t* tile) {
    double lower[3], upper[3];
    for (int world_axis = 0; world_axis < 3; ++world_axis) {
        const int slice_axis = slice_error_slice_axis_for_world(
            &o->view, world_axis);
        if (slice_axis < 0) {
            lower[world_axis] = upper[world_axis] =
                o->view.plane.origin[world_axis];
            continue;
        }
        const double sign = slice_error_slice_axis_sign(
            &o->view, world_axis);
        lower[world_axis] = o->view.plane.origin[world_axis] +
            sign * (sign > 0.0 ? tile->uv_min[slice_axis]
                               : tile->uv_max[slice_axis]);
        upper[world_axis] = o->view.plane.origin[world_axis] +
            sign * (sign > 0.0 ? tile->uv_max[slice_axis]
                               : tile->uv_min[slice_axis]);
    }
    const double lo[3] = {box->min_x, box->min_y, box->min_z};
    const double hi[3] = {box->max_x, box->max_y, box->max_z};
    for (int axis = 0; axis < 3; ++axis) {
        if (!isfinite(lower[axis]) || !isfinite(upper[axis])) return 0;
        const double margin = 64.0 * DBL_EPSILON *
            fmax(1.0, fmax(fabs(lower[axis]), fabs(upper[axis])));
        if (hi[axis] < lower[axis] - margin ||
            lo[axis] > upper[axis] + margin) return 1;
    }
    return 0;
}

static int slice_error_plane_index(const slice_error_axis_plane_t* planes,
                                   size_t count, uint32_t primitive_id) {
    for (size_t i = 0; i < count; ++i)
        if (planes[i].primitive_id == primitive_id) return (int)i;
    return -1;
}

static int slice_error_collect_primitives(const alea_system_t* sys,
                                          alea_node_id_t id, size_t depth,
                                          slice_error_axis_plane_t* planes,
                                          size_t* count, size_t max_count,
                                          size_t* work) {
    if (!*work || depth > 128 || id >= sys->nodes.count) return 0;
    --*work;
    const alea_node_t* node = &sys->nodes.data[id];
    const alea_operation_t op = ALEA_GET_OPERATION(node);
    if (op == ALEA_OP_PRIMITIVE) {
        const uint32_t primitive_id = node->primitive.primitive_id;
        if (primitive_id >= sys->primitives.count) return 0;
        if (slice_error_plane_index(planes, *count, primitive_id) >= 0)
            return 1;
        if (*count == max_count) return 0;
        planes[*count].primitive_id = primitive_id;
        ++*count;
        return 1;
    }
    if (op != ALEA_OP_UNION && op != ALEA_OP_INTERSECTION &&
        op != ALEA_OP_DIFFERENCE && op != ALEA_OP_COMPLEMENT) return 0;
    if (!slice_error_collect_primitives(sys, node->operation.left, depth + 1,
                                        planes, count, max_count, work)) return 0;
    return op == ALEA_OP_COMPLEMENT ||
        slice_error_collect_primitives(sys, node->operation.right, depth + 1,
                                       planes, count, max_count, work);
}

typedef struct {
    double coordinate;
    double uncertainty;
    int surface_id;
    uint32_t primitive_id;
} slice_error_grid_line_t;

typedef struct {
    alea_point_coverage_kind_t kind;
    size_t owner_count;
    int owner_cell_ids[ALEA_SLICE_ERROR_OWNER_CAPACITY];
} slice_error_face_t;

static int slice_error_grid_line_compare(const void* aa, const void* bb) {
    const slice_error_grid_line_t* a = aa;
    const slice_error_grid_line_t* b = bb;
    return (a->coordinate > b->coordinate) -
           (a->coordinate < b->coordinate);
}

static int slice_error_axis_view_supported(const alea_slice_view_t* view) {
    const alea_slice_plane_t* p = &view->plane;
    int axes[3] = {-1, -1, -1};
    const double* vectors[3] = {p->u_axis, p->v_axis, p->normal};
    for (int vector = 0; vector < 3; ++vector) {
        for (int world_axis = 0; world_axis < 3; ++world_axis) {
            const double component = vectors[vector][world_axis];
            if (component == 0.0) continue;
            if (fabs(component) != 1.0 || axes[vector] >= 0)
                return 0;
            axes[vector] = world_axis;
        }
        if (axes[vector] < 0) return 0;
    }
    return axes[0] != axes[1] && axes[0] != axes[2] &&
           axes[1] != axes[2];
}

static int slice_error_slice_axis_for_world(const alea_slice_view_t* view,
                                           int world_axis) {
    return view->plane.u_axis[world_axis] != 0.0 ? 0
        : view->plane.v_axis[world_axis] != 0.0 ? 1 : -1;
}

static double slice_error_slice_axis_sign(const alea_slice_view_t* view,
                                          int world_axis) {
    return view->plane.u_axis[world_axis] != 0.0
        ? view->plane.u_axis[world_axis]
        : view->plane.v_axis[world_axis];
}

/* Return the sign of a sphere predicate only when its whole tile is
 * separated from the sphere. Expanded coordinate ranges cover rounding in
 * the world-to-slice map; close and tangent cases remain unresolved. */
static int slice_error_sphere_tile_sign(
    const alea_sphere_data_t* sphere, const alea_slice_view_t* view,
    const alea_transition_slice_critical_tile_t* tile, int* sign) {
    const double center[3] = {sphere->center_x, sphere->center_y,
                              sphere->center_z};
    if (!isfinite(sphere->radius) || !(sphere->radius > 0.0)) return 0;
    long double min_squared = 0.0L, max_squared = 0.0L;
    for (int axis = 0; axis < 3; ++axis) {
        if (!isfinite(center[axis])) return 0;
        const int slice_axis = slice_error_slice_axis_for_world(view, axis);
        const long double origin = view->plane.origin[axis];
        long double lo = origin, hi = origin;
        long double extent = 0.0L;
        if (slice_axis >= 0) {
            const long double direction =
                slice_error_slice_axis_sign(view, axis);
            const long double a = tile->uv_min[slice_axis];
            const long double b = tile->uv_max[slice_axis];
            lo = origin + direction * (direction > 0.0L ? a : b);
            hi = origin + direction * (direction > 0.0L ? b : a);
            extent = fmaxl(fabsl(a), fabsl(b));
        }
        const long double scale = fmaxl(1.0L, fabsl(origin) +
            extent + fabsl((long double)center[axis]) +
            (long double)sphere->radius);
        /* Keep the proof inside the finite range of the double-precision
         * primitive evaluator, including its sum of three squared terms. */
        if (!(scale < 1.0e150L)) return 0;
        const long double margin = 128.0L * DBL_EPSILON * scale;
        lo -= margin;
        hi += margin;
        if (!isfinite(lo) || !isfinite(hi) || lo > hi) return 0;
        const long double c = center[axis];
        const long double nearest = c < lo ? lo - c : c > hi ? c - hi : 0.0L;
        const long double farthest = fmaxl(fabsl(lo - c), fabsl(hi - c));
        min_squared += nearest * nearest;
        max_squared += farthest * farthest;
    }
    const long double radius_squared =
        (long double)sphere->radius * (long double)sphere->radius;
    if (!isfinite(min_squared) || !isfinite(max_squared) ||
        !isfinite(radius_squared)) return 0;
    const long double tolerance = 128.0L * DBL_EPSILON *
        fmaxl(1.0L, fmaxl(max_squared, radius_squared));
    if (!isfinite(tolerance)) return 0;
    if (min_squared > radius_squared + tolerance) {
        *sign = 1;
        return 1;
    }
    if (max_squared + tolerance < radius_squared) {
        *sign = -1;
        return 1;
    }
    return 0;
}

/* A one-cell sphere section forms exactly two faces, even when the circle
 * crosses the page edge. No polygon approximation is used in the ownership
 * proof. Tangencies and corner coincidences remain unresolved. */
static int slice_error_classify_isolated_circle(
    const alea_system_t* sys, const alea_slice_view_t* view,
    const alea_transition_slice_critical_tile_t* tile, size_t cell_index,
    size_t contextual_bytes, size_t output_limit,
    alea_slice_error_page_t* out) {
    const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
    if (cell->fill_universe > 0 || cell->lat_type != 0 ||
        cell->original_root_node_id != ALEA_NODE_ID_INVALID ||
        cell->root_node_id >= sys->nodes.count) return 0;
    const alea_node_t* node = &sys->nodes.data[cell->root_node_id];
    if (ALEA_GET_OPERATION(node) != ALEA_OP_PRIMITIVE ||
        node->primitive.primitive_id >= sys->primitives.count) return 0;
    const uint32_t primitive_id = node->primitive.primitive_id;
    const alea_primitive_entry_t* primitive =
        &sys->primitives.data[primitive_id];
    alea_sphere_data_t sphere;
    if (primitive->type == ALEA_PRIMITIVE_SPHERE &&
        primitive->payload_index < sys->primitive_spheres.count) {
        sphere = sys->primitive_spheres.data[primitive->payload_index];
    } else if (primitive->type == ALEA_PRIMITIVE_SPH &&
               primitive->payload_index < sys->primitive_sphs.count) {
        const alea_sph_data_t* sph =
            &sys->primitive_sphs.data[primitive->payload_index];
        sphere = (alea_sphere_data_t){sph->center_x, sph->center_y,
                                      sph->center_z, sph->radius};
    } else return 0;
    const double center[3] = {sphere.center_x, sphere.center_y,
                              sphere.center_z};
    if (!(sphere.radius > 0.0) || !isfinite(sphere.radius)) return 0;
    int surface_id = 0;
    for (size_t i = 0; i < sys->surfaces.count; ++i) {
        const alea_surface_entry_t* surface = &sys->surfaces.data[i];
        if (surface->primitive_id != primitive_id) continue;
        if (surface_id || surface->mc_surface_id <= 0) return 0;
        surface_id = surface->mc_surface_id;
    }
    if (!surface_id) return 0;
    long double center_uv[2] = {0.0L, 0.0L};
    long double normal_distance = 0.0L;
    long double scale = fmaxl(1.0L, sphere.radius);
    for (int axis = 0; axis < 3; ++axis) {
        if (!isfinite(center[axis])) return 0;
        const long double origin = view->plane.origin[axis];
        const long double coordinate = center[axis];
        scale = fmaxl(scale, fabsl(origin) + fabsl(coordinate));
        const int slice_axis = slice_error_slice_axis_for_world(view, axis);
        if (slice_axis < 0)
            normal_distance = coordinate - origin;
        else
            center_uv[slice_axis] = (coordinate - origin) /
                slice_error_slice_axis_sign(view, axis);
    }
    for (int axis = 0; axis < 2; ++axis)
        scale = fmaxl(scale, fmaxl(fabsl((long double)tile->uv_min[axis]),
                                    fabsl((long double)tile->uv_max[axis])));
    if (!(scale < 1.0e150L)) return 0;
    const long double r2 = (long double)sphere.radius * sphere.radius;
    const long double h2 = normal_distance * normal_distance;
    const long double difference = r2 - h2;
    const long double squared_error = 256.0L * DBL_EPSILON *
        fmaxl(1.0L, r2 + h2);
    if (!(difference > squared_error) || !isfinite(difference)) return 0;
    const long double radius = sqrtl(difference);
    const long double uncertainty = 1024.0L * DBL_EPSILON *
        scale * (1.0L + (long double)sphere.radius / radius);
    if (!isfinite(radius) || !isfinite(uncertainty) ||
        !(uncertainty < radius / 16.0L)) return 0;
    long double closest_squared = 0.0L, farthest_squared = 0.0L;
    const long double boundary_margin = 8.0L * uncertainty;
    for (int axis = 0; axis < 2; ++axis) {
        const long double lo = tile->uv_min[axis];
        const long double hi = tile->uv_max[axis];
        const long double c = center_uv[axis];
        if (!isfinite(c)) return 0;
        const long double nearest = c < lo ? lo - c : c > hi ? c - hi : 0.0L;
        const long double farthest = fmaxl(fabsl(lo - c), fabsl(hi - c));
        closest_squared += nearest * nearest;
        farthest_squared += farthest * farthest;
        /* Edge tangencies are numerically ill-conditioned for clipping. */
        if (fabsl(fabsl(lo - c) - radius) <= boundary_margin ||
            fabsl(fabsl(hi - c) - radius) <= boundary_margin) return 0;
    }
    const long double circle_r2 = radius * radius;
    const long double boundary_squared_margin =
        4.0L * radius * boundary_margin +
        boundary_margin * boundary_margin;
    if (!(closest_squared + boundary_squared_margin < circle_r2) ||
        !(farthest_squared > circle_r2 + boundary_squared_margin)) return 0;
    for (int x = 0; x < 2; ++x)
        for (int y = 0; y < 2; ++y) {
            const long double dx = (long double)(x ? tile->uv_max[0]
                                                   : tile->uv_min[0]) -
                center_uv[0];
            const long double dy = (long double)(y ? tile->uv_max[1]
                                                   : tile->uv_min[1]) -
                center_uv[1];
            if (!(fabsl(dx * dx + dy * dy - circle_r2) >
                  boundary_squared_margin)) return 0;
        }
    if (contextual_bytes > output_limit ||
        sizeof(alea_slice_error_circle_t) > output_limit - contextual_bytes)
        return 2;
    const int flip = (node->primitive.sense > 0) !=
                     (node->primitive.inverted != 0);
    const int inside_owned = !flip;
    const double u = (double)center_uv[0], v = (double)center_uv[1];
    double center_world[3], corner_world[3];
    for (int axis = 0; axis < 3; ++axis) {
        center_world[axis] = view->plane.origin[axis] +
            view->plane.u_axis[axis] * u + view->plane.v_axis[axis] * v;
        corner_world[axis] = view->plane.origin[axis] +
            view->plane.u_axis[axis] * tile->uv_min[0] +
            view->plane.v_axis[axis] * tile->uv_min[1];
    }
    const long double corner_dx = (long double)tile->uv_min[0] - center_uv[0];
    const long double corner_dy = (long double)tile->uv_min[1] - center_uv[1];
    const int corner_sphere_inside =
        corner_dx * corner_dx + corner_dy * corner_dy < circle_r2;
    const int corner_owned = corner_sphere_inside == inside_owned;
    if (alea_point_inside(sys, cell->root_node_id,
            center_world[0], center_world[1], center_world[2]) !=
            inside_owned ||
        alea_point_inside(sys, cell->root_node_id,
            corner_world[0], corner_world[1], corner_world[2]) !=
            corner_owned) return 0;
    out->circles = calloc(1, sizeof(*out->circles));
    if (!out->circles) return -1;
    out->circle_count = 1;
    alea_slice_error_circle_t* circle = out->circles;
    circle->surface_id = surface_id;
    circle->primitive_id = primitive_id;
    circle->center_uv[0] = u;
    circle->center_uv[1] = v;
    circle->radius = (double)radius;
    circle->geometry_uncertainty = (double)uncertainty;
    circle->start_angle = 0.0;
    circle->end_angle = 6.283185307179586476925286766559;
    circle->inside_kind = inside_owned ? ALEA_POINT_COVERAGE_UNIQUE
                                       : ALEA_POINT_COVERAGE_GAP;
    circle->outside_kind = inside_owned ? ALEA_POINT_COVERAGE_GAP
                                        : ALEA_POINT_COVERAGE_UNIQUE;
    if (inside_owned) {
        circle->inside_owner_count = 1;
        circle->inside_owner_cell_ids[0] = cell->mc_cell_id;
    } else {
        circle->outside_owner_count = 1;
        circle->outside_owner_cell_ids[0] = cell->mc_cell_id;
    }
    return 1;
}

/* Two inside-sphere cells with complete circles admit bounded proofs for
 * strictly nested or transverse pairs. Tangent, coincident, and clipped pairs
 * remain unresolved. */
static int slice_error_classify_crossing_circles(
    const alea_system_t* sys, const alea_slice_view_t* view,
    const alea_transition_slice_critical_tile_t* tile,
    const size_t cells[2], size_t scratch_available, size_t* scratch_used,
    size_t contextual_bytes, size_t output_limit,
    alea_slice_error_page_t* out) {
    const long double tau = 6.283185307179586476925286766559005768L;
    alea_slice_error_page_t individual[2];
    memset(individual, 0, sizeof(individual));
    *scratch_used = 0;
    for (int i = 0; i < 2; ++i) {
        const alea_cell_entry_t* cell = &sys->cells.data[cells[i]];
        if (cell->fill_universe > 0 || cell->lat_type != 0 ||
            cell->original_root_node_id != ALEA_NODE_ID_INVALID ||
            cell->root_node_id >= sys->nodes.count) return 0;
        const alea_node_t* node = &sys->nodes.data[cell->root_node_id];
        if (ALEA_GET_OPERATION(node) != ALEA_OP_PRIMITIVE ||
            node->primitive.primitive_id >= sys->primitives.count) return 0;
        const alea_primitive_entry_t* primitive =
            &sys->primitives.data[node->primitive.primitive_id];
        if (primitive->type != ALEA_PRIMITIVE_SPHERE &&
            primitive->type != ALEA_PRIMITIVE_SPH) return 0;
    }
    if (scratch_available < 2 * sizeof(alea_slice_error_circle_t)) return 3;
    int status = 0;
    for (int i = 0; i < 2; ++i) {
        const int result = slice_error_classify_isolated_circle(
            sys, view, tile, cells[i], 0, SIZE_MAX, &individual[i]);
        if (result < 0) { status = -1; goto done; }
        if (result != 1) goto done;
        *scratch_used = (size_t)(i + 1) *
            sizeof(alea_slice_error_circle_t);
        const alea_slice_error_circle_t* circle = individual[i].circles;
        if (circle->inside_kind != ALEA_POINT_COVERAGE_UNIQUE ||
            circle->outside_kind != ALEA_POINT_COVERAGE_GAP ||
            circle->inside_owner_count != 1 ||
            circle->outside_owner_count != 0) goto done;
        for (int axis = 0; axis < 2; ++axis) {
            const long double center = circle->center_uv[axis];
            const long double radius = circle->radius;
            const long double error = 8.0L * circle->geometry_uncertainty;
            if (!(center - radius - error > tile->uv_min[axis]) ||
                !(center + radius + error < tile->uv_max[axis])) goto done;
        }
    }
    const alea_slice_error_circle_t* a = individual[0].circles;
    const alea_slice_error_circle_t* b = individual[1].circles;
    if (a->primitive_id == b->primitive_id ||
        a->surface_id == b->surface_id ||
        a->inside_owner_cell_ids[0] == b->inside_owner_cell_ids[0]) goto done;
    const long double dx = (long double)b->center_uv[0] - a->center_uv[0];
    const long double dy = (long double)b->center_uv[1] - a->center_uv[1];
    const long double distance = sqrtl(dx * dx + dy * dy);
    const long double r[2] = {a->radius, b->radius};
    const long double scale = fmaxl(1.0L,
        fmaxl(distance, fmaxl(r[0], r[1])));
    const long double margin = 64.0L *
        ((long double)a->geometry_uncertainty +
         (long double)b->geometry_uncertainty) +
        512.0L * DBL_EPSILON * scale;
    if (!isfinite(distance) || !isfinite(margin)) goto done;
    const int outer = r[0] >= r[1] ? 0 : 1;
    const int inner = 1 - outer;
    const long double nesting_clearance = r[outer] - r[inner] - distance;
    if (nesting_clearance > 64.0L * margin) {
        const long double probe = 16.0L * margin;
        if (!(probe < r[inner] / 8.0L) ||
            !(probe < nesting_clearance / 4.0L)) goto done;
        for (int boundary = 0; boundary < 2; ++boundary) {
            const alea_slice_error_circle_t* circle =
                individual[boundary].circles;
            for (int side = 0; side < 2; ++side) {
                const double u = circle->center_uv[0] +
                    circle->radius + (side ? (double)probe : -(double)probe);
                const double v = circle->center_uv[1];
                double world[3];
                for (int axis = 0; axis < 3; ++axis)
                    world[axis] = view->plane.origin[axis] +
                        view->plane.u_axis[axis] * u +
                        view->plane.v_axis[axis] * v;
                const int outer_expected = boundary == outer ? !side : 1;
                const int inner_expected = boundary == inner && !side;
                if (alea_point_inside(sys,
                        sys->cells.data[cells[outer]].root_node_id,
                        world[0], world[1], world[2]) != outer_expected ||
                    alea_point_inside(sys,
                        sys->cells.data[cells[inner]].root_node_id,
                        world[0], world[1], world[2]) != inner_expected)
                    goto done;
            }
        }
        if (contextual_bytes > output_limit ||
            2 * sizeof(alea_slice_error_circle_t) >
                output_limit - contextual_bytes) {
            status = 2;
            goto done;
        }
        out->circles = calloc(2, sizeof(*out->circles));
        if (!out->circles) { status = -1; goto done; }
        out->circles[0] = *individual[outer].circles;
        out->circles[1] = *individual[inner].circles;
        alea_slice_error_circle_t* inner_circle = &out->circles[1];
        inner_circle->inside_kind = ALEA_POINT_COVERAGE_OVERLAP;
        inner_circle->inside_owner_count = 2;
        inner_circle->inside_owner_cell_ids[0] = a->inside_owner_cell_ids[0];
        inner_circle->inside_owner_cell_ids[1] = b->inside_owner_cell_ids[0];
        inner_circle->outside_kind = ALEA_POINT_COVERAGE_UNIQUE;
        inner_circle->outside_owner_count = 1;
        inner_circle->outside_owner_cell_ids[0] =
            individual[outer].circles->inside_owner_cell_ids[0];
        out->circle_count = 2;
        status = 1;
        goto done;
    }
    if (!(distance > fabsl(r[0] - r[1]) + margin) ||
        !(distance + margin < r[0] + r[1])) goto done;
    const long double along =
        (r[0] * r[0] - r[1] * r[1] + distance * distance) /
        (2.0L * distance);
    const long double height_squared = r[0] * r[0] - along * along;
    if (!(height_squared > 16.0L * margin * scale) ||
        !isfinite(height_squared)) goto done;
    const long double height = sqrtl(height_squared);
    const long double intersection_uncertainty =
        64.0L * margin * (1.0L + scale / height);
    if (!isfinite(intersection_uncertainty) ||
        !(intersection_uncertainty < height / 16.0L)) goto done;
    const long double unit_x = dx / distance, unit_y = dy / distance;
    const long double base_x = (long double)a->center_uv[0] + along * unit_x;
    const long double base_y = (long double)a->center_uv[1] + along * unit_y;
    long double crossing[2][2] = {
        {base_x - height * unit_y, base_y + height * unit_x},
        {base_x + height * unit_y, base_y - height * unit_x}
    };
    double angle[2][2];
    for (int point = 0; point < 2; ++point) {
        for (int axis = 0; axis < 2; ++axis)
            if (!(crossing[point][axis] >
                      (long double)tile->uv_min[axis] +
                          intersection_uncertainty) ||
                !(crossing[point][axis] <
                      (long double)tile->uv_max[axis] -
                          intersection_uncertainty)) goto done;
        for (int i = 0; i < 2; ++i) {
            const alea_slice_error_circle_t* circle =
                individual[i].circles;
            const long double px = crossing[point][0] - circle->center_uv[0];
            const long double py = crossing[point][1] - circle->center_uv[1];
            if (!(fabsl(px * px + py * py - r[i] * r[i]) <
                  8.0L * margin * scale)) goto done;
            long double theta = atan2l(py, px);
            if (theta < 0.0L) theta += tau;
            angle[i][point] = (double)theta;
        }
    }
    alea_slice_error_circle_t arcs[4];
    for (int i = 0; i < 2; ++i) {
        const int other = 1 - i;
        const alea_slice_error_circle_t* own = individual[i].circles;
        const alea_slice_error_circle_t* neighbor =
            individual[other].circles;
        double first = fmin(angle[i][0], angle[i][1]);
        double second = fmax(angle[i][0], angle[i][1]);
        const long double angular_margin = 8.0L * margin / r[i];
        if (!((long double)second - first > angular_margin) ||
            !(tau - ((long double)second - first) > angular_margin)) goto done;
        for (int part = 0; part < 2; ++part) {
            const long double start = part == 0 ? first : second;
            const long double end = part == 0 ? second : first + tau;
            const long double midpoint = (start + end) * 0.5L;
            const long double ux = cosl(midpoint), uy = sinl(midpoint);
            const long double px = (long double)own->center_uv[0] + r[i] * ux;
            const long double py = (long double)own->center_uv[1] + r[i] * uy;
            const long double ox = px - neighbor->center_uv[0];
            const long double oy = py - neighbor->center_uv[1];
            const long double other_distance = sqrtl(ox * ox + oy * oy);
            const int other_inside = other_distance < r[other];
            const long double separation = fabsl(other_distance - r[other]);
            const long double probe = 16.0L * margin;
            if (!(separation > 4.0L * probe) ||
                !(probe < r[i] / 8.0L) || !isfinite(separation)) goto done;
            for (int side = 0; side < 2; ++side) {
                const long double radial = r[i] + (side ? probe : -probe);
                const long double u = (long double)own->center_uv[0] +
                    radial * ux;
                const long double v = (long double)own->center_uv[1] +
                    radial * uy;
                double world[3];
                for (int axis = 0; axis < 3; ++axis)
                    world[axis] = view->plane.origin[axis] +
                        view->plane.u_axis[axis] * (double)u +
                        view->plane.v_axis[axis] * (double)v;
                if (alea_point_inside(sys, sys->cells.data[cells[i]].root_node_id,
                        world[0], world[1], world[2]) != !side ||
                    alea_point_inside(sys,
                        sys->cells.data[cells[other]].root_node_id,
                        world[0], world[1], world[2]) != other_inside)
                    goto done;
            }
            alea_slice_error_circle_t* arc = &arcs[2 * i + part];
            *arc = *own;
            arc->start_angle = (double)start;
            arc->end_angle = (double)end;
            arc->geometry_uncertainty = (double)fmaxl(
                (long double)own->geometry_uncertainty,
                intersection_uncertainty);
            memset(arc->inside_owner_cell_ids, 0,
                   sizeof(arc->inside_owner_cell_ids));
            memset(arc->outside_owner_cell_ids, 0,
                   sizeof(arc->outside_owner_cell_ids));
            if (other_inside) {
                arc->inside_kind = ALEA_POINT_COVERAGE_OVERLAP;
                arc->inside_owner_count = 2;
                arc->inside_owner_cell_ids[0] = a->inside_owner_cell_ids[0];
                arc->inside_owner_cell_ids[1] = b->inside_owner_cell_ids[0];
                arc->outside_kind = ALEA_POINT_COVERAGE_UNIQUE;
                arc->outside_owner_count = 1;
                arc->outside_owner_cell_ids[0] =
                    neighbor->inside_owner_cell_ids[0];
            } else {
                arc->inside_kind = ALEA_POINT_COVERAGE_UNIQUE;
                arc->inside_owner_count = 1;
                arc->inside_owner_cell_ids[0] =
                    own->inside_owner_cell_ids[0];
                arc->outside_kind = ALEA_POINT_COVERAGE_GAP;
                arc->outside_owner_count = 0;
            }
        }
    }
    if (contextual_bytes > output_limit ||
        4 * sizeof(alea_slice_error_circle_t) >
            output_limit - contextual_bytes) {
        status = 2;
        goto done;
    }
    out->circles = calloc(4, sizeof(*out->circles));
    if (!out->circles) { status = -1; goto done; }
    memcpy(out->circles, arcs, sizeof(arcs));
    out->circle_count = 4;
    status = 1;
done:
    free(individual[0].circles);
    free(individual[1].circles);
    return status;
}

/* This proof path accepts root cells formed from planes aligned to world
 * axes. A plane parallel to the slice has a constant sign when separated
 * from it. Each primitive sign is fixed in an open grid face. */
static int slice_error_axis_node_inside(
    const alea_system_t* sys, alea_node_id_t id,
    const slice_error_axis_plane_t* planes, size_t plane_count,
    const slice_error_grid_line_t* x, size_t xi,
    const slice_error_grid_line_t* y, size_t yi,
    size_t depth, size_t* work_remaining, int* out_inside) {
    if (*work_remaining == 0) return 0;
    --*work_remaining;
    if (id >= sys->nodes.count || depth > 128) return 0;
    const alea_node_t* node = &sys->nodes.data[id];
    const alea_operation_t op = ALEA_GET_OPERATION(node);
    if (op == ALEA_OP_PRIMITIVE) {
        const uint32_t primitive_id = node->primitive.primitive_id;
        if (primitive_id >= sys->primitives.count) return 0;
        const int plane_index = slice_error_plane_index(
            planes, plane_count, primitive_id);
        if (plane_index < 0) return 0;
        const slice_error_axis_plane_t* plane = &planes[plane_index];
        int raw_negative;
        if (plane->axis == 2) {
            raw_negative = plane->coefficient_sign < 0;
        } else {
            const slice_error_grid_line_t* lines = plane->axis == 0 ? x : y;
            const size_t index = plane->axis == 0 ? xi : yi;
            int side;
            if (lines[index + 1].coordinate <= plane->coordinate)
                side = -1;
            else if (lines[index].coordinate >= plane->coordinate)
                side = 1;
            else
                return 0;
            raw_negative = side * plane->coefficient_sign < 0;
        }
        const int flip = (node->primitive.sense > 0) !=
                         (node->primitive.inverted != 0);
        *out_inside = raw_negative != flip;
        return 1;
    }
    if (op != ALEA_OP_UNION && op != ALEA_OP_INTERSECTION &&
        op != ALEA_OP_DIFFERENCE && op != ALEA_OP_COMPLEMENT) return 0;
    int left = 0, right = 0;
    if (!slice_error_axis_node_inside(sys, node->operation.left, planes,
                                      plane_count,
                                      x, xi, y, yi, depth + 1,
                                      work_remaining, &left)) return 0;
    if (op == ALEA_OP_COMPLEMENT) {
        *out_inside = !left;
        return 1;
    }
    if (!slice_error_axis_node_inside(sys, node->operation.right, planes,
                                      plane_count,
                                      x, xi, y, yi, depth + 1,
                                      work_remaining, &right)) return 0;
    *out_inside = op == ALEA_OP_UNION ? left || right
        : op == ALEA_OP_INTERSECTION ? left && right
        : left && !right;
    return 1;
}

/* Bound a page in an occurrence's local frame. The extra rounding envelope
 * covers both the view-to-world and inverse-fill affine evaluations. */
static int slice_error_occurrence_tile_box(
    const alea_slice_view_t* view,
    const alea_transition_slice_critical_tile_t* tile,
    const alea_matrix_t* transform, alea_bbox_t* box) {
    long double low[3] = {LDBL_MAX, LDBL_MAX, LDBL_MAX};
    long double high[3] = {-LDBL_MAX, -LDBL_MAX, -LDBL_MAX};
    long double scale[3] = {0, 0, 0};
    if (transform && !transform->has_inverse) return 0;
    for (int iu = 0; iu < 2; ++iu)
        for (int iv = 0; iv < 2; ++iv) {
            const long double u = tile->uv_min[0] * (1 - iu) +
                                  tile->uv_max[0] * iu;
            const long double v = tile->uv_min[1] * (1 - iv) +
                                  tile->uv_max[1] * iv;
            long double world[3], world_scale[3];
            for (int axis = 0; axis < 3; ++axis)
                world[axis] = (long double)view->plane.origin[axis] +
                    (long double)view->plane.u_axis[axis] * u +
                    (long double)view->plane.v_axis[axis] * v;
            for (int axis = 0; axis < 3; ++axis)
                world_scale[axis] =
                    fabsl((long double)view->plane.origin[axis]) +
                    fabsl((long double)view->plane.u_axis[axis] * u) +
                    fabsl((long double)view->plane.v_axis[axis] * v);
            for (int axis = 0; axis < 3; ++axis) {
                long double value = world[axis];
                long double magnitude = world_scale[axis];
                if (transform) {
                    const double* row = &transform->inv[4 * axis];
                    value = (long double)row[0] * world[0] +
                        (long double)row[1] * world[1] +
                        (long double)row[2] * world[2] + row[3];
                    magnitude = fabsl((long double)row[0] * world[0]) +
                        fabsl((long double)row[1] * world[1]) +
                        fabsl((long double)row[2] * world[2]) +
                        fabsl((long double)row[3]) +
                        fabsl((long double)row[0]) * world_scale[0] *
                            DBL_EPSILON +
                        fabsl((long double)row[1]) * world_scale[1] *
                            DBL_EPSILON +
                        fabsl((long double)row[2]) * world_scale[2] *
                            DBL_EPSILON;
                }
                if (!isfinite(value) || fabsl(value) > DBL_MAX / 2 ||
                    !isfinite(magnitude)) return 0;
                if (value < low[axis]) low[axis] = value;
                if (value > high[axis]) high[axis] = value;
                if (magnitude > scale[axis]) scale[axis] = magnitude;
            }
        }
    double* bounds[6] = {&box->min_x, &box->max_x,
                         &box->min_y, &box->max_y,
                         &box->min_z, &box->max_z};
    for (int axis = 0; axis < 3; ++axis) {
        const long double margin = 1024.0L * DBL_EPSILON *
            fmaxl(1.0L, scale[axis]);
        if (!isfinite(margin) || margin > DBL_MAX / 4) return 0;
        *bounds[2 * axis] = nextafter((double)(low[axis] - margin),
                                     -INFINITY);
        *bounds[2 * axis + 1] = nextafter((double)(high[axis] + margin),
                                         INFINITY);
    }
    return 1;
}

/* -1: outside, +1: inside, 0: not constant, -2: malformed or work limit. */
static int slice_error_occurrence_node_constant(
    const alea_system_t* sys, alea_node_id_t id, const alea_bbox_t* box,
    size_t depth, size_t* work) {
    if (!*work || depth > 128 || id >= sys->nodes.count) return -2;
    --*work;
    const alea_node_t* node = &sys->nodes.data[id];
    const alea_operation_t op = ALEA_GET_OPERATION(node);
    if (op == ALEA_OP_PRIMITIVE) {
        const uint32_t id = node->primitive.primitive_id;
        if (id >= sys->primitives.count) return -2;
        const alea_primitive_entry_t* primitive = &sys->primitives.data[id];
        const alea_primitive_type_t type = primitive->type;
        if (node->primitive.prim_type != type ||
            (type != ALEA_PRIMITIVE_PLANE &&
             type != ALEA_PRIMITIVE_SPHERE &&
             type != ALEA_PRIMITIVE_SPH &&
             type != ALEA_PRIMITIVE_CYLINDER_X &&
             type != ALEA_PRIMITIVE_CYLINDER_Y &&
             type != ALEA_PRIMITIVE_CYLINDER_Z &&
             type != ALEA_PRIMITIVE_RPP)) return 0;
        alea_primitive_data_t data;
        if (!alea_primitive_copy_data(sys, id, &data)) return -2;
        const alea_interval_t interval =
            alea_primitive_interval_eval(type, &data, box);
        long double scale = 1.0L;
        if (type == ALEA_PRIMITIVE_PLANE) {
            const alea_plane_data_t* p = &data.plane;
            scale += fabsl(p->a) * fmaxl(fabsl(box->min_x), fabsl(box->max_x)) +
                fabsl(p->b) * fmaxl(fabsl(box->min_y), fabsl(box->max_y)) +
                fabsl(p->c) * fmaxl(fabsl(box->min_z), fabsl(box->max_z)) +
                fabsl(p->d);
        } else if (type == ALEA_PRIMITIVE_RPP) {
            const alea_box_data_t* b = &data.box;
            if (fabsl(b->min_x) + fabsl(b->max_x) > DBL_MAX / 2 ||
                fabsl(b->min_y) + fabsl(b->max_y) > DBL_MAX / 2 ||
                fabsl(b->min_z) + fabsl(b->max_z) > DBL_MAX / 2)
                return 0;
            scale += fabsl(b->min_x) + fabsl(b->max_x) +
                fabsl(b->min_y) + fabsl(b->max_y) +
                fabsl(b->min_z) + fabsl(b->max_z) +
                fabsl(box->min_x) + fabsl(box->max_x) +
                fabsl(box->min_y) + fabsl(box->max_y) +
                fabsl(box->min_z) + fabsl(box->max_z);
        } else {
            double cx = 0.0, cy = 0.0, cz = 0.0, radius = 0.0;
            if (type == ALEA_PRIMITIVE_SPHERE) {
                cx = data.sphere.center_x;
                cy = data.sphere.center_y;
                cz = data.sphere.center_z;
                radius = data.sphere.radius;
            } else if (type == ALEA_PRIMITIVE_SPH) {
                cx = data.sph.center_x;
                cy = data.sph.center_y;
                cz = data.sph.center_z;
                radius = data.sph.radius;
            } else if (type == ALEA_PRIMITIVE_CYLINDER_X) {
                cy = data.cyl_x.center_y;
                cz = data.cyl_x.center_z;
                radius = data.cyl_x.radius;
            } else if (type == ALEA_PRIMITIVE_CYLINDER_Y) {
                cx = data.cyl_y.center_x;
                cz = data.cyl_y.center_z;
                radius = data.cyl_y.radius;
            } else {
                cx = data.cyl_z.center_x;
                cy = data.cyl_z.center_y;
                radius = data.cyl_z.radius;
            }
            const long double axes[3] = {
                fmaxl(fabsl(box->min_x - cx), fabsl(box->max_x - cx)),
                fmaxl(fabsl(box->min_y - cy), fabsl(box->max_y - cy)),
                fmaxl(fabsl(box->min_z - cz), fabsl(box->max_z - cz))};
            scale += axes[0] * axes[0] + axes[1] * axes[1] +
                     axes[2] * axes[2] + (long double)radius * radius;
        }
        const long double margin = 1024.0L * DBL_EPSILON * scale;
        if (!isfinite(interval.min) || !isfinite(interval.max) ||
            !isfinite(margin)) return 0;
        int inside = interval.max < -margin ? 1
            : interval.min > margin ? -1 : 0;
        if (inside && ((node->primitive.sense > 0) !=
                       (node->primitive.inverted != 0)))
            inside = -inside;
        return inside;
    }
    if (op != ALEA_OP_UNION && op != ALEA_OP_INTERSECTION &&
        op != ALEA_OP_DIFFERENCE && op != ALEA_OP_COMPLEMENT) return -2;
    const int left = slice_error_occurrence_node_constant(
        sys, node->operation.left, box, depth + 1, work);
    if (left == -2) return -2;
    if (op == ALEA_OP_COMPLEMENT) return -left;
    const int right = slice_error_occurrence_node_constant(
        sys, node->operation.right, box, depth + 1, work);
    if (right == -2) return -2;
    if (op == ALEA_OP_UNION)
        return left == 1 || right == 1 ? 1
            : left == -1 && right == -1 ? -1 : 0;
    if (op == ALEA_OP_INTERSECTION)
        return left == -1 || right == -1 ? -1
            : left == 1 && right == 1 ? 1 : 0;
    return left == -1 || right == 1 ? -1
        : left == 1 && right == -1 ? 1 : 0;
}

/* Certify that a rectangle remains in one concrete rectangular lattice
 * element, or that it misses the element. A crossing is indeterminate here. */
static int slice_error_rect_element_relation(
    const alea_cell_entry_t* cell,
    const alea_hier_spatial_chain_hit_t* hit, size_t level,
    const alea_bbox_t* box) {
    if (cell->lat_type != 1) return 0;
    const int64_t dims[3] = {
        (int64_t)cell->lat_fill_dims[1] - cell->lat_fill_dims[0] + 1,
        (int64_t)cell->lat_fill_dims[3] - cell->lat_fill_dims[2] + 1,
        (int64_t)cell->lat_fill_dims[5] - cell->lat_fill_dims[4] + 1};
    const int indices[3] = {hit->ancestor_lattice_i[level],
                            hit->ancestor_lattice_j[level],
                            hit->ancestor_lattice_k[level]};
    const double box_min[3] = {box->min_x, box->min_y, box->min_z};
    const double box_max[3] = {box->max_x, box->max_y, box->max_z};
    for (int axis = 0; axis < 3; ++axis) {
        if (dims[axis] <= 0) return 0;
        const double pitch = cell->lat_pitch[axis];
        const int tiled = axis < 2 || dims[axis] > 1 ||
            (cell->lat_fill_repeating && pitch > 0.0);
        if (!tiled && pitch <= 0.0) continue;
        if (!(pitch > 0.0) || !isfinite(pitch)) return 0;
        const long double offset = (long double)indices[axis] -
            (long double)cell->lat_fill_dims[2 * axis];
        const long double lo = (long double)cell->lat_lower_left[axis] +
            offset * pitch;
        const long double hi = lo + pitch;
        if (!isfinite(lo) || !isfinite(hi)) return 0;
        if ((long double)box_max[axis] < lo ||
            (long double)box_min[axis] > hi) return -1;
        if (!((long double)box_min[axis] > lo &&
              (long double)box_max[axis] < hi)) return 0;
    }
    return 1;
}

typedef struct {
    const alea_slice_error_query_t* query;
    const alea_transition_slice_critical_tile_t* tile;
    alea_slice_error_unresolved_reason_t reason;
    size_t work;
    int owner_ids[ALEA_SLICE_ERROR_OWNER_CAPACITY];
    size_t owner_count;
} slice_error_constant_occurrence_visit_t;

static int slice_error_visit_constant_occurrence(
    const alea_hier_spatial_chain_hit_t* hit, void* userdata) {
    slice_error_constant_occurrence_visit_t* visit = userdata;
    if (!hit->hit.is_terminal) return 0;
    if (hit->chain_truncated ||
        hit->hit.cell_index >= visit->query->sys->cells.count) {
        visit->reason = ALEA_SLICE_ERROR_UNRESOLVED_OCCURRENCE;
        return 1;
    }
    int contained = 1;
    for (size_t level = 0; level <= hit->ancestor_count; ++level) {
        const uint32_t cell_index = level == hit->ancestor_count
            ? hit->hit.cell_index : hit->ancestor_cell_indices[level];
        if (cell_index >= visit->query->sys->cells.count) {
            visit->reason = ALEA_SLICE_ERROR_UNRESOLVED_OCCURRENCE;
            return 1;
        }
        const alea_cell_entry_t* cell =
            &visit->query->sys->cells.data[cell_index];
        if (cell->original_root_node_id != ALEA_NODE_ID_INVALID) {
            visit->reason = ALEA_SLICE_ERROR_UNRESOLVED_OCCURRENCE;
            return 1;
        }
        alea_bbox_t local_box;
        const alea_matrix_t* transform = level == hit->ancestor_count
            ? &hit->hit.transform : &hit->ancestor_transforms[level];
        if (!slice_error_occurrence_tile_box(
                &visit->query->options.view, visit->tile,
                transform, &local_box)) {
            visit->reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 1;
        }
        if (level < hit->ancestor_count &&
            hit->ancestor_is_lattice[level]) {
            const int relation = slice_error_rect_element_relation(
                cell, hit, level, &local_box);
            if (relation == 0) {
                visit->reason = ALEA_SLICE_ERROR_UNRESOLVED_OCCURRENCE;
                return 1;
            }
            if (relation < 0) { contained = 0; break; }
            continue;
        }
        if (cell->lat_type != 0) {
            visit->reason = ALEA_SLICE_ERROR_UNRESOLVED_OCCURRENCE;
            return 1;
        }
        const int relation = slice_error_occurrence_node_constant(
            visit->query->sys, cell->root_node_id, &local_box, 0,
            &visit->work);
        if (relation == -2) {
            visit->reason = visit->work
                ? ALEA_SLICE_ERROR_UNRESOLVED_PRIMITIVE
                : ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
            return 1;
        }
        if (relation == 0) {
            visit->reason = ALEA_SLICE_ERROR_UNRESOLVED_PLANAR_ARRANGEMENT;
            return 1;
        }
        if (relation < 0) { contained = 0; break; }
    }
    if (!contained) return 0;
    if (visit->owner_count == ALEA_SLICE_ERROR_OWNER_CAPACITY) {
        visit->reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
        return 1;
    }
    size_t position = visit->owner_count;
    while (position > 0 &&
           visit->owner_ids[position - 1] > hit->hit.cell_id) {
        visit->owner_ids[position] = visit->owner_ids[position - 1];
        position--;
    }
    visit->owner_ids[position] = hit->hit.cell_id;
    visit->owner_count++;
    return 0;
}

/* Certify entire pages whose relevant occurrence CSG is constant. The chain
 * query must be complete; a single varying primitive leaves the page for the
 * boundary classifier. Lattice element extents need their own proof. */
static int slice_error_classify_constant_occurrences(
    const alea_slice_error_query_t* query,
    const alea_transition_slice_critical_tile_t* tile,
    size_t contextual_bytes, alea_slice_error_page_t* out,
    alea_slice_error_unresolved_reason_t* reason,
    int* output_omitted, size_t* peak_scratch,
    int* classified_owner_ids, size_t* classified_owner_count) {
    const alea_slice_error_query_options_t* options = &query->options;
    const size_t capacity =
        options->scan_options.max_exhaustive_occurrence_hits;
    if (!capacity) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
        return 0;
    }
    *peak_scratch = sizeof(slice_error_constant_occurrence_visit_t);
    alea_bbox_t world_box;
    if (!slice_error_occurrence_tile_box(&options->view, tile, NULL,
                                         &world_box)) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
        return 0;
    }
    slice_error_constant_occurrence_visit_t visit = {
        .query = query, .tile = tile,
        .reason = ALEA_SLICE_ERROR_UNRESOLVED_OCCURRENCE,
        .work = options->scan_options.max_active_boundary_tests};
    alea_hier_region_chain_status_t traversal =
        ALEA_HIER_REGION_CHAIN_UNSUPPORTED;
    size_t hit_count = 0;
    if (alea_hier_spatial_visit_region_chain_bounded(
            query->sys, &world_box, slice_error_visit_constant_occurrence,
            &visit, capacity, &hit_count, &traversal) != 0)
        return -1;
    if (traversal != ALEA_HIER_REGION_CHAIN_COMPLETE) {
        *reason = traversal == ALEA_HIER_REGION_CHAIN_VISITOR_STOPPED
            ? visit.reason
            : traversal == ALEA_HIER_REGION_CHAIN_MAX_HITS
                ? ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT
                : ALEA_SLICE_ERROR_UNRESOLVED_OCCURRENCE;
        return 0;
    }
    const size_t owner_count = visit.owner_count;
    const int* owner_ids = visit.owner_ids;
    if (out && owner_count != 1) {
        if (contextual_bytes > options->scan_options.max_output_bytes ||
            sizeof(alea_slice_error_region_t) >
                options->scan_options.max_output_bytes - contextual_bytes) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
            *output_omitted = 1;
            return 0;
        }
        out->regions = calloc(1, sizeof(*out->regions));
        if (!out->regions) return -1;
        out->region_count = 1;
        alea_slice_error_region_t* region = &out->regions[0];
        memcpy(region->uv_min, tile->uv_min, sizeof(tile->uv_min));
        memcpy(region->uv_max, tile->uv_max, sizeof(tile->uv_max));
        region->kind = owner_count == 0 ? ALEA_POINT_COVERAGE_GAP
                                         : ALEA_POINT_COVERAGE_OVERLAP;
        region->owner_count = owner_count;
        memcpy(region->owner_cell_ids, owner_ids,
               owner_count * sizeof(*owner_ids));
    }
    if (classified_owner_count) {
        *classified_owner_count = owner_count;
        if (classified_owner_ids)
            memcpy(classified_owner_ids, owner_ids,
                   owner_count * sizeof(*owner_ids));
    }
    *reason = ALEA_SLICE_ERROR_RESOLVED;
    return 1;
}

static int slice_error_axis_grid_valid(
    const slice_error_grid_line_t* lines, size_t count) {
    for (size_t i = 1; i < count; ++i) {
        const double a = lines[i - 1].coordinate;
        const double b = lines[i].coordinate;
        const double margin = lines[i - 1].uncertainty +
            lines[i].uncertainty +
            4.0 * DBL_EPSILON * fmax(1.0, fmax(fabs(a), fabs(b)));
        if (!(b > a) || !(b - a > margin)) return 0;
    }
    return 1;
}

static int slice_error_axis_is_defect(alea_point_coverage_kind_t kind) {
    return kind == ALEA_POINT_COVERAGE_GAP ||
           kind == ALEA_POINT_COVERAGE_OVERLAP;
}

typedef struct {
    size_t count;
    double uv[ALEA_SLICE_ERROR_POLYGON_CAPACITY][2];
    double uncertainty[ALEA_SLICE_ERROR_POLYGON_CAPACITY][2];
} slice_error_polygon_t;

static int slice_error_oblique_node_inside(
    const alea_system_t* sys, alea_node_id_t id,
    const uint32_t* primitive_ids, const int* raw_negative,
    size_t plane_count, const slice_error_axis_plane_t* constant_planes,
    size_t constant_count, size_t depth, size_t* work, int* inside) {
    if (!*work || depth > 128 || id >= sys->nodes.count) return 0;
    --*work;
    const alea_node_t* node = &sys->nodes.data[id];
    const alea_operation_t op = ALEA_GET_OPERATION(node);
    if (op == ALEA_OP_PRIMITIVE) {
        size_t plane = 0;
        while (plane < plane_count &&
               primitive_ids[plane] != node->primitive.primitive_id)
            ++plane;
        int negative;
        if (plane < plane_count) {
            negative = raw_negative[plane];
        } else {
            size_t constant = 0;
            while (constant < constant_count &&
                   (constant_planes[constant].axis != 2 ||
                    constant_planes[constant].primitive_id !=
                        node->primitive.primitive_id))
                ++constant;
            if (constant == constant_count) return 0;
            negative = constant_planes[constant].coefficient_sign < 0;
        }
        const int flip = (node->primitive.sense > 0) !=
                         (node->primitive.inverted != 0);
        *inside = negative != flip;
        return 1;
    }
    if (op != ALEA_OP_UNION && op != ALEA_OP_INTERSECTION &&
        op != ALEA_OP_DIFFERENCE && op != ALEA_OP_COMPLEMENT) return 0;
    int left = 0, right = 0;
    if (!slice_error_oblique_node_inside(
            sys, node->operation.left, primitive_ids, raw_negative,
            plane_count, constant_planes, constant_count,
            depth + 1, work, &left)) return 0;
    if (op == ALEA_OP_COMPLEMENT) {
        *inside = !left;
        return 1;
    }
    if (!slice_error_oblique_node_inside(
            sys, node->operation.right, primitive_ids, raw_negative,
            plane_count, constant_planes, constant_count,
            depth + 1, work, &right)) return 0;
    *inside = op == ALEA_OP_UNION ? left || right
        : op == ALEA_OP_INTERSECTION ? left && right
        : left && !right;
    return 1;
}

static int slice_error_polygon_add(slice_error_polygon_t* polygon,
                                   double u, double v,
                                   double uncertainty) {
    if (polygon->count == ALEA_SLICE_ERROR_POLYGON_CAPACITY ||
        !isfinite(u) || !isfinite(v) || !isfinite(uncertainty)) return 0;
    const size_t index = polygon->count++;
    polygon->uv[index][0] = u;
    polygon->uv[index][1] = v;
    polygon->uncertainty[index][0] = uncertainty;
    polygon->uncertainty[index][1] = uncertainty;
    return 1;
}

static void slice_error_project_plane(const alea_plane_data_t* plane,
                                      const alea_slice_plane_t* frame,
                                      long double* a, long double* b,
                                      long double* c) {
    const long double normal[3] = {
        plane->a, plane->b, plane->c
    };
    *a = *b = 0.0L;
    *c = plane->d;
    for (int axis = 0; axis < 3; ++axis) {
        *a += normal[axis] * frame->u_axis[axis];
        *b += normal[axis] * frame->v_axis[axis];
        *c += normal[axis] * frame->origin[axis];
    }
}

static int slice_error_world_point(const alea_slice_plane_t* frame,
                                   double u, double v, double world[3]) {
    for (int axis = 0; axis < 3; ++axis) {
        world[axis] = frame->origin[axis] +
            frame->u_axis[axis] * u + frame->v_axis[axis] * v;
        if (!isfinite(world[axis])) return 0;
    }
    return 1;
}

static void slice_error_mixed_owners(
    const alea_cell_entry_t* plane_cell,
    const alea_cell_entry_t* sphere_cell,
    int plane_owned, int sphere_owned,
    alea_point_coverage_kind_t* kind, size_t* count, int ids[2]) {
    *count = 0;
    if (plane_owned) ids[(*count)++] = plane_cell->mc_cell_id;
    if (sphere_owned) ids[(*count)++] = sphere_cell->mc_cell_id;
    *kind = *count == 0 ? ALEA_POINT_COVERAGE_GAP
        : *count == 1 ? ALEA_POINT_COVERAGE_UNIQUE
        : ALEA_POINT_COVERAGE_OVERLAP;
}

/* A line and a circle have at most two intersections. Both crossings must be
 * safely inside the page. The circular arcs may cross the page edge and are
 * clipped by the consumer to the receipt's core rectangle. */
static int slice_error_classify_mixed_line_circle(
    const alea_system_t* sys, const alea_slice_view_t* view,
    const alea_transition_slice_critical_tile_t* tile,
    const size_t cells[2], size_t scratch_available, size_t* scratch_used,
    size_t contextual_bytes, size_t output_limit,
    alea_slice_error_page_t* out) {
    const long double tau = 6.283185307179586476925286766559005768L;
    *scratch_used = 0;
    const alea_cell_entry_t* cell[2] = {
        &sys->cells.data[cells[0]], &sys->cells.data[cells[1]]
    };
    const alea_node_t* node[2];
    const alea_primitive_entry_t* primitive[2];
    int sphere_index = -1, plane_index = -1;
    for (int i = 0; i < 2; ++i) {
        if (cell[i]->fill_universe > 0 || cell[i]->lat_type != 0 ||
            cell[i]->original_root_node_id != ALEA_NODE_ID_INVALID ||
            cell[i]->root_node_id >= sys->nodes.count) return 0;
        node[i] = &sys->nodes.data[cell[i]->root_node_id];
        if (ALEA_GET_OPERATION(node[i]) != ALEA_OP_PRIMITIVE ||
            node[i]->primitive.primitive_id >= sys->primitives.count)
            return 0;
        primitive[i] = &sys->primitives.data[
            node[i]->primitive.primitive_id];
        if (primitive[i]->type == ALEA_PRIMITIVE_SPHERE ||
            primitive[i]->type == ALEA_PRIMITIVE_SPH) {
            if (sphere_index >= 0) return 0;
            sphere_index = i;
        } else if (primitive[i]->type == ALEA_PRIMITIVE_PLANE &&
                   primitive[i]->payload_index <
                       sys->primitive_planes.count) {
            if (plane_index >= 0) return 0;
            plane_index = i;
        } else return 0;
    }
    if (sphere_index < 0 || plane_index < 0) return 0;
    if (scratch_available < sizeof(alea_slice_error_circle_t)) return 3;
    const alea_cell_entry_t* scell = cell[sphere_index];
    const alea_cell_entry_t* pcell = cell[plane_index];
    alea_slice_error_page_t source = {0};
    int status = slice_error_classify_isolated_circle(
        sys, view, tile, cells[sphere_index], 0, SIZE_MAX, &source);
    if (status != 1) return status < 0 ? -1 : 0;
    status = 0;
    *scratch_used = sizeof(alea_slice_error_circle_t);
    const alea_slice_error_circle_t* circle = source.circles;
    long double a, b, c;
    const alea_plane_data_t* p = &sys->primitive_planes.data[
        primitive[plane_index]->payload_index];
    slice_error_project_plane(p, &view->plane, &a, &b, &c);
    const long double norm = sqrtl(a * a + b * b);
    if (!(norm > 0.0L) || !isfinite(norm)) goto done;
    a /= norm; b /= norm; c /= norm;
    int surface_id = 0;
    for (size_t i = 0; i < sys->surfaces.count; ++i) {
        if (sys->surfaces.data[i].primitive_id !=
            node[plane_index]->primitive.primitive_id) continue;
        if (surface_id || sys->surfaces.data[i].mc_surface_id <= 0)
            goto done;
        surface_id = sys->surfaces.data[i].mc_surface_id;
    }
    if (!surface_id) goto done;
    const long double scale = 1.0L + fabsl(c) +
        fmaxl(fabsl(tile->uv_min[0]), fabsl(tile->uv_max[0])) +
        fmaxl(fabsl(tile->uv_min[1]), fabsl(tile->uv_max[1])) +
        circle->radius;
    const long double margin = 4096.0L * DBL_EPSILON * scale +
        32.0L * circle->geometry_uncertainty;
    if (!isfinite(margin) || !(margin < circle->radius / 64.0L))
        goto done;
    const long double corner[4][2] = {
        {tile->uv_min[0], tile->uv_min[1]},
        {tile->uv_max[0], tile->uv_min[1]},
        {tile->uv_max[0], tile->uv_max[1]},
        {tile->uv_min[0], tile->uv_max[1]}
    };
    long double endpoint[2][2];
    int endpoint_count = 0;
    int first_corner_negative = -1;
    for (int i = 0; i < 4; ++i) {
        const int j = (i + 1) % 4;
        const long double s0 = a * corner[i][0] + b * corner[i][1] + c;
        const long double s1 = a * corner[j][0] + b * corner[j][1] + c;
        if (!isfinite(s0) || !isfinite(s1) ||
            fabsl(s0) <= margin || fabsl(s1) <= margin) goto done;
        if (first_corner_negative < 0)
            first_corner_negative = s0 < 0.0L;
        if ((s0 < 0.0L) == (s1 < 0.0L)) continue;
        if (endpoint_count == 2) goto done;
        const long double t = s0 / (s0 - s1);
        for (int axis = 0; axis < 2; ++axis)
            endpoint[endpoint_count][axis] = corner[i][axis] +
                t * (corner[j][axis] - corner[i][axis]);
        ++endpoint_count;
    }
    if (endpoint_count != 0 && endpoint_count != 2) goto done;
    const int line_visible = endpoint_count == 2;
    if (!line_visible) {
        for (int i = 0; i < 4; ++i)
            if ((a * corner[i][0] + b * corner[i][1] + c < 0.0L) !=
                first_corner_negative) goto done;
        endpoint[0][0] = -c * a + b;
        endpoint[0][1] = -c * b - a;
        endpoint[1][0] = -c * a - b;
        endpoint[1][1] = -c * b + a;
    }
    const long double center_u = circle->center_uv[0];
    const long double center_v = circle->center_uv[1];
    const long double signed_distance = a * center_u + b * center_v + c;
    const long double radius = circle->radius;
    if (!isfinite(signed_distance) ||
        !(fabsl(fabsl(signed_distance) - radius) >
          64.0L * margin)) goto done;
    const int intersection_count =
        fabsl(signed_distance) < radius ? 2 : 0;
    long double crossing[2][2] = {{0.0L, 0.0L}, {0.0L, 0.0L}};
    long double crossing_error = 64.0L * margin;
    if (intersection_count) {
        const long double height2 = radius * radius -
            signed_distance * signed_distance;
        const long double height = sqrtl(height2);
        if (!(height > 64.0L * margin) || !isfinite(height))
            goto done;
        crossing_error *= 1.0L + radius / height;
        if (!(crossing_error < height / 16.0L) ||
            !isfinite(crossing_error)) goto done;
        const long double foot[2] = {
            center_u - signed_distance * a,
            center_v - signed_distance * b
        };
        const long double tangent[2] = {-b, a};
        for (int i = 0; i < 2; ++i)
            for (int axis = 0; axis < 2; ++axis)
                crossing[i][axis] = foot[axis] +
                    (i ? height : -height) * tangent[axis];
    }
    const long double direction[2] = {
        endpoint[1][0] - endpoint[0][0],
        endpoint[1][1] - endpoint[0][1]
    };
    const long double length2 = direction[0] * direction[0] +
        direction[1] * direction[1];
    if (!(length2 > 0.0L) || !isfinite(length2)) goto done;
    long double crossing_t[2] = {0.0L, 0.0L};
    long double angle[2] = {0.0L, 0.0L};
    int visible_crossing_count = 0;
    long double visible_t[2];
    long double visible_crossing[2][2];
    for (int i = 0; i < intersection_count; ++i) {
        crossing_t[i] = ((crossing[i][0] - endpoint[0][0]) *
                         direction[0] +
                         (crossing[i][1] - endpoint[0][1]) *
                         direction[1]) / length2;
        const long double parameter_margin =
            crossing_error / sqrtl(length2);
        if (!isfinite(crossing_t[i]) ||
            (line_visible &&
             (fabsl(crossing_t[i]) <= parameter_margin ||
              fabsl(crossing_t[i] - 1.0L) <= parameter_margin)))
            goto done;
        angle[i] = atan2l(crossing[i][1] - center_v,
                          crossing[i][0] - center_u);
        if (angle[i] < 0.0L) angle[i] += tau;
        if (line_visible && crossing_t[i] > 0.0L &&
            crossing_t[i] < 1.0L) {
            for (int axis = 0; axis < 2; ++axis)
                if (!(crossing[i][axis] > tile->uv_min[axis] +
                          crossing_error) ||
                    !(crossing[i][axis] < tile->uv_max[axis] -
                          crossing_error)) goto done;
            visible_t[visible_crossing_count] = crossing_t[i];
            memcpy(visible_crossing[visible_crossing_count], crossing[i],
                   sizeof(crossing[i]));
            ++visible_crossing_count;
        }
    }
    if (visible_crossing_count == 2 && visible_t[0] > visible_t[1]) {
        const long double t = visible_t[0];
        visible_t[0] = visible_t[1]; visible_t[1] = t;
        for (int axis = 0; axis < 2; ++axis) {
            const long double x = visible_crossing[0][axis];
            visible_crossing[0][axis] = visible_crossing[1][axis];
            visible_crossing[1][axis] = x;
        }
    }
    if (intersection_count && angle[0] > angle[1]) {
        const long double t = angle[0];
        angle[0] = angle[1]; angle[1] = t;
    }
    if (visible_crossing_count == 2 &&
        !(visible_t[1] - visible_t[0] >
          16.0L * crossing_error / sqrtl(length2))) goto done;
    if (intersection_count &&
        (!(angle[1] - angle[0] >
           16.0L * crossing_error / radius) ||
         !(tau - angle[1] + angle[0] >
           16.0L * crossing_error / radius))) goto done;
    const size_t arc_count = intersection_count ? 2u : 1u;
    const size_t interval_count = line_visible
        ? (size_t)visible_crossing_count + 1u : 0u;
    if (contextual_bytes > output_limit ||
        arc_count * sizeof(alea_slice_error_circle_t) +
        interval_count * sizeof(alea_slice_error_interval_t) >
            output_limit - contextual_bytes) {
        status = 2;
        goto done;
    }
    out->circles = calloc(arc_count, sizeof(*out->circles));
    out->intervals = interval_count
        ? calloc(interval_count, sizeof(*out->intervals)) : NULL;
    if (!out->circles || (interval_count && !out->intervals)) {
        status = -1;
        goto done;
    }
    const int plane_negative_owned =
        !((node[plane_index]->primitive.sense > 0) !=
          (node[plane_index]->primitive.inverted != 0));
    const int sphere_inside_owned =
        circle->inside_owner_count == 1;
    for (size_t part = 0; part < arc_count; ++part) {
        const long double start = !intersection_count ? 0.0L
            : part == 0 ? angle[0] : angle[1];
        const long double end = !intersection_count ? tau
            : part == 0 ? angle[1] : angle[0] + tau;
        const long double middle = (start + end) * 0.5L;
        const long double u = center_u + radius * cosl(middle);
        const long double v = center_v + radius * sinl(middle);
        const long double line_distance = a * u + b * v + c;
        if (!(fabsl(line_distance) > 16.0L * margin)) goto done;
        const int plane_owned = (line_distance < 0.0L) ==
            plane_negative_owned;
        alea_slice_error_circle_t* arc = &out->circles[part];
        *arc = *circle;
        arc->start_angle = (double)start;
        arc->end_angle = (double)end;
        arc->geometry_uncertainty = (double)crossing_error;
        int ids[2];
        slice_error_mixed_owners(pcell, scell, plane_owned,
            sphere_inside_owned, &arc->inside_kind,
            &arc->inside_owner_count, ids);
        memcpy(arc->inside_owner_cell_ids, ids,
               arc->inside_owner_count * sizeof(int));
        slice_error_mixed_owners(pcell, scell, plane_owned,
            !sphere_inside_owned, &arc->outside_kind,
            &arc->outside_owner_count, ids);
        memcpy(arc->outside_owner_cell_ids, ids,
               arc->outside_owner_count * sizeof(int));
        for (int side = 0; side < 2; ++side) {
            const long double radial = radius +
                (side ? 32.0L * margin : -32.0L * margin);
            double world[3];
            if (!slice_error_world_point(&view->plane,
                    (double)(center_u + radial * cosl(middle)),
                    (double)(center_v + radial * sinl(middle)),
                    world) ||
                alea_point_inside(sys, scell->root_node_id,
                    world[0], world[1], world[2]) !=
                    (side ? !sphere_inside_owned : sphere_inside_owned) ||
                alea_point_inside(sys, pcell->root_node_id,
                    world[0], world[1], world[2]) != plane_owned)
                goto done;
        }
    }
    long double t[4] = {0.0L, 0.0L, 0.0L, 0.0L};
    for (int i = 0; i < visible_crossing_count; ++i)
        t[i + 1] = visible_t[i];
    if (interval_count) t[interval_count] = 1.0L;
    for (size_t part = 0; part < interval_count; ++part) {
        const long double middle = (t[part] + t[part + 1]) * 0.5L;
        const long double u = endpoint[0][0] + middle * direction[0];
        const long double v = endpoint[0][1] + middle * direction[1];
        const long double sphere_distance =
            sqrtl((u - center_u) * (u - center_u) +
                  (v - center_v) * (v - center_v));
        if (!(fabsl(sphere_distance - radius) > 16.0L * margin))
            goto done;
        const int sphere_owned = (sphere_distance < radius) ==
            sphere_inside_owned;
        alea_slice_error_interval_t* interval = &out->intervals[part];
        interval->evidence_scope =
            ALEA_SLICE_BOUNDARY_EVIDENCE_VERIFIED_INTERVAL;
        interval->surface_id = surface_id;
        interval->primitive_id = node[plane_index]->primitive.primitive_id;
        interval->axis = -1;
        const long double* first = part == 0 ? endpoint[0] :
            visible_crossing[part - 1];
        const long double* last = part + 1 == interval_count
            ? endpoint[1] : visible_crossing[part];
        for (int axis = 0; axis < 2; ++axis) {
            interval->uv_start[axis] = (double)first[axis];
            interval->uv_end[axis] = (double)last[axis];
        }
        interval->endpoint_uncertainty[0] =
            (double)(part == 0 ? margin : crossing_error);
        interval->endpoint_uncertainty[1] =
            (double)(part + 1 == interval_count ? margin : crossing_error);
        int ids[2];
        slice_error_mixed_owners(pcell, scell, plane_negative_owned,
            sphere_owned, &interval->negative_side_kind,
            &interval->negative_owner_count, ids);
        memcpy(interval->negative_owner_cell_ids, ids,
               interval->negative_owner_count * sizeof(int));
        slice_error_mixed_owners(pcell, scell, !plane_negative_owned,
            sphere_owned, &interval->positive_side_kind,
            &interval->positive_owner_count, ids);
        memcpy(interval->positive_owner_cell_ids, ids,
               interval->positive_owner_count * sizeof(int));
        for (int side = 0; side < 2; ++side) {
            double world[3];
            const long double offset = side ? 32.0L * margin
                                            : -32.0L * margin;
            if (!slice_error_world_point(&view->plane,
                    (double)(u + offset * a),
                    (double)(v + offset * b), world) ||
                alea_point_inside(sys, pcell->root_node_id,
                    world[0], world[1], world[2]) !=
                    (side ? !plane_negative_owned :
                            plane_negative_owned) ||
                alea_point_inside(sys, scell->root_node_id,
                    world[0], world[1], world[2]) != sphere_owned)
                goto done;
        }
    }
    out->circle_count = arc_count;
    out->interval_count = interval_count;
    status = 1;
done:
    if (status != 1) {
        free(out->circles); out->circles = NULL;
        free(out->intervals); out->intervals = NULL;
    }
    free(source.circles);
    return status;
}

/* A single oblique line partitions a rectangular core into two convex faces.
 * This deliberately excludes line corners and multiple-line arrangements. */
static int slice_error_classify_single_oblique_tile(
    const alea_slice_error_query_t* query,
    const alea_transition_slice_critical_tile_t* tile,
    const size_t* cells, size_t cell_count,
    const slice_error_axis_plane_t* selected,
    const slice_error_axis_plane_t* constant_planes,
    size_t constant_count,
    size_t contextual_bytes, size_t* work,
    alea_slice_error_page_t* out,
    alea_slice_error_unresolved_reason_t* reason,
    int* output_omitted) {
    const alea_slice_error_query_options_t* o = &query->options;
    const alea_system_t* sys = query->sys;
    if (selected->surface_id <= 0) return 0;
    const alea_primitive_entry_t* primitive =
        &sys->primitives.data[selected->primitive_id];
    if (primitive->type != ALEA_PRIMITIVE_PLANE ||
        primitive->payload_index >= sys->primitive_planes.count)
        return 0;
    const alea_plane_data_t* p =
        &sys->primitive_planes.data[primitive->payload_index];
    if (!isfinite(p->a) || !isfinite(p->b) ||
        !isfinite(p->c) || !isfinite(p->d)) return 0;
    const alea_slice_plane_t* frame = &o->view.plane;
    long double a, b, c;
    slice_error_project_plane(p, frame, &a, &b, &c);
    if (a == 0.0L && b == 0.0L) return 0;
    const long double norm = fmaxl(fabsl(a), fabsl(b));
    if (!(norm > 0.0L) || !isfinite(norm)) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
        return 0;
    }
    a /= norm; b /= norm; c /= norm;
    const double corner[4][2] = {
        {tile->uv_min[0], tile->uv_min[1]},
        {tile->uv_max[0], tile->uv_min[1]},
        {tile->uv_max[0], tile->uv_max[1]},
        {tile->uv_min[0], tile->uv_max[1]}
    };
    long double value[4];
    const long double scale = 1.0L + fabsl(c) +
        fabsl(a) * fmaxl(fabsl(tile->uv_min[0]),
                         fabsl(tile->uv_max[0])) +
        fabsl(b) * fmaxl(fabsl(tile->uv_min[1]),
                         fabsl(tile->uv_max[1]));
    const long double tolerance = scale *
        (128.0L * LDBL_EPSILON + 16.0L * DBL_EPSILON);
    if (!isfinite(tolerance)) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
        return 0;
    }
    for (int i = 0; i < 4; ++i) {
        value[i] = a * corner[i][0] + b * corner[i][1] + c;
        if (!isfinite(value[i]) || fabsl(value[i]) <= tolerance) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
    }
    slice_error_polygon_t polygon[2] = {{0}, {0}};
    slice_error_polygon_t crossings = {0};
    for (int i = 0; i < 4; ++i) {
        const int next = (i + 1) & 3;
        const int side = value[i] > 0.0L;
        if (!slice_error_polygon_add(&polygon[side],
                corner[i][0], corner[i][1], 0.0)) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
            return 0;
        }
        if ((value[i] > 0.0L) == (value[next] > 0.0L)) continue;
        const long double t = value[i] / (value[i] - value[next]);
        const long double u = (long double)corner[i][0] + t *
            ((long double)corner[next][0] - corner[i][0]);
        const long double v = (long double)corner[i][1] + t *
            ((long double)corner[next][1] - corner[i][1]);
        const double du = (double)u, dv = (double)v;
        const double uncertainty = (double)(
            4.0L * tolerance / fabsl(value[i] - value[next]) *
            fmaxl(fabsl((long double)corner[next][0] - corner[i][0]),
                  fabsl((long double)corner[next][1] - corner[i][1])) +
            8.0L * DBL_EPSILON *
                fmaxl(1.0L, fmaxl(fabsl(u), fabsl(v))));
        const long double edge_length = fmaxl(
            fabsl((long double)corner[next][0] - corner[i][0]),
            fabsl((long double)corner[next][1] - corner[i][1]));
        const long double corner_distance =
            fminl(t, 1.0L - t) * edge_length;
        if (!(t > 0.0L && t < 1.0L) || !isfinite(du) ||
            !isfinite(dv) || !isfinite(uncertainty) ||
            !((long double)uncertainty * 4.0L < corner_distance) ||
            !slice_error_polygon_add(&crossings, du, dv, uncertainty) ||
            !slice_error_polygon_add(&polygon[0], du, dv, uncertainty) ||
            !slice_error_polygon_add(&polygon[1], du, dv, uncertainty)) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
    }
    if (crossings.count != 0 && crossings.count != 2) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
        return 0;
    }
    slice_error_face_t face[2] = {{0}, {0}};
    for (int side = 0; side < 2; ++side) {
        if (polygon[side].count == 0) continue;
        if (polygon[side].count < 3) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
        double u = 0.0, v = 0.0;
        for (size_t j = 0; j < polygon[side].count; ++j) {
            u += polygon[side].uv[j][0] / polygon[side].count;
            v += polygon[side].uv[j][1] / polygon[side].count;
        }
        const long double residual = a * (long double)u +
                                     b * (long double)v + c;
        if (!isfinite(u) || !isfinite(v) ||
            (side ? residual <= tolerance : residual >= -tolerance)) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
        double world[3];
        if (!slice_error_world_point(frame, u, v, world)) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
        for (size_t ci = 0; ci < cell_count; ++ci) {
            const alea_cell_entry_t* cell = &sys->cells.data[cells[ci]];
            int inside = 0;
            const uint32_t primitive_id = selected->primitive_id;
            const int raw_negative = side == 0;
            if (!slice_error_oblique_node_inside(sys, cell->root_node_id,
                    &primitive_id, &raw_negative, 1,
                    constant_planes, constant_count, 0, work, &inside)) {
                *reason = *work == 0
                    ? ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT
                    : ALEA_SLICE_ERROR_UNRESOLVED_UNSUPPORTED_GEOMETRY;
                return 0;
            }
            if (inside != alea_point_inside(sys, cell->root_node_id,
                    world[0], world[1], world[2])) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                return 0;
            }
            if (inside) {
                if (face[side].owner_count == ALEA_SLICE_ERROR_OWNER_CAPACITY) {
                    *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
                    return 0;
                }
                face[side].owner_cell_ids[face[side].owner_count++] =
                    cell->mc_cell_id;
            }
        }
        face[side].kind = face[side].owner_count == 0
            ? ALEA_POINT_COVERAGE_GAP
            : face[side].owner_count == 1
                ? ALEA_POINT_COVERAGE_UNIQUE
                : ALEA_POINT_COVERAGE_OVERLAP;
    }
    size_t region_count = 0;
    for (int side = 0; side < 2; ++side)
        if (polygon[side].count && slice_error_axis_is_defect(face[side].kind))
            ++region_count;
    const size_t interval_count = crossings.count == 2 &&
        slice_error_axis_is_defect(face[0].kind) !=
        slice_error_axis_is_defect(face[1].kind) ? 1u : 0u;
    const size_t output_limit = o->scan_options.max_output_bytes;
    if (contextual_bytes > output_limit ||
        region_count > (output_limit - contextual_bytes) /
            sizeof(alea_slice_error_region_t) ||
        interval_count > (output_limit - contextual_bytes -
            region_count * sizeof(alea_slice_error_region_t)) /
            sizeof(alea_slice_error_interval_t)) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
        *output_omitted = 1;
        return 0;
    }
    out->regions = region_count
        ? calloc(region_count, sizeof(*out->regions)) : NULL;
    out->intervals = interval_count
        ? calloc(interval_count, sizeof(*out->intervals)) : NULL;
    if ((region_count && !out->regions) ||
        (interval_count && !out->intervals)) return -1;
    for (int side = 0; side < 2; ++side) {
        if (!polygon[side].count ||
            !slice_error_axis_is_defect(face[side].kind)) continue;
        alea_slice_error_region_t* region =
            &out->regions[out->region_count++];
        region->kind = face[side].kind;
        region->owner_count = face[side].owner_count;
        memcpy(region->owner_cell_ids, face[side].owner_cell_ids,
               region->owner_count * sizeof(int));
        region->polygon_vertex_count = polygon[side].count;
        region->uv_min[0] = region->uv_min[1] = INFINITY;
        region->uv_max[0] = region->uv_max[1] = -INFINITY;
        for (size_t j = 0; j < polygon[side].count; ++j) {
            for (int axis = 0; axis < 2; ++axis) {
                const double coordinate = polygon[side].uv[j][axis];
                region->polygon_uv[j][axis] = coordinate;
                region->polygon_uv_uncertainty[j][axis] =
                    polygon[side].uncertainty[j][axis];
                region->uv_min[axis] = fmin(region->uv_min[axis], coordinate);
                region->uv_max[axis] = fmax(region->uv_max[axis], coordinate);
                if (region->uv_min[axis] == coordinate)
                    region->uv_min_uncertainty[axis] = fmax(
                        region->uv_min_uncertainty[axis],
                        polygon[side].uncertainty[j][axis]);
                if (region->uv_max[axis] == coordinate)
                    region->uv_max_uncertainty[axis] = fmax(
                        region->uv_max_uncertainty[axis],
                        polygon[side].uncertainty[j][axis]);
            }
        }
    }
    if (interval_count) {
        alea_slice_error_interval_t* interval = &out->intervals[0];
        interval->evidence_scope =
            ALEA_SLICE_BOUNDARY_EVIDENCE_VERIFIED_INTERVAL;
        interval->surface_id = selected->surface_id;
        interval->primitive_id = selected->primitive_id;
        interval->axis = -1;
        const int reverse = crossings.uv[0][0] > crossings.uv[1][0] ||
            (crossings.uv[0][0] == crossings.uv[1][0] &&
             crossings.uv[0][1] > crossings.uv[1][1]);
        const int first = reverse ? 1 : 0, last = reverse ? 0 : 1;
        memcpy(interval->uv_start, crossings.uv[first],
               sizeof(interval->uv_start));
        memcpy(interval->uv_end, crossings.uv[last],
               sizeof(interval->uv_end));
        interval->endpoint_uncertainty[0] =
            crossings.uncertainty[first][0];
        interval->endpoint_uncertainty[1] =
            crossings.uncertainty[last][0];
        interval->negative_side_kind = face[0].kind;
        interval->positive_side_kind = face[1].kind;
        interval->negative_owner_count = face[0].owner_count;
        interval->positive_owner_count = face[1].owner_count;
        memcpy(interval->negative_owner_cell_ids, face[0].owner_cell_ids,
               face[0].owner_count * sizeof(int));
        memcpy(interval->positive_owner_cell_ids, face[1].owner_cell_ids,
               face[1].owner_count * sizeof(int));
        out->interval_count = 1;
    }
    *reason = ALEA_SLICE_ERROR_RESOLVED;
    return 1;
}

typedef struct {
    int plane;
    int negative_face, positive_face;
    double uv[2][2];
    double uncertainty[2];
} slice_error_oblique_segment_t;

static int slice_error_publish_oblique_faces(
    const alea_slice_error_query_options_t* options,
    const slice_error_axis_plane_t* selected,
    const slice_error_polygon_t* polygon,
    const slice_error_face_t* face, size_t face_count,
    const slice_error_oblique_segment_t* segments, size_t segment_count,
    size_t contextual_bytes, alea_slice_error_page_t* out,
    alea_slice_error_unresolved_reason_t* reason, int* output_omitted) {
    size_t region_count = 0, interval_count = 0;
    for (size_t i = 0; i < segment_count; ++i)
        if (segments[i].negative_face < 0 ||
            segments[i].positive_face < 0 ||
            (size_t)segments[i].negative_face >= face_count ||
            (size_t)segments[i].positive_face >= face_count ||
            !polygon[segments[i].negative_face].count ||
            !polygon[segments[i].positive_face].count) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
    for (size_t i = 0; i < face_count; ++i)
        if (polygon[i].count && slice_error_axis_is_defect(face[i].kind))
            ++region_count;
    for (size_t i = 0; i < segment_count; ++i)
        if (slice_error_axis_is_defect(
                face[segments[i].negative_face].kind) !=
            slice_error_axis_is_defect(
                face[segments[i].positive_face].kind)) ++interval_count;
    const size_t output_limit = options->scan_options.max_output_bytes;
    if (contextual_bytes > output_limit ||
        region_count > (output_limit - contextual_bytes) /
            sizeof(alea_slice_error_region_t) ||
        interval_count > (output_limit - contextual_bytes -
            region_count * sizeof(alea_slice_error_region_t)) /
            sizeof(alea_slice_error_interval_t)) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
        *output_omitted = 1;
        return 0;
    }
    out->regions = region_count
        ? calloc(region_count, sizeof(*out->regions)) : NULL;
    out->intervals = interval_count
        ? calloc(interval_count, sizeof(*out->intervals)) : NULL;
    if ((region_count && !out->regions) ||
        (interval_count && !out->intervals)) return -1;
    for (size_t i = 0; i < face_count; ++i) {
        if (!polygon[i].count ||
            !slice_error_axis_is_defect(face[i].kind)) continue;
        alea_slice_error_region_t* region =
            &out->regions[out->region_count++];
        region->kind = face[i].kind;
        region->owner_count = face[i].owner_count;
        memcpy(region->owner_cell_ids, face[i].owner_cell_ids,
               region->owner_count * sizeof(int));
        region->polygon_vertex_count = polygon[i].count;
        region->uv_min[0] = region->uv_min[1] = INFINITY;
        region->uv_max[0] = region->uv_max[1] = -INFINITY;
        for (size_t j = 0; j < polygon[i].count; ++j)
            for (int axis = 0; axis < 2; ++axis) {
                const double coordinate = polygon[i].uv[j][axis];
                const double uncertainty = polygon[i].uncertainty[j][axis];
                region->polygon_uv[j][axis] = coordinate;
                region->polygon_uv_uncertainty[j][axis] = uncertainty;
                region->uv_min[axis] =
                    fmin(region->uv_min[axis], coordinate);
                region->uv_max[axis] =
                    fmax(region->uv_max[axis], coordinate);
                if (region->uv_min[axis] == coordinate)
                    region->uv_min_uncertainty[axis] = fmax(
                        region->uv_min_uncertainty[axis], uncertainty);
                if (region->uv_max[axis] == coordinate)
                    region->uv_max_uncertainty[axis] = fmax(
                        region->uv_max_uncertainty[axis], uncertainty);
            }
    }
    for (size_t i = 0; i < segment_count; ++i) {
        const slice_error_oblique_segment_t* segment = &segments[i];
        const slice_error_face_t* negative =
            &face[segment->negative_face];
        const slice_error_face_t* positive =
            &face[segment->positive_face];
        if (slice_error_axis_is_defect(negative->kind) ==
            slice_error_axis_is_defect(positive->kind)) continue;
        alea_slice_error_interval_t* interval =
            &out->intervals[out->interval_count++];
        interval->evidence_scope =
            ALEA_SLICE_BOUNDARY_EVIDENCE_VERIFIED_INTERVAL;
        interval->surface_id = selected[segment->plane].surface_id;
        interval->primitive_id = selected[segment->plane].primitive_id;
        interval->axis = -1;
        const int reverse = segment->uv[0][0] > segment->uv[1][0] ||
            (segment->uv[0][0] == segment->uv[1][0] &&
             segment->uv[0][1] > segment->uv[1][1]);
        const int first = reverse ? 1 : 0, last = reverse ? 0 : 1;
        memcpy(interval->uv_start, segment->uv[first],
               sizeof(interval->uv_start));
        memcpy(interval->uv_end, segment->uv[last],
               sizeof(interval->uv_end));
        interval->endpoint_uncertainty[0] =
            segment->uncertainty[first];
        interval->endpoint_uncertainty[1] =
            segment->uncertainty[last];
        interval->negative_side_kind = negative->kind;
        interval->positive_side_kind = positive->kind;
        interval->negative_owner_count = negative->owner_count;
        interval->positive_owner_count = positive->owner_count;
        memcpy(interval->negative_owner_cell_ids, negative->owner_cell_ids,
               negative->owner_count * sizeof(int));
        memcpy(interval->positive_owner_cell_ids, positive->owner_cell_ids,
               positive->owner_count * sizeof(int));
    }
    *reason = ALEA_SLICE_ERROR_RESOLVED;
    return 1;
}

/* Two parallel oblique planes divide the core into at most three convex
 * bands. The primitive signs are fixed in each band, including a thin band
 * between the lines. Nonparallel or nearly coincident lines remain
 * unresolved. */
static int slice_error_classify_parallel_oblique_tile(
    const alea_slice_error_query_t* query,
    const alea_transition_slice_critical_tile_t* tile,
    const size_t* cells, size_t cell_count,
    const slice_error_axis_plane_t selected[2],
    const slice_error_axis_plane_t* constant_planes,
    size_t constant_count,
    size_t contextual_bytes, size_t* work,
    alea_slice_error_page_t* out,
    alea_slice_error_unresolved_reason_t* reason,
    int* output_omitted) {
    const alea_slice_error_query_options_t* o = &query->options;
    const alea_system_t* sys = query->sys;
    const alea_slice_plane_t* frame = &o->view.plane;
    long double coefficient[2][2], threshold[2];
    double normal[2][3];
    int orientation[2] = {1, 0};
    uint32_t primitive_ids[2];
    for (int i = 0; i < 2; ++i) {
        if (selected[i].surface_id <= 0) return 0;
        primitive_ids[i] = selected[i].primitive_id;
        const alea_primitive_entry_t* primitive =
            &sys->primitives.data[primitive_ids[i]];
        if (primitive->type != ALEA_PRIMITIVE_PLANE ||
            primitive->payload_index >= sys->primitive_planes.count)
            return 0;
        const alea_plane_data_t* p =
            &sys->primitive_planes.data[primitive->payload_index];
        if (!isfinite(p->a) || !isfinite(p->b) || !isfinite(p->c) ||
            !isfinite(p->d)) return 0;
        normal[i][0] = p->a;
        normal[i][1] = p->b;
        normal[i][2] = p->c;
        long double a, b, c;
        slice_error_project_plane(p, frame, &a, &b, &c);
        const long double norm = fmaxl(fabsl(a), fabsl(b));
        if (!(norm > 0.0L) || !isfinite(norm)) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
        coefficient[i][0] = a / norm;
        coefficient[i][1] = b / norm;
        threshold[i] = -c / norm;
        if (!isfinite(threshold[i])) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
    }
    /* Exact matching normals avoid treating two nearly parallel lines as
     * disjoint bands when their intersection is hidden inside the tile. */
    if (normal[1][0] == normal[0][0] &&
        normal[1][1] == normal[0][1] &&
        normal[1][2] == normal[0][2])
        orientation[1] = 1;
    else if (normal[1][0] == -normal[0][0] &&
             normal[1][1] == -normal[0][1] &&
             normal[1][2] == -normal[0][2]) {
        orientation[1] = -1;
        threshold[1] = -threshold[1];
    } else return 0;
    const int line[2] = {
        threshold[0] < threshold[1] ? 0 : 1,
        threshold[0] < threshold[1] ? 1 : 0
    };
    const long double bound[2] = {
        threshold[line[0]], threshold[line[1]]
    };
    const double corner[4][2] = {
        {tile->uv_min[0], tile->uv_min[1]},
        {tile->uv_max[0], tile->uv_min[1]},
        {tile->uv_max[0], tile->uv_max[1]},
        {tile->uv_min[0], tile->uv_max[1]}
    };
    long double value[4];
    const long double a = coefficient[0][0], b = coefficient[0][1];
    const long double scale = 1.0L +
        fmaxl(fabsl(bound[0]), fabsl(bound[1])) +
        fabsl(a) * fmaxl(fabsl(tile->uv_min[0]),
                         fabsl(tile->uv_max[0])) +
        fabsl(b) * fmaxl(fabsl(tile->uv_min[1]),
                         fabsl(tile->uv_max[1]));
    const long double numerical_tolerance = scale *
        (128.0L * LDBL_EPSILON + 16.0L * DBL_EPSILON);
    const long double tolerance = numerical_tolerance;
    const long double separation_tolerance = 8.0L * numerical_tolerance;
    if (!isfinite(tolerance) || !(bound[1] - bound[0] >
                                  separation_tolerance)) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
        return 0;
    }
    for (int i = 0; i < 4; ++i) {
        value[i] = a * corner[i][0] + b * corner[i][1];
        if (!isfinite(value[i]) ||
            fabsl(value[i] - bound[0]) <= tolerance ||
            fabsl(value[i] - bound[1]) <= tolerance) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
    }
    slice_error_polygon_t polygon[3] = {{0}, {0}, {0}};
    slice_error_polygon_t crossings[2] = {{0}, {0}};
    for (int i = 0; i < 4; ++i) {
        const int next = (i + 1) & 3;
        const int band = value[i] < bound[0] ? 0
            : value[i] < bound[1] ? 1 : 2;
        if (!slice_error_polygon_add(&polygon[band],
                corner[i][0], corner[i][1], 0.0)) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
            return 0;
        }
        const int increasing = value[next] > value[i];
        for (int step = 0; step < 2; ++step) {
            const int k = increasing ? step : 1 - step;
            if ((value[i] < bound[k]) == (value[next] < bound[k]))
                continue;
            const long double t =
                (bound[k] - value[i]) / (value[next] - value[i]);
            const long double du_edge =
                (long double)corner[next][0] - corner[i][0];
            const long double dv_edge =
                (long double)corner[next][1] - corner[i][1];
            const long double u = (long double)corner[i][0] + t * du_edge;
            const long double v = (long double)corner[i][1] + t * dv_edge;
            const double du = (double)u, dv = (double)v;
            const long double edge_length =
                fmaxl(fabsl(du_edge), fabsl(dv_edge));
            const double uncertainty = (double)(
                4.0L * tolerance / fabsl(value[next] - value[i]) *
                    edge_length +
                8.0L * DBL_EPSILON *
                    fmaxl(1.0L, fmaxl(fabsl(u), fabsl(v))));
            const long double corner_distance =
                fminl(t, 1.0L - t) * edge_length;
            if (!(t > 0.0L && t < 1.0L) || !isfinite(du) ||
                !isfinite(dv) || !isfinite(uncertainty) ||
                !((long double)uncertainty * 4.0L < corner_distance) ||
                !slice_error_polygon_add(&crossings[k], du, dv,
                                         uncertainty) ||
                !slice_error_polygon_add(&polygon[k], du, dv,
                                         uncertainty) ||
                !slice_error_polygon_add(&polygon[k + 1], du, dv,
                                         uncertainty)) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                return 0;
            }
        }
    }
    for (int k = 0; k < 2; ++k)
        if (crossings[k].count != 0 && crossings[k].count != 2) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
    slice_error_face_t face[3] = {{0}, {0}, {0}};
    for (int band = 0; band < 3; ++band) {
        if (!polygon[band].count) continue;
        if (polygon[band].count < 3) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
        double u = 0.0, v = 0.0;
        for (size_t j = 0; j < polygon[band].count; ++j) {
            u += polygon[band].uv[j][0] / polygon[band].count;
            v += polygon[band].uv[j][1] / polygon[band].count;
        }
        const long double s = a * (long double)u + b * (long double)v;
        if (!isfinite(u) || !isfinite(v) ||
            (band > 0 && !(s - bound[band - 1] > tolerance)) ||
            (band < 2 && !(bound[band] - s > tolerance))) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
        double world[3];
        if (!slice_error_world_point(frame, u, v, world)) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
        int raw_negative[2];
        for (int i = 0; i < 2; ++i)
            raw_negative[i] = orientation[i] > 0
                ? s < threshold[i] : s > threshold[i];
        for (size_t ci = 0; ci < cell_count; ++ci) {
            const alea_cell_entry_t* cell = &sys->cells.data[cells[ci]];
            int inside = 0;
            if (!slice_error_oblique_node_inside(sys,
                    cell->root_node_id, primitive_ids, raw_negative,
                    2, constant_planes, constant_count,
                    0, work, &inside)) {
                *reason = *work == 0
                    ? ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT
                    : ALEA_SLICE_ERROR_UNRESOLVED_UNSUPPORTED_GEOMETRY;
                return 0;
            }
            if (inside != alea_point_inside(sys, cell->root_node_id,
                    world[0], world[1], world[2])) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                return 0;
            }
            if (inside) {
                if (face[band].owner_count ==
                    ALEA_SLICE_ERROR_OWNER_CAPACITY) {
                    *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
                    return 0;
                }
                face[band].owner_cell_ids[face[band].owner_count++] =
                    cell->mc_cell_id;
            }
        }
        face[band].kind = face[band].owner_count == 0
            ? ALEA_POINT_COVERAGE_GAP
            : face[band].owner_count == 1
                ? ALEA_POINT_COVERAGE_UNIQUE
                : ALEA_POINT_COVERAGE_OVERLAP;
    }
    slice_error_oblique_segment_t segments[2] = {{0}, {0}};
    size_t segment_count = 0;
    for (int k = 0; k < 2; ++k) {
        if (crossings[k].count != 2) continue;
        const int plane = line[k];
        slice_error_oblique_segment_t* segment =
            &segments[segment_count++];
        segment->plane = plane;
        segment->negative_face = orientation[plane] > 0 ? k : k + 1;
        segment->positive_face = orientation[plane] > 0 ? k + 1 : k;
        memcpy(segment->uv, crossings[k].uv, sizeof(segment->uv));
        segment->uncertainty[0] = crossings[k].uncertainty[0][0];
        segment->uncertainty[1] = crossings[k].uncertainty[1][0];
    }
    return slice_error_publish_oblique_faces(
        o, selected, polygon, face, 3, segments, segment_count,
        contextual_bytes, out, reason, output_omitted);
}

/* Two nonparallel slice lines partition a rectangular core into at most four convex
 * faces. An interior crossing is a vertex of all four faces. A crossing well
 * outside the core has no face vertex. Boundary or poorly conditioned
 * crossings are withheld. */
static int slice_error_classify_crossing_oblique_tile(
    const alea_slice_error_query_t* query,
    const alea_transition_slice_critical_tile_t* tile,
    const size_t* cells, size_t cell_count,
    const slice_error_axis_plane_t selected[2],
    const slice_error_axis_plane_t* constant_planes,
    size_t constant_count,
    size_t contextual_bytes, size_t* work,
    alea_slice_error_page_t* out,
    alea_slice_error_unresolved_reason_t* reason,
    int* output_omitted) {
    const alea_slice_error_query_options_t* o = &query->options;
    const alea_system_t* sys = query->sys;
    const alea_slice_plane_t* frame = &o->view.plane;
    uint32_t primitive_ids[2];
    long double a[2], b[2], c[2], tolerance[2];
    const double corner[4][2] = {
        {tile->uv_min[0], tile->uv_min[1]},
        {tile->uv_max[0], tile->uv_min[1]},
        {tile->uv_max[0], tile->uv_max[1]},
        {tile->uv_min[0], tile->uv_max[1]}
    };
    long double value[2][4];
    for (int i = 0; i < 2; ++i) {
        if (selected[i].surface_id <= 0) return 0;
        primitive_ids[i] = selected[i].primitive_id;
        const alea_primitive_entry_t* primitive =
            &sys->primitives.data[primitive_ids[i]];
        if (primitive->type != ALEA_PRIMITIVE_PLANE ||
            primitive->payload_index >= sys->primitive_planes.count)
            return 0;
        const alea_plane_data_t* p =
            &sys->primitive_planes.data[primitive->payload_index];
        if (!isfinite(p->a) || !isfinite(p->b) || !isfinite(p->c) ||
            !isfinite(p->d)) return 0;
        slice_error_project_plane(p, frame, &a[i], &b[i], &c[i]);
        const long double norm = fmaxl(fabsl(a[i]), fabsl(b[i]));
        if (!(norm > 0.0L) || !isfinite(norm)) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
        a[i] /= norm; b[i] /= norm; c[i] /= norm;
        const long double scale = 1.0L + fabsl(c[i]) +
            fabsl(a[i]) * fmaxl(fabsl(tile->uv_min[0]),
                                fabsl(tile->uv_max[0])) +
            fabsl(b[i]) * fmaxl(fabsl(tile->uv_min[1]),
                                fabsl(tile->uv_max[1]));
        tolerance[i] = scale *
            (128.0L * LDBL_EPSILON + 16.0L * DBL_EPSILON);
        if (!isfinite(tolerance[i])) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
        for (int j = 0; j < 4; ++j) {
            value[i][j] = a[i] * corner[j][0] +
                          b[i] * corner[j][1] + c[i];
            if (!isfinite(value[i][j]) ||
                fabsl(value[i][j]) <= tolerance[i]) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                return 0;
            }
        }
    }
    const long double det = a[0] * b[1] - a[1] * b[0];
    if (!isfinite(det) || fabsl(det) <=
        128.0L * LDBL_EPSILON) return 0;
    const long double center_u = (b[0] * c[1] - b[1] * c[0]) / det;
    const long double center_v = (a[1] * c[0] - a[0] * c[1]) / det;
    const double crossing[2] = {(double)center_u, (double)center_v};
    const double center_uncertainty = (double)(
        32.0L * (tolerance[0] + tolerance[1]) / fabsl(det) +
        8.0L * DBL_EPSILON *
            fmaxl(1.0L, fmaxl(fabsl(center_u), fabsl(center_v))));
    if (!isfinite(center_u) || !isfinite(center_v) ||
        !isfinite(crossing[0]) || !isfinite(crossing[1]) ||
        !isfinite(center_uncertainty)) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
        return 0;
    }
    int center_inside = 1, center_outside = 0;
    for (int axis = 0; axis < 2; ++axis) {
        if (!(crossing[axis] - tile->uv_min[axis] >
              4.0 * center_uncertainty &&
              tile->uv_max[axis] - crossing[axis] >
              4.0 * center_uncertainty)) center_inside = 0;
        if (crossing[axis] < tile->uv_min[axis] -
                4.0 * center_uncertainty ||
            crossing[axis] > tile->uv_max[axis] +
                4.0 * center_uncertainty) center_outside = 1;
    }
    if (!center_inside && !center_outside) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
        return 0;
    }
    slice_error_polygon_t polygon[4] = {{0}, {0}, {0}, {0}};
    slice_error_polygon_t edge_hit[2] = {{0}, {0}};
    for (int j = 0; j < 4; ++j) {
        const int mask = (value[0][j] > 0.0L ? 1 : 0) |
                         (value[1][j] > 0.0L ? 2 : 0);
        if (!slice_error_polygon_add(&polygon[mask],
                corner[j][0], corner[j][1], 0.0)) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
            return 0;
        }
        const int next = (j + 1) & 3;
        for (int i = 0; i < 2; ++i) {
            if ((value[i][j] > 0.0L) == (value[i][next] > 0.0L))
                continue;
            const long double t = value[i][j] /
                (value[i][j] - value[i][next]);
            const long double du_edge =
                (long double)corner[next][0] - corner[j][0];
            const long double dv_edge =
                (long double)corner[next][1] - corner[j][1];
            const long double u = corner[j][0] + t * du_edge;
            const long double v = corner[j][1] + t * dv_edge;
            const double du = (double)u, dv = (double)v;
            const long double edge_length =
                fmaxl(fabsl(du_edge), fabsl(dv_edge));
            const double uncertainty = (double)(
                4.0L * tolerance[i] /
                    fabsl(value[i][j] - value[i][next]) * edge_length +
                8.0L * DBL_EPSILON *
                    fmaxl(1.0L, fmaxl(fabsl(u), fabsl(v))));
            const int other = 1 - i;
            const long double other_value =
                a[other] * u + b[other] * v + c[other];
            if (!(t > 0.0L && t < 1.0L) || !isfinite(du) ||
                !isfinite(dv) || !isfinite(uncertainty) ||
                !((long double)uncertainty * 4.0L <
                  fminl(t, 1.0L - t) * edge_length) ||
                !(fabsl(other_value) > 4.0L *
                    (tolerance[other] +
                     (fabsl(a[other]) + fabsl(b[other])) *
                        uncertainty))) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                return 0;
            }
            const int other_bit = other_value > 0.0L ? 1 : 0;
            const int negative_face = i == 0 ? other_bit << 1 : other_bit;
            const int positive_face = negative_face | (1 << i);
            if (!slice_error_polygon_add(&edge_hit[i], du, dv,
                                         uncertainty) ||
                !slice_error_polygon_add(&polygon[negative_face], du, dv,
                                         uncertainty) ||
                !slice_error_polygon_add(&polygon[positive_face], du, dv,
                                         uncertainty)) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
                return 0;
            }
        }
    }
    for (int i = 0; i < 2; ++i) {
        if (edge_hit[i].count != 0 && edge_hit[i].count != 2) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
        if (center_inside && edge_hit[i].count != 2) return 0;
        if (center_inside)
            for (int j = 0; j < 2; ++j) {
                const double separation = fmax(
                    fabs(crossing[0] - edge_hit[i].uv[j][0]),
                    fabs(crossing[1] - edge_hit[i].uv[j][1]));
                if (!(separation > 4.0 * (center_uncertainty +
                            edge_hit[i].uncertainty[j][0]))) {
                    *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                    return 0;
                }
            }
    }
    slice_error_face_t face[4] = {{0}, {0}, {0}, {0}};
    for (int mask = 0; mask < 4; ++mask) {
        if (center_inside && !slice_error_polygon_add(&polygon[mask],
                crossing[0], crossing[1], center_uncertainty)) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
        if (!polygon[mask].count) continue;
        if (polygon[mask].count < 3) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
        double centroid_u = 0.0, centroid_v = 0.0;
        for (size_t j = 0; j < polygon[mask].count; ++j) {
            centroid_u += polygon[mask].uv[j][0] / polygon[mask].count;
            centroid_v += polygon[mask].uv[j][1] / polygon[mask].count;
        }
        /* Every point belongs to this convex face. Angular order gives its
         * boundary, including a face with no rectangle corner. */
        for (size_t j = 1; j < polygon[mask].count; ++j) {
            const double u = polygon[mask].uv[j][0];
            const double v = polygon[mask].uv[j][1];
            const double uncertainty = polygon[mask].uncertainty[j][0];
            const double angle = atan2(v - centroid_v, u - centroid_u);
            size_t k = j;
            while (k > 0 && atan2(
                    polygon[mask].uv[k - 1][1] - centroid_v,
                    polygon[mask].uv[k - 1][0] - centroid_u) > angle) {
                memcpy(polygon[mask].uv[k], polygon[mask].uv[k - 1],
                       sizeof(polygon[mask].uv[k]));
                memcpy(polygon[mask].uncertainty[k],
                       polygon[mask].uncertainty[k - 1],
                       sizeof(polygon[mask].uncertainty[k]));
                --k;
            }
            polygon[mask].uv[k][0] = u;
            polygon[mask].uv[k][1] = v;
            polygon[mask].uncertainty[k][0] = uncertainty;
            polygon[mask].uncertainty[k][1] = uncertainty;
        }
        for (size_t j = 0; j < polygon[mask].count; ++j)
            for (size_t k = j + 1; k < polygon[mask].count; ++k) {
                const double separation = fmax(
                    fabs(polygon[mask].uv[j][0] - polygon[mask].uv[k][0]),
                    fabs(polygon[mask].uv[j][1] - polygon[mask].uv[k][1]));
                if (!(separation > 4.0 *
                    (polygon[mask].uncertainty[j][0] +
                     polygon[mask].uncertainty[k][0]))) {
                    *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                    return 0;
                }
            }
        const long double residual[2] = {
            a[0] * centroid_u + b[0] * centroid_v + c[0],
            a[1] * centroid_u + b[1] * centroid_v + c[1]
        };
        for (int i = 0; i < 2; ++i)
            if (((mask & (1 << i)) != 0)
                    ? !(residual[i] > tolerance[i])
                    : !(residual[i] < -tolerance[i])) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                return 0;
            }
        double world[3];
        if (!slice_error_world_point(frame, centroid_u, centroid_v, world)) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
        const int raw_negative[2] = {
            (mask & 1) == 0, (mask & 2) == 0
        };
        for (size_t ci = 0; ci < cell_count; ++ci) {
            const alea_cell_entry_t* cell = &sys->cells.data[cells[ci]];
            int inside = 0;
            if (!slice_error_oblique_node_inside(sys,
                    cell->root_node_id, primitive_ids, raw_negative,
                    2, constant_planes, constant_count,
                    0, work, &inside)) {
                *reason = *work == 0
                    ? ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT
                    : ALEA_SLICE_ERROR_UNRESOLVED_UNSUPPORTED_GEOMETRY;
                return 0;
            }
            if (inside != alea_point_inside(sys, cell->root_node_id,
                    world[0], world[1], world[2])) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                return 0;
            }
            if (inside) {
                if (face[mask].owner_count ==
                    ALEA_SLICE_ERROR_OWNER_CAPACITY) {
                    *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
                    return 0;
                }
                face[mask].owner_cell_ids[face[mask].owner_count++] =
                    cell->mc_cell_id;
            }
        }
        face[mask].kind = face[mask].owner_count == 0
            ? ALEA_POINT_COVERAGE_GAP
            : face[mask].owner_count == 1
                ? ALEA_POINT_COVERAGE_UNIQUE
                : ALEA_POINT_COVERAGE_OVERLAP;
    }
    slice_error_oblique_segment_t segments[4] = {{0}, {0}, {0}, {0}};
    size_t segment_count = 0;
    for (int i = 0; i < 2; ++i) {
        const int other = 1 - i;
        const int parts = center_inside ? 2
            : edge_hit[i].count == 2 ? 1 : 0;
        for (int j = 0; j < parts; ++j) {
            const double* start = edge_hit[i].uv[j];
            const double* end = center_inside
                ? crossing : edge_hit[i].uv[1];
            const long double mid_u =
                ((long double)start[0] + end[0]) * 0.5L;
            const long double mid_v =
                ((long double)start[1] + end[1]) * 0.5L;
            const long double other_value =
                a[other] * mid_u + b[other] * mid_v + c[other];
            if (!(fabsl(other_value) > 4.0L * tolerance[other])) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                return 0;
            }
            const int other_bit = other_value > 0.0L ? 1 : 0;
            slice_error_oblique_segment_t* segment =
                &segments[segment_count++];
            segment->plane = i;
            segment->negative_face = i == 0
                ? other_bit << 1 : other_bit;
            segment->positive_face =
                segment->negative_face | (1 << i);
            memcpy(segment->uv[0], start,
                   sizeof(segment->uv[0]));
            memcpy(segment->uv[1], end,
                   sizeof(segment->uv[1]));
            segment->uncertainty[0] =
                edge_hit[i].uncertainty[j][0];
            segment->uncertainty[1] = center_inside
                ? center_uncertainty : edge_hit[i].uncertainty[1][0];
        }
    }
    return slice_error_publish_oblique_faces(
        o, selected, polygon, face, 4, segments, segment_count,
        contextual_bytes, out, reason, output_omitted);
}

typedef struct {
    long double uv[2];
    long double uncertainty;
} slice_error_arr_vertex_t;

typedef struct {
    long double a, b, c, tolerance;
    uint32_t primitive_id;
    int surface_id;
} slice_error_arr_line_t;

typedef struct {
    slice_error_arr_vertex_t* vertex;
    size_t count;
    uint64_t positive_mask;
    slice_error_face_t ownership;
} slice_error_arr_face_t;

typedef struct {
    size_t face, neighbor, line, vertex;
} slice_error_arr_edge_t;

static long double slice_error_arr_value(
    const slice_error_arr_line_t* line,
    const slice_error_arr_vertex_t* vertex) {
    return line->a * vertex->uv[0] +
           line->b * vertex->uv[1] + line->c;
}

static long double slice_error_arr_tolerance(
    const slice_error_arr_line_t* line,
    const slice_error_arr_vertex_t* vertex) {
    return line->tolerance +
        (fabsl(line->a) + fabsl(line->b)) * vertex->uncertainty;
}

/* Input vertices are strictly separated from the clipping line. */
static int slice_error_arr_clip(
    const slice_error_arr_vertex_t* input, size_t count,
    const slice_error_arr_line_t* line, int keep_positive,
    slice_error_arr_vertex_t* output, size_t capacity,
    size_t* output_count, size_t* work) {
    *output_count = 0;
    for (size_t i = 0; i < count; ++i) {
        if (!*work) return 0;
        --*work;
        const slice_error_arr_vertex_t* first = &input[i];
        const slice_error_arr_vertex_t* second = &input[(i + 1) % count];
        const long double f0 = slice_error_arr_value(line, first);
        const long double f1 = slice_error_arr_value(line, second);
        const int inside = keep_positive ? f0 > 0.0L : f0 < 0.0L;
        const int next_inside = keep_positive ? f1 > 0.0L : f1 < 0.0L;
        if (inside) {
            if (*output_count == capacity) return 0;
            output[(*output_count)++] = *first;
        }
        if (inside == next_inside) continue;
        const long double fraction = f0 / (f0 - f1);
        const long double length = fmaxl(
            fabsl(second->uv[0] - first->uv[0]),
            fabsl(second->uv[1] - first->uv[1]));
        const long double error =
            4.0L * (slice_error_arr_tolerance(line, first) +
                    slice_error_arr_tolerance(line, second)) /
                fabsl(f0 - f1) * length +
            first->uncertainty + second->uncertainty;
        if (!(fraction > 0.0L && fraction < 1.0L) ||
            !isfinite(error) ||
            !(4.0L * error < fminl(fraction, 1.0L - fraction) * length) ||
            *output_count == capacity) return 0;
        slice_error_arr_vertex_t* crossing = &output[(*output_count)++];
        for (int axis = 0; axis < 2; ++axis)
            crossing->uv[axis] = first->uv[axis] + fraction *
                (second->uv[axis] - first->uv[axis]);
        crossing->uncertainty = error +
            16.0L * DBL_EPSILON * fmaxl(1.0L,
                fmaxl(fabsl(crossing->uv[0]),
                      fabsl(crossing->uv[1])));
        if (!isfinite(crossing->uv[0]) ||
            !isfinite(crossing->uv[1]) ||
            !isfinite(crossing->uncertainty)) return 0;
    }
    return *output_count >= 3;
}

static int slice_error_arr_tile_edge(
    const slice_error_arr_vertex_t* first,
    const slice_error_arr_vertex_t* second,
    const alea_transition_slice_critical_tile_t* tile) {
    for (int axis = 0; axis < 2; ++axis)
        if ((first->uv[axis] == tile->uv_min[axis] &&
             second->uv[axis] == tile->uv_min[axis]) ||
            (first->uv[axis] == tile->uv_max[axis] &&
             second->uv[axis] == tile->uv_max[axis])) return 1;
    return 0;
}

static int slice_error_arr_same_vertex(
    const slice_error_arr_vertex_t* first,
    const slice_error_arr_vertex_t* second) {
    const long double separation = fmaxl(
        fabsl(first->uv[0] - second->uv[0]),
        fabsl(first->uv[1] - second->uv[1]));
    const long double tolerance = 4.0L *
        (first->uncertainty + second->uncertainty) +
        16.0L * DBL_EPSILON * fmaxl(1.0L,
            fmaxl(fabsl(first->uv[0]), fabsl(first->uv[1])));
    return separation <= tolerance;
}

/* Split every convex face by every projected plane. The line and face counts
 * are limited by the caller's scratch and work budgets before publication. */
static int slice_error_classify_line_arrangement(
    const alea_slice_error_query_t* query,
    const alea_transition_slice_critical_tile_t* tile,
    const size_t* cells, size_t cell_count,
    const slice_error_axis_plane_t* planes, size_t primitive_count,
    size_t contextual_bytes, size_t* work,
    alea_slice_error_page_t* out,
    alea_slice_error_unresolved_reason_t* reason,
    int* output_omitted, size_t* peak_scratch) {
    const alea_system_t* sys = query->sys;
    const alea_slice_error_query_options_t* options = &query->options;
    const alea_slice_plane_t* frame = &options->view.plane;
    slice_error_arr_line_t line[63] = {{0}};
    uint32_t primitive_ids[63];
    size_t line_count = 0;
    *reason = ALEA_SLICE_ERROR_UNRESOLVED_PLANAR_ARRANGEMENT;
    for (size_t pi = 0; pi < primitive_count; ++pi) {
        if (planes[pi].axis == 2) continue;
        if (line_count == 63) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
            return 0;
        }
        const alea_primitive_entry_t* primitive =
            &sys->primitives.data[planes[pi].primitive_id];
        const alea_plane_data_t* p =
            &sys->primitive_planes.data[primitive->payload_index];
        slice_error_arr_line_t* current = &line[line_count];
        slice_error_project_plane(p, frame,
            &current->a, &current->b, &current->c);
        const long double norm = fmaxl(fabsl(current->a),
                                       fabsl(current->b));
        if (!(norm > 0.0L) || !isfinite(norm)) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
        current->a /= norm;
        current->b /= norm;
        current->c /= norm;
        const long double scale = 1.0L + fabsl(current->c) +
            fabsl(current->a) * fmaxl(fabsl(tile->uv_min[0]),
                                      fabsl(tile->uv_max[0])) +
            fabsl(current->b) * fmaxl(fabsl(tile->uv_min[1]),
                                      fabsl(tile->uv_max[1]));
        current->tolerance = scale *
            (128.0L * LDBL_EPSILON + 16.0L * DBL_EPSILON);
        if (!isfinite(current->c) ||
            !isfinite(current->tolerance)) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            return 0;
        }
        current->primitive_id = planes[pi].primitive_id;
        current->surface_id = planes[pi].surface_id;
        primitive_ids[line_count++] = planes[pi].primitive_id;
    }
    if (line_count < 3) return 0;
    const size_t max_faces = 1 + line_count * (line_count + 1) / 2;
    const size_t max_vertices = line_count + 4;
    const size_t edge_capacity = max_faces * max_vertices;
    const size_t arrangement_bytes =
        max_faces * sizeof(slice_error_arr_face_t) +
        edge_capacity * (sizeof(slice_error_arr_vertex_t) +
                         sizeof(slice_error_arr_edge_t)) +
        2 * max_vertices * sizeof(slice_error_arr_vertex_t) +
        sizeof(line) + sizeof(primitive_ids);
    const size_t scratch_limit =
        options->scan_options.max_critical_scratch_bytes;
    if (*peak_scratch > scratch_limit ||
        arrangement_bytes > scratch_limit - *peak_scratch ||
        line_count > *work / (max_faces * max_vertices)) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
        return 0;
    }
    *peak_scratch += arrangement_bytes;
    slice_error_arr_face_t* face = calloc(max_faces, sizeof(*face));
    slice_error_arr_vertex_t* vertex = calloc(edge_capacity,
                                               sizeof(*vertex));
    slice_error_arr_vertex_t* negative = calloc(max_vertices,
                                                 sizeof(*negative));
    slice_error_arr_vertex_t* positive = calloc(max_vertices,
                                                 sizeof(*positive));
    slice_error_arr_edge_t* edge = calloc(edge_capacity, sizeof(*edge));
    if (!face || !vertex || !negative || !positive || !edge) {
        free(face); free(vertex); free(negative); free(positive); free(edge);
        return -1;
    }
    int status = 0;
    size_t face_count = 1, edge_count = 0;
    for (size_t i = 0; i < max_faces; ++i)
        face[i].vertex = vertex + i * max_vertices;
    face[0].count = 4;
    const double corner[4][2] = {
        {tile->uv_min[0], tile->uv_min[1]},
        {tile->uv_max[0], tile->uv_min[1]},
        {tile->uv_max[0], tile->uv_max[1]},
        {tile->uv_min[0], tile->uv_max[1]}
    };
    for (int i = 0; i < 4; ++i) {
        face[0].vertex[i].uv[0] = corner[i][0];
        face[0].vertex[i].uv[1] = corner[i][1];
    }
    for (size_t li = 0; li < line_count; ++li) {
        const size_t existing_faces = face_count;
        for (size_t fi = 0; fi < existing_faces; ++fi) {
            slice_error_arr_face_t* current = &face[fi];
            int seen_negative = 0, seen_positive = 0;
            for (size_t vi = 0; vi < current->count; ++vi) {
                if (!*work) {
                    *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
                    goto done;
                }
                --*work;
                const long double value = slice_error_arr_value(
                    &line[li], &current->vertex[vi]);
                if (!isfinite(value) ||
                    !(fabsl(value) > slice_error_arr_tolerance(
                        &line[li], &current->vertex[vi]))) {
                    *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                    goto done;
                }
                seen_negative |= value < 0.0L;
                seen_positive |= value > 0.0L;
            }
            if (!seen_negative) {
                current->positive_mask |= UINT64_C(1) << li;
                continue;
            }
            if (!seen_positive) continue;
            if (face_count == max_faces) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
                goto done;
            }
            size_t negative_count, positive_count;
            if (!slice_error_arr_clip(current->vertex, current->count,
                    &line[li], 0, negative, max_vertices,
                    &negative_count, work) ||
                !slice_error_arr_clip(current->vertex, current->count,
                    &line[li], 1, positive, max_vertices,
                    &positive_count, work)) {
                *reason = *work
                    ? ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL
                    : ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
                goto done;
            }
            slice_error_arr_face_t* added = &face[face_count++];
            added->count = positive_count;
            added->positive_mask =
                current->positive_mask | (UINT64_C(1) << li);
            memcpy(added->vertex, positive,
                   positive_count * sizeof(*positive));
            current->count = negative_count;
            memcpy(current->vertex, negative,
                   negative_count * sizeof(*negative));
        }
    }
    for (size_t fi = 0; fi < face_count; ++fi) {
        slice_error_arr_face_t* current = &face[fi];
        long double u = 0.0L, v = 0.0L;
        for (size_t vi = 0; vi < current->count; ++vi) {
            u += current->vertex[vi].uv[0] / current->count;
            v += current->vertex[vi].uv[1] / current->count;
            if (!isfinite(current->vertex[vi].uv[0]) ||
                !isfinite(current->vertex[vi].uv[1]) ||
                !isfinite(current->vertex[vi].uncertainty) ||
                fabsl(current->vertex[vi].uv[0]) > DBL_MAX ||
                fabsl(current->vertex[vi].uv[1]) > DBL_MAX ||
                current->vertex[vi].uncertainty +
                    32.0L * DBL_EPSILON * fmaxl(1.0L,
                        fmaxl(fabsl(current->vertex[vi].uv[0]),
                              fabsl(current->vertex[vi].uv[1]))) >
                    DBL_MAX) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                goto done;
            }
        }
        if (!isfinite(u) || !isfinite(v) ||
            fabsl(u) > DBL_MAX || fabsl(v) > DBL_MAX) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            goto done;
        }
        int raw_negative[63];
        for (size_t li = 0; li < line_count; ++li) {
            const long double residual = line[li].a * u +
                line[li].b * v + line[li].c;
            const int positive_side =
                (current->positive_mask & (UINT64_C(1) << li)) != 0;
            if ((positive_side &&
                 !(residual > 4.0L * line[li].tolerance)) ||
                (!positive_side &&
                 !(residual < -4.0L * line[li].tolerance))) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                goto done;
            }
            raw_negative[li] = !positive_side;
        }
        double world[3];
        if (!slice_error_world_point(frame, (double)u, (double)v, world)) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            goto done;
        }
        for (size_t ci = 0; ci < cell_count; ++ci) {
            const alea_cell_entry_t* cell = &sys->cells.data[cells[ci]];
            int inside = 0;
            if (!slice_error_oblique_node_inside(sys,
                    cell->root_node_id, primitive_ids, raw_negative,
                    line_count, planes, primitive_count, 0, work, &inside)) {
                *reason = *work
                    ? ALEA_SLICE_ERROR_UNRESOLVED_PLANAR_ARRANGEMENT
                    : ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
                goto done;
            }
            if (inside != alea_point_inside(sys, cell->root_node_id,
                    world[0], world[1], world[2])) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                goto done;
            }
            if (!inside) continue;
            slice_error_face_t* ownership = &current->ownership;
            if (ownership->owner_count == ALEA_SLICE_ERROR_OWNER_CAPACITY) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
                goto done;
            }
            ownership->owner_cell_ids[ownership->owner_count++] =
                cell->mc_cell_id;
        }
        current->ownership.kind = current->ownership.owner_count == 0
            ? ALEA_POINT_COVERAGE_GAP
            : current->ownership.owner_count == 1
                ? ALEA_POINT_COVERAGE_UNIQUE
                : ALEA_POINT_COVERAGE_OVERLAP;
    }
    size_t region_count = 0;
    for (size_t fi = 0; fi < face_count; ++fi)
        if (slice_error_axis_is_defect(face[fi].ownership.kind)) {
            const size_t pieces = face[fi].count <=
                ALEA_SLICE_ERROR_POLYGON_CAPACITY
                ? 1 : face[fi].count - 2;
            if (region_count > SIZE_MAX - pieces) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
                goto done;
            }
            region_count += pieces;
        }
    for (size_t fi = 0; fi < face_count; ++fi) {
        const slice_error_arr_face_t* current = &face[fi];
        for (size_t vi = 0; vi < current->count; ++vi) {
            if (!*work) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
                goto done;
            }
            --*work;
            const slice_error_arr_vertex_t* first = &current->vertex[vi];
            const slice_error_arr_vertex_t* second =
                &current->vertex[(vi + 1) % current->count];
            if (slice_error_arr_tile_edge(first, second, tile)) continue;
            size_t matching_line = SIZE_MAX;
            for (size_t li = 0; li < line_count; ++li) {
                if (!*work) {
                    *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
                    goto done;
                }
                --*work;
                if (fabsl(slice_error_arr_value(&line[li], first)) >
                        4.0L * slice_error_arr_tolerance(&line[li], first) ||
                    fabsl(slice_error_arr_value(&line[li], second)) >
                        4.0L * slice_error_arr_tolerance(&line[li], second))
                    continue;
                if (matching_line != SIZE_MAX) {
                    *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                    goto done;
                }
                matching_line = li;
            }
            if (matching_line == SIZE_MAX) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                goto done;
            }
            /* Publish each interface from its negative-sign face once. */
            if (current->positive_mask & (UINT64_C(1) << matching_line))
                continue;
            const uint64_t neighbor_mask = current->positive_mask ^
                (UINT64_C(1) << matching_line);
            size_t neighbor = 0;
            while (neighbor < face_count &&
                   face[neighbor].positive_mask != neighbor_mask) {
                if (!*work) {
                    *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
                    goto done;
                }
                --*work;
                ++neighbor;
            }
            if (neighbor == face_count) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                goto done;
            }
            int shared = 0;
            for (size_t other = 0; other < face[neighbor].count; ++other) {
                const slice_error_arr_vertex_t* a =
                    &face[neighbor].vertex[other];
                const slice_error_arr_vertex_t* b =
                    &face[neighbor].vertex[(other + 1) %
                        face[neighbor].count];
                if (slice_error_arr_same_vertex(first, b) &&
                    slice_error_arr_same_vertex(second, a)) {
                    shared = 1;
                    break;
                }
            }
            if (!shared) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                goto done;
            }
            if (slice_error_axis_is_defect(current->ownership.kind) ==
                slice_error_axis_is_defect(face[neighbor].ownership.kind))
                continue;
            if (edge_count == edge_capacity) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
                goto done;
            }
            edge[edge_count++] = (slice_error_arr_edge_t){
                fi, neighbor, matching_line, vi
            };
        }
    }
    const size_t output_limit = options->scan_options.max_output_bytes;
    if (contextual_bytes > output_limit ||
        region_count > (output_limit - contextual_bytes) /
            sizeof(alea_slice_error_region_t) ||
        edge_count > (output_limit - contextual_bytes -
            region_count * sizeof(alea_slice_error_region_t)) /
            sizeof(alea_slice_error_interval_t)) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
        *output_omitted = 1;
        goto done;
    }
    out->regions = region_count
        ? calloc(region_count, sizeof(*out->regions)) : NULL;
    out->intervals = edge_count
        ? calloc(edge_count, sizeof(*out->intervals)) : NULL;
    if ((region_count && !out->regions) ||
        (edge_count && !out->intervals)) {
        status = -1;
        goto done;
    }
    for (size_t fi = 0; fi < face_count; ++fi) {
        const slice_error_arr_face_t* current = &face[fi];
        if (!slice_error_axis_is_defect(current->ownership.kind)) continue;
        const size_t pieces = current->count <=
            ALEA_SLICE_ERROR_POLYGON_CAPACITY
            ? 1 : current->count - 2;
        for (size_t piece = 0; piece < pieces; ++piece) {
            alea_slice_error_region_t* region =
                &out->regions[out->region_count++];
            region->kind = current->ownership.kind;
            region->owner_count = current->ownership.owner_count;
            memcpy(region->owner_cell_ids,
                current->ownership.owner_cell_ids,
                region->owner_count * sizeof(int));
            region->polygon_vertex_count = pieces == 1
                ? current->count : 3;
            region->uv_min[0] = region->uv_min[1] = INFINITY;
            region->uv_max[0] = region->uv_max[1] = -INFINITY;
            for (size_t vi = 0; vi < region->polygon_vertex_count; ++vi) {
                const size_t source = pieces == 1 ? vi
                    : vi == 0 ? 0 : piece + vi;
                const slice_error_arr_vertex_t* point =
                    &current->vertex[source];
                for (int axis = 0; axis < 2; ++axis) {
                    const double value = (double)point->uv[axis];
                    const double uncertainty = (double)(
                        point->uncertainty +
                        fabsl(point->uv[axis] - (long double)value));
                    region->polygon_uv[vi][axis] = value;
                    region->polygon_uv_uncertainty[vi][axis] =
                        uncertainty;
                    region->uv_min[axis] =
                        fmin(region->uv_min[axis], value);
                    region->uv_max[axis] =
                        fmax(region->uv_max[axis], value);
                    region->uv_min_uncertainty[axis] = fmax(
                        region->uv_min_uncertainty[axis], uncertainty);
                    region->uv_max_uncertainty[axis] = fmax(
                        region->uv_max_uncertainty[axis], uncertainty);
                }
            }
        }
    }
    for (size_t ei = 0; ei < edge_count; ++ei) {
        const slice_error_arr_edge_t* selected = &edge[ei];
        const slice_error_arr_face_t* negative_face =
            &face[selected->face];
        const slice_error_arr_face_t* positive_face =
            &face[selected->neighbor];
        const slice_error_arr_vertex_t* endpoints[2] = {
            &negative_face->vertex[selected->vertex],
            &negative_face->vertex[(selected->vertex + 1) %
                                   negative_face->count]
        };
        const int reverse = endpoints[0]->uv[0] > endpoints[1]->uv[0] ||
            (endpoints[0]->uv[0] == endpoints[1]->uv[0] &&
             endpoints[0]->uv[1] > endpoints[1]->uv[1]);
        alea_slice_error_interval_t* interval =
            &out->intervals[out->interval_count++];
        interval->evidence_scope =
            ALEA_SLICE_BOUNDARY_EVIDENCE_VERIFIED_INTERVAL;
        interval->surface_id = line[selected->line].surface_id;
        interval->primitive_id = line[selected->line].primitive_id;
        interval->axis = -1;
        for (int end = 0; end < 2; ++end) {
            const slice_error_arr_vertex_t* point =
                endpoints[reverse ? 1 - end : end];
            for (int axis = 0; axis < 2; ++axis)
                (end ? interval->uv_end : interval->uv_start)[axis] =
                    (double)point->uv[axis];
            interval->endpoint_uncertainty[end] = (double)(
                point->uncertainty +
                16.0L * DBL_EPSILON * fmaxl(1.0L,
                    fmaxl(fabsl(point->uv[0]),
                          fabsl(point->uv[1]))));
        }
        interval->negative_side_kind = negative_face->ownership.kind;
        interval->positive_side_kind = positive_face->ownership.kind;
        interval->negative_owner_count =
            negative_face->ownership.owner_count;
        interval->positive_owner_count =
            positive_face->ownership.owner_count;
        memcpy(interval->negative_owner_cell_ids,
            negative_face->ownership.owner_cell_ids,
            interval->negative_owner_count * sizeof(int));
        memcpy(interval->positive_owner_cell_ids,
            positive_face->ownership.owner_cell_ids,
            interval->positive_owner_count * sizeof(int));
    }
    *reason = ALEA_SLICE_ERROR_RESOLVED;
    status = 1;
done:
    free(face); free(vertex); free(negative); free(positive); free(edge);
    return status;
}

typedef struct {
    const alea_system_t* sys;
    const alea_slice_view_t* view;
    const alea_transition_slice_critical_tile_t* tile;
    int axis;
    double coordinate;
    double uncertainty;
    int surface_id;
    uint32_t primitive_id;
    size_t* work;
} slice_error_occurrence_line_t;

static int slice_error_offer_occurrence_line(
    slice_error_occurrence_line_t* line, const alea_plane_data_t* world,
    int surface_id, uint32_t primitive_id) {
    long double u, v, offset;
    slice_error_project_plane(world, &line->view->plane, &u, &v, &offset);
    if (u == 0.0L && v == 0.0L) {
        const long double scale = fabsl(world->a) + fabsl(world->b) +
            fabsl(world->c) + fabsl(world->d) +
            fabsl((long double)world->a * line->view->plane.origin[0]) +
            fabsl((long double)world->b * line->view->plane.origin[1]) +
            fabsl((long double)world->c * line->view->plane.origin[2]);
        return fabsl(offset) >
            4096.0L * DBL_EPSILON * fmaxl(1.0L, scale);
    }
    const int axis = u != 0.0L && v == 0.0L ? 0
        : v != 0.0L && u == 0.0L ? 1 : -1;
    if (axis < 0 || !isfinite(offset)) return 0;
    const long double slope = axis == 0 ? u : v;
    const long double coordinate = -offset / slope;
    if (!isfinite(coordinate) || fabsl(coordinate) > DBL_MAX / 2)
        return 0;
    const double rounded = (double)coordinate;
    const long double scale = fabsl(world->a) + fabsl(world->b) +
        fabsl(world->c) + fabsl(world->d) +
        fabsl((long double)world->a * line->view->plane.origin[0]) +
        fabsl((long double)world->b * line->view->plane.origin[1]) +
        fabsl((long double)world->c * line->view->plane.origin[2]) +
        fabsl(slope * coordinate);
    const long double error = fabsl(coordinate - rounded) +
        4096.0L * DBL_EPSILON * fmaxl(1.0L, scale) / fabsl(slope);
    if (!isfinite(error) || error > DBL_MAX / 4) return 0;
    const double uncertainty = (double)error;
    if (!(rounded > line->tile->uv_min[axis] + uncertainty &&
          rounded < line->tile->uv_max[axis] - uncertainty)) return 1;
    if (line->axis >= 0) {
        if (line->axis != axis || line->coordinate != rounded) return 0;
        /* A lattice seam can be discovered from both adjacent elements. */
        if (line->surface_id == 0 && surface_id == 0) return 1;
        return line->surface_id == surface_id &&
               line->primitive_id == primitive_id;
    }
    line->axis = axis;
    line->coordinate = rounded;
    line->uncertainty = uncertainty;
    line->surface_id = surface_id;
    line->primitive_id = primitive_id;
    return 1;
}

static int slice_error_offer_local_axis_line(
    slice_error_occurrence_line_t* line, const alea_matrix_t* transform,
    int local_axis, double coordinate) {
    if (!transform || !transform->has_inverse || local_axis < 0 ||
        local_axis > 2 || !isfinite(coordinate)) return 0;
    const double* row = &transform->inv[4 * local_axis];
    const alea_plane_data_t world = {
        row[0], row[1], row[2], row[3] - coordinate};
    return slice_error_offer_occurrence_line(
        line, &world, 0, UINT32_MAX);
}

/* Collect every boundary-capable primitive, including inactive CSG branches.
 * Extra lines can only make this narrow proof decline a page. */
static int slice_error_collect_occurrence_line_node(
    slice_error_occurrence_line_t* line, alea_node_id_t id,
    const alea_matrix_t* transform, const alea_bbox_t* local_box,
    size_t depth) {
    if (!*line->work || depth > 128 || id >= line->sys->nodes.count)
        return 0;
    --*line->work;
    const alea_node_t* node = &line->sys->nodes.data[id];
    const alea_operation_t op = ALEA_GET_OPERATION(node);
    if (op != ALEA_OP_PRIMITIVE) {
        if (op != ALEA_OP_UNION && op != ALEA_OP_INTERSECTION &&
            op != ALEA_OP_DIFFERENCE && op != ALEA_OP_COMPLEMENT)
            return 0;
        if (!slice_error_collect_occurrence_line_node(
                line, node->operation.left, transform, local_box,
                depth + 1)) return 0;
        return op == ALEA_OP_COMPLEMENT ||
            slice_error_collect_occurrence_line_node(
                line, node->operation.right, transform, local_box,
                depth + 1);
    }
    size_t check_work = 2;
    const int constant = slice_error_occurrence_node_constant(
        line->sys, id, local_box, 0, &check_work);
    if (constant == 1 || constant == -1) return 1;
    const uint32_t primitive_id = node->primitive.primitive_id;
    if (primitive_id >= line->sys->primitives.count ||
        line->sys->primitives.data[primitive_id].type !=
            ALEA_PRIMITIVE_PLANE || !transform->has_inverse)
        return 0;
    alea_primitive_data_t data;
    if (!alea_primitive_copy_data(line->sys, primitive_id, &data)) return 0;
    const alea_plane_data_t* p = &data.plane;
    const double* m = transform->inv;
    const long double a = (long double)p->a * m[0] +
        (long double)p->b * m[4] + (long double)p->c * m[8];
    const long double b = (long double)p->a * m[1] +
        (long double)p->b * m[5] + (long double)p->c * m[9];
    const long double c = (long double)p->a * m[2] +
        (long double)p->b * m[6] + (long double)p->c * m[10];
    const long double d = (long double)p->a * m[3] +
        (long double)p->b * m[7] + (long double)p->c * m[11] + p->d;
    if (!isfinite(a) || !isfinite(b) || !isfinite(c) || !isfinite(d) ||
        fabsl(a) > DBL_MAX / 2 || fabsl(b) > DBL_MAX / 2 ||
        fabsl(c) > DBL_MAX / 2 || fabsl(d) > DBL_MAX / 2)
        return 0;
    const alea_plane_data_t world = {
        (double)a, (double)b, (double)c, (double)d};
    int surface_id = 0;
    for (size_t i = 0; i < line->sys->surfaces.count; ++i) {
        const alea_surface_entry_t* surface = &line->sys->surfaces.data[i];
        if (surface->primitive_id != primitive_id) continue;
        if (surface_id && surface_id != surface->mc_surface_id) return 0;
        surface_id = surface->mc_surface_id;
    }
    if (!surface_id) return 0;
    return slice_error_offer_occurrence_line(
        line, &world, surface_id, primitive_id);
}

typedef struct {
    const alea_slice_error_query_t* query;
    const alea_transition_slice_critical_tile_t* tile;
    slice_error_occurrence_line_t line;
    int complete;
} slice_error_line_occurrence_visit_t;

static int slice_error_visit_line_occurrence(
    const alea_hier_spatial_chain_hit_t* hit, void* userdata) {
    slice_error_line_occurrence_visit_t* visit = userdata;
    if (!hit->hit.is_terminal) return 0;
    if (hit->chain_truncated) {
        visit->complete = 0;
        return 1;
    }
    for (size_t level = 0; level <= hit->ancestor_count; ++level) {
        if (level < hit->ancestor_count &&
            hit->ancestor_is_lattice[level]) {
            const uint32_t lattice_index =
                hit->ancestor_cell_indices[level];
            if (lattice_index >= visit->query->sys->cells.count) {
                visit->complete = 0;
                return 1;
            }
            const alea_cell_entry_t* lattice =
                &visit->query->sys->cells.data[lattice_index];
            if (lattice->lat_type != 1) {
                visit->complete = 0;
                return 1;
            }
            const int indices[3] = {
                hit->ancestor_lattice_i[level],
                hit->ancestor_lattice_j[level],
                hit->ancestor_lattice_k[level]};
            for (int axis = 0; axis < 3; ++axis) {
                const int64_t dim =
                    (int64_t)lattice->lat_fill_dims[2 * axis + 1] -
                    lattice->lat_fill_dims[2 * axis] + 1;
                const double pitch = lattice->lat_pitch[axis];
                const int tiled = axis < 2 || dim > 1 ||
                    (lattice->lat_fill_repeating && pitch > 0.0);
                if (!tiled && pitch <= 0.0) continue;
                if (dim <= 0 || !(pitch > 0.0) || !isfinite(pitch)) {
                    visit->complete = 0;
                    return 1;
                }
                const long double offset =
                    (long double)indices[axis] -
                    lattice->lat_fill_dims[2 * axis];
                const long double lo =
                    (long double)lattice->lat_lower_left[axis] +
                    offset * pitch;
                const long double hi_bound = lo + pitch;
                if (!isfinite(lo) || !isfinite(hi_bound) ||
                    fabsl(lo) > DBL_MAX / 2 ||
                    fabsl(hi_bound) > DBL_MAX / 2 ||
                    !slice_error_offer_local_axis_line(
                        &visit->line, &hit->ancestor_transforms[level],
                        axis, (double)lo) ||
                    !slice_error_offer_local_axis_line(
                        &visit->line, &hit->ancestor_transforms[level],
                        axis, (double)hi_bound)) {
                    visit->complete = 0;
                    return 1;
                }
            }
            continue;
        }
        const uint32_t cell_index = level == hit->ancestor_count
            ? hit->hit.cell_index : hit->ancestor_cell_indices[level];
        if (cell_index >= visit->query->sys->cells.count) {
            visit->complete = 0;
            return 1;
        }
        const alea_cell_entry_t* cell =
            &visit->query->sys->cells.data[cell_index];
        if (cell->original_root_node_id != ALEA_NODE_ID_INVALID ||
            cell->lat_type != 0) {
            visit->complete = 0;
            return 1;
        }
        const alea_matrix_t* transform = level == hit->ancestor_count
            ? &hit->hit.transform : &hit->ancestor_transforms[level];
        alea_bbox_t local_box;
        if (!slice_error_occurrence_tile_box(
                &visit->query->options.view, visit->tile,
                transform, &local_box) ||
            !slice_error_collect_occurrence_line_node(
                &visit->line, cell->root_node_id,
                transform, &local_box, 0)) {
            visit->complete = 0;
            return 1;
        }
    }
    return 0;
}

static int slice_error_classify_single_occurrence_line(
    const alea_slice_error_query_t* query,
    const alea_transition_slice_critical_tile_t* tile,
    size_t contextual_bytes, alea_slice_error_page_t* out,
    alea_slice_error_unresolved_reason_t* reason,
    int* output_omitted, size_t* peak_scratch) {
    const alea_slice_error_query_options_t* options = &query->options;
    const size_t capacity =
        options->scan_options.max_exhaustive_occurrence_hits;
    if (!capacity) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
        return 0;
    }
    *peak_scratch = sizeof(slice_error_line_occurrence_visit_t);
    alea_bbox_t world_box;
    if (!slice_error_occurrence_tile_box(&options->view, tile, NULL,
                                         &world_box)) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
        return 0;
    }
    size_t work = options->scan_options.max_active_boundary_tests;
    slice_error_line_occurrence_visit_t visit = {
        .query = query, .tile = tile,
        .line = {.sys = query->sys, .view = &options->view, .tile = tile,
                 .axis = -1, .work = &work},
        .complete = 1};
    alea_hier_region_chain_status_t traversal =
        ALEA_HIER_REGION_CHAIN_UNSUPPORTED;
    size_t hit_count = 0;
    if (alea_hier_spatial_visit_region_chain_bounded(
            query->sys, &world_box, slice_error_visit_line_occurrence,
            &visit, capacity, &hit_count, &traversal) != 0)
        return -1;
    if (traversal != ALEA_HIER_REGION_CHAIN_COMPLETE) {
        *reason = traversal == ALEA_HIER_REGION_CHAIN_MAX_HITS ||
                (traversal == ALEA_HIER_REGION_CHAIN_VISITOR_STOPPED && !work)
            ? ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT
            : ALEA_SLICE_ERROR_UNRESOLVED_OCCURRENCE;
        return 0;
    }
    slice_error_occurrence_line_t line = visit.line;
    if (!visit.complete || line.axis < 0) {
        *reason = work ? ALEA_SLICE_ERROR_UNRESOLVED_OCCURRENCE
                       : ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
        return 0;
    }
    const int axis = line.axis;
    const double separation = fmax(4.0 * line.uncertainty,
        4096.0 * DBL_EPSILON * fmax(1.0, fabs(line.coordinate)));
    alea_transition_slice_critical_tile_t side[2] = {*tile, *tile};
    side[0].uv_max[axis] = line.coordinate - separation;
    side[1].uv_min[axis] = line.coordinate + separation;
    if (!(side[0].uv_max[axis] > tile->uv_min[axis] &&
          side[1].uv_min[axis] < tile->uv_max[axis])) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
        return 0;
    }
    int owners[2][ALEA_SLICE_ERROR_OWNER_CAPACITY] = {{0}};
    size_t counts[2] = {0, 0};
    for (int side_index = 0; side_index < 2; ++side_index) {
        size_t side_scratch = 0;
        const int classified = slice_error_classify_constant_occurrences(
            query, &side[side_index], 0, NULL, reason,
            output_omitted, &side_scratch,
            owners[side_index], &counts[side_index]);
        if (*peak_scratch < side_scratch) *peak_scratch = side_scratch;
        if (classified != 1) return classified;
    }
    const alea_point_coverage_kind_t kinds[2] = {
        counts[0] == 0 ? ALEA_POINT_COVERAGE_GAP
            : counts[0] == 1 ? ALEA_POINT_COVERAGE_UNIQUE
            : ALEA_POINT_COVERAGE_OVERLAP,
        counts[1] == 0 ? ALEA_POINT_COVERAGE_GAP
            : counts[1] == 1 ? ALEA_POINT_COVERAGE_UNIQUE
            : ALEA_POINT_COVERAGE_OVERLAP};
    const size_t region_count =
        slice_error_axis_is_defect(kinds[0]) +
        slice_error_axis_is_defect(kinds[1]);
    const size_t interval_count =
        slice_error_axis_is_defect(kinds[0]) !=
        slice_error_axis_is_defect(kinds[1]);
    const size_t output_limit = options->scan_options.max_output_bytes;
    if (contextual_bytes > output_limit ||
        region_count * sizeof(*out->regions) >
            output_limit - contextual_bytes ||
        interval_count * sizeof(*out->intervals) >
            output_limit - contextual_bytes -
                region_count * sizeof(*out->regions)) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
        *output_omitted = 1;
        return 0;
    }
    out->regions = region_count ? calloc(region_count, sizeof(*out->regions))
                                : NULL;
    out->intervals = interval_count
        ? calloc(interval_count, sizeof(*out->intervals)) : NULL;
    if ((region_count && !out->regions) ||
        (interval_count && !out->intervals)) return -1;
    for (int s = 0; s < 2; ++s) {
        if (!slice_error_axis_is_defect(kinds[s])) continue;
        alea_slice_error_region_t* region =
            &out->regions[out->region_count++];
        memcpy(region->uv_min, tile->uv_min, sizeof(tile->uv_min));
        memcpy(region->uv_max, tile->uv_max, sizeof(tile->uv_max));
        if (s == 0) region->uv_max[axis] = line.coordinate;
        else region->uv_min[axis] = line.coordinate;
        if (s == 0) region->uv_max_uncertainty[axis] = line.uncertainty;
        else region->uv_min_uncertainty[axis] = line.uncertainty;
        region->kind = kinds[s];
        region->owner_count = counts[s];
        memcpy(region->owner_cell_ids, owners[s],
               counts[s] * sizeof(*owners[s]));
    }
    if (interval_count) {
        alea_slice_error_interval_t* interval = &out->intervals[0];
        out->interval_count = 1;
        interval->evidence_scope =
            ALEA_SLICE_BOUNDARY_EVIDENCE_VERIFIED_INTERVAL;
        interval->surface_id = line.surface_id;
        interval->primitive_id = line.primitive_id;
        interval->axis = axis;
        interval->uv_start[axis] = interval->uv_end[axis] =
            line.coordinate;
        interval->uv_start[1 - axis] = tile->uv_min[1 - axis];
        interval->uv_end[1 - axis] = tile->uv_max[1 - axis];
        interval->endpoint_uncertainty[0] = line.uncertainty;
        interval->endpoint_uncertainty[1] = line.uncertainty;
        interval->negative_side_kind = kinds[0];
        interval->positive_side_kind = kinds[1];
        interval->negative_owner_count = counts[0];
        interval->positive_owner_count = counts[1];
        memcpy(interval->negative_owner_cell_ids, owners[0],
               counts[0] * sizeof(*owners[0]));
        memcpy(interval->positive_owner_cell_ids, owners[1],
               counts[1] * sizeof(*owners[1]));
    }
    *reason = ALEA_SLICE_ERROR_RESOLVED;
    return 1;
}

/* Returns 1 for a complete proof, 0 with a whole-tile unresolved reason for
 * unsupported/numerically uncertain work, or -1 on allocation failure. */
static int slice_error_classify_axis_tile(
    const alea_slice_error_query_t* query,
    const alea_transition_slice_critical_tile_t* tile,
    size_t contextual_bytes, alea_slice_error_page_t* out,
    alea_slice_error_unresolved_reason_t* reason,
    int* output_omitted, size_t* peak_scratch) {
    *output_omitted = 0;
    *peak_scratch = 0;
    alea_system_t* sys = query->sys;
    const alea_slice_error_query_options_t* options = &query->options;
    size_t primitive_count = 0;
    size_t cell_count = 0;
    if (!slice_error_axis_view_supported(&options->view)) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_SLICE_FRAME;
        return 0;
    }
    const size_t candidate_limit = options->scan_options.max_curves_per_tile;
    const size_t scratch_limit =
        options->scan_options.max_critical_scratch_bytes;
    size_t cell_capacity = sys->cells.count;
    if (cell_capacity > candidate_limit) cell_capacity = candidate_limit;
    if (cell_capacity > scratch_limit / sizeof(size_t))
        cell_capacity = scratch_limit / sizeof(size_t);
    size_t* cells = cell_capacity
        ? calloc(cell_capacity, sizeof(*cells)) : NULL;
    if (cell_capacity && !cells) return -1;
    *peak_scratch = cell_capacity * sizeof(*cells);
    int status = 0;
    *reason = ALEA_SLICE_ERROR_UNRESOLVED_PLANAR_ARRANGEMENT;
    size_t discovery_work = options->scan_options.max_active_boundary_tests;
    const int x_slice_axis = slice_error_slice_axis_for_world(
        &options->view, 0);
    const double x_slice_sign = slice_error_slice_axis_sign(
        &options->view, 0);
    const double world_max_x = options->view.plane.origin[0] +
        (x_slice_axis < 0 ? 0.0 :
         x_slice_sign * (x_slice_sign > 0.0
             ? tile->uv_max[x_slice_axis] : tile->uv_min[x_slice_axis]));
    const double x_margin = 64.0 * DBL_EPSILON *
        fmax(1.0, fabs(world_max_x));
    const double x_cutoff = world_max_x + x_margin;
    for (size_t i = 0; i < query->indexed_count; ++i) {
        const slice_error_index_cell_t* indexed =
            &query->indexed_cells[i];
        if (isfinite(x_cutoff) && indexed->box.min_x > x_cutoff)
            break;
        if (slice_error_box_outside_tile(&indexed->box, options, tile))
            continue;
        if (cell_count == cell_capacity) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
            goto selection_done;
        }
        cells[cell_count++] = indexed->cell_index;
    }
    for (size_t i = 0; i < query->uncertain_count; ++i) {
        if (cell_count == cell_capacity) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
            goto selection_done;
        }
        cells[cell_count++] = query->uncertain_cells[i];
    }
    for (size_t ci = query->uncached_cell_begin;
         ci < sys->cells.count; ++ci) {
        const alea_cell_entry_t* cell = &sys->cells.data[ci];
        if (cell->universe_id != 0) continue;
        if (cell->original_root_node_id == ALEA_NODE_ID_INVALID &&
            slice_error_bbox_tree_safe(sys, cell->root_node_id, 0,
                                       &discovery_work)) {
            alea_bbox_t box = alea_get_bbox(sys, cell->root_node_id);
            if (box.min_x <= box.max_x && box.min_y <= box.max_y &&
                box.min_z <= box.max_z) {
                box.min_x = slice_error_trusted_lower(box.min_x);
                box.max_x = slice_error_trusted_upper(box.max_x);
                box.min_y = slice_error_trusted_lower(box.min_y);
                box.max_y = slice_error_trusted_upper(box.max_y);
                box.min_z = slice_error_trusted_lower(box.min_z);
                box.max_z = slice_error_trusted_upper(box.max_z);
                if (slice_error_box_outside_tile(&box, options, tile))
                    continue;
            }
        }
        if (cell_count == cell_capacity) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
            goto selection_done;
        }
        cells[cell_count++] = ci;
    }
    if (cell_count > 1)
        qsort(cells, cell_count, sizeof(*cells), slice_error_size_compare);
    if (cell_count == 0) {
        const size_t output_limit = options->scan_options.max_output_bytes;
        if (contextual_bytes > output_limit ||
            sizeof(alea_slice_error_region_t) >
                output_limit - contextual_bytes) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
            *output_omitted = 1;
            goto selection_done;
        }
        out->regions = calloc(1, sizeof(*out->regions));
        if (!out->regions) { status = -1; goto selection_done; }
        out->region_count = 1;
        memcpy(out->regions[0].uv_min, tile->uv_min,
               sizeof(tile->uv_min));
        memcpy(out->regions[0].uv_max, tile->uv_max,
               sizeof(tile->uv_max));
        out->regions[0].kind = ALEA_POINT_COVERAGE_GAP;
        *reason = ALEA_SLICE_ERROR_RESOLVED;
        status = 1;
        goto selection_done;
    }
    if (cell_count == 1) {
        status = slice_error_classify_isolated_circle(
            sys, &options->view, tile, cells[0], contextual_bytes,
            options->scan_options.max_output_bytes, out);
        if (status == 1) {
            *reason = ALEA_SLICE_ERROR_RESOLVED;
            goto selection_done;
        }
        if (status == 2) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
            *output_omitted = 1;
            status = 0;
            goto selection_done;
        }
        if (status < 0) goto selection_done;
    }
    if (cell_count == 2) {
        size_t mixed_scratch = 0;
        status = slice_error_classify_mixed_line_circle(
            sys, &options->view, tile, cells,
            scratch_limit - *peak_scratch, &mixed_scratch, contextual_bytes,
            options->scan_options.max_output_bytes, out);
        *peak_scratch += mixed_scratch;
        if (status == 1) {
            *reason = ALEA_SLICE_ERROR_RESOLVED;
            goto selection_done;
        }
        if (status == 2) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
            *output_omitted = 1;
            status = 0;
            goto selection_done;
        }
        if (status == 3) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
            status = 0;
            goto selection_done;
        }
        if (status < 0) goto selection_done;
        if (mixed_scratch) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_PLANAR_ARRANGEMENT;
            goto selection_done;
        }
        size_t pair_scratch = 0;
        status = slice_error_classify_crossing_circles(
            sys, &options->view, tile, cells,
            scratch_limit - *peak_scratch, &pair_scratch, contextual_bytes,
            options->scan_options.max_output_bytes, out);
        *peak_scratch += pair_scratch;
        if (status == 1) {
            *reason = ALEA_SLICE_ERROR_RESOLVED;
            goto selection_done;
        }
        if (status == 2) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
            *output_omitted = 1;
            status = 0;
            goto selection_done;
        }
        if (status == 3) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
            status = 0;
            goto selection_done;
        }
        if (status < 0) goto selection_done;
    }
    size_t plane_capacity = sys->primitives.count;
    if (plane_capacity > candidate_limit) plane_capacity = candidate_limit;
    const size_t remaining_scratch = scratch_limit -
        cell_capacity * sizeof(*cells);
    const size_t per_plane = sizeof(slice_error_axis_plane_t) +
        2 * sizeof(slice_error_grid_line_t);
    if (remaining_scratch < 4 * sizeof(slice_error_grid_line_t) ||
        plane_capacity > (remaining_scratch -
            4 * sizeof(slice_error_grid_line_t)) / per_plane)
        plane_capacity = remaining_scratch <
            4 * sizeof(slice_error_grid_line_t) ? 0 :
            (remaining_scratch - 4 * sizeof(slice_error_grid_line_t)) /
                per_plane;
    if (!plane_capacity) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
        goto selection_done;
    }
    slice_error_axis_plane_t* planes =
        calloc(plane_capacity, sizeof(*planes));
    slice_error_grid_line_t* x =
        calloc(plane_capacity + 2, sizeof(*x));
    slice_error_grid_line_t* y =
        calloc(plane_capacity + 2, sizeof(*y));
    if (!planes || !x || !y) {
        free(planes); free(x); free(y);
        status = -1;
        goto selection_done;
    }
    const size_t plane_peak_scratch = cell_capacity * sizeof(*cells) +
        plane_capacity * sizeof(*planes) +
        2 * (plane_capacity + 2) * sizeof(*x);
    if (*peak_scratch < plane_peak_scratch)
        *peak_scratch = plane_peak_scratch;
    size_t nx = 1, ny = 1;
    for (size_t i = 0; i < cell_count; ++i) {
        const alea_cell_entry_t* cell = &sys->cells.data[cells[i]];
        if (cell->fill_universe > 0 || cell->lat_type != 0 ||
            cell->original_root_node_id != ALEA_NODE_ID_INVALID) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_OCCURRENCE;
            goto done;
        }
        if (!slice_error_collect_primitives(
                sys, cell->root_node_id, 0, planes, &primitive_count,
                plane_capacity, &discovery_work)) {
            *reason = !discovery_work || primitive_count == plane_capacity
                ? ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT
                : ALEA_SLICE_ERROR_UNRESOLVED_UNSUPPORTED_GEOMETRY;
            goto done;
        }
    }
    for (size_t si = 0; si < sys->surfaces.count; ++si) {
        const alea_surface_entry_t* surface = &sys->surfaces.data[si];
        const int index = slice_error_plane_index(
            planes, primitive_count, surface->primitive_id);
        if (index < 0) continue;
        slice_error_axis_plane_t* plane = &planes[index];
        if (plane->surface_id != 0) {
            /* Primitive deduplication defines the system's geometry
             * equivalence. Multiple card IDs for that primitive do not add a
             * second physical boundary to the ownership arrangement. */
            continue;
        }
        plane->surface_id = surface->mc_surface_id;
    }
    slice_error_axis_plane_t selected[2];
    size_t line_count = 0, variable_count = 0;
    int oblique = 0, tilted = 0, all_supported = 1;
    for (size_t pi = 0; pi < primitive_count; ++pi) {
        const alea_primitive_entry_t* primitive =
            &sys->primitives.data[planes[pi].primitive_id];
        if (primitive->type == ALEA_PRIMITIVE_SPHERE ||
            primitive->type == ALEA_PRIMITIVE_SPH) {
            alea_sphere_data_t sphere;
            if (primitive->type == ALEA_PRIMITIVE_SPHERE) {
                if (primitive->payload_index >= sys->primitive_spheres.count) {
                    *reason = ALEA_SLICE_ERROR_UNRESOLVED_PRIMITIVE;
                    goto done;
                }
                sphere = sys->primitive_spheres.data[primitive->payload_index];
            } else {
                if (primitive->payload_index >= sys->primitive_sphs.count) {
                    *reason = ALEA_SLICE_ERROR_UNRESOLVED_PRIMITIVE;
                    goto done;
                }
                const alea_sph_data_t* sph =
                    &sys->primitive_sphs.data[primitive->payload_index];
                sphere = (alea_sphere_data_t){sph->center_x, sph->center_y,
                                              sph->center_z, sph->radius};
            }
            int sign;
            if (!slice_error_sphere_tile_sign(
                    &sphere, &options->view, tile, &sign)) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_PRIMITIVE;
                goto done;
            }
            planes[pi].axis = 2;
            planes[pi].coefficient_sign = sign;
            continue;
        }
        if (primitive->type != ALEA_PRIMITIVE_PLANE ||
            primitive->payload_index >= sys->primitive_planes.count ||
            planes[pi].surface_id <= 0) {
            all_supported = 0;
            break;
        }
        const alea_plane_data_t* p =
            &sys->primitive_planes.data[primitive->payload_index];
        if (!isfinite(p->a) || !isfinite(p->b) ||
            !isfinite(p->c) || !isfinite(p->d)) {
            all_supported = 0;
            break;
        }
        long double u, v, offset;
        slice_error_project_plane(p, &options->view.plane,
                                  &u, &v, &offset);
        if (!isfinite(u) || !isfinite(v) || !isfinite(offset)) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            goto done;
        }
        if (u == 0.0L && v == 0.0L) {
            const long double magnitude = fabsl((long double)p->d) +
                fabsl((long double)p->a * options->view.plane.origin[0]) +
                fabsl((long double)p->b * options->view.plane.origin[1]) +
                fabsl((long double)p->c * options->view.plane.origin[2]);
            const long double uncertainty = 32.0L * DBL_EPSILON *
                fmaxl(1.0L, magnitude);
            if (!(fabsl(offset) > uncertainty)) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                goto done;
            }
            planes[pi].axis = 2;
            planes[pi].coefficient_sign = offset > 0.0L ? 1 : -1;
            continue;
        }
        ++variable_count;
        if (line_count < 2) selected[line_count++] = planes[pi];
        oblique += u != 0.0L && v != 0.0L;
        tilted += ((p->a != 0.0) + (p->b != 0.0) +
                   (p->c != 0.0)) > 1;
    }
    if (all_supported && variable_count == 1 && (oblique || tilted)) {
        status = slice_error_classify_single_oblique_tile(
            query, tile, cells, cell_count, &selected[0],
            planes, primitive_count, contextual_bytes,
            &discovery_work, out, reason, output_omitted);
        goto done;
    }
    if (all_supported && variable_count == 2 && (oblique || tilted)) {
        if (oblique == 2) {
            status = slice_error_classify_parallel_oblique_tile(
                query, tile, cells, cell_count, selected,
                planes, primitive_count, contextual_bytes,
                &discovery_work, out, reason, output_omitted);
        }
        if (oblique < 2 ||
            (status == 0 && *reason ==
             ALEA_SLICE_ERROR_UNRESOLVED_PLANAR_ARRANGEMENT))
            status = slice_error_classify_crossing_oblique_tile(
                query, tile, cells, cell_count, selected,
                planes, primitive_count, contextual_bytes,
                &discovery_work, out, reason, output_omitted);
        goto done;
    }
    if (all_supported && variable_count >= 3 && (oblique || tilted)) {
        status = slice_error_classify_line_arrangement(
            query, tile, cells, cell_count, planes, primitive_count,
            contextual_bytes, &discovery_work, out, reason,
            output_omitted, peak_scratch);
        goto done;
    }
    for (size_t pi = 0; pi < primitive_count; ++pi) {
        const uint32_t primitive_id = planes[pi].primitive_id;
        const alea_primitive_entry_t* primitive =
            &sys->primitives.data[primitive_id];
        if (planes[pi].axis == 2 &&
            (primitive->type == ALEA_PRIMITIVE_SPHERE ||
             primitive->type == ALEA_PRIMITIVE_SPH)) continue;
        if (primitive->type != ALEA_PRIMITIVE_PLANE ||
            primitive->payload_index >= sys->primitive_planes.count ||
            planes[pi].surface_id <= 0) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_PRIMITIVE;
            goto done;
        }
        const alea_plane_data_t* p =
            &sys->primitive_planes.data[primitive->payload_index];
        if (!isfinite(p->a) || !isfinite(p->b) ||
            !isfinite(p->c) || !isfinite(p->d))
            goto done;
        const int world_axis = p->a != 0.0 && p->b == 0.0 && p->c == 0.0 ? 0
            : p->b != 0.0 && p->a == 0.0 && p->c == 0.0 ? 1
            : p->c != 0.0 && p->a == 0.0 && p->b == 0.0 ? 2 : -1;
        if (world_axis < 0) goto done;
        const int axis = slice_error_slice_axis_for_world(
            &options->view, world_axis);
        const double coefficient = world_axis == 0 ? p->a
            : world_axis == 1 ? p->b : p->c;
        if (axis < 0) {
            const long double term = (long double)coefficient *
                (long double)options->view.plane.origin[world_axis];
            const long double value = term + (long double)p->d;
            const long double uncertainty =
                32.0L * DBL_EPSILON *
                fmaxl(1.0L, fabsl(term) + fabsl((long double)p->d));
            if (!isfinite(value) ||
                !(fabsl(value) > uncertainty)) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                goto done;
            }
            planes[pi].axis = 2;
            planes[pi].coefficient_sign = value > 0.0L ? 1 : -1;
            continue;
        }
        const double slice_sign = slice_error_slice_axis_sign(
            &options->view, world_axis);
        const long double quotient =
            -(long double)p->d / (long double)coefficient;
        const long double origin =
            (long double)options->view.plane.origin[world_axis];
        const long double coordinate =
            (quotient - origin) / (long double)slice_sign;
        if (!isfinite(coordinate) || fabsl(coordinate) > DBL_MAX) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            goto done;
        }
        const double coord = (double)coordinate;
        const long double uncertainty_ld =
            64.0L * LDBL_EPSILON *
                fmaxl(1.0L, fmaxl(fabsl(quotient), fabsl(origin))) +
            fabsl(coordinate - (long double)coord) +
            8.0L * DBL_EPSILON * fmaxl(1.0L, fabsl(coordinate));
        if (!isfinite(uncertainty_ld) || uncertainty_ld > DBL_MAX) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            goto done;
        }
        const double uncertainty = (double)uncertainty_ld;
        planes[pi] = (slice_error_axis_plane_t){
            .axis = axis,
            .coefficient_sign = coefficient * slice_sign > 0.0 ? 1 : -1,
            .coordinate = coord,
            .uncertainty = uncertainty,
            .surface_id = planes[pi].surface_id,
            .primitive_id = primitive_id
        };
        const double lo = tile->uv_min[axis], hi = tile->uv_max[axis];
        if ((fabsl(coordinate - (long double)lo) <= uncertainty &&
             coordinate != (long double)lo) ||
            (fabsl(coordinate - (long double)hi) <= uncertainty &&
             coordinate != (long double)hi)) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            goto done;
        }
    }
    double support_min[2] = {tile->uv_min[0], tile->uv_min[1]};
    for (int axis = 0; axis < 2; ++axis) {
        const double lo = tile->uv_min[axis];
        if (lo == options->required_uv_min[axis]) continue;
        int seam_plane = 0;
        double previous = options->required_uv_min[axis];
        double previous_uncertainty = 0.0;
        for (size_t pi = 0; pi < primitive_count; ++pi) {
            const slice_error_axis_plane_t* plane = &planes[pi];
            if (plane->axis == 2) continue;
            if (plane->axis != axis) continue;
            if (plane->coordinate == lo) seam_plane = 1;
            if (plane->coordinate < lo &&
                plane->coordinate > previous) {
                previous = plane->coordinate;
                previous_uncertainty = plane->uncertainty;
            }
        }
        if (!seam_plane) continue;
        support_min[axis] = previous + (lo - previous) * 0.5;
        if (!(support_min[axis] > previous && support_min[axis] < lo) ||
            support_min[axis] - previous <= previous_uncertainty +
                4.0 * DBL_EPSILON * fmax(1.0, fabs(previous))) {
            *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
            goto done;
        }
    }
    x[0].coordinate = support_min[0];
    y[0].coordinate = support_min[1];
    for (size_t pi = 0; pi < primitive_count; ++pi) {
        const slice_error_axis_plane_t* plane = &planes[pi];
        if (plane->axis == 2) continue;
        const int axis = plane->axis;
        const double coord = plane->coordinate;
        if (coord > support_min[axis] &&
            coord < tile->uv_max[axis]) {
            slice_error_grid_line_t line = {
                .coordinate = coord,
                .uncertainty = plane->uncertainty,
                .surface_id = plane->surface_id,
                .primitive_id = plane->primitive_id
            };
            if (axis == 0) x[nx++] = line;
            else y[ny++] = line;
        }
    }
    x[nx++].coordinate = tile->uv_max[0];
    y[ny++].coordinate = tile->uv_max[1];
    qsort(x, nx, sizeof(*x), slice_error_grid_line_compare);
    qsort(y, ny, sizeof(*y), slice_error_grid_line_compare);
    if (!slice_error_axis_grid_valid(x, nx) ||
        !slice_error_axis_grid_valid(y, ny)) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
        goto done;
    }
    const size_t x_faces = nx - 1, y_faces = ny - 1;
    if (x_faces > SIZE_MAX / y_faces ||
        x_faces * y_faces > options->scan_options.max_active_boundary_tests ||
        (nx - 2) > options->scan_options.max_curve_pairs /
            (ny - 2 ? ny - 2 : 1)) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
        goto done;
    }
    const size_t face_count = x_faces * y_faces;
    if (face_count > discovery_work / cell_count) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
        goto done;
    }
    const size_t scratch_bytes = cell_capacity * sizeof(*cells) +
        plane_capacity * sizeof(*planes) +
        2 * (plane_capacity + 2) * sizeof(*x);
    if (face_count > SIZE_MAX / sizeof(slice_error_face_t) ||
        scratch_bytes > options->scan_options.max_critical_scratch_bytes ||
        face_count * sizeof(slice_error_face_t) >
            options->scan_options.max_critical_scratch_bytes - scratch_bytes) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
        goto done;
    }
    *peak_scratch = scratch_bytes +
        face_count * sizeof(slice_error_face_t);
    slice_error_face_t* faces = calloc(face_count, sizeof(*faces));
    if (!faces) { status = -1; goto done; }
    size_t work_remaining = discovery_work;
    for (size_t yi = 0; yi < y_faces; ++yi) {
        for (size_t xi = 0; xi < x_faces; ++xi) {
            slice_error_face_t* face = &faces[yi * x_faces + xi];
            const double u = x[xi].coordinate +
                0.5 * (x[xi + 1].coordinate - x[xi].coordinate);
            const double v = y[yi].coordinate +
                0.5 * (y[yi + 1].coordinate - y[yi].coordinate);
            double world[3];
            for (int axis = 0; axis < 3; ++axis)
                world[axis] = options->view.plane.origin[axis] +
                    options->view.plane.u_axis[axis] * u +
                    options->view.plane.v_axis[axis] * v;
            if (!(u > x[xi].coordinate && u < x[xi + 1].coordinate) ||
                !(v > y[yi].coordinate && v < y[yi + 1].coordinate) ||
                !isfinite(world[0]) || !isfinite(world[1]) ||
                !isfinite(world[2])) {
                *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                free(faces); goto done;
            }
            for (size_t ci = 0; ci < cell_count; ++ci) {
                const alea_cell_entry_t* cell = &sys->cells.data[cells[ci]];
                int inside = 0;
                if (!slice_error_axis_node_inside(
                        sys, cell->root_node_id, planes, primitive_count,
                        x, xi, y, yi, 0, &work_remaining, &inside)) {
                    *reason = work_remaining == 0
                        ? ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT
                        : ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                    free(faces); goto done;
                }
                if (inside != alea_point_inside(
                        sys, cell->root_node_id,
                        world[0], world[1], world[2])) {
                    *reason = ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
                    free(faces); goto done;
                }
                if (inside) {
                    if (face->owner_count == ALEA_SLICE_ERROR_OWNER_CAPACITY) {
                        *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
                        free(faces); goto done;
                    }
                    face->owner_cell_ids[face->owner_count++] =
                        cell->mc_cell_id;
                }
            }
            face->kind = face->owner_count == 0
                ? ALEA_POINT_COVERAGE_GAP
                : face->owner_count == 1
                    ? ALEA_POINT_COVERAGE_UNIQUE
                    : ALEA_POINT_COVERAGE_OVERLAP;
        }
    }
    size_t region_count = 0, interval_count = 0;
    for (size_t yi = 0; yi < y_faces; ++yi)
        for (size_t xi = 0; xi < x_faces; ++xi)
            if (x[xi + 1].coordinate > tile->uv_min[0] &&
                x[xi].coordinate < tile->uv_max[0] &&
                y[yi + 1].coordinate > tile->uv_min[1] &&
                y[yi].coordinate < tile->uv_max[1])
                region_count += slice_error_axis_is_defect(
                    faces[yi * x_faces + xi].kind);
    for (size_t yi = 0; yi < y_faces; ++yi)
        for (size_t xi = 1; xi < x_faces; ++xi) {
            if (x[xi].coordinate < tile->uv_min[0] ||
                x[xi].coordinate >= tile->uv_max[0] ||
                y[yi + 1].coordinate <= tile->uv_min[1] ||
                y[yi].coordinate >= tile->uv_max[1]) continue;
            const slice_error_face_t* a = &faces[yi * x_faces + xi - 1];
            const slice_error_face_t* b = &faces[yi * x_faces + xi];
            interval_count += slice_error_axis_is_defect(a->kind) !=
                              slice_error_axis_is_defect(b->kind);
        }
    for (size_t yi = 1; yi < y_faces; ++yi)
        for (size_t xi = 0; xi < x_faces; ++xi) {
            if (y[yi].coordinate < tile->uv_min[1] ||
                y[yi].coordinate >= tile->uv_max[1] ||
                x[xi + 1].coordinate <= tile->uv_min[0] ||
                x[xi].coordinate >= tile->uv_max[0]) continue;
            const slice_error_face_t* a = &faces[(yi - 1) * x_faces + xi];
            const slice_error_face_t* b = &faces[yi * x_faces + xi];
            interval_count += slice_error_axis_is_defect(a->kind) !=
                              slice_error_axis_is_defect(b->kind);
        }
    const size_t output_limit = options->scan_options.max_output_bytes;
    if (contextual_bytes > output_limit ||
        region_count > (output_limit - contextual_bytes) /
            sizeof(alea_slice_error_region_t) ||
        interval_count > (output_limit - contextual_bytes -
            region_count * sizeof(alea_slice_error_region_t)) /
            sizeof(alea_slice_error_interval_t)) {
        *reason = ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
        *output_omitted = 1;
        free(faces); goto done;
    }
    out->regions = region_count
        ? calloc(region_count, sizeof(*out->regions)) : NULL;
    out->intervals = interval_count
        ? calloc(interval_count, sizeof(*out->intervals)) : NULL;
    if ((region_count && !out->regions) ||
        (interval_count && !out->intervals)) {
        free(faces); status = -1; goto done;
    }
    for (size_t yi = 0; yi < y_faces; ++yi)
        for (size_t xi = 0; xi < x_faces; ++xi) {
            const slice_error_face_t* face = &faces[yi * x_faces + xi];
            if (!slice_error_axis_is_defect(face->kind) ||
                x[xi + 1].coordinate <= tile->uv_min[0] ||
                x[xi].coordinate >= tile->uv_max[0] ||
                y[yi + 1].coordinate <= tile->uv_min[1] ||
                y[yi].coordinate >= tile->uv_max[1]) continue;
            alea_slice_error_region_t* region =
                &out->regions[out->region_count++];
            region->uv_min[0] = fmax(x[xi].coordinate, tile->uv_min[0]);
            region->uv_max[0] = fmin(x[xi + 1].coordinate, tile->uv_max[0]);
            region->uv_min[1] = fmax(y[yi].coordinate, tile->uv_min[1]);
            region->uv_max[1] = fmin(y[yi + 1].coordinate, tile->uv_max[1]);
            region->uv_min_uncertainty[0] =
                region->uv_min[0] == x[xi].coordinate
                    ? x[xi].uncertainty : 0.0;
            region->uv_max_uncertainty[0] =
                region->uv_max[0] == x[xi + 1].coordinate
                    ? x[xi + 1].uncertainty : 0.0;
            region->uv_min_uncertainty[1] =
                region->uv_min[1] == y[yi].coordinate
                    ? y[yi].uncertainty : 0.0;
            region->uv_max_uncertainty[1] =
                region->uv_max[1] == y[yi + 1].coordinate
                    ? y[yi + 1].uncertainty : 0.0;
            region->kind = face->kind;
            region->owner_count = face->owner_count;
            memcpy(region->owner_cell_ids, face->owner_cell_ids,
                   face->owner_count * sizeof(int));
        }
    for (int axis = 0; axis < 2; ++axis) {
        const size_t major_faces = axis == 0 ? x_faces : y_faces;
        const size_t minor_faces = axis == 0 ? y_faces : x_faces;
        const slice_error_grid_line_t* major = axis == 0 ? x : y;
        const slice_error_grid_line_t* minor = axis == 0 ? y : x;
        for (size_t mi = 1; mi < major_faces; ++mi)
            for (size_t ni = 0; ni < minor_faces; ++ni) {
                if (major[mi].coordinate < tile->uv_min[axis] ||
                    major[mi].coordinate >= tile->uv_max[axis] ||
                    minor[ni + 1].coordinate <= tile->uv_min[1 - axis] ||
                    minor[ni].coordinate >= tile->uv_max[1 - axis]) continue;
                const size_t negative_index = axis == 0
                    ? ni * x_faces + mi - 1
                    : (mi - 1) * x_faces + ni;
                const size_t positive_index = axis == 0
                    ? ni * x_faces + mi
                    : mi * x_faces + ni;
                const slice_error_face_t* negative =
                    &faces[negative_index];
                const slice_error_face_t* positive =
                    &faces[positive_index];
                if (slice_error_axis_is_defect(negative->kind) ==
                    slice_error_axis_is_defect(positive->kind)) continue;
                alea_slice_error_interval_t* interval =
                    &out->intervals[out->interval_count++];
                interval->evidence_scope =
                    ALEA_SLICE_BOUNDARY_EVIDENCE_VERIFIED_INTERVAL;
                interval->surface_id = major[mi].surface_id;
                interval->primitive_id = major[mi].primitive_id;
                interval->axis = axis;
                interval->uv_start[axis] = major[mi].coordinate;
                interval->uv_end[axis] = major[mi].coordinate;
                interval->uv_start[1 - axis] = fmax(
                    minor[ni].coordinate, tile->uv_min[1 - axis]);
                interval->uv_end[1 - axis] = fmin(
                    minor[ni + 1].coordinate, tile->uv_max[1 - axis]);
                interval->endpoint_uncertainty[0] = fmax(
                    major[mi].uncertainty, minor[ni].uncertainty);
                interval->endpoint_uncertainty[1] = fmax(
                    major[mi].uncertainty, minor[ni + 1].uncertainty);
                interval->negative_side_kind = negative->kind;
                interval->positive_side_kind = positive->kind;
                interval->negative_owner_count = negative->owner_count;
                interval->positive_owner_count = positive->owner_count;
                memcpy(interval->negative_owner_cell_ids,
                       negative->owner_cell_ids,
                       negative->owner_count * sizeof(int));
                memcpy(interval->positive_owner_cell_ids,
                       positive->owner_cell_ids,
                       positive->owner_count * sizeof(int));
            }
    }
    free(faces);
    *reason = ALEA_SLICE_ERROR_RESOLVED;
    status = 1;
done:
    free(planes); free(x); free(y);
    if (status != 1) {
        free(out->intervals); out->intervals = NULL;
        free(out->circles); out->circles = NULL;
        free(out->regions); out->regions = NULL;
        out->interval_count = out->circle_count = out->region_count = 0;
    }
selection_done:
    free(cells);
    return status;
}

typedef struct {
    alea_system_t* sys;
    const alea_slice_view_t* view;
    const alea_transition_slice_critical_tile_t* tile;
    size_t candidate_count;
    int failed;
} slice_error_density_count_t;

static int slice_error_count_candidate(
    const alea_spatial_hit_t* candidate, void* userdata) {
    (void)candidate;
    (void)userdata;
    return 0;
}

static int slice_error_count_occurrence(
    const alea_hier_spatial_chain_hit_t* occurrence, void* userdata) {
    slice_error_density_count_t* count = userdata;
    alea_bbox_t local_box;
    if (!slice_error_occurrence_tile_box(
            count->view, count->tile, &occurrence->hit.transform,
            &local_box)) {
        count->failed = 1;
        return 1;
    }
    size_t visited = 0;
    if (alea_hier_spatial_visit_universe_region(
            count->sys, occurrence->hit.universe_id, &local_box,
            slice_error_count_candidate, NULL, &visited) != 0 ||
        count->candidate_count > SIZE_MAX - visited) {
        count->failed = 1;
        return 1;
    }
    count->candidate_count += visited;
    return 0;
}

static void slice_error_uniform_critical_tiles(
    const alea_transition_slice_critical_tile_t* page, size_t divisions,
    alea_transition_slice_critical_tile_t* tiles) {
    for (size_t j = 0; j < divisions; ++j) {
        for (size_t i = 0; i < divisions; ++i) {
            alea_transition_slice_critical_tile_t* part =
                &tiles[j * divisions + i];
            const double fu0 = (double)i / (double)divisions;
            const double fu1 = (double)(i + 1u) / (double)divisions;
            const double fv0 = (double)j / (double)divisions;
            const double fv1 = (double)(j + 1u) / (double)divisions;
            part->uv_min[0] = page->uv_min[0] +
                (page->uv_max[0] - page->uv_min[0]) * fu0;
            part->uv_max[0] = page->uv_min[0] +
                (page->uv_max[0] - page->uv_min[0]) * fu1;
            part->uv_min[1] = page->uv_min[1] +
                (page->uv_max[1] - page->uv_min[1]) * fv0;
            part->uv_max[1] = page->uv_min[1] +
                (page->uv_max[1] - page->uv_min[1]) * fv1;
        }
    }
}

static size_t slice_error_page_output_bytes(
    const alea_slice_error_page_t* page) {
    if (!page) return 0;
    if (page->finding_count > SIZE_MAX / sizeof(*page->findings) ||
        page->witness_count > SIZE_MAX / sizeof(*page->witnesses) ||
        page->interval_count > SIZE_MAX / sizeof(*page->intervals) ||
        page->circle_count > SIZE_MAX / sizeof(*page->circles) ||
        page->region_count > SIZE_MAX / sizeof(*page->regions))
        return SIZE_MAX;
    size_t bytes = page->finding_count * sizeof(*page->findings);
    const size_t witness_bytes =
        page->witness_count * sizeof(*page->witnesses);
    const size_t interval_bytes =
        page->interval_count * sizeof(*page->intervals);
    const size_t circle_bytes = page->circle_count * sizeof(*page->circles);
    const size_t region_bytes = page->region_count * sizeof(*page->regions);
    if (bytes > SIZE_MAX - witness_bytes) return SIZE_MAX;
    bytes += witness_bytes;
    if (bytes > SIZE_MAX - interval_bytes) return SIZE_MAX;
    bytes += interval_bytes;
    if (bytes > SIZE_MAX - circle_bytes) return SIZE_MAX;
    bytes += circle_bytes;
    if (bytes > SIZE_MAX - region_bytes) return SIZE_MAX;
    return bytes + region_bytes;
}

/* Append a classifier's independently allocated products. Return zero when
 * the caller's output budget cannot retain the complete set. */
static int slice_error_append_products(
    alea_slice_error_page_t* destination,
    const alea_slice_error_page_t* source, size_t output_limit) {
    const size_t current = slice_error_page_output_bytes(destination);
    const size_t added = slice_error_page_output_bytes(source);
    if (current == SIZE_MAX || added == SIZE_MAX || current > output_limit ||
        added > output_limit - current)
        return 0;
#define APPEND_SLICE_PRODUCTS(member, count_member) do {                  \
        if (source->count_member) {                                       \
            if (destination->count_member >                               \
                    SIZE_MAX - source->count_member) return -1;            \
            const size_t total = destination->count_member +              \
                source->count_member;                                     \
            if (total > SIZE_MAX / sizeof(*destination->member))          \
                return -1;                                                \
            void* next = realloc(destination->member,                     \
                                 total * sizeof(*destination->member));    \
            if (!next) return -1;                                         \
            destination->member = next;                                   \
            memcpy(destination->member + destination->count_member,       \
                   source->member,                                        \
                   source->count_member * sizeof(*source->member));       \
            destination->count_member = total;                            \
        }                                                                 \
    } while (0)
    APPEND_SLICE_PRODUCTS(intervals, interval_count);
    APPEND_SLICE_PRODUCTS(circles, circle_count);
    APPEND_SLICE_PRODUCTS(regions, region_count);
#undef APPEND_SLICE_PRODUCTS
    return 1;
}

static int slice_error_append_unresolved_region(
    alea_slice_error_page_t* page,
    const alea_transition_slice_critical_tile_t* tile,
    size_t output_limit) {
    alea_slice_error_page_t unresolved = {0};
    alea_slice_error_region_t region = {0};
    memcpy(region.uv_min, tile->uv_min, sizeof(region.uv_min));
    memcpy(region.uv_max, tile->uv_max, sizeof(region.uv_max));
    region.kind = ALEA_POINT_COVERAGE_UNRESOLVED;
    unresolved.regions = &region;
    unresolved.region_count = 1;
    return slice_error_append_products(page, &unresolved, output_limit);
}

static alea_slice_error_unresolved_reason_t
slice_error_reason_from_critical_stop(
    alea_transition_slice_critical_stop_reason_t reason) {
    switch (reason) {
    case ALEA_TRANSITION_SLICE_CRITICAL_NONE:
        return ALEA_SLICE_ERROR_UNRESOLVED_CLASSIFIER_PENDING;
    case ALEA_TRANSITION_SLICE_CRITICAL_UNSUPPORTED_CURVE:
        return ALEA_SLICE_ERROR_UNRESOLVED_PRIMITIVE;
    case ALEA_TRANSITION_SLICE_CRITICAL_UNSUPPORTED_OCCURRENCE_TRAVERSAL:
    case ALEA_TRANSITION_SLICE_CRITICAL_CHAIN_TRUNCATED:
        return ALEA_SLICE_ERROR_UNRESOLVED_OCCURRENCE;
    case ALEA_TRANSITION_SLICE_CRITICAL_NUMERICAL_UNRESOLVED:
        return ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL;
    default:
        return ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT;
    }
}

static int slice_error_classify_verified_tile(
    const alea_slice_error_query_t* query,
    const alea_transition_slice_critical_tile_t* tile,
    size_t contextual_bytes, alea_slice_error_page_t* out,
    alea_slice_error_unresolved_reason_t* reason,
    int* output_omitted, size_t* peak_scratch) {
    int verified = 0;
    *output_omitted = 0;
    *peak_scratch = 0;
    if (query->has_hierarchy &&
        slice_error_axis_view_supported(&query->options.view))
        verified = slice_error_classify_constant_occurrences(
            query, tile, contextual_bytes, out, reason,
            output_omitted, peak_scratch, NULL, NULL);
    if (query->has_hierarchy && verified == 0 && !*output_omitted &&
        slice_error_axis_view_supported(&query->options.view)) {
        size_t line_peak = 0;
        verified = slice_error_classify_single_occurrence_line(
            query, tile, contextual_bytes, out, reason,
            output_omitted, &line_peak);
        if (*peak_scratch < line_peak) *peak_scratch = line_peak;
    }
    if (verified == 0 && !*output_omitted) {
        const alea_slice_error_unresolved_reason_t occurrence_reason =
            *reason;
        size_t axis_peak = 0;
        verified = slice_error_classify_axis_tile(
            query, tile, contextual_bytes, out, reason,
            output_omitted, &axis_peak);
        if (*peak_scratch < axis_peak) *peak_scratch = axis_peak;
        if (verified == 0 &&
            *reason == ALEA_SLICE_ERROR_UNRESOLVED_OCCURRENCE &&
            (occurrence_reason ==
                 ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT ||
             occurrence_reason == ALEA_SLICE_ERROR_UNRESOLVED_NUMERICAL))
            *reason = occurrence_reason;
        else if (verified == 0 &&
                 *reason == ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT &&
                 occurrence_reason !=
                     ALEA_SLICE_ERROR_UNRESOLVED_CLASSIFIER_PENDING &&
                 occurrence_reason !=
                     ALEA_SLICE_ERROR_UNRESOLVED_CANDIDATE_LIMIT)
            *reason = occurrence_reason;
    }
    return verified;
}

int alea_slice_error_query_run_page(alea_slice_error_query_t* query,
                                    size_t page_index,
                                    alea_slice_error_page_t* page) {
    if (!query || !page || page_index >= query->page_count ||
        alea_system_geometry_generation(query->sys) !=
            query->geometry_generation) return -1;
    const double page_started = alea_monotonic_seconds();
    const alea_slice_error_query_options_t* options = &query->options;
    const size_t col = page_index % options->tile_columns;
    const size_t row = page_index / options->tile_columns;
    alea_transition_slice_critical_tile_t tile = {0};
    for (size_t axis = 0; axis < 2; ++axis) {
        const size_t index = axis == 0 ? col : row;
        const size_t count = axis == 0
            ? options->tile_columns : options->tile_rows;
        const double lo = options->required_uv_min[axis];
        const double hi = options->required_uv_max[axis];
        const double width = hi - lo;
        if (!isfinite(width)) return -1;
        tile.uv_min[axis] = index == 0 ? lo : lo + width *
            ((double)index / (double)count);
        tile.uv_max[axis] = index + 1 == count ? hi : lo + width *
            ((double)(index + 1) / (double)count);
        if (!(tile.uv_max[axis] > tile.uv_min[axis])) return -1;
    }
    if (alea_raycast_ensure_hier_caches(query->sys) != 0) return -1;
    alea_transition_slice_critical_tile_t* critical_tiles = &tile;
    size_t critical_tile_count = 1;
    alea_transition_slice_critical_tile_t* allocated_tiles = NULL;
    size_t critical_divisions = 1;
    size_t max_critical_divisions = 1;
    while (max_critical_divisions + 1u <=
               options->scan_options.max_critical_tiles /
                   (max_critical_divisions + 1u))
        max_critical_divisions++;
    if (query->has_hierarchy) {
        alea_bbox_t world_box;
        if (!slice_error_occurrence_tile_box(
                &options->view, &tile, NULL, &world_box))
            return -1;
        size_t occurrence_count = 0;
        slice_error_density_count_t density = {
            .sys = query->sys, .view = &options->view, .tile = &tile};
        alea_hier_region_chain_status_t occurrence_status =
            ALEA_HIER_REGION_CHAIN_UNSUPPORTED;
        if (alea_hier_spatial_visit_region_occurrences_bounded(
                query->sys, &world_box, slice_error_count_occurrence, &density,
                options->scan_options.max_exhaustive_occurrence_hits,
                &occurrence_count, &occurrence_status) != 0)
            return -1;
        if (density.failed) return -1;
        if (occurrence_status == ALEA_HIER_REGION_CHAIN_COMPLETE) {
            size_t target = options->scan_options.max_curves_per_tile / 8u;
            if (target < 64u) target = 64u;
            const size_t density_count = density.candidate_count >
                    occurrence_count
                ? density.candidate_count : occurrence_count;
            const size_t tile_goal = density_count > target
                ? (density_count + target - 1u) / target : 1u;
            while (critical_divisions < max_critical_divisions &&
                   critical_divisions * critical_divisions < tile_goal)
                critical_divisions++;
            if (critical_divisions > 1) {
                critical_tile_count =
                    critical_divisions * critical_divisions;
                allocated_tiles = calloc(
                    critical_tile_count, sizeof(*allocated_tiles));
                if (!allocated_tiles) return -1;
                slice_error_uniform_critical_tiles(
                    &tile, critical_divisions, allocated_tiles);
                critical_tiles = allocated_tiles;
            }
        }
    }
    alea_transition_slice_stats_t stats = {0};
    stats.critical_stop_reason = ALEA_TRANSITION_SLICE_CRITICAL_NONE;
    alea_slice_error_page_t candidate = {0};
    slice_error_finding_sink_t sink = {
        .page = &candidate,
        .max_findings = options->scan_options.max_critical_findings,
        .max_output_bytes = options->scan_options.max_output_bytes
    };
    alea_slice_view_t analysis_view = options->view;
    analysis_view.u_min = options->required_uv_min[0];
    analysis_view.u_max = options->required_uv_max[0];
    analysis_view.v_min = options->required_uv_min[1];
    analysis_view.v_max = options->required_uv_max[1];
    int critical_rc = 0;
    alea_transition_slice_critical_stop_reason_t* tile_stop_reasons = NULL;
    uint8_t* localized_tiles = NULL;
    slice_error_numerical_region_sink_t numerical_sink = {0};
    for (;;) {
        tile_stop_reasons = calloc(
            critical_tile_count, sizeof(*tile_stop_reasons));
        localized_tiles = calloc(critical_tile_count, sizeof(*localized_tiles));
        if (!tile_stop_reasons || !localized_tiles) {
            free(tile_stop_reasons);
            free(localized_tiles);
            free(allocated_tiles);
            free(candidate.findings);
            free(candidate.regions);
            return -1;
        }
        numerical_sink = (slice_error_numerical_region_sink_t){
            .page = &candidate,
            .localized_tiles = localized_tiles,
            .tile_count = critical_tile_count,
            .max_output_bytes = options->scan_options.max_output_bytes
        };
        critical_rc = alea_transition_slice_enumerate_critical_tiles_with_regions(
            query->sys, &analysis_view, &options->scan_options,
            critical_tiles, critical_tile_count,
            slice_error_retain_finding, &sink, &stats,
            tile_stop_reasons, critical_tile_count,
            slice_error_retain_numerical_region, &numerical_sink);
        if (critical_rc != 0 ||
            stats.critical_stop_reason !=
                ALEA_TRANSITION_SLICE_CRITICAL_MAX_CURVES ||
            critical_divisions >= max_critical_divisions)
            break;

        free(candidate.findings);
        free(candidate.regions);
        memset(&candidate, 0, sizeof(candidate));
        sink.page = &candidate;
        sink.omitted_findings = 0;
        memset(&stats, 0, sizeof(stats));
        stats.critical_stop_reason = ALEA_TRANSITION_SLICE_CRITICAL_NONE;
        free(tile_stop_reasons);
        tile_stop_reasons = NULL;
        free(localized_tiles);
        localized_tiles = NULL;
        free(allocated_tiles);
        allocated_tiles = NULL;
        size_t next = critical_divisions <= max_critical_divisions / 2u
            ? critical_divisions * 2u : max_critical_divisions;
        if (next <= critical_divisions) break;
        critical_divisions = next;
        critical_tile_count = critical_divisions * critical_divisions;
        allocated_tiles = calloc(
            critical_tile_count, sizeof(*allocated_tiles));
        if (!allocated_tiles) return -1;
        slice_error_uniform_critical_tiles(
            &tile, critical_divisions, allocated_tiles);
        critical_tiles = allocated_tiles;
    }
    if (critical_rc != 0 ||
        alea_system_geometry_generation(query->sys) !=
            query->geometry_generation) {
        free(tile_stop_reasons);
        free(localized_tiles);
        free(allocated_tiles);
        free(candidate.findings);
        free(candidate.regions);
        return -1;
    }

    slice_error_confirmation_t confirmation = {0};
    (void)slice_error_confirmation_init(
        &confirmation, options->scan_options.max_coverage_hits);
    const double confirmation_started = alea_monotonic_seconds();
    const size_t interior_axis_count =
        options->scan_options.coverage_uniform_probes_per_ray;
    size_t interior_probe_count = 0;
    for (size_t row = 0; row < interior_axis_count; ++row) {
        for (size_t column = 0; column < interior_axis_count; ++column) {
            double uv[2] = {
                tile.uv_min[0] + (tile.uv_max[0] - tile.uv_min[0]) *
                    ((double)column + 0.5) / (double)interior_axis_count,
                tile.uv_min[1] + (tile.uv_max[1] - tile.uv_min[1]) *
                    ((double)row + 0.5) / (double)interior_axis_count};
            double world[3];
            for (int axis = 0; axis < 3; ++axis)
                world[axis] = options->view.plane.origin[axis] +
                    uv[0] * options->view.plane.u_axis[axis] +
                    uv[1] * options->view.plane.v_axis[axis];
            (void)slice_error_confirm_point(
                query, &candidate, &confirmation, uv, world,
                ALEA_SLICE_ERROR_WITNESS_INTERIOR_PROBE, NULL);
            interior_probe_count++;
        }
    }
    for (size_t i = 0; i < candidate.finding_count; ++i) {
        const alea_transition_slice_critical_finding_t* finding =
            &candidate.findings[i];
        if (finding->transition.kind != ALEA_TRANSITION_GAP &&
            finding->transition.kind != ALEA_TRANSITION_OVERLAP) continue;
        double uv[2] = {finding->uv[0], finding->uv[1]};
        double world[3];
        memcpy(world, finding->world_point, sizeof(world));
        int displaced = 0;
        for (int axis = 0; axis < 3; ++axis)
            displaced |= finding->transition.before_point[axis] !=
                         finding->transition.after_point[axis];
        if (displaced) {
            const double distance = finding->transition.probe_distance;
            for (int axis = 0; axis < 3; ++axis)
                world[axis] += distance * finding->direction[axis];
            for (int axis = 0; axis < 3; ++axis) {
                uv[0] += distance * finding->direction[axis] *
                         options->view.plane.u_axis[axis];
                uv[1] += distance * finding->direction[axis] *
                         options->view.plane.v_axis[axis];
            }
        }
        if (uv[0] < tile.uv_min[0] || uv[0] > tile.uv_max[0] ||
            uv[1] < tile.uv_min[1] || uv[1] > tile.uv_max[1]) continue;
        (void)slice_error_confirm_point(
            query, &candidate, &confirmation, uv, world,
            ALEA_SLICE_ERROR_WITNESS_BOUNDARY_PROBE, finding);
    }
    const double confirmation_seconds =
        alea_monotonic_seconds() - confirmation_started;
    alea_slice_error_page_receipt_t receipt = {0};
    receipt.query_id = query->query_id;
    receipt.geometry_generation = query->geometry_generation;
    receipt.boundary_analysis_policy_version =
        ALEA_SLICE_ERROR_BOUNDARY_ANALYSIS_POLICY_VERSION;
    receipt.close_crossing_relative_tolerance = query->options.scan_options
        .critical_relative_distance_tolerance;
    receipt.page_index = page_index;
    memcpy(receipt.core_uv_min, tile.uv_min, sizeof(tile.uv_min));
    memcpy(receipt.core_uv_max, tile.uv_max, sizeof(tile.uv_max));
    receipt.occurrence_paths = stats.critical_occurrence_paths;
    receipt.candidate_curves = stats.critical_curves;
    receipt.candidate_pairs_tested = stats.critical_curve_pairs_tested;
    receipt.peak_scratch_bytes = stats.peak_critical_scratch_bytes;
    receipt.close_crossing_observations =
        stats.critical_close_crossing_observations;
    receipt.symbolic_one_sided_intervals =
        stats.critical_symbolic_one_sided_intervals;
    receipt.numerical_unsafe_probe_intervals =
        stats.critical_numerical_unsafe_probe_intervals;
    receipt.numerical_unrepresentable_probe_intervals =
        stats.critical_numerical_unrepresentable_probe_intervals;
    receipt.numerical_inconsistent_probe_intervals =
        stats.critical_numerical_inconsistent_probe_intervals;
    receipt.query_index_bytes = query->index_bytes;
    receipt.contextual_finding_count = candidate.finding_count;
    receipt.omitted_contextual_findings = sink.omitted_findings;
    receipt.omitted_contextual_boundary_evidence =
        stats.omitted_critical_boundary_evidence;
    receipt.interior_probe_count = interior_probe_count;
    receipt.confirmation_attempt_count = confirmation.attempts;
    receipt.confirmation_failure_count = confirmation.failures;
    receipt.confirmed_gap_witness_count = confirmation.confirmed_gaps;
    receipt.confirmed_overlap_witness_count =
        confirmation.confirmed_overlaps;
    receipt.omitted_confirmed_witnesses = confirmation.omitted;
    receipt.confirmation_seconds = confirmation_seconds;
    receipt.scan_stop_reason = stats.critical_stop_reason;
    receipt.unresolved_reason = slice_error_reason_from_critical_stop(
        stats.critical_stop_reason);
    receipt.requested_work_complete = stats.critical_stop_reason ==
        ALEA_TRANSITION_SLICE_CRITICAL_NONE && !sink.omitted_findings;
    receipt.scope_classified = 0;
    receipt.output_complete = sink.omitted_findings == 0 &&
        stats.omitted_critical_boundary_evidence == 0 &&
        numerical_sink.omitted_regions == 0 && confirmation.omitted == 0;
    if (receipt.requested_work_complete) {
        size_t classified_tiles = 0;
        for (size_t ti = 0; ti < critical_tile_count; ++ti) {
            alea_slice_error_page_t classified = {0};
            alea_slice_error_unresolved_reason_t tile_reason =
                ALEA_SLICE_ERROR_UNRESOLVED_CLASSIFIER_PENDING;
            int output_omitted = 0;
            size_t classifier_peak = 0;
            const int verified = slice_error_classify_verified_tile(
                query, &critical_tiles[ti],
                slice_error_page_output_bytes(&candidate), &classified,
                &tile_reason, &output_omitted, &classifier_peak);
            if (receipt.peak_scratch_bytes < classifier_peak)
                receipt.peak_scratch_bytes = classifier_peak;
            if (verified < 0) {
                free(classified.intervals);
                free(classified.circles);
                free(classified.regions);
                goto page_failed;
            }
            if (verified > 0) {
                const int appended = slice_error_append_products(
                    &candidate, &classified,
                    options->scan_options.max_output_bytes);
                if (appended < 0) {
                    free(classified.intervals);
                    free(classified.circles);
                    free(classified.regions);
                    goto page_failed;
                }
                if (appended == 0) receipt.output_complete = 0;
                classified_tiles++;
            } else {
                const int appended = slice_error_append_unresolved_region(
                    &candidate, &critical_tiles[ti],
                    options->scan_options.max_output_bytes);
                if (appended < 0) {
                    free(classified.intervals);
                    free(classified.circles);
                    free(classified.regions);
                    goto page_failed;
                }
                if (appended == 0 || output_omitted)
                    receipt.output_complete = 0;
                if (receipt.unresolved_reason ==
                        ALEA_SLICE_ERROR_UNRESOLVED_CLASSIFIER_PENDING)
                    receipt.unresolved_reason = tile_reason;
            }
            free(classified.intervals);
            free(classified.circles);
            free(classified.regions);
        }
        receipt.scope_classified = classified_tiles == critical_tile_count;
        if (receipt.scope_classified)
            receipt.unresolved_reason = ALEA_SLICE_ERROR_RESOLVED;
    } else {
        size_t classified_tiles = 0;
        for (size_t ti = 0; ti < critical_tile_count; ++ti) {
            if (tile_stop_reasons[ti] !=
                    ALEA_TRANSITION_SLICE_CRITICAL_NONE) {
                if (tile_stop_reasons[ti] !=
                        ALEA_TRANSITION_SLICE_CRITICAL_NUMERICAL_UNRESOLVED ||
                    !localized_tiles[ti]) {
                    const int appended = slice_error_append_unresolved_region(
                        &candidate, &critical_tiles[ti],
                        options->scan_options.max_output_bytes);
                    if (appended < 0) goto page_failed;
                    if (appended == 0) receipt.output_complete = 0;
                }
                continue;
            }
            alea_slice_error_page_t classified = {0};
            alea_slice_error_unresolved_reason_t tile_reason =
                ALEA_SLICE_ERROR_UNRESOLVED_CLASSIFIER_PENDING;
            int output_omitted = 0;
            size_t classifier_peak = 0;
            const int verified = slice_error_classify_verified_tile(
                query, &critical_tiles[ti],
                slice_error_page_output_bytes(&candidate), &classified,
                &tile_reason, &output_omitted, &classifier_peak);
            if (receipt.peak_scratch_bytes < classifier_peak)
                receipt.peak_scratch_bytes = classifier_peak;
            if (verified < 0) {
                free(classified.intervals);
                free(classified.circles);
                free(classified.regions);
                goto page_failed;
            }
            if (verified > 0) {
                const int appended = slice_error_append_products(
                    &candidate, &classified,
                    options->scan_options.max_output_bytes);
                if (appended < 0) {
                    free(classified.intervals);
                    free(classified.circles);
                    free(classified.regions);
                    goto page_failed;
                }
                if (appended == 0) receipt.output_complete = 0;
                else classified_tiles++;
            } else {
                const int appended = slice_error_append_unresolved_region(
                    &candidate, &critical_tiles[ti],
                    options->scan_options.max_output_bytes);
                if (appended < 0) {
                    free(classified.intervals);
                    free(classified.circles);
                    free(classified.regions);
                    goto page_failed;
                }
                if (appended == 0 || output_omitted)
                    receipt.output_complete = 0;
                if (receipt.unresolved_reason ==
                        ALEA_SLICE_ERROR_UNRESOLVED_CLASSIFIER_PENDING)
                    receipt.unresolved_reason = tile_reason;
            }
            free(classified.intervals);
            free(classified.circles);
            free(classified.regions);
        }
        /* A partial page intentionally remains unclassified as a whole. The
         * retained verified products and UNRESOLVED rectangles identify what
         * is actionable and what still needs refinement. */
        (void)classified_tiles;
    }
    free(tile_stop_reasons);
    free(localized_tiles);
    free(allocated_tiles);
    receipt.verified_interval_count = candidate.interval_count;
    receipt.verified_circle_count = candidate.circle_count;
    receipt.region_count = candidate.region_count;
    receipt.elapsed_seconds = alea_monotonic_seconds() - page_started;
    candidate.receipt = receipt;
    candidate.populated = 1;
    free(page->findings);
    free(page->witnesses);
    free(page->intervals);
    free(page->circles);
    free(page->regions);
    *page = candidate;
    slice_error_confirmation_free(&confirmation);
    return 0;

page_failed:
    slice_error_confirmation_free(&confirmation);
    free(tile_stop_reasons);
    free(localized_tiles);
    free(allocated_tiles);
    free(candidate.findings);
    free(candidate.witnesses);
    free(candidate.intervals);
    free(candidate.circles);
    free(candidate.regions);
    return -1;
}

typedef struct {
    alea_slice_error_query_t* query;
    const size_t* page_indices;
    alea_slice_error_page_t* const* pages;
    int* statuses;
} slice_error_batch_context_t;

static int slice_error_batch_run(
    void* opaque, size_t worker_index, size_t begin, size_t end) {
    (void)worker_index;
    slice_error_batch_context_t* context = opaque;
    for (size_t i = begin; i < end; ++i)
        context->statuses[i] = alea_slice_error_query_run_page(
            context->query, context->page_indices[i], context->pages[i]);
    return 0;
}

int alea_slice_error_query_run_pages(
    alea_slice_error_query_t* query, const size_t* page_indices,
    size_t page_count, size_t requested_workers,
    uint64_t max_parallel_scratch_bytes,
    alea_slice_error_page_t* const* pages,
    alea_transition_slice_batch_stats_t* out_stats) {
    if (!query || (!page_indices && page_count) || (!pages && page_count) ||
        !out_stats) return -1;
    memset(out_stats, 0, sizeof(*out_stats));
    out_stats->page_count = page_count;
    out_stats->requested_workers = requested_workers;
    if (!page_count) return 0;
    for (size_t i = 0; i < page_count; ++i) {
        if (!pages[i] || page_indices[i] >= query->page_count) return -1;
        for (size_t j = 0; j < i; ++j)
            if (pages[i] == pages[j]) return -1;
    }
    if (alea_system_geometry_generation(query->sys) !=
            query->geometry_generation ||
        alea_raycast_ensure_hier_caches(query->sys) != 0)
        return -1;

    uint64_t worker_scratch =
        query->options.scan_options.max_critical_scratch_bytes;
    if (!worker_scratch) worker_scratch = 1;
    out_stats->reserved_scratch_bytes_per_worker = worker_scratch;
    size_t workers = alea_parallel_effective_workers(
        page_count, 1, requested_workers);
    if (!max_parallel_scratch_bytes) {
        workers = 1;
    } else {
        uint64_t budget_workers =
            max_parallel_scratch_bytes / worker_scratch;
        if (!budget_workers) budget_workers = 1;
        if ((uint64_t)workers > budget_workers)
            workers = (size_t)budget_workers;
    }
    if (!workers) workers = 1;
    if (workers > UINT64_MAX / worker_scratch) return -1;
    out_stats->reserved_parallel_scratch_bytes =
        (uint64_t)workers * worker_scratch;

    int* statuses = calloc(page_count, sizeof(*statuses));
    if (!statuses) return -1;
    slice_error_batch_context_t context = {
        .query = query, .page_indices = page_indices,
        .pages = pages, .statuses = statuses};
    size_t actual_workers = 1;
    const alea_parallel_status_t parallel_status = alea_parallel_for(
        page_count, 1, workers, ALEA_PARALLEL_STATIC_BLOCK,
        slice_error_batch_run, &context, &actual_workers);
    if (parallel_status != ALEA_PARALLEL_OK) {
        free(statuses);
        return -1;
    }
    out_stats->actual_workers = actual_workers;
    int failed = 0;
    for (size_t i = 0; i < page_count; ++i) {
        if (statuses[i] == 0) out_stats->completed_page_count++;
        else failed = 1;
    }
    free(statuses);
    return failed ? -1 : 0;
}

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
    options->max_curves_per_tile = 2048;
    options->max_critical_points = 2048;
    options->max_active_boundary_tests = 100000;
    options->max_critical_probes = 4096;
    options->max_critical_findings = 1024;
    options->max_critical_boundary_evidence = 1024;
    options->max_curve_pairs = 100000;
    options->max_critical_sector_witnesses = 8192;
    options->critical_relative_distance_tolerance = 1e-6;
    options->occurrence_discovery = ALEA_TRANSITION_SLICE_OCCURRENCE_SAMPLED;
    /* Occurrence receipts are streamed, so this bounds traversal work without
     * increasing the fixed critical scratch allocation. */
    options->max_exhaustive_occurrence_hits = 65536;
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
    case ALEA_TRANSITION_SLICE_CRITICAL_MAX_OCCURRENCE_HITS:
        return "max_exhaustive_occurrence_hits";
    case ALEA_TRANSITION_SLICE_CRITICAL_UNSUPPORTED_OCCURRENCE_TRAVERSAL:
        return "unsupported_occurrence_traversal";
    case ALEA_TRANSITION_SLICE_CRITICAL_NUMERICAL_UNRESOLVED:
        return "numerical_unresolved";
    }
    return "unknown";
}

const char* alea_slice_error_numerical_cause_name(
    alea_slice_error_numerical_cause_t cause) {
    switch (cause) {
    case ALEA_SLICE_ERROR_NUMERICAL_NONE: return "none";
    case ALEA_SLICE_ERROR_NUMERICAL_CLOSE_BREAKPOINTS:
        return "close_breakpoints";
    case ALEA_SLICE_ERROR_NUMERICAL_INVALID_PROBE_RADIUS:
        return "invalid_probe_radius";
    case ALEA_SLICE_ERROR_NUMERICAL_UNREPRESENTABLE_PROBES:
        return "unrepresentable_probes";
    case ALEA_SLICE_ERROR_NUMERICAL_INCONSISTENT_CLASSIFICATION:
        return "inconsistent_classification";
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

static int transition_slice_scan_ray_reuse(
    alea_system_t* sys, const alea_slice_view_t* view,
    const alea_transition_slice_options_t* options,
    alea_transition_slice_orientation_t orientation, size_t ray_index,
    size_t base_ray_index, uint32_t refinement_depth,
    double transverse_coordinate, alea_raycast_result_t* scratch,
    alea_ray_boundary_event_result_t* events,
    transition_slice_coverage_scratch_t* coverage_scratch,
    alea_transition_workspace_t* transition_workspace,
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
            if (alea_check_selected_boundary_event_transition_reuse_nocache(
                    sys, event, &transition_options, &transition,
                    transition_workspace) != 0)
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
    if (options->critical_full_view) {
        if (seed_count == SIZE_MAX) return -1;
        seed_count++;
    }
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
    if (options->critical_full_view) {
        transition_slice_tile_seed_t* item = &seeds[seed++];
        item->uv_min[0] = view->u_min;
        item->uv_min[1] = view->v_min;
        item->uv_max[0] = view->u_max;
        item->uv_max[1] = view->v_max;
        item->kind = ALEA_TRANSITION_SLICE_TILE_SOURCE_FULL_VIEW;
        item->priority = -1;
    }
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
    alea_transition_workspace_t* transition_workspace,
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
        int rc = transition_slice_scan_ray_reuse(
            sys, view, options, orientation, ray_index, base, 0, coordinate,
            scratch, events, coverage_scratch, transition_workspace, result,
            &rows[base]);
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
            int rc = transition_slice_scan_ray_reuse(
                sys, view, options, orientation, ray_index, SIZE_MAX,
                depth + 1, coordinate, scratch, events, coverage_scratch,
                transition_workspace, result, &next[output]);
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
        !isfinite(options.critical_relative_distance_tolerance) ||
        options.critical_relative_distance_tolerance < 0.0 ||
        options.occurrence_discovery <
            ALEA_TRANSITION_SLICE_OCCURRENCE_SAMPLED ||
        options.occurrence_discovery >
            ALEA_TRANSITION_SLICE_OCCURRENCE_EXHAUSTIVE ||
        (options.occurrence_discovery ==
             ALEA_TRANSITION_SLICE_OCCURRENCE_EXHAUSTIVE &&
         options.max_exhaustive_occurrence_hits == 0) ||
        (options.refine_signals &
         ~(ALEA_TRANSITION_SLICE_REFINE_SIGNATURE |
           ALEA_TRANSITION_SLICE_REFINE_FINDING)))
        return -1;

    alea_transition_slice_result_t candidate = {0};
    candidate.stats.occurrence_discovery = options.occurrence_discovery;
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
        alea_transition_workspace_t transition_workspace;
        memset(&coverage_scratch, 0, sizeof(coverage_scratch));
        const int coverage_enabled =
            options.coverage_uniform_probes_per_ray != 0 ||
            options.coverage_probe_selected_intervals;
        size_t coverage_scratch_bytes = coverage_enabled
            ? transition_slice_coverage_scratch_bytes(
                options.max_coverage_hits) : 0;
        alea_raycast_result_init(&scratch);
        alea_ray_boundary_event_result_init(&events);
        alea_transition_workspace_init(&transition_workspace);
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
                &scratch, &events, &coverage_scratch, &transition_workspace,
                coverage_scratch_bytes, &candidate);
            if (rc < 0) failed = 1;
            if (rc > 0) { stop = 1; scan_incomplete = 1; }
            if (failed) break;
        }
        alea_ray_boundary_event_result_free(&events);
        alea_raycast_result_free(&scratch);
        transition_slice_coverage_scratch_free(&coverage_scratch);
        alea_transition_workspace_free(&transition_workspace);
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
    candidate.stats.occurrence_enumeration_complete =
        options.occurrence_discovery ==
            ALEA_TRANSITION_SLICE_OCCURRENCE_EXHAUSTIVE &&
        candidate.stats.critical_complete;
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

typedef struct {
    alea_system_t* sys;
    const alea_slice_view_t* views;
    const alea_transition_slice_options_t* options;
    alea_transition_slice_result_t* const* results;
    int* statuses;
} alea_transition_slice_batch_context_t;

static int alea_transition_slice_batch_run(
    void* opaque, size_t worker_index, size_t begin, size_t end) {
    (void)worker_index;
    alea_transition_slice_batch_context_t* context = opaque;
    for (size_t page = begin; page < end; page++) {
        context->statuses[page] = alea_transition_slice_screen(
            context->sys, &context->views[page], context->options,
            context->results[page]);
    }
    return 0;
}

int alea_transition_slice_screen_batch(
    alea_system_t* sys, const alea_slice_view_t* views, size_t page_count,
    const alea_transition_slice_options_t* input,
    size_t requested_workers, uint64_t max_parallel_scratch_bytes,
    alea_transition_slice_result_t* const* results,
    alea_transition_slice_batch_stats_t* out_stats) {
    if (!sys || (!views && page_count) || (!results && page_count) ||
            !out_stats) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "transition slice batch argument is null");
        return -1;
    }
    memset(out_stats, 0, sizeof(*out_stats));
    out_stats->page_count = page_count;
    out_stats->requested_workers = requested_workers;
    if (page_count == 0) return 0;
    for (size_t page = 0; page < page_count; page++) {
        if (!results[page]) {
            alea_set_error_detail(ALEA_ERR_NULL_ARG,
                                  "transition slice batch result is null");
            return -1;
        }
    }

    alea_transition_slice_options_t defaults, options;
    alea_transition_slice_options_init(&defaults);
    options = defaults;
    if (input) {
        if (input->struct_size < sizeof(input->struct_size)) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "transition slice options are too small");
            return -1;
        }
        size_t bytes = input->struct_size < sizeof(options)
            ? input->struct_size : sizeof(options);
        memcpy(&options, input, bytes);
    }
    uint64_t worker_scratch = options.max_scratch_bytes;
    if (worker_scratch > UINT64_MAX - options.max_critical_scratch_bytes) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "transition slice worker scratch overflows");
        return -1;
    }
    worker_scratch += options.max_critical_scratch_bytes;
    if (worker_scratch == 0) worker_scratch = 1;
    out_stats->reserved_scratch_bytes_per_worker = worker_scratch;

    size_t workers = alea_parallel_effective_workers(page_count, 1,
                                                     requested_workers);
    if (max_parallel_scratch_bytes == 0) {
        workers = 1;
    } else {
        uint64_t budget_workers = max_parallel_scratch_bytes / worker_scratch;
        if (budget_workers == 0) budget_workers = 1;
        if ((uint64_t)workers > budget_workers)
            workers = (size_t)budget_workers;
    }
    if (workers > UINT64_MAX / worker_scratch) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "transition slice parallel scratch overflows");
        return -1;
    }
    out_stats->reserved_parallel_scratch_bytes =
        (uint64_t)workers * worker_scratch;

    if (alea_raycast_ensure_hier_caches(sys) != 0) return -1;
    int* statuses = calloc(page_count, sizeof(*statuses));
    if (!statuses) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "transition slice batch status allocation failed");
        return -1;
    }
    size_t actual_workers = 1;
    alea_transition_slice_batch_context_t context = {
        sys, views, &options, results, statuses
    };
    alea_parallel_status_t parallel_status = alea_parallel_for(
        page_count, 1, workers, ALEA_PARALLEL_STATIC_BLOCK,
        alea_transition_slice_batch_run, &context, &actual_workers);
    if (parallel_status != ALEA_PARALLEL_OK) {
        free(statuses);
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              alea_parallel_status_string(parallel_status));
        return -1;
    }
    out_stats->actual_workers = actual_workers;
    int failed = 0;
    for (size_t page = 0; page < page_count; page++) {
        if (statuses[page] == 0) out_stats->completed_page_count++;
        else failed = 1;
    }
    free(statuses);
    return failed ? -1 : 0;
}
