// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "geo_validator.h"
#include "raycast/ray_bbox.h"

#include "alea_slice.h"
#include "core/alea_cell.h"
#include "core/alea_system.h"
#include "core/alea_universe.h"
#include "primitives/bbox.h"
#include "raycast/ray_epsilon.h"
#include "raycast/raycast.h"
#include "util/math.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define VALIDATOR_HIT_CAP 64
#define VALIDATOR_DEFAULT_RAYS 128
#define VALIDATOR_DEFAULT_MAX_ERRORS 1024
#define VALIDATOR_DEFAULT_MAX_CROSSINGS 100000

typedef enum {
    COVERAGE_NONE = 0,
    COVERAGE_ONE = 1,
    COVERAGE_MULTI = 2
} coverage_class_t;

typedef struct {
    coverage_class_t klass;
    int primary_cell_id;
    int primary_cell_idx;
    uint64_t primary_occurrence_key;
    /* Root-to-primary chain for the selected concrete occurrence. Adjacency is
     * defined on the crossed primitive's owning cell, which may be an
     * enclosing fill container while point coverage projects to a child. */
    int primary_ancestor_cell_indices[VALIDATOR_HIT_CAP];
    size_t primary_ancestor_count;
    int secondary_cell_id;
    uint64_t secondary_occurrence_key;
    int universe_id;
    int depth;
    int count_at_depth;
    int target_depth;
    int truncated;
} point_coverage_t;

static int cell_references_primitive(const alea_system_t* sys,
                                     const alea_cell_entry_t* cell,
                                     uint32_t primitive_id);

static int point_coverage_primary_has_cell_index(const point_coverage_t* cov,
                                                  int cell_index) {
    if (!cov || cell_index < 0) return 0;
    for (size_t index = 0; index < cov->primary_ancestor_count; index++)
        if (cov->primary_ancestor_cell_indices[index] == cell_index)
            return 1;
    return 0;
}

static int point_coverage_primary_ancestor_references_primitive(
    const alea_system_t* sys, const point_coverage_t* cov,
    uint32_t primitive_id) {
    if (!sys || !cov || primitive_id == ALEA_PRIMITIVE_ID_INVALID) return 0;
    for (size_t index = 0; index < cov->primary_ancestor_count; index++) {
        int cell_index = cov->primary_ancestor_cell_indices[index];
        if (cell_index < 0 ||
            (size_t)cell_index >= alea_vec_count(&sys->cells)) continue;
        if (cell_references_primitive(sys, &sys->cells.data[cell_index],
                                      primitive_id))
            return 1;
    }
    return 0;
}

static void point_coverage_set_primary_ancestors(
    point_coverage_t* out, const alea_cell_hit_t* hits,
    const uint64_t* occurrence_keys, const uint64_t* parent_occurrence_keys,
    size_t hit_count) {
    if (!out || !hits || !occurrence_keys || !parent_occurrence_keys) return;
    out->primary_ancestor_count = 0;
    uint64_t key = out->primary_occurrence_key;
    while (key != 0 && out->primary_ancestor_count < VALIDATOR_HIT_CAP) {
        size_t index = 0;
        while (index < hit_count && occurrence_keys[index] != key)
            index++;
        if (index == hit_count) break;
        out->primary_ancestor_cell_indices[out->primary_ancestor_count++] =
            hits[index].cell_index;
        key = parent_occurrence_keys[index];
    }
}

static void copy3(double dst[3], const double src[3]) {
    dst[0] = src[0];
    dst[1] = src[1];
    dst[2] = src[2];
}

static void init_geom_error(alea_geom_error_t* err) {
    memset(err, 0, sizeof(*err));
    err->source = ALEA_GEOM_EVENT_SOURCE_UNKNOWN;
    err->curve_index = SIZE_MAX;
    err->primitive_id = ALEA_PRIMITIVE_ID_INVALID;
}

typedef struct {
    int type, source, previous_cell_id, found_cell_id, secondary_cell_id;
    int expected_neighbor_cell_id, universe_depth;
    uint32_t primitive_id;
} error_signature_t;

typedef struct {
    error_signature_t key;
    size_t count;
    int used;
} error_signature_entry_t;

typedef struct {
    error_signature_entry_t* entries;
    size_t capacity, count;
} error_signature_table_t;

static uint64_t signature_hash_mix(uint64_t hash, uint64_t value) {
    hash ^= value;
    return hash * UINT64_C(1099511628211);
}

static error_signature_t error_signature(const alea_geom_error_t* error) {
    return (error_signature_t){
        .type = (int)error->type, .source = (int)error->source,
        .previous_cell_id = error->previous_cell_id,
        .found_cell_id = error->found_cell_id,
        .secondary_cell_id = error->secondary_cell_id,
        .expected_neighbor_cell_id = error->expected_neighbor_cell_id,
        .universe_depth = error->universe_depth,
        .primitive_id = error->primitive_id,
    };
}

static int signatures_equal(error_signature_t a, error_signature_t b) {
    return a.type == b.type && a.source == b.source &&
        a.previous_cell_id == b.previous_cell_id &&
        a.found_cell_id == b.found_cell_id &&
        a.secondary_cell_id == b.secondary_cell_id &&
        a.expected_neighbor_cell_id == b.expected_neighbor_cell_id &&
        a.universe_depth == b.universe_depth && a.primitive_id == b.primitive_id;
}

static uint64_t signature_hash(error_signature_t key) {
    uint64_t hash = UINT64_C(1469598103934665603);
    hash = signature_hash_mix(hash, (uint32_t)key.type);
    hash = signature_hash_mix(hash, (uint32_t)key.source);
    hash = signature_hash_mix(hash, (uint32_t)key.previous_cell_id);
    hash = signature_hash_mix(hash, (uint32_t)key.found_cell_id);
    hash = signature_hash_mix(hash, (uint32_t)key.secondary_cell_id);
    hash = signature_hash_mix(hash, (uint32_t)key.expected_neighbor_cell_id);
    hash = signature_hash_mix(hash, (uint32_t)key.universe_depth);
    return signature_hash_mix(hash, key.primitive_id);
}

static int signature_table_grow(error_signature_table_t* table) {
    size_t capacity = table->capacity ? table->capacity * 2 : 16;
    if (capacity < table->capacity ||
        capacity > SIZE_MAX / sizeof(*table->entries)) return -1;
    error_signature_entry_t* entries = calloc(capacity, sizeof(*entries));
    if (!entries) return -1;
    for (size_t i = 0; i < table->capacity; i++) {
        error_signature_entry_t entry = table->entries[i];
        if (!entry.used) continue;
        size_t index = (size_t)signature_hash(entry.key) & (capacity - 1);
        while (entries[index].used) index = (index + 1) & (capacity - 1);
        entries[index] = entry;
    }
    free(table->entries);
    table->entries = entries;
    table->capacity = capacity;
    return 0;
}

static error_signature_entry_t* signature_table_entry(
    alea_geom_validator_result_t* result, const alea_geom_error_t* error) {
    error_signature_table_t* table = result->signature_table;
    if (!table) {
        table = calloc(1, sizeof(*table));
        if (!table) return NULL;
        result->signature_table = table;
    }
    if (table->capacity == 0 || table->count * 10 >= table->capacity * 7) {
        if (signature_table_grow(table) != 0) return NULL;
    }
    error_signature_t key = error_signature(error);
    size_t index = (size_t)signature_hash(key) & (table->capacity - 1);
    while (table->entries[index].used &&
           !signatures_equal(table->entries[index].key, key))
        index = (index + 1) & (table->capacity - 1);
    if (!table->entries[index].used) {
        table->entries[index].used = 1;
        table->entries[index].key = key;
        table->count++;
    }
    return &table->entries[index];
}

void alea_geom_validator_options_init(alea_geom_validator_options_t* options) {
    if (!options) return;
    memset(options, 0, sizeof(*options));
    options->flags = ALEA_GEOM_VALIDATE_RAYS |
                     ALEA_GEOM_VALIDATE_STRICT_ADJACENCY;
    options->universe_depth = -1;
    options->max_errors = VALIDATOR_DEFAULT_MAX_ERRORS;
    options->max_samples_per_signature = 4;
    options->max_crossings = VALIDATOR_DEFAULT_MAX_CROSSINGS;
    options->sample_offset = SURFACE_SAMPLE_OFFSET;
    options->t_max = 0.0;
    options->seed = 42u;
    options->ray_count = VALIDATOR_DEFAULT_RAYS;
}

void alea_geom_validator_result_init(alea_geom_validator_result_t* result) {
    if (!result) return;
    memset(result, 0, sizeof(*result));
}

void alea_geom_validator_result_free(alea_geom_validator_result_t* result) {
    if (!result) return;
    free(result->errors);
    if (result->signature_table) {
        error_signature_table_t* table = result->signature_table;
        free(table->entries);
        free(table);
    }
    memset(result, 0, sizeof(*result));
}

const char* alea_geom_error_type_name(alea_geom_error_type_t type) {
    switch (type) {
        case ALEA_GEOM_ERR_UNDEFINED_AFTER_CROSSING:
            return "undefined_after_crossing";
        case ALEA_GEOM_ERR_OVERLAP_AFTER_CROSSING:
            return "overlap_after_crossing";
        case ALEA_GEOM_ERR_NON_ADJACENT_TRANSITION:
            return "non_adjacent_transition";
        case ALEA_GEOM_ERR_MISSING_NEIGHBOR:
            return "missing_neighbor";
        case ALEA_GEOM_ERR_AMBIGUOUS_BOUNDARY:
            return "ambiguous_boundary";
        case ALEA_GEOM_ERR_INTERIOR_GAP:
            return "interior_gap";
        default:
            return "unknown";
    }
}

size_t alea_geom_validator_error_count(const alea_geom_validator_result_t* result) {
    return result ? result->error_count : 0;
}

int alea_geom_validator_error_get(const alea_geom_validator_result_t* result,
                                  size_t index,
                                  alea_geom_error_t* out_error) {
    if (!result || !out_error || index >= result->error_count)
        return -1;
    *out_error = result->errors[index];
    return 0;
}

static int append_error(alea_geom_validator_result_t* result,
                        const alea_geom_validator_options_t* options,
                        const alea_geom_error_t* error) {
    error_signature_entry_t* signature_entry = NULL;
    if (options->max_samples_per_signature > 0) {
        signature_entry = signature_table_entry(result, error);
        if (!signature_entry) {
            alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                                  "geometry validator: out of memory tracking signatures");
            return -1;
        }
        if (signature_entry->count >= options->max_samples_per_signature) {
            result->suppressed_samples++;
            return 0;
        }
    }
    size_t max_errors = options->max_errors;
    if (max_errors == 0) max_errors = VALIDATOR_DEFAULT_MAX_ERRORS;
    if (result->error_count >= max_errors) {
        result->truncated = 1;
        return 0;
    }

    if (result->error_count == result->error_capacity) {
        size_t new_cap = result->error_capacity ? result->error_capacity * 2 : 32;
        if (new_cap > max_errors) new_cap = max_errors;
        alea_geom_error_t* new_errors =
            realloc(result->errors, new_cap * sizeof(*new_errors));
        if (!new_errors) {
            alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                                  "geometry validator: out of memory appending error");
            return -1;
        }
        result->errors = new_errors;
        result->error_capacity = new_cap;
    }

    result->errors[result->error_count++] = *error;
    if (signature_entry) signature_entry->count++;
    return 0;
}

typedef struct {
    const alea_ray_t* ray;
    const alea_geom_validator_options_t* options;
    alea_geom_validator_result_t* result;
    int report_gaps;
} ray_coverage_finding_context_t;

typedef struct {
    double t_enter;
    double t_exit;
    alea_ray_coverage_kind_t kind;
    size_t owner_count;
    alea_ray_coverage_owner_t owners[VALIDATOR_HIT_CAP];
} ray_coverage_trace_interval_t;

typedef struct {
    ray_coverage_finding_context_t findings;
    ray_coverage_trace_interval_t* intervals;
    size_t count;
    size_t capacity;
    int failed;
} ray_coverage_trace_t;

/* Coverage is the diagnostic oracle for conditions a selected owner trace can
 * hide. Gaps are emitted only by an explicit domain-aware sweep: a legacy
 * validator ray may legitimately start or end outside the model. */
static int append_ray_coverage_finding(
    void* context, const alea_ray_coverage_interval_t* interval) {
    ray_coverage_finding_context_t* ctx = context;
    if (interval->kind == ALEA_RAY_COVERAGE_TRUNCATED ||
        interval->kind == ALEA_RAY_COVERAGE_UNRESOLVED) {
        /* Complete coverage is no longer knowable after owner saturation or
         * an unresolved ownership chain. Preserve earlier findings, but do
         * not continue and present a partial diagnostic as exhaustive. */
        ctx->result->truncated = 1;
        return 1;
    }
    const int is_gap = interval->kind == ALEA_RAY_COVERAGE_GAP &&
        ctx->report_gaps;
    if (interval->kind != ALEA_RAY_COVERAGE_OVERLAP &&
        interval->kind != ALEA_RAY_COVERAGE_UNDEFINED_FILL && !is_gap)
        return 0;

    alea_geom_error_t error;
    init_geom_error(&error);
    error.type = is_gap ? ALEA_GEOM_ERR_INTERIOR_GAP
        : interval->kind == ALEA_RAY_COVERAGE_OVERLAP
            ? ALEA_GEOM_ERR_OVERLAP_AFTER_CROSSING
            : ALEA_GEOM_ERR_UNDEFINED_AFTER_CROSSING;
    error.source = ALEA_GEOM_EVENT_SOURCE_RAY;
    error.t = interval->t_enter;
    error.offset = interval->t_exit - interval->t_enter;
    error.found_cell_count = (int)interval->owner_count;
    if (interval->owner_count > 0) {
        error.found_cell_id = interval->owners[0].cell_id;
        error.universe_id = interval->owners[0].universe_id;
        error.universe_depth = interval->owners[0].depth;
    }
    if (interval->owner_count > 1)
        error.secondary_cell_id = interval->owners[1].cell_id;
    alea_ray_point_at(ctx->ray, interval->t_enter,
                      &error.crossing_point[0], &error.crossing_point[1],
                      &error.crossing_point[2]);
    const double sample_t = interval->t_enter +
        0.381966011250105 * (interval->t_exit - interval->t_enter);
    alea_ray_point_at(ctx->ray, sample_t,
                      &error.sample_point[0], &error.sample_point[1],
                      &error.sample_point[2]);
    error.direction[0] = ctx->ray->dx;
    error.direction[1] = ctx->ray->dy;
    error.direction[2] = ctx->ray->dz;
    if (append_error(ctx->result, ctx->options, &error) != 0)
        return -1;
    return ctx->result->truncated ? 1 : 0;
}

static int append_ray_coverage_gap_finding(
    void* context, const alea_ray_coverage_interval_t* interval) {
    if (interval->kind != ALEA_RAY_COVERAGE_GAP) return 0;
    return append_ray_coverage_finding(context, interval);
}

static void ray_coverage_trace_free(ray_coverage_trace_t* trace) {
    if (!trace) return;
    free(trace->intervals);
    memset(trace, 0, sizeof(*trace));
}

/* Retain the complete, occurrence-sensitive interval result for transition
 * classification as well as emitting its ray-wide findings.  Callback owners
 * are scratch-backed, so the validator must materialize this small trace. */
static int append_ray_coverage_trace(
    void* context, const alea_ray_coverage_interval_t* interval) {
    ray_coverage_trace_t* trace = context;
    const int finding_rc = append_ray_coverage_finding(
        &trace->findings, interval);
    if (finding_rc != 0) return finding_rc;
    if (interval->owner_count > VALIDATOR_HIT_CAP) {
        trace->failed = 1;
        return 1;
    }
    if (trace->count == trace->capacity) {
        const size_t next = trace->capacity ? trace->capacity * 2 : 16;
        ray_coverage_trace_interval_t* grown = realloc(
            trace->intervals, next * sizeof(*grown));
        if (!grown) {
            trace->failed = 1;
            return 1;
        }
        trace->intervals = grown;
        trace->capacity = next;
    }
    ray_coverage_trace_interval_t* saved = &trace->intervals[trace->count++];
    saved->t_enter = interval->t_enter;
    saved->t_exit = interval->t_exit;
    saved->kind = interval->kind;
    saved->owner_count = interval->owner_count;
    memcpy(saved->owners, interval->owners,
           interval->owner_count * sizeof(*interval->owners));
    return 0;
}

/* Match on the canonical (deduplicated) primitive identity, mirroring MCNP's
 * surface equivalence.  Two distinct surface cards that fold onto the same
 * primitive must be treated as the same boundary, otherwise we would invent
 * transition errors that MCNP never sees.  mc_surface_id is report-only. */
static int cell_references_primitive(const alea_system_t* sys,
                                     const alea_cell_entry_t* cell,
                                     uint32_t primitive_id) {
    if (!sys || !cell || primitive_id == ALEA_PRIMITIVE_ID_INVALID ||
        !cell->surface_indices || cell->surface_index_count == 0) {
        return 0;
    }

    for (size_t i = 0; i < cell->surface_index_count; i++) {
        uint32_t surf_idx = cell->surface_indices[i];
        if (surf_idx >= alea_vec_count(&sys->surfaces)) continue;
        if (sys->surfaces.data[surf_idx].primitive_id == primitive_id)
            return 1;
    }
    return 0;
}

static int cell_can_be_implicit_neighbor_by_primitive(
    const alea_system_t* sys,
    const alea_cell_entry_t* cell,
    uint32_t primitive_id) {
    return cell_references_primitive(sys, cell, primitive_id);
}

static int point_inside_cell(const alea_system_t* sys, int cell_idx,
                             const double p[3]);

/* Neighbor lookup keyed on canonical primitive identity instead of
 * mc_surface_id.  alea_find_neighbor_cell() matches on mc_surface_id, which
 * would miss deduplicated equivalent surfaces. */
static int find_neighbor_by_primitive(const alea_system_t* sys,
                                      uint32_t cell_index,
                                      uint32_t primitive_id) {
    if (!sys || primitive_id == ALEA_PRIMITIVE_ID_INVALID ||
        cell_index >= alea_vec_count(&sys->cells) ||
        !sys->cell_adjacency_built) {
        return -1;
    }

    const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
    for (size_t i = 0; i < cell->neighbor_count; i++) {
        uint32_t surf_idx = cell->neighbors[i].surface_index;
        if (surf_idx < alea_vec_count(&sys->surfaces) &&
            sys->surfaces.data[surf_idx].primitive_id == primitive_id) {
            return (int)cell->neighbors[i].neighbor_index;
        }
    }
    return -1;
}

/* A primitive can bound more than two cells: for example, a cylindrical
 * boundary may be partitioned into angular sectors on its outside.  The
 * adjacency table retains all candidates, so resolve the sampled side rather
 * than treating insertion order as geometry. */
static int find_neighbor_by_primitive_at_point(const alea_system_t* sys,
                                               uint32_t cell_index,
                                               uint32_t primitive_id,
                                               const double point[3]) {
    if (!sys || !point || primitive_id == ALEA_PRIMITIVE_ID_INVALID ||
        cell_index >= alea_vec_count(&sys->cells) ||
        !sys->cell_adjacency_built) return -1;

    const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
    for (size_t i = 0; i < cell->neighbor_count; i++) {
        uint32_t surf_idx = cell->neighbors[i].surface_index;
        uint32_t neighbor_idx = cell->neighbors[i].neighbor_index;
        if (surf_idx >= alea_vec_count(&sys->surfaces) ||
            sys->surfaces.data[surf_idx].primitive_id != primitive_id ||
            neighbor_idx >= alea_vec_count(&sys->cells)) continue;
        if (point_inside_cell(sys, (int)neighbor_idx, point))
            return (int)neighbor_idx;
    }
    return -1;
}

static int find_neighbor_by_primitive_in_coverage(
    const alea_system_t* sys, uint32_t cell_index, uint32_t primitive_id,
    const point_coverage_t* cov, const double point[3]) {
    if (!sys || !cov || primitive_id == ALEA_PRIMITIVE_ID_INVALID ||
        cell_index >= alea_vec_count(&sys->cells) ||
        !sys->cell_adjacency_built) return -1;

    const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
    for (size_t i = 0; i < cell->neighbor_count; i++) {
        uint32_t surf_idx = cell->neighbors[i].surface_index;
        uint32_t neighbor_idx = cell->neighbors[i].neighbor_index;
        if (surf_idx >= alea_vec_count(&sys->surfaces) ||
            sys->surfaces.data[surf_idx].primitive_id != primitive_id ||
            neighbor_idx >= alea_vec_count(&sys->cells)) continue;
        if (cov->primary_cell_idx == (int)neighbor_idx ||
            point_coverage_primary_has_cell_index(cov, (int)neighbor_idx) ||
            (alea_cell_entry_is_container(&sys->cells.data[neighbor_idx]) &&
             point_inside_cell(sys, (int)neighbor_idx, point)))
            return (int)neighbor_idx;
    }
    return -1;
}

static int find_point_coverage(alea_system_t* sys,
                               double x, double y, double z,
                               int universe_depth,
                               point_coverage_t* out) {
    if (!sys || !out) return -1;
    memset(out, 0, sizeof(*out));
    out->klass = COVERAGE_NONE;
    out->primary_cell_id = -1;
    out->primary_cell_idx = -1;
    out->primary_occurrence_key = 0;
    out->secondary_cell_id = -1;
    out->secondary_occurrence_key = 0;
    out->universe_id = 0;
    out->depth = universe_depth;
    out->target_depth = universe_depth;

    alea_cell_hit_t hits[VALIDATOR_HIT_CAP];
    uint64_t occurrence_keys[VALIDATOR_HIT_CAP];
    uint64_t parent_occurrence_keys[VALIDATOR_HIT_CAP];
    int n = alea_find_all_cells_at_point_coverage_chain_recursive(
        sys, x, y, z, hits, occurrence_keys, parent_occurrence_keys,
        VALIDATOR_HIT_CAP);
    if (n < 0) return -1;
    if (n >= VALIDATOR_HIT_CAP) out->truncated = 1;
    if (n == 0) return 0;

    int target_depth = universe_depth;
    if (target_depth < 0) {
        target_depth = hits[0].depth;
        for (int i = 1; i < n; i++) {
            if (hits[i].depth > target_depth) target_depth = hits[i].depth;
        }
    }
    out->target_depth = target_depth;

    int count = 0;
    for (int i = 0; i < n; i++) {
        if (hits[i].depth != target_depth) continue;
        if (count == 0) {
            out->primary_cell_id = hits[i].cell_id;
            out->primary_cell_idx = hits[i].cell_index;
            out->primary_occurrence_key = occurrence_keys[i];
            out->universe_id = hits[i].universe_id;
            out->depth = hits[i].depth;
        } else if (count == 1) {
            out->secondary_cell_id = hits[i].cell_id;
            out->secondary_occurrence_key = occurrence_keys[i];
        }
        count++;
    }

    out->count_at_depth = count;
    if (count == 1) out->klass = COVERAGE_ONE;
    else if (count > 1) out->klass = COVERAGE_MULTI;
    if (out->primary_cell_idx >= 0)
        point_coverage_set_primary_ancestors(
            out, hits, occurrence_keys, parent_occurrence_keys, (size_t)n);
    return 0;
}

/* Convert a complete coverage interval into the same target-depth view used
 * by sampled transition validation.  This avoids treating an enclosing
 * container and its filled child as a false overlap. */
static int coverage_interval_to_point(
    const ray_coverage_trace_interval_t* interval, int universe_depth,
    point_coverage_t* out) {
    if (!interval || !out) return 0;
    memset(out, 0, sizeof(*out));
    out->klass = COVERAGE_NONE;
    out->primary_cell_id = -1;
    out->primary_cell_idx = -1;
    out->secondary_cell_id = -1;
    out->universe_id = 0;
    out->depth = universe_depth;
    out->target_depth = universe_depth;
    out->truncated = interval->kind == ALEA_RAY_COVERAGE_TRUNCATED ||
        interval->kind == ALEA_RAY_COVERAGE_UNRESOLVED;
    if (interval->kind == ALEA_RAY_COVERAGE_GAP ||
        interval->kind == ALEA_RAY_COVERAGE_UNRESOLVED ||
        interval->kind == ALEA_RAY_COVERAGE_TRUNCATED ||
        interval->owner_count == 0) {
        return interval->kind != ALEA_RAY_COVERAGE_TRUNCATED &&
            interval->kind != ALEA_RAY_COVERAGE_UNRESOLVED;
    }

    int target_depth = universe_depth;
    if (target_depth < 0) {
        target_depth = interval->owners[0].depth;
        for (size_t i = 1; i < interval->owner_count; i++)
            if (interval->owners[i].depth > target_depth)
                target_depth = interval->owners[i].depth;
    }
    out->target_depth = target_depth;
    int count = 0;
    for (size_t i = 0; i < interval->owner_count; i++) {
        const alea_ray_coverage_owner_t* owner = &interval->owners[i];
        if (owner->depth != target_depth) continue;
        if (count == 0) {
            out->primary_cell_id = owner->cell_id;
            out->primary_cell_idx = owner->cell_index;
            out->primary_occurrence_key = owner->occurrence_key;
            out->universe_id = owner->universe_id;
            out->depth = owner->depth;
        } else if (count == 1) {
            out->secondary_cell_id = owner->cell_id;
            out->secondary_occurrence_key = owner->occurrence_key;
        }
        count++;
    }
    out->count_at_depth = count;
    if (count == 1) out->klass = COVERAGE_ONE;
    else if (count > 1) out->klass = COVERAGE_MULTI;
    if (out->primary_cell_idx >= 0) {
        uint64_t key = out->primary_occurrence_key;
        while (key != 0 && out->primary_ancestor_count < interval->owner_count) {
            size_t index = 0;
            while (index < interval->owner_count &&
                   interval->owners[index].occurrence_key != key)
                index++;
            if (index == interval->owner_count) break;
            out->primary_ancestor_cell_indices[out->primary_ancestor_count++] =
                interval->owners[index].cell_index;
            key = interval->owners[index].parent_occurrence_key;
        }
    }
    return 1;
}

/* Return the interval immediately before or after t.  A merged interval may
 * straddle an irrelevant/tangent primitive hit, which correctly makes both
 * sides resolve to the same complete owner set. */
static int ray_coverage_trace_at(const ray_coverage_trace_t* trace,
                                 double t, int after, int universe_depth,
                                 point_coverage_t* out) {
    if (!trace || !out) return 0;
    for (size_t i = 0; i < trace->count; i++) {
        const ray_coverage_trace_interval_t* interval = &trace->intervals[i];
        if (after) {
            if (interval->t_enter <= t + DEDUP_EPSILON &&
                interval->t_exit > t + DEDUP_EPSILON)
                return coverage_interval_to_point(interval, universe_depth, out);
        } else if (interval->t_enter < t - DEDUP_EPSILON &&
                   interval->t_exit >= t - DEDUP_EPSILON) {
            return coverage_interval_to_point(interval, universe_depth, out);
        }
    }
    return 0;
}

static int same_coverage(const point_coverage_t* a, const point_coverage_t* b) {
    if (a->klass != b->klass) return 0;
    if (a->klass == COVERAGE_ONE &&
        a->primary_occurrence_key != b->primary_occurrence_key) {
        return 0;
    }
    if (a->klass == COVERAGE_MULTI &&
        (a->primary_occurrence_key != b->primary_occurrence_key ||
         a->secondary_occurrence_key != b->secondary_occurrence_key)) {
        return 0;
    }
    return 1;
}

static int sample_coverage_ladder(alea_system_t* sys,
                                  const double crossing_point[3],
                                  const double direction[3],
                                  const alea_geom_validator_options_t* options,
                                  point_coverage_t* out,
                                  double sample_point[3],
                                  double* out_offset,
                                  int* out_ambiguous,
                                  uint32_t* out_flags,
                                  alea_geom_validator_result_t* result) {
    static const double factors[] = { 1.0, 2.0, 5.0, 10.0 };
    double base = options->sample_offset > 0.0
        ? options->sample_offset
        : SURFACE_SAMPLE_OFFSET;

    point_coverage_t first;
    int have_first = 0;
    *out_ambiguous = 0;
    *out_offset = base;
    if (out_flags) *out_flags = 0;

    for (size_t i = 0; i < sizeof(factors) / sizeof(factors[0]); i++) {
        double offset = base * factors[i];
        double p[3] = {
            crossing_point[0] + offset * direction[0],
            crossing_point[1] + offset * direction[1],
            crossing_point[2] + offset * direction[2]
        };
        point_coverage_t cur;
        if (find_point_coverage(sys, p[0], p[1], p[2],
                                options->universe_depth, &cur) != 0) {
            return -1;
        }
        result->exact_queries++;
        if (cur.truncated && out_flags)
            *out_flags |= ALEA_GEOM_EVENT_TRUNCATED_COVERAGE;

        if (!have_first) {
            first = cur;
            copy3(sample_point, p);
            *out_offset = offset;
            have_first = 1;
            continue;
        }
        if (!same_coverage(&first, &cur)) {
            if (i > 1) break;
            *out_ambiguous = 1;
            break;
        }
    }

    *out = first;
    return 0;
}

static int validate_initial_point(alea_system_t* sys,
                                  const alea_ray_t* ray,
                                  const alea_geom_validator_options_t* options,
                                  alea_geom_validator_result_t* result,
                                  point_coverage_t* out_initial) {
    double offset = options->sample_offset > 0.0
        ? options->sample_offset
        : SURFACE_SAMPLE_OFFSET;
    double p[3] = {
        ray->ox + offset * ray->dx,
        ray->oy + offset * ray->dy,
        ray->oz + offset * ray->dz
    };
    point_coverage_t cov;
    if (find_point_coverage(sys, p[0], p[1], p[2],
                            options->universe_depth, &cov) != 0) {
        return -1;
    }
    result->exact_queries++;
    if (out_initial) *out_initial = cov;

    if (!(options->flags & ALEA_GEOM_VALIDATE_STRICT_ADJACENCY) ||
        cov.klass != COVERAGE_MULTI) {
        return 0;
    }

    alea_geom_error_t err;
    init_geom_error(&err);
    err.type = ALEA_GEOM_ERR_OVERLAP_AFTER_CROSSING;
    err.source = ALEA_GEOM_EVENT_SOURCE_INITIAL_POINT;
    err.previous_cell_id = -1;
    err.found_cell_id = cov.primary_cell_id;
    err.expected_neighbor_cell_id = -1;
    err.secondary_cell_id = cov.secondary_cell_id;
    err.found_cell_count = cov.count_at_depth;
    err.surface_id = -1;
    err.universe_id = cov.universe_id;
    err.universe_depth = cov.depth;
    err.crossing_point[0] = ray->ox;
    err.crossing_point[1] = ray->oy;
    err.crossing_point[2] = ray->oz;
    copy3(err.sample_point, p);
    err.direction[0] = ray->dx;
    err.direction[1] = ray->dy;
    err.direction[2] = ray->dz;
    err.t = 0.0;
    err.offset = offset;
    err.flags = ALEA_GEOM_EVENT_INITIAL_POINT;
    if (cov.truncated) err.flags |= ALEA_GEOM_EVENT_TRUNCATED_COVERAGE;
    return append_error(result, options, &err);
}

static int point_inside_cell(const alea_system_t* sys, int cell_idx,
                             const double p[3]) {
    if (!sys || cell_idx < 0 ||
        (size_t)cell_idx >= alea_vec_count(&sys->cells)) {
        return 0;
    }
    const alea_cell_entry_t* cell = &sys->cells.data[cell_idx];
    if (cell->root_node_id == ALEA_NODE_ID_INVALID ||
        cell->root_node_id >= alea_vec_count(&sys->nodes)) {
        return 0;
    }
    return alea_point_inside(sys, cell->root_node_id, p[0], p[1], p[2]) ? 1 : 0;
}

static void coverage_from_cell(const alea_system_t* sys, int cell_idx,
                               point_coverage_t* out) {
    memset(out, 0, sizeof(*out));
    out->klass = COVERAGE_NONE;
    out->primary_cell_id = -1;
    out->primary_cell_idx = -1;
    out->primary_occurrence_key = 0;
    out->secondary_cell_id = -1;
    out->secondary_occurrence_key = 0;
    out->universe_id = 0;
    out->depth = 0;
    out->target_depth = 0;
    if (!sys || cell_idx < 0 ||
        (size_t)cell_idx >= alea_vec_count(&sys->cells)) {
        return;
    }
    const alea_cell_entry_t* cell = &sys->cells.data[cell_idx];
    out->klass = COVERAGE_ONE;
    out->primary_cell_id = cell->mc_cell_id;
    out->primary_cell_idx = cell_idx;
    out->primary_ancestor_cell_indices[0] = cell_idx;
    out->primary_ancestor_count = 1;
    out->universe_id = cell->universe_id;
    out->count_at_depth = 1;
}

/* Shared classifier: given the previous cell and the exact coverage on the
 * "after" side of a boundary, decide which (if any) structured event to emit.
 * Keyed on canonical primitive_id for matching; surface_id is report-only.
 * Used by both the ray-driven and surface/slice-driven drivers. */
static int classify_transition(alea_system_t* sys,
                               int previous_cell_idx,
                               int surface_id,
                               uint32_t primitive_id,
                               const point_coverage_t* cov,
                               int ambiguous,
                               uint32_t event_flags,
                               const double crossing_point[3],
                               const double sample_point[3],
                               const double direction[3],
                               double t,
                               double offset,
                               alea_geom_event_source_t source,
                               size_t curve_index,
                               uint32_t component_index,
                               const double uv[2],
                               const alea_geom_validator_options_t* options,
                               alea_geom_validator_result_t* result) {
    int previous_cell_id = -1;
    int expected_neighbor_idx = -1;
    int expected_neighbor_id = -1;
    int previous_references_surface = 0;
    uint32_t flags = event_flags;

    if (previous_cell_idx >= 0 &&
        (size_t)previous_cell_idx < alea_vec_count(&sys->cells)) {
        const alea_cell_entry_t* prev = &sys->cells.data[previous_cell_idx];
        previous_cell_id = prev->mc_cell_id;
        previous_references_surface =
            cell_references_primitive(sys, prev, primitive_id);
        if (previous_references_surface) {
            expected_neighbor_idx =
                find_neighbor_by_primitive(sys, (uint32_t)previous_cell_idx,
                                           primitive_id);
            if (expected_neighbor_idx >= 0 &&
                (size_t)expected_neighbor_idx < alea_vec_count(&sys->cells)) {
                expected_neighbor_id =
                    sys->cells.data[expected_neighbor_idx].mc_cell_id;
                result->adjacency_hits++;
            } else {
                flags |= ALEA_GEOM_EVENT_MISSING_ADJACENCY;
            }
        } else {
            flags |= ALEA_GEOM_EVENT_PREVIOUS_NO_SURFACE;
        }
    }

    alea_geom_error_t err;
    init_geom_error(&err);
    err.source = source;
    err.previous_cell_id = previous_cell_id;
    err.found_cell_id = cov->primary_cell_id;
    err.expected_neighbor_cell_id = expected_neighbor_id;
    err.secondary_cell_id = cov->secondary_cell_id;
    err.found_cell_count = cov->count_at_depth;
    err.surface_id = surface_id;
    err.primitive_id = primitive_id;
    err.universe_id = cov->universe_id;
    err.universe_depth = cov->depth;
    copy3(err.crossing_point, crossing_point);
    copy3(err.sample_point, sample_point);
    copy3(err.direction, direction);
    err.t = t;
    err.offset = offset;
    err.curve_index = curve_index;
    err.component_index = component_index;
    if (uv) {
        err.uv[0] = uv[0];
        err.uv[1] = uv[1];
    }
    err.flags = flags;

    if (ambiguous || (event_flags & ALEA_GEOM_EVENT_COINCIDENT_SURFACES)) {
        err.type = ALEA_GEOM_ERR_AMBIGUOUS_BOUNDARY;
        result->ambiguous_crossings++;
        return append_error(result, options, &err);
    }

    if (cov->klass == COVERAGE_NONE) {
        if ((options->flags & ALEA_GEOM_VALIDATE_ALLOW_EXTERIOR_VOID) &&
            expected_neighbor_idx < 0) {
            return 0;
        }
        err.type = ALEA_GEOM_ERR_UNDEFINED_AFTER_CROSSING;
        if (event_flags & ALEA_GEOM_EVENT_VIEWPORT_EDGE) {
            err.type = ALEA_GEOM_ERR_AMBIGUOUS_BOUNDARY;
            result->ambiguous_crossings++;
        }
        return append_error(result, options, &err);
    }

    if (cov->klass == COVERAGE_MULTI) {
        /* Ray validation has already emitted the complete extent of this
         * overlap from the coverage sweep.  A sampled boundary can only add a
         * duplicate point finding; keep it for slice validation, where no
         * ray-wide coverage interval exists. */
        if (source == ALEA_GEOM_EVENT_SOURCE_RAY)
            return 0;
        err.type = (event_flags & ALEA_GEOM_EVENT_VIEWPORT_EDGE)
            ? ALEA_GEOM_ERR_AMBIGUOUS_BOUNDARY
            : ALEA_GEOM_ERR_OVERLAP_AFTER_CROSSING;
        if (err.type == ALEA_GEOM_ERR_AMBIGUOUS_BOUNDARY)
            result->ambiguous_crossings++;
        return append_error(result, options, &err);
    }

    if (expected_neighbor_idx >= 0) {
        if (cov->primary_cell_idx != expected_neighbor_idx) {
            if (find_neighbor_by_primitive_in_coverage(
                    sys, (uint32_t)previous_cell_idx, primitive_id, cov,
                    sample_point) >= 0) {
                return 0;
            }
            /* A boundary on a parent fill cell legitimately enters one of
             * that fill's terminal descendants. ``primary_cell_idx`` is the
             * projected/deepest owner, so compare its concrete ancestor chain
             * before claiming a non-adjacent transition. */
            if (cov->klass == COVERAGE_ONE &&
                point_coverage_primary_has_cell_index(cov,
                                                      expected_neighbor_idx)) {
                return 0;
            }
            /* Some legacy recursive point paths publish only the terminal
             * entry despite resolving through a root-level fill. Retain a
             * conservative fallback for that representation: the expected
             * neighbor must itself be a fill/lattice container and must
             * contain the sampled point. Ordinary sibling cells still use the
             * exact occurrence-chain comparison above. */
            if (cov->klass == COVERAGE_ONE &&
                (size_t)expected_neighbor_idx < alea_vec_count(&sys->cells) &&
                alea_cell_entry_is_container(
                    &sys->cells.data[expected_neighbor_idx]) &&
                point_inside_cell(sys, expected_neighbor_idx, sample_point)) {
                return 0;
            }
            if (cov->klass == COVERAGE_ONE &&
                cov->primary_cell_idx >= 0 &&
                (size_t)cov->primary_cell_idx < alea_vec_count(&sys->cells) &&
                cell_can_be_implicit_neighbor_by_primitive(
                    sys, &sys->cells.data[cov->primary_cell_idx],
                    primitive_id)) {
                return 0;
            }
            err.type = (event_flags & ALEA_GEOM_EVENT_VIEWPORT_EDGE)
                ? ALEA_GEOM_ERR_AMBIGUOUS_BOUNDARY
                : ALEA_GEOM_ERR_NON_ADJACENT_TRANSITION;
            if (err.type == ALEA_GEOM_ERR_AMBIGUOUS_BOUNDARY)
                result->ambiguous_crossings++;
            return append_error(result, options, &err);
        }
        return 0;
    }

    if (previous_references_surface) {
        /* Dense shared surfaces are deliberately omitted from the adjacency
         * clique to keep memory bounded.  The resolved occurrence chain still
         * proves a valid transition when one of its enclosing fill/lattice
         * owners references the crossed primitive. */
        if (cov->klass == COVERAGE_ONE &&
            point_coverage_primary_ancestor_references_primitive(
                sys, cov, primitive_id)) {
            return 0;
        }
        if (cov->klass == COVERAGE_ONE &&
            cov->primary_cell_idx >= 0 &&
            (size_t)cov->primary_cell_idx < alea_vec_count(&sys->cells) &&
            cell_can_be_implicit_neighbor_by_primitive(
                sys, &sys->cells.data[cov->primary_cell_idx],
                primitive_id)) {
            return 0;
        }
        if (cov->klass == COVERAGE_ONE)
            err.flags |= ALEA_GEOM_EVENT_FOUND_WITHOUT_ADJACENCY;
        err.type = (event_flags & ALEA_GEOM_EVENT_VIEWPORT_EDGE)
            ? ALEA_GEOM_ERR_AMBIGUOUS_BOUNDARY
            : ALEA_GEOM_ERR_MISSING_NEIGHBOR;
        if (err.type == ALEA_GEOM_ERR_AMBIGUOUS_BOUNDARY)
            result->ambiguous_crossings++;
        return append_error(result, options, &err);
    }

    return 0;
}

static int validate_crossing(alea_system_t* sys,
                             int previous_cell_idx,
                             int surface_id,
                             uint32_t primitive_id,
                             const double crossing_point[3],
                             const double direction[3],
                             double t,
                             uint32_t event_flags,
                             const alea_geom_validator_options_t* options,
                             alea_geom_validator_result_t* result,
                             const point_coverage_t* exact_after,
                             point_coverage_t* out_after) {
    if (surface_id <= 0) return 0;

    point_coverage_t cov;
    double sample_point[3];
    double offset = 0.0;
    int ambiguous = 0;
    uint32_t flags = event_flags;
    if (sample_coverage_ladder(sys, crossing_point, direction, options,
                               &cov, sample_point, &offset, &ambiguous,
                               &flags, result) != 0) {
        return -1;
    }
    if (exact_after) cov = *exact_after;
    if (out_after) *out_after = cov;

    /* The ladder resets flags; restore caller flags (e.g. COINCIDENT_SURFACES,
     * carried-over truncation) so classification still sees them. */
    flags |= event_flags;

    return classify_transition(sys, previous_cell_idx, surface_id, primitive_id,
                               &cov, ambiguous, flags, crossing_point,
                               sample_point, direction, t, offset,
                               ALEA_GEOM_EVENT_SOURCE_RAY, SIZE_MAX, 0, NULL,
                               options, result);
}

static int validate_crossing_fast(alea_system_t* sys,
                                  int previous_cell_idx,
                                  int surface_id,
                                  uint32_t primitive_id,
                                  const double crossing_point[3],
                                  const double direction[3],
                                  double t,
                                  uint32_t event_flags,
                                  const alea_geom_validator_options_t* options,
                                  alea_geom_validator_result_t* result,
                                  const point_coverage_t* exact_after,
                                  point_coverage_t* out_after) {
    if (surface_id <= 0) return 0;

    double offset = options->sample_offset > 0.0
        ? options->sample_offset
        : SURFACE_SAMPLE_OFFSET;
    double sample_point[3] = {
        crossing_point[0] + offset * direction[0],
        crossing_point[1] + offset * direction[1],
        crossing_point[2] + offset * direction[2]
    };

    int previous_cell_id = -1;
    int expected_neighbor_idx = -1;
    int expected_neighbor_id = -1;
    int previous_references_surface = 0;
    uint32_t flags = event_flags;

    if (previous_cell_idx >= 0 &&
        (size_t)previous_cell_idx < alea_vec_count(&sys->cells)) {
        const alea_cell_entry_t* prev = &sys->cells.data[previous_cell_idx];
        previous_cell_id = prev->mc_cell_id;
        previous_references_surface =
            cell_references_primitive(sys, prev, primitive_id);
        if (previous_references_surface) {
            expected_neighbor_idx =
                find_neighbor_by_primitive(sys, (uint32_t)previous_cell_idx,
                                           primitive_id);
            int sampled_neighbor_idx =
                find_neighbor_by_primitive_at_point(
                    sys, (uint32_t)previous_cell_idx, primitive_id,
                    sample_point);
            if (sampled_neighbor_idx >= 0)
                expected_neighbor_idx = sampled_neighbor_idx;
            if (expected_neighbor_idx >= 0 &&
                (size_t)expected_neighbor_idx < alea_vec_count(&sys->cells)) {
                expected_neighbor_id =
                    sys->cells.data[expected_neighbor_idx].mc_cell_id;
                result->adjacency_hits++;
            } else {
                flags |= ALEA_GEOM_EVENT_MISSING_ADJACENCY;
            }
        } else {
            flags |= ALEA_GEOM_EVENT_PREVIOUS_NO_SURFACE;
        }
    }

    if (expected_neighbor_idx >= 0 &&
        point_inside_cell(sys, expected_neighbor_idx, sample_point)) {
        if (out_after) coverage_from_cell(sys, expected_neighbor_idx, out_after);
        return 0;
    }

    if (expected_neighbor_idx < 0 &&
        (options->flags & ALEA_GEOM_VALIDATE_ALLOW_EXTERIOR_VOID)) {
        if (out_after) {
            memset(out_after, 0, sizeof(*out_after));
            out_after->klass = COVERAGE_NONE;
            out_after->primary_cell_id = -1;
            out_after->primary_cell_idx = -1;
            out_after->secondary_cell_id = -1;
        }
        return 0;
    }

    if (expected_neighbor_idx < 0 && previous_references_surface) {
        alea_geom_error_t err;
        init_geom_error(&err);
        err.type = ALEA_GEOM_ERR_MISSING_NEIGHBOR;
        err.source = ALEA_GEOM_EVENT_SOURCE_RAY;
        err.previous_cell_id = previous_cell_id;
        err.found_cell_id = -1;
        err.expected_neighbor_cell_id = expected_neighbor_id;
        err.secondary_cell_id = -1;
        err.found_cell_count = 0;
        err.surface_id = surface_id;
        err.primitive_id = primitive_id;
        copy3(err.crossing_point, crossing_point);
        copy3(err.sample_point, sample_point);
        copy3(err.direction, direction);
        err.t = t;
        err.offset = offset;
        err.flags = flags;
        return append_error(result, options, &err);
    }

    return validate_crossing(sys, previous_cell_idx, surface_id, primitive_id,
                             crossing_point, direction, t, event_flags,
                             options, result, exact_after, out_after);
}

static uint64_t lcg_next(uint64_t* state) {
    *state = *state * UINT64_C(6364136223846793005) +
             UINT64_C(1442695040888963407);
    return *state;
}

static double lcg_double(uint64_t* state) {
    return (double)(lcg_next(state) >> 11) * (1.0 / 9007199254740992.0);
}

static alea_bbox_t validator_bounds(alea_system_t* sys) {
    alea_bbox_t bbox = alea_bbox_empty();
    const alea_universe_t* root = alea_get_universe(sys, 0);
    if (root) bbox = root->bbox;

    if (!alea_bbox_is_valid(&bbox) ||
        bbox.min_x <= -1.0e9 || bbox.max_x >= 1.0e9 ||
        bbox.min_y <= -1.0e9 || bbox.max_y >= 1.0e9 ||
        bbox.min_z <= -1.0e9 || bbox.max_z >= 1.0e9) {
        bbox.min_x = bbox.min_y = bbox.min_z = -10.0;
        bbox.max_x = bbox.max_y = bbox.max_z = 10.0;
    }
    return bbox;
}

static void generate_validation_ray(uint64_t* rng,
                                    const alea_bbox_t* bounds,
                                    alea_ray_t* ray) {
    double cx = 0.5 * (bounds->min_x + bounds->max_x);
    double cy = 0.5 * (bounds->min_y + bounds->max_y);
    double cz = 0.5 * (bounds->min_z + bounds->max_z);
    double bx = bounds->max_x - bounds->min_x;
    double by = bounds->max_y - bounds->min_y;
    double bz = bounds->max_z - bounds->min_z;
    double radius = 0.5 * sqrt(bx * bx + by * by + bz * bz);
    if (radius <= 0.0 || !isfinite(radius)) radius = 10.0;
    radius *= 1.25;

    double z = 2.0 * lcg_double(rng) - 1.0;
    double phi = 2.0 * M_PI * lcg_double(rng);
    double r = sqrt(fmax(0.0, 1.0 - z * z));
    double dx = r * cos(phi);
    double dy = r * sin(phi);
    double dz = z;

    double ox = cx - dx * radius;
    double oy = cy - dy * radius;
    double oz = cz - dz * radius;
    alea_ray_init(ray, ox, oy, oz, dx, dy, dz);
}

static int prepare_validator(alea_system_t* sys,
                             alea_geom_validator_options_t* local_options,
                             const alea_geom_validator_options_t* options) {
    if (options) {
        *local_options = *options;
    } else {
        alea_geom_validator_options_init(local_options);
    }
    if (local_options->ray_count <= 0)
        local_options->ray_count = VALIDATOR_DEFAULT_RAYS;
    if (local_options->sample_offset <= 0.0)
        local_options->sample_offset = SURFACE_SAMPLE_OFFSET;
    if (local_options->flags & ALEA_GEOM_VALIDATE_DOMAIN_BOUNDS) {
        const double* b = local_options->validation_bounds;
        if (!isfinite(b[0]) || !isfinite(b[1]) || !isfinite(b[2]) ||
            !isfinite(b[3]) || !isfinite(b[4]) || !isfinite(b[5]) ||
            b[0] > b[1] || b[2] > b[3] || b[4] > b[5]) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "geometry validator: invalid validation_bounds");
            return -1;
        }
    }

    return alea_system_prepare_query_caches(sys, ALEA_CACHE_RAYCAST);
}

/* Convert the optional closed world-space validation AABB into the scalar
 * coverage sweep's ray-parameter domain.  A ray that misses the domain has no
 * coverage obligation and therefore produces no domain-gap findings. */
static int validator_coverage_domain_for_ray(
    const alea_geom_validator_options_t* options, const alea_ray_t* ray,
    double t_max, alea_ray_coverage_domain_t* out_domain) {
    if (!(options->flags & ALEA_GEOM_VALIDATE_DOMAIN_BOUNDS))
        return 0;
    const double* b = options->validation_bounds;
    const alea_bbox_t bounds = {
        .min_x = b[0], .max_x = b[1],
        .min_y = b[2], .max_y = b[3],
        .min_z = b[4], .max_z = b[5]
    };
    double t_enter = 0.0, t_exit = 0.0;
    if (!ray_bbox_slab_enter_exit(ray, &bounds, 0.0, t_max,
                                  &t_enter, &t_exit))
        return 1;
    *out_domain = (alea_ray_coverage_domain_t){
        .t_min = t_enter,
        .t_max = t_exit,
        .has_domain = 1,
        .report_allowed_exterior = 0
    };
    return 2;
}

static int validate_one_ray(alea_system_t* sys,
                            const alea_ray_t* ray,
                            double t_max,
                            const alea_geom_validator_options_t* options,
                            alea_geom_validator_result_t* result) {
    alea_raycast_result_t coverage_scratch;
    alea_raycast_result_init(&coverage_scratch);
    ray_coverage_trace_t coverage_trace = {
        .findings = { .ray = ray, .options = options, .result = result }
    };
    const int coverage_rc = alea_ray_coverage_sweep_reuse_nocache(
        sys, ray, t_max, &coverage_scratch,
        append_ray_coverage_trace, &coverage_trace);
    alea_raycast_result_free(&coverage_scratch);
    if (coverage_rc < 0 || coverage_trace.failed) {
        ray_coverage_trace_free(&coverage_trace);
        return -1;
    }
    if (result->truncated) {
        ray_coverage_trace_free(&coverage_trace);
        return 0;
    }

    alea_ray_coverage_domain_t domain;
    const int domain_rc = validator_coverage_domain_for_ray(
        options, ray, t_max, &domain);
    if (domain_rc == 1) {
        ray_coverage_trace_free(&coverage_trace);
        return 0;
    }
    if (domain_rc == 2) {
        ray_coverage_finding_context_t gap_findings = {
            .ray = ray, .options = options, .result = result, .report_gaps = 1
        };
        alea_raycast_result_t domain_scratch;
        alea_raycast_result_init(&domain_scratch);
        const int gap_rc = alea_ray_coverage_sweep_domain_reuse_nocache(
            sys, ray, t_max, &domain, &domain_scratch,
            append_ray_coverage_gap_finding, &gap_findings);
        alea_raycast_result_free(&domain_scratch);
        if (gap_rc < 0) {
            ray_coverage_trace_free(&coverage_trace);
            return -1;
        }
        if (result->truncated) {
            ray_coverage_trace_free(&coverage_trace);
            return 0;
        }
    }

    point_coverage_t previous_cov;
    if (validate_initial_point(sys, ray, options, result, &previous_cov) != 0) {
        ray_coverage_trace_free(&coverage_trace);
        return -1;
    }

    alea_raycast_result_t ray_result;
    alea_raycast_result_init(&ray_result);
    /* Use the diagnostic breakpoint engine rather than bare root surfaces:
     * fills and lattice transitions must participate in the same crossing
     * vocabulary as boundary provenance and coverage validation. */
    int rc = alea_raycast_global_breakpoints_reuse_nocache(
        sys, ray, t_max, &ray_result);
    if (rc != 0) {
        alea_raycast_result_free(&ray_result);
        ray_coverage_trace_free(&coverage_trace);
        return -1;
    }

    size_t max_crossings = options->max_crossings;
    if (max_crossings == 0) max_crossings = VALIDATOR_DEFAULT_MAX_CROSSINGS;

    for (size_t i = 0; i < ray_result.hits.count;) {
        if (result->crossings_checked >= max_crossings) {
            result->truncated = 1;
            break;
        }
        if (result->truncated) break;

        /* One diagnostic crossing is one coincident distance group, not one
         * raw mathematical hit.  Surface cards folded onto a primitive and
         * genuinely coincident primitives must therefore share the same
         * before/after coverage sample and yield at most one finding. */
        const double t = ray_result.hits.data[i].t;
        size_t group_end = i + 1;
        while (group_end < ray_result.hits.count &&
               fabs(ray_result.hits.data[group_end].t - t) <= DEDUP_EPSILON) {
            group_end++;
        }

        size_t physical_index = SIZE_MAX;
        uint32_t physical_primitive = ALEA_PRIMITIVE_ID_INVALID;
        int has_synthetic = 0;
        uint32_t event_flags = previous_cov.truncated
            ? ALEA_GEOM_EVENT_TRUNCATED_COVERAGE : 0;
        for (size_t j = i; j < group_end; j++) {
            const alea_ray_hit_t* candidate = &ray_result.hits.data[j];
            if (candidate->surface_id == 0) {
                has_synthetic = 1;
            } else if (candidate->surface_id > 0) {
                if (physical_index == SIZE_MAX) {
                    physical_index = j;
                    physical_primitive = candidate->primitive_id;
                } else if (candidate->primitive_id != physical_primitive) {
                    event_flags |= ALEA_GEOM_EVENT_COINCIDENT_SURFACES;
                }
            }
        }

        /* A pure lattice DDA boundary has no physical surface to validate,
         * but it still changes the concrete occurrence that later physical
         * crossings must use as their predecessor.  A synthetic event that
         * coincides with a physical one is covered by the physical crossing;
         * otherwise it refreshes coverage without an adjacency verdict. */
        if (physical_index == SIZE_MAX) {
            if (!has_synthetic) {
                i = group_end;
                continue;
            }
            double crossing[3];
            alea_ray_point_at(ray, t, &crossing[0], &crossing[1],
                              &crossing[2]);
            double direction[3] = { ray->dx, ray->dy, ray->dz };
            alea_geom_validator_options_t depth_options;
            const alea_geom_validator_options_t* sampling_options = options;
            if (options->universe_depth < 0 &&
                previous_cov.klass == COVERAGE_ONE) {
                depth_options = *options;
                depth_options.universe_depth = previous_cov.depth;
                sampling_options = &depth_options;
            }
            point_coverage_t after_cov;
            double sample_point[3];
            double offset = 0.0;
            int ambiguous = 0;
            uint32_t flags = 0;
            if (sample_coverage_ladder(sys, crossing, direction,
                                       sampling_options, &after_cov,
                                       sample_point, &offset, &ambiguous,
                                       &flags, result) != 0) {
                alea_raycast_result_free(&ray_result);
                ray_coverage_trace_free(&coverage_trace);
                return -1;
            }
            point_coverage_t exact_after;
            if (ray_coverage_trace_at(&coverage_trace, t, 1,
                                      sampling_options->universe_depth,
                                      &exact_after) &&
                exact_after.klass != COVERAGE_NONE) {
                after_cov = exact_after;
            }
            previous_cov = after_cov;
            result->crossings_checked++;
            i = group_end;
            continue;
        }

        const alea_ray_hit_t* hit = &ray_result.hits.data[physical_index];

        point_coverage_t exact_before;
        if (ray_coverage_trace_at(&coverage_trace, t, 0,
                                  options->universe_depth, &exact_before)) {
            previous_cov = exact_before;
        }
        int previous_cell_idx = (previous_cov.klass == COVERAGE_ONE)
            ? previous_cov.primary_cell_idx
            : -1;

        double crossing[3];
        alea_ray_point_at(ray, t, &crossing[0], &crossing[1], &crossing[2]);
        double direction[3] = { ray->dx, ray->dy, ray->dz };
        point_coverage_t after_cov;
        alea_geom_validator_options_t depth_options;
        const alea_geom_validator_options_t* crossing_options = options;
        if (options->universe_depth < 0 && previous_cov.klass == COVERAGE_ONE) {
            depth_options = *options;
            depth_options.universe_depth = previous_cov.depth;
            crossing_options = &depth_options;
        }
        point_coverage_t exact_after;
        const point_coverage_t* exact_after_ptr = NULL;
        if (ray_coverage_trace_at(&coverage_trace, t, 1,
                                  crossing_options->universe_depth,
                                  &exact_after) &&
            exact_after.klass != COVERAGE_NONE) {
            exact_after_ptr = &exact_after;
        }

        if (crossing_options->flags & ALEA_GEOM_VALIDATE_STRICT_ADJACENCY) {
            rc = validate_crossing(sys, previous_cell_idx, hit->surface_id,
                                   hit->primitive_id, crossing, direction,
                                   hit->t, event_flags,
                                   crossing_options, result, exact_after_ptr,
                                   &after_cov);
        } else {
            rc = validate_crossing_fast(sys, previous_cell_idx, hit->surface_id,
                                        hit->primitive_id, crossing, direction,
                                        hit->t, event_flags,
                                        crossing_options, result, exact_after_ptr,
                                        &after_cov);
        }
        if (rc != 0) {
            alea_raycast_result_free(&ray_result);
            ray_coverage_trace_free(&coverage_trace);
            return -1;
        }
        previous_cov = after_cov;
        result->crossings_checked++;
        i = group_end;
    }

    alea_raycast_result_free(&ray_result);
    ray_coverage_trace_free(&coverage_trace);
    return 0;
}

int alea_validate_geometry(alea_system_t* sys,
                           const alea_geom_validator_options_t* options,
                           alea_geom_validator_result_t* result) {
    if (!sys || !result) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "alea_validate_geometry: NULL argument");
        return -1;
    }

    alea_geom_validator_options_t local_options;
    if (prepare_validator(sys, &local_options, options) != 0)
        return -1;
    if (!(local_options.flags & ALEA_GEOM_VALIDATE_RAYS))
        return 0;

    alea_bbox_t bounds = validator_bounds(sys);
    double bx = bounds.max_x - bounds.min_x;
    double by = bounds.max_y - bounds.min_y;
    double bz = bounds.max_z - bounds.min_z;
    double diag = sqrt(bx * bx + by * by + bz * bz);
    if (diag <= 0.0 || !isfinite(diag)) diag = 20.0;
    double t_max = local_options.t_max > 0.0
        ? local_options.t_max
        : 3.0 * diag;

    uint64_t rng = local_options.seed ? local_options.seed : UINT64_C(42);
    for (int i = 0; i < local_options.ray_count; i++) {
        if (g_alea_interrupted) {
            alea_set_error_detail(ALEA_ERR_INTERRUPTED,
                                  "geometry validation interrupted");
            return -1;
        }
        if (result->truncated) break;
        alea_ray_t ray;
        generate_validation_ray(&rng, &bounds, &ray);
        if (validate_one_ray(sys, &ray, t_max, &local_options, result) != 0)
            return -1;
    }

    return 0;
}

int alea_validate_geometry_ray(alea_system_t* sys,
                               const alea_geom_validator_options_t* options,
                               double ox, double oy, double oz,
                               double dx, double dy, double dz,
                               double t_max,
                               alea_geom_validator_result_t* result) {
    if (!sys || !result) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "alea_validate_geometry_ray: NULL argument");
        return -1;
    }

    alea_geom_validator_options_t local_options;
    if (prepare_validator(sys, &local_options, options) != 0)
        return -1;

    alea_ray_t ray;
    if (alea_ray_init(&ray, ox, oy, oz, dx, dy, dz) != 0) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "alea_validate_geometry_ray: zero-length direction");
        return -1;
    }

    double effective_t_max = t_max > 0.0 ? t_max : local_options.t_max;
    if (effective_t_max <= 0.0) {
        alea_bbox_t bounds = validator_bounds(sys);
        double bx = bounds.max_x - bounds.min_x;
        double by = bounds.max_y - bounds.min_y;
        double bz = bounds.max_z - bounds.min_z;
        double diag = sqrt(bx * bx + by * by + bz * bz);
        if (diag <= 0.0 || !isfinite(diag)) diag = 20.0;
        effective_t_max = 3.0 * diag;
    }

    return validate_one_ray(sys, &ray, effective_t_max, &local_options, result);
}

/* ------------------------------------------------------------------------- *
 * Surface/slice-driven validation
 * ------------------------------------------------------------------------- */

/* Validate one boundary sample: world point `p` with in-plane normal mapped to
 * the 3D direction `dir`.  Evaluates exact coverage on both sides and reuses the
 * shared classifier.  The side carrying a single cell is treated as "previous"
 * so the adjacency comparison engages; orientation is chosen so a MULTI/NONE
 * side is surfaced as the "after" coverage. */
static int validate_surface_sample(alea_system_t* sys,
                                   const double p[3],
                                   const double dir[3],
                                   int surface_id,
                                   uint32_t primitive_id,
                                   size_t curve_index,
                                   uint32_t component_index,
                                   double curve_t,
                                   const double uv[2],
                                   uint32_t event_flags,
                                   const alea_geom_validator_options_t* options,
                                   alea_geom_validator_result_t* result) {
    double neg_dir[3] = { -dir[0], -dir[1], -dir[2] };

    point_coverage_t cov_plus, cov_minus;
    double sp_plus[3], sp_minus[3];
    double off_plus = 0.0, off_minus = 0.0;
    int amb_plus = 0, amb_minus = 0;
    uint32_t flags_plus = 0, flags_minus = 0;

    if (sample_coverage_ladder(sys, p, dir, options, &cov_plus, sp_plus,
                               &off_plus, &amb_plus, &flags_plus, result) != 0)
        return -1;
    if (sample_coverage_ladder(sys, p, neg_dir, options, &cov_minus, sp_minus,
                               &off_minus, &amb_minus, &flags_minus, result) != 0)
        return -1;

    int ambiguous = amb_plus || amb_minus;
    uint32_t flags = flags_plus | flags_minus | event_flags;

    if (!ambiguous &&
        cov_plus.klass == COVERAGE_ONE &&
        cov_minus.klass == COVERAGE_ONE &&
        cov_plus.primary_occurrence_key == cov_minus.primary_occurrence_key) {
        return 0;
    }

    /* Choose orientation: prefer a single-cell side as "previous". */
    int prev_idx;
    const point_coverage_t* after;
    const double* after_sp;
    const double* after_dir;
    double after_off;

    if (cov_minus.klass == COVERAGE_ONE) {
        prev_idx = cov_minus.primary_cell_idx;
        after = &cov_plus;  after_sp = sp_plus;  after_dir = dir;  after_off = off_plus;
    } else if (cov_plus.klass == COVERAGE_ONE) {
        prev_idx = cov_plus.primary_cell_idx;
        after = &cov_minus; after_sp = sp_minus; after_dir = neg_dir; after_off = off_minus;
    } else {
        /* Neither side is a clean single cell: surface the more severe side. */
        prev_idx = -1;
        if (cov_minus.klass == COVERAGE_MULTI) {
            after = &cov_minus; after_sp = sp_minus; after_dir = neg_dir; after_off = off_minus;
        } else {
            after = &cov_plus;  after_sp = sp_plus;  after_dir = dir;  after_off = off_plus;
        }
    }

    return classify_transition(sys, prev_idx, surface_id, primitive_id,
                               after, ambiguous, flags, p, after_sp, after_dir,
                               curve_t, after_off,
                               ALEA_GEOM_EVENT_SOURCE_SLICE_CURVE,
                               curve_index, component_index, uv,
                               options, result);
}

static int slice_sample_is_on_viewport_edge(const alea_slice_view_t* view,
                                            double u, double v) {
    double u_span = fabs(view->u_max - view->u_min);
    double v_span = fabs(view->v_max - view->v_min);
    double scale = fmax(1.0, fmax(u_span, v_span));
    /* Curve clipping and plane transforms accumulate more than a few ulps on
     * kilometre-scale views. This remains a sub-pixel geometric tolerance,
     * not a probe distance. */
    double tolerance = fmax(64.0 * DBL_EPSILON * scale, 1.0e-9 * scale);
    return fabs(u - view->u_min) <= tolerance ||
           fabs(u - view->u_max) <= tolerance ||
           fabs(v - view->v_min) <= tolerance ||
           fabs(v - view->v_max) <= tolerance;
}

int alea_validate_geometry_slice(alea_system_t* sys,
                                 const alea_slice_view_t* view,
                                 const alea_slice_curves_t* curves,
                                 const alea_geom_validator_options_t* options,
                                 alea_geom_validator_result_t* result) {
    if (!sys || !view || !curves || !result) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG,
                              "alea_validate_geometry_slice: NULL argument");
        return -1;
    }

    alea_geom_validator_options_t local;
    if (prepare_validator(sys, &local, options) != 0)
        return -1;

    size_t ncurves = alea_slice_curves_count(curves);
    double vp_w = view->u_max - view->u_min;
    double vp_h = view->v_max - view->v_min;
    double vp_diag = sqrt(vp_w * vp_w + vp_h * vp_h);
    double sample_spacing = vp_diag / 200.0;
    if (!(sample_spacing > 0.0) || !isfinite(sample_spacing))
        sample_spacing = 0.05;

    size_t max_crossings = local.max_crossings;
    if (max_crossings == 0) max_crossings = VALIDATOR_DEFAULT_MAX_CROSSINGS;

    const double* origin = view->plane.origin;
    const double* u_axis = view->plane.u_axis;
    const double* v_axis = view->plane.v_axis;

    int rc = 0;
    for (size_t ci = 0; ci < ncurves; ci++) {
        if (g_alea_interrupted) {
            alea_set_error_detail(ALEA_ERR_INTERRUPTED,
                                  "geometry slice validation interrupted");
            rc = -1;
            break;
        }
        if (result->truncated) break;

        alea_curve_t c;
        if (alea_slice_curves_get(curves, ci, &c) != 0) continue;
        if (c.type == ALEA_CURVE_NONE || c.type == ALEA_CURVE_POINT) continue;
        if (c.surface_id <= 0) continue;  /* synthetic boundary, not physical */

        double t_lo, t_hi;
        alea_slice_curve_param_range(curves, ci, view, &t_lo, &t_hi);
        if (!(t_hi > t_lo)) continue;

        /* Approximate arc length to pick a sample count (mirrors slice checker). */
        double arc = t_hi - t_lo;
        if (c.type == ALEA_CURVE_CIRCLE || c.type == ALEA_CURVE_ARC) {
            arc = (t_hi - t_lo) * c.data.circle.radius;
        } else if (c.type == ALEA_CURVE_ELLIPSE || c.type == ALEA_CURVE_ELLIPSE_ARC) {
            arc = (t_hi - t_lo) * 0.5 * (c.data.ellipse.semi_a + c.data.ellipse.semi_b);
        }
        int n_samples = (int)(arc / sample_spacing);
        if (n_samples < 2) n_samples = 2;
        if (n_samples > 2000) n_samples = 2000;

        double dt_finite = (t_hi - t_lo) * 1e-6;
        if (dt_finite < 1e-15) dt_finite = 1e-15;

        for (int i = 0; i <= n_samples; i++) {
            if (result->crossings_checked >= max_crossings) {
                result->truncated = 1;
                break;
            }
            if (result->truncated) break;

            double t = t_lo + (t_hi - t_lo) * (double)i / (double)n_samples;
            double u, v;
            if (alea_slice_curve_eval(curves, ci, t, &u, &v) != 0) continue;
            if (u < view->u_min || u > view->u_max ||
                v < view->v_min || v > view->v_max) continue;

            /* Normal from a finite-difference tangent in plane coordinates. */
            double t_back = t - dt_finite, t_fwd = t + dt_finite;
            if (t_back < t_lo) t_back = t_lo;
            if (t_fwd > t_hi) t_fwd = t_hi;
            double u1, v1, u2, v2;
            if (alea_slice_curve_eval(curves, ci, t_back, &u1, &v1) != 0 ||
                alea_slice_curve_eval(curves, ci, t_fwd, &u2, &v2) != 0) continue;
            double tang_u = u2 - u1, tang_v = v2 - v1;
            double tlen = sqrt(tang_u * tang_u + tang_v * tang_v);
            if (tlen < 1e-20) continue;
            double nu = -tang_v / tlen, nv = tang_u / tlen;

            double pw[3] = {
                origin[0] + u * u_axis[0] + v * v_axis[0],
                origin[1] + u * u_axis[1] + v * v_axis[1],
                origin[2] + u * u_axis[2] + v * v_axis[2]
            };
            double dirw[3] = {
                nu * u_axis[0] + nv * v_axis[0],
                nu * u_axis[1] + nv * v_axis[1],
                nu * u_axis[2] + nv * v_axis[2]
            };

            double uv[2] = { u, v };
            uint32_t event_flags = slice_sample_is_on_viewport_edge(
                view, u, v) ? ALEA_GEOM_EVENT_VIEWPORT_EDGE : 0;
            if (validate_surface_sample(sys, pw, dirw, c.surface_id,
                                        c.primitive_id, ci, 0, t, uv, event_flags,
                                        &local, result) != 0) {
                rc = -1;
                break;
            }
            result->crossings_checked++;
        }
        if (rc != 0) break;
    }
    return rc;
}
