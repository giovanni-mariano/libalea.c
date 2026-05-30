// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "geo_validator.h"

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
    int secondary_cell_id;
    int universe_id;
    int depth;
    int count_at_depth;
    int target_depth;
    int truncated;
} point_coverage_t;

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

void alea_geom_validator_options_init(alea_geom_validator_options_t* options) {
    if (!options) return;
    memset(options, 0, sizeof(*options));
    options->flags = ALEA_GEOM_VALIDATE_RAYS |
                     ALEA_GEOM_VALIDATE_STRICT_ADJACENCY;
    options->universe_depth = -1;
    options->max_errors = VALIDATOR_DEFAULT_MAX_ERRORS;
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

static int find_point_coverage(alea_system_t* sys,
                               double x, double y, double z,
                               int universe_depth,
                               point_coverage_t* out) {
    if (!sys || !out) return -1;
    memset(out, 0, sizeof(*out));
    out->klass = COVERAGE_NONE;
    out->primary_cell_id = -1;
    out->primary_cell_idx = -1;
    out->secondary_cell_id = -1;
    out->universe_id = 0;
    out->depth = universe_depth;
    out->target_depth = universe_depth;

    alea_cell_hit_t hits[VALIDATOR_HIT_CAP];
    int n = alea_find_all_cells_at_point_recursive(sys, x, y, z,
                                                   hits, VALIDATOR_HIT_CAP);
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
            out->universe_id = hits[i].universe_id;
            out->depth = hits[i].depth;
        } else if (count == 1) {
            out->secondary_cell_id = hits[i].cell_id;
        }
        count++;
    }

    out->count_at_depth = count;
    if (count == 1) out->klass = COVERAGE_ONE;
    else if (count > 1) out->klass = COVERAGE_MULTI;
    return 0;
}

static int same_coverage(const point_coverage_t* a, const point_coverage_t* b) {
    if (a->klass != b->klass) return 0;
    if (a->klass == COVERAGE_ONE &&
        a->primary_cell_idx != b->primary_cell_idx) {
        return 0;
    }
    if (a->klass == COVERAGE_MULTI &&
        (a->primary_cell_id != b->primary_cell_id ||
         a->secondary_cell_id != b->secondary_cell_id)) {
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

static int point_inside_cell(alea_system_t* sys, int cell_idx,
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
    out->secondary_cell_id = -1;
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
        return append_error(result, options, &err);
    }

    if (cov->klass == COVERAGE_MULTI) {
        err.type = ALEA_GEOM_ERR_OVERLAP_AFTER_CROSSING;
        return append_error(result, options, &err);
    }

    if (expected_neighbor_idx >= 0) {
        if (cov->primary_cell_idx != expected_neighbor_idx) {
            err.type = ALEA_GEOM_ERR_NON_ADJACENT_TRANSITION;
            return append_error(result, options, &err);
        }
        return 0;
    }

    if (previous_references_surface) {
        if (cov->klass == COVERAGE_ONE)
            err.flags |= ALEA_GEOM_EVENT_FOUND_WITHOUT_ADJACENCY;
        err.type = ALEA_GEOM_ERR_MISSING_NEIGHBOR;
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
                             options, result, out_after);
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

    unsigned cache_flags = alea_system_spatial_mode_prefers_hier(sys)
        ? ALEA_CACHE_RAYCAST_HIER
        : ALEA_CACHE_RAYCAST;
    if (local_options->flags & ALEA_GEOM_VALIDATE_HIERARCHICAL)
        cache_flags = ALEA_CACHE_RAYCAST_HIER;
    return alea_system_prepare_query_caches(sys, cache_flags);
}

static int validate_one_ray(alea_system_t* sys,
                            const alea_ray_t* ray,
                            double t_max,
                            const alea_geom_validator_options_t* options,
                            alea_geom_validator_result_t* result) {
    point_coverage_t previous_cov;
    if (validate_initial_point(sys, ray, options, result, &previous_cov) != 0)
        return -1;

    alea_raycast_result_t ray_result;
    alea_raycast_result_init(&ray_result);
    int rc = alea_raycast_surfaces_nocache(sys, ray, 0.0, t_max, &ray_result);
    if (rc != 0) {
        alea_raycast_result_free(&ray_result);
        return -1;
    }

    size_t max_crossings = options->max_crossings;
    if (max_crossings == 0) max_crossings = VALIDATOR_DEFAULT_MAX_CROSSINGS;

    for (size_t i = 0; i < ray_result.hits.count; i++) {
        if (result->crossings_checked >= max_crossings) {
            result->truncated = 1;
            break;
        }
        if (result->truncated) break;

        alea_ray_hit_t* hit = &ray_result.hits.data[i];
        if (hit->surface_id <= 0) continue;

        /* Collapse duplicate hits from deduplicated surface cards: the same
         * physical primitive can appear once per mc_surface_id at the same t.
         * Treat them as a single crossing so we neither double-classify nor
         * report a spurious coincident-surface ambiguity. */
        if (i > 0 &&
            fabs(ray_result.hits.data[i - 1].t - hit->t) <= DEDUP_EPSILON &&
            ray_result.hits.data[i - 1].primitive_id == hit->primitive_id) {
            continue;
        }

        int previous_cell_idx = (previous_cov.klass == COVERAGE_ONE)
            ? previous_cov.primary_cell_idx
            : -1;

        uint32_t event_flags = 0;
        if (previous_cov.truncated) event_flags |= ALEA_GEOM_EVENT_TRUNCATED_COVERAGE;

        /* Coincident-surface ambiguity is keyed on canonical primitive identity:
         * two hits at the same t are only genuinely coincident if they are
         * distinct primitives, not duplicate cards of one primitive. */
        for (size_t j = i + 1; j < ray_result.hits.count; j++) {
            if (fabs(ray_result.hits.data[j].t - hit->t) > DEDUP_EPSILON)
                break;
            if (ray_result.hits.data[j].primitive_id != hit->primitive_id) {
                event_flags |= ALEA_GEOM_EVENT_COINCIDENT_SURFACES;
                break;
            }
        }
        if (i > 0 &&
            fabs(ray_result.hits.data[i - 1].t - hit->t) <= DEDUP_EPSILON &&
            ray_result.hits.data[i - 1].primitive_id != hit->primitive_id) {
            event_flags |= ALEA_GEOM_EVENT_COINCIDENT_SURFACES;
        }

        double crossing[3];
        alea_ray_point_at(ray, hit->t, &crossing[0], &crossing[1], &crossing[2]);
        double direction[3] = { ray->dx, ray->dy, ray->dz };
        point_coverage_t after_cov;
        alea_geom_validator_options_t depth_options;
        const alea_geom_validator_options_t* crossing_options = options;
        if (options->universe_depth < 0 && previous_cov.klass == COVERAGE_ONE) {
            depth_options = *options;
            depth_options.universe_depth = previous_cov.depth;
            crossing_options = &depth_options;
        }

        if (crossing_options->flags & ALEA_GEOM_VALIDATE_STRICT_ADJACENCY) {
            rc = validate_crossing(sys, previous_cell_idx, hit->surface_id,
                                   hit->primitive_id, crossing, direction,
                                   hit->t, event_flags,
                                   crossing_options, result, &after_cov);
        } else {
            rc = validate_crossing_fast(sys, previous_cell_idx, hit->surface_id,
                                        hit->primitive_id, crossing, direction,
                                        hit->t, event_flags,
                                        crossing_options, result, &after_cov);
        }
        if (rc != 0) {
            alea_raycast_result_free(&ray_result);
            return -1;
        }
        previous_cov = after_cov;
        result->crossings_checked++;
    }

    alea_raycast_result_free(&ray_result);
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
    uint32_t flags = flags_plus | flags_minus;

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
            if (validate_surface_sample(sys, pw, dirw, c.surface_id,
                                        c.primitive_id, ci, 0, t, uv,
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
