// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_geo_validator.h"

#include "raycast/raycast.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

struct alea_ray_slice_validation_result {
    size_t row_count;
    size_t interval_count;
    uint32_t fields;
    uint32_t executed_trace_mask;
    uint32_t reused_trace_mask;
    alea_ray_slice_refinement_status_t refinement_status;
    uint64_t* row_offsets;
    double* row_transverse_coordinates;
    uint8_t* row_direction_tags;
    size_t* row_base_indices;
    double* u_enter;
    double* u_exit;
    uint32_t* diagnostic_flags;
    int32_t* fast_forward_cell_ids;
    int32_t* fast_reverse_cell_ids;
    uint64_t* fast_forward_occurrence_keys;
    uint64_t* fast_reverse_occurrence_keys;
    uint32_t* coverage_owner_counts;
    int32_t* u_enter_forward_surface_ids;
    int32_t* u_enter_reverse_surface_ids;
    int32_t* u_exit_forward_surface_ids;
    int32_t* u_exit_reverse_surface_ids;
    uint32_t* u_enter_provenance_flags;
    uint32_t* u_exit_provenance_flags;
};

typedef struct {
    double* u_enter;
    double* u_exit;
    uint32_t* flags;
    int32_t* forward_cells;
    int32_t* reverse_cells;
    uint64_t* forward_keys;
    uint64_t* reverse_keys;
    uint32_t* coverage_owners;
    size_t count;
    size_t capacity;
} validation_row_t;

typedef struct {
    const alea_raycast_batch_result_t* result;
    const uint64_t* offsets;
    const double* t_enter;
    const double* t_exit;
    const int32_t* cells;
    const uint64_t* occurrence_keys;
    int reverse;
    double u_min;
    double u_max;
} trace_view_t;

static void validation_result_free_buffers(
    alea_ray_slice_validation_result_t* result) {
    if (!result) return;
    free(result->row_offsets);
    free(result->row_transverse_coordinates);
    free(result->row_direction_tags);
    free(result->row_base_indices);
    free(result->u_enter);
    free(result->u_exit);
    free(result->diagnostic_flags);
    free(result->fast_forward_cell_ids);
    free(result->fast_reverse_cell_ids);
    free(result->fast_forward_occurrence_keys);
    free(result->fast_reverse_occurrence_keys);
    free(result->coverage_owner_counts);
    free(result->u_enter_forward_surface_ids);
    free(result->u_enter_reverse_surface_ids);
    free(result->u_exit_forward_surface_ids);
    free(result->u_exit_reverse_surface_ids);
    free(result->u_enter_provenance_flags);
    free(result->u_exit_provenance_flags);
    memset(result, 0, sizeof(*result));
}

static void validation_row_free(validation_row_t* row) {
    if (!row) return;
    free(row->u_enter);
    free(row->u_exit);
    free(row->flags);
    free(row->forward_cells);
    free(row->reverse_cells);
    free(row->forward_keys);
    free(row->reverse_keys);
    free(row->coverage_owners);
    memset(row, 0, sizeof(*row));
}

static int validation_row_reserve(validation_row_t* row, size_t needed) {
    if (needed <= row->capacity) return 0;
    size_t capacity = row->capacity ? row->capacity * 2 : 8;
    if (capacity < needed) capacity = needed;
    if (capacity > SIZE_MAX / sizeof(*row->u_enter)) return -1;

#define GROW_ROW_FIELD(field) do { \
    void* p = realloc(row->field, capacity * sizeof(*row->field)); \
    if (!p) return -1; \
    row->field = p; \
} while (0)
    GROW_ROW_FIELD(u_enter);
    GROW_ROW_FIELD(u_exit);
    GROW_ROW_FIELD(flags);
    GROW_ROW_FIELD(forward_cells);
    GROW_ROW_FIELD(reverse_cells);
    GROW_ROW_FIELD(forward_keys);
    GROW_ROW_FIELD(reverse_keys);
    GROW_ROW_FIELD(coverage_owners);
#undef GROW_ROW_FIELD
    row->capacity = capacity;
    return 0;
}

static int validation_row_append(validation_row_t* row,
                                 double u_enter, double u_exit,
                                 uint32_t flags, int32_t forward_cell,
                                 int32_t reverse_cell, uint64_t forward_key,
                                 uint64_t reverse_key,
                                 uint32_t coverage_owners) {
    if (validation_row_reserve(row, row->count + 1) != 0) return -1;
    size_t i = row->count++;
    row->u_enter[i] = u_enter;
    row->u_exit[i] = u_exit;
    row->flags[i] = flags;
    row->forward_cells[i] = forward_cell;
    row->reverse_cells[i] = reverse_cell;
    row->forward_keys[i] = forward_key;
    row->reverse_keys[i] = reverse_key;
    row->coverage_owners[i] = coverage_owners;
    return 0;
}

/* Map one complete-coverage classification onto the published diagnostic
 * flags.  A unique chain contributes nothing: it is the expected answer. */
static uint32_t coverage_kind_flag(uint8_t kind) {
    switch ((alea_ray_coverage_kind_t)kind) {
    case ALEA_RAY_COVERAGE_GAP:
        return ALEA_RAY_SLICE_DIAG_COVERAGE_GAP;
    case ALEA_RAY_COVERAGE_OVERLAP:
        return ALEA_RAY_SLICE_DIAG_COVERAGE_OVERLAP;
    case ALEA_RAY_COVERAGE_UNDEFINED_FILL:
        return ALEA_RAY_SLICE_DIAG_COVERAGE_UNDEFINED_FILL;
    case ALEA_RAY_COVERAGE_UNRESOLVED:
        return ALEA_RAY_SLICE_DIAG_COVERAGE_UNRESOLVED;
    case ALEA_RAY_COVERAGE_TRUNCATED:
        return ALEA_RAY_SLICE_DIAG_COVERAGE_TRUNCATED;
    case ALEA_RAY_COVERAGE_ALLOWED_EXTERIOR:
        return ALEA_RAY_SLICE_DIAG_COVERAGE_ALLOWED_EXTERIOR;
    case ALEA_RAY_COVERAGE_UNIQUE:
        break;
    }
    return 0;
}

/* Coverage rows are traced in view U coordinates: the row ray starts at u_min
 * and advances along +u_axis, so u == u_min + t.  This matches the forward
 * slice trace, which shifts its own t values by u_min. */
static int coverage_row_span(const alea_ray_coverage_slice_result_t* coverage,
                             size_t row, size_t* begin, size_t* end) {
    if (!coverage || row >= coverage->row_count) return -1;
    *begin = coverage->row_offsets[row];
    *end = coverage->row_offsets[row + 1];
    return (*begin <= *end && *end <= coverage->interval_count) ? 0 : -1;
}

static void coverage_flags_at_u(
    const alea_ray_coverage_slice_result_t* coverage, size_t begin, size_t end,
    double u_min, double u, uint32_t coverage_flags, uint32_t* out_flags,
    uint32_t* out_owners) {
    *out_flags = 0;
    *out_owners = 0;
    for (size_t i = begin; i < end; i++) {
        const double enter = u_min + coverage->t_enter[i];
        const double exit = u_min + coverage->t_exit[i];
        if (u <= enter || u >= exit) continue;
        uint32_t flag = coverage_kind_flag(coverage->kinds[i]);
        /* Under the uniform policy no domain was imposed, so the sweep's
         * unowned intervals are exterior by the caller's declaration. */
        if (flag == ALEA_RAY_SLICE_DIAG_COVERAGE_GAP &&
            (coverage_flags & ALEA_RAY_SLICE_COVERAGE_UNOWNED_IS_EXTERIOR))
            flag = ALEA_RAY_SLICE_DIAG_COVERAGE_ALLOWED_EXTERIOR;
        if (flag == ALEA_RAY_SLICE_DIAG_COVERAGE_ALLOWED_EXTERIOR &&
            !(coverage_flags & ALEA_RAY_SLICE_COVERAGE_REPORT_EXTERIOR))
            flag = 0;
        *out_flags = flag;
        const size_t owners =
            coverage->owner_offsets[i + 1] - coverage->owner_offsets[i];
        *out_owners = owners > UINT32_MAX ? UINT32_MAX : (uint32_t)owners;
        return;
    }
}

static int add_coverage_boundaries(
    const alea_ray_coverage_slice_result_t* coverage, size_t begin, size_t end,
    double u_min, double u_max, double* boundaries, size_t* count,
    size_t capacity) {
    if (end - begin > (capacity - *count) / 2) return -1;
    for (size_t i = begin; i < end; i++) {
        double enter = u_min + coverage->t_enter[i];
        double exit = u_min + coverage->t_exit[i];
        if (enter < u_min) enter = u_min;
        if (enter > u_max) enter = u_max;
        if (exit < u_min) exit = u_min;
        if (exit > u_max) exit = u_max;
        boundaries[(*count)++] = enter;
        boundaries[(*count)++] = exit;
    }
    return 0;
}

/* Build one coverage row per validation row, matching the forward slice row
 * centres.  The domain is an explicit caller decision, never inferred from the
 * viewport extent. */
static int build_coverage_rows(const alea_slice_view_t* view, size_t row_count,
                               const alea_ray_slice_validation_options_t* options,
                               alea_ray_coverage_row_t** out_rows) {
    alea_ray_coverage_row_t* rows = calloc(row_count, sizeof(*rows));
    if (!rows) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate coverage rows");
        return -1;
    }
    const double span = view->u_max - view->u_min;
    const int has_domain =
        (options->coverage_flags & ALEA_RAY_SLICE_COVERAGE_HAS_DOMAIN) != 0;
    for (size_t row = 0; row < row_count; row++) {
        const double v = view->v_min + ((double)row + 0.5) *
            (view->v_max - view->v_min) / (double)row_count;
        double origin[3];
        for (size_t axis = 0; axis < 3; axis++)
            origin[axis] = view->plane.origin[axis] +
                view->u_min * view->plane.u_axis[axis] +
                v * view->plane.v_axis[axis];
        if (alea_ray_init(&rows[row].ray, origin[0], origin[1], origin[2],
                          view->plane.u_axis[0], view->plane.u_axis[1],
                          view->plane.u_axis[2]) != 0) {
            free(rows);
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "coverage row direction is degenerate");
            return -1;
        }
        rows[row].t_max = span;
        rows[row].direction_tag = 0;
        rows[row].transverse_coordinate = v;
        /* Without an explicit domain the caller has declared all unowned space
         * exterior, so no closed domain is imposed and nothing is a gap. */
        rows[row].use_domain = (uint8_t)has_domain;
        rows[row].domain.has_domain = (uint8_t)has_domain;
        rows[row].domain.t_min = options->coverage_domain_u_min - view->u_min;
        rows[row].domain.t_max = options->coverage_domain_u_max - view->u_min;
        rows[row].domain.report_allowed_exterior = (uint8_t)
            ((options->coverage_flags &
              ALEA_RAY_SLICE_COVERAGE_REPORT_EXTERIOR) != 0);
    }
    *out_rows = rows;
    return 0;
}

static int compare_double(const void* a, const void* b) {
    double da = *(const double*)a;
    double db = *(const double*)b;
    return (da > db) - (da < db);
}

static int view_owner_at(const trace_view_t* trace, size_t row, double u,
                         double tolerance, int32_t* out_cell,
                         uint64_t* out_key) {
    size_t begin = (size_t)trace->offsets[row];
    size_t end = (size_t)trace->offsets[row + 1];
    *out_cell = -1;
    *out_key = UINT64_MAX;
    for (size_t i = begin; i < end; i++) {
        double enter = trace->t_enter[i];
        double exit = trace->t_exit[i];
        if (trace->reverse) {
            double r_enter = trace->u_max - exit;
            double r_exit = trace->u_max - enter;
            enter = r_enter;
            exit = r_exit;
        }
        if (u > enter - tolerance && u < exit + tolerance) {
            *out_cell = trace->cells[i];
            if (trace->occurrence_keys && trace->occurrence_keys[i] != 0)
                *out_key = trace->occurrence_keys[i];
            return 0;
        }
    }
    return 0;
}

static int add_trace_boundaries(const trace_view_t* trace, size_t row,
                                double* boundaries, size_t* count,
                                size_t capacity) {
    size_t begin = (size_t)trace->offsets[row];
    size_t end = (size_t)trace->offsets[row + 1];
    if (end - begin > (capacity - *count) / 2) return -1;
    for (size_t i = begin; i < end; i++) {
        double enter = trace->t_enter[i];
        double exit = trace->t_exit[i];
        if (trace->reverse) {
            double r_enter = trace->u_max - exit;
            double r_exit = trace->u_max - enter;
            enter = r_enter;
            exit = r_exit;
        }
        if (enter < trace->u_min) enter = trace->u_min;
        if (enter > trace->u_max) enter = trace->u_max;
        if (exit < trace->u_min) exit = trace->u_min;
        if (exit > trace->u_max) exit = trace->u_max;
        boundaries[(*count)++] = enter;
        boundaries[(*count)++] = exit;
    }
    return 0;
}

static uint64_t selected_key(const trace_view_t* trace, size_t row, double u,
                             double tolerance, int32_t* cell) {
    uint64_t key;
    (void)view_owner_at(trace, row, u, tolerance, cell, &key);
    return key;
}

/* Find the canonical event evidence nearest a normalized U boundary.  The
 * event cache uses pixel-centre endpoints, so its local tolerance also covers
 * the small traversal epsilon used for synthetic lattice crossings. */
static int event_evidence_at_u(
    const alea_slice_directional_event_cache_t* cache, int width,
    size_t row, int reverse, double u, const alea_slice_view_t* view,
    double tolerance, int32_t* out_surface_id, uint32_t* out_flags) {
    const alea_ray_boundary_event_t* events = NULL;
    size_t count = 0;
    *out_surface_id = -1;
    *out_flags = 0;
    if (alea_slice_directional_event_cache_line_events(
            cache, ALEA_SLICE_EDGE_RIGHT, reverse, row, &events, &count) != 0)
        return -1;
    const double du = (view->u_max - view->u_min) / (double)width;
    const double origin_u = reverse ? view->u_max - 0.5 * du
                                    : view->u_min + 0.5 * du;
    const double event_tolerance = fmax(tolerance, fabs(du) * 1e-7);
    for (size_t i = 0; i < count; i++) {
        const double event_u = reverse ? origin_u - events[i].t
                                       : origin_u + events[i].t;
        if (fabs(event_u - u) > event_tolerance) continue;
        if (events[i].kind == ALEA_RAY_BOUNDARY_EVENT_SYNTHETIC_LATTICE) {
            *out_flags |= ALEA_RAY_SLICE_BOUNDARY_PROVENANCE_SYNTHETIC;
        } else if (events[i].kind == ALEA_RAY_BOUNDARY_EVENT_UNRESOLVED) {
            *out_flags |= ALEA_RAY_SLICE_BOUNDARY_PROVENANCE_UNRESOLVED;
        } else if (events[i].surface_id > 0) {
            if (*out_surface_id > 0)
                *out_flags |= ALEA_RAY_SLICE_BOUNDARY_PROVENANCE_COINCIDENT;
            if (*out_surface_id < 0 || events[i].surface_id < *out_surface_id)
                *out_surface_id = events[i].surface_id;
        }
    }
    return 0;
}

static int populate_boundary_provenance(
    alea_ray_slice_validation_result_t* candidate,
    const size_t* row_base_indices,
    const alea_slice_directional_event_cache_t* cache, int width,
    const alea_slice_view_t* view, double tolerance) {
    for (size_t row = 0; row < candidate->row_count; row++) {
        /* The event cache is defined on the requested viewport rows.  A refined
         * row has no cached line, so it publishes no boundary provenance rather
         * than borrowing a neighbour's. */
        const size_t base = row_base_indices[row];
        if (base == SIZE_MAX) {
            for (size_t i = candidate->row_offsets[row];
                 i < candidate->row_offsets[row + 1]; i++) {
                candidate->u_enter_forward_surface_ids[i] = -1;
                candidate->u_enter_reverse_surface_ids[i] = -1;
                candidate->u_exit_forward_surface_ids[i] = -1;
                candidate->u_exit_reverse_surface_ids[i] = -1;
                candidate->u_enter_provenance_flags[i] = 0;
                candidate->u_exit_provenance_flags[i] = 0;
            }
            continue;
        }
        for (size_t i = candidate->row_offsets[row];
             i < candidate->row_offsets[row + 1]; i++) {
            uint32_t forward_flags, reverse_flags;
            if (event_evidence_at_u(cache, width, base, 0, candidate->u_enter[i],
                                    view, tolerance,
                                    &candidate->u_enter_forward_surface_ids[i],
                                    &forward_flags) != 0 ||
                event_evidence_at_u(cache, width, base, 1, candidate->u_enter[i],
                                    view, tolerance,
                                    &candidate->u_enter_reverse_surface_ids[i],
                                    &reverse_flags) != 0)
                return -1;
            candidate->u_enter_provenance_flags[i] = forward_flags | reverse_flags;
            if (event_evidence_at_u(cache, width, base, 0, candidate->u_exit[i],
                                    view, tolerance,
                                    &candidate->u_exit_forward_surface_ids[i],
                                    &forward_flags) != 0 ||
                event_evidence_at_u(cache, width, base, 1, candidate->u_exit[i],
                                    view, tolerance,
                                    &candidate->u_exit_reverse_surface_ids[i],
                                    &reverse_flags) != 0)
                return -1;
            candidate->u_exit_provenance_flags[i] = forward_flags | reverse_flags;
        }
    }
    return 0;
}

void alea_ray_slice_validation_options_init(
    alea_ray_slice_validation_options_t* options) {
    if (!options) return;
    memset(options, 0, sizeof(*options));
    options->struct_size = sizeof(*options);
    options->projected_depth = -1;
}

alea_ray_slice_validation_result_t*
alea_ray_slice_validation_result_create(void) {
    return calloc(1, sizeof(alea_ray_slice_validation_result_t));
}

void alea_ray_slice_validation_result_destroy(
    alea_ray_slice_validation_result_t* result) {
    if (!result) return;
    validation_result_free_buffers(result);
    free(result);
}

size_t alea_ray_slice_validation_row_count(
    const alea_ray_slice_validation_result_t* result) {
    return result ? result->row_count : 0;
}
size_t alea_ray_slice_validation_interval_count(
    const alea_ray_slice_validation_result_t* result) {
    return result ? result->interval_count : 0;
}
uint32_t alea_ray_slice_validation_fields(
    const alea_ray_slice_validation_result_t* result) {
    return result ? result->fields : 0;
}
uint32_t alea_ray_slice_validation_executed_trace_mask(
    const alea_ray_slice_validation_result_t* result) {
    return result ? result->executed_trace_mask : 0;
}
uint32_t alea_ray_slice_validation_reused_trace_mask(
    const alea_ray_slice_validation_result_t* result) {
    return result ? result->reused_trace_mask : 0;
}
const uint64_t* alea_ray_slice_validation_row_offsets(
    const alea_ray_slice_validation_result_t* result) { return result ? result->row_offsets : NULL; }
const double* alea_ray_slice_validation_u_enter(
    const alea_ray_slice_validation_result_t* result) { return result ? result->u_enter : NULL; }
const double* alea_ray_slice_validation_u_exit(
    const alea_ray_slice_validation_result_t* result) { return result ? result->u_exit : NULL; }
const uint32_t* alea_ray_slice_validation_diagnostic_flags(
    const alea_ray_slice_validation_result_t* result) { return result ? result->diagnostic_flags : NULL; }
const int32_t* alea_ray_slice_validation_fast_forward_cell_ids(
    const alea_ray_slice_validation_result_t* result) { return result ? result->fast_forward_cell_ids : NULL; }
const int32_t* alea_ray_slice_validation_fast_reverse_cell_ids(
    const alea_ray_slice_validation_result_t* result) { return result ? result->fast_reverse_cell_ids : NULL; }
const uint64_t* alea_ray_slice_validation_fast_forward_occurrence_keys(
    const alea_ray_slice_validation_result_t* result) { return result ? result->fast_forward_occurrence_keys : NULL; }
const uint64_t* alea_ray_slice_validation_fast_reverse_occurrence_keys(
    const alea_ray_slice_validation_result_t* result) { return result ? result->fast_reverse_occurrence_keys : NULL; }
const uint32_t* alea_ray_slice_validation_coverage_owner_counts(
    const alea_ray_slice_validation_result_t* result) { return result ? result->coverage_owner_counts : NULL; }
const double* alea_ray_slice_validation_row_transverse_coordinates(
    const alea_ray_slice_validation_result_t* result) { return result ? result->row_transverse_coordinates : NULL; }
const uint8_t* alea_ray_slice_validation_row_direction_tags(
    const alea_ray_slice_validation_result_t* result) { return result ? result->row_direction_tags : NULL; }
const size_t* alea_ray_slice_validation_row_base_indices(
    const alea_ray_slice_validation_result_t* result) { return result ? result->row_base_indices : NULL; }
alea_ray_slice_refinement_status_t alea_ray_slice_validation_refinement_status(
    const alea_ray_slice_validation_result_t* result) {
    return result ? result->refinement_status
                  : ALEA_RAY_SLICE_REFINEMENT_NOT_REQUESTED;
}
const int32_t* alea_ray_slice_validation_u_enter_forward_surface_ids(
    const alea_ray_slice_validation_result_t* result) { return result ? result->u_enter_forward_surface_ids : NULL; }
const int32_t* alea_ray_slice_validation_u_enter_reverse_surface_ids(
    const alea_ray_slice_validation_result_t* result) { return result ? result->u_enter_reverse_surface_ids : NULL; }
const int32_t* alea_ray_slice_validation_u_exit_forward_surface_ids(
    const alea_ray_slice_validation_result_t* result) { return result ? result->u_exit_forward_surface_ids : NULL; }
const int32_t* alea_ray_slice_validation_u_exit_reverse_surface_ids(
    const alea_ray_slice_validation_result_t* result) { return result ? result->u_exit_reverse_surface_ids : NULL; }
const uint32_t* alea_ray_slice_validation_u_enter_provenance_flags(
    const alea_ray_slice_validation_result_t* result) { return result ? result->u_enter_provenance_flags : NULL; }
const uint32_t* alea_ray_slice_validation_u_exit_provenance_flags(
    const alea_ray_slice_validation_result_t* result) { return result ? result->u_exit_provenance_flags : NULL; }
const int32_t* alea_ray_slice_validation_u_enter_forward_surface_ids_internal(
    const alea_ray_slice_validation_result_t* result) { return result ? result->u_enter_forward_surface_ids : NULL; }
const int32_t* alea_ray_slice_validation_u_enter_reverse_surface_ids_internal(
    const alea_ray_slice_validation_result_t* result) { return result ? result->u_enter_reverse_surface_ids : NULL; }
const int32_t* alea_ray_slice_validation_u_exit_forward_surface_ids_internal(
    const alea_ray_slice_validation_result_t* result) { return result ? result->u_exit_forward_surface_ids : NULL; }
const int32_t* alea_ray_slice_validation_u_exit_reverse_surface_ids_internal(
    const alea_ray_slice_validation_result_t* result) { return result ? result->u_exit_reverse_surface_ids : NULL; }
const uint32_t* alea_ray_slice_validation_u_enter_provenance_flags_internal(
    const alea_ray_slice_validation_result_t* result) { return result ? result->u_enter_provenance_flags : NULL; }
const uint32_t* alea_ray_slice_validation_u_exit_provenance_flags_internal(
    const alea_ray_slice_validation_result_t* result) { return result ? result->u_exit_provenance_flags : NULL; }

static int allocate_candidate(alea_ray_slice_validation_result_t* candidate,
                              validation_row_t* rows, size_t row_count,
                              const alea_ray_slice_validation_options_t* options,
                              int include_provenance, int include_coverage,
                              const double* row_coordinates,
                              const uint8_t* row_tags,
                              const size_t* row_base_indices) {
    const int include_rows = row_coordinates != NULL;
    size_t total = 0, bytes = 0;
    for (size_t row = 0; row < row_count; row++) {
        if (rows[row].count > SIZE_MAX - total) return -1;
        total += rows[row].count;
    }
    if (options->max_output_intervals && total > options->max_output_intervals) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW, "ray slice validation output interval limit exceeded");
        return -1;
    }
#define ADD_BYTES(count, type) do { \
    if ((count) > SIZE_MAX / sizeof(type) || bytes > SIZE_MAX - (count) * sizeof(type)) return -1; \
    bytes += (count) * sizeof(type); \
} while (0)
    ADD_BYTES(row_count + 1, uint64_t);
    ADD_BYTES(total, double); ADD_BYTES(total, double); ADD_BYTES(total, uint32_t);
    ADD_BYTES(total, int32_t); ADD_BYTES(total, int32_t);
    ADD_BYTES(total, uint64_t); ADD_BYTES(total, uint64_t);
    if (include_coverage) ADD_BYTES(total, uint32_t);
    if (include_rows) {
        ADD_BYTES(row_count, double); ADD_BYTES(row_count, uint8_t);
        ADD_BYTES(row_count, size_t);
    }
    if (include_provenance) {
        ADD_BYTES(total, int32_t); ADD_BYTES(total, int32_t);
        ADD_BYTES(total, int32_t); ADD_BYTES(total, int32_t);
        ADD_BYTES(total, uint32_t); ADD_BYTES(total, uint32_t);
    }
#undef ADD_BYTES
    if (options->max_output_bytes && bytes > options->max_output_bytes) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW, "ray slice validation output byte limit exceeded");
        return -1;
    }
    candidate->row_offsets = calloc(row_count + 1, sizeof(*candidate->row_offsets));
    candidate->u_enter = malloc(total ? total * sizeof(*candidate->u_enter) : 1);
    candidate->u_exit = malloc(total ? total * sizeof(*candidate->u_exit) : 1);
    candidate->diagnostic_flags = malloc(total ? total * sizeof(*candidate->diagnostic_flags) : 1);
    candidate->fast_forward_cell_ids = malloc(total ? total * sizeof(*candidate->fast_forward_cell_ids) : 1);
    candidate->fast_reverse_cell_ids = malloc(total ? total * sizeof(*candidate->fast_reverse_cell_ids) : 1);
    candidate->fast_forward_occurrence_keys = malloc(total ? total * sizeof(*candidate->fast_forward_occurrence_keys) : 1);
    candidate->fast_reverse_occurrence_keys = malloc(total ? total * sizeof(*candidate->fast_reverse_occurrence_keys) : 1);
    if (!candidate->row_offsets || !candidate->u_enter || !candidate->u_exit ||
        !candidate->diagnostic_flags || !candidate->fast_forward_cell_ids ||
        !candidate->fast_reverse_cell_ids || !candidate->fast_forward_occurrence_keys ||
        !candidate->fast_reverse_occurrence_keys) return -1;
    if (include_coverage) {
        candidate->coverage_owner_counts = malloc(total ? total * sizeof(*candidate->coverage_owner_counts) : 1);
        if (!candidate->coverage_owner_counts) return -1;
    }
    if (include_rows) {
        candidate->row_transverse_coordinates = malloc(row_count ? row_count * sizeof(*candidate->row_transverse_coordinates) : 1);
        candidate->row_direction_tags = malloc(row_count ? row_count * sizeof(*candidate->row_direction_tags) : 1);
        candidate->row_base_indices = malloc(row_count ? row_count * sizeof(*candidate->row_base_indices) : 1);
        if (!candidate->row_transverse_coordinates ||
            !candidate->row_direction_tags || !candidate->row_base_indices)
            return -1;
        memcpy(candidate->row_transverse_coordinates, row_coordinates,
               row_count * sizeof(*row_coordinates));
        memcpy(candidate->row_direction_tags, row_tags,
               row_count * sizeof(*row_tags));
        memcpy(candidate->row_base_indices, row_base_indices,
               row_count * sizeof(*row_base_indices));
    }
    if (include_provenance) {
        candidate->u_enter_forward_surface_ids = malloc(total ? total * sizeof(*candidate->u_enter_forward_surface_ids) : 1);
        candidate->u_enter_reverse_surface_ids = malloc(total ? total * sizeof(*candidate->u_enter_reverse_surface_ids) : 1);
        candidate->u_exit_forward_surface_ids = malloc(total ? total * sizeof(*candidate->u_exit_forward_surface_ids) : 1);
        candidate->u_exit_reverse_surface_ids = malloc(total ? total * sizeof(*candidate->u_exit_reverse_surface_ids) : 1);
        candidate->u_enter_provenance_flags = malloc(total ? total * sizeof(*candidate->u_enter_provenance_flags) : 1);
        candidate->u_exit_provenance_flags = malloc(total ? total * sizeof(*candidate->u_exit_provenance_flags) : 1);
        if (!candidate->u_enter_forward_surface_ids ||
            !candidate->u_enter_reverse_surface_ids ||
            !candidate->u_exit_forward_surface_ids ||
            !candidate->u_exit_reverse_surface_ids ||
            !candidate->u_enter_provenance_flags ||
            !candidate->u_exit_provenance_flags) return -1;
    }
    candidate->row_count = row_count;
    candidate->interval_count = total;
    candidate->fields = ALEA_RAY_SLICE_VALIDATION_FIELD_FAST_FORWARD |
                        ALEA_RAY_SLICE_VALIDATION_FIELD_FAST_REVERSE |
                        ALEA_RAY_SLICE_VALIDATION_FIELD_OCCURRENCE_KEYS;
    if (include_coverage)
        candidate->fields |= ALEA_RAY_SLICE_VALIDATION_FIELD_COVERAGE;
    if (include_rows)
        candidate->fields |= ALEA_RAY_SLICE_VALIDATION_FIELD_ADAPTIVE_ROWS;
    size_t dst = 0;
    for (size_t row = 0; row < row_count; row++) {
        candidate->row_offsets[row] = dst;
        for (size_t i = 0; i < rows[row].count; i++, dst++) {
            candidate->u_enter[dst] = rows[row].u_enter[i];
            candidate->u_exit[dst] = rows[row].u_exit[i];
            candidate->diagnostic_flags[dst] = rows[row].flags[i];
            candidate->fast_forward_cell_ids[dst] = rows[row].forward_cells[i];
            candidate->fast_reverse_cell_ids[dst] = rows[row].reverse_cells[i];
            candidate->fast_forward_occurrence_keys[dst] = rows[row].forward_keys[i];
            candidate->fast_reverse_occurrence_keys[dst] = rows[row].reverse_keys[i];
            if (include_coverage)
                candidate->coverage_owner_counts[dst] = rows[row].coverage_owners[i];
        }
    }
    candidate->row_offsets[row_count] = dst;
    return 0;
}

int alea_validate_ray_slice_compact_with_event_cache(
    alea_system_t* sys, const alea_slice_view_t* view, size_t row_count,
    const alea_ray_slice_validation_options_t* validation_options,
    const alea_raycast_batch_options_t* render_options,
    alea_raycast_batch_result_t* inout_fast_forward,
    const alea_slice_directional_event_cache_t* event_cache,
    alea_ray_slice_validation_result_t* out_validation) {
    alea_ray_slice_validation_options_t defaults, options;
    alea_raycast_batch_options_t forward_options, reverse_options;
    alea_raycast_batch_result_t *fresh_forward = NULL, *fresh_reverse = NULL;
    const alea_raycast_batch_result_t* reverse = NULL;
    const alea_raycast_batch_result_t* forward = NULL;
    validation_row_t* rows = NULL;
    alea_ray_slice_validation_result_t candidate = {0};
    double *origins = NULL, *directions = NULL;
    alea_ray_coverage_row_t* coverage_rows = NULL;
    alea_ray_coverage_slice_result_t coverage_result;
    alea_ray_coverage_executor_t coverage_executor;
    const alea_ray_coverage_slice_result_t* coverage = NULL;
    size_t* row_base_indices = NULL;
    uint8_t* row_tags = NULL;
    const double* row_coordinates = NULL;
    size_t out_row_count = row_count;
    int coverage_initialized = 0;
    int coverage_executor_initialized = 0;
    int event_cache_width = 0, event_cache_height = 0;
    int cache_hit = 0, rc = -1;

    if (!sys || !view || !out_validation) {
        alea_set_error_detail(ALEA_ERR_NULL_ARG, "ray slice validation requires system, view, and result");
        return -1;
    }
    if (row_count == 0) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG, "ray slice validation row_count must be positive");
        return -1;
    }
    if (validation_options) {
        if (validation_options->struct_size < sizeof(*validation_options)) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG, "ray slice validation options struct_size is too small");
            return -1;
        }
        options = *validation_options;
    } else {
        alea_ray_slice_validation_options_init(&defaults);
        options = defaults;
    }
    if (options.checks & ~(uint32_t)(ALEA_RAY_SLICE_VALIDATE_FAST_BIDIRECTIONAL |
                                     ALEA_RAY_SLICE_VALIDATE_COVERAGE) ||
        options.flags & ~ALEA_RAY_SLICE_VALIDATION_INCLUDE_AGREEMENTS ||
        options.coverage_flags & ~(uint32_t)(ALEA_RAY_SLICE_COVERAGE_HAS_DOMAIN |
                                             ALEA_RAY_SLICE_COVERAGE_UNOWNED_IS_EXTERIOR |
                                             ALEA_RAY_SLICE_COVERAGE_REPORT_EXTERIOR)) {
        alea_set_error_detail(ALEA_ERR_UNSUPPORTED, "unsupported ray slice validation option");
        return -1;
    }
    if (options.checks & ALEA_RAY_SLICE_VALIDATE_COVERAGE) {
        /* Unowned space cannot be classified without an explicit decision
         * about where ownership is required.  The viewport extent is an
         * observation window, not a validation domain. */
        const int has_domain =
            (options.coverage_flags & ALEA_RAY_SLICE_COVERAGE_HAS_DOMAIN) != 0;
        const int uniform =
            (options.coverage_flags & ALEA_RAY_SLICE_COVERAGE_UNOWNED_IS_EXTERIOR) != 0;
        if (has_domain == uniform) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "coverage validation requires exactly one of an explicit "
                                  "domain or a uniform unowned-space policy");
            return -1;
        }
        if (has_domain &&
            !(options.coverage_domain_u_min < options.coverage_domain_u_max)) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "coverage validation domain must be a non-empty U range");
            return -1;
        }
        if (options.coverage_refine_signals &
            ~(uint32_t)(ALEA_RAY_SLICE_REFINE_SIGNATURE |
                        ALEA_RAY_SLICE_REFINE_DISPLACEMENT |
                        ALEA_RAY_SLICE_REFINE_DENSITY |
                        ALEA_RAY_SLICE_REFINE_FINDING)) {
            alea_set_error_detail(ALEA_ERR_UNSUPPORTED,
                                  "unsupported coverage refinement signal");
            return -1;
        }
    } else if (options.coverage_max_refinement_depth != 0) {
        /* Refinement drives coverage probes; without the coverage check there
         * is nothing for a refined row to measure. */
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "adaptive refinement requires the coverage check");
        return -1;
    }
    if (options.checks == 0 && !inout_fast_forward) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG, "empty validation requires a forward cache destination");
        return -1;
    }
    if (render_options && render_options->struct_size < sizeof(*render_options)) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG, "render options struct_size is too small");
        return -1;
    }
    if (event_cache) {
        if (row_count > INT_MAX) {
            alea_set_error_detail(ALEA_ERR_OVERFLOW,
                                  "validation rows exceed event-cache range");
            return -1;
        }
        if (alea_slice_directional_event_cache_dimensions(
                event_cache, &event_cache_width, &event_cache_height) != 0 ||
            event_cache_height != (int)row_count ||
            !alea_slice_directional_event_cache_matches(
                event_cache, sys, view, event_cache_width, event_cache_height)) {
            alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                  "event cache does not match validation view or rows");
            return -1;
        }
    }

    /* The directional cache owns both U traces.  Prefer a caller-supplied
     * compatible forward trace (it may carry richer requested fields), then
     * fall back to the cache's ownership trace. */
    cache_hit = inout_fast_forward &&
        alea_raycast_batch_result_matches_fast_slice_cache(
            inout_fast_forward, sys, view, row_count, render_options,
            options.projected_depth);
    if (cache_hit && options.max_trace_intervals &&
        alea_raycast_batch_segment_count(inout_fast_forward) > options.max_trace_intervals)
        cache_hit = 0;

    if (!cache_hit && event_cache && options.projected_depth == -1) {
        const alea_raycast_batch_result_t* cached =
            alea_slice_directional_event_cache_ownership_trace(
                event_cache, ALEA_SLICE_EDGE_RIGHT, 0, -1,
                ALEA_RAY_BATCH_PROJECTED_OWNER);
        if (cached && (!options.max_trace_intervals ||
                       alea_raycast_batch_segment_count(cached) <= options.max_trace_intervals)) {
            forward = cached;
            cache_hit = 1;
            candidate.reused_trace_mask |= ALEA_RAY_SLICE_TRACE_FAST_FORWARD;
        }
    }

    if (cache_hit) {
        if (!forward) {
            forward = inout_fast_forward;
            candidate.reused_trace_mask |= ALEA_RAY_SLICE_TRACE_FAST_FORWARD;
        }
    } else {
        memset(&forward_options, 0, sizeof(forward_options));
        if (render_options) forward_options = *render_options;
        forward_options.struct_size = sizeof(forward_options);
        forward_options.fields |= ALEA_RAY_BATCH_PROJECTED_OWNER;
        forward_options.projected_depth = options.projected_depth;
        if (options.max_trace_intervals &&
            (!forward_options.max_segments || forward_options.max_segments > options.max_trace_intervals))
            forward_options.max_segments = options.max_trace_intervals;
        fresh_forward = alea_raycast_batch_result_create();
        if (!fresh_forward || alea_trace_ray_slice_compact(sys, view, row_count,
                                                           &forward_options, fresh_forward) != 0)
            goto cleanup;
        forward = fresh_forward;
        candidate.executed_trace_mask |= ALEA_RAY_SLICE_TRACE_FAST_FORWARD;
    }

    if (options.max_trace_intervals &&
        alea_raycast_batch_segment_count(forward) > options.max_trace_intervals) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW, "ray slice validation trace interval limit exceeded");
        goto cleanup;
    }

    if (options.checks & ALEA_RAY_SLICE_VALIDATE_FAST_BIDIRECTIONAL) {
        size_t forward_count = alea_raycast_batch_segment_count(forward);
        uint64_t remaining = 0;
        if (options.max_trace_intervals) {
            if (forward_count >= options.max_trace_intervals) {
                alea_set_error_detail(ALEA_ERR_OVERFLOW, "ray slice validation trace interval limit exceeded");
                goto cleanup;
            }
            remaining = options.max_trace_intervals - forward_count;
        }
        if (row_count > SIZE_MAX / (3 * sizeof(double))) {
            alea_set_error_detail(ALEA_ERR_OVERFLOW, "ray slice validation row storage overflows");
            goto cleanup;
        }
        origins = malloc(row_count * 3 * sizeof(*origins));
        directions = malloc(row_count * 3 * sizeof(*directions));
        if (!origins || !directions) {
            alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "failed to allocate reverse slice rows");
            goto cleanup;
        }
        for (size_t row = 0; row < row_count; row++) {
            double v = view->v_min + ((double)row + 0.5) *
                       (view->v_max - view->v_min) / (double)row_count;
            for (size_t axis = 0; axis < 3; axis++) {
                origins[row * 3 + axis] = view->plane.origin[axis] +
                    view->u_max * view->plane.u_axis[axis] + v * view->plane.v_axis[axis];
                directions[row * 3 + axis] = -view->plane.u_axis[axis];
            }
        }
        memset(&reverse_options, 0, sizeof(reverse_options));
        reverse_options.struct_size = sizeof(reverse_options);
        reverse_options.fields = ALEA_RAY_BATCH_PROJECTED_OWNER;
        reverse_options.projected_depth = options.projected_depth;
        reverse_options.max_segments = remaining;
        if (event_cache && options.projected_depth == -1) {
            reverse = alea_slice_directional_event_cache_ownership_trace(
                event_cache, ALEA_SLICE_EDGE_RIGHT, 1, -1,
                ALEA_RAY_BATCH_PROJECTED_OWNER);
            if (reverse && (!options.max_trace_intervals ||
                            alea_raycast_batch_segment_count(reverse) <= remaining))
                candidate.reused_trace_mask |= ALEA_RAY_SLICE_TRACE_FAST_REVERSE;
            else
                reverse = NULL;
        }
        if (!reverse) {
            fresh_reverse = alea_raycast_batch_result_create();
            if (!fresh_reverse || alea_raycast_hier_batch(sys, origins, directions, row_count,
                                                          view->u_max - view->u_min,
                                                          &reverse_options, fresh_reverse) != 0)
                goto cleanup;
            reverse = fresh_reverse;
            candidate.executed_trace_mask |= ALEA_RAY_SLICE_TRACE_FAST_REVERSE;
        }
    }

    if (options.checks & ALEA_RAY_SLICE_VALIDATE_COVERAGE) {
        alea_ray_coverage_slice_limits_t coverage_limits;
        alea_ray_coverage_slice_limits_init(&coverage_limits);
        coverage_limits.max_owners = (size_t)options.max_coverage_owners;
        coverage_limits.max_intervals = (size_t)options.max_trace_intervals;
        coverage_limits.max_rows = (size_t)options.max_coverage_rows;
        if (build_coverage_rows(view, row_count, &options, &coverage_rows) != 0)
            goto cleanup;
        alea_ray_coverage_slice_result_init(&coverage_result);
        alea_ray_coverage_executor_init(&coverage_executor);
        size_t coverage_workers = 1;
#ifdef _OPENMP
        coverage_workers = (size_t)omp_get_max_threads();
#endif
        if (alea_ray_coverage_executor_prepare(&coverage_executor,
                                               coverage_workers) != 0)
            goto cleanup;
        coverage_executor_initialized = 1;
        coverage_initialized = 1;
        if (options.coverage_max_refinement_depth == 0) {
            if (alea_ray_coverage_slice_build_executor_nocache(
                    sys, coverage_rows, row_count, &coverage_limits,
                    &coverage_executor, &coverage_result) != 0) {
                if (alea_error_code() == ALEA_OK)
                    alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                                          "coverage slice classification failed");
                goto cleanup;
            }
        } else {
            alea_ray_coverage_refinement_policy_t policy;
            alea_ray_coverage_refinement_policy_init(&policy);
            if (options.coverage_refine_signals != 0)
                policy.signals = options.coverage_refine_signals;
            policy.endpoint_displacement = options.coverage_endpoint_displacement;
            policy.crossing_density = (size_t)options.coverage_crossing_density;
            policy.min_transverse_spacing = options.coverage_min_transverse_spacing;
            if (alea_ray_coverage_slice_build_adaptive_policy_executor_nocache(
                sys, coverage_rows, row_count,
                options.coverage_max_refinement_depth, &policy,
                &coverage_limits, &coverage_executor,
                &coverage_result) != 0) {
                if (alea_error_code() == ALEA_OK)
                    alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                                          "coverage refinement policy is inconsistent");
                goto cleanup;
            }
            switch (coverage_result.refinement_status) {
            case ALEA_RAY_COVERAGE_REFINEMENT_COMPLETE:
                candidate.refinement_status = ALEA_RAY_SLICE_REFINEMENT_CONVERGED;
                break;
            case ALEA_RAY_COVERAGE_REFINEMENT_MAX_DEPTH:
                candidate.refinement_status = ALEA_RAY_SLICE_REFINEMENT_MAX_DEPTH;
                break;
            case ALEA_RAY_COVERAGE_REFINEMENT_MAX_ROWS:
                candidate.refinement_status = ALEA_RAY_SLICE_REFINEMENT_MAX_ROWS;
                break;
            case ALEA_RAY_COVERAGE_REFINEMENT_MIN_SPACING:
                candidate.refinement_status = ALEA_RAY_SLICE_REFINEMENT_MIN_SPACING;
                break;
            }
        }
        coverage = &coverage_result;
        out_row_count = coverage_result.row_count;
        if (out_row_count < row_count) {
            alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                                  "coverage refinement dropped requested rows");
            goto cleanup;
        }
    }

    /* Map published rows back to the requested viewport rows.  Refined rows sit
     * strictly between requested transverse coordinates, so a single ordered
     * walk separates them; a refined row keeps SIZE_MAX and carries coverage
     * evidence only. */
    row_base_indices = malloc(out_row_count * sizeof(*row_base_indices));
    if (!row_base_indices) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate validation row provenance");
        goto cleanup;
    }
    if (out_row_count == row_count) {
        for (size_t row = 0; row < row_count; row++) row_base_indices[row] = row;
    } else {
        size_t base = 0;
        for (size_t row = 0; row < out_row_count; row++) {
            if (base < row_count &&
                coverage->row_transverse_coordinates[row] ==
                    coverage_rows[base].transverse_coordinate)
                row_base_indices[row] = base++;
            else
                row_base_indices[row] = SIZE_MAX;
        }
        if (base != row_count) {
            alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                                  "coverage refinement lost a requested row coordinate");
            goto cleanup;
        }
    }

    rows = calloc(out_row_count, sizeof(*rows));
    if (!rows) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "failed to allocate validation rows");
        goto cleanup;
    }
    if (reverse || coverage) {
        trace_view_t fwd = {
            .result = forward,
            .offsets = alea_raycast_batch_ray_offsets(forward),
            .t_enter = alea_raycast_batch_t_enter(forward),
            .t_exit = alea_raycast_batch_t_exit(forward),
            .cells = options.projected_depth >= 0 ?
                alea_raycast_batch_projected_cell_ids(forward) :
                alea_raycast_batch_cell_ids(forward),
            .occurrence_keys = alea_raycast_batch_projected_occurrence_keys(forward),
            .u_min = view->u_min, .u_max = view->u_max
        };
        trace_view_t rev = {
            .result = reverse,
            .offsets = alea_raycast_batch_ray_offsets(reverse),
            .t_enter = alea_raycast_batch_t_enter(reverse),
            .t_exit = alea_raycast_batch_t_exit(reverse),
            .cells = options.projected_depth >= 0 ?
                alea_raycast_batch_projected_cell_ids(reverse) :
                alea_raycast_batch_cell_ids(reverse),
            .occurrence_keys = alea_raycast_batch_projected_occurrence_keys(reverse),
            .reverse = 1, .u_min = view->u_min, .u_max = view->u_max
        };
        if (!fwd.offsets || !fwd.t_enter || !fwd.t_exit || !fwd.cells ||
            (reverse && (!rev.offsets || !rev.t_enter || !rev.t_exit || !rev.cells))) {
            alea_set_error_detail(ALEA_ERR_INVALID_STATE, "ray slice validation trace lacks required fields");
            goto cleanup;
        }
        double tolerance = options.absolute_tolerance > 0.0 ? options.absolute_tolerance : 1e-9;
        double relative = options.relative_tolerance > 0.0 ? options.relative_tolerance : 1e-12;
        double span = view->u_max - view->u_min;
        if (relative * fmax(1.0, span) > tolerance)
            tolerance = relative * fmax(1.0, span);
        for (size_t row = 0; row < out_row_count; row++) {
            /* Selected traces exist only for the requested viewport rows; a
             * refined row contributes coverage evidence alone. */
            const size_t base = row_base_indices[row];
            const int has_selected = base != SIZE_MAX;
            size_t fc = has_selected ?
                (size_t)(fwd.offsets[base + 1] - fwd.offsets[base]) : 0;
            size_t rcnt = (has_selected && reverse) ?
                (size_t)(rev.offsets[base + 1] - rev.offsets[base]) : 0;
            /* Coverage breakpoints partition the row alongside the selected
             * traces: a coverage boundary need not be a selected transition. */
            size_t coverage_begin = 0, coverage_end = 0;
            if (coverage &&
                coverage_row_span(coverage, row, &coverage_begin,
                                  &coverage_end) != 0) {
                alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                                      "coverage slice row is malformed");
                goto cleanup;
            }
            const size_t ccnt = coverage_end - coverage_begin;
            if (fc > (SIZE_MAX - 2) / 2 || rcnt > (SIZE_MAX - 2 - 2 * fc) / 2 ||
                ccnt > (SIZE_MAX - 2 - 2 * (fc + rcnt)) / 2) {
                alea_set_error_detail(ALEA_ERR_OVERFLOW, "ray slice validation boundaries overflow");
                goto cleanup;
            }
            size_t cap = 2 + 2 * (fc + rcnt + ccnt);
            double* boundaries = malloc(cap * sizeof(*boundaries));
            if (!boundaries) {
                alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "failed to allocate validation boundaries");
                goto cleanup;
            }
            size_t count = 0;
            boundaries[count++] = view->u_min;
            boundaries[count++] = view->u_max;
            if ((has_selected &&
                 add_trace_boundaries(&fwd, base, boundaries, &count, cap) != 0) ||
                (has_selected && reverse &&
                 add_trace_boundaries(&rev, base, boundaries, &count, cap) != 0) ||
                (coverage &&
                 add_coverage_boundaries(coverage, coverage_begin, coverage_end,
                                         view->u_min, view->u_max, boundaries,
                                         &count, cap) != 0)) {
                free(boundaries);
                alea_set_error_detail(ALEA_ERR_OVERFLOW, "ray slice validation boundary storage overflow");
                goto cleanup;
            }
            qsort(boundaries, count, sizeof(*boundaries), compare_double);
            size_t unique = 0;
            for (size_t i = 0; i < count; i++) {
                if (unique == 0 || boundaries[i] - boundaries[unique - 1] > tolerance)
                    boundaries[unique++] = boundaries[i];
            }
            for (size_t i = 0; i + 1 < unique; i++) {
                double a = boundaries[i], b = boundaries[i + 1];
                if (b - a <= tolerance) continue;
                double mid = a + 0.5 * (b - a);
                int32_t fcell = -1, rcell = -1;
                uint64_t fkey = UINT64_MAX, rkey = UINT64_MAX;
                uint32_t flags = 0;
                if (reverse && has_selected) {
                    fkey = selected_key(&fwd, base, mid, tolerance, &fcell);
                    rkey = selected_key(&rev, base, mid, tolerance, &rcell);
                    const int mismatch = (fkey != UINT64_MAX && rkey != UINT64_MAX) ?
                        fkey != rkey : fcell != rcell;
                    /* Directional disagreement is trace-consistency evidence
                     * and stays in its own flag; it is never merged into a
                     * coverage finding. */
                    if (mismatch)
                        flags |= ALEA_RAY_SLICE_DIAG_FAST_DIRECTION_MISMATCH;
                }
                uint32_t coverage_owners = 0;
                if (coverage) {
                    uint32_t coverage_diag = 0;
                    coverage_flags_at_u(coverage, coverage_begin, coverage_end,
                                        view->u_min, mid, options.coverage_flags,
                                        &coverage_diag, &coverage_owners);
                    flags |= coverage_diag;
                }
                if (!flags && !(options.flags & ALEA_RAY_SLICE_VALIDATION_INCLUDE_AGREEMENTS))
                    continue;
                if (rows[row].count > 0) {
                    size_t prev = rows[row].count - 1;
                    if (fabs(rows[row].u_exit[prev] - a) <= tolerance &&
                        rows[row].flags[prev] == flags &&
                        rows[row].forward_cells[prev] == fcell &&
                        rows[row].reverse_cells[prev] == rcell &&
                        rows[row].forward_keys[prev] == fkey &&
                        rows[row].reverse_keys[prev] == rkey &&
                        rows[row].coverage_owners[prev] == coverage_owners) {
                        rows[row].u_exit[prev] = b;
                        continue;
                    }
                }
                if (validation_row_append(&rows[row], a, b, flags, fcell, rcell,
                                          fkey, rkey, coverage_owners) != 0) {
                    free(boundaries);
                    alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "failed to append validation interval");
                    goto cleanup;
                }
            }
            free(boundaries);
        }
    }
    /* Row provenance is published whenever refinement was requested, so the
     * caller never has to infer which published rows it asked for. */
    if (options.coverage_max_refinement_depth != 0) {
        row_tags = calloc(out_row_count, sizeof(*row_tags));
        if (!row_tags) {
            alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                                  "failed to allocate validation row tags");
            goto cleanup;
        }
        for (size_t row = 0; row < out_row_count; row++)
            row_tags[row] = coverage->row_direction_tags[row];
        row_coordinates = coverage->row_transverse_coordinates;
    }
    if (allocate_candidate(&candidate, rows, out_row_count, &options,
                           event_cache != NULL, coverage != NULL,
                           row_coordinates, row_tags, row_base_indices) != 0) {
        if (alea_error_code() == ALEA_OK)
            alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "failed to allocate validation result");
        goto cleanup;
    }
    if (event_cache) {
        double tolerance = options.absolute_tolerance > 0.0 ?
            options.absolute_tolerance : 1e-9;
        if (populate_boundary_provenance(&candidate, row_base_indices,
                                         event_cache,
                                         event_cache_width, view,
                                         tolerance) != 0) {
            alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                                  "failed to read matching event cache");
            goto cleanup;
        }
    }
    validation_result_free_buffers(out_validation);
    *out_validation = candidate;
    memset(&candidate, 0, sizeof(candidate));
    if (fresh_forward && inout_fast_forward) {
        alea_raycast_batch_result_swap_internal(inout_fast_forward, fresh_forward);
    }
    rc = 0;

cleanup:
    free(origins);
    free(directions);
    free(coverage_rows);
    free(row_base_indices);
    free(row_tags);
    if (coverage_initialized) {
        alea_ray_coverage_slice_result_free(&coverage_result);
    }
    if (coverage_executor_initialized)
        alea_ray_coverage_executor_free(&coverage_executor);
    if (rows) {
        for (size_t row = 0; row < out_row_count; row++) validation_row_free(&rows[row]);
        free(rows);
    }
    validation_result_free_buffers(&candidate);
    alea_raycast_batch_result_destroy(fresh_reverse);
    alea_raycast_batch_result_destroy(fresh_forward);
    return rc;
}

int alea_validate_ray_slice_compact(
    alea_system_t* sys, const alea_slice_view_t* view, size_t row_count,
    const alea_ray_slice_validation_options_t* validation_options,
    const alea_raycast_batch_options_t* render_options,
    alea_raycast_batch_result_t* inout_fast_forward,
    alea_ray_slice_validation_result_t* out_validation) {
    return alea_validate_ray_slice_compact_with_event_cache(
        sys, view, row_count, validation_options, render_options,
        inout_fast_forward, NULL, out_validation);
}

int alea_validate_ray_slice_compact_with_directional_cache(
    alea_system_t* sys, const alea_slice_view_t* view, size_t row_count,
    const alea_ray_slice_validation_options_t* validation_options,
    const alea_raycast_batch_options_t* render_options,
    alea_raycast_batch_result_t* inout_fast_forward,
    const alea_slice_directional_trace_cache_t* cache,
    alea_ray_slice_validation_result_t* out_validation) {
    return alea_validate_ray_slice_compact_with_event_cache(
        sys, view, row_count, validation_options, render_options,
        inout_fast_forward, cache, out_validation);
}
