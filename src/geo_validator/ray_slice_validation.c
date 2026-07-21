// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_geo_validator.h"

#include "raycast/raycast.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

struct alea_ray_slice_validation_result {
    size_t row_count;
    size_t interval_count;
    uint32_t fields;
    uint32_t executed_trace_mask;
    uint32_t reused_trace_mask;
    uint64_t* row_offsets;
    double* u_enter;
    double* u_exit;
    uint32_t* diagnostic_flags;
    int32_t* fast_forward_cell_ids;
    int32_t* fast_reverse_cell_ids;
    uint64_t* fast_forward_occurrence_keys;
    uint64_t* fast_reverse_occurrence_keys;
};

typedef struct {
    double* u_enter;
    double* u_exit;
    uint32_t* flags;
    int32_t* forward_cells;
    int32_t* reverse_cells;
    uint64_t* forward_keys;
    uint64_t* reverse_keys;
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
    free(result->u_enter);
    free(result->u_exit);
    free(result->diagnostic_flags);
    free(result->fast_forward_cell_ids);
    free(result->fast_reverse_cell_ids);
    free(result->fast_forward_occurrence_keys);
    free(result->fast_reverse_occurrence_keys);
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
#undef GROW_ROW_FIELD
    row->capacity = capacity;
    return 0;
}

static int validation_row_append(validation_row_t* row,
                                 double u_enter, double u_exit,
                                 uint32_t flags, int32_t forward_cell,
                                 int32_t reverse_cell, uint64_t forward_key,
                                 uint64_t reverse_key) {
    if (validation_row_reserve(row, row->count + 1) != 0) return -1;
    size_t i = row->count++;
    row->u_enter[i] = u_enter;
    row->u_exit[i] = u_exit;
    row->flags[i] = flags;
    row->forward_cells[i] = forward_cell;
    row->reverse_cells[i] = reverse_cell;
    row->forward_keys[i] = forward_key;
    row->reverse_keys[i] = reverse_key;
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

static int allocate_candidate(alea_ray_slice_validation_result_t* candidate,
                              validation_row_t* rows, size_t row_count,
                              const alea_ray_slice_validation_options_t* options) {
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
    candidate->row_count = row_count;
    candidate->interval_count = total;
    candidate->fields = ALEA_RAY_SLICE_VALIDATION_FIELD_FAST_FORWARD |
                        ALEA_RAY_SLICE_VALIDATION_FIELD_FAST_REVERSE |
                        ALEA_RAY_SLICE_VALIDATION_FIELD_OCCURRENCE_KEYS;
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
        }
    }
    candidate->row_offsets[row_count] = dst;
    return 0;
}

int alea_validate_ray_slice_compact(
    alea_system_t* sys, const alea_slice_view_t* view, size_t row_count,
    const alea_ray_slice_validation_options_t* validation_options,
    const alea_raycast_batch_options_t* render_options,
    alea_raycast_batch_result_t* inout_fast_forward,
    alea_ray_slice_validation_result_t* out_validation) {
    alea_ray_slice_validation_options_t defaults, options;
    alea_raycast_batch_options_t forward_options, reverse_options;
    alea_raycast_batch_result_t *fresh_forward = NULL, *reverse = NULL;
    const alea_raycast_batch_result_t* forward;
    validation_row_t* rows = NULL;
    alea_ray_slice_validation_result_t candidate = {0};
    double *origins = NULL, *directions = NULL;
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
    if (options.checks & ~ALEA_RAY_SLICE_VALIDATE_FAST_BIDIRECTIONAL ||
        options.flags & ~ALEA_RAY_SLICE_VALIDATION_INCLUDE_AGREEMENTS) {
        alea_set_error_detail(ALEA_ERR_UNSUPPORTED, "unsupported ray slice validation option");
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

    cache_hit = inout_fast_forward &&
        alea_raycast_batch_result_matches_fast_slice_cache(
            inout_fast_forward, sys, view, row_count, render_options,
            options.projected_depth);
    if (cache_hit && options.max_trace_intervals &&
        alea_raycast_batch_segment_count(inout_fast_forward) > options.max_trace_intervals)
        cache_hit = 0;

    if (cache_hit) {
        forward = inout_fast_forward;
        candidate.reused_trace_mask |= ALEA_RAY_SLICE_TRACE_FAST_FORWARD;
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
        reverse = alea_raycast_batch_result_create();
        if (!reverse || alea_raycast_hier_batch(sys, origins, directions, row_count,
                                                view->u_max - view->u_min,
                                                &reverse_options, reverse) != 0)
            goto cleanup;
        candidate.executed_trace_mask |= ALEA_RAY_SLICE_TRACE_FAST_REVERSE;
    }

    rows = calloc(row_count, sizeof(*rows));
    if (!rows) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "failed to allocate validation rows");
        goto cleanup;
    }
    if (reverse) {
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
            !rev.offsets || !rev.t_enter || !rev.t_exit || !rev.cells) {
            alea_set_error_detail(ALEA_ERR_INVALID_STATE, "ray slice validation trace lacks required fields");
            goto cleanup;
        }
        double tolerance = options.absolute_tolerance > 0.0 ? options.absolute_tolerance : 1e-9;
        double relative = options.relative_tolerance > 0.0 ? options.relative_tolerance : 1e-12;
        double span = view->u_max - view->u_min;
        if (relative * fmax(1.0, span) > tolerance)
            tolerance = relative * fmax(1.0, span);
        for (size_t row = 0; row < row_count; row++) {
            size_t fc = (size_t)(fwd.offsets[row + 1] - fwd.offsets[row]);
            size_t rcnt = (size_t)(rev.offsets[row + 1] - rev.offsets[row]);
            if (fc > (SIZE_MAX - 2) / 2 || rcnt > (SIZE_MAX - 2 - 2 * fc) / 2) {
                alea_set_error_detail(ALEA_ERR_OVERFLOW, "ray slice validation boundaries overflow");
                goto cleanup;
            }
            size_t cap = 2 + 2 * (fc + rcnt);
            double* boundaries = malloc(cap * sizeof(*boundaries));
            if (!boundaries) {
                alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "failed to allocate validation boundaries");
                goto cleanup;
            }
            size_t count = 0;
            boundaries[count++] = view->u_min;
            boundaries[count++] = view->u_max;
            if (add_trace_boundaries(&fwd, row, boundaries, &count, cap) != 0 ||
                add_trace_boundaries(&rev, row, boundaries, &count, cap) != 0) {
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
                int32_t fcell, rcell;
                uint64_t fkey = selected_key(&fwd, row, mid, tolerance, &fcell);
                uint64_t rkey = selected_key(&rev, row, mid, tolerance, &rcell);
                int mismatch = (fkey != UINT64_MAX && rkey != UINT64_MAX) ?
                    fkey != rkey : fcell != rcell;
                uint32_t flags = mismatch ? ALEA_RAY_SLICE_DIAG_FAST_DIRECTION_MISMATCH : 0;
                if (!flags && !(options.flags & ALEA_RAY_SLICE_VALIDATION_INCLUDE_AGREEMENTS))
                    continue;
                if (rows[row].count > 0) {
                    size_t prev = rows[row].count - 1;
                    if (fabs(rows[row].u_exit[prev] - a) <= tolerance &&
                        rows[row].flags[prev] == flags &&
                        rows[row].forward_cells[prev] == fcell &&
                        rows[row].reverse_cells[prev] == rcell &&
                        rows[row].forward_keys[prev] == fkey &&
                        rows[row].reverse_keys[prev] == rkey) {
                        rows[row].u_exit[prev] = b;
                        continue;
                    }
                }
                if (validation_row_append(&rows[row], a, b, flags, fcell, rcell,
                                          fkey, rkey) != 0) {
                    free(boundaries);
                    alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "failed to append validation interval");
                    goto cleanup;
                }
            }
            free(boundaries);
        }
    }
    if (allocate_candidate(&candidate, rows, row_count, &options) != 0) {
        if (alea_error_code() == ALEA_OK)
            alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY, "failed to allocate validation result");
        goto cleanup;
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
    if (rows) {
        for (size_t row = 0; row < row_count; row++) validation_row_free(&rows[row]);
        free(rows);
    }
    validation_result_free_buffers(&candidate);
    alea_raycast_batch_result_destroy(reverse);
    alea_raycast_batch_result_destroy(fresh_forward);
    return rc;
}
