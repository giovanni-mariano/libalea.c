// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file ray_coverage.c
 * @brief Complete-coverage diagnostic sweep for ray geometry queries.
 *
 * This intentionally does not use selected-owner segments.  A selected
 * trace is allowed to hide overlapping claimants; coverage answers whether
 * every open interval is uniquely owned.
 */

#include "raycast.h"
#include "core/alea_universe.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Keep one sentinel slot beyond the published diagnostic budget.  Recursive
 * all-owner resolution stops at its supplied capacity, so this is how the
 * sweep distinguishes a complete 32-owner set from a truncated one. */
enum {
    ALEA_RAY_COVERAGE_OWNER_BUDGET = 32,
    ALEA_RAY_COVERAGE_MAX_OWNERS = ALEA_RAY_COVERAGE_OWNER_BUDGET + 1
};
#define ALEA_RAY_COVERAGE_MIN_INTERVAL 1e-9

static int coverage_owner_sets_equal(const uint64_t* a,
                                     const uint64_t* a_parents, int a_count,
                                     const uint64_t* b,
                                     const uint64_t* b_parents, int b_count) {
    if (a_count != b_count) return 0;
    /* Cell ID is presentation data, not ownership identity.  The recursive
     * coverage resolver derives each key from the concrete fill/lattice path,
     * so repeated or transformed occurrences cannot collapse here. */
    for (int i = 0; i < a_count; i++) {
        if (a[i] != b[i] || a_parents[i] != b_parents[i])
            return 0;
    }
    return 1;
}

static int coverage_owner_index(const alea_ray_coverage_owner_t* owners,
                                size_t owner_count, uint64_t occurrence_key) {
    for (size_t i = 0; i < owner_count; i++)
        if (owners[i].occurrence_key == occurrence_key)
            return (int)i;
    return -1;
}

static alea_ray_coverage_kind_t coverage_kind(
    const alea_ray_coverage_owner_t* owners, size_t owner_count) {
    if (owner_count == 0) return ALEA_RAY_COVERAGE_GAP;

    size_t root_count = 0;
    size_t child_counts[ALEA_RAY_COVERAGE_MAX_OWNERS] = {0};
    int deepest = 0;
    for (size_t i = 0; i < owner_count; i++) {
        if (owners[i].depth > owners[deepest].depth)
            deepest = (int)i;
        if (owners[i].parent_occurrence_key == 0) {
            root_count++;
            continue;
        }
        const int parent = coverage_owner_index(
            owners, owner_count, owners[i].parent_occurrence_key);
        if (parent < 0 || owners[i].depth != owners[parent].depth + 1)
            return ALEA_RAY_COVERAGE_UNRESOLVED;
        child_counts[parent]++;
        if (child_counts[parent] > 1)
            return ALEA_RAY_COVERAGE_OVERLAP;
    }
    if (root_count != 1)
        return root_count > 1 ? ALEA_RAY_COVERAGE_OVERLAP
                              : ALEA_RAY_COVERAGE_UNRESOLVED;
    return owners[deepest].resolution_flags & ALEA_RESOLVE_UNDEFINED_FILL
        ? ALEA_RAY_COVERAGE_UNDEFINED_FILL : ALEA_RAY_COVERAGE_UNIQUE;
}

static void coverage_make_legacy_finding(const alea_ray_coverage_owner_t* owners,
                                         size_t owner_count,
                                         alea_ray_coverage_kind_t kind,
                                         double t_enter, double t_exit,
                                         alea_ray_interval_finding_t* out) {
    *out = (alea_ray_interval_finding_t){
        .t_enter = t_enter, .t_exit = t_exit,
        .kind = ALEA_INTERVAL_OK, .cell_id = -1,
        .overlap_cell_id = -1, .depth = -1
    };
    if (kind == ALEA_RAY_COVERAGE_GAP) {
        out->kind = ALEA_INTERVAL_GAP;
        return;
    }
    if (kind == ALEA_RAY_COVERAGE_UNRESOLVED) {
        out->kind = ALEA_INTERVAL_UNRESOLVED;
        return;
    }
    if (kind == ALEA_RAY_COVERAGE_TRUNCATED) {
        out->kind = ALEA_INTERVAL_TRUNCATED;
        return;
    }
    if (kind == ALEA_RAY_COVERAGE_OVERLAP) {
        out->kind = ALEA_INTERVAL_OVERLAP;
        for (size_t i = 0; i < owner_count; i++) {
            for (size_t j = i + 1; j < owner_count; j++) {
                if (owners[i].parent_occurrence_key ==
                    owners[j].parent_occurrence_key) {
                    out->cell_id = owners[i].cell_id;
                    out->overlap_cell_id = owners[j].cell_id;
                    out->depth = owners[i].depth;
                    return;
                }
            }
            if (owners[i].parent_occurrence_key == 0) {
                out->cell_id = owners[i].cell_id;
                for (size_t j = i + 1; j < owner_count; j++) {
                    if (owners[j].parent_occurrence_key == 0) {
                        out->overlap_cell_id = owners[j].cell_id;
                        out->depth = owners[i].depth;
                        return;
                    }
                }
            }
        }
        return;
    }
    size_t deepest = 0;
    for (size_t i = 1; i < owner_count; i++)
        if (owners[i].depth > owners[deepest].depth)
            deepest = i;
    out->cell_id = owners[deepest].cell_id;
    out->depth = owners[deepest].depth;
    if (kind == ALEA_RAY_COVERAGE_UNDEFINED_FILL)
        out->kind = ALEA_INTERVAL_UNDEFINED_FILL;
}

int alea_ray_coverage_sweep_domain_reuse_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    const alea_ray_coverage_domain_t* domain,
    alea_raycast_result_t* breakpoint_scratch,
    alea_ray_coverage_interval_callback_t callback, void* context) {
    if (!sys || !ray || !breakpoint_scratch || !callback || t_max <= 0.0)
        return -1;
    if (domain && domain->has_domain && domain->t_max < domain->t_min)
        return -1;
    if (alea_raycast_global_breakpoints_reuse_nocache(
            sys, ray, t_max, breakpoint_scratch) != 0)
        return -1;

    const int has_domain = domain && domain->has_domain;
    const double domain_t_min = has_domain
        ? (domain->t_min > 0.0 ? domain->t_min : 0.0) : 0.0;
    const double domain_t_max = has_domain
        ? (domain->t_max < t_max ? domain->t_max : t_max) : t_max;
    const int report_allowed_exterior = domain &&
        domain->report_allowed_exterior;

    uint64_t previous_keys[ALEA_RAY_COVERAGE_MAX_OWNERS];
    uint64_t previous_parent_keys[ALEA_RAY_COVERAGE_MAX_OWNERS];
    int previous_count = -1;
    alea_ray_coverage_owner_t current_owners[ALEA_RAY_COVERAGE_MAX_OWNERS];
    alea_ray_coverage_owner_t resolved_owners[ALEA_RAY_COVERAGE_MAX_OWNERS];
    alea_ray_coverage_interval_t current = {0};
    int have_current = 0;
    int total = 0;
    double t_previous = 0.0;

    size_t hit_index = 0;
    while (t_previous < t_max) {
        while (hit_index < breakpoint_scratch->hits.count &&
               breakpoint_scratch->hits.data[hit_index].t <=
                   t_previous + ALEA_RAY_COVERAGE_MIN_INTERVAL)
            hit_index++;
        double t_current = t_max;
        if (hit_index < breakpoint_scratch->hits.count &&
            breakpoint_scratch->hits.data[hit_index].t < t_current)
            t_current = breakpoint_scratch->hits.data[hit_index].t;
        if (has_domain && domain_t_min > t_previous +
                              ALEA_RAY_COVERAGE_MIN_INTERVAL &&
            domain_t_min < t_current)
            t_current = domain_t_min;
        if (has_domain && domain_t_max > t_previous +
                              ALEA_RAY_COVERAGE_MIN_INTERVAL &&
            domain_t_max < t_current)
            t_current = domain_t_max;
        if (t_current - t_previous <= ALEA_RAY_COVERAGE_MIN_INTERVAL) {
            t_previous = t_current;
            continue;
        }
        const double sample_t = t_previous + 0.381966011250105 *
            (t_current - t_previous);
        double x, y, z;
        alea_ray_point_at(ray, sample_t, &x, &y, &z);
        alea_cell_hit_t hits[ALEA_RAY_COVERAGE_MAX_OWNERS];
        uint64_t occurrence_keys[ALEA_RAY_COVERAGE_MAX_OWNERS];
        uint64_t parent_occurrence_keys[ALEA_RAY_COVERAGE_MAX_OWNERS];
        const int hit_count =
            alea_find_all_cells_at_point_coverage_chain_recursive(
                sys, x, y, z, hits, occurrence_keys, parent_occurrence_keys,
                ALEA_RAY_COVERAGE_MAX_OWNERS);
        if (hit_count < 0) return -1;
        const int owners_truncated =
            hit_count >= ALEA_RAY_COVERAGE_MAX_OWNERS;
        const int retained_count = owners_truncated
            ? ALEA_RAY_COVERAGE_OWNER_BUDGET : hit_count;
        for (int owner = 0; owner < retained_count; owner++) {
            resolved_owners[owner] = (alea_ray_coverage_owner_t){
                .cell_id = hits[owner].cell_id,
                .cell_index = hits[owner].cell_index,
                .material_id = hits[owner].material_id,
                .universe_id = hits[owner].universe_id,
                .fill_universe = hits[owner].fill_universe,
                .depth = hits[owner].depth,
                .occurrence_key = occurrence_keys[owner],
                .parent_occurrence_key = parent_occurrence_keys[owner],
                .resolution_flags = hits[owner].resolution_flags
            };
        }
        alea_ray_coverage_kind_t interval_kind = owners_truncated
            ? ALEA_RAY_COVERAGE_TRUNCATED
            : coverage_kind(resolved_owners, (size_t)retained_count);
        const int in_domain = !has_domain ||
            (sample_t >= domain_t_min && sample_t <= domain_t_max);
        if (hit_count == 0 && !in_domain) {
            if (!report_allowed_exterior) {
                if (have_current) {
                    total++;
                    if (callback(context, &current) != 0) return total;
                    have_current = 0;
                }
                t_previous = t_current;
                continue;
            }
            interval_kind = ALEA_RAY_COVERAGE_ALLOWED_EXTERIOR;
        }
        /* A saturated owner buffer proves that the complete set is not
         * available, but it is still useful diagnostic output.  Keep the
         * retained prefix and publish an explicit truncated interval instead
         * of failing the whole sweep. */
        if (!owners_truncated && have_current && coverage_owner_sets_equal(
                occurrence_keys, parent_occurrence_keys, hit_count,
                previous_keys, previous_parent_keys, previous_count) &&
            current.kind == interval_kind) {
            current.t_exit = t_current;
        } else {
            if (have_current) {
                total++;
                if (callback(context, &current) != 0) return total;
            }
            memcpy(current_owners, resolved_owners,
                   (size_t)retained_count * sizeof(*current_owners));
            current = (alea_ray_coverage_interval_t){
                .t_enter = t_previous, .t_exit = t_current,
                .kind = interval_kind,
                .owners = current_owners,
                .owner_count = (size_t)retained_count,
                .owner_count_lower_bound = owners_truncated
                    ? ALEA_RAY_COVERAGE_MAX_OWNERS : (size_t)retained_count
            };
            have_current = 1;
            previous_count = retained_count;
            memcpy(previous_keys, occurrence_keys,
                   (size_t)retained_count * sizeof(*occurrence_keys));
            memcpy(previous_parent_keys, parent_occurrence_keys,
                   (size_t)retained_count * sizeof(*parent_occurrence_keys));
        }
        t_previous = t_current;
        if (t_previous >= t_max) break;
    }
    if (have_current) {
        total++;
        if (callback(context, &current) != 0) return total;
    }
    return total;
}

int alea_ray_coverage_sweep_reuse_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    alea_raycast_result_t* breakpoint_scratch,
    alea_ray_coverage_interval_callback_t callback, void* context) {
    return alea_ray_coverage_sweep_domain_reuse_nocache(
        sys, ray, t_max, NULL, breakpoint_scratch, callback, context);
}

typedef struct {
    alea_ray_coverage_row_interval_callback_t callback;
    void* context;
    size_t row_index;
    int callback_failed;
} coverage_row_callback_context_t;

static int coverage_row_interval_callback(
    void* context, const alea_ray_coverage_interval_t* interval) {
    coverage_row_callback_context_t* row_context = context;
    if (row_context->callback(row_context->context, row_context->row_index,
                              interval) != 0) {
        row_context->callback_failed = 1;
        return 1;
    }
    return 0;
}

int alea_ray_coverage_rows_serial_reuse_nocache(
    alea_system_t* sys, const alea_ray_coverage_row_t* rows, size_t row_count,
    alea_raycast_result_t* breakpoint_scratch,
    alea_ray_coverage_row_interval_callback_t callback, void* context) {
    if (!sys || (!rows && row_count != 0) || !breakpoint_scratch || !callback)
        return -1;
    for (size_t row_index = 0; row_index < row_count; row_index++) {
        const alea_ray_coverage_row_t* row = &rows[row_index];
        if (row->t_max <= 0.0) return -1;
        coverage_row_callback_context_t row_context = {
            .callback = callback,
            .context = context,
            .row_index = row_index,
            .callback_failed = 0
        };
        const int rc = row->use_domain
            ? alea_ray_coverage_sweep_domain_reuse_nocache(
                  sys, &row->ray, row->t_max, &row->domain,
                  breakpoint_scratch, coverage_row_interval_callback,
                  &row_context)
            : alea_ray_coverage_sweep_reuse_nocache(
                  sys, &row->ray, row->t_max, breakpoint_scratch,
                  coverage_row_interval_callback, &row_context);
        if (rc < 0 || row_context.callback_failed) return -1;
    }
    return 0;
}

void alea_ray_coverage_slice_limits_init(
    alea_ray_coverage_slice_limits_t* limits) {
    if (limits) memset(limits, 0, sizeof(*limits));
}

void alea_ray_coverage_slice_result_init(
    alea_ray_coverage_slice_result_t* result) {
    if (result) memset(result, 0, sizeof(*result));
}

void alea_ray_coverage_slice_result_free(
    alea_ray_coverage_slice_result_t* result) {
    if (!result) return;
    free(result->row_offsets);
    free(result->row_direction_tags);
    free(result->row_transverse_coordinates);
    free(result->t_enter);
    free(result->t_exit);
    free(result->kinds);
    free(result->owner_offsets);
    free(result->owner_count_lower_bounds);
    free(result->owner_cell_ids);
    free(result->owner_cell_indices);
    free(result->owner_material_ids);
    free(result->owner_universe_ids);
    free(result->owner_fill_universes);
    free(result->owner_depths);
    free(result->owner_occurrence_keys);
    free(result->owner_parent_occurrence_keys);
    free(result->owner_resolution_flags);
    memset(result, 0, sizeof(*result));
}

static int coverage_slice_add_bytes(size_t* total, size_t count,
                                    size_t element_size) {
    if (count != 0 && element_size > (SIZE_MAX - *total) / count)
        return -1;
    *total += count * element_size;
    return 0;
}

static int coverage_slice_published_bytes(size_t rows, size_t intervals,
                                          size_t owners, size_t* bytes) {
    size_t total = 0;
    if (rows == SIZE_MAX || intervals == SIZE_MAX ||
        coverage_slice_add_bytes(&total, rows + 1, sizeof(size_t)) ||
        coverage_slice_add_bytes(&total, rows, sizeof(uint8_t)) ||
        coverage_slice_add_bytes(&total, rows, sizeof(double)) ||
        coverage_slice_add_bytes(&total, intervals, sizeof(double)) ||
        coverage_slice_add_bytes(&total, intervals, sizeof(double)) ||
        coverage_slice_add_bytes(&total, intervals, sizeof(uint8_t)) ||
        coverage_slice_add_bytes(&total, intervals + 1, sizeof(size_t)) ||
        coverage_slice_add_bytes(&total, intervals, sizeof(size_t)) ||
        coverage_slice_add_bytes(&total, owners, sizeof(int)) ||
        coverage_slice_add_bytes(&total, owners, sizeof(int)) ||
        coverage_slice_add_bytes(&total, owners, sizeof(int)) ||
        coverage_slice_add_bytes(&total, owners, sizeof(int)) ||
        coverage_slice_add_bytes(&total, owners, sizeof(int)) ||
        coverage_slice_add_bytes(&total, owners, sizeof(int)) ||
        coverage_slice_add_bytes(&total, owners, sizeof(uint64_t)) ||
        coverage_slice_add_bytes(&total, owners, sizeof(uint64_t)) ||
        coverage_slice_add_bytes(&total, owners, sizeof(uint8_t)))
        return -1;
    *bytes = total;
    return 0;
}

#define COVERAGE_SLICE_REALLOC(result, member, count) do { \
    void* coverage_slice_reallocated = realloc((result)->member, \
        (count) * sizeof(*(result)->member)); \
    if (!coverage_slice_reallocated) return -1; \
    (result)->member = coverage_slice_reallocated; \
} while (0)

static int coverage_slice_reserve_rows(alea_ray_coverage_slice_result_t* result,
                                       size_t row_count) {
    COVERAGE_SLICE_REALLOC(result, row_offsets, row_count + 1);
    if (row_count == 0) return 0;
    COVERAGE_SLICE_REALLOC(result, row_direction_tags, row_count);
    COVERAGE_SLICE_REALLOC(result, row_transverse_coordinates, row_count);
    return 0;
}

static int coverage_slice_reserve_intervals(
    alea_ray_coverage_slice_result_t* result, size_t interval_count) {
    if (interval_count != 0) {
        COVERAGE_SLICE_REALLOC(result, t_enter, interval_count);
        COVERAGE_SLICE_REALLOC(result, t_exit, interval_count);
        COVERAGE_SLICE_REALLOC(result, kinds, interval_count);
        COVERAGE_SLICE_REALLOC(result, owner_count_lower_bounds, interval_count);
    }
    COVERAGE_SLICE_REALLOC(result, owner_offsets, interval_count + 1);
    return 0;
}

static int coverage_slice_reserve_owners(
    alea_ray_coverage_slice_result_t* result, size_t owner_count) {
    if (owner_count == 0) return 0;
    COVERAGE_SLICE_REALLOC(result, owner_cell_ids, owner_count);
    COVERAGE_SLICE_REALLOC(result, owner_cell_indices, owner_count);
    COVERAGE_SLICE_REALLOC(result, owner_material_ids, owner_count);
    COVERAGE_SLICE_REALLOC(result, owner_universe_ids, owner_count);
    COVERAGE_SLICE_REALLOC(result, owner_fill_universes, owner_count);
    COVERAGE_SLICE_REALLOC(result, owner_depths, owner_count);
    COVERAGE_SLICE_REALLOC(result, owner_occurrence_keys, owner_count);
    COVERAGE_SLICE_REALLOC(result, owner_parent_occurrence_keys, owner_count);
    COVERAGE_SLICE_REALLOC(result, owner_resolution_flags, owner_count);
    return 0;
}

typedef struct {
    alea_ray_coverage_slice_result_t* result;
    const alea_ray_coverage_slice_limits_t* limits;
    size_t next_row_offset;
    int failed;
    /* Set only by the parallel executor.  Local arena counts remain useful
     * for compaction, but resource limits describe one published operation,
     * not an independently budgeted worker shard. */
    struct coverage_executor_counters* executor_counters;
} coverage_slice_builder_t;

typedef struct coverage_executor_counters {
    atomic_size_t intervals;
    atomic_size_t owners;
    atomic_size_t bytes;
} coverage_executor_counters_t;

/* Atomically reserve a portion of an operation-wide budget.  Reservations are
 * deliberately not rolled back on a later allocation failure: that operation
 * is discarded transactionally, and retaining the reservation prevents other
 * workers from growing staging while failure is propagating. */
static int coverage_executor_reserve(atomic_size_t* counter, size_t amount,
                                     size_t limit) {
    size_t current = atomic_load(counter);
    for (;;) {
        if (amount > SIZE_MAX - current ||
            (limit != 0 && amount > limit - current))
            return -1;
        if (atomic_compare_exchange_weak(counter, &current, current + amount))
            return 0;
    }
}

static int coverage_slice_append(void* context, size_t row_index,
                                 const alea_ray_coverage_interval_t* interval) {
    coverage_slice_builder_t* builder = context;
    alea_ray_coverage_slice_result_t* result = builder->result;
    const size_t interval_count = result->interval_count;
    const size_t owner_count = result->owner_count;
    size_t bytes;
    if (interval_count == SIZE_MAX ||
        interval->owner_count > SIZE_MAX - owner_count ||
        coverage_slice_published_bytes(result->row_count, interval_count + 1,
                                       owner_count + interval->owner_count,
                                       &bytes)) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "coverage slice output limit exceeded");
        builder->failed = 1;
        return -1;
    }
    if (builder->executor_counters) {
        size_t empty_bytes, appended_bytes;
        if (coverage_slice_published_bytes(0, 0, 0, &empty_bytes) ||
            coverage_slice_published_bytes(0, 1, interval->owner_count,
                                           &appended_bytes) ||
            appended_bytes < empty_bytes ||
            coverage_executor_reserve(&builder->executor_counters->intervals,
                                      1, builder->limits->max_intervals) ||
            coverage_executor_reserve(&builder->executor_counters->owners,
                                      interval->owner_count,
                                      builder->limits->max_owners) ||
            coverage_executor_reserve(&builder->executor_counters->bytes,
                                      appended_bytes - empty_bytes,
                                      builder->limits->max_bytes)) {
            alea_set_error_detail(ALEA_ERR_OVERFLOW,
                                  "coverage slice output limit exceeded");
            builder->failed = 1;
            return -1;
        }
    } else if ((builder->limits->max_intervals != 0 &&
                interval_count >= builder->limits->max_intervals) ||
               (builder->limits->max_owners != 0 &&
                interval->owner_count > builder->limits->max_owners - owner_count) ||
               (builder->limits->max_bytes != 0 &&
                bytes > builder->limits->max_bytes)) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "coverage slice output limit exceeded");
        builder->failed = 1;
        return -1;
    }
    while (builder->next_row_offset <= row_index)
        result->row_offsets[builder->next_row_offset++] = interval_count;

    if (coverage_slice_reserve_intervals(result, interval_count + 1) != 0) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate coverage slice intervals");
        builder->failed = 1;
        return -1;
    }
    if (coverage_slice_reserve_owners(
            result, owner_count + interval->owner_count) != 0) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate coverage slice owners");
        builder->failed = 1;
        return -1;
    }
    result->t_enter[interval_count] = interval->t_enter;
    result->t_exit[interval_count] = interval->t_exit;
    result->kinds[interval_count] = (uint8_t)interval->kind;
    result->owner_offsets[interval_count] = owner_count;
    result->owner_offsets[interval_count + 1] = owner_count + interval->owner_count;
    result->owner_count_lower_bounds[interval_count] =
        interval->owner_count_lower_bound;
    for (size_t owner = 0; owner < interval->owner_count; owner++) {
        const size_t index = owner_count + owner;
        const alea_ray_coverage_owner_t* source = &interval->owners[owner];
        result->owner_cell_ids[index] = source->cell_id;
        result->owner_cell_indices[index] = source->cell_index;
        result->owner_material_ids[index] = source->material_id;
        result->owner_universe_ids[index] = source->universe_id;
        result->owner_fill_universes[index] = source->fill_universe;
        result->owner_depths[index] = source->depth;
        result->owner_occurrence_keys[index] = source->occurrence_key;
        result->owner_parent_occurrence_keys[index] = source->parent_occurrence_key;
        result->owner_resolution_flags[index] = source->resolution_flags;
    }
    result->interval_count++;
    result->owner_count += interval->owner_count;
    return 0;
}

int alea_ray_coverage_slice_build_serial_nocache(
    alea_system_t* sys, const alea_ray_coverage_row_t* rows, size_t row_count,
    const alea_ray_coverage_slice_limits_t* limits,
    alea_raycast_result_t* breakpoint_scratch,
    alea_ray_coverage_slice_result_t* result) {
    if (!sys || (!rows && row_count != 0) || !breakpoint_scratch || !result)
        return -1;
    const alea_ray_coverage_slice_limits_t unlimited = {0};
    if (!limits) limits = &unlimited;
    size_t row_bytes;
    if ((limits->max_rows != 0 && row_count > limits->max_rows) ||
        coverage_slice_published_bytes(row_count, 0, 0, &row_bytes) != 0 ||
        (limits->max_bytes != 0 && row_bytes > limits->max_bytes)) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "coverage slice row limit exceeded");
        return -1;
    }

    alea_ray_coverage_slice_result_t candidate;
    alea_ray_coverage_slice_result_init(&candidate);
    candidate.row_count = row_count;
    if (coverage_slice_reserve_rows(&candidate, row_count) != 0) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate coverage slice rows");
        alea_ray_coverage_slice_result_free(&candidate);
        return -1;
    }
    if (coverage_slice_reserve_intervals(&candidate, 0) != 0) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate coverage slice offsets");
        alea_ray_coverage_slice_result_free(&candidate);
        return -1;
    }
    for (size_t row = 0; row < row_count; row++) {
        candidate.row_direction_tags[row] = rows[row].direction_tag;
        candidate.row_transverse_coordinates[row] = rows[row].transverse_coordinate;
    }
    coverage_slice_builder_t builder = {
        .result = &candidate, .limits = limits, .next_row_offset = 0
    };
    if (alea_ray_coverage_rows_serial_reuse_nocache(
            sys, rows, row_count, breakpoint_scratch, coverage_slice_append,
            &builder) != 0 || builder.failed) {
        alea_ray_coverage_slice_result_free(&candidate);
        return -1;
    }
    while (builder.next_row_offset <= row_count)
        candidate.row_offsets[builder.next_row_offset++] = candidate.interval_count;
    alea_ray_coverage_slice_result_free(result);
    *result = candidate;
    return 0;
}

static int coverage_slice_row_interval_range(
    const alea_ray_coverage_slice_result_t* result, size_t row,
    size_t* begin, size_t* end) {
    if (!result || row >= result->row_count || !result->row_offsets ||
        !result->owner_offsets ||
        (result->interval_count != 0 &&
         (!result->kinds || !result->owner_count_lower_bounds)) ||
        (result->owner_count != 0 &&
         (!result->owner_occurrence_keys ||
          !result->owner_parent_occurrence_keys ||
          !result->owner_resolution_flags)))
        return -1;
    *begin = result->row_offsets[row];
    *end = result->row_offsets[row + 1];
    if (*begin > *end || *end > result->interval_count ||
        result->owner_offsets[0] != 0 ||
        result->owner_offsets[result->interval_count] != result->owner_count)
        return -1;
    return 0;
}

int alea_ray_coverage_slice_rows_same_signature(
    const alea_ray_coverage_slice_result_t* result,
    size_t first_row, size_t second_row) {
    size_t first_begin, first_end, second_begin, second_end;
    if (coverage_slice_row_interval_range(result, first_row,
                                          &first_begin, &first_end) != 0 ||
        coverage_slice_row_interval_range(result, second_row,
                                          &second_begin, &second_end) != 0)
        return -1;
    if (first_end - first_begin != second_end - second_begin) return 0;
    for (size_t offset = 0; offset < first_end - first_begin; offset++) {
        const size_t first = first_begin + offset;
        const size_t second = second_begin + offset;
        const size_t first_owners = result->owner_offsets[first + 1] -
            result->owner_offsets[first];
        const size_t second_owners = result->owner_offsets[second + 1] -
            result->owner_offsets[second];
        if (result->owner_offsets[first] > result->owner_offsets[first + 1] ||
            result->owner_offsets[second] > result->owner_offsets[second + 1] ||
            result->owner_offsets[first + 1] > result->owner_count ||
            result->owner_offsets[second + 1] > result->owner_count)
            return -1;
        if (result->kinds[first] != result->kinds[second] ||
            result->owner_count_lower_bounds[first] !=
                result->owner_count_lower_bounds[second] ||
            first_owners != second_owners)
            return 0;
        for (size_t owner = 0; owner < first_owners; owner++) {
            const size_t first_owner = result->owner_offsets[first] + owner;
            const size_t second_owner = result->owner_offsets[second] + owner;
            if (result->owner_occurrence_keys[first_owner] !=
                    result->owner_occurrence_keys[second_owner] ||
                result->owner_parent_occurrence_keys[first_owner] !=
                    result->owner_parent_occurrence_keys[second_owner] ||
                result->owner_resolution_flags[first_owner] !=
                    result->owner_resolution_flags[second_owner])
                return 0;
        }
    }
    return 1;
}

void alea_ray_coverage_refinement_policy_init(
    alea_ray_coverage_refinement_policy_t* policy) {
    if (!policy) return;
    memset(policy, 0, sizeof(*policy));
    policy->signals = ALEA_RAY_COVERAGE_REFINE_SIGNATURE;
}

/* A coverage kind that a probe should look at more closely.  Allowed exterior
 * is a configured domain answer, not a defect. */
static int coverage_kind_is_finding(uint8_t kind) {
    return kind != ALEA_RAY_COVERAGE_UNIQUE &&
           kind != ALEA_RAY_COVERAGE_ALLOWED_EXTERIOR;
}

static int coverage_row_has_finding(
    const alea_ray_coverage_slice_result_t* result, size_t begin, size_t end) {
    for (size_t i = begin; i < end; i++)
        if (coverage_kind_is_finding(result->kinds[i])) return 1;
    return 0;
}

/* Paired endpoint displacement between two rows of equal interval count.  The
 * signature deliberately ignores endpoints, so this is the only signal that
 * sees a boundary sliding while owner identity stays put. */
static int coverage_rows_endpoints_displaced(
    const alea_ray_coverage_slice_result_t* result, size_t first_begin,
    size_t second_begin, size_t count, double tolerance) {
    if (!result->t_enter || !result->t_exit) return -1;
    for (size_t offset = 0; offset < count; offset++) {
        const size_t first = first_begin + offset;
        const size_t second = second_begin + offset;
        if (fabs(result->t_enter[first] - result->t_enter[second]) > tolerance ||
            fabs(result->t_exit[first] - result->t_exit[second]) > tolerance)
            return 1;
    }
    return 0;
}

int alea_ray_coverage_slice_mark_refinement_boundaries_policy(
    const alea_ray_coverage_slice_result_t* result,
    const alea_ray_coverage_refinement_policy_t* policy,
    uint8_t* out_refine_between, size_t* out_spacing_limited) {
    alea_ray_coverage_refinement_policy_t defaults;
    if (!policy) {
        alea_ray_coverage_refinement_policy_init(&defaults);
        policy = &defaults;
    }
    if (!result || (result->row_count > 1 && !out_refine_between) ||
        (result->row_count != 0 && (!result->row_direction_tags ||
                                    !result->row_transverse_coordinates)))
        return -1;
    if (policy->signals & ~(uint32_t)(ALEA_RAY_COVERAGE_REFINE_SIGNATURE |
                                      ALEA_RAY_COVERAGE_REFINE_DISPLACEMENT |
                                      ALEA_RAY_COVERAGE_REFINE_DENSITY |
                                      ALEA_RAY_COVERAGE_REFINE_FINDING))
        return -1;
    if (policy->min_transverse_spacing < 0.0 ||
        policy->endpoint_displacement < 0.0 ||
        ((policy->signals & ALEA_RAY_COVERAGE_REFINE_DISPLACEMENT) &&
         !(policy->endpoint_displacement > 0.0)) ||
        ((policy->signals & ALEA_RAY_COVERAGE_REFINE_DENSITY) &&
         policy->crossing_density == 0))
        return -1;
    if (out_spacing_limited) *out_spacing_limited = 0;

    int marked = 0;
    for (size_t row = 0; row + 1 < result->row_count; row++) {
        out_refine_between[row] = 0;
        if (result->row_direction_tags[row] != result->row_direction_tags[row + 1])
            continue;
        const double first_coordinate = result->row_transverse_coordinates[row];
        const double second_coordinate =
            result->row_transverse_coordinates[row + 1];
        if (first_coordinate >= second_coordinate) return -1;
        size_t first_begin, first_end, second_begin, second_end;
        if (coverage_slice_row_interval_range(result, row, &first_begin,
                                              &first_end) != 0 ||
            coverage_slice_row_interval_range(result, row + 1, &second_begin,
                                              &second_end) != 0)
            return -1;
        const size_t first_count = first_end - first_begin;
        const size_t second_count = second_end - second_begin;

        int select = 0;
        if (policy->signals & ALEA_RAY_COVERAGE_REFINE_SIGNATURE) {
            const int same = alea_ray_coverage_slice_rows_same_signature(
                result, row, row + 1);
            if (same < 0) return -1;
            if (!same) select = 1;
        }
        if (!select && (policy->signals & ALEA_RAY_COVERAGE_REFINE_DENSITY) &&
            (first_count >= policy->crossing_density ||
             second_count >= policy->crossing_density))
            select = 1;
        if (!select && (policy->signals & ALEA_RAY_COVERAGE_REFINE_FINDING) &&
            (coverage_row_has_finding(result, first_begin, first_end) ||
             coverage_row_has_finding(result, second_begin, second_end)))
            select = 1;
        if (!select && (policy->signals & ALEA_RAY_COVERAGE_REFINE_DISPLACEMENT) &&
            first_count == second_count) {
            const int displaced = coverage_rows_endpoints_displaced(
                result, first_begin, second_begin, first_count,
                policy->endpoint_displacement);
            if (displaced < 0) return -1;
            if (displaced) select = 1;
        }
        if (!select) continue;

        /* Splitting halves the gap; refuse to generate rows the caller has
         * declared too close to distinguish. */
        if (policy->min_transverse_spacing > 0.0 &&
            0.5 * (second_coordinate - first_coordinate) <
                policy->min_transverse_spacing) {
            if (out_spacing_limited) (*out_spacing_limited)++;
            continue;
        }
        out_refine_between[row] = 1;
        marked++;
    }
    return marked;
}

int alea_ray_coverage_slice_mark_refinement_boundaries(
    const alea_ray_coverage_slice_result_t* result,
    uint8_t* out_refine_between) {
    return alea_ray_coverage_slice_mark_refinement_boundaries_policy(
        result, NULL, out_refine_between, NULL);
}

static int coverage_rows_midpoint(const alea_ray_coverage_row_t* first,
                                  const alea_ray_coverage_row_t* second,
                                  alea_ray_coverage_row_t* midpoint) {
    if (first->direction_tag != second->direction_tag ||
        first->transverse_coordinate >= second->transverse_coordinate ||
        first->ray.dx != second->ray.dx || first->ray.dy != second->ray.dy ||
        first->ray.dz != second->ray.dz ||
        first->use_domain != second->use_domain ||
        first->domain.has_domain != second->domain.has_domain ||
        first->domain.report_allowed_exterior !=
            second->domain.report_allowed_exterior)
        return -1;
    *midpoint = *first;
    midpoint->transverse_coordinate = 0.5 *
        (first->transverse_coordinate + second->transverse_coordinate);
    midpoint->t_max = 0.5 * (first->t_max + second->t_max);
    midpoint->domain.t_min = 0.5 *
        (first->domain.t_min + second->domain.t_min);
    midpoint->domain.t_max = 0.5 *
        (first->domain.t_max + second->domain.t_max);
    return alea_ray_init(&midpoint->ray,
                         0.5 * (first->ray.ox + second->ray.ox),
                         0.5 * (first->ray.oy + second->ray.oy),
                         0.5 * (first->ray.oz + second->ray.oz),
                         first->ray.dx, first->ray.dy, first->ray.dz);
}

int alea_ray_coverage_rows_refine_midpoints(
    const alea_ray_coverage_row_t* rows, size_t row_count,
    const uint8_t* refine_between, size_t max_rows,
    alea_ray_coverage_row_t* out_rows, size_t out_capacity,
    size_t* out_row_count) {
    if ((!rows && row_count != 0) || !out_row_count ||
        (row_count > 1 && !refine_between) ||
        (row_count != 0 && out_rows == rows))
        return -1;
    size_t marked = 0;
    for (size_t row = 0; row + 1 < row_count; row++) {
        if (!refine_between[row]) continue;
        alea_ray_coverage_row_t midpoint;
        if (coverage_rows_midpoint(&rows[row], &rows[row + 1],
                                   &midpoint) != 0)
            return -1;
        if (marked == SIZE_MAX - row_count) return -1;
        marked++;
    }
    const size_t refined_count = row_count + marked;
    if ((max_rows != 0 && refined_count > max_rows) ||
        refined_count > out_capacity ||
        (refined_count != 0 && !out_rows)) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "coverage refinement row limit exceeded");
        return -1;
    }
    size_t output = 0;
    for (size_t row = 0; row < row_count; row++) {
        out_rows[output++] = rows[row];
        if (row + 1 < row_count && refine_between[row]) {
            if (coverage_rows_midpoint(&rows[row], &rows[row + 1],
                                       &out_rows[output]) != 0)
                return -1;
            output++;
        }
    }
    *out_row_count = refined_count;
    return 0;
}

int alea_ray_coverage_slice_build_adaptive_policy_serial_nocache(
    alea_system_t* sys, const alea_ray_coverage_row_t* initial_rows,
    size_t initial_row_count, size_t max_refinement_depth,
    const alea_ray_coverage_refinement_policy_t* policy,
    const alea_ray_coverage_slice_limits_t* limits,
    alea_raycast_result_t* breakpoint_scratch,
    alea_ray_coverage_slice_result_t* result) {
    if (!sys || (!initial_rows && initial_row_count != 0) ||
        !breakpoint_scratch || !result)
        return -1;
    const alea_ray_coverage_slice_limits_t unlimited = {0};
    if (!limits) limits = &unlimited;
    const alea_ray_coverage_row_t* current_rows = initial_rows;
    size_t current_row_count = initial_row_count;
    alea_ray_coverage_row_t* owned_rows = NULL;
    alea_ray_coverage_slice_result_t current;
    alea_ray_coverage_slice_result_init(&current);

    for (size_t depth = 0;; depth++) {
        if (alea_ray_coverage_slice_build_serial_nocache(
                sys, current_rows, current_row_count, limits,
                breakpoint_scratch, &current) != 0)
            goto fail;
        uint8_t* refine_between = NULL;
        if (current_row_count > 1) {
            refine_between = calloc(current_row_count - 1,
                                    sizeof(*refine_between));
            if (!refine_between) {
                alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                                      "failed to allocate coverage refinement marks");
                goto fail;
            }
        }
        size_t spacing_limited = 0;
        const int marked = alea_ray_coverage_slice_mark_refinement_boundaries_policy(
            &current, policy, refine_between, &spacing_limited);
        if (marked < 0) {
            free(refine_between);
            goto fail;
        }
        if (marked == 0 || depth >= max_refinement_depth) {
            if (marked != 0)
                current.refinement_status =
                    ALEA_RAY_COVERAGE_REFINEMENT_MAX_DEPTH;
            else if (spacing_limited != 0)
                current.refinement_status =
                    ALEA_RAY_COVERAGE_REFINEMENT_MIN_SPACING;
            else
                current.refinement_status =
                    ALEA_RAY_COVERAGE_REFINEMENT_COMPLETE;
            free(refine_between);
            free(owned_rows);
            alea_ray_coverage_slice_result_free(result);
            *result = current;
            return 0;
        }
        if ((size_t)marked > SIZE_MAX - current_row_count ||
            (limits->max_rows != 0 &&
             current_row_count + (size_t)marked > limits->max_rows)) {
            current.refinement_status = ALEA_RAY_COVERAGE_REFINEMENT_MAX_ROWS;
            free(refine_between);
            free(owned_rows);
            alea_ray_coverage_slice_result_free(result);
            *result = current;
            return 0;
        }
        const size_t next_row_count = current_row_count + (size_t)marked;
        if (next_row_count > SIZE_MAX / sizeof(*owned_rows)) {
            free(refine_between);
            alea_set_error_detail(ALEA_ERR_OVERFLOW,
                                  "coverage refinement row storage overflows");
            goto fail;
        }
        alea_ray_coverage_row_t* next_rows = malloc(
            next_row_count * sizeof(*next_rows));
        if (!next_rows) {
            free(refine_between);
            alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                                  "failed to allocate coverage refinement rows");
            goto fail;
        }
        size_t produced = 0;
        const int refine_rc = alea_ray_coverage_rows_refine_midpoints(
            current_rows, current_row_count, refine_between, limits->max_rows,
            next_rows, next_row_count, &produced);
        free(refine_between);
        if (refine_rc != 0 || produced != next_row_count) {
            free(next_rows);
            goto fail;
        }
        free(owned_rows);
        owned_rows = next_rows;
        current_rows = owned_rows;
        current_row_count = next_row_count;
    }

fail:
    free(owned_rows);
    alea_ray_coverage_slice_result_free(&current);
    return -1;
}

int alea_ray_coverage_slice_build_adaptive_serial_nocache(
    alea_system_t* sys, const alea_ray_coverage_row_t* initial_rows,
    size_t initial_row_count, size_t max_refinement_depth,
    const alea_ray_coverage_slice_limits_t* limits,
    alea_raycast_result_t* breakpoint_scratch,
    alea_ray_coverage_slice_result_t* result) {
    return alea_ray_coverage_slice_build_adaptive_policy_serial_nocache(
        sys, initial_rows, initial_row_count, max_refinement_depth, NULL,
        limits, breakpoint_scratch, result);
}

void alea_ray_coverage_executor_init(alea_ray_coverage_executor_t* executor) {
    if (executor) memset(executor, 0, sizeof(*executor));
}

void alea_ray_coverage_executor_free(alea_ray_coverage_executor_t* executor) {
    if (!executor) return;
    for (size_t worker = 0; worker < executor->worker_count; worker++) {
        alea_raycast_result_free(&executor->workers[worker].breakpoint_scratch);
        alea_ray_coverage_slice_result_free(&executor->workers[worker].arena);
    }
    free(executor->workers);
    memset(executor, 0, sizeof(*executor));
}

int alea_ray_coverage_executor_prepare(alea_ray_coverage_executor_t* executor,
                                        size_t worker_count) {
    if (!executor || worker_count == 0 ||
        worker_count > SIZE_MAX / sizeof(*executor->workers)) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "invalid coverage executor worker count");
        return -1;
    }
    alea_ray_coverage_worker_scratch_t* workers = calloc(
        worker_count, sizeof(*workers));
    if (!workers) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate coverage executor workers");
        return -1;
    }
    for (size_t worker = 0; worker < worker_count; worker++)
        alea_raycast_result_init(&workers[worker].breakpoint_scratch);
    alea_ray_coverage_executor_free(executor);
    executor->workers = workers;
    executor->worker_count = worker_count;
    return 0;
}

alea_ray_coverage_worker_scratch_t*
alea_ray_coverage_executor_worker_for_row(
    alea_ray_coverage_executor_t* executor, size_t row_index) {
    if (!executor || !executor->workers || executor->worker_count == 0)
        return NULL;
    return &executor->workers[row_index % executor->worker_count];
}

/* Keep arena allocations owned by the executor across operations.  Counts and
 * row provenance are overwritten before use; published result ownership never
 * aliases this staging storage. */
static int coverage_executor_reset_arena(
    alea_ray_coverage_worker_scratch_t* worker, size_t row_count) {
    alea_ray_coverage_slice_result_t* arena = &worker->arena;
    arena->row_count = row_count;
    arena->interval_count = 0;
    arena->owner_count = 0;
    arena->refinement_status = ALEA_RAY_COVERAGE_REFINEMENT_COMPLETE;
    if (coverage_slice_reserve_rows(arena, row_count) != 0 ||
        coverage_slice_reserve_intervals(arena, 0) != 0)
        return -1;
    return 0;
}

static int coverage_executor_build_worker(
    alea_system_t* sys, const alea_ray_coverage_row_t* rows, size_t row_count,
    const alea_ray_coverage_slice_limits_t* limits,
    alea_ray_coverage_executor_t* executor, size_t worker_index,
    coverage_executor_counters_t* counters) {
    alea_ray_coverage_worker_scratch_t* worker = &executor->workers[worker_index];
    const size_t workers = executor->worker_count;
    const size_t owned_rows = row_count > worker_index
        ? 1 + (row_count - 1 - worker_index) / workers : 0;
    if (coverage_executor_reset_arena(worker, owned_rows) != 0) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate coverage worker arena");
        return -1;
    }
    coverage_slice_builder_t builder = {
        .result = &worker->arena, .limits = limits, .next_row_offset = 0,
        .executor_counters = counters
    };
    size_t local_row = 0;
    for (size_t row = worker_index; row < row_count; row += workers, local_row++) {
        worker->arena.row_direction_tags[local_row] = rows[row].direction_tag;
        worker->arena.row_transverse_coordinates[local_row] =
            rows[row].transverse_coordinate;
        if (rows[row].t_max <= 0.0) return -1;
        coverage_row_callback_context_t context = {
            .callback = coverage_slice_append, .context = &builder,
            .row_index = local_row, .callback_failed = 0
        };
        const int rc = rows[row].use_domain
            ? alea_ray_coverage_sweep_domain_reuse_nocache(
                  sys, &rows[row].ray, rows[row].t_max, &rows[row].domain,
                  &worker->breakpoint_scratch, coverage_row_interval_callback,
                  &context)
            : alea_ray_coverage_sweep_reuse_nocache(
                  sys, &rows[row].ray, rows[row].t_max,
                  &worker->breakpoint_scratch, coverage_row_interval_callback,
                  &context);
        if (rc < 0 || context.callback_failed || builder.failed) return -1;
    }
    while (builder.next_row_offset <= owned_rows)
        worker->arena.row_offsets[builder.next_row_offset++] =
            worker->arena.interval_count;
    return 0;
}

static int coverage_executor_compact(
    const alea_ray_coverage_row_t* rows, size_t row_count,
    const alea_ray_coverage_slice_limits_t* limits,
    alea_ray_coverage_executor_t* executor,
    alea_ray_coverage_slice_result_t* result) {
    alea_ray_coverage_slice_result_t candidate;
    alea_ray_coverage_slice_result_init(&candidate);
    candidate.row_count = row_count;
    if (coverage_slice_reserve_rows(&candidate, row_count) != 0 ||
        coverage_slice_reserve_intervals(&candidate, 0) != 0) {
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                              "failed to allocate compacted coverage slice");
        goto fail;
    }
    coverage_slice_builder_t builder = {
        .result = &candidate, .limits = limits, .next_row_offset = 0
    };
    for (size_t row = 0; row < row_count; row++) {
        const size_t worker_index = row % executor->worker_count;
        const size_t local_row = row / executor->worker_count;
        const alea_ray_coverage_slice_result_t* arena =
            &executor->workers[worker_index].arena;
        const size_t begin = arena->row_offsets[local_row];
        const size_t end = arena->row_offsets[local_row + 1];
        candidate.row_direction_tags[row] = rows[row].direction_tag;
        candidate.row_transverse_coordinates[row] = rows[row].transverse_coordinate;
        for (size_t interval_index = begin; interval_index < end; interval_index++) {
            const size_t owner_begin = arena->owner_offsets[interval_index];
            const size_t owner_end = arena->owner_offsets[interval_index + 1];
            if (owner_end < owner_begin || owner_end - owner_begin >
                ALEA_RAY_COVERAGE_OWNER_BUDGET) goto fail;
            alea_ray_coverage_owner_t owners[ALEA_RAY_COVERAGE_OWNER_BUDGET];
            for (size_t owner = 0; owner < owner_end - owner_begin; owner++) {
                const size_t source = owner_begin + owner;
                owners[owner] = (alea_ray_coverage_owner_t){
                    .cell_id = arena->owner_cell_ids[source],
                    .cell_index = arena->owner_cell_indices[source],
                    .material_id = arena->owner_material_ids[source],
                    .universe_id = arena->owner_universe_ids[source],
                    .fill_universe = arena->owner_fill_universes[source],
                    .depth = arena->owner_depths[source],
                    .occurrence_key = arena->owner_occurrence_keys[source],
                    .parent_occurrence_key = arena->owner_parent_occurrence_keys[source],
                    .resolution_flags = arena->owner_resolution_flags[source]
                };
            }
            const alea_ray_coverage_interval_t interval = {
                .t_enter = arena->t_enter[interval_index],
                .t_exit = arena->t_exit[interval_index],
                .kind = (alea_ray_coverage_kind_t)arena->kinds[interval_index],
                .owners = owners, .owner_count = owner_end - owner_begin,
                .owner_count_lower_bound =
                    arena->owner_count_lower_bounds[interval_index]
            };
            if (coverage_slice_append(&builder, row, &interval) != 0) goto fail;
        }
    }
    while (builder.next_row_offset <= row_count)
        candidate.row_offsets[builder.next_row_offset++] = candidate.interval_count;
    alea_ray_coverage_slice_result_free(result);
    *result = candidate;
    return 0;
fail:
    alea_ray_coverage_slice_result_free(&candidate);
    return -1;
}

int alea_ray_coverage_slice_build_executor_nocache(
    alea_system_t* sys, const alea_ray_coverage_row_t* rows, size_t row_count,
    const alea_ray_coverage_slice_limits_t* limits,
    alea_ray_coverage_executor_t* executor,
    alea_ray_coverage_slice_result_t* result) {
    if (!sys || (!rows && row_count != 0) || !executor ||
        !executor->workers || executor->worker_count == 0 || !result)
        return -1;
    const alea_ray_coverage_slice_limits_t unlimited = {0};
    if (!limits) limits = &unlimited;
    size_t row_bytes;
    if ((limits->max_rows != 0 && row_count > limits->max_rows) ||
        coverage_slice_published_bytes(row_count, 0, 0, &row_bytes) != 0 ||
        (limits->max_bytes != 0 && row_bytes > limits->max_bytes)) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "coverage slice row limit exceeded");
        return -1;
    }
    int failed = 0;
    coverage_executor_counters_t counters;
    atomic_init(&counters.intervals, 0);
    atomic_init(&counters.owners, 0);
    atomic_init(&counters.bytes, row_bytes);
#ifdef _OPENMP
    if (executor->worker_count <= (size_t)INT_MAX) {
        /* Parallelize worker arenas, not individual rows: a worker index owns
         * its scratch for the entire operation even if the runtime chooses
         * fewer threads than requested. */
        #pragma omp parallel for schedule(static) shared(failed)
        for (size_t worker = 0; worker < executor->worker_count; worker++) {
            int rc = coverage_executor_build_worker(
                sys, rows, row_count, limits, executor, worker, &counters);
            if (rc != 0) {
                #pragma omp atomic write
                failed = 1;
            }
        }
    } else
#endif
    {
        for (size_t worker = 0; worker < executor->worker_count; worker++)
            if (coverage_executor_build_worker(sys, rows, row_count, limits,
                                               executor, worker, &counters) != 0) {
                failed = 1;
                break;
            }
    }
    if (failed || alea_interrupted()) {
        if (alea_interrupted())
            alea_set_error_detail(ALEA_ERR_INTERRUPTED,
                                  "coverage slice execution interrupted");
        return -1;
    }
    return coverage_executor_compact(rows, row_count, limits, executor, result);
}

int alea_ray_coverage_slice_build_adaptive_policy_executor_nocache(
    alea_system_t* sys, const alea_ray_coverage_row_t* initial_rows,
    size_t initial_row_count, size_t max_refinement_depth,
    const alea_ray_coverage_refinement_policy_t* policy,
    const alea_ray_coverage_slice_limits_t* limits,
    alea_ray_coverage_executor_t* executor,
    alea_ray_coverage_slice_result_t* result) {
    if (!sys || (!initial_rows && initial_row_count != 0) || !executor || !result)
        return -1;
    const alea_ray_coverage_slice_limits_t unlimited = {0};
    if (!limits) limits = &unlimited;
    const alea_ray_coverage_row_t* current_rows = initial_rows;
    size_t current_row_count = initial_row_count;
    alea_ray_coverage_row_t* owned_rows = NULL;
    alea_ray_coverage_slice_result_t current;
    alea_ray_coverage_slice_result_init(&current);
    for (size_t depth = 0;; depth++) {
        if (alea_ray_coverage_slice_build_executor_nocache(
                sys, current_rows, current_row_count, limits, executor,
                &current) != 0) goto fail;
        uint8_t* marks = current_row_count > 1
            ? calloc(current_row_count - 1, sizeof(*marks)) : NULL;
        if (current_row_count > 1 && !marks) {
            alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                                  "failed to allocate coverage refinement marks");
            goto fail;
        }
        size_t spacing_limited = 0;
        const int marked = alea_ray_coverage_slice_mark_refinement_boundaries_policy(
            &current, policy, marks, &spacing_limited);
        if (marked < 0) { free(marks); goto fail; }
        if (marked == 0 || depth >= max_refinement_depth) {
            current.refinement_status = marked != 0
                ? ALEA_RAY_COVERAGE_REFINEMENT_MAX_DEPTH
                : spacing_limited != 0 ? ALEA_RAY_COVERAGE_REFINEMENT_MIN_SPACING
                                       : ALEA_RAY_COVERAGE_REFINEMENT_COMPLETE;
            free(marks); free(owned_rows);
            alea_ray_coverage_slice_result_free(result);
            *result = current;
            return 0;
        }
        if ((size_t)marked > SIZE_MAX - current_row_count ||
            (limits->max_rows != 0 &&
             current_row_count + (size_t)marked > limits->max_rows)) {
            current.refinement_status = ALEA_RAY_COVERAGE_REFINEMENT_MAX_ROWS;
            free(marks); free(owned_rows);
            alea_ray_coverage_slice_result_free(result);
            *result = current;
            return 0;
        }
        const size_t next_count = current_row_count + (size_t)marked;
        if (next_count > SIZE_MAX / sizeof(*owned_rows)) { free(marks); goto fail; }
        alea_ray_coverage_row_t* next = malloc(next_count * sizeof(*next));
        if (!next) { free(marks); goto fail; }
        size_t produced = 0;
        const int rc = alea_ray_coverage_rows_refine_midpoints(
            current_rows, current_row_count, marks, limits->max_rows, next,
            next_count, &produced);
        free(marks);
        if (rc != 0 || produced != next_count) { free(next); goto fail; }
        free(owned_rows);
        owned_rows = next;
        current_rows = owned_rows;
        current_row_count = next_count;
    }
fail:
    free(owned_rows);
    alea_ray_coverage_slice_result_free(&current);
    return -1;
}

#undef COVERAGE_SLICE_REALLOC

typedef struct {
    alea_ray_interval_finding_t* out;
    size_t max_out;
    size_t total;
} legacy_coverage_builder_t;

static int coverage_append_legacy(void* context,
                                  const alea_ray_coverage_interval_t* interval) {
    legacy_coverage_builder_t* builder = context;
    if (builder->total < builder->max_out && builder->out) {
        coverage_make_legacy_finding(interval->owners, interval->owner_count,
                                     interval->kind,
                                     interval->t_enter, interval->t_exit,
                                     &builder->out[builder->total]);
    }
    builder->total++;
    return 0;
}

int alea_ray_coverage_classify_reuse_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    alea_raycast_result_t* breakpoint_scratch,
    alea_ray_interval_finding_t* out, size_t max_out) {
    legacy_coverage_builder_t builder = {
        .out = out, .max_out = max_out, .total = 0
    };
    const int result = alea_ray_coverage_sweep_reuse_nocache(
        sys, ray, t_max, breakpoint_scratch, coverage_append_legacy, &builder);
    if (result < 0 || builder.total > (size_t)INT_MAX) return -1;
    return (int)builder.total;
}
