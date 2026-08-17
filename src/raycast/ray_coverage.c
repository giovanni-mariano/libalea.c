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
