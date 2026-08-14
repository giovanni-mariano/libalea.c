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

static int coverage_owner_sets_equal(const uint64_t* a, int a_count,
                                     const uint64_t* b, int b_count) {
    if (a_count != b_count) return 0;
    /* Cell ID is presentation data, not ownership identity.  The recursive
     * coverage resolver derives each key from the concrete fill/lattice path,
     * so repeated or transformed occurrences cannot collapse here. */
    for (int i = 0; i < a_count; i++) {
        if (a[i] != b[i])
            return 0;
    }
    return 1;
}

static void coverage_make_legacy_finding(const alea_ray_coverage_owner_t* owners,
                                         size_t owner_count,
                                         double t_enter, double t_exit,
                                         alea_ray_interval_finding_t* out) {
    *out = (alea_ray_interval_finding_t){
        .t_enter = t_enter, .t_exit = t_exit,
        .kind = ALEA_INTERVAL_OK, .cell_id = -1,
        .overlap_cell_id = -1, .depth = -1
    };
    if (owner_count == 0) {
        out->kind = ALEA_INTERVAL_GAP;
        return;
    }

    for (size_t i = 0; i < owner_count; i++) {
        for (size_t j = i + 1; j < owner_count; j++) {
            if (owners[i].depth == owners[j].depth) {
                out->kind = ALEA_INTERVAL_OVERLAP;
                out->cell_id = owners[i].cell_id;
                out->overlap_cell_id = owners[j].cell_id;
                out->depth = owners[i].depth;
                return;
            }
        }
    }

    size_t deepest = 0;
    while (deepest + 1 < owner_count &&
           owners[deepest + 1].depth == owners[deepest].depth + 1)
        deepest++;
    out->cell_id = owners[deepest].cell_id;
    out->depth = owners[deepest].depth;
    if (owners[deepest].resolution_flags & ALEA_RESOLVE_UNDEFINED_FILL)
        out->kind = ALEA_INTERVAL_UNDEFINED_FILL;
}

static alea_ray_coverage_kind_t coverage_kind(
    const alea_ray_coverage_owner_t* owners, size_t owner_count) {
    if (owner_count == 0) return ALEA_RAY_COVERAGE_GAP;
    for (size_t i = 0; i < owner_count; i++)
        for (size_t j = i + 1; j < owner_count; j++)
            if (owners[i].depth == owners[j].depth)
                return ALEA_RAY_COVERAGE_OVERLAP;
    size_t deepest = 0;
    while (deepest + 1 < owner_count &&
           owners[deepest + 1].depth == owners[deepest].depth + 1)
        deepest++;
    return owners[deepest].resolution_flags & ALEA_RESOLVE_UNDEFINED_FILL
        ? ALEA_RAY_COVERAGE_UNDEFINED_FILL : ALEA_RAY_COVERAGE_UNIQUE;
}

int alea_ray_coverage_sweep_reuse_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    alea_raycast_result_t* breakpoint_scratch,
    alea_ray_coverage_interval_callback_t callback, void* context) {
    if (!sys || !ray || !breakpoint_scratch || !callback || t_max <= 0.0)
        return -1;
    if (alea_raycast_global_breakpoints_reuse_nocache(
            sys, ray, t_max, breakpoint_scratch) != 0)
        return -1;

    uint64_t previous_keys[ALEA_RAY_COVERAGE_MAX_OWNERS];
    int previous_count = -1;
    alea_ray_coverage_owner_t current_owners[ALEA_RAY_COVERAGE_MAX_OWNERS];
    alea_ray_coverage_interval_t current = {0};
    int have_current = 0;
    int total = 0;
    double t_previous = 0.0;

    for (size_t i = 0; i <= breakpoint_scratch->hits.count; i++) {
        double t_current = i < breakpoint_scratch->hits.count
            ? breakpoint_scratch->hits.data[i].t : t_max;
        if (t_current > t_max) t_current = t_max;
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
        const int hit_count = alea_find_all_cells_at_point_coverage_recursive(
            sys, x, y, z, hits, occurrence_keys, ALEA_RAY_COVERAGE_MAX_OWNERS);
        if (hit_count < 0) return -1;
        if (hit_count >= ALEA_RAY_COVERAGE_MAX_OWNERS) {
            alea_set_error_detail(ALEA_ERR_OVERFLOW,
                                  "ray coverage owner budget exceeded");
            return -1;
        }
        if (have_current && coverage_owner_sets_equal(
                occurrence_keys, hit_count, previous_keys, previous_count)) {
            current.t_exit = t_current;
        } else {
            if (have_current) {
                total++;
                if (callback(context, &current) != 0) return total;
            }
            for (int owner = 0; owner < hit_count; owner++) {
                current_owners[owner] = (alea_ray_coverage_owner_t){
                    .cell_id = hits[owner].cell_id,
                    .cell_index = hits[owner].cell_index,
                    .material_id = hits[owner].material_id,
                    .universe_id = hits[owner].universe_id,
                    .fill_universe = hits[owner].fill_universe,
                    .depth = hits[owner].depth,
                    .occurrence_key = occurrence_keys[owner],
                    .resolution_flags = hits[owner].resolution_flags
                };
            }
            current = (alea_ray_coverage_interval_t){
                .t_enter = t_previous, .t_exit = t_current,
                .kind = coverage_kind(current_owners, (size_t)hit_count),
                .owners = current_owners, .owner_count = (size_t)hit_count
            };
            have_current = 1;
            previous_count = hit_count;
            memcpy(previous_keys, occurrence_keys,
                   (size_t)hit_count * sizeof(*occurrence_keys));
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
