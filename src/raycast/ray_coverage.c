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

static int coverage_owner_sets_equal(const alea_cell_hit_t* a, int a_count,
                                     const alea_cell_hit_t* b, int b_count) {
    if (a_count != b_count) return 0;
    /* Cell ID is presentation data, not ownership identity.  Preserve the
     * concrete cell/universe/depth tuple so duplicated user IDs cannot make
     * two different coverage sets merge.  Full occurrence keys will replace
     * this tuple when they are promoted from hierarchy scratch state. */
    for (int i = 0; i < a_count; i++) {
        if (a[i].cell_index != b[i].cell_index ||
            a[i].universe_id != b[i].universe_id ||
            a[i].fill_universe != b[i].fill_universe ||
            a[i].depth != b[i].depth ||
            a[i].resolution_flags != b[i].resolution_flags)
            return 0;
    }
    return 1;
}

static void coverage_make_legacy_finding(const alea_cell_hit_t* hits, int n,
                                         double t_enter, double t_exit,
                                         alea_ray_interval_finding_t* out) {
    *out = (alea_ray_interval_finding_t){
        .t_enter = t_enter, .t_exit = t_exit,
        .kind = ALEA_INTERVAL_OK, .cell_id = -1,
        .overlap_cell_id = -1, .depth = -1
    };
    if (n == 0) {
        out->kind = ALEA_INTERVAL_GAP;
        return;
    }

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (hits[i].depth == hits[j].depth) {
                out->kind = ALEA_INTERVAL_OVERLAP;
                out->cell_id = hits[i].cell_id;
                out->overlap_cell_id = hits[j].cell_id;
                out->depth = hits[i].depth;
                return;
            }
        }
    }

    int deepest = 0;
    while (deepest + 1 < n &&
           hits[deepest + 1].depth == hits[deepest].depth + 1)
        deepest++;
    out->cell_id = hits[deepest].cell_id;
    out->depth = hits[deepest].depth;
    if (hits[deepest].resolution_flags & ALEA_RESOLVE_UNDEFINED_FILL)
        out->kind = ALEA_INTERVAL_UNDEFINED_FILL;
}

int alea_ray_coverage_classify_reuse_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    alea_raycast_result_t* breakpoint_scratch,
    alea_ray_interval_finding_t* out, size_t max_out) {
    if (!sys || !ray || !breakpoint_scratch || t_max <= 0.0) return -1;
    if (alea_raycast_global_breakpoints_reuse_nocache(
            sys, ray, t_max, breakpoint_scratch) != 0)
        return -1;

    alea_cell_hit_t previous[ALEA_RAY_COVERAGE_MAX_OWNERS];
    int previous_count = -1;
    alea_ray_interval_finding_t current = {0};
    int have_current = 0;
    size_t total = 0;
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
        const int hit_count = alea_find_all_cells_at_point_recursive(
            sys, x, y, z, hits, ALEA_RAY_COVERAGE_MAX_OWNERS);
        if (hit_count < 0)
            return -1;
        if (hit_count >= ALEA_RAY_COVERAGE_MAX_OWNERS) {
            /* The legacy classifier has no truncated finding ABI.  Failing
             * is therefore the only honest result: reporting no overlap here
             * would turn an incomplete diagnostic into apparent success. */
            alea_set_error_detail(ALEA_ERR_OVERFLOW,
                                  "ray coverage owner budget exceeded");
            return -1;
        }

        if (have_current && coverage_owner_sets_equal(
                hits, hit_count, previous, previous_count)) {
            current.t_exit = t_current;
        } else {
            if (have_current) {
                if (total < max_out && out) out[total] = current;
                total++;
            }
            coverage_make_legacy_finding(hits, hit_count, t_previous,
                                         t_current, &current);
            have_current = 1;
            previous_count = hit_count;
            memcpy(previous, hits, (size_t)hit_count * sizeof(*hits));
        }
        t_previous = t_current;
        if (t_previous >= t_max) break;
    }

    if (have_current) {
        if (total < max_out && out) out[total] = current;
        total++;
    }
    return total > (size_t)INT_MAX ? -1 : (int)total;
}
