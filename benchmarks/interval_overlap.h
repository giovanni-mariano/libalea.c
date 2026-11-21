// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_BENCH_INTERVAL_OVERLAP_H
#define ALEA_BENCH_INTERVAL_OVERLAP_H

#include "alea_raycast.h"
#include "core/alea_system.h"
#include "raycast/raycast.h"
#include <math.h>

/* Diagnostic scoring for finite prepared/public ray replays. */
typedef struct {
    int cell, material;
    double density;
    uint8_t flags;
} replay_owner_t;

typedef struct {
    alea_system_t* sys;
    replay_owner_t prepared, legacy;
    double covered;
} overlap_check_t;

static bool replay_owner_matches(const overlap_check_t* check, replay_owner_t label,
                                 const alea_ray_coverage_owner_t* owner) {
    return owner && !owner->resolution_flags && !label.flags &&
           label.cell == owner->cell_id && label.material == owner->material_id &&
           label.density == check->sys->cells.data[owner->cell_index].density;
}

/* The strict sweep supplies the whole occurrence tree, including containers.
 * Select the earliest sibling at each level, never a global minimum cell ID.
 * The legacy result must also identify a terminal claimant of this interval. */
static int check_overlap_interval(void* context, const alea_ray_coverage_interval_t* span) {
    overlap_check_t* check = context;
    if (span->kind != ALEA_RAY_COVERAGE_OVERLAP || span->t_enter != check->covered) return 1;
    uint64_t parent = 0;
    const alea_ray_coverage_owner_t* selected = NULL;
    for (size_t depth = 0; depth <= span->owner_count; ++depth) {
        const alea_ray_coverage_owner_t* first = NULL;
        for (size_t i = 0; i < span->owner_count; ++i) {
            const alea_ray_coverage_owner_t* o = &span->owners[i];
            if (o->parent_occurrence_key == parent && (!first || o->cell_index < first->cell_index))
                first = o;
        }
        if (!first) break;
        if (depth == span->owner_count) return 1;
        selected = first;
        parent = first->occurrence_key;
    }
    if (!replay_owner_matches(check, check->prepared, selected)) return 1;
    for (size_t i = 0; i < span->owner_count; ++i) {
        const alea_ray_coverage_owner_t* o = &span->owners[i];
        if (!replay_owner_matches(check, check->legacy, o)) continue;
        size_t j = 0;
        while (j < span->owner_count && span->owners[j].parent_occurrence_key != o->occurrence_key) ++j;
        if (j == span->owner_count) {
            check->covered = span->t_exit;
            return 0;
        }
    }
    return 1;
}

/* Overlay the two partitions of the original ray. The ordinary tolerance
 * test has already failed. Account for every disagreement before excusing
 * the ray, including stretches before or after the first shifted row. */
static bool overlap_only_ray(alea_system_t* sys, const alea_raycast_batch_result_t* batch,
                             size_t begin, size_t end, const alea_raycast_result_t* prepared,
                             double distance, alea_raycast_result_t* scratch) {
    size_t i = begin, j = 0;
    double t = 0;
    bool found = false;
    while (i < end && j < prepared->segments.count) {
        const alea_ray_segment_t* s = &prepared->segments.data[j];
        double enter = alea_raycast_batch_t_enter(batch)[i];
        double exit = alea_raycast_batch_t_exit(batch)[i];
        double hi = fmin(exit, s->t_exit);
        if (!isfinite(exit) || !isfinite(s->t_exit) || !(hi > t) ||
            enter != (i == begin ? 0 : alea_raycast_batch_t_exit(batch)[i-1]) ||
            s->t_enter != (j == 0 ? 0 : prepared->segments.data[j-1].t_exit)) return false;
        replay_owner_t p = {s->cell_id, s->material_id, s->density, s->resolution_flags};
        replay_owner_t h = {alea_raycast_batch_cell_ids(batch)[i], alea_raycast_batch_material_ids(batch)[i],
                            alea_raycast_batch_densities(batch)[i], alea_raycast_batch_resolution_flags(batch)[i]};
        if (p.cell == h.cell) {
            if (p.material != h.material || p.density != h.density || p.flags != h.flags) return false;
        } else {
            overlap_check_t check = {.sys = sys, .prepared = p, .legacy = h, .covered = t};
            int rc = alea_ray_coverage_sweep_strict_reuse_nocache(sys, &prepared->ray, t, hi,
                8192, NULL, scratch, check_overlap_interval, &check);
            if (rc <= 0 || check.covered != hi) return false;
            found = true;
        }
        /* A spurious split between identical cells is a separate structural
         * mismatch, even if this ray also traverses a genuine overlap. */
        bool public_split = exit == hi && i+1 < end && alea_raycast_batch_cell_ids(batch)[i+1] == h.cell;
        bool prepared_split = s->t_exit == hi && j+1 < prepared->segments.count &&
                              prepared->segments.data[j+1].cell_id == p.cell;
        if (public_split != prepared_split) return false;
        if (exit == hi) ++i;
        if (s->t_exit == hi) ++j;
        t = hi;
    }
    return found && i == end && j == prepared->segments.count && t == distance;
}

#endif
