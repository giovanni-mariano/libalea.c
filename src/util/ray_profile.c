// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "ray_profile.h"
#ifdef ALEA_RAY_PROFILE
#include <stdlib.h>
#include <string.h>

ALEA_THREAD_LOCAL ray_profile_state alea_ray_profile;

void alea_ray_profile_free(void) {
    free(alea_ray_profile.roots);
    free(alea_ray_profile.lists);
    free(alea_ray_profile.recoveries);
    memset(&alea_ray_profile, 0, sizeof(alea_ray_profile));
}

int alea_ray_profile_start(void) {
    alea_ray_profile_free();
    alea_ray_profile.roots = calloc(RAY_PROFILE_ROOTS, sizeof(ray_profile_root));
    alea_ray_profile.lists = calloc(RAY_PROFILE_LISTS, sizeof(ray_profile_list));
    alea_ray_profile.recoveries = calloc(RAY_PROFILE_RECOVERIES, sizeof(ray_profile_recovery_event));
    if (!alea_ray_profile.roots || !alea_ray_profile.lists || !alea_ray_profile.recoveries) {
        alea_ray_profile_free();
        return -1;
    }
    alea_ray_profile.enabled = 1;
    alea_ray_profile.last = alea_monotonic_seconds();
    return 0;
}

int alea_ray_profile_enter(int stage) {
    if (!alea_ray_profile.enabled) return -1;
    int previous = alea_ray_profile.stage;
    ++alea_ray_profile.calls[stage];
    if (previous != stage) {
        double now = alea_monotonic_seconds();
        alea_ray_profile.seconds[previous] += now - alea_ray_profile.last;
        alea_ray_profile.last = now;
        alea_ray_profile.stage = stage;
    }
    return previous;
}

void alea_ray_profile_leave(int previous) {
    if (previous < 0 || previous == alea_ray_profile.stage) return;
    double now = alea_monotonic_seconds();
    alea_ray_profile.seconds[alea_ray_profile.stage] += now - alea_ray_profile.last;
    alea_ray_profile.last = now;
    alea_ray_profile.stage = previous;
}

void alea_ray_profile_stop(void) {
    if (!alea_ray_profile.enabled) return;
    double now = alea_monotonic_seconds();
    alea_ray_profile.seconds[alea_ray_profile.stage] += now - alea_ray_profile.last;
    alea_ray_profile.enabled = 0;
}

void alea_ray_profile_root_record(uint32_t root, double seconds) {
    if (!alea_ray_profile.enabled) return;
    unsigned slot = (root * UINT32_C(2654435761)) & (RAY_PROFILE_ROOTS - 1);
    /* Bound profiler overhead even when the table is saturated. */
    for (unsigned n = 0; n < 32; ++n, slot = (slot + 1) & (RAY_PROFILE_ROOTS - 1)) {
        ray_profile_root *r = &alea_ray_profile.roots[slot];
        if (!r->used || r->root == root) {
            r->used = 1;
            r->root = root;
            ++r->calls;
            r->seconds += seconds;
            return;
        }
    }
    ++alea_ray_profile.dropped_roots;
}

void alea_ray_profile_list_record(const uint32_t *indices, size_t count, uint64_t tests) {
    if (!alea_ray_profile.enabled) return;
    unsigned slot = ((uintptr_t)indices >> 4) & (RAY_PROFILE_LISTS - 1);
    for (unsigned n = 0; n < 32; ++n, slot = (slot + 1) & (RAY_PROFILE_LISTS - 1)) {
        ray_profile_list *r = &alea_ray_profile.lists[slot];
        if (!r->calls || (r->indices == indices && r->count == count)) {
            r->indices = indices;
            r->count = count;
            ++r->calls;
            r->tests += tests;
            return;
        }
    }
    ++alea_ray_profile.dropped_lists;
}
#endif
