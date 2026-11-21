// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_RAY_PROFILE_H
#define ALEA_RAY_PROFILE_H
/* Private, compile-time instrumentation. No code/counters in ordinary builds. */
#ifdef ALEA_RAY_PROFILE
#include "compat.h"
#include <stdint.h>
enum {
    RAY_OTHER, RAY_HIERARCHY, RAY_MEMBERSHIP, RAY_INTERSECTION,
    RAY_RECOVERY, RAY_RECOVERY_SCAN, RAY_RECOVERY_CACHE, RAY_STAGES
};
typedef struct { uint32_t root, used; uint64_t calls; double seconds; } ray_profile_root;
typedef struct {
    const uint32_t *indices;
    size_t count;
    uint64_t calls, tests;
} ray_profile_list;
typedef enum {
    RAY_OWNER_NEIGHBOR, RAY_OWNER_CACHED_PATH, RAY_OWNER_PATH_DESCENT,
    RAY_OWNER_ROOT, RAY_OWNER_FULL
} ray_profile_owner_t;

typedef struct {
    size_t ray;
    int attempt, path_length, root_cell, owner_cell, owner_universe;
    ray_profile_owner_t source;
    uint64_t occurrence;
    double current_t, proposed_t, recovered_t;
    int proposed_surface, recovered_surface;
} ray_profile_recovery_event;

/* Deepest-owner descent is recursive. Self time and membership exclude
 * child frames. Recovery time includes work below the boundary retry. */
typedef struct {
    uint64_t lookups, membership_calls, enumerations, candidates, recoveries;
    double lookup_seconds, membership_seconds, enumeration_seconds, recovery_seconds;
} ray_profile_universe;

#ifndef RAY_PROFILE_RECOVERIES
#define RAY_PROFILE_RECOVERIES 1024u
#endif
#define RAY_PROFILE_LISTS 4096u
#define RAY_PROFILE_TYPES 32u
typedef struct {
    double seconds[RAY_STAGES], last;
    uint64_t calls[RAY_STAGES], csg_nodes, csg_primitives, dropped_roots;
    int enabled, stage;
    ray_profile_universe *universes, *active_universe;
    double universe_start, membership_start;
    uint64_t membership_calls_start;
    ray_profile_root *roots;
    ray_profile_list *lists;
    ray_profile_recovery_event *recoveries;
    size_t ray_ordinal, recovery_count;
    uint64_t dropped_recoveries;
    uint64_t dropped_lists, recovery_cache_hits, plane_index_groups, plane_index_pruned;
    uint64_t recovery_root_builds, recovery_root_tests, recovery_root_queries;
    uint64_t prepared_card_checks_reused, cached_empty_windows;
    size_t peak_root_bytes;
    uint64_t indexed_root_lists, roots_sorted;
    uint64_t fill_passes, fill_improved, fill_children, fill_duplicates;
    uint64_t fill_bbox_rejected, fill_descents;
    uint64_t cached_root_candidates, root_tests_by_type[RAY_PROFILE_TYPES];
    uint64_t roots_retained_by_type[RAY_PROFILE_TYPES];
    uint64_t scan_tests[RAY_PROFILE_TYPES], scan_window_roots[RAY_PROFILE_TYPES];
} ray_profile_state;
extern ALEA_THREAD_LOCAL ray_profile_state alea_ray_profile;
#define RAY_PROFILE_ROOTS 131072u
int alea_ray_profile_start(size_t universe_count);
void alea_ray_profile_switch_universe(ray_profile_universe* next);
void alea_ray_profile_stop(void);
void alea_ray_profile_free(void);
int alea_ray_profile_enter(int stage);
void alea_ray_profile_leave(int previous);
void alea_ray_profile_root_record(uint32_t root, double seconds);
void alea_ray_profile_list_record(const uint32_t *indices, size_t count, uint64_t tests);
#define RAY_PROFILE_COUNT(field) do { if (alea_ray_profile.enabled) ++alea_ray_profile.field; } while(0)
#define RAY_PROFILE_ADD(field,n) do { if (alea_ray_profile.enabled) alea_ray_profile.field += (n); } while(0)
/* Keep error/early-return paths balanced without compiler-specific cleanup. */
#define RAY_PROFILE_FUNCTION(type,name,params,args,stage) \
    static type name##_profile_impl params; \
    type name params { \
        int previous = alea_ray_profile_enter(stage); \
        type result = name##_profile_impl args; \
        alea_ray_profile_leave(previous); \
        return result; \
    } \
    static type name##_profile_impl params
#else
#define RAY_PROFILE_COUNT(field) ((void)0)
#define RAY_PROFILE_ADD(field,n) ((void)0)
#define RAY_PROFILE_FUNCTION(type,name,params,args,stage) type name params
#endif
#endif
