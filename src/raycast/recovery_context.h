// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_RECOVERY_CONTEXT_H
#define ALEA_RECOVERY_CONTEXT_H

#include "raycast.h"
#include "util/ray_profile.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Optional scratch owned by one ray/walk. Empty-window certificates belong to
 * one accumulating nearest recovery and must be reset between recoveries.
 * Polynomial roots survive those resets, keyed by the exact occurrence ray.
 * Explicit arrays avoid hashing/comparing C struct padding. */
typedef struct {
    const uint32_t* cards;
    size_t count;
    double values[33]; /* ray including reciprocals, transform and inverse */
    unsigned flags;
} recovery_context_key_t;
typedef struct {
    double t;
    uint32_t ordinal;
    uint32_t root;
} recovery_root_t;
#define RECOVERY_ROOT_BUDGET (8u * 1024u * 1024u)
typedef struct {
    recovery_context_key_t key;
    double lo, hi;
    bool complete;
    bool roots_attempted;
    bool roots_complete;
    bool roots_binned;
    uint64_t sorted_bins;
    double roots_lo, roots_hi;
    recovery_root_t* roots;
    size_t root_count;
    uint32_t* exceptional;
    size_t exceptional_count;
} recovery_context_entry_t;
typedef struct {
    size_t capacity;
    size_t root_bytes;
    recovery_context_entry_t entries[];
} recovery_context_t;

static inline recovery_context_t* recovery_context_create(size_t capacity) {
    if (!capacity || (capacity & (capacity-1)) ||
        capacity > (1024*1024-sizeof(recovery_context_t))/sizeof(recovery_context_entry_t))
        return NULL;
    recovery_context_t* table = calloc(1, sizeof(*table)+capacity*sizeof(table->entries[0]));
    if (table) table->capacity = capacity;
    return table;
}

static inline void recovery_context_destroy(recovery_context_t* table) {
    if (!table) return;
    for (size_t i=0;i<table->capacity;i++) {
        free(table->entries[i].roots);
        free(table->entries[i].exceptional);
    }
    free(table);
}

static inline int recovery_root_compare(const void* a, const void* b) {
    const recovery_root_t *x=a, *y=b;
    if (x->t != y->t) return x->t < y->t ? -1 : 1;
    if (x->ordinal != y->ordinal) return x->ordinal < y->ordinal ? -1 : 1;
    return (x->root > y->root) - (x->root < y->root);
}

static inline const recovery_root_t* recovery_sorted_root_first(
    const recovery_root_t* roots, size_t count, double lo, double hi, double excluded) {
    size_t first=0, end=count;
    while (first<end) {
        size_t mid=first+(end-first)/2;
        double t=roots[mid].t;
        if (t<lo || t<=excluded) first=mid+1;
        else end=mid;
    }
    return first<count && roots[first].t<=hi ? &roots[first] : NULL;
}

/* Only the bin assignment is approximate. It is monotone for the finite,
 * nonnegative cache interval. Build and query use the same function, and all
 * comparisons against the actual query window use the original double roots. */
#define RECOVERY_ROOT_BINS 64u
static inline unsigned recovery_root_bin(double t, double lo, double hi) {
    if (t <= lo) return 0;
    if (t >= hi) return RECOVERY_ROOT_BINS-1;
    unsigned bin = (unsigned)(((t-lo)/(hi-lo))*RECOVERY_ROOT_BINS);
    return bin < RECOVERY_ROOT_BINS ? bin : RECOVERY_ROOT_BINS-1;
}

/* A stable scatter retains card/root encounter order within each bin.
 * Store the offsets immediately after the roots in one owned allocation. */
static inline recovery_root_t* recovery_bin_roots(
    const recovery_root_t* roots, size_t count, double lo, double hi,
    size_t budget, size_t* bytes) {
    uint32_t offsets[RECOVERY_ROOT_BINS+1] = {0}, cursor[RECOVERY_ROOT_BINS];
    if (count < RECOVERY_ROOT_BINS || budget < sizeof(offsets) ||
        count > (budget-sizeof(offsets))/sizeof(*roots) || count > UINT32_MAX) return NULL;
    *bytes = count*sizeof(*roots)+sizeof(offsets);
    recovery_root_t* indexed = malloc(*bytes);
    if (!indexed) return NULL;
    for (size_t i=0; i<count; ++i) ++offsets[1+recovery_root_bin(roots[i].t,lo,hi)];
    for (unsigned i=0; i<RECOVERY_ROOT_BINS; ++i) offsets[i+1] += offsets[i];
    memcpy(cursor,offsets,sizeof(cursor));
    for (size_t i=0; i<count; ++i) {
        unsigned bin = recovery_root_bin(roots[i].t,lo,hi);
        indexed[cursor[bin]++] = roots[i];
    }
    memcpy(indexed+count,offsets,sizeof(offsets));
    return indexed;
}

/* Walk bins in increasing distance. Small bins retain the scan's first-card
 * tie order. Sort a busy bin once, only when a lookup actually reaches it. */
static inline const recovery_root_t* recovery_root_first(
    recovery_context_entry_t* entry, double lo, double hi, double current) {
    double excluded = nextafter(current,INFINITY);
    if (!entry->roots_binned)
        return recovery_sorted_root_first(entry->roots,entry->root_count,lo,hi,excluded);
    const uint32_t* offsets = (const uint32_t*)(entry->roots+entry->root_count);
    unsigned first = recovery_root_bin(lo,entry->roots_lo,entry->roots_hi);
    unsigned last = recovery_root_bin(hi,entry->roots_lo,entry->roots_hi);
    for (unsigned bin=first; bin<=last; ++bin) {
        recovery_root_t* roots = entry->roots+offsets[bin];
        size_t count = offsets[bin+1]-offsets[bin];
        const recovery_root_t* best = NULL;
        if (count <= 16) {
            for (size_t i=0; i<count; ++i) {
                double t = roots[i].t;
                if (t>=lo && t<=hi && t>excluded && (!best || t<best->t)) best = &roots[i];
            }
        } else {
            uint64_t bit = UINT64_C(1)<<bin;
            if (!(entry->sorted_bins & bit)) {
                qsort(roots,count,sizeof(*roots),recovery_root_compare);
                RAY_PROFILE_ADD(roots_sorted,count);
                entry->sorted_bins |= bit;
            }
            best = recovery_sorted_root_first(roots,count,lo,hi,excluded);
        }
        if (best) return best;
    }
    return NULL;
}

static inline bool recovery_context_key(recovery_context_key_t* key,
    const uint32_t* cards, size_t count, const alea_ray_t* ray,
    const alea_matrix_t* transform) {
    *key = (recovery_context_key_t){.cards=cards, .count=count};
    const double values[] = {ray->ox,ray->oy,ray->oz,ray->dx,ray->dy,ray->dz,
                            ray->inv_dx,ray->inv_dy,ray->inv_dz};
    for (size_t i=0;i<9;i++) {
        /* Infinite reciprocals of zero directions are ordinary ray data. */
        if (i<6 ? !isfinite(values[i]) : isnan(values[i])) return false;
        key->values[i]=values[i];
    }
    if (transform) {
        key->flags=1u | (transform->has_inverse ? 2u : 0u);
        for (size_t i=0;i<12;i++) {
            if (!isfinite(transform->m[i])) return false;
            key->values[9+i]=transform->m[i];
            if (transform->has_inverse) {
                if (!isfinite(transform->inv[i])) return false;
                key->values[21+i]=transform->inv[i];
            }
        }
    }
    return true;
}

/* Return an existing key or an unused slot. Probe effort is bounded; collisions
 * and saturation fall back to the scan. A wider window never counts as a hit. */
static inline recovery_context_entry_t* recovery_context_find(
    recovery_context_t* table, const recovery_context_key_t* key) {
    if (!table) return NULL;
    uint64_t hash=(uint64_t)(uintptr_t)key->cards ^ (uint64_t)key->count ^ key->flags;
    const unsigned char* bytes=(const unsigned char*)key->values;
    for(size_t i=0;i<sizeof(key->values);i++)
        hash=(hash^bytes[i])*UINT64_C(1099511628211);
    for(size_t probe=0;probe<16 && probe<table->capacity;probe++) {
        recovery_context_entry_t* entry=&table->entries[(hash+probe)&(table->capacity-1)];
        if (!entry->complete || (entry->key.cards==key->cards &&
            entry->key.count==key->count && entry->key.flags==key->flags &&
            !memcmp(entry->key.values,key->values,sizeof(key->values)))) return entry;
    }
    return NULL;
}
static inline bool recovery_context_contains(const recovery_context_entry_t* entry,
                                            double lo, double hi) {
    return entry && entry->complete && isfinite(lo) && isfinite(hi) &&
           lo>=0 && hi>=lo && lo>=entry->lo && hi<=entry->hi;
}
#endif
