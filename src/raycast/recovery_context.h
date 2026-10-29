// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_RECOVERY_CONTEXT_H
#define ALEA_RECOVERY_CONTEXT_H

#include "raycast.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Optional scratch owned by a single nearest recovery and accumulating result.
 * Explicit arrays avoid hashing/comparing C struct padding. */
typedef struct {
    const uint32_t* cards;
    size_t count;
    double values[33]; /* ray including reciprocals, transform and inverse */
    unsigned flags;
} recovery_context_key_t;
typedef struct {
    recovery_context_key_t key;
    double lo, hi;
    bool complete;
} recovery_context_entry_t;
typedef struct {
    size_t capacity;
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
