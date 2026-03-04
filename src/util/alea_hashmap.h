// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_hashmap.h
 * @brief Header-only macro-based generic hash map (open addressing, linear probing)
 *
 * Uses backward-shift deletion (no tombstones). Grows at load > 0.7.
 *
 * Usage:
 *   // Define hash and equality functions
 *   static inline uint64_t my_hash(int key) { return (uint64_t)key * 0x9E3779B97F4A7C15ULL; }
 *   static inline bool my_eq(int a, int b) { return a == b; }
 *
 *   // Generate the type and functions (empty_key is a sentinel never used as real key)
 *   ALEA_HASHMAP_DEFINE(my_map, int, double, my_hash, my_eq, INT_MIN)
 *
 *   // Use it
 *   my_map_t map = my_map_create(16);
 *   my_map_put(&map, 42, 3.14);
 *   double* val = my_map_get(&map, 42);  // returns pointer or NULL
 *   my_map_destroy(&map);
 */

#ifndef ALEA_HASHMAP_H
#define ALEA_HASHMAP_H

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

/**
 * Generate a typed hash map with inline functions.
 *
 * @param name      Prefix for generated types/functions
 * @param key_type  Key type
 * @param val_type  Value type
 * @param hash_fn   uint64_t hash_fn(key_type key)
 * @param eq_fn     bool eq_fn(key_type a, key_type b)
 * @param empty_key Sentinel value for empty slots (must never be a valid key)
 */
#define ALEA_HASHMAP_DEFINE(name, key_type, val_type, hash_fn, eq_fn, empty_key) \
\
typedef struct { key_type key; val_type value; } name##_entry_t; \
\
typedef struct { \
    name##_entry_t* entries; \
    size_t capacity; \
    size_t count; \
} name##_t; \
\
static inline name##_t name##_create(size_t initial_cap) { \
    name##_t map; \
    size_t _cap = initial_cap < 8 ? 8 : initial_cap; \
    /* Round up to next power of two */ \
    _cap--; _cap |= _cap >> 1; _cap |= _cap >> 2; \
    _cap |= _cap >> 4; _cap |= _cap >> 8; _cap |= _cap >> 16; \
    _cap |= _cap >> 32; _cap++; \
    map.capacity = _cap; \
    map.count = 0; \
    map.entries = (name##_entry_t*)malloc(map.capacity * sizeof(name##_entry_t)); \
    if (map.entries) { \
        for (size_t _i = 0; _i < map.capacity; _i++) { \
            map.entries[_i].key = (empty_key); \
        } \
    } else { \
        map.capacity = 0; \
    } \
    return map; \
} \
\
static inline void name##_destroy(name##_t* map) { \
    free(map->entries); \
    map->entries = NULL; \
    map->capacity = 0; \
    map->count = 0; \
} \
\
static inline void name##_clear(name##_t* map) { \
    for (size_t _i = 0; _i < map->capacity; _i++) { \
        map->entries[_i].key = (empty_key); \
    } \
    map->count = 0; \
} \
\
static inline val_type* name##_get(const name##_t* map, key_type key) { \
    if (!map->entries || map->count == 0) return NULL; \
    size_t mask = map->capacity - 1; \
    size_t idx = (size_t)(hash_fn(key)) & mask; \
    for (size_t _probe = 0; _probe < map->capacity; _probe++) { \
        name##_entry_t* e = &map->entries[idx]; \
        if (eq_fn(e->key, (empty_key))) return NULL; \
        if (eq_fn(e->key, key)) return &e->value; \
        idx = (idx + 1) & mask; \
    } \
    return NULL; \
} \
\
static inline bool name##_resize(name##_t* map, size_t new_cap) { \
    name##_entry_t* old = map->entries; \
    size_t old_cap = map->capacity; \
    map->entries = (name##_entry_t*)malloc(new_cap * sizeof(name##_entry_t)); \
    if (!map->entries) { map->entries = old; return false; } \
    map->capacity = new_cap; \
    map->count = 0; \
    for (size_t _i = 0; _i < new_cap; _i++) { \
        map->entries[_i].key = (empty_key); \
    } \
    /* Re-insert old entries */ \
    size_t mask = new_cap - 1; \
    for (size_t _i = 0; _i < old_cap; _i++) { \
        if (!eq_fn(old[_i].key, (empty_key))) { \
            size_t idx = (size_t)(hash_fn(old[_i].key)) & mask; \
            while (!eq_fn(map->entries[idx].key, (empty_key))) { \
                idx = (idx + 1) & mask; \
            } \
            map->entries[idx] = old[_i]; \
            map->count++; \
        } \
    } \
    free(old); \
    return true; \
} \
\
static inline bool name##_put(name##_t* map, key_type key, val_type value) { \
    if (!map->entries) return false; \
    /* Grow if load > 70% */ \
    if (map->count * 10 >= map->capacity * 7) { \
        if (!name##_resize(map, map->capacity * 2)) return false; \
    } \
    size_t mask = map->capacity - 1; \
    size_t idx = (size_t)(hash_fn(key)) & mask; \
    for (;;) { \
        name##_entry_t* e = &map->entries[idx]; \
        if (eq_fn(e->key, (empty_key))) { \
            e->key = key; \
            e->value = value; \
            map->count++; \
            return true; \
        } \
        if (eq_fn(e->key, key)) { \
            e->value = value; /* Update existing */ \
            return true; \
        } \
        idx = (idx + 1) & mask; \
    } \
} \
\
static inline bool name##_remove(name##_t* map, key_type key) { \
    if (!map->entries || map->count == 0) return false; \
    size_t mask = map->capacity - 1; \
    size_t idx = (size_t)(hash_fn(key)) & mask; \
    /* Find the entry */ \
    for (size_t _probe = 0; _probe < map->capacity; _probe++) { \
        if (eq_fn(map->entries[idx].key, (empty_key))) return false; \
        if (eq_fn(map->entries[idx].key, key)) break; \
        idx = (idx + 1) & mask; \
        if (_probe == map->capacity - 1) return false; \
    } \
    /* Backward-shift deletion */ \
    map->entries[idx].key = (empty_key); \
    map->count--; \
    size_t next = (idx + 1) & mask; \
    while (!eq_fn(map->entries[next].key, (empty_key))) { \
        size_t natural = (size_t)(hash_fn(map->entries[next].key)) & mask; \
        /* Check if 'next' is displaced past 'idx' */ \
        bool displaced; \
        if (next >= idx) { \
            displaced = (natural <= idx || natural > next); \
        } else { /* wrapped */ \
            displaced = (natural <= idx && natural > next); \
        } \
        if (displaced) { \
            map->entries[idx] = map->entries[next]; \
            map->entries[next].key = (empty_key); \
            idx = next; \
        } \
        next = (next + 1) & mask; \
    } \
    return true; \
} \
\
/* Sentinel to allow trailing backslash in macro */ \
typedef name##_t name##_hashmap_defined_t

#endif /* ALEA_HASHMAP_H */
