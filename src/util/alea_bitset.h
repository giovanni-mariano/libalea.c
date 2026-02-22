// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_bitset.h
 * @brief Header-only compact bitset backed by uint64_t words
 *
 * Replaces bool* calloc patterns with 8x less memory.
 *
 * Usage:
 *   alea_bitset_t bs = alea_bitset_create(1000);
 *   alea_bitset_set(&bs, 42);
 *   if (alea_bitset_test(&bs, 42)) { ... }
 *   alea_bitset_destroy(&bs);
 */

#ifndef ALEA_BITSET_H
#define ALEA_BITSET_H

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    uint64_t* words;
    size_t    n_words;
} alea_bitset_t;

/** Create a bitset with capacity for n_bits, all cleared. Returns zeroed struct on OOM. */
static inline alea_bitset_t alea_bitset_create(size_t n_bits) {
    alea_bitset_t bs;
    bs.n_words = (n_bits + 63) / 64;
    if (bs.n_words == 0) bs.n_words = 1;
    bs.words = (uint64_t*)calloc(bs.n_words, sizeof(uint64_t));
    if (!bs.words) {
        bs.n_words = 0;
    }
    return bs;
}

/** Free bitset memory. */
static inline void alea_bitset_destroy(alea_bitset_t* bs) {
    free(bs->words);
    bs->words = NULL;
    bs->n_words = 0;
}

/** Set bit i. */
static inline void alea_bitset_set(alea_bitset_t* bs, size_t i) {
    bs->words[i >> 6] |= (1ULL << (i & 63));
}

/** Clear bit i. */
static inline void alea_bitset_clear(alea_bitset_t* bs, size_t i) {
    bs->words[i >> 6] &= ~(1ULL << (i & 63));
}

/** Test bit i. */
static inline bool alea_bitset_test(const alea_bitset_t* bs, size_t i) {
    return (bs->words[i >> 6] & (1ULL << (i & 63))) != 0;
}

/** Clear all bits. */
static inline void alea_bitset_clear_all(alea_bitset_t* bs) {
    memset(bs->words, 0, bs->n_words * sizeof(uint64_t));
}

/** Count set bits (popcount). */
static inline size_t alea_bitset_popcount(const alea_bitset_t* bs) {
    size_t count = 0;
    for (size_t i = 0; i < bs->n_words; i++) {
        count += (size_t)__builtin_popcountll(bs->words[i]);
    }
    return count;
}

#endif /* ALEA_BITSET_H */
