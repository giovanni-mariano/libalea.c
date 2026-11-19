// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#ifndef ALEA_FINGERPRINT_H
#define ALEA_FINGERPRINT_H
#include "alea.h"
#ifdef __cplusplus
extern "C" {
#endif
#define ALEA_FINGERPRINT_VERSION 2u
typedef struct { unsigned char bytes[32]; } alea_digest_t;
typedef struct { int id; alea_digest_t geometry,content; } alea_fingerprint_t;
typedef struct {
    unsigned version;
    size_t cell_count,universe_count;
    alea_fingerprint_t* cells; /* Model cell order. */
    alea_fingerprint_t* universes; /* Sorted by universe number. */
} alea_fingerprints_t;
/** Exact canonical fingerprints, independent of IDs and commutative term order.
 * Equality covers the supported canonical rewrites, not arbitrary CSG equivalence.
 * Bounds, comments, names and wrapper metadata are excluded. Inputs must not be
 * mutated concurrently. Cyclic fills other than lattice self-fill are errors. */
alea_fingerprints_t* alea_fingerprints(const alea_system_t* sys);
void alea_fingerprints_destroy(alea_fingerprints_t* fingerprints);
void alea_digest_hex(const alea_digest_t* digest,char hex[65]);
typedef enum {
    ALEA_DIFF_UNCHANGED, ALEA_DIFF_GEOMETRY, ALEA_DIFF_CONTENT,
    ALEA_DIFF_RENUMBERED, ALEA_DIFF_REMOVED, ALEA_DIFF_ADDED,
    ALEA_DIFF_AMBIGUOUS
} alea_diff_kind_t;
typedef struct { int universe; int old_id,new_id; alea_diff_kind_t kind; } alea_model_change_t;
typedef struct { size_t count; alea_model_change_t* changes; } alea_model_diff_t;
/** Includes unchanged entries. Ambiguous removed matches have new_id=-1;
 * their unmatched destinations are also reported as added entries. */
alea_model_diff_t* alea_fingerprint_diff(const alea_fingerprints_t* before,
                                       const alea_fingerprints_t* after);
void alea_model_diff_destroy(alea_model_diff_t* diff);
typedef enum {
    ALEA_FINGERPRINT_LEGACY = 0,
    ALEA_FINGERPRINT_SEMANTIC = 1
} alea_fingerprint_profile_t;
typedef struct {
    size_t struct_size;
    alea_fingerprint_profile_t profile;
} alea_fingerprint_options_t;
typedef struct alea_fingerprint_set alea_fingerprint_set_t;
/** Opt-in semantic comparison normalizes material/mixture fractions and merges
 * duplicate components. No numeric quantization is performed. This is not an
 * acceleration reuse key. The source remains unchanged. */
alea_fingerprint_set_t* alea_fingerprint_set_create(const alea_system_t* sys,
    const alea_fingerprint_options_t* options);
/** Borrowed, read-only values, valid until set destruction. */
const alea_fingerprints_t* alea_fingerprint_set_values(const alea_fingerprint_set_t* set);
alea_fingerprint_profile_t alea_fingerprint_set_profile(const alea_fingerprint_set_t* set);
void alea_fingerprint_set_destroy(alea_fingerprint_set_t* set);
/** Rejects incompatible profiles. */
alea_model_diff_t* alea_fingerprint_set_diff(const alea_fingerprint_set_t* before,
    const alea_fingerprint_set_t* after);
#ifdef __cplusplus
}
#endif
#endif
