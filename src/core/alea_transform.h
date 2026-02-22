// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_TRANSFORM_H
#define ALEA_TRANSFORM_H

/**
 * @file alea_transform.h
 * @brief Transform types and operations
 *
 * Transform struct (TRn cards + inline FILL), transform vec, and
 * all transform-related API functions.
 */

#include "alea_types.h"
#include "util/alea_vec.h"

#ifdef __cplusplus
extern "C" {
#endif

// ============================================================================
// TRANSFORM TYPE
// ============================================================================

/**
 * @brief Transform information
 *
 * Stores both TRn card transforms and inline transforms from FILL.
 * Inline transforms get auto-assigned IDs starting from a high number
 * to avoid collision with user-defined TRn cards.
 *
 * The transform stores both original values (for export) and computed
 * direction cosines (for use in calculations). Format:
 *   data[0-2]: ox, oy, oz (translation - same in both forms)
 *   data[3-11]: b1-b9 (original values - may be degrees or cosines)
 *   cosines[3-11]: b1-b9 as direction cosines (always ready to use)
 */
typedef struct {
    int transform_id;           // MCNP TRn number
    double data[12];            // Original values (angles if degrees=true, cosines if degrees=false)
    double cosines[12];         // Computed direction cosines (always usable as-is)
    int value_count;            // 3 (translation only) or 12 (full transform)
    bool degrees;               // Original data[3-11] are angles in degrees (vs direction cosines)
    int from_inline;            // 1 = from inline FILL, 0 = from TRn card
} alea_transform_t;

ALEA_VEC_DEFINE(alea_transform_vec, alea_transform_t);

// ============================================================================
// API - TRANSFORM OPERATIONS
// ============================================================================

/**
 * @brief Add a transform from TRn card
 *
 * @param sys CSG system
 * @param transform_id TRn number (e.g., 1 for TR1)
 * @param data Transform data (3 or 12 values)
 * @param value_count Number of values (3 for translation, 12 for full)
 * @param degrees 1 if angles in degrees (*TRn), 0 for direction cosines
 * @return Result with 0 on success, or error code
 */
int alea_add_transform(alea_system_t* sys, int transform_id,
                      const double* data, int value_count, int degrees);

/**
 * @brief Add an inline transform (from FILL) and assign an ID
 *
 * @param sys CSG system
 * @param data Transform data (3 or 12 values)
 * @param value_count Number of values (3 for translation, 12 for full)
 * @param degrees 1 if angles in degrees (*FILL), 0 for direction cosines
 * @return Assigned transform ID, or -1 on error
 */
int alea_add_inline_transform(alea_system_t* sys, const double* data,
                             int value_count, int degrees);

/**
 * @brief Get a transform by ID
 *
 * @param sys CSG system
 * @param transform_id Transform ID
 * @return Pointer to transform, or NULL if not found
 */
const alea_transform_t* alea_get_transform(const alea_system_t* sys, int transform_id);

/**
 * @brief Apply transform to a point (auxiliary -> main coordinates)
 *
 * @param tr Transform to apply
 * @param x X coordinate (modified in place)
 * @param y Y coordinate (modified in place)
 * @param z Z coordinate (modified in place)
 */
void alea_transform_point(const alea_transform_t* tr, double* x, double* y, double* z);

/**
 * @brief Apply transform to a direction vector (rotation only, no translation)
 *
 * @param tr Transform to apply
 * @param vx X component (modified in place)
 * @param vy Y component (modified in place)
 * @param vz Z component (modified in place)
 */
void alea_transform_vector(const alea_transform_t* tr, double* vx, double* vy, double* vz);

/**
 * @brief Apply transform to primitive, converting to global coordinates
 *
 * Transforms the primitive geometry from auxiliary (local) coordinates
 * to main (global) coordinates. May change the primitive type if the
 * transform converts an axis-aligned primitive to a general form.
 *
 * @param tr Transform to apply
 * @param in_type Input primitive type
 * @param in_data Input primitive data
 * @param out_type Output primitive type (may differ from input)
 * @param out_data Output primitive data (in global coordinates)
 * @return true on success, false if transform cannot be applied
 */
bool alea_apply_transform_to_primitive(const alea_transform_t* tr,
                                      alea_primitive_type_t in_type,
                                      const alea_primitive_data_t* in_data,
                                      alea_primitive_type_t* out_type,
                                      alea_primitive_data_t* out_data);

/**
 * @brief Apply inverse transform to primitive (for export roundtrip)
 *
 * Converts primitive from global coordinates back to local coordinates
 * by applying the inverse of the given transform.
 *
 * @param tr Transform whose inverse should be applied
 * @param in_type Input primitive type (in global coordinates)
 * @param in_data Input primitive data
 * @param out_type Output primitive type (may differ from input)
 * @param out_data Output primitive data (in local coordinates)
 * @return true on success, false if inverse transform cannot be applied
 */
bool alea_apply_inverse_transform_to_primitive(const alea_transform_t* tr,
                                              alea_primitive_type_t in_type,
                                              const alea_primitive_data_t* in_data,
                                              alea_primitive_type_t* out_type,
                                              alea_primitive_data_t* out_data);

/**
 * @brief Check if transform is translation-only (no rotation)
 *
 * @param tr Transform to check
 * @return true if translation only, false if includes rotation
 */
bool alea_transform_is_translation_only(const alea_transform_t* tr);

/**
 * @brief Finalize inline transform ID counter after loading TRn cards
 *
 * Scans existing transforms and sets next_inline_transform_id = max(transform_ids) + 1.
 * Call this after loading all TRn cards and before cell conversion begins.
 *
 * @param sys CSG system
 */
void alea_finalize_transform_ids(alea_system_t* sys);

#ifdef __cplusplus
}
#endif

#endif /* ALEA_TRANSFORM_H */
