// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef PRIMITIVE_DESC_H
#define PRIMITIVE_DESC_H

/**
 * @file primitive_desc.h
 * @brief Primitive descriptor table for data-driven dispatch
 *
 * Instead of switch statements scattered throughout the codebase,
 * we use a descriptor table that maps each primitive type to its
 * operations. This makes adding new primitives straightforward and
 * keeps all primitive-specific code in one place.
 */

#include "alea_types.h"


/* ============================================================================
 * INTERVAL ARITHMETIC TYPES
 * ============================================================================ */

/**
 * @brief Interval for bounding function values over a region
 */
typedef struct {
    double min;
    double max;
} alea_interval_t;

/**
 * @brief Relation of bounding box to surface
 */
typedef enum {
    ALEA_RELATION_NEGATIVE,   /* f(x) < 0 everywhere in box (inside surface) */
    ALEA_RELATION_POSITIVE,   /* f(x) > 0 everywhere in box (outside surface) */
    ALEA_RELATION_INTERSECT   /* f(x) changes sign in box (crosses surface) */
} alea_box_relation_t;

/* Function pointer types for primitive operations */
typedef double (*prim_eval_fn)(const alea_primitive_data_t* data, double x, double y, double z);
typedef alea_bbox_t (*prim_bbox_fn)(const alea_primitive_data_t* data);
typedef alea_interval_t (*prim_interval_eval_fn)(const alea_primitive_data_t* data, const alea_bbox_t* box);

/**
 * @brief Transform function pointer type
 *
 * Transforms a primitive by a 3x4 affine matrix.
 * May change the primitive type (e.g., axis-aligned cylinder -> quadric).
 *
 * @param in_data   Input primitive data
 * @param mat       3x4 affine matrix as 12 doubles, row-major: [R00,R01,R02,Tx, R10,R11,R12,Ty, R20,R21,R22,Tz]
 * @param out_type  Output primitive type (may differ from input)
 * @param out_data  Output primitive data
 * @return true if transform succeeded, false on error
 */
typedef bool (*prim_transform_fn)(
    const alea_primitive_data_t* in_data,
    const double* mat,
    alea_primitive_type_t* out_type,
    alea_primitive_data_t* out_data
);

/**
 * @brief Descriptor for a primitive type
 *
 * Contains all operations that can be performed on a primitive.
 * NULL pointers indicate the operation is not supported.
 */
typedef struct {
    const char* name;                   /* Human-readable name ("plane", "sphere", etc.) */
    prim_eval_fn eval;                  /* Signed distance evaluation */
    prim_bbox_fn bbox;                  /* Bounding box computation */
    prim_interval_eval_fn interval_eval;/* Interval evaluation over bounding box */
    prim_transform_fn transform;        /* Affine transformation */
} alea_primitive_desc_t;

/**
 * @brief Get descriptor for a primitive type
 * @param type Primitive type enum value
 * @return Pointer to descriptor, or NULL if type is invalid
 */
const alea_primitive_desc_t* alea_primitive_get_desc(alea_primitive_type_t type);

/**
 * @brief Get primitive type name
 * @param type Primitive type enum value
 * @return Name string or "unknown"
 */
const char* alea_primitive_type_name(alea_primitive_type_t type);

/**
 * @brief Transform a primitive by an affine matrix
 *
 * @param in_type   Input primitive type
 * @param in_data   Input primitive data
 * @param mat       3x4 affine matrix (12 doubles, row-major)
 * @param out_type  Output primitive type (may change, e.g. cylinder -> quadric)
 * @param out_data  Output primitive data
 * @return true on success, false if transform not supported for this type
 */
bool alea_primitive_transform(
    alea_primitive_type_t in_type,
    const alea_primitive_data_t* in_data,
    const double* mat,
    alea_primitive_type_t* out_type,
    alea_primitive_data_t* out_data
);

/* ============================================================================
 * INTERVAL ARITHMETIC OPERATIONS
 * ============================================================================ */

/**
 * @brief Create an interval
 */
static inline alea_interval_t alea_iv_make(double min, double max) {
    return (alea_interval_t){min, max};
}

/**
 * @brief Add two intervals: [a,b] + [c,d] = [a+c, b+d]
 */
static inline alea_interval_t alea_iv_add(alea_interval_t a, alea_interval_t b) {
    return (alea_interval_t){a.min + b.min, a.max + b.max};
}

/**
 * @brief Subtract intervals: [a,b] - [c,d] = [a-d, b-c]
 */
static inline alea_interval_t alea_iv_sub(alea_interval_t a, alea_interval_t b) {
    return (alea_interval_t){a.min - b.max, a.max - b.min};
}

/**
 * @brief Multiply interval by scalar
 */
static inline alea_interval_t alea_iv_mul_scalar(alea_interval_t a, double s) {
    if (s >= 0) return (alea_interval_t){a.min * s, a.max * s};
    return (alea_interval_t){a.max * s, a.min * s};
}

/**
 * @brief Multiply two intervals
 */
static inline alea_interval_t alea_iv_mul(alea_interval_t a, alea_interval_t b) {
    double p1 = a.min * b.min;
    double p2 = a.min * b.max;
    double p3 = a.max * b.min;
    double p4 = a.max * b.max;

    double min_v = p1;
    if (p2 < min_v) min_v = p2;
    if (p3 < min_v) min_v = p3;
    if (p4 < min_v) min_v = p4;

    double max_v = p1;
    if (p2 > max_v) max_v = p2;
    if (p3 > max_v) max_v = p3;
    if (p4 > max_v) max_v = p4;

    return (alea_interval_t){min_v, max_v};
}

/**
 * @brief Square an interval with special handling for zero crossing
 */
static inline alea_interval_t alea_iv_sqr(alea_interval_t a) {
    double s1 = a.min * a.min;
    double s2 = a.max * a.max;
    /* If range crosses zero, min square is 0 */
    if (a.min <= 0 && a.max >= 0) {
        return (alea_interval_t){0.0, (s1 > s2) ? s1 : s2};
    }
    return (alea_interval_t){(s1 < s2) ? s1 : s2, (s1 > s2) ? s1 : s2};
}

/**
 * @brief Evaluate interval of primitive surface function over bounding box
 * @param type Primitive type
 * @param data Primitive parameters
 * @param box Bounding box to evaluate over
 * @return Interval bounding the function values
 */
alea_interval_t alea_primitive_interval_eval(
    alea_primitive_type_t type,
    const alea_primitive_data_t* data,
    const alea_bbox_t* box
);


#endif /* PRIMITIVE_DESC_H */
