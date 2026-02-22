// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef PRIMITIVE_EVAL_H
#define PRIMITIVE_EVAL_H

#include "alea_types.h"
#include <stddef.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file primitive_eval.h
 * @brief Signed distance evaluation for all primitive types
 * 
 * This module provides efficient evaluation of signed distance functions
 * for geometric primitives used in CSG operations.
 * 
 * Signed distance convention:
 * - Negative values: point is inside the primitive
 * - Positive values: point is outside the primitive
 * - Zero: point is exactly on the surface
 * 
 * This convention allows CSG operations to work naturally:
 * - Union: min(d1, d2)
 * - Intersection: max(d1, d2)
 * - Difference: max(d1, -d2)
 */

/**
 * Evaluate signed distance for a primitive at a point
 * 
 * @param type Type of primitive (PLANE, SPHERE, etc.)
 * @param data Primitive data (must match type)
 * @param x X coordinate
 * @param y Y coordinate
 * @param z Z coordinate
 * @return Signed distance (negative = inside, positive = outside)
 * 
 * Example:
 * ```c
 * alea_primitive_data_t sphere;
 * sphere.sphere.center_x = 0.0;
 * sphere.sphere.center_y = 0.0;
 * sphere.sphere.center_z = 0.0;
 * sphere.sphere.radius = 5.0;
 * 
 * double d = alea_primitive_eval(ALEA_PRIMITIVE_SPHERE, &sphere, 1.0, 2.0, 3.0);
 * if (d < 0) {
 *     printf("Point is inside sphere\n");
 * }
 * ```
 */
double alea_primitive_eval(
    alea_primitive_type_t type,
    const alea_primitive_data_t* data,
    double x, double y, double z
);

/**
 * Evaluate primitive with sense (orientation)
 * 
 * @param type Type of primitive
 * @param data Primitive data
 * @param sense +1 for normal, -1 for inverted
 * @param x X coordinate
 * @param y Y coordinate
 * @param z Z coordinate
 * @return Signed distance with sense applied
 * 
 * Sense allows primitives to be "flipped" without modifying data.
 * For example, a plane with sense=-1 represents the opposite half-space.
 */
double alea_primitive_eval_with_sense(
    alea_primitive_type_t type,
    const alea_primitive_data_t* data,
    int8_t sense,
    double x, double y, double z
);

/**
 * Check if a point is inside a primitive (convenience function)
 * 
 * @param type Type of primitive
 * @param data Primitive data
 * @param x X coordinate
 * @param y Y coordinate
 * @param z Z coordinate
 * @return true if point is inside (distance < 0)
 */
static inline bool alea_primitive_contains(
    alea_primitive_type_t type,
    const alea_primitive_data_t* data,
    double x, double y, double z
) {
    return alea_primitive_eval(type, data, x, y, z) < 0.0;
}

#ifdef __cplusplus
}
#endif

#endif // PRIMITIVE_EVAL_H