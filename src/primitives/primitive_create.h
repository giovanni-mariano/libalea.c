// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file primitive_create.h
 * @brief Primitive factory functions
 *
 * Updated to use alea_world_t unified context.
 *
 * This module provides TWO levels of primitive creation:
 *
 * LEVEL 1: Data Creation (alea_make_*)
 *   - Creates primitive data structures only
 *   - No pool allocation
 *   - Returns alea_primitive_data_t
 *   - Use when you want to manually manage nodes
 *
 * LEVEL 2: Node Creation (alea_create_*_node)
 *   - Creates complete CSG primitive nodes
 *   - Stores primitives in the primitive pool (with deduplication)
 *   - Allocates nodes from the node pool
 *   - Sets up all node properties (material, bbox, etc.)
 *   - Returns alea_node_id_t
 *   - Convenience functions for common use cases
 */

#ifndef PRIMITIVE_CREATE_H
#define PRIMITIVE_CREATE_H

#include "core/alea_system.h"
#include "bbox.h"

#ifdef __cplusplus
extern "C" {
#endif

// ============================================================================
// LEVEL 1: PRIMITIVE DATA CREATION
// ============================================================================
// These functions create primitive data structures only (no node allocation)
// These are UNCHANGED - they just create data

/**
 * Create a plane: ax + by + cz + d = 0
 *
 * @param a X coefficient
 * @param b Y coefficient
 * @param c Z coefficient
 * @param d Constant term
 * @return Plane primitive data
 *
 * Note: Coefficients should typically be normalized (a²+b²+c²=1)
 * for proper distance calculations
 */
alea_primitive_data_t alea_make_plane(double a, double b, double c, double d);

/**
 * Create a sphere: (x-cx)² + (y-cy)² + (z-cz)² = r²
 *
 * @param cx Center X
 * @param cy Center Y
 * @param cz Center Z
 * @param radius Sphere radius (must be > 0)
 * @return Sphere primitive data
 */
alea_primitive_data_t alea_make_sphere(double cx, double cy, double cz, double radius);

/**
 * Create an axis-aligned box
 *
 * @param x_min Minimum X
 * @param x_max Maximum X
 * @param y_min Minimum Y
 * @param y_max Maximum Y
 * @param z_min Minimum Z
 * @param z_max Maximum Z
 * @return Box primitive data
 *
 * Note: Automatically swaps min/max if given in wrong order
 */
alea_primitive_data_t alea_make_box(
    double x_min, double x_max,
    double y_min, double y_max,
    double z_min, double z_max
);

// ============================================================================
// CYLINDERS (ALL THREE AXES)
// ============================================================================

/**
 * Create a cylinder along X-axis (infinite axial extent)
 *
 * @param center_y Center Y coordinate
 * @param center_z Center Z coordinate
 * @param radius Cylinder radius
 * @return Cylinder primitive data
 */
alea_primitive_data_t alea_make_cylinder_x(
    double center_y, double center_z, double radius
);

/**
 * Create a cylinder along Y-axis (infinite axial extent)
 */
alea_primitive_data_t alea_make_cylinder_y(
    double center_x, double center_z, double radius
);

/**
 * Create a cylinder along Z-axis (infinite axial extent)
 */
alea_primitive_data_t alea_make_cylinder_z(
    double center_x, double center_y, double radius
);

// ============================================================================
// CONES (ALL THREE AXES)
// ============================================================================

/**
 * Create a cone along X-axis (infinite axial extent)
 *
 * @param apex_x Apex X coordinate
 * @param apex_y Apex Y coordinate
 * @param apex_z Apex Z coordinate
 * @param tan_angle Tangent of half-angle
 * @return Cone primitive data
 *
 * Note: tan_angle is converted to tan²(angle) internally
 */
alea_primitive_data_t alea_make_cone_x(
    double apex_x, double apex_y, double apex_z,
    double tan_angle
);

/**
 * Create a cone along Y-axis (infinite axial extent)
 */
alea_primitive_data_t alea_make_cone_y(
    double apex_x, double apex_y, double apex_z,
    double tan_angle
);

/**
 * Create a cone along Z-axis (infinite axial extent)
 */
alea_primitive_data_t alea_make_cone_z(
    double apex_x, double apex_y, double apex_z,
    double tan_angle
);

// ============================================================================
// ADVANCED PRIMITIVES
// ============================================================================

/**
 * Create a general quadric surface
 * Ax² + By² + Cz² + Dxy + Eyz + Fxz + Gx + Hy + Iz + J = 0
 *
 * @param coeffs Array of 10 coefficients [A, B, C, D, E, F, G, H, I, J]
 * @return Quadric primitive data
 */
alea_primitive_data_t alea_make_quadric(const double coeffs[10]);

/**
 * Create a torus along X-axis
 *
 * @param center_x Center X coordinate
 * @param center_y Center Y coordinate
 * @param center_z Center Z coordinate
 * @param major_radius Major radius (center to tube center)
 * @param minor_radius Minor radius (tube radius)
 * @param axial_semiwidth_B Axial semiwidth parameter
 * @return Torus primitive data
 */
alea_primitive_data_t alea_make_torus_x(
    double center_x, double center_y, double center_z,
    double major_radius, double minor_radius,
    double axial_semiwidth_B
);

/**
 * Create a torus along Y-axis
 */
alea_primitive_data_t alea_make_torus_y(
    double center_x, double center_y, double center_z,
    double major_radius, double minor_radius,
    double axial_semiwidth_B
);

/**
 * Create a torus along Z-axis
 */
alea_primitive_data_t alea_make_torus_z(
    double center_x, double center_y, double center_z,
    double major_radius, double minor_radius,
    double axial_semiwidth_B
);

#ifdef __cplusplus
}
#endif

#endif // PRIMITIVE_CREATE_H
