// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file primitive_create.c
 * @brief Implementation of primitive factory functions
 * 
 * 
 * This module provides TWO levels of primitive creation:
 * 1. Data creation: alea_make_*() - Creates primitive data structures
 * 2. Node creation: alea_create_*_node() - Creates complete CSG nodes with pool storage
 */

#include "primitive_create.h"
#include "bbox.h"
#include "util/math.h"
#include <string.h>

// ============================================================================
// LEVEL 1: PRIMITIVE DATA CREATION
// ============================================================================
// These functions create primitive data structures (no nodes, no pools)
// These remain unchanged - they just create data

alea_primitive_data_t alea_make_plane(double a, double b, double c, double d) {
    alea_primitive_data_t data;
    data.plane.a = a;
    data.plane.b = b;
    data.plane.c = c;
    data.plane.d = d;
    return data;
}

alea_primitive_data_t alea_make_sphere(double cx, double cy, double cz, double radius) {
    alea_primitive_data_t data;
    data.sphere.center_x = cx;
    data.sphere.center_y = cy;
    data.sphere.center_z = cz;
    data.sphere.radius = radius;
    return data;
}

alea_primitive_data_t alea_make_box(
    double x_min, double x_max,
    double y_min, double y_max,
    double z_min, double z_max
) {
    alea_primitive_data_t data;
    
    // Ensure min <= max (swap if needed)
    if (x_min > x_max) { double t = x_min; x_min = x_max; x_max = t; }
    if (y_min > y_max) { double t = y_min; y_min = y_max; y_max = t; }
    if (z_min > z_max) { double t = z_min; z_min = z_max; z_max = t; }
    
    data.box.min_x = x_min;
    data.box.max_x = x_max;
    data.box.min_y = y_min;
    data.box.max_y = y_max;
    data.box.min_z = z_min;
    data.box.max_z = z_max;
    
    return data;
}

// ============================================================================
// CYLINDERS
// ============================================================================

alea_primitive_data_t alea_make_cylinder_x(
    double center_y, double center_z, double radius
) {
    alea_primitive_data_t data;
    memset(&data, 0, sizeof(data));
    data.cyl_x.center_y = center_y;
    data.cyl_x.center_z = center_z;
    data.cyl_x.radius = radius;
    return data;
}

alea_primitive_data_t alea_make_cylinder_y(
    double center_x, double center_z, double radius
) {
    alea_primitive_data_t data;
    memset(&data, 0, sizeof(data));
    data.cyl_y.center_x = center_x;
    data.cyl_y.center_z = center_z;
    data.cyl_y.radius = radius;
    return data;
}

alea_primitive_data_t alea_make_cylinder_z(
    double center_x, double center_y, double radius
) {
    alea_primitive_data_t data;
    memset(&data, 0, sizeof(data));
    data.cyl_z.center_x = center_x;
    data.cyl_z.center_y = center_y;
    data.cyl_z.radius = radius;
    return data;
}

// ============================================================================
// CONES
// ============================================================================

alea_primitive_data_t alea_make_cone_x(
    double apex_x, double apex_y, double apex_z,
    double tan_angle
) {
    alea_primitive_data_t data;
    memset(&data, 0, sizeof(data));
    data.cone_x.apex_x = apex_x;
    data.cone_x.apex_y = apex_y;
    data.cone_x.apex_z = apex_z;
    data.cone_x.tan_angle_sq = tan_angle * tan_angle;
    return data;
}

alea_primitive_data_t alea_make_cone_y(
    double apex_x, double apex_y, double apex_z,
    double tan_angle
) {
    alea_primitive_data_t data;
    memset(&data, 0, sizeof(data));
    data.cone_y.apex_x = apex_x;
    data.cone_y.apex_y = apex_y;
    data.cone_y.apex_z = apex_z;
    data.cone_y.tan_angle_sq = tan_angle * tan_angle;
    return data;
}

alea_primitive_data_t alea_make_cone_z(
    double apex_x, double apex_y, double apex_z,
    double tan_angle
) {
    alea_primitive_data_t data;
    memset(&data, 0, sizeof(data));
    data.cone_z.apex_x = apex_x;
    data.cone_z.apex_y = apex_y;
    data.cone_z.apex_z = apex_z;
    data.cone_z.tan_angle_sq = tan_angle * tan_angle;
    return data;
}

// ============================================================================
// ADVANCED PRIMITIVES
// ============================================================================

alea_primitive_data_t alea_make_quadric(const double coeffs[10]) {
    alea_primitive_data_t data;
    
    if (coeffs) {
        for (int i = 0; i < 10; i++) {
            data.quadric.coeffs[i] = coeffs[i];
        }
    } else {
        // Initialize to zero if NULL
        memset(data.quadric.coeffs, 0, sizeof(data.quadric.coeffs));
    }
    
    return data;
}

alea_primitive_data_t alea_make_torus_x(
    double center_x, double center_y, double center_z,
    double major_radius, double minor_radius,
    double axial_semiwidth_B
) {
    alea_primitive_data_t data;
    
    data.torus.axis = ALEA_AXIS_X;
    data.torus.center_x = center_x;    
    data.torus.center_y = center_y;
    data.torus.center_z = center_z;
    data.torus.major_radius = major_radius;
    data.torus.minor_radius = minor_radius;
    data.torus.axial_semiwidth_B = axial_semiwidth_B;
    
    return data;
}

alea_primitive_data_t alea_make_torus_y(
    double center_x, double center_y, double center_z,
    double major_radius, double minor_radius,
    double axial_semiwidth_B
) {
    alea_primitive_data_t data;
    
    data.torus.axis = ALEA_AXIS_Y;
    data.torus.center_x = center_x;    
    data.torus.center_y = center_y;
    data.torus.center_z = center_z;
    data.torus.major_radius = major_radius;
    data.torus.minor_radius = minor_radius;
    data.torus.axial_semiwidth_B = axial_semiwidth_B;
    
    return data;
}

alea_primitive_data_t alea_make_torus_z(
    double center_x, double center_y, double center_z,
    double major_radius, double minor_radius,
    double axial_semiwidth_B
) {
    alea_primitive_data_t data;
    
    data.torus.axis = ALEA_AXIS_Z;
    data.torus.center_x = center_x;    
    data.torus.center_y = center_y;
    data.torus.center_z = center_z;
    data.torus.major_radius = major_radius;
    data.torus.minor_radius = minor_radius;
    data.torus.axial_semiwidth_B = axial_semiwidth_B;
    
    return data;
}



