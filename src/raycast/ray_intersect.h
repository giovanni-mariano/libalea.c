// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef RAY_INTERSECT_H
#define RAY_INTERSECT_H

#include "raycast.h"
#include "alea_types.h"

/**
 * @file ray_intersect.h
 * @brief Ray-primitive intersection routines (internal)
 *
 * Each function returns the number of intersections found (0, 1, or 2)
 * and fills t_out array with distances. Caller must provide t_out[2].
 *
 * For infinite surfaces (planes, infinite cylinders), returns 1 intersection.
 * For closed surfaces (spheres, bounded cylinders), returns 0 or 2.
 */

/* Plane: ax + by + cz + d = 0 */
int ray_intersect_plane(const alea_ray_t* ray,
                        const alea_plane_data_t* plane,
                        double* restrict t_out,
                        double* restrict nx, double* restrict ny, double* restrict nz);

/* Sphere */
int ray_intersect_sphere(const alea_ray_t* ray,
                         const alea_sphere_data_t* sphere,
                         double* restrict t_out);

/* Axis-aligned cylinders (infinite) */
int ray_intersect_cylinder_x(const alea_ray_t* ray,
                             const alea_cylinder_x_data_t* cyl,
                             double* restrict t_out);

int ray_intersect_cylinder_y(const alea_ray_t* ray,
                             const alea_cylinder_y_data_t* cyl,
                             double* restrict t_out);

int ray_intersect_cylinder_z(const alea_ray_t* ray,
                             const alea_cylinder_z_data_t* cyl,
                             double* restrict t_out);

/* Axis-aligned cones (infinite, both sheets) */
int ray_intersect_cone_x(const alea_ray_t* ray,
                         const alea_cone_x_data_t* cone,
                         double* restrict t_out);

int ray_intersect_cone_y(const alea_ray_t* ray,
                         const alea_cone_y_data_t* cone,
                         double* restrict t_out);

int ray_intersect_cone_z(const alea_ray_t* ray,
                         const alea_cone_z_data_t* cone,
                         double* restrict t_out);

/* Axis-aligned box */
int ray_intersect_box(const alea_ray_t* ray,
                      const alea_box_data_t* box,
                      double* restrict t_out);

/* General quadric */
int ray_intersect_quadric(const alea_ray_t* ray,
                          const alea_quadric_data_t* quadric,
                          double* restrict t_out);

/* Torus (expensive!) */
int ray_intersect_torus(const alea_ray_t* ray,
                        const alea_torus_data_t* torus,
                        double* restrict t_out);

/* RCC macrobody (finite cylinder, arbitrary orientation) */
int ray_intersect_rcc(const alea_ray_t* ray,
                      const alea_rcc_data_t* rcc,
                      double* restrict t_out);

/* TRC macrobody (truncated cone) */
int ray_intersect_trc(const alea_ray_t* ray,
                      const alea_trc_data_t* trc,
                      double* restrict t_out);

/**
 * @brief Dispatch to appropriate intersection function by primitive type
 *
 * @param ray Ray to test
 * @param type Primitive type
 * @param data Primitive data
 * @param t_out Output array for distances (size >= 2)
 * @return Number of intersections found
 */
int ray_intersect_primitive(const alea_ray_t* ray,
                            alea_primitive_type_t type,
                            const alea_primitive_data_t* data,
                            double* restrict t_out);

/**
 * @brief Compute surface normal at intersection point
 *
 * @param type Primitive type
 * @param data Primitive data
 * @param px, py, pz Point on surface
 * @param nx, ny, nz Output normal (normalized, pointing outward)
 */
void primitive_normal_at(alea_primitive_type_t type,
                         const alea_primitive_data_t* data,
                         double px, double py, double pz,
                         double* nx, double* ny, double* nz);

#endif /* RAY_INTERSECT_H */
