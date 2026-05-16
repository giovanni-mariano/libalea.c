// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef BBOX_H
#define BBOX_H

#include "core/alea_system.h"
#include <stdbool.h>

/* Forward decl to avoid pulling alea_universe.h transitively */
struct alea_matrix;


/**
 * @file bbox.h
 * @brief Bounding box operations for spatial queries
 * 
 * Provides efficient bounding box computation and intersection tests
 * for geometric primitives. Bounding boxes are used for:
 * 
 * 1. Early rejection - skip evaluation if point is outside bbox
 * 2. Spatial partitioning - organize geometry hierarchically
 * 3. Collision detection - quick overlap tests
 */

/**
 * Compute bounding box for a primitive
 * 
 * @param type Type of primitive
 * @param data Primitive data
 * @return Bounding box enclosing the primitive
 * 
 * Note: For infinite primitives (e.g., unbounded planes), returns
 * a conservatively large bounding box.
 * 
 * Example:
 * ```c
 * alea_primitive_data_t sphere;
 * sphere.sphere.center_x = 0.0;
 * sphere.sphere.center_y = 0.0;
 * sphere.sphere.center_z = 0.0;
 * sphere.sphere.radius = 5.0;
 * 
 * alea_bbox_t bbox = alea_primitive_bbox(ALEA_PRIMITIVE_SPHERE, &sphere);
 * // bbox = {-5, 5, -5, 5, -5, 5}
 * ```
 */
alea_bbox_t alea_primitive_bbox(
    alea_primitive_type_t type,
    const alea_primitive_data_t* data
);

/**
 * Recursively compute bounding box for a CSG tree node.
 *
 * @param sys  CSG system containing the tree
 * @param id   Node ID of the CSG tree root
 * @return Bounding box enclosing the CSG subtree
 */
alea_bbox_t alea_get_bbox(alea_system_t* sys,  
                        alea_node_id_t id);

/**
 * Compute union of two bounding boxes
 * Returns smallest box containing both input boxes
 * 
 * @param a First bounding box
 * @param b Second bounding box
 * @return Union (smallest enclosing box)
 */
alea_bbox_t alea_bbox_union(const alea_bbox_t* a, const alea_bbox_t* b);

/**
 * Compute intersection of two bounding boxes
 * Returns overlap region (empty if boxes don't intersect)
 * 
 * @param a First bounding box
 * @param b Second bounding box
 * @return Intersection (overlap region)
 */
alea_bbox_t alea_bbox_intersection(const alea_bbox_t* a, const alea_bbox_t* b);

/**
 * Test if two bounding boxes intersect
 * 
 * @param a First bounding box
 * @param b Second bounding box
 * @return true if boxes overlap
 * 
 * Fast separating axis test - O(6) comparisons
 */
bool alea_bbox_intersects(const alea_bbox_t* a, const alea_bbox_t* b);

/**
 * Test if a bounding box contains a point
 * 
 * @param bbox Bounding box
 * @param x X coordinate
 * @param y Y coordinate
 * @param z Z coordinate
 * @return true if point is inside or on boundary
 * 
 * Useful for early rejection in ray tracing and point queries
 */
bool alea_bbox_contains_point(const alea_bbox_t* bbox, double x, double y, double z);

/**
 * Check if one bounding box entirely contains another
 *
 * @param outer The potentially containing bbox
 * @param inner The potentially contained bbox
 * @return true if inner is entirely within outer (or equal)
 */
bool alea_bbox_contains(const alea_bbox_t* outer, const alea_bbox_t* inner);

/**
 * Get center point of bounding box
 * 
 * @param bbox Bounding box
 * @param out_x Output X coordinate
 * @param out_y Output Y coordinate
 * @param out_z Output Z coordinate
 */
void alea_bbox_center(const alea_bbox_t* bbox, double* out_x, double* out_y, double* out_z);

/**
 * Get dimensions (width, height, depth) of bounding box
 * 
 * @param bbox Bounding box
 * @param out_width Output width (x dimension)
 * @param out_height Output height (y dimension)
 * @param out_depth Output depth (z dimension)
 */
void alea_bbox_dimensions(const alea_bbox_t* bbox, double* out_width, double* out_height, double* out_depth);

/**
 * Get volume of bounding box
 * 
 * @param bbox Bounding box
 * @return Volume (width × height × depth)
 */
double alea_bbox_volume(const alea_bbox_t* bbox);

/**
 * Get surface area of bounding box
 * 
 * @param bbox Bounding box
 * @return Surface area
 */
double alea_bbox_surface_area(const alea_bbox_t* bbox);

/**
 * Check if bounding box is valid (min <= max for all dimensions)
 * 
 * @param bbox Bounding box to validate
 * @return true if valid
 */
bool alea_bbox_is_valid(const alea_bbox_t* bbox);

/**
 * Create an empty bounding box (min > max)
 * Useful as initial value for computing unions
 *
 * @return Empty bounding box
 */
alea_bbox_t alea_bbox_empty(void);

/**
 * Transform an axis-aligned bounding box by a matrix and re-AABB the result.
 *
 * Conservative under rotation: returns the smallest AABB enclosing the
 * 8 transformed corners (computed in O(1) from min/max contributions, not
 * by enumerating corners). Safe for cull tests — false positives only.
 *
 * @param bbox Bounding box in source frame
 * @param mat  Transform matrix (row-major 3x4)
 * @return Transformed bounding box
 */
alea_bbox_t alea_bbox_transform(const alea_bbox_t* bbox,
                                const struct alea_matrix* mat);

/**
 * Create an infinite bounding box (covers all space)
 * Useful for unbounded primitives
 *
 * @return Infinite bounding box
 */
alea_bbox_t alea_bbox_infinite(void);

/**
 * Compute bounding box for a half-space (primitive + sense)
 *
 * Unlike alea_primitive_bbox which returns the bbox of the surface,
 * this function returns the bbox of the half-space selected by sense.
 *
 * For bounded primitives (sphere, box, rcc):
 * - Negative sense (inside): primitive's bbox
 * - Positive sense (outside): infinite bbox
 *
 * For axis-aligned planes:
 * - Positive sense: half-space on positive side of plane
 * - Negative sense: half-space on negative side of plane
 *
 * @param type Primitive type
 * @param data Primitive data
 * @param sense +1 for positive half-space, -1 for negative
 * @return Bounding box of the selected half-space
 */
alea_bbox_t alea_halfspace_bbox(
    alea_primitive_type_t type,
    const alea_primitive_data_t* data,
    int8_t sense
);

/**
 * Tighten a cell's bounding box using interval arithmetic.
 *
 * For each axis, binary-searches from both sides to discard regions
 * that are entirely outside (POSITIVE relation). Converges in
 * O(log(extent/tol)) steps per axis.
 *
 * @param sys   CSG system
 * @param root  Root node of the cell's CSG tree
 * @param start Initial (conservative) bounding box
 * @param tol   Spatial tolerance for convergence
 * @param out   Output tightened bounding box
 */
void alea_tighten_tree_bbox(const alea_system_t* sys,
                            alea_node_id_t root,
                            const alea_bbox_t* start,
                            double tol,
                            alea_bbox_t* out);

/**
 * Compute a tight bounding sphere for the entire model.
 *
 * Tightens each bounded cell's bbox via interval arithmetic, unions them,
 * and returns the enclosing sphere.
 *
 * @param sys      CSG system
 * @param tol      Spatial tolerance for tightening (e.g. 1.0)
 * @param cx,cy,cz Output center
 * @param radius   Output radius
 * @return 0 on success, -1 on error
 */
int alea_compute_bounding_sphere(alea_system_t* sys,
                                 double tol,
                                 double* cx, double* cy, double* cz,
                                 double* radius);

/**
 * Tighten all cell bboxes in the system in place.
 *
 * Iterates over every bounded cell, applies alea_tighten_tree_bbox(),
 * and writes the result back to the root node's bbox.
 *
 * @param sys  CSG system (modified: node bboxes updated)
 * @param tol  Spatial tolerance for tightening
 * @return Number of cells tightened, or -1 on error
 */
int alea_tighten_all_bboxes(alea_system_t* sys, double tol);

/**
 * Numerically tighten a cell's bounding box.
 *
 * Uses two complementary strategies:
 * 1. LP vertex enumeration — for pure plane intersections, enumerate
 *    triple-plane intersection vertices via Cramer's rule.
 * 2. Octree recursive bisection — for general cells, recursively
 *    subdivide a large search box and discard empty sub-boxes.
 *
 * After either method produces a finite bbox, a final binary-search
 * pass (alea_tighten_tree_bbox) refines to the given tolerance.
 *
 * @param sys   CSG system
 * @param root  Root node of the cell's CSG tree
 * @param tol   Spatial tolerance for final binary-search refinement
 * @param out   Output tightened bounding box
 * @return 0 on success (finite bbox found), -1 on failure
 */
int alea_tighten_bbox_numerical(const alea_system_t* sys,
                               alea_node_id_t root,
                               double tol,
                               alea_bbox_t* out);


#endif // BBOX_H