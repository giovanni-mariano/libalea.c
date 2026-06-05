// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_simplify.c
 * @brief CSG tree simplification implementation
 * 
 * 
 * Simplifications performed:
 * 1. NNF conversion (push complements to leaves by flipping sense)
 * 2. Contradiction detection (A ∩ ¬A = empty)
 * 3. Tautology detection (A ∪ ¬A = universe)  
 * 4. Idempotent reduction (A ∩ A = A, A ∪ A = A)
 * 5. Statistical simplification (remove always-true/always-false terms)
 */

#include "alea_simplify.h"
#include "core/alea_system.h"
#include "alea_ops.h"
#include "alea_cell_complement.h"
#include "primitives/bbox.h"
#include "alea_eval.h"
#include "primitives/primitive_eval.h"
#include "primitives/bbox.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "util/alea_log.h"
#include "util/alea_bitset.h"

// ============================================================================
// INTERVAL ARITHMETIC ENGINE
// ============================================================================

typedef struct {
    double lo;
    double hi;
} interval_t;

static inline interval_t iv_make(double lo, double hi) {
    return (interval_t){lo, hi};
}

static inline interval_t iv_add(interval_t a, interval_t b) {
    return (interval_t){a.lo + b.lo, a.hi + b.hi};
}

static inline interval_t iv_sub(interval_t a, interval_t b) {
    return (interval_t){a.lo - b.hi, a.hi - b.lo};
}

static inline interval_t iv_mul_scalar(interval_t a, double s) {
    if (s >= 0) return (interval_t){a.lo * s, a.hi * s};
    return (interval_t){a.hi * s, a.lo * s};
}

// Full interval multiplication (for cross terms like x*y)
static inline interval_t iv_mul(interval_t a, interval_t b) {
    double p1 = a.lo * b.lo;
    double p2 = a.lo * b.hi;
    double p3 = a.hi * b.lo;
    double p4 = a.hi * b.hi;

    double lo = p1;
    if (p2 < lo) lo = p2;
    if (p3 < lo) lo = p3;
    if (p4 < lo) lo = p4;

    double hi = p1;
    if (p2 > hi) hi = p2;
    if (p3 > hi) hi = p3;
    if (p4 > hi) hi = p4;

    return (interval_t){lo, hi};
}

// Interval square - handles zero-crossing correctly
static inline interval_t iv_sqr(interval_t a) {
    double s1 = a.lo * a.lo;
    double s2 = a.hi * a.hi;
    // If range crosses zero, minimum squared is 0
    if (a.lo <= 0 && a.hi >= 0) {
        return (interval_t){0.0, (s1 > s2) ? s1 : s2};
    }
    return (interval_t){(s1 < s2) ? s1 : s2, (s1 > s2) ? s1 : s2};
}

// Result of checking box against surface
typedef enum {
    BOX_STRICTLY_NEGATIVE,  // f(x) < 0 everywhere in box (inside surface)
    BOX_STRICTLY_POSITIVE,  // f(x) > 0 everywhere in box (outside surface)
    BOX_INTERSECTS          // f(x) changes sign in box
} box_surface_relation_t;

// Tolerance for interval comparisons
#define IV_TOL 1e-9

// ============================================================================
// SURFACE-SPECIFIC INTERVAL SOLVERS
// ============================================================================

/**
 * Check plane: ax + by + cz + d = 0
 * Positive half-space: ax + by + cz + d > 0
 */
static box_surface_relation_t check_plane(
    const alea_plane_data_t* plane,
    const alea_bbox_t* box
) {
    interval_t x = iv_make(box->min_x, box->max_x);
    interval_t y = iv_make(box->min_y, box->max_y);
    interval_t z = iv_make(box->min_z, box->max_z);

    interval_t result = iv_make(plane->d, plane->d);
    result = iv_add(result, iv_mul_scalar(x, plane->a));
    result = iv_add(result, iv_mul_scalar(y, plane->b));
    result = iv_add(result, iv_mul_scalar(z, plane->c));

    if (result.hi < -IV_TOL) return BOX_STRICTLY_NEGATIVE;
    if (result.lo > IV_TOL)  return BOX_STRICTLY_POSITIVE;
    return BOX_INTERSECTS;
}

/**
 * Check sphere: (x-cx)² + (y-cy)² + (z-cz)² - r² = 0
 * Inside sphere: d² < r² (negative)
 * Outside sphere: d² > r² (positive)
 */
static box_surface_relation_t check_sphere(
    const alea_sphere_data_t* sphere,
    const alea_bbox_t* box
) {
    interval_t dx = iv_sub(iv_make(box->min_x, box->max_x),
                           iv_make(sphere->center_x, sphere->center_x));
    interval_t dy = iv_sub(iv_make(box->min_y, box->max_y),
                           iv_make(sphere->center_y, sphere->center_y));
    interval_t dz = iv_sub(iv_make(box->min_z, box->max_z),
                           iv_make(sphere->center_z, sphere->center_z));

    interval_t d2 = iv_add(iv_add(iv_sqr(dx), iv_sqr(dy)), iv_sqr(dz));
    double r2 = sphere->radius * sphere->radius;

    // d² - r²: negative means inside, positive means outside
    if (d2.hi < r2 - IV_TOL) return BOX_STRICTLY_NEGATIVE;
    if (d2.lo > r2 + IV_TOL) return BOX_STRICTLY_POSITIVE;
    return BOX_INTERSECTS;
}

/**
 * Check cylinder along Z-axis: (x-cx)² + (y-cy)² - r² = 0
 */
static box_surface_relation_t check_cylinder_z(
    const alea_cylinder_z_data_t* cyl,
    const alea_bbox_t* box
) {
    interval_t dx = iv_sub(iv_make(box->min_x, box->max_x),
                           iv_make(cyl->center_x, cyl->center_x));
    interval_t dy = iv_sub(iv_make(box->min_y, box->max_y),
                           iv_make(cyl->center_y, cyl->center_y));

    interval_t d2 = iv_add(iv_sqr(dx), iv_sqr(dy));
    double r2 = cyl->radius * cyl->radius;

    if (d2.hi < r2 - IV_TOL) return BOX_STRICTLY_NEGATIVE;
    if (d2.lo > r2 + IV_TOL) return BOX_STRICTLY_POSITIVE;
    return BOX_INTERSECTS;
}

/**
 * Check cylinder along X-axis: (y-cy)² + (z-cz)² - r² = 0
 */
static box_surface_relation_t check_cylinder_x(
    const alea_cylinder_x_data_t* cyl,
    const alea_bbox_t* box
) {
    interval_t dy = iv_sub(iv_make(box->min_y, box->max_y),
                           iv_make(cyl->center_y, cyl->center_y));
    interval_t dz = iv_sub(iv_make(box->min_z, box->max_z),
                           iv_make(cyl->center_z, cyl->center_z));

    interval_t d2 = iv_add(iv_sqr(dy), iv_sqr(dz));
    double r2 = cyl->radius * cyl->radius;

    if (d2.hi < r2 - IV_TOL) return BOX_STRICTLY_NEGATIVE;
    if (d2.lo > r2 + IV_TOL) return BOX_STRICTLY_POSITIVE;
    return BOX_INTERSECTS;
}

/**
 * Check cylinder along Y-axis: (x-cx)² + (z-cz)² - r² = 0
 */
static box_surface_relation_t check_cylinder_y(
    const alea_cylinder_y_data_t* cyl,
    const alea_bbox_t* box
) {
    interval_t dx = iv_sub(iv_make(box->min_x, box->max_x),
                           iv_make(cyl->center_x, cyl->center_x));
    interval_t dz = iv_sub(iv_make(box->min_z, box->max_z),
                           iv_make(cyl->center_z, cyl->center_z));

    interval_t d2 = iv_add(iv_sqr(dx), iv_sqr(dz));
    double r2 = cyl->radius * cyl->radius;

    if (d2.hi < r2 - IV_TOL) return BOX_STRICTLY_NEGATIVE;
    if (d2.lo > r2 + IV_TOL) return BOX_STRICTLY_POSITIVE;
    return BOX_INTERSECTS;
}

/**
 * Check cone along Z-axis: (x-ax)² + (y-ay)² - t²(z-az)² = 0
 * where t² = tan²(half-angle)
 */
static box_surface_relation_t check_cone_z(
    const alea_cone_z_data_t* cone,
    const alea_bbox_t* box
) {
    interval_t dx = iv_sub(iv_make(box->min_x, box->max_x),
                           iv_make(cone->apex_x, cone->apex_x));
    interval_t dy = iv_sub(iv_make(box->min_y, box->max_y),
                           iv_make(cone->apex_y, cone->apex_y));
    interval_t dz = iv_sub(iv_make(box->min_z, box->max_z),
                           iv_make(cone->apex_z, cone->apex_z));

    interval_t r2 = iv_add(iv_sqr(dx), iv_sqr(dy));
    interval_t h2 = iv_mul_scalar(iv_sqr(dz), cone->tan_angle_sq);

    // r² - t²z²: negative means inside cone, positive means outside
    interval_t result = iv_sub(r2, h2);

    if (result.hi < -IV_TOL) return BOX_STRICTLY_NEGATIVE;
    if (result.lo > IV_TOL)  return BOX_STRICTLY_POSITIVE;
    return BOX_INTERSECTS;
}

/**
 * Check cone along X-axis: (y-ay)² + (z-az)² - t²(x-ax)² = 0
 */
static box_surface_relation_t check_cone_x(
    const alea_cone_x_data_t* cone,
    const alea_bbox_t* box
) {
    interval_t dx = iv_sub(iv_make(box->min_x, box->max_x),
                           iv_make(cone->apex_x, cone->apex_x));
    interval_t dy = iv_sub(iv_make(box->min_y, box->max_y),
                           iv_make(cone->apex_y, cone->apex_y));
    interval_t dz = iv_sub(iv_make(box->min_z, box->max_z),
                           iv_make(cone->apex_z, cone->apex_z));

    interval_t r2 = iv_add(iv_sqr(dy), iv_sqr(dz));
    interval_t h2 = iv_mul_scalar(iv_sqr(dx), cone->tan_angle_sq);

    interval_t result = iv_sub(r2, h2);

    if (result.hi < -IV_TOL) return BOX_STRICTLY_NEGATIVE;
    if (result.lo > IV_TOL)  return BOX_STRICTLY_POSITIVE;
    return BOX_INTERSECTS;
}

/**
 * Check cone along Y-axis: (x-ax)² + (z-az)² - t²(y-ay)² = 0
 */
static box_surface_relation_t check_cone_y(
    const alea_cone_y_data_t* cone,
    const alea_bbox_t* box
) {
    interval_t dx = iv_sub(iv_make(box->min_x, box->max_x),
                           iv_make(cone->apex_x, cone->apex_x));
    interval_t dy = iv_sub(iv_make(box->min_y, box->max_y),
                           iv_make(cone->apex_y, cone->apex_y));
    interval_t dz = iv_sub(iv_make(box->min_z, box->max_z),
                           iv_make(cone->apex_z, cone->apex_z));

    interval_t r2 = iv_add(iv_sqr(dx), iv_sqr(dz));
    interval_t h2 = iv_mul_scalar(iv_sqr(dy), cone->tan_angle_sq);

    interval_t result = iv_sub(r2, h2);

    if (result.hi < -IV_TOL) return BOX_STRICTLY_NEGATIVE;
    if (result.lo > IV_TOL)  return BOX_STRICTLY_POSITIVE;
    return BOX_INTERSECTS;
}

/**
 * Check general quadric: Ax² + By² + Cz² + Dxy + Eyz + Fzx + Gx + Hy + Iz + J = 0
 * coeffs[10] = {A, B, C, D, E, F, G, H, I, J}
 */
static box_surface_relation_t check_quadric(
    const alea_quadric_data_t* gq,
    const alea_bbox_t* box
) {
    interval_t x = iv_make(box->min_x, box->max_x);
    interval_t y = iv_make(box->min_y, box->max_y);
    interval_t z = iv_make(box->min_z, box->max_z);

    // Start with constant term J
    interval_t sum = iv_make(gq->coeffs[9], gq->coeffs[9]);

    // Squared terms: Ax², By², Cz²
    sum = iv_add(sum, iv_mul_scalar(iv_sqr(x), gq->coeffs[0]));
    sum = iv_add(sum, iv_mul_scalar(iv_sqr(y), gq->coeffs[1]));
    sum = iv_add(sum, iv_mul_scalar(iv_sqr(z), gq->coeffs[2]));

    // Linear terms: Gx, Hy, Iz
    sum = iv_add(sum, iv_mul_scalar(x, gq->coeffs[6]));
    sum = iv_add(sum, iv_mul_scalar(y, gq->coeffs[7]));
    sum = iv_add(sum, iv_mul_scalar(z, gq->coeffs[8]));

    // Cross terms: Dxy, Eyz, Fzx
    if (fabs(gq->coeffs[3]) > IV_TOL)
        sum = iv_add(sum, iv_mul_scalar(iv_mul(x, y), gq->coeffs[3]));
    if (fabs(gq->coeffs[4]) > IV_TOL)
        sum = iv_add(sum, iv_mul_scalar(iv_mul(y, z), gq->coeffs[4]));
    if (fabs(gq->coeffs[5]) > IV_TOL)
        sum = iv_add(sum, iv_mul_scalar(iv_mul(z, x), gq->coeffs[5]));

    if (sum.hi < -IV_TOL) return BOX_STRICTLY_NEGATIVE;
    if (sum.lo > IV_TOL)  return BOX_STRICTLY_POSITIVE;
    return BOX_INTERSECTS;
}

/**
 * Check torus (Z-axis): distance to major circle skeleton
 *
 * For a torus centered at (cx, cy, cz) with axis along Z:
 * The skeleton is a circle of radius R in the XY plane at z=cz.
 * A point is inside if its distance to the skeleton is less than r.
 */
static box_surface_relation_t check_torus_z(
    const alea_torus_data_t* torus,
    const alea_bbox_t* box
) {
    interval_t dx = iv_sub(iv_make(box->min_x, box->max_x),
                           iv_make(torus->center_x, torus->center_x));
    interval_t dy = iv_sub(iv_make(box->min_y, box->max_y),
                           iv_make(torus->center_y, torus->center_y));
    interval_t dz = iv_sub(iv_make(box->min_z, box->max_z),
                           iv_make(torus->center_z, torus->center_z));

    // Distance squared from axis in XY plane
    interval_t rho_sq = iv_add(iv_sqr(dx), iv_sqr(dy));

    // Distance interval from axis
    double rho_lo = sqrt(rho_sq.lo);
    double rho_hi = sqrt(rho_sq.hi);

    // Distance to major circle in radial direction
    double dR_lo, dR_hi;
    double R = torus->major_radius;

    if (rho_lo > R) {
        dR_lo = rho_lo - R;
        dR_hi = rho_hi - R;
    } else if (rho_hi < R) {
        dR_lo = R - rho_hi;
        dR_hi = R - rho_lo;
    } else {
        // Range crosses R - minimum distance is 0
        dR_lo = 0.0;
        dR_hi = fmax(rho_hi - R, R - rho_lo);
    }

    // Total distance squared to skeleton circle
    interval_t h_sq = iv_sqr(dz);
    double tot_lo = dR_lo * dR_lo + h_sq.lo;
    double tot_hi = dR_hi * dR_hi + h_sq.hi;

    double r2 = torus->minor_radius * torus->minor_radius;

    // distance² - r²: negative means inside tube, positive means outside
    if (tot_hi < r2 - IV_TOL) return BOX_STRICTLY_NEGATIVE;
    if (tot_lo > r2 + IV_TOL) return BOX_STRICTLY_POSITIVE;
    return BOX_INTERSECTS;
}

/**
 * Check torus (X-axis): similar to Z but rotated
 */
static box_surface_relation_t check_torus_x(
    const alea_torus_data_t* torus,
    const alea_bbox_t* box
) {
    interval_t dx = iv_sub(iv_make(box->min_x, box->max_x),
                           iv_make(torus->center_x, torus->center_x));
    interval_t dy = iv_sub(iv_make(box->min_y, box->max_y),
                           iv_make(torus->center_y, torus->center_y));
    interval_t dz = iv_sub(iv_make(box->min_z, box->max_z),
                           iv_make(torus->center_z, torus->center_z));

    // For X-axis torus, skeleton is in YZ plane
    interval_t rho_sq = iv_add(iv_sqr(dy), iv_sqr(dz));
    double rho_lo = sqrt(rho_sq.lo);
    double rho_hi = sqrt(rho_sq.hi);

    double R = torus->major_radius;
    double dR_lo, dR_hi;

    if (rho_lo > R) {
        dR_lo = rho_lo - R;
        dR_hi = rho_hi - R;
    } else if (rho_hi < R) {
        dR_lo = R - rho_hi;
        dR_hi = R - rho_lo;
    } else {
        dR_lo = 0.0;
        dR_hi = fmax(rho_hi - R, R - rho_lo);
    }

    interval_t h_sq = iv_sqr(dx);
    double tot_lo = dR_lo * dR_lo + h_sq.lo;
    double tot_hi = dR_hi * dR_hi + h_sq.hi;

    double r2 = torus->minor_radius * torus->minor_radius;

    if (tot_hi < r2 - IV_TOL) return BOX_STRICTLY_NEGATIVE;
    if (tot_lo > r2 + IV_TOL) return BOX_STRICTLY_POSITIVE;
    return BOX_INTERSECTS;
}

/**
 * Check torus (Y-axis)
 */
static box_surface_relation_t check_torus_y(
    const alea_torus_data_t* torus,
    const alea_bbox_t* box
) {
    interval_t dx = iv_sub(iv_make(box->min_x, box->max_x),
                           iv_make(torus->center_x, torus->center_x));
    interval_t dy = iv_sub(iv_make(box->min_y, box->max_y),
                           iv_make(torus->center_y, torus->center_y));
    interval_t dz = iv_sub(iv_make(box->min_z, box->max_z),
                           iv_make(torus->center_z, torus->center_z));

    // For Y-axis torus, skeleton is in XZ plane
    interval_t rho_sq = iv_add(iv_sqr(dx), iv_sqr(dz));
    double rho_lo = sqrt(rho_sq.lo);
    double rho_hi = sqrt(rho_sq.hi);

    double R = torus->major_radius;
    double dR_lo, dR_hi;

    if (rho_lo > R) {
        dR_lo = rho_lo - R;
        dR_hi = rho_hi - R;
    } else if (rho_hi < R) {
        dR_lo = R - rho_hi;
        dR_hi = R - rho_lo;
    } else {
        dR_lo = 0.0;
        dR_hi = fmax(rho_hi - R, R - rho_lo);
    }

    interval_t h_sq = iv_sqr(dy);
    double tot_lo = dR_lo * dR_lo + h_sq.lo;
    double tot_hi = dR_hi * dR_hi + h_sq.hi;

    double r2 = torus->minor_radius * torus->minor_radius;

    if (tot_hi < r2 - IV_TOL) return BOX_STRICTLY_NEGATIVE;
    if (tot_lo > r2 + IV_TOL) return BOX_STRICTLY_POSITIVE;
    return BOX_INTERSECTS;
}

// ============================================================================
// DISPATCHER: CHECK BOX-SURFACE RELATION
// ============================================================================

/**
 * Dispatch to appropriate checker based on primitive type
 */
static box_surface_relation_t check_box_surface_relation(
    const alea_system_t* sys,
    alea_node_id_t prim_node_id,
    const alea_bbox_t* box
) {
    const alea_node_t* node = &sys->nodes.data[prim_node_id];
    alea_primitive_id_t prim_id = node->primitive.primitive_id;
    const alea_primitive_entry_t* prim = &sys->primitives.data[prim_id];

    switch (prim->type) {
        case ALEA_PRIMITIVE_PLANE:
            return check_plane(&prim->data.plane, box);

        case ALEA_PRIMITIVE_SPHERE:
            return check_sphere(&prim->data.sphere, box);

        case ALEA_PRIMITIVE_CYLINDER_X:
            return check_cylinder_x(&prim->data.cyl_x, box);

        case ALEA_PRIMITIVE_CYLINDER_Y:
            return check_cylinder_y(&prim->data.cyl_y, box);

        case ALEA_PRIMITIVE_CYLINDER_Z:
            return check_cylinder_z(&prim->data.cyl_z, box);

        case ALEA_PRIMITIVE_CONE_X:
            return check_cone_x(&prim->data.cone_x, box);

        case ALEA_PRIMITIVE_CONE_Y:
            return check_cone_y(&prim->data.cone_y, box);

        case ALEA_PRIMITIVE_CONE_Z:
            return check_cone_z(&prim->data.cone_z, box);

        case ALEA_PRIMITIVE_QUADRIC:
            return check_quadric(&prim->data.quadric, box);

        case ALEA_PRIMITIVE_TORUS_X:
            return check_torus_x(&prim->data.torus, box);

        case ALEA_PRIMITIVE_TORUS_Y:
            return check_torus_y(&prim->data.torus, box);

        case ALEA_PRIMITIVE_TORUS_Z:
            return check_torus_z(&prim->data.torus, box);

        case ALEA_PRIMITIVE_RPP:
        case ALEA_PRIMITIVE_RCC:
            // Box and RCC are compound - conservative: assume intersection
            return BOX_INTERSECTS;

        default:
            return BOX_INTERSECTS;
    }
}

// ============================================================================
// HELPER: GET OPPOSITE SENSE NODE
// ============================================================================

/**
 * Given a primitive node, find the node with opposite sense.
 * Uses the surface_entry which already has both pos_node and neg_node.
 * Matches by sense value, not by node ID, to handle flattened systems
 * where multiple nodes may reference the same primitive.
 */
static alea_node_id_t get_opposite_sense_node(
    const alea_system_t* sys,
    alea_node_id_t prim_node_id
) {
    if (prim_node_id == ALEA_NODE_ID_INVALID || prim_node_id >= alea_vec_count(&sys->nodes)) {
        return ALEA_NODE_ID_INVALID;
    }

    const alea_node_t* node = &sys->nodes.data[prim_node_id];
    if (ALEA_GET_OPERATION(node) != ALEA_OP_PRIMITIVE) {
        return ALEA_NODE_ID_INVALID;
    }

    int mcnp_id = node->primitive.mc_surface_id;
    int8_t sense = node->primitive.sense;

    // Find surface entry by MCNP ID
    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        if (sys->surfaces.data[i].mc_surface_id == mcnp_id) {
            // Return the opposite sense node based on current node's sense
            if (sense > 0) {
                return sys->surfaces.data[i].neg_node;
            } else {
                return sys->surfaces.data[i].pos_node;
            }
        }
    }

    return ALEA_NODE_ID_INVALID;
}

// ============================================================================
// NNF CONVERSION (push complements to primitives)
// ============================================================================

/**
 * Convert tree to Negation Normal Form.
 * 
 * The key insight: we DON'T create new primitive nodes. When we need to
 * negate a primitive, we look up the surface_entry and return the 
 * pre-existing node with opposite sense.
 * 
 * @param sys CSG system
 * @param node_id Current node
 * @param negated Whether we're under an odd number of complements
 * @param stats Statistics (can be NULL)
 * @return Node ID of NNF-converted tree
 */
static alea_node_id_t to_nnf(
    alea_system_t* sys,
    alea_node_id_t node_id,
    bool negated,
    alea_simplify_stats_t* stats
) {
    if (node_id == ALEA_NODE_ID_INVALID || node_id >= alea_vec_count(&sys->nodes)) {
        return ALEA_NODE_ID_INVALID;
    }
    
    const alea_node_t* node = &sys->nodes.data[node_id];
    alea_operation_t op = ALEA_GET_OPERATION(node);
    
    switch (op) {
        case ALEA_OP_PRIMITIVE: {
            if (negated) {
                // Need to flip sense - get the opposite sense node
                alea_node_id_t opposite = get_opposite_sense_node(sys, node_id);
                if (opposite != ALEA_NODE_ID_INVALID) {
                    if (stats) stats->complements_eliminated++;
                    return opposite;
                }
                // Fallback: shouldn't happen if surfaces are set up correctly
                ALEA_LOG_WARN("couldn't find opposite sense for surface %d",
                        node->primitive.mc_surface_id);
                return node_id;
            }
            return node_id;
        }
        
        case ALEA_OP_COMPLEMENT: {
            alea_node_id_t child_id = node->operation.left;
            
            // Check for unresolved cell complement
            if (child_id == ALEA_NODE_ID_INVALID) {
                return node_id;  // Can't simplify unresolved reference
            }
            
            if (stats) stats->complements_eliminated++;
            if (negated && stats) stats->double_negations++;
            
            // Complement flips the negation
            return to_nnf(sys, child_id, !negated, stats);
        }
        
        case ALEA_OP_UNION:
        case ALEA_OP_INTERSECTION: {
            // De Morgan: ¬(A ∪ B) = ¬A ∩ ¬B, ¬(A ∩ B) = ¬A ∪ ¬B
            alea_operation_t result_op = op;
            if (negated) {
                result_op = (op == ALEA_OP_UNION) ? ALEA_OP_INTERSECTION : ALEA_OP_UNION;
            }
            
            // Save child IDs before recursing (reallocation safety)
            alea_node_id_t left_id = node->operation.left;
            alea_node_id_t right_id = node->operation.right;
            
            alea_node_id_t new_left = to_nnf(sys, left_id, negated, stats);
            alea_node_id_t new_right = to_nnf(sys, right_id, negated, stats);
            
            if (new_left == ALEA_NODE_ID_INVALID || new_right == ALEA_NODE_ID_INVALID) {
                return ALEA_NODE_ID_INVALID;
            }
            
            // If nothing changed and not negated, return original
            if (!negated && new_left == left_id && new_right == right_id) {
                return node_id;
            }
            
            // Create new operation node
            if (result_op == ALEA_OP_UNION) {
                return alea_create_union(sys, new_left, new_right);
            } else {
                return alea_create_intersection(sys, new_left, new_right);
            }
        }
        
        case ALEA_OP_DIFFERENCE: {
            // A - B = A ∩ ¬B
            // ¬(A - B) = ¬A ∪ B
            alea_node_id_t left_id = node->operation.left;
            alea_node_id_t right_id = node->operation.right;
            
            alea_node_id_t new_left, new_right;
            alea_operation_t result_op;
            
            if (negated) {
                // ¬(A - B) = ¬A ∪ B
                new_left = to_nnf(sys, left_id, true, stats);
                new_right = to_nnf(sys, right_id, false, stats);
                result_op = ALEA_OP_UNION;
            } else {
                // A - B = A ∩ ¬B
                new_left = to_nnf(sys, left_id, false, stats);
                new_right = to_nnf(sys, right_id, true, stats);
                result_op = ALEA_OP_INTERSECTION;
            }
            
            if (new_left == ALEA_NODE_ID_INVALID || new_right == ALEA_NODE_ID_INVALID) {
                return ALEA_NODE_ID_INVALID;
            }
            
            if (result_op == ALEA_OP_UNION) {
                return alea_create_union(sys, new_left, new_right);
            } else {
                return alea_create_intersection(sys, new_left, new_right);
            }
        }
        
        default:
            return ALEA_NODE_ID_INVALID;
    }
}

// ============================================================================
// TERM COLLECTION FOR FLATTENING
// ============================================================================

ALEA_VEC_DEFINE(node_id_vec, alea_node_id_t);

typedef struct {
    node_id_vec_t vec;
} node_list_t;

static node_list_t* node_list_create(size_t initial_cap) {
    node_list_t* list = malloc(sizeof(node_list_t));
    if (!list) return NULL;
    alea_vec_init(&list->vec);
    int res = alea_vec_reserve(&list->vec, initial_cap, alea_node_id_t);
    if (res != 0) { free(list); return NULL; }
    return list;
}

static void node_list_destroy(node_list_t* list) {
    if (!list) return;
    alea_vec_free(&list->vec);
    free(list);
}

static bool node_list_add(node_list_t* list, alea_node_id_t node_id) {
    int res = alea_vec_push(&list->vec, node_id, alea_node_id_t);
    return res == 0;
}

/**
 * Collect all terms from a chain of same-type operations.
 * Flattens (A ∩ B) ∩ C into [A, B, C].
 */
static void collect_terms(
    const alea_system_t* sys,
    alea_node_id_t node_id,
    alea_operation_t target_op,
    node_list_t* list
) {
    if (node_id == ALEA_NODE_ID_INVALID || node_id >= alea_vec_count(&sys->nodes)) return;
    
    const alea_node_t* node = &sys->nodes.data[node_id];
    alea_operation_t op = ALEA_GET_OPERATION(node);
    
    if (op == target_op) {
        collect_terms(sys, node->operation.left, target_op, list);
        collect_terms(sys, node->operation.right, target_op, list);
    } else {
        node_list_add(list, node_id);
    }
}

// ============================================================================
// CONTRADICTION AND TAUTOLOGY DETECTION
// ============================================================================

/**
 * Check if two nodes are complements of each other.
 * They must be primitives with the same primitive_id but opposite sense.
 */
static bool nodes_are_complements(const alea_system_t* sys, alea_node_id_t a, alea_node_id_t b) {
    if (a == ALEA_NODE_ID_INVALID || b == ALEA_NODE_ID_INVALID) return false;
    if (a >= alea_vec_count(&sys->nodes) || b >= alea_vec_count(&sys->nodes)) return false;
    
    const alea_node_t* na = &sys->nodes.data[a];
    const alea_node_t* nb = &sys->nodes.data[b];
    
    if (ALEA_GET_OPERATION(na) != ALEA_OP_PRIMITIVE) return false;
    if (ALEA_GET_OPERATION(nb) != ALEA_OP_PRIMITIVE) return false;
    
    return (na->primitive.primitive_id == nb->primitive.primitive_id) &&
           (na->primitive.sense == -nb->primitive.sense);
}

/**
 * Check if a list of terms (for intersection) contains a contradiction.
 * Returns true if A and ¬A both appear.
 */
static bool has_contradiction(const alea_system_t* sys, const node_list_t* list) {
    for (size_t i = 0; i < list->vec.count; i++) {
        for (size_t j = i + 1; j < list->vec.count; j++) {
            if (nodes_are_complements(sys, list->vec.data[i], list->vec.data[j])) {
                return true;
            }
        }
    }
    return false;
}

/**
 * Check if a list of terms (for union) contains a tautology.
 * Returns true if A and ¬A both appear.
 */
static bool has_tautology(const alea_system_t* sys, const node_list_t* list) {
    return has_contradiction(sys, list);
}

// ============================================================================
// SEMANTIC DEDUPLICATION (same primitive_id + sense = duplicate)
// ============================================================================

typedef struct {
    alea_node_id_t node_id;
    alea_primitive_id_t primitive_id;
    int8_t sense;
    bool is_primitive;
} dedup_key_t;

// ============================================================================
// LITERAL SET DATA STRUCTURES (for union optimization)
// ============================================================================

typedef struct {
    alea_primitive_id_t prim_id;
    int8_t eff_sense;
} literal_t;

typedef struct {
    literal_t* items;
    size_t count;
} literal_set_t;

static int compare_literals(const void* a_ptr, const void* b_ptr) {
    const literal_t* a = (const literal_t*)a_ptr;
    const literal_t* b = (const literal_t*)b_ptr;
    if (a->prim_id < b->prim_id) return -1;
    if (a->prim_id > b->prim_id) return 1;
    if (a->eff_sense < b->eff_sense) return -1;
    if (a->eff_sense > b->eff_sense) return 1;
    return 0;
}

/**
 * Extract sorted literal set from a node.
 * - ALEA_OP_PRIMITIVE → single-element set
 * - ALEA_OP_INTERSECTION → collect all primitive children (returns false if any non-primitive)
 * - Other → false
 */
static bool extract_literal_set(const alea_system_t* sys, alea_node_id_t node_id, literal_set_t* out) {
    out->items = NULL;
    out->count = 0;

    if (node_id == ALEA_NODE_ID_INVALID || node_id >= alea_vec_count(&sys->nodes))
        return false;

    const alea_node_t* node = &sys->nodes.data[node_id];
    alea_operation_t op = ALEA_GET_OPERATION(node);

    if (op == ALEA_OP_PRIMITIVE) {
        out->items = malloc(sizeof(literal_t));
        if (!out->items) return false;
        int8_t eff = node->primitive.sense;
        if (node->primitive.inverted) eff = -eff;
        out->items[0].prim_id = node->primitive.primitive_id;
        out->items[0].eff_sense = eff;
        out->count = 1;
        return true;
    }

    if (op == ALEA_OP_INTERSECTION) {
        node_list_t* terms = node_list_create(16);
        if (!terms) return false;
        collect_terms(sys, node_id, ALEA_OP_INTERSECTION, terms);

        out->items = malloc(terms->vec.count * sizeof(literal_t));
        if (!out->items) { node_list_destroy(terms); return false; }

        for (size_t i = 0; i < terms->vec.count; i++) {
            const alea_node_t* tn = &sys->nodes.data[terms->vec.data[i]];
            if (ALEA_GET_OPERATION(tn) != ALEA_OP_PRIMITIVE) {
                free(out->items);
                out->items = NULL;
                out->count = 0;
                node_list_destroy(terms);
                return false;
            }
            int8_t eff = tn->primitive.sense;
            if (tn->primitive.inverted) eff = -eff;
            out->items[i].prim_id = tn->primitive.primitive_id;
            out->items[i].eff_sense = eff;
        }
        out->count = terms->vec.count;
        node_list_destroy(terms);
        qsort(out->items, out->count, sizeof(literal_t), compare_literals);
        return true;
    }

    return false;
}

static void literal_set_free(literal_set_t* set) {
    free(set->items);
    set->items = NULL;
    set->count = 0;
}

/**
 * Check if set A is a subset of set B (both sorted).
 * O(|A| + |B|) merge.
 */
static bool literal_set_is_subset(const literal_set_t* a, const literal_set_t* b) {
    if (a->count > b->count) return false;
    size_t ia = 0, ib = 0;
    while (ia < a->count && ib < b->count) {
        int c = compare_literals(&a->items[ia], &b->items[ib]);
        if (c == 0) { ia++; ib++; }
        else if (c > 0) { ib++; }
        else { return false; }  /* element in A not found in B */
    }
    return ia == a->count;
}

/**
 * Compute intersection of two sorted literal sets.
 * Result is sorted. Returns true if non-empty intersection.
 */
static bool literal_set_intersect(const literal_set_t* a, const literal_set_t* b, literal_set_t* result) {
    size_t max_count = (a->count < b->count) ? a->count : b->count;
    result->items = malloc(max_count * sizeof(literal_t));
    if (!result->items) { result->count = 0; return false; }
    result->count = 0;

    size_t ia = 0, ib = 0;
    while (ia < a->count && ib < b->count) {
        int c = compare_literals(&a->items[ia], &b->items[ib]);
        if (c == 0) {
            result->items[result->count++] = a->items[ia];
            ia++; ib++;
        } else if (c < 0) { ia++; }
        else { ib++; }
    }
    return result->count > 0;
}

/**
 * Compute set difference: elements in A that are not in B (both sorted).
 */
static void literal_set_difference(const literal_set_t* a, const literal_set_t* b, literal_set_t* result) {
    result->items = malloc(a->count * sizeof(literal_t));
    if (!result->items) { result->count = 0; return; }
    result->count = 0;

    size_t ia = 0, ib = 0;
    while (ia < a->count) {
        if (ib >= b->count) {
            result->items[result->count++] = a->items[ia++];
        } else {
            int c = compare_literals(&a->items[ia], &b->items[ib]);
            if (c == 0) { ia++; ib++; }
            else if (c < 0) { result->items[result->count++] = a->items[ia++]; }
            else { ib++; }
        }
    }
}

// ============================================================================
// DEDUP KEYS
// ============================================================================

static int compare_dedup_keys(const void* a_ptr, const void* b_ptr) {
    const dedup_key_t* a = (const dedup_key_t*)a_ptr;
    const dedup_key_t* b = (const dedup_key_t*)b_ptr;

    // Group: all primitives first, then non-primitives
    if (a->is_primitive != b->is_primitive) {
        return a->is_primitive ? -1 : 1;
    }

    if (a->is_primitive) {
        // Primitives: sort by (primitive_id, sense)
        if (a->primitive_id < b->primitive_id) return -1;
        if (a->primitive_id > b->primitive_id) return 1;
        if (a->sense < b->sense) return -1;
        if (a->sense > b->sense) return 1;
        return 0;
    } else {
        // Non-primitives: sort by node_id
        if (a->node_id < b->node_id) return -1;
        if (a->node_id > b->node_id) return 1;
        return 0;
    }
}

/**
 * Remove semantically duplicate terms from list.
 *
 * After flatten, different node IDs can reference the same (primitive_id, sense).
 * This compares by (primitive_id, sense) for primitives, by node_id for non-primitives.
 *
 * A secondary pass also deduplicates by (mc_surface_id, sense) to catch cases
 * where primitive dedup didn't merge identical geometry (e.g., floating-point
 * differences from transform application).
 *
 * Returns number of duplicates removed.
 */
static size_t remove_semantic_duplicates(const alea_system_t* sys, node_list_t* list) {
    if (list->vec.count < 2) return 0;

    dedup_key_t* keys = malloc(list->vec.count * sizeof(dedup_key_t));
    if (!keys) return 0;

    for (size_t i = 0; i < list->vec.count; i++) {
        alea_node_id_t nid = list->vec.data[i];
        const alea_node_t* node = &sys->nodes.data[nid];
        keys[i].node_id = nid;

        if (ALEA_GET_OPERATION(node) == ALEA_OP_PRIMITIVE) {
            keys[i].is_primitive = true;
            keys[i].primitive_id = node->primitive.primitive_id;
            // Use effective sense (accounting for inverted flag) to match exporter behavior
            int8_t eff_sense = node->primitive.sense;
            if (node->primitive.inverted) eff_sense = -eff_sense;
            keys[i].sense = eff_sense;
        } else {
            keys[i].is_primitive = false;
            keys[i].primitive_id = 0;
            keys[i].sense = 0;
        }
    }

    qsort(keys, list->vec.count, sizeof(dedup_key_t), compare_dedup_keys);

    // After sorting, duplicates are adjacent
    list->vec.data[0] = keys[0].node_id;
    size_t write = 1;
    size_t removed = 0;

    for (size_t read = 1; read < list->vec.count; read++) {
        bool dup = false;
        if (keys[read].is_primitive && keys[read - 1].is_primitive) {
            dup = (keys[read].primitive_id == keys[read - 1].primitive_id &&
                   keys[read].sense == keys[read - 1].sense);
        } else if (!keys[read].is_primitive && !keys[read - 1].is_primitive) {
            dup = (keys[read].node_id == keys[read - 1].node_id);
        }

        if (!dup) {
            list->vec.data[write++] = keys[read].node_id;
        } else {
            removed++;
        }
    }

    list->vec.count = write;
    free(keys);

    // Secondary pass: geometric identity dedup.
    // Compares actual primitive coefficients (via alea_primitives_equal) to catch
    // different primitive_ids that represent the same surface. This handles cases
    // where transform roundtrips or separate creation paths produced duplicate
    // primitives that the hash-based dedup at creation time didn't merge.
    if (list->vec.count >= 2) {
        size_t write2 = 0;
        for (size_t i = 0; i < list->vec.count; i++) {
            alea_node_id_t nid = list->vec.data[i];
            const alea_node_t* node = &sys->nodes.data[nid];

            if (ALEA_GET_OPERATION(node) != ALEA_OP_PRIMITIVE) {
                list->vec.data[write2++] = nid;
                continue;
            }

            uint32_t prim_id = node->primitive.primitive_id;
            const alea_primitive_entry_t* prim = &sys->primitives.data[prim_id];
            int8_t eff_sense = node->primitive.sense;
            if (node->primitive.inverted) eff_sense = -eff_sense;

            bool dup = false;
            for (size_t j = 0; j < write2; j++) {
                const alea_node_t* kept = &sys->nodes.data[list->vec.data[j]];
                if (ALEA_GET_OPERATION(kept) != ALEA_OP_PRIMITIVE) continue;

                int8_t k_sense = kept->primitive.sense;
                if (kept->primitive.inverted) k_sense = -k_sense;

                uint32_t k_prim_id = kept->primitive.primitive_id;
                if (k_prim_id == prim_id) {
                    // Same primitive — already handled by primary pass
                    if (k_sense == eff_sense) { dup = true; break; }
                    continue;
                }

                const alea_primitive_entry_t* k_prim = &sys->primitives.data[k_prim_id];
                int8_t match_inv = 0;
                if (alea_primitives_equal(prim->type, &prim->data,
                                         k_prim->type, &k_prim->data,
                                         &sys->config, &match_inv)) {
                    // Geometrically identical surfaces
                    int8_t adjusted_sense = eff_sense;
                    if (match_inv) adjusted_sense = -adjusted_sense;

                    if (adjusted_sense == k_sense) {
                        dup = true;  // Same half-space → idempotent
                        break;
                    }
                    // Opposite effective sense on same surface → contradiction
                    // (handled elsewhere via complement detection)
                }
            }

            if (!dup) {
                list->vec.data[write2++] = nid;
            } else {
                removed++;
            }
        }
        list->vec.count = write2;
    }

    return removed;
}

// ============================================================================
// PAIRWISE CONTAINMENT CHECKS
// ============================================================================

/**
 * Check if two parallel planes with same sense have a containment relationship.
 *
 * For same-direction normals and same sense:
 * - sense < 0 (below plane): tighter = larger normalized d
 * - sense > 0 (above plane): tighter = smaller normalized d
 *
 * Returns: 1 if B is redundant, -1 if A is redundant, 0 if no subsumption.
 */
static int check_parallel_planes(
    const alea_plane_data_t* pa, int8_t sense_a,
    const alea_plane_data_t* pb, int8_t sense_b
) {
    if (sense_a != sense_b) return 0;

    double mag_a = sqrt(pa->a * pa->a + pa->b * pa->b + pa->c * pa->c);
    double mag_b = sqrt(pb->a * pb->a + pb->b * pb->b + pb->c * pb->c);
    if (mag_a < 1e-12 || mag_b < 1e-12) return 0;

    // Check normals are parallel: |cross product|² ≈ 0
    double cx = pa->b * pb->c - pa->c * pb->b;
    double cy = pa->c * pb->a - pa->a * pb->c;
    double cz = pa->a * pb->b - pa->b * pb->a;
    double cross_sq = cx * cx + cy * cy + cz * cz;
    double threshold = 1e-16 * mag_a * mag_a * mag_b * mag_b;

    if (cross_sq > threshold) return 0;

    // Check same direction (canonicalized planes should have same direction)
    double dot = pa->a * pb->a + pa->b * pb->b + pa->c * pb->c;
    if (dot < 0) return 0;

    // Normalized offsets along the normal direction
    double da = pa->d / mag_a;
    double db = pb->d / mag_b;

    double diff = da - db;

    if (sense_a < 0) {
        // Region: n̂·x < -d/|n|. Tighter = larger d/|n|.
        // Equal offsets (same plane) → B is redundant (idempotent).
        return (diff >= -1e-10) ? 1 : -1;
    } else {
        // Region: n̂·x > -d/|n|. Tighter = smaller d/|n|.
        return (diff <= 1e-10) ? 1 : -1;
    }
}

/**
 * Check coaxial cylinders (same center in perpendicular plane) for containment.
 *
 * Returns: 1 if B is redundant, -1 if A is redundant, 0 if no subsumption.
 */
static int check_coaxial_cylinders(
    double cx_a, double cy_a, double r_a, int8_t sense_a,
    double cx_b, double cy_b, double r_b, int8_t sense_b
) {
    if (sense_a != sense_b) return 0;
    if (fabs(cx_a - cx_b) > 1e-10 || fabs(cy_a - cy_b) > 1e-10) return 0;
    if (fabs(r_a - r_b) < 1e-10) return 0;

    if (sense_a < 0) {
        // Inside: smaller radius is tighter, larger is redundant
        return (r_a < r_b) ? 1 : -1;
    } else {
        // Outside: larger radius is tighter, smaller is redundant
        return (r_a > r_b) ? 1 : -1;
    }
}

/**
 * Check sphere containment.
 * If one sphere is entirely inside the other, the outer one is redundant (for same sense).
 *
 * Returns: 1 if B is redundant, -1 if A is redundant, 0 if no subsumption.
 */
static int check_sphere_containment(
    const alea_sphere_data_t* sa, int8_t sense_a,
    const alea_sphere_data_t* sb, int8_t sense_b
) {
    if (sense_a != sense_b) return 0;

    double dx = sa->center_x - sb->center_x;
    double dy = sa->center_y - sb->center_y;
    double dz = sa->center_z - sb->center_z;
    double dist = sqrt(dx * dx + dy * dy + dz * dz);

    if (sense_a < 0) {
        // Inside: if SA ⊂ SB (ball A inside ball B), SB is redundant
        if (dist + sa->radius <= sb->radius + 1e-10) return 1;
        if (dist + sb->radius <= sa->radius + 1e-10) return -1;
    } else {
        // Outside: if SA ⊂ SB, then outside(SB) ⊂ outside(SA), so A is redundant
        if (dist + sa->radius <= sb->radius + 1e-10) return -1;
        if (dist + sb->radius <= sa->radius + 1e-10) return 1;
    }

    return 0;
}

/**
 * Dispatch to specialized pairwise containment check based on primitive type.
 *
 * Returns: 1 if B is redundant, -1 if A is redundant, 0 if no subsumption detected.
 */
static int check_typed_pairwise_subsumption(
    const alea_primitive_entry_t* prim_a, int8_t sense_a,
    const alea_primitive_entry_t* prim_b, int8_t sense_b
) {
    if (prim_a->type != prim_b->type) return 0;

    switch (prim_a->type) {
        case ALEA_PRIMITIVE_PLANE:
            return check_parallel_planes(
                &prim_a->data.plane, sense_a,
                &prim_b->data.plane, sense_b);

        case ALEA_PRIMITIVE_SPHERE:
            return check_sphere_containment(
                &prim_a->data.sphere, sense_a,
                &prim_b->data.sphere, sense_b);

        case ALEA_PRIMITIVE_CYLINDER_Z:
            return check_coaxial_cylinders(
                prim_a->data.cyl_z.center_x, prim_a->data.cyl_z.center_y,
                prim_a->data.cyl_z.radius, sense_a,
                prim_b->data.cyl_z.center_x, prim_b->data.cyl_z.center_y,
                prim_b->data.cyl_z.radius, sense_b);

        case ALEA_PRIMITIVE_CYLINDER_X:
            return check_coaxial_cylinders(
                prim_a->data.cyl_x.center_y, prim_a->data.cyl_x.center_z,
                prim_a->data.cyl_x.radius, sense_a,
                prim_b->data.cyl_x.center_y, prim_b->data.cyl_x.center_z,
                prim_b->data.cyl_x.radius, sense_b);

        case ALEA_PRIMITIVE_CYLINDER_Y:
            return check_coaxial_cylinders(
                prim_a->data.cyl_y.center_x, prim_a->data.cyl_y.center_z,
                prim_a->data.cyl_y.radius, sense_a,
                prim_b->data.cyl_y.center_x, prim_b->data.cyl_y.center_z,
                prim_b->data.cyl_y.radius, sense_b);

        default:
            return 0;
    }
}

/**
 * General pairwise subsumption check using interval arithmetic.
 *
 * Uses one primitive's bbox to evaluate the other's surface equation.
 * If B is satisfied everywhere in A's bbox, then A implies B and B is redundant.
 *
 * Returns: 1 if B is redundant, -1 if A is redundant, 0 if no subsumption detected.
 */
static int check_general_pairwise_subsumption(
    const alea_system_t* sys,
    alea_node_id_t nid_a, const alea_node_t* na,
    alea_node_id_t nid_b, const alea_node_t* nb
) {
    const double UNBOUNDED_THRESHOLD = 9e5;

    // Try: does A imply B? (evaluate B's surface over A's bbox)
    {
        const alea_bbox_t ba_v = alea_node_bbox_get(&na->bbox);
        const alea_bbox_t* ba = &ba_v;
        if (ba->min_x <= ba->max_x &&
            ba->max_x - ba->min_x < UNBOUNDED_THRESHOLD &&
            ba->max_y - ba->min_y < UNBOUNDED_THRESHOLD &&
            ba->max_z - ba->min_z < UNBOUNDED_THRESHOLD) {

            box_surface_relation_t rel = check_box_surface_relation(sys, nid_b, ba);
            int8_t sense_b = nb->primitive.sense;
            if (nb->primitive.inverted) sense_b = -sense_b;

            if ((sense_b < 0 && rel == BOX_STRICTLY_NEGATIVE) ||
                (sense_b > 0 && rel == BOX_STRICTLY_POSITIVE)) {
                return 1;  // B satisfied everywhere in A's bbox → B redundant
            }
        }
    }

    // Try: does B imply A? (evaluate A's surface over B's bbox)
    {
        const alea_bbox_t bb_v = alea_node_bbox_get(&nb->bbox);
        const alea_bbox_t* bb = &bb_v;
        if (bb->min_x <= bb->max_x &&
            bb->max_x - bb->min_x < UNBOUNDED_THRESHOLD &&
            bb->max_y - bb->min_y < UNBOUNDED_THRESHOLD &&
            bb->max_z - bb->min_z < UNBOUNDED_THRESHOLD) {

            box_surface_relation_t rel = check_box_surface_relation(sys, nid_a, bb);
            int8_t sense_a = na->primitive.sense;
            if (na->primitive.inverted) sense_a = -sense_a;

            if ((sense_a < 0 && rel == BOX_STRICTLY_NEGATIVE) ||
                (sense_a > 0 && rel == BOX_STRICTLY_POSITIVE)) {
                return -1;  // A satisfied everywhere in B's bbox → A redundant
            }
        }
    }

    return 0;
}

/**
 * Remove terms subsumed by other terms in an intersection.
 *
 * Uses O(n²) pairwise checks: specialized (planes, cylinders, spheres)
 * and general (interval arithmetic over bboxes).
 *
 * Returns number of terms removed.
 */
static size_t remove_subsumed_terms(
    const alea_system_t* sys,
    node_list_t* terms,
    alea_simplify_stats_t* stats
) {
    if (terms->vec.count < 2) return 0;

    alea_bitset_t redundant = alea_bitset_create(terms->vec.count);
    if (!redundant.words) return 0;

    size_t removed = 0;

    for (size_t i = 0; i < terms->vec.count; i++) {
        if (alea_bitset_test(&redundant, i)) continue;

        alea_node_id_t nid_i = terms->vec.data[i];
        const alea_node_t* ni = &sys->nodes.data[nid_i];
        if (ALEA_GET_OPERATION(ni) != ALEA_OP_PRIMITIVE) continue;

        alea_primitive_id_t pid_i = ni->primitive.primitive_id;
        int8_t sense_i = ni->primitive.sense;
        if (ni->primitive.inverted) sense_i = -sense_i;
        const alea_primitive_entry_t* prim_i = &sys->primitives.data[pid_i];

        for (size_t j = i + 1; j < terms->vec.count; j++) {
            if (alea_bitset_test(&redundant, j)) continue;

            alea_node_id_t nid_j = terms->vec.data[j];
            const alea_node_t* nj = &sys->nodes.data[nid_j];
            if (ALEA_GET_OPERATION(nj) != ALEA_OP_PRIMITIVE) continue;

            alea_primitive_id_t pid_j = nj->primitive.primitive_id;
            int8_t sense_j = nj->primitive.sense;
            if (nj->primitive.inverted) sense_j = -sense_j;
            const alea_primitive_entry_t* prim_j = &sys->primitives.data[pid_j];

            // Try specialized check first
            int result = check_typed_pairwise_subsumption(
                prim_i, sense_i, prim_j, sense_j);

            // Fall back to general interval arithmetic
            if (result == 0) {
                result = check_general_pairwise_subsumption(
                    sys, nid_i, ni, nid_j, nj);
            }

            if (result == 1) {
                // B (j) is redundant
                alea_bitset_set(&redundant, j);
                removed++;
            } else if (result == -1) {
                // A (i) is redundant
                alea_bitset_set(&redundant, i);
                removed++;
                break;
            }
        }
    }

    if (removed > 0) {
        size_t write = 0;
        for (size_t i = 0; i < terms->vec.count; i++) {
            if (!alea_bitset_test(&redundant, i)) {
                terms->vec.data[write++] = terms->vec.data[i];
            }
        }
        terms->vec.count = write;
        if (stats) stats->absorption_reductions += removed;
    }

    alea_bitset_destroy(&redundant);
    return removed;
}

// Forward declarations for union optimization
static alea_node_id_t build_balanced_tree(
    alea_system_t* sys, const alea_node_id_t* nodes, size_t count, alea_operation_t op);
static bool cell_is_provably_empty(
    const alea_system_t* sys, alea_node_id_t root, const alea_bbox_t* box, int depth);

// ============================================================================
// UNION OPTIMIZATION PHASES
// ============================================================================

/**
 * Phase 1: Algebraic Absorption for unions.
 *
 * If literals(Ti) ⊆ literals(Tj), Ti has fewer constraints → Ti ⊇ Tj
 * → Tj is redundant in the union. O(k²·m).
 *
 * Returns number of branches removed.
 */
static size_t remove_absorbed_union_branches(
    const alea_system_t* sys,
    node_list_t* terms,
    alea_simplify_stats_t* stats
) {
    if (terms->vec.count < 2) return 0;

    /* Extract literal sets for all branches */
    literal_set_t* sets = calloc(terms->vec.count, sizeof(literal_set_t));
    alea_bitset_t valid = alea_bitset_create(terms->vec.count);
    if (!sets || !valid.words) { free(sets); alea_bitset_destroy(&valid); return 0; }

    for (size_t i = 0; i < terms->vec.count; i++) {
        if (extract_literal_set(sys, terms->vec.data[i], &sets[i]))
            alea_bitset_set(&valid, i);
    }

    alea_bitset_t redundant = alea_bitset_create(terms->vec.count);
    if (!redundant.words) {
        for (size_t i = 0; i < terms->vec.count; i++) if (alea_bitset_test(&valid, i)) literal_set_free(&sets[i]);
        free(sets); alea_bitset_destroy(&valid);
        return 0;
    }

    size_t removed = 0;

    for (size_t i = 0; i < terms->vec.count; i++) {
        if (!alea_bitset_test(&valid, i) || alea_bitset_test(&redundant, i)) continue;
        for (size_t j = 0; j < terms->vec.count; j++) {
            if (i == j || !alea_bitset_test(&valid, j) || alea_bitset_test(&redundant, j)) continue;

            /* If Ti ⊆ Tj, then Ti (fewer constraints) absorbs Tj */
            if (literal_set_is_subset(&sets[i], &sets[j])) {
                alea_bitset_set(&redundant, j);
                removed++;
            }
        }
    }

    if (removed > 0) {
        size_t write = 0;
        for (size_t i = 0; i < terms->vec.count; i++) {
            if (!alea_bitset_test(&redundant, i)) {
                terms->vec.data[write++] = terms->vec.data[i];
            }
        }
        terms->vec.count = write;
        if (stats) stats->union_branches_absorbed += removed;
    }

    /* Clean up: free all literal sets (including redundant ones) */
    size_t orig_count = terms->vec.count + removed;
    for (size_t i = 0; i < orig_count; i++) {
        if (alea_bitset_test(&valid, i)) literal_set_free(&sets[i]);
    }
    free(sets);
    alea_bitset_destroy(&valid);
    alea_bitset_destroy(&redundant);
    return removed;
}

/**
 * Phase 2: Common Factor Extraction for unions.
 *
 * (A∩B∩C) ∪ (A∩B∩D) → A ∩ B ∩ (C ∪ D)
 *
 * Requires ALL branches to have valid literal sets.
 * Returns new node or ALEA_NODE_ID_INVALID if factoring didn't apply.
 */
static alea_node_id_t factor_common_union_literals(
    alea_system_t* sys,
    node_list_t* terms,
    alea_simplify_stats_t* stats
) {
    if (terms->vec.count < 2) return ALEA_NODE_ID_INVALID;

    /* Extract literal sets for all branches */
    literal_set_t* sets = calloc(terms->vec.count, sizeof(literal_set_t));
    if (!sets) return ALEA_NODE_ID_INVALID;

    for (size_t i = 0; i < terms->vec.count; i++) {
        if (!extract_literal_set(sys, terms->vec.data[i], &sets[i])) {
            /* Not all branches are pure literal sets — can't factor */
            for (size_t j = 0; j <= i; j++) literal_set_free(&sets[j]);
            free(sets);
            return ALEA_NODE_ID_INVALID;
        }
    }

    /* Compute intersection of all literal sets */
    literal_set_t common = {0};
    /* Start with first set, intersect with each subsequent */
    common.items = malloc(sets[0].count * sizeof(literal_t));
    if (!common.items) {
        for (size_t i = 0; i < terms->vec.count; i++) literal_set_free(&sets[i]);
        free(sets);
        return ALEA_NODE_ID_INVALID;
    }
    memcpy(common.items, sets[0].items, sets[0].count * sizeof(literal_t));
    common.count = sets[0].count;

    for (size_t i = 1; i < terms->vec.count && common.count > 0; i++) {
        literal_set_t new_common;
        literal_set_intersect(&common, &sets[i], &new_common);
        free(common.items);
        common = new_common;
    }

    if (common.count == 0) {
        /* No common factors */
        literal_set_free(&common);
        for (size_t i = 0; i < terms->vec.count; i++) literal_set_free(&sets[i]);
        free(sets);
        return ALEA_NODE_ID_INVALID;
    }

    /* Build the common-factor intersection node */
    /* For each common literal, find its node in the first branch */
    alea_node_id_t* common_nodes = malloc(common.count * sizeof(alea_node_id_t));
    if (!common_nodes) {
        literal_set_free(&common);
        for (size_t i = 0; i < terms->vec.count; i++) literal_set_free(&sets[i]);
        free(sets);
        return ALEA_NODE_ID_INVALID;
    }

    /* Find existing nodes for common literals from first branch */
    {
        node_list_t* first_terms = node_list_create(16);
        if (!first_terms) {
            free(common_nodes);
            literal_set_free(&common);
            for (size_t i = 0; i < terms->vec.count; i++) literal_set_free(&sets[i]);
            free(sets);
            return ALEA_NODE_ID_INVALID;
        }

        const alea_node_t* first_node = &sys->nodes.data[terms->vec.data[0]];
        alea_operation_t first_op = ALEA_GET_OPERATION(first_node);
        if (first_op == ALEA_OP_INTERSECTION) {
            collect_terms(sys, terms->vec.data[0], ALEA_OP_INTERSECTION, first_terms);
        } else {
            /* Single primitive */
            node_list_add(first_terms, terms->vec.data[0]);
        }

        for (size_t c = 0; c < common.count; c++) {
            common_nodes[c] = ALEA_NODE_ID_INVALID;
            for (size_t t = 0; t < first_terms->vec.count; t++) {
                const alea_node_t* tn = &sys->nodes.data[first_terms->vec.data[t]];
                if (ALEA_GET_OPERATION(tn) != ALEA_OP_PRIMITIVE) continue;
                int8_t eff = tn->primitive.sense;
                if (tn->primitive.inverted) eff = -eff;
                if (tn->primitive.primitive_id == common.items[c].prim_id &&
                    eff == common.items[c].eff_sense) {
                    common_nodes[c] = first_terms->vec.data[t];
                    break;
                }
            }
        }
        node_list_destroy(first_terms);
    }

    /* Build remainder union: for each branch, remove common literals */
    alea_node_id_t* remainder_nodes = malloc(terms->vec.count * sizeof(alea_node_id_t));
    if (!remainder_nodes) {
        free(common_nodes);
        literal_set_free(&common);
        for (size_t i = 0; i < terms->vec.count; i++) literal_set_free(&sets[i]);
        free(sets);
        return ALEA_NODE_ID_INVALID;
    }

    bool has_empty_remainder = false;
    size_t remainder_count = 0;

    for (size_t i = 0; i < terms->vec.count; i++) {
        literal_set_t diff;
        literal_set_difference(&sets[i], &common, &diff);

        if (diff.count == 0) {
            /* Branch is entirely common factors — union collapses */
            has_empty_remainder = true;
            literal_set_free(&diff);
            break;
        }

        /* Build node for remainder literals */
        /* Find nodes from branch i for each diff literal */
        node_list_t* branch_terms = node_list_create(16);
        if (!branch_terms) {
            literal_set_free(&diff);
            free(remainder_nodes);
            free(common_nodes);
            literal_set_free(&common);
            for (size_t j = 0; j < terms->vec.count; j++) literal_set_free(&sets[j]);
            free(sets);
            return ALEA_NODE_ID_INVALID;
        }

        const alea_node_t* bnode = &sys->nodes.data[terms->vec.data[i]];
        alea_operation_t bop = ALEA_GET_OPERATION(bnode);
        if (bop == ALEA_OP_INTERSECTION) {
            collect_terms(sys, terms->vec.data[i], ALEA_OP_INTERSECTION, branch_terms);
        } else {
            node_list_add(branch_terms, terms->vec.data[i]);
        }

        alea_node_id_t* diff_nodes = malloc(diff.count * sizeof(alea_node_id_t));
        if (!diff_nodes) {
            node_list_destroy(branch_terms);
            literal_set_free(&diff);
            free(remainder_nodes);
            free(common_nodes);
            literal_set_free(&common);
            for (size_t j = 0; j < terms->vec.count; j++) literal_set_free(&sets[j]);
            free(sets);
            return ALEA_NODE_ID_INVALID;
        }

        for (size_t d = 0; d < diff.count; d++) {
            diff_nodes[d] = ALEA_NODE_ID_INVALID;
            for (size_t t = 0; t < branch_terms->vec.count; t++) {
                const alea_node_t* tn = &sys->nodes.data[branch_terms->vec.data[t]];
                if (ALEA_GET_OPERATION(tn) != ALEA_OP_PRIMITIVE) continue;
                int8_t eff = tn->primitive.sense;
                if (tn->primitive.inverted) eff = -eff;
                if (tn->primitive.primitive_id == diff.items[d].prim_id &&
                    eff == diff.items[d].eff_sense) {
                    diff_nodes[d] = branch_terms->vec.data[t];
                    break;
                }
            }
        }

        remainder_nodes[remainder_count++] = build_balanced_tree(sys, diff_nodes, diff.count, ALEA_OP_INTERSECTION);
        free(diff_nodes);
        node_list_destroy(branch_terms);
        literal_set_free(&diff);
    }

    /* Save common count before cleanup */
    size_t common_count = common.count;

    /* Clean up literal sets */
    literal_set_free(&common);
    for (size_t i = 0; i < terms->vec.count; i++) literal_set_free(&sets[i]);
    free(sets);

    alea_node_id_t result;

    if (has_empty_remainder) {
        /* One branch was entirely common factors → union collapses to just the common intersection */
        result = build_balanced_tree(sys, common_nodes, common_count, ALEA_OP_INTERSECTION);
    } else {
        /* Build: common_factors ∩ (remainder_union) */
        alea_node_id_t remainder_union = build_balanced_tree(sys, remainder_nodes, remainder_count, ALEA_OP_UNION);
        alea_node_id_t common_inter = build_balanced_tree(sys, common_nodes, common_count, ALEA_OP_INTERSECTION);
        result = alea_create_intersection(sys, common_inter, remainder_union);
    }

    free(common_nodes);
    free(remainder_nodes);

    if (stats) stats->union_common_factors++;
    return result;
}

/**
 * Phase 3: Geometric Branch Subsumption for unions.
 *
 * For each branch Ti, check if it's geometrically subsumed by the union of others.
 * Build difference(Ti, union(others)) and test with cell_is_provably_empty().
 *
 * Returns number of branches removed.
 */
static size_t remove_geometrically_subsumed_branches(
    alea_system_t* sys,
    node_list_t* terms,
    alea_simplify_stats_t* stats
) {
    if (terms->vec.count < 2) return 0;

    alea_bitset_t redundant = alea_bitset_create(terms->vec.count);
    if (!redundant.words) return 0;

    size_t removed = 0;

    for (size_t i = 0; i < terms->vec.count; i++) {
        if (alea_bitset_test(&redundant, i)) continue;

        /* Build union of all other non-redundant branches */
        node_list_t* others = node_list_create(terms->vec.count);
        if (!others) continue;

        for (size_t j = 0; j < terms->vec.count; j++) {
            if (j != i && !alea_bitset_test(&redundant, j)) {
                node_list_add(others, terms->vec.data[j]);
            }
        }

        if (others->vec.count == 0) {
            node_list_destroy(others);
            continue;
        }

        alea_node_id_t others_union = build_balanced_tree(
            sys, others->vec.data, others->vec.count, ALEA_OP_UNION);
        node_list_destroy(others);

        /* Build difference: Ti - union(others) = Ti ∩ ¬union(others) */
        alea_node_id_t ti = terms->vec.data[i];
        alea_node_id_t diff = alea_create_difference(sys, ti, others_union);

        /* Get bbox for Ti */
        const alea_node_t* ti_node = &sys->nodes.data[ti];
        const alea_bbox_t ti_bbox_v = alea_node_bbox_get(&ti_node->bbox);
        const alea_bbox_t* ti_bbox = &ti_bbox_v;

        /* Skip if bbox is degenerate */
        if (ti_bbox->min_x > ti_bbox->max_x) continue;

        if (cell_is_provably_empty(sys, diff, ti_bbox, 3)) {
            alea_bitset_set(&redundant, i);
            removed++;
        }
    }

    if (removed > 0) {
        size_t write = 0;
        for (size_t i = 0; i < terms->vec.count; i++) {
            if (!alea_bitset_test(&redundant, i)) {
                terms->vec.data[write++] = terms->vec.data[i];
            }
        }
        terms->vec.count = write;
        if (stats) stats->union_branches_subsumed += removed;
    }

    alea_bitset_destroy(&redundant);
    return removed;
}

// ============================================================================
// BUILD BALANCED TREE FROM TERM LIST
// ============================================================================

static alea_node_id_t build_balanced_tree(
    alea_system_t* sys,
    const alea_node_id_t* nodes,
    size_t count,
    alea_operation_t op
) {
    if (count == 0) return ALEA_NODE_ID_INVALID;
    if (count == 1) return nodes[0];
    if (count == 2) {
        return (op == ALEA_OP_INTERSECTION) ?
               alea_create_intersection(sys, nodes[0], nodes[1]) :
               alea_create_union(sys, nodes[0], nodes[1]);
    }
    
    size_t mid = count / 2;
    alea_node_id_t left = build_balanced_tree(sys, nodes, mid, op);
    alea_node_id_t right = build_balanced_tree(sys, nodes + mid, count - mid, op);
    
    return (op == ALEA_OP_INTERSECTION) ?
           alea_create_intersection(sys, left, right) :
           alea_create_union(sys, left, right);
}

// ============================================================================
// FLATTEN AND OPTIMIZE
// ============================================================================

#define RESULT_EMPTY    ((alea_node_id_t)(ALEA_NODE_ID_INVALID - 1))
#define RESULT_UNIVERSE ((alea_node_id_t)(ALEA_NODE_ID_INVALID - 2))

/**
 * Check if a primitive term is geometrically redundant given a tight bounding box.
 *
 * A term is redundant in an intersection if:
 * - We require the positive half-space (sense > 0) AND box is strictly positive
 * - We require the negative half-space (sense < 0) AND box is strictly negative
 *
 * Returns: 1 = redundant (can remove), 0 = keep, -1 = contradiction (empty)
 */
static int check_geometric_redundancy(
    const alea_system_t* sys,
    alea_node_id_t prim_node_id,
    const alea_bbox_t* tight_box
) {
    const alea_node_t* node = &sys->nodes.data[prim_node_id];

    if (ALEA_GET_OPERATION(node) != ALEA_OP_PRIMITIVE) {
        return 0;  // Not a primitive, can't check
    }

    box_surface_relation_t rel = check_box_surface_relation(sys, prim_node_id, tight_box);
    int8_t sense = node->primitive.sense;
    if (node->primitive.inverted) sense = -sense;

    if (sense > 0) {
        // Require positive half-space (outside surface)
        if (rel == BOX_STRICTLY_POSITIVE) {
            return 1;  // Redundant: box already entirely in positive half-space
        }
        if (rel == BOX_STRICTLY_NEGATIVE) {
            return -1; // Contradiction: box entirely in negative, but we need positive
        }
    } else if (sense < 0) {
        // Require negative half-space (inside surface)
        if (rel == BOX_STRICTLY_NEGATIVE) {
            return 1;  // Redundant: box already entirely in negative half-space
        }
        if (rel == BOX_STRICTLY_POSITIVE) {
            return -1; // Contradiction: box entirely in positive, but we need negative
        }
    }

    return 0;  // BOX_INTERSECTS or sense == 0 - must keep
}

/**
 * Compute the tight intersection bounding box from primitive half-spaces.
 *
 * Both interior AND exterior half-spaces contribute useful one-sided bounds
 * (e.g., z>5 gives min_z=5 even though max_z is unbounded).  Each individual
 * bound is checked against BOUND_INF (1e9) to filter out unbounded sides.
 * Node bboxes already account for the inverted flag (effective sense).
 */
static bool compute_intersection_bbox(
    const alea_system_t* sys,
    const node_list_t* terms,
    alea_bbox_t* out_box
) {
    /* Threshold for "effectively infinite" — BBOX_LARGE is 1e10, so any
     * individual bound with absolute value beyond this is unbounded. */
    const double BOUND_INF = 1e9;

    // Start with infinite box
    out_box->min_x = -DBL_MAX;
    out_box->max_x = DBL_MAX;
    out_box->min_y = -DBL_MAX;
    out_box->max_y = DBL_MAX;
    out_box->min_z = -DBL_MAX;
    out_box->max_z = DBL_MAX;

    bool found_bound = false;

    for (size_t i = 0; i < terms->vec.count; i++) {
        alea_node_id_t tid = terms->vec.data[i];
        if (tid >= alea_vec_count(&sys->nodes)) continue;

        const alea_node_t* node = &sys->nodes.data[tid];
        if (ALEA_GET_OPERATION(node) != ALEA_OP_PRIMITIVE) continue;

        /* Use node's precomputed halfspace bbox.  This already accounts for
         * sense: e.g. plane z=5, sense>0 → bbox min_z=5; sense<0 → max_z=5.
         * Both interior AND exterior half-spaces can provide useful one-sided
         * bounds (planes always do, exterior of closed surfaces don't). */
        const alea_bbox_t nb_v = alea_node_bbox_get(&node->bbox);
        const alea_bbox_t* nb = &nb_v;
        if (nb->min_x > nb->max_x) continue;  /* degenerate */

        /* Check each bound independently.  A half-space typically provides
         * only ONE finite bound per axis (e.g. z>5 gives min_z=5 but
         * max_z=BBOX_LARGE).  The old code checked axis *extent* and
         * discarded the whole axis if one side was unbounded. */
        if (nb->min_x > -BOUND_INF && nb->min_x > out_box->min_x) {
            out_box->min_x = nb->min_x; found_bound = true;
        }
        if (nb->max_x < BOUND_INF && nb->max_x < out_box->max_x) {
            out_box->max_x = nb->max_x; found_bound = true;
        }
        if (nb->min_y > -BOUND_INF && nb->min_y > out_box->min_y) {
            out_box->min_y = nb->min_y; found_bound = true;
        }
        if (nb->max_y < BOUND_INF && nb->max_y < out_box->max_y) {
            out_box->max_y = nb->max_y; found_bound = true;
        }
        if (nb->min_z > -BOUND_INF && nb->min_z > out_box->min_z) {
            out_box->min_z = nb->min_z; found_bound = true;
        }
        if (nb->max_z < BOUND_INF && nb->max_z < out_box->max_z) {
            out_box->max_z = nb->max_z; found_bound = true;
        }
    }

    // Check if bbox is empty (contradiction)
    if (out_box->min_x > out_box->max_x ||
        out_box->min_y > out_box->max_y ||
        out_box->min_z > out_box->max_z) {
        return false;  // Empty intersection
    }

    // If no finite bounds found, use a large but finite default
    if (!found_bound) {
        double DEFAULT_BOUND = 1e6;
        out_box->min_x = -DEFAULT_BOUND;
        out_box->max_x = DEFAULT_BOUND;
        out_box->min_y = -DEFAULT_BOUND;
        out_box->max_y = DEFAULT_BOUND;
        out_box->min_z = -DEFAULT_BOUND;
        out_box->max_z = DEFAULT_BOUND;
    }

    return true;
}

// Forward declaration for mutual recursion
static alea_node_id_t simplify_recursive(
    alea_system_t* sys, alea_node_id_t node_id, alea_simplify_stats_t* stats);

static alea_node_id_t flatten_and_optimize(
    alea_system_t* sys,
    alea_node_id_t node_id,
    alea_simplify_stats_t* stats
) {
    if (node_id == ALEA_NODE_ID_INVALID || node_id >= alea_vec_count(&sys->nodes)) {
        return ALEA_NODE_ID_INVALID;
    }

    const alea_node_t* node = &sys->nodes.data[node_id];
    alea_operation_t op = ALEA_GET_OPERATION(node);

    if (op != ALEA_OP_INTERSECTION && op != ALEA_OP_UNION) {
        return node_id;
    }

    node_list_t* terms = node_list_create(32);
    if (!terms) return node_id;

    collect_terms(sys, node_id, op, terms);

    // Recursively simplify each collected term (handles nested union/intersection)
    {
        size_t write = 0;
        for (size_t i = 0; i < terms->vec.count; i++) {
            alea_node_id_t simplified = simplify_recursive(sys, terms->vec.data[i], stats);

            if (simplified == RESULT_EMPTY) {
                if (op == ALEA_OP_INTERSECTION) {
                    // Empty term in intersection → entire intersection is empty
                    node_list_destroy(terms);
                    return RESULT_EMPTY;
                }
                // Empty term in union → skip it
                continue;
            } else if (simplified == RESULT_UNIVERSE) {
                if (op == ALEA_OP_UNION) {
                    // Universe term in union → entire union is universe
                    node_list_destroy(terms);
                    return RESULT_UNIVERSE;
                }
                // Universe term in intersection → skip it
                continue;
            }

            terms->vec.data[write++] = simplified;
        }
        terms->vec.count = write;

        if (terms->vec.count == 0) {
            node_list_destroy(terms);
            return (op == ALEA_OP_INTERSECTION) ? RESULT_UNIVERSE : RESULT_EMPTY;
        }
    }

    // Check contradiction/tautology (algebraic)
    if (op == ALEA_OP_INTERSECTION && has_contradiction(sys, terms)) {
        if (stats) stats->contradictions_found++;
        node_list_destroy(terms);
        return RESULT_EMPTY;
    }

    if (op == ALEA_OP_UNION && has_tautology(sys, terms)) {
        if (stats) stats->tautologies_found++;
        node_list_destroy(terms);
        return RESULT_UNIVERSE;
    }

    // Remove semantic duplicates (idempotent rule: A ∩ A = A)
    // Uses (primitive_id, sense) comparison instead of node_id to catch
    // duplicates introduced by flatten (different node IDs, same surface+sense)
    size_t dups = remove_semantic_duplicates(sys, terms);
    if (stats && dups > 0) stats->idempotent_reductions += dups;

    // ========================================================================
    // UNION OPTIMIZATION (algebraic absorption, factoring, geometric subsumption)
    // ========================================================================
    if (op == ALEA_OP_UNION && terms->vec.count > 1) {
        remove_absorbed_union_branches(sys, terms, stats);
        if (terms->vec.count > 1) {
            alea_node_id_t factored = factor_common_union_literals(sys, terms, stats);
            if (factored != ALEA_NODE_ID_INVALID) {
                node_list_destroy(terms);
                return simplify_recursive(sys, factored, stats);
            }
        }
        if (terms->vec.count > 1) {
            remove_geometrically_subsumed_branches(sys, terms, stats);
        }
    }

    // ========================================================================
    // ITERATIVE GEOMETRIC REDUNDANCY PRUNING (only for intersections)
    // ========================================================================
    #define MAX_SIMPLIFY_PASSES 8
    if (op == ALEA_OP_INTERSECTION && terms->vec.count > 1) {
        for (int pass = 0; pass < MAX_SIMPLIFY_PASSES; pass++) {
            size_t before = terms->vec.count;

            // Step 1: Recompute tight intersection bbox
            alea_bbox_t tight_box;
            if (!compute_intersection_bbox(sys, terms, &tight_box)) {
                if (stats) stats->empty_cells_removed++;
                node_list_destroy(terms);
                return RESULT_EMPTY;
            }

            // Step 2: Check each term for geometric redundancy against bbox
            size_t write = 0;
            for (size_t i = 0; i < terms->vec.count; i++) {
                alea_node_id_t tid = terms->vec.data[i];
                const alea_node_t* tn = &sys->nodes.data[tid];

                if (ALEA_GET_OPERATION(tn) == ALEA_OP_PRIMITIVE) {
                    int result = check_geometric_redundancy(sys, tid, &tight_box);

                    if (result == -1) {
                        if (stats) stats->empty_cells_removed++;
                        node_list_destroy(terms);
                        return RESULT_EMPTY;
                    }

                    if (result == 1) {
                        if (stats) stats->absorption_reductions++;
                        continue;
                    }
                }

                terms->vec.data[write++] = tid;
            }
            terms->vec.count = write;

            // Step 3: Pairwise containment pruning
            remove_subsumed_terms(sys, terms, stats);

            if (terms->vec.count >= before) break;  // Fixed point reached
        }
    }
    #undef MAX_SIMPLIFY_PASSES

    // ========================================================================
    // BUILD RESULT
    // ========================================================================
    if (terms->vec.count == 0) {
        node_list_destroy(terms);
        return (op == ALEA_OP_INTERSECTION) ? RESULT_UNIVERSE : RESULT_EMPTY;
    }

    if (terms->vec.count == 1) {
        alea_node_id_t result = terms->vec.data[0];
        node_list_destroy(terms);
        return result;
    }

    alea_node_id_t result = build_balanced_tree(sys, terms->vec.data, terms->vec.count, op);
    node_list_destroy(terms);

    return result;
}

// ============================================================================
// RECURSIVE SIMPLIFICATION
// ============================================================================

static alea_node_id_t simplify_recursive(
    alea_system_t* sys,
    alea_node_id_t node_id,
    alea_simplify_stats_t* stats
) {
    if (node_id == ALEA_NODE_ID_INVALID || node_id >= alea_vec_count(&sys->nodes)) {
        return ALEA_NODE_ID_INVALID;
    }
    
    const alea_node_t* node = &sys->nodes.data[node_id];
    alea_operation_t op = ALEA_GET_OPERATION(node);
    
    if (op == ALEA_OP_PRIMITIVE) {
        return node_id;
    }
    
    if (op == ALEA_OP_INTERSECTION || op == ALEA_OP_UNION) {
        return flatten_and_optimize(sys, node_id, stats);
    }
    
    if (op == ALEA_OP_COMPLEMENT) {
        alea_node_id_t child = simplify_recursive(sys, node->operation.left, stats);
        if (child == ALEA_NODE_ID_INVALID || child == RESULT_EMPTY) {
            return RESULT_UNIVERSE;
        }
        if (child == RESULT_UNIVERSE) {
            return RESULT_EMPTY;
        }
        return alea_create_complement(sys, child);
    }
    
    if (op == ALEA_OP_DIFFERENCE) {
        alea_node_id_t left = simplify_recursive(sys, node->operation.left, stats);
        alea_node_id_t right = simplify_recursive(sys, node->operation.right, stats);
        
        if (left == ALEA_NODE_ID_INVALID || left == RESULT_EMPTY) {
            return RESULT_EMPTY;
        }
        if (right == ALEA_NODE_ID_INVALID || right == RESULT_EMPTY) {
            return left;
        }
        if (right == RESULT_UNIVERSE) {
            return RESULT_EMPTY;
        }
        if (left == RESULT_UNIVERSE) {
            return alea_create_complement(sys, right);
        }
        
        return alea_create_difference(sys, left, right);
    }
    
    return node_id;
}

// ============================================================================
// PUBLIC API
// ============================================================================

alea_node_id_t alea_tree_simplify(
    alea_system_t* sys,
    alea_node_id_t root_id,
    alea_simplify_stats_t* stats
) {
    if (!sys || root_id == ALEA_NODE_ID_INVALID) {
        return ALEA_NODE_ID_INVALID;
    }
    
    if (stats) memset(stats, 0, sizeof(*stats));
    
    if (!alea_vec_empty(&sys->cell_refs) && !alea_cell_complements_resolved(sys)) {
        alea_resolve_cell_complements(sys);
    }
    
    size_t nodes_before = alea_vec_count(&sys->nodes);
    
    // Step 1: Convert to NNF
    alea_node_id_t nnf = to_nnf(sys, root_id, false, stats);
    if (nnf == ALEA_NODE_ID_INVALID) {
        return ALEA_NODE_ID_INVALID;
    }
    
    // Step 2: Flatten and optimize
    alea_node_id_t result = simplify_recursive(sys, nnf, stats);
    
    if (result == RESULT_EMPTY) {
        result = ALEA_NODE_ID_INVALID;
    } else if (result == RESULT_UNIVERSE) {
        ALEA_LOG_WARN("cell simplified to universe");
        result = nnf;
    }
    
    if (stats) {
        stats->nodes_before = nodes_before;
        stats->nodes_after = alea_vec_count(&sys->nodes);
    }
    
    return result;
}

void alea_simplify_all_cells(
    alea_system_t* sys,
    alea_simplify_stats_t* stats
) {
    if (!sys) return;
    
    alea_simplify_stats_t total = {0};
    
    alea_resolve_cell_complements(sys);
    
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        alea_cell_entry_t* cell = &sys->cells.data[i];
        if (cell->root_node_id == ALEA_NODE_ID_INVALID) continue;
        
        alea_simplify_stats_t cell_stats = {0};
        alea_node_id_t new_root = alea_tree_simplify(sys, cell->root_node_id, &cell_stats);
        
        if (new_root == ALEA_NODE_ID_INVALID) {
            total.empty_cells_removed++;
            ALEA_LOG_INFO("Cell %d simplified to empty set", cell->mc_cell_id);
        }
        
        cell->root_node_id = new_root;
        
        total.complements_eliminated += cell_stats.complements_eliminated;
        total.double_negations += cell_stats.double_negations;
        total.idempotent_reductions += cell_stats.idempotent_reductions;
        total.absorption_reductions += cell_stats.absorption_reductions;
        total.contradictions_found += cell_stats.contradictions_found;
        total.tautologies_found += cell_stats.tautologies_found;
        total.union_branches_absorbed += cell_stats.union_branches_absorbed;
        total.union_common_factors += cell_stats.union_common_factors;
        total.union_branches_subsumed += cell_stats.union_branches_subsumed;
    }
    
    total.nodes_after = alea_vec_count(&sys->nodes);
    
    if (stats) {
        *stats = total;
    }
}

/**
 * @brief Prove a cell is empty using recursive bbox subdivision.
 *
 * Evaluates the full CSG tree with interval arithmetic over the cell's
 * bounding box.  If the result is "entirely outside" the cell is provably
 * empty.  When the initial check is inconclusive the box is bisected along
 * its longest axis and both halves are tested recursively, up to @p depth
 * levels of subdivision.
 */
static bool cell_is_provably_empty(
    const alea_system_t* sys,
    alea_node_id_t root,
    const alea_bbox_t* box,
    int depth
) {
    alea_box_relation_t rel = alea_tree_box_relation(sys, root, box);
    if (rel == ALEA_RELATION_POSITIVE) return true;   /* all outside */
    if (rel == ALEA_RELATION_NEGATIVE) return false;   /* all inside  */
    if (depth <= 0) return false;                      /* inconclusive */

    /* Bisect along longest axis */
    double dx = box->max_x - box->min_x;
    double dy = box->max_y - box->min_y;
    double dz = box->max_z - box->min_z;

    alea_bbox_t a = *box, b = *box;
    if (dx >= dy && dx >= dz) {
        double mid = (box->min_x + box->max_x) * 0.5;
        a.max_x = mid;  b.min_x = mid;
    } else if (dy >= dz) {
        double mid = (box->min_y + box->max_y) * 0.5;
        a.max_y = mid;  b.min_y = mid;
    } else {
        double mid = (box->min_z + box->max_z) * 0.5;
        a.max_z = mid;  b.min_z = mid;
    }

    return cell_is_provably_empty(sys, root, &a, depth - 1) &&
           cell_is_provably_empty(sys, root, &b, depth - 1);
}

void alea_flatten_all_cells(
    alea_system_t* sys,
    alea_simplify_stats_t* stats
) {
    alea_simplify_all_cells(sys, stats);

    /* Second pass: use full-tree interval arithmetic with bbox subdivision
     * to detect empty cells that the per-primitive simplifier missed.
     * Typical case: fill-universe cells entirely outside the parent boundary
     * where the contradiction only becomes apparent through composed intervals. */
    size_t interval_detected = 0;
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        alea_cell_entry_t* cell = &sys->cells.data[i];
        if (cell->root_node_id == ALEA_NODE_ID_INVALID) continue;

        const alea_node_t* root_node = &sys->nodes.data[cell->root_node_id];
        const alea_bbox_t bbox_v = alea_node_bbox_get(&root_node->bbox);
        const alea_bbox_t* bbox = &bbox_v;

        /* Skip cells with degenerate bbox */
        if (bbox->min_x > bbox->max_x) continue;

        if (cell_is_provably_empty(sys, cell->root_node_id, bbox, 4)) {
            ALEA_LOG_INFO("Cell %d is provably empty (interval subdivision)",
                         cell->mc_cell_id);
            cell->root_node_id = ALEA_NODE_ID_INVALID;
            interval_detected++;
        }
    }
    if (stats) stats->empty_cells_removed += interval_detected;

    /* Compact: remove all cells with root_node_id == ALEA_NODE_ID_INVALID.
     * After flatten these are fill-universe cells entirely outside the parent
     * cell boundary — their intersection with the parent is geometrically empty. */
    size_t write = 0;
    size_t removed = 0;
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        alea_cell_entry_t* cell = &sys->cells.data[i];
        if (cell->root_node_id == ALEA_NODE_ID_INVALID) {
            free(cell->surface_indices);
            free(cell->neighbors);
            free(cell->lat_fill);
            free(cell->comments);
            free(cell->inline_comment);
            removed++;
        } else {
            if (write != i) {
                sys->cells.data[write] = sys->cells.data[i];
            }
            write++;
        }
    }
    sys->cells.count = write;
    alea_system_invalidate_query_caches(sys, ALEA_CACHE_ALL);

    if (removed > 0) {
        ALEA_LOG_INFO("Removed %zu empty cells (%zu by simplifier, %zu by interval analysis)",
                     removed, removed - interval_detected, interval_detected);
    }
}

// ============================================================================
// CELL SPLITTING
// ============================================================================

int alea_split_union_cells(alea_system_t* sys) {
    if (!sys) return -1;

    size_t original_count = alea_vec_count(&sys->cells);
    int new_cells_created = 0;

    /* Phase 1: For each cell with a top-level union, create new cells for branches
     * and mark the original for removal. */
    for (size_t i = 0; i < original_count; i++) {
        alea_cell_entry_t* cell = &sys->cells.data[i];
        if (cell->root_node_id == ALEA_NODE_ID_INVALID) continue;

        const alea_node_t* root_node = &sys->nodes.data[cell->root_node_id];
        if (ALEA_GET_OPERATION(root_node) != ALEA_OP_UNION) continue;

        /* Collect union branches */
        node_list_t* branches = node_list_create(16);
        if (!branches) return -1;

        collect_terms(sys, cell->root_node_id, ALEA_OP_UNION, branches);

        if (branches->vec.count <= 1) {
            node_list_destroy(branches);
            continue;
        }

        /* Save cell properties before creating new cells (vector may realloc) */
        int material_index = cell->material_index;
        double density = cell->density;
        unsigned int is_mass_density = cell->is_mass_density;
        int universe_id = cell->universe_id;
        double temperature = cell->temperature;
        unsigned int has_temperature = cell->has_temperature;
        size_t original_idx = i;

        /* Create one new cell per branch */
        for (size_t b = 0; b < branches->vec.count; b++) {
            int idx = alea_add_cell(sys, 0, branches->vec.data[b],
                                   material_index, density, universe_id);
            if (idx >= 0) {
                alea_cell_entry_t* new_cell = &sys->cells.data[idx];
                new_cell->is_mass_density = is_mass_density;
                new_cell->temperature = temperature;
                new_cell->has_temperature = has_temperature;
                /* on_cell_added callback was already fired by alea_add_cell;
                 * notify copy so module can copy params from original cell */
                if (sys->on_cell_copied)
                    sys->on_cell_copied(sys->cell_hook_userdata, (size_t)idx, original_idx);
                new_cells_created++;
            }
        }

        node_list_destroy(branches);

        /* Mark original cell for removal.
         * Re-fetch pointer since alea_add_cell may have reallocated. */
        sys->cells.data[i].root_node_id = ALEA_NODE_ID_INVALID;
    }

    if (new_cells_created == 0) return 0;

    /* Phase 2: Compact — remove marked cells */
    size_t write = 0;
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        alea_cell_entry_t* cell = &sys->cells.data[i];
        if (cell->root_node_id == ALEA_NODE_ID_INVALID) {
            free(cell->surface_indices);
            free(cell->neighbors);
            free(cell->lat_fill);
            free(cell->comments);
            free(cell->inline_comment);
        } else {
            if (write != i) {
                sys->cells.data[write] = sys->cells.data[i];
            }
            write++;
        }
    }
    sys->cells.count = write;
    alea_system_invalidate_query_caches(sys, ALEA_CACHE_ALL);

    ALEA_LOG_INFO("Split union cells: %d new cells created", new_cells_created);
    return new_cells_created;
}

alea_node_id_t alea_tree_to_nnf(
    alea_system_t* sys,
    alea_node_id_t root_id
) {
    if (!sys || root_id == ALEA_NODE_ID_INVALID) {
        return ALEA_NODE_ID_INVALID;
    }
    return to_nnf(sys, root_id, false, NULL);
}

bool alea_tree_is_nnf(
    const alea_system_t* sys,
    alea_node_id_t root_id
) {
    if (!sys || root_id == ALEA_NODE_ID_INVALID || root_id >= alea_vec_count(&sys->nodes)) {
        return true;
    }
    
    const alea_node_t* node = &sys->nodes.data[root_id];
    alea_operation_t op = ALEA_GET_OPERATION(node);
    
    if (op == ALEA_OP_COMPLEMENT) {
        return false;
    }
    
    if (op == ALEA_OP_PRIMITIVE) {
        return true;
    }
    
    return alea_tree_is_nnf(sys, node->operation.left) &&
           alea_tree_is_nnf(sys, node->operation.right);
}

void alea_simplify_stats_print(const alea_simplify_stats_t* stats) {
    if (!stats) return;

    ALEA_LOG_INFO("Simplification statistics:");
    ALEA_LOG_INFO("  Nodes before:              %zu", stats->nodes_before);
    ALEA_LOG_INFO("  Nodes after:               %zu", stats->nodes_after);
    ALEA_LOG_INFO("  Complements eliminated:    %zu", stats->complements_eliminated);
    ALEA_LOG_INFO("  Double negations:          %zu", stats->double_negations);
    ALEA_LOG_INFO("  Idempotent reductions:     %zu", stats->idempotent_reductions);
    ALEA_LOG_INFO("  Redundant terms removed:   %zu", stats->absorption_reductions);
    ALEA_LOG_INFO("  Contradictions (A∩¬A):     %zu", stats->contradictions_found);
    ALEA_LOG_INFO("  Tautologies (A∪¬A):        %zu", stats->tautologies_found);
    ALEA_LOG_INFO("  Union branches absorbed:   %zu", stats->union_branches_absorbed);
    ALEA_LOG_INFO("  Union common factors:      %zu", stats->union_common_factors);
    ALEA_LOG_INFO("  Union branches subsumed:   %zu", stats->union_branches_subsumed);
    ALEA_LOG_INFO("  Empty cells:               %zu", stats->empty_cells_removed);
}

// ============================================================================
// DE MORGAN COMPLEMENT (creates new tree without COMPLEMENT nodes)
// ============================================================================

/**
 * @brief Apply De Morgan's law to create the complement of a tree
 *
 * Instead of wrapping in a COMPLEMENT node, this recursively transforms:
 *   #(A ∩ B) → #A ∪ #B
 *   #(A ∪ B) → #A ∩ #B
 *   #(primitive with sense S) → new primitive node with sense -S
 *
 * This produces a flattened tree without COMPLEMENT nodes.
 * Works on any tree - does NOT require primitives to have MCNP surface IDs.
 *
 * @param sys CSG system
 * @param node_id Node to complement
 * @return New node ID representing the complement, or ALEA_NODE_ID_INVALID on error
 */
alea_node_id_t alea_apply_demorgan(alea_system_t* sys, alea_node_id_t node_id) {
    if (!sys || node_id == ALEA_NODE_ID_INVALID || node_id >= alea_vec_count(&sys->nodes)) {
        return ALEA_NODE_ID_INVALID;
    }

    const alea_node_t* node = &sys->nodes.data[node_id];
    alea_operation_t op = ALEA_GET_OPERATION(node);

    if (op == ALEA_OP_PRIMITIVE) {
        /* Create new node with flipped sense */
        const alea_node_t* src = &sys->nodes.data[node_id];
        int8_t flipped_sense = -src->primitive.sense;
        return alea_add_primitive_node(sys, src->primitive.primitive_id, flipped_sense,
                                      src->primitive.inverted, 0);
    }

    if (op == ALEA_OP_INTERSECTION) {
        /* #(A ∩ B) → #A ∪ #B */
        alea_node_id_t left = node->operation.left;
        alea_node_id_t right = node->operation.right;

        alea_node_id_t new_left = alea_apply_demorgan(sys, left);
        alea_node_id_t new_right = alea_apply_demorgan(sys, right);

        if (new_left == ALEA_NODE_ID_INVALID || new_right == ALEA_NODE_ID_INVALID) {
            return ALEA_NODE_ID_INVALID;
        }

        return alea_create_union(sys, new_left, new_right);
    }

    if (op == ALEA_OP_UNION) {
        /* #(A ∪ B) → #A ∩ #B */
        alea_node_id_t left = node->operation.left;
        alea_node_id_t right = node->operation.right;

        alea_node_id_t new_left = alea_apply_demorgan(sys, left);
        alea_node_id_t new_right = alea_apply_demorgan(sys, right);

        if (new_left == ALEA_NODE_ID_INVALID || new_right == ALEA_NODE_ID_INVALID) {
            return ALEA_NODE_ID_INVALID;
        }

        return alea_create_intersection(sys, new_left, new_right);
    }

    if (op == ALEA_OP_COMPLEMENT) {
        /* ##A → A (double negation elimination) */
        return node->operation.left;
    }

    if (op == ALEA_OP_DIFFERENCE) {
        /* #(A - B) = #(A ∩ #B) = #A ∪ B */
        alea_node_id_t left = node->operation.left;
        alea_node_id_t right = node->operation.right;

        alea_node_id_t new_left = alea_apply_demorgan(sys, left);
        /* right is already negated in difference, so we just use it as-is */

        if (new_left == ALEA_NODE_ID_INVALID) {
            return ALEA_NODE_ID_INVALID;
        }

        return alea_create_union(sys, new_left, right);
    }

    return ALEA_NODE_ID_INVALID;
}
