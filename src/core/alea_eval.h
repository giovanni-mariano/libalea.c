// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_EVAL_H
#define ALEA_EVAL_H

/**
 * @file alea_eval.h
 * @brief CSG tree evaluation functions
 * 
 * This module provides functions for evaluating CSG trees at points,
 * including signed distance evaluation, containment tests, and material
 * queries.
 * 
 * The CSG tree is evaluated recursively using boolean operations:
 * - Union: min(left, right)
 * - Intersection: max(left, right)
 * - Difference: max(left, -right)
 */

#include <stddef.h>
#include <stdbool.h>
#include "alea_types.h"


// Forward declaration
typedef struct alea_system alea_system_t;

// Sample point structure (defined in alea_types.h or here if not)
#ifndef ALEA_SAMPLE_POINT_T_DEFINED
#define ALEA_SAMPLE_POINT_T_DEFINED
typedef struct {
    double x, y, z;           // Point coordinates
    bool inside;              // Result: inside solid?
    alea_material_id_t material; // Material at this point
} alea_sample_point_t;
#endif

// Performance counters (defined in alea_types.h or here if not)
#ifndef ALEA_PERF_COUNTERS_T_DEFINED
#define ALEA_PERF_COUNTERS_T_DEFINED
typedef struct {
    size_t primitive_evaluations;  // Number of primitive evals
    size_t boolean_operations;      // Number of boolean ops
    size_t bbox_tests;              // Number of bbox culls
    size_t cache_hits;              // Cache hits (if caching enabled)
    size_t cache_misses;            // Cache misses
} alea_perf_counters_t;
#endif

// ============================================================================
// POINT EVALUATION
// ============================================================================

/**
 * @brief Evaluate signed distance at a point
 * 
 * Returns the signed distance from the point to the surface:
 * - Negative: inside solid
 * - Zero: on surface
 * - Positive: outside solid
 * 
 * Uses early termination and bounding box culling for performance.
 * 
 * @param sys CSG system containing the geometry
 * @param node_id Root node of CSG tree to evaluate
 * @param x X coordinate
 * @param y Y coordinate
 * @param z Z coordinate
 * @return Signed distance value
 */
double alea_evaluate_point(
    const alea_system_t* sys,
    alea_node_id_t node_id,
    double x, double y, double z
);

/**
 * @brief Test if a point is inside a solid
 * 
 * Convenience function equivalent to: alea_evaluate_point() <= 0.0
 * 
 * @param sys CSG system containing the geometry
 * @param node_id Root node of CSG tree
 * @param x X coordinate
 * @param y Y coordinate
 * @param z Z coordinate
 * @return true if point is inside or on surface
 */
bool alea_contains_point(
    const alea_system_t* sys,
    alea_node_id_t node_id,
    double x, double y, double z
);

/**
 * @brief Get material ID at a point
 * 
 * Returns the material at the given point, or ALEA_MATERIAL_NONE
 * if the point is outside all geometry.
 * 
 * For boolean operations:
 * - Union: returns first material found
 * - Intersection: returns left operand's material
 * - Difference: returns left material if not subtracted
 * 
 * @param sys CSG system containing the geometry
 * @param node_id Root node of CSG tree
 * @param x X coordinate
 * @param y Y coordinate
 * @param z Z coordinate
 * @return Material ID or ALEA_MATERIAL_NONE
 */
alea_material_id_t alea_get_material_at(
    const alea_system_t* sys,
    alea_node_id_t node_id,
    double x, double y, double z
);

// ============================================================================
// BATCH EVALUATION
// ============================================================================

/**
 * @brief Evaluate multiple points in batch
 * 
 * More efficient than individual calls for large point sets.
 * Updates the 'inside' and 'material' fields of each sample point.
 * 
 * Typical use case: Monte Carlo volume estimation
 * 
 * @param sys CSG system containing the geometry
 * @param node_id Root node of CSG tree
 * @param points Array of sample points (x,y,z fields must be set)
 * @param count Number of points
 */
void alea_evaluate_points_batch(
    const alea_system_t* sys,
    alea_node_id_t node_id,
    alea_sample_point_t* points,
    size_t count
);

// ============================================================================
// INTERVAL EVALUATION (EXPERIMENTAL)
// ============================================================================

/*
 * Interval-based CSG tree evaluation for region culling.
 * Enable with ALEA_USE_INTERVAL_CULLING=1 at compile time.
 *
 * This evaluates signed distance bounds over a bounding box region,
 * allowing quick rejection of regions that are entirely inside or outside.
 */

#include "primitives/primitive_desc.h"  /* For alea_interval_t, alea_box_relation_t */

/**
 * @brief Evaluate CSG tree with interval arithmetic
 *
 * Returns an interval [min, max] bounding the signed distance function
 * over the given bounding box region.
 *
 * @param sys CSG system
 * @param node_id Root of CSG tree
 * @param box Query bounding box
 * @return Interval bounding the SDF over the box
 */
alea_interval_t alea_evaluate_interval(
    const alea_system_t* sys,
    alea_node_id_t node_id,
    const alea_bbox_t* box
);

/**
 * @brief Classify a bounding box against a CSG tree
 *
 * @param sys CSG system
 * @param node_id Root of CSG tree
 * @param box Query bounding box
 * @return ALEA_RELATION_NEGATIVE (all inside), ALEA_RELATION_POSITIVE (all outside),
 *         or ALEA_RELATION_INTERSECT (crosses surface)
 */
alea_box_relation_t alea_tree_box_relation(
    const alea_system_t* sys,
    alea_node_id_t node_id,
    const alea_bbox_t* box
);

// ============================================================================
// PERFORMANCE MONITORING
// ============================================================================

/**
 * @brief Reset performance counters
 * 
 * Call this before a benchmark to start with clean counters.
 */
void alea_perf_reset(void);

/**
 * @brief Get current performance counters
 * 
 * @return Copy of current counters
 */
alea_perf_counters_t alea_perf_get(void);

/**
 * @brief Print performance statistics to stdout
 * 
 * Useful for debugging and optimization.
 */
void alea_perf_print(void);


#endif // ALEA_EVAL_H