// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_OPS_H
#define ALEA_OPS_H

/**
 * @file alea_ops.h
 * @brief CSG boolean operations
 * 
 * This module provides creation of CSG boolean operation nodes:
 * - Union (OR)
 * - Intersection (AND)
 * - Difference (subtraction)
 * 
 * Each operation creates a new node in the CSG tree and
 * computes the appropriate bounding box.
 */

#include "core/alea_system.h"

#ifdef __cplusplus
extern "C" {
#endif

// ============================================================================
// BOOLEAN OPERATIONS
// ============================================================================

/**
 * @brief Create a union node (A ∪ B)
 * 
 * The result contains points that are in either A or B (or both).
 * Signed distance: min(dist_A, dist_B)
 * 
 * Bounding box: union of child bounding boxes
 * Material priority: left, then right
 * 
 * @param left Left operand node ID
 * @param right Right operand node ID
 * @return New union node ID, or ALEA_NODE_ID_INVALID on error
 */
alea_node_id_t alea_create_union(
    alea_system_t* sys,
    alea_node_id_t left,
    alea_node_id_t right
);

/**
 * @brief Create an intersection node (A ∩ B)
 * 
 * The result contains only points that are in both A and B.
 * Signed distance: max(dist_A, dist_B)
 * 
 * Bounding box: intersection of child bounding boxes
 * Material: left operand's material
 * 
 * @param left Left operand node ID
 * @param right Right operand node ID
 * @return New intersection node ID, or ALEA_NODE_ID_INVALID on error
 */
alea_node_id_t alea_create_intersection(
    alea_system_t* sys,
    alea_node_id_t left,
    alea_node_id_t right
);

/**
 * @brief Create a difference node (A - B)
 * 
 * The result contains points that are in A but not in B.
 * Signed distance: max(dist_A, -dist_B)
 * 
 * Bounding box: left operand's bounding box (conservative)
 * Material: left operand's material (where not subtracted)
 * 
 * @param left Left operand node ID (what to keep)
 * @param right Right operand node ID (what to subtract)
 * @return New difference node ID, or ALEA_NODE_ID_INVALID on error
 */
alea_node_id_t alea_create_difference(
    alea_system_t* sys,
    alea_node_id_t left,
    alea_node_id_t right
);

// ============================================================================
// CONVENIENCE FUNCTIONS
// ============================================================================

/**
 * @brief Create n-ary union (A ∪ B ∪ C ∪ ...)
 * 
 * Builds a left-balanced tree of union operations.
 * More efficient than manually creating nested unions.
 * 
 * @param node_ids Array of node IDs to union
 * @param count Number of nodes
 * @return Root of union tree, or ALEA_NODE_ID_INVALID on error
 */
alea_node_id_t alea_create_union_many(
    alea_system_t* sys,
    const alea_node_id_t* node_ids,
    size_t count
);

/**
 * @brief Create n-ary intersection (A ∩ B ∩ C ∩ ...)
 * 
 * Builds a left-balanced tree of intersection operations.
 * 
 * @param node_ids Array of node IDs to intersect
 * @param count Number of nodes
 * @return Root of intersection tree, or ALEA_NODE_ID_INVALID on error
 */
alea_node_id_t alea_create_intersection_many(
    alea_system_t* sys,
    const alea_node_id_t* node_ids,
    size_t count
);

/**
 * @brief Create symmetric difference (A ⊕ B = (A - B) ∪ (B - A))
 * 
 * The result contains points that are in A or B but not both.
 * Useful for finding non-overlapping regions.
 * 
 * @param left Left operand node ID
 * @param right Right operand node ID
 * @return New symmetric difference node ID, or ALEA_NODE_ID_INVALID on error
 */
alea_node_id_t alea_create_symmetric_difference(
    alea_system_t* sys,
    alea_node_id_t left,
    alea_node_id_t right
);

/**
 * @brief Create complement node (NOT operation)
 * 
 * The complement of a region A is everything outside A.
 * 
 * @param operand Operand node ID
 * @return New complement node ID, or ALEA_NODE_ID_INVALID on error
 */
alea_node_id_t alea_create_complement(
    alea_system_t* sys,
    alea_node_id_t operand
);

#ifdef __cplusplus
}
#endif

#endif // ALEA_OPS_H