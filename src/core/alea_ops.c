// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_ops.c
 * @brief CSG boolean operations implementation
 * 
 * Creates CSG tree nodes for boolean operations with
 * appropriate bounding box computation.
 */

#include "alea_ops.h"
#include "../primitives/bbox.h"
#include <stddef.h>

// ============================================================================
// BINARY OPERATIONS
// ============================================================================

alea_node_id_t alea_create_union(
    alea_system_t* sys,
    alea_node_id_t left,
    alea_node_id_t right
) {
    // Validate inputs
    if (left == ALEA_NODE_ID_INVALID || right == ALEA_NODE_ID_INVALID) {
        return ALEA_NODE_ID_INVALID;
    }
    
    // Allocate new node
    alea_node_id_t node_id = alea_alloc_node(sys);
    if (node_id == ALEA_NODE_ID_INVALID) {
        return ALEA_NODE_ID_INVALID;
    }

    // Get node pointer for setup
    alea_node_t* node = alea_get_node(sys, node_id);

    // Set operation type
    ALEA_SET_OPERATION(node, ALEA_OP_UNION);
    
    // Operations don't have materials (children do)
    node->material_id = ALEA_MATERIAL_NONE;
    
    // Set operands
    node->operation.left = left;
    node->operation.right = right;
    
    // Compute bounding box: union of child boxes
    const alea_node_t* left_node = alea_get_node(sys, left);
    const alea_node_t* right_node = alea_get_node(sys, right);
    
    node->bbox = alea_bbox_union(&left_node->bbox, &right_node->bbox);
    
    return node_id;
}

alea_node_id_t alea_create_intersection(
    alea_system_t* sys,
    alea_node_id_t left,
    alea_node_id_t right
) {
    // Validate inputs
    if (left == ALEA_NODE_ID_INVALID || right == ALEA_NODE_ID_INVALID) {
        return ALEA_NODE_ID_INVALID;
    }
    
    // Allocate new node
    alea_node_id_t node_id = alea_alloc_node(sys);
    if (node_id == ALEA_NODE_ID_INVALID) {
        return ALEA_NODE_ID_INVALID;
    }

    // Get node pointer for setup
    alea_node_t* node = alea_get_node(sys, node_id);

    // Set operation type
    ALEA_SET_OPERATION(node, ALEA_OP_INTERSECTION);
    
    // Operations don't have materials
    node->material_id = ALEA_MATERIAL_NONE;
    
    // Set operands
    node->operation.left = left;
    node->operation.right = right;
    
    // Compute bounding box: intersection of child boxes
    const alea_node_t* left_node = alea_get_node(sys, left);
    const alea_node_t* right_node = alea_get_node(sys, right);
    
    node->bbox = alea_bbox_intersection(&left_node->bbox, &right_node->bbox);
    
    return node_id;
}

alea_node_id_t alea_create_difference(
    alea_system_t* sys,
    alea_node_id_t left,
    alea_node_id_t right
) {
    // Validate inputs
    if (left == ALEA_NODE_ID_INVALID || right == ALEA_NODE_ID_INVALID) {
        return ALEA_NODE_ID_INVALID;
    }
    
    // Allocate new node
    alea_node_id_t node_id = alea_alloc_node(sys);
    if (node_id == ALEA_NODE_ID_INVALID) {
        return ALEA_NODE_ID_INVALID;
    }

    // Get node pointer for setup
    alea_node_t* node = alea_get_node(sys, node_id);

    // Set operation type
    ALEA_SET_OPERATION(node, ALEA_OP_DIFFERENCE);
    
    // Operations don't have materials
    node->material_id = ALEA_MATERIAL_NONE;
    
    // Set operands
    node->operation.left = left;
    node->operation.right = right;
    
    // Compute bounding box: conservative - use left bbox
    // (actual result may be smaller, but this is safe)
    const alea_node_t* left_node = alea_get_node(sys, left);
    node->bbox = left_node->bbox;
    
    return node_id;
}

// ============================================================================
// N-ARY OPERATIONS
// ============================================================================

/* Recursive helper: build balanced binary tree of unions */
static alea_node_id_t union_balanced(
    alea_system_t* sys,
    const alea_node_id_t* ids,
    size_t lo,
    size_t hi
) {
    if (hi - lo == 1) return ids[lo];
    size_t mid = lo + (hi - lo) / 2;
    alea_node_id_t left = union_balanced(sys, ids, lo, mid);
    if (left == ALEA_NODE_ID_INVALID) return ALEA_NODE_ID_INVALID;
    alea_node_id_t right = union_balanced(sys, ids, mid, hi);
    if (right == ALEA_NODE_ID_INVALID) return ALEA_NODE_ID_INVALID;
    return alea_create_union(sys, left, right);
}

alea_node_id_t alea_create_union_many(
    alea_system_t* sys,
    const alea_node_id_t* node_ids,
    size_t count
) {
    if (count == 0 || !node_ids) {
        return ALEA_NODE_ID_INVALID;
    }

    if (count == 1) {
        return node_ids[0];
    }

    // Build balanced binary tree
    // For [A, B, C, D]: (A ∪ B) ∪ (C ∪ D)
    return union_balanced(sys, node_ids, 0, count);
}

/* Recursive helper: build balanced binary tree of intersections */
static alea_node_id_t intersection_balanced(
    alea_system_t* sys,
    const alea_node_id_t* ids,
    size_t lo,
    size_t hi
) {
    if (hi - lo == 1) return ids[lo];
    size_t mid = lo + (hi - lo) / 2;
    alea_node_id_t left = intersection_balanced(sys, ids, lo, mid);
    if (left == ALEA_NODE_ID_INVALID) return ALEA_NODE_ID_INVALID;
    alea_node_id_t right = intersection_balanced(sys, ids, mid, hi);
    if (right == ALEA_NODE_ID_INVALID) return ALEA_NODE_ID_INVALID;
    return alea_create_intersection(sys, left, right);
}

alea_node_id_t alea_create_intersection_many(
    alea_system_t* sys,
    const alea_node_id_t* node_ids,
    size_t count
) {
    if (count == 0 || !node_ids) {
        return ALEA_NODE_ID_INVALID;
    }

    if (count == 1) {
        return node_ids[0];
    }

    // Build balanced binary tree
    return intersection_balanced(sys, node_ids, 0, count);
}

alea_node_id_t alea_create_symmetric_difference(
    alea_system_t* sys,
    alea_node_id_t left,
    alea_node_id_t right
) {
    // Symmetric difference: (A - B) ∪ (B - A)
    
    // Create A - B
    alea_node_id_t a_minus_b = alea_create_difference(sys, left, right);
    if (a_minus_b == ALEA_NODE_ID_INVALID) {
        return ALEA_NODE_ID_INVALID;
    }
    
    // Create B - A
    alea_node_id_t b_minus_a = alea_create_difference(sys, right, left);
    if (b_minus_a == ALEA_NODE_ID_INVALID) {
        return ALEA_NODE_ID_INVALID;
    }
    
    // Union them
    return alea_create_union(sys, a_minus_b, b_minus_a);
}

// ============================================================================
// UNARY OPERATIONS
// ============================================================================


alea_node_id_t alea_create_complement(
    alea_system_t* sys,
    alea_node_id_t operand
) {
    // Validate input
    if (operand == ALEA_NODE_ID_INVALID) {
        return ALEA_NODE_ID_INVALID;
    }
    
    // Allocate new node
    alea_node_id_t node_id = alea_alloc_node(sys);
    if (node_id == ALEA_NODE_ID_INVALID) {
        return ALEA_NODE_ID_INVALID;
    }

    // Get node pointer for setup
    alea_node_t* node = alea_get_node(sys, node_id);

    // Set operation type
    ALEA_SET_OPERATION(node, ALEA_OP_COMPLEMENT);
    
    // Complement doesn't have material (it's "everything else")
    node->material_id = ALEA_MATERIAL_NONE;
    
    // Store operand in left child, right is unused
    node->operation.left = operand;
    node->operation.right = ALEA_NODE_ID_INVALID;
    
    // Bounding box: complement of finite region is infinite
    // Use large conservative bounds
    node->bbox = (alea_bbox_t){
        -1.0e6, 1.0e6,
        -1.0e6, 1.0e6,
        -1.0e6, 1.0e6
    };
    
    return node_id;
}