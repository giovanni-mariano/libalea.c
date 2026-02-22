// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_eval.c
 * @brief CSG tree evaluation implementation
 * 
 * Evaluates CSG trees using:
 * - Recursive tree traversal
 * - Early termination for boolean operations
 * - Bounding box culling
 * - Primitive evaluation dispatch
 */

#include "alea_eval.h"
#include "core/alea_system.h"
#include "primitives/primitive_eval.h"
#include "primitives/bbox.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "util/alea_log.h"

// ============================================================================
// PERFORMANCE COUNTERS (THREAD-LOCAL)
// ============================================================================

// Thread-local performance counters for safe parallel execution
static _Thread_local alea_perf_counters_t g_alea_perf = {0};

// ============================================================================
// HELPER: EVALUATE PRIMITIVE NODE
// ============================================================================

/* Debug trace flag - linked from alea_universe.c */
extern _Thread_local int g_debug_point_trace;
static _Thread_local int g_trace_cell_id = 0;  /* Cell ID to trace primitives for */

static inline double eval_primitive_node(
    const alea_system_t* sys,
    const alea_node_t* node,
    double x, double y, double z
) {
    g_alea_perf.primitive_evaluations++;

    // Get primitive data directly from system array
    if (node->primitive.primitive_id >= alea_vec_count(&sys->primitives)) {
        ALEA_LOG_ERROR("Invalid primitive ID %u",
                node->primitive.primitive_id);
        return 1.0; // Outside if invalid
    }

    const alea_primitive_entry_t* prim = &sys->primitives.data[node->primitive.primitive_id];

    // Evaluate the primitive
    double distance = alea_primitive_eval(
        prim->type,
        &prim->data,
        x, y, z
    );

    // Apply inverted flag (primitive was normalized/canonicalized)
    // and sense (orientation): XOR them since both flip the sign
    bool flip = (node->primitive.sense > 0) != (node->primitive.inverted != 0);
    double final_dist = flip ? -distance : distance;

    // Debug trace for specific cell
    if (g_debug_point_trace && g_trace_cell_id != 0 && final_dist <= 0) {
        ALEA_LOG_DEBUG("surf=%d prim=%u type=%d sense=%d inverted=%d raw=%.4f final=%.4f",
               node->primitive.mcnp_surface_id, node->primitive.primitive_id,
               prim->type, node->primitive.sense, node->primitive.inverted,
               distance, final_dist);
    }

    return final_dist;
}

// ============================================================================
// POINT EVALUATION (RECURSIVE)
// ============================================================================

static _Thread_local int g_eval_trace_depth = 0;

double alea_evaluate_point(
    const alea_system_t* sys,
    alea_node_id_t node_id,
    double x, double y, double z
) {
    if (!sys || node_id == ALEA_NODE_ID_INVALID || node_id >= alea_vec_count(&sys->nodes)) {
        return 1.0; // Outside
    }

    // Get node directly from system array
    const alea_node_t* node = &sys->nodes.data[node_id];


    // Early bounding box test
    if (!alea_bbox_contains_point(&node->bbox, x, y, z)) {
        //g_alea_perf.bbox_tests++;
        if (g_debug_point_trace && g_eval_trace_depth == 0) {
            ALEA_LOG_DEBUG("node %u: bbox reject (bbox=[%.2f,%.2f]x[%.2f,%.2f]x[%.2f,%.2f])",
                   node_id, node->bbox.min_x, node->bbox.max_x,
                   node->bbox.min_y, node->bbox.max_y,
                   node->bbox.min_z, node->bbox.max_z);
        }
        return 1.0;
    }
    
    alea_operation_t op = ALEA_GET_OPERATION(node);
    
    // Base case: primitive
    if (op == ALEA_OP_PRIMITIVE) {
        return eval_primitive_node(sys, node, x, y, z);
    }
    
    // Boolean operations
    g_alea_perf.boolean_operations++;

    // Handle unary complement first (only uses left child)
    if (op == ALEA_OP_COMPLEMENT) {
        double child_val = alea_evaluate_point(sys, node->operation.left, x, y, z);
        return -child_val;  // Complement: negate the signed distance
    }
    
    // Evaluate left child
    double left_val = alea_evaluate_point(sys, node->operation.left, x, y, z);
    
    // Early termination optimizations
    switch (op) {
        case ALEA_OP_UNION:
            if (left_val <= 0.0) return left_val;
            break;
        case ALEA_OP_INTERSECTION:
            if (left_val > 0.0) return left_val;
            break;
        case ALEA_OP_DIFFERENCE:
            if (left_val > 0.0) return left_val;
            break;
        default:
            break;
    }
    
    // Evaluate right child
    double right_val = alea_evaluate_point(sys, node->operation.right, x, y, z);
    
    // Combine using boolean operation
    switch (op) {
        case ALEA_OP_UNION:
            return (left_val < right_val) ? left_val : right_val;
        case ALEA_OP_INTERSECTION:
            return (left_val > right_val) ? left_val : right_val;
        case ALEA_OP_DIFFERENCE:
            return (left_val > -right_val) ? left_val : -right_val;
        default:
            return 1.0;
    }
}

bool alea_contains_point(
    const alea_system_t* sys,
    alea_node_id_t node_id,
    double x, double y, double z
) {
    return alea_evaluate_point(sys, node_id, x, y, z) <= 0.0;
}

// ============================================================================
// MATERIAL IDENTIFICATION
// ============================================================================

alea_material_id_t alea_get_material_at(
    const alea_system_t* sys,
    alea_node_id_t node_id,
    double x, double y, double z
) {
    if (!sys || !alea_contains_point(sys, node_id, x, y, z)) {
        return ALEA_MATERIAL_NONE;
    }
    
    if (node_id >= alea_vec_count(&sys->nodes)) {
        return ALEA_MATERIAL_NONE;
    }
    
    const alea_node_t* node = &sys->nodes.data[node_id];
    
    alea_operation_t op = ALEA_GET_OPERATION(node);
    
    // Base case: primitive has material
    if (op == ALEA_OP_PRIMITIVE) {
        return node->material_id;
    }
    
    // Boolean operations
    switch (op) {

        case ALEA_OP_COMPLEMENT: {
            // Complement region typically has no material (void/outside)
            // If you're inside the complement, you're outside the original
            return ALEA_MATERIAL_NONE;
        }
        
        case ALEA_OP_UNION: {
            alea_material_id_t left_mat = alea_get_material_at(
                sys, node->operation.left, x, y, z
            );
            if (left_mat != ALEA_MATERIAL_NONE) return left_mat;
            return alea_get_material_at(sys, node->operation.right, x, y, z);
        }
        
        case ALEA_OP_INTERSECTION: {
            alea_material_id_t left_mat = alea_get_material_at(
                sys, node->operation.left, x, y, z
            );
            alea_material_id_t right_mat = alea_get_material_at(
                sys, node->operation.right, x, y, z
            );
            
            if (left_mat != ALEA_MATERIAL_NONE && right_mat != ALEA_MATERIAL_NONE) {
                return left_mat;
            }
            return ALEA_MATERIAL_NONE;
        }
        
        case ALEA_OP_DIFFERENCE: {
            alea_material_id_t left_mat = alea_get_material_at(
                sys, node->operation.left, x, y, z
            );
            if (left_mat == ALEA_MATERIAL_NONE) {
                return ALEA_MATERIAL_NONE;
            }
            
            if (alea_contains_point(sys, node->operation.right, x, y, z)) {
                return ALEA_MATERIAL_NONE;
            }
            
            return left_mat;
        }
        
        default:
            return ALEA_MATERIAL_NONE;
    }
}

// ============================================================================
// BATCH EVALUATION
// ============================================================================

void alea_evaluate_points_batch(
    const alea_system_t* sys,
    alea_node_id_t node_id,
    alea_sample_point_t* points,
    size_t count
) {
    for (size_t i = 0; i < count; i++) {
        double eval = alea_evaluate_point(
            sys, node_id,
            points[i].x, points[i].y, points[i].z
        );
        
        points[i].inside = (eval <= 0.0);
        
        if (points[i].inside) {
            points[i].material = alea_get_material_at(
                sys, node_id,
                points[i].x, points[i].y, points[i].z
            );
        } else {
            points[i].material = ALEA_MATERIAL_NONE;
        }
    }
}

// ============================================================================
// PERFORMANCE MONITORING
// ============================================================================

void alea_perf_reset(void) {
    memset(&g_alea_perf, 0, sizeof(g_alea_perf));
}

alea_perf_counters_t alea_perf_get(void) {
    return g_alea_perf;
}

void alea_perf_print(void) {
    ALEA_LOG_INFO("CSG Performance Counters:");
    ALEA_LOG_INFO("  Primitive evaluations: %zu", g_alea_perf.primitive_evaluations);
    ALEA_LOG_INFO("  Boolean operations: %zu", g_alea_perf.boolean_operations);
    ALEA_LOG_INFO("  Bounding box tests: %zu", g_alea_perf.bbox_tests);
    ALEA_LOG_INFO("  Cache hits: %zu", g_alea_perf.cache_hits);
    ALEA_LOG_INFO("  Cache misses: %zu", g_alea_perf.cache_misses);
}

// ============================================================================
// INTERVAL EVALUATION (EXPERIMENTAL)
// ============================================================================
//
// Define ALEA_USE_INTERVAL_CULLING=1 to enable interval-based culling.
// This is experimental - remove this section if it doesn't improve performance.
//

#include "primitives/primitive_desc.h"

/**
 * Interval operations for CSG boolean combinations:
 * - UNION: min(f1, f2)      -> [min(a.min, b.min), min(a.max, b.max)]
 * - INTERSECTION: max(f1, f2) -> [max(a.min, b.min), max(a.max, b.max)]
 * - COMPLEMENT: -f          -> [-a.max, -a.min]
 */
static inline alea_interval_t iv_min(alea_interval_t a, alea_interval_t b) {
    return alea_iv_make(
        (a.min < b.min) ? a.min : b.min,
        (a.max < b.max) ? a.max : b.max
    );
}

static inline alea_interval_t iv_max(alea_interval_t a, alea_interval_t b) {
    return alea_iv_make(
        (a.min > b.min) ? a.min : b.min,
        (a.max > b.max) ? a.max : b.max
    );
}

static inline alea_interval_t iv_negate(alea_interval_t a) {
    return alea_iv_make(-a.max, -a.min);
}

/* Evaluate a primitive node with interval arithmetic */
static alea_interval_t eval_primitive_interval(
    const alea_system_t* sys,
    const alea_node_t* node,
    const alea_bbox_t* box
) {
    if (node->primitive.primitive_id >= alea_vec_count(&sys->primitives)) {
        return alea_iv_make(-1e30, 1e30);  /* Unknown - assume anything */
    }

    const alea_primitive_entry_t* prim = &sys->primitives.data[node->primitive.primitive_id];

    alea_interval_t iv = alea_primitive_interval_eval(prim->type, &prim->data, box);

    /* Apply inverted flag and sense (same logic as point evaluation) */
    bool flip = (node->primitive.sense > 0) != (node->primitive.inverted != 0);
    if (flip) {
        iv = iv_negate(iv);
    }

    return iv;
}

alea_interval_t alea_evaluate_interval(
    const alea_system_t* sys,
    alea_node_id_t node_id,
    const alea_bbox_t* box
) {
    /* Default: assume anything is possible */
    alea_interval_t unknown = alea_iv_make(-1e30, 1e30);

    if (!sys || node_id == ALEA_NODE_ID_INVALID || node_id >= alea_vec_count(&sys->nodes)) {
        return alea_iv_make(1.0, 1e30);  /* Outside */
    }

    const alea_node_t* node = &sys->nodes.data[node_id];

    /* Early bbox test: if query box doesn't intersect node bbox, definitely outside */
    if (box->max_x < node->bbox.min_x || box->min_x > node->bbox.max_x ||
        box->max_y < node->bbox.min_y || box->min_y > node->bbox.max_y ||
        box->max_z < node->bbox.min_z || box->min_z > node->bbox.max_z) {
        return alea_iv_make(1.0, 1e30);  /* Positive = outside */
    }

    alea_operation_t op = ALEA_GET_OPERATION(node);

    /* Base case: primitive */
    if (op == ALEA_OP_PRIMITIVE) {
        return eval_primitive_interval(sys, node, box);
    }

    /* Unary complement */
    if (op == ALEA_OP_COMPLEMENT) {
        alea_interval_t child = alea_evaluate_interval(sys, node->operation.left, box);
        return iv_negate(child);
    }

    /* Binary operations */
    alea_interval_t left = alea_evaluate_interval(sys, node->operation.left, box);
    alea_interval_t right = alea_evaluate_interval(sys, node->operation.right, box);

    switch (op) {
        case ALEA_OP_UNION:
            return iv_min(left, right);

        case ALEA_OP_INTERSECTION:
            return iv_max(left, right);

        case ALEA_OP_DIFFERENCE:
            /* difference(A, B) = intersection(A, complement(B)) */
            return iv_max(left, iv_negate(right));

        default:
            return unknown;
    }
}

alea_box_relation_t alea_tree_box_relation(
    const alea_system_t* sys,
    alea_node_id_t node_id,
    const alea_bbox_t* box
) {
    alea_interval_t iv = alea_evaluate_interval(sys, node_id, box);

    if (iv.max < 0.0) {
        return ALEA_RELATION_NEGATIVE;  /* Entirely inside */
    } else if (iv.min > 0.0) {
        return ALEA_RELATION_POSITIVE;  /* Entirely outside */
    }
    return ALEA_RELATION_INTERSECT;     /* Crosses surface */
}