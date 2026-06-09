// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_cell_complement.c
 * @brief Cell complement (#cell) resolution
 * 
 * MCNP allows referencing cells in geometry expressions:
 *   3 0  1 #2 #7    means: inside surface 1 AND NOT(cell 2) AND NOT(cell 7)
 * 
 * The challenge is that when parsing cell 3, cells 2 and 7 may not
 * have been converted yet. 
 * 
 * Solution: Two-pass approach
 *   Pass 1: Store unresolved cell references during parsing
 *   Pass 2: After all cells converted, resolve references
 * 
 * Integration:
 *   1. Add alea_cell_ref_t array to alea_system_t
 *   2. During parsing, when #cell_id is seen, create placeholder node
 *      and record it in the cell_refs array
 *   3. After mcnp_convert_file() converts all cells, call
 *      alea_resolve_cell_complements() to fix up the references
 */
#include "core/alea_cell_complement.h"
#include "core/alea_system.h"
#include "core/alea_ops.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "util/alea_log.h"

/* ============================================================================
 * CELL REFERENCE TRACKING
 * ============================================================================ */


/**
 * @brief Storage for unresolved cell references
 * 
 * Add this to alea_system_t:
 *   alea_cell_ref_t* cell_refs;
 *   size_t cell_ref_count;
 *   size_t cell_ref_capacity;
 */

/* ============================================================================
 * DURING PARSING: Record cell complement reference
 * ============================================================================ */

/**
 * @brief Create placeholder for cell complement and record for later resolution
 * 
 * Called from geom_parser.c when #cell_id is encountered.
 * Creates a complement node with invalid child (to be fixed later).
 * 
 * @param sys CSG system
 * @param referenced_cell_id The cell being referenced (#N)
 * @param referencing_cell_id The cell containing this reference
 * @return Placeholder node ID, or ALEA_NODE_ID_INVALID on error
 */
alea_node_id_t alea_record_cell_complement(alea_system_t* sys,
                                          int referenced_cell_id,
                                          int referencing_cell_id) {
    if (!sys) return ALEA_NODE_ID_INVALID;
    
    /* Create a placeholder complement node
     * We'll create it as complement of a "universe" node that we'll replace later.
     * For now, use ALEA_NODE_ID_INVALID as child - evaluation will fail until resolved.
     */
    alea_node_id_t placeholder = alea_alloc_node(sys);
    if (placeholder == ALEA_NODE_ID_INVALID) return ALEA_NODE_ID_INVALID;

    alea_node_t* node = &sys->nodes.data[placeholder];
    ALEA_SET_OPERATION(node, ALEA_OP_COMPLEMENT);
    node->operation.left = ALEA_NODE_ID_INVALID;  /* Will be resolved later */
    node->operation.right = ALEA_NODE_ID_INVALID;
    node->material_id = ALEA_MATERIAL_NONE;

    /* Set bbox to infinite - will be updated when resolved */
    node->bbox.min_x = -1e30; node->bbox.max_x = 1e30;
    node->bbox.min_y = -1e30; node->bbox.max_y = 1e30;
    node->bbox.min_z = -1e30; node->bbox.max_z = 1e30;

    /* Record the reference for later resolution */
    alea_cell_ref_t* ref = alea_vec_push_uninit(&sys->cell_refs, alea_cell_ref_t);
    if (!ref) {
        /* Allocation failed - can't easily undo node allocation, but return error */
        return ALEA_NODE_ID_INVALID;
    }
    ref->placeholder_node = placeholder;
    ref->referenced_cell_id = referenced_cell_id;
    ref->referencing_cell_id = referencing_cell_id;
    
    ALEA_LOG_INFO("  Recorded cell complement: #%d (placeholder node %u)",
           referenced_cell_id, placeholder);
    
    return placeholder;
}

/* ============================================================================
 * AFTER ALL CELLS PARSED: Resolve references
 * ============================================================================ */

/**
 * @brief Resolve all cell complement references
 * 
 * Call this after all cells have been converted. Walks through
 * recorded cell references and updates placeholder nodes to point
 * to the actual cell geometry.
 * 
 * @param sys CSG system with converted cells
 * @return Number of resolved references, or -1 on error
 */
int alea_resolve_cell_complements(alea_system_t* sys) {
    if (!sys) return -1;
    
    if (alea_vec_empty(&sys->cell_refs)) {
        return 0;  /* Nothing to resolve */
    }

    int resolved = 0;
    int errors = 0;

    for (size_t i = 0; i < alea_vec_count(&sys->cell_refs); i++) {
        alea_cell_ref_t* ref = &sys->cell_refs.data[i];
        
        /* Find the referenced cell's root node. The cell index is maintained
         * during load, so avoid refs * cells linear scans on large models. */
        alea_node_id_t cell_root = ALEA_NODE_ID_INVALID;
        int cell_idx = alea_find_cell_by_id(sys, ref->referenced_cell_id);
        if (cell_idx >= 0) {
            cell_root = sys->cells.data[cell_idx].root_node_id;
        }
        
        if (cell_root == ALEA_NODE_ID_INVALID) {
            ALEA_LOG_ERROR("Cell #%d referenced by cell %d not found",
                    ref->referenced_cell_id, ref->referencing_cell_id);
            errors++;
            continue;
        }
        
        /* Update the placeholder node to point to the cell's geometry */
        alea_node_t* node = &sys->nodes.data[ref->placeholder_node];
        node->operation.left = cell_root;
        
        /* Update bounding box (complement has infinite extent conceptually,
         * but we use the original cell's bbox for culling) */
        if (cell_root < alea_vec_count(&sys->nodes)) {
            node->bbox = sys->nodes.data[cell_root].bbox;
            /* Expand slightly since complement is "everything outside" */
            double expand = 1e10;
            node->bbox.min_x = -expand;
            node->bbox.max_x = expand;
            node->bbox.min_y = -expand;
            node->bbox.max_y = expand;
            node->bbox.min_z = -expand;
            node->bbox.max_z = expand;
        }
        
        ALEA_LOG_INFO("Resolved cell complement: #%d -> node %u (complement of node %u)",
               ref->referenced_cell_id, ref->placeholder_node, cell_root);
        
        resolved++;
    }
    
    if (errors > 0) {
        ALEA_LOG_ERROR("%d cell complement references could not be resolved", errors);
        return -1;
    }
    
    if (resolved > 0) {
        ALEA_LOG_INFO("Resolved %d cell complement references", resolved);
    }
    
    return resolved;
}

/**
 * @brief Check if all cell complements are resolved
 * 
 * @param sys CSG system
 * @return true if all resolved, false if any remain invalid
 */
bool alea_cell_complements_resolved(const alea_system_t* sys) {
    if (!sys) return true;

    for (size_t i = 0; i < alea_vec_count(&sys->cell_refs); i++) {
        const alea_cell_ref_t* ref = &sys->cell_refs.data[i];
        if (ref->placeholder_node >= alea_vec_count(&sys->nodes)) continue;

        const alea_node_t* node = &sys->nodes.data[ref->placeholder_node];
        if (node->operation.left == ALEA_NODE_ID_INVALID) {
            return false;
        }
    }

    return true;
}

/**
 * @brief Free cell reference tracking data
 *
 * Call from alea_system_destroy()
 */
void alea_free_cell_refs(alea_system_t* sys) {
    if (sys) {
        alea_vec_free(&sys->cell_refs);
    }
}
