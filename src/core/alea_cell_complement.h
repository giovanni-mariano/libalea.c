// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_CELL_COMPLEMENT_H
#define ALEA_CELL_COMPLEMENT_H

/**
 * @file alea_cell_complement.h
 * @brief Cell complement (#cell) resolution
 * 
 * Handles MCNP cell complement operator, e.g.:
 *   3 0  1 #2 #7
 * Means: inside surface 1 AND NOT(cell 2 geometry) AND NOT(cell 7 geometry)
 */

#include "alea_types.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Forward declaration */
typedef struct alea_system alea_system_t;

/**
 * @brief Record of an unresolved cell complement reference
 */
typedef struct {
    alea_node_id_t placeholder_node;  /* Node that needs to be updated */
    int referenced_cell_id;          /* MCNP cell ID being referenced */
    int referencing_cell_id;         /* MCNP cell ID containing the reference */
} alea_cell_ref_t;

/**
 * @brief Create placeholder for cell complement during parsing
 * 
 * Called from geom_parser.c when #cell_id is encountered.
 * Creates a complement node and records it for later resolution.
 * 
 * @param sys CSG system
 * @param referenced_cell_id The cell being referenced (#N)
 * @param referencing_cell_id The cell containing this reference (for error msgs)
 * @return Placeholder node ID, or ALEA_NODE_ID_INVALID on error
 */
alea_node_id_t alea_record_cell_complement(alea_system_t* sys,
                                          int referenced_cell_id,
                                          int referencing_cell_id);

/**
 * @brief Resolve all cell complement references
 * 
 * Call AFTER all cells have been converted. Updates placeholder nodes
 * to point to the actual referenced cell geometry trees.
 * 
 * @param sys CSG system with converted cells
 * @return Number of resolved references, or -1 on error
 */
int alea_resolve_cell_complements(alea_system_t* sys);

/**
 * @brief Check if all cell complements are resolved
 * 
 * @param sys CSG system
 * @return true if all resolved, false if any remain invalid
 */
bool alea_cell_complements_resolved(const alea_system_t* sys);

/**
 * @brief Free cell reference tracking data
 * 
 * Called from alea_system_destroy()
 */
void alea_free_cell_refs(alea_system_t* sys);

#ifdef __cplusplus
}
#endif

#endif /* ALEA_CELL_COMPLEMENT_H */