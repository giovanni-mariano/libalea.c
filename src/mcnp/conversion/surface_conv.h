// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef SURFACE_CONV_H
#define SURFACE_CONV_H

#include "core/alea_system.h"
#include "mcnp/parser/mcnp_parser.h"
#include <stdbool.h>

/**
 * @file surface_conv.h
 * @brief MCNP surface to CSG primitive conversion - Updated for flattened architecture
 */

// ============================================================================
// SURFACE CONVERSION - Simplified interface
// ============================================================================

typedef enum {
    SURFACE_CONV_SUCCESS = 0,
    SURFACE_CONV_INVALID_PARAMS,
    SURFACE_CONV_UNKNOWN_MNEMONIC,
    SURFACE_CONV_PARSE_ERROR,
    SURFACE_CONV_MEMORY_ERROR
} surface_conversion_result_t;

/**
 * @brief Build fast lookup table from MCNP surface IDs to CSG node IDs
 * 
 * This function:
 * 1. Finds the maximum MCNP surface ID in the surface map
 * 2. Allocates an array indexed by surface ID
 * 3. Populates it with corresponding CSG node IDs
 * 
 * Must be called after all surfaces are converted and before cell conversion.
 * 
 * @param sys CSG system with populated surface_map
 * @return 0 on success, -1 on allocation failure
 */
int alea_build_surface_lookup(alea_system_t* sys);

/**
 * @brief Convert an MCNP surface to a CSG primitive
 * 
 * This function:
 * 1. Parses the MCNP surface definition
 * 2. Creates/finds a primitive (with automatic deduplication)
 * 3. Stores the surface-to-primitive mapping
 * 
 * @param sys CSG system
 * @param mcnp_surf MCNP surface to convert
 * @return 0 on success, -1 on failure
 */
int alea_convert_surface(alea_system_t* sys, const mcnp_surface_t* mcnp_surf);

/**
 * @brief Find primitive ID for an MCNP surface ID
 * 
 * @param sys CSG system
 * @param mcnp_surface_id MCNP surface ID (can be positive or negative)
 * @return Primitive ID, or UINT32_MAX if not found
 */
uint32_t alea_find_surface_primitive(const alea_system_t* sys, int mcnp_surface_id);

/**
 * @brief Determine sense from MCNP surface ID
 * 
 * In MCNP:
 *   Positive surface ID = outside  (positive sense)
 *   Negative surface ID = inside (negative sense)
 * 
 * @param mcnp_surface_id MCNP surface ID
 * @return true for outside (positive), false for inside (negative)
 */
bool alea_surface_sense(int mcnp_surface_id);

/**
 * Get the appropriate sense node for a signed surface reference.
 * 
 * @param sys System
 * @param mcnp_surface_id Signed surface ID:
 *        +5 → positive sense (outside)
 *        -5 → negative sense (inside)
 * @return Node ID, or ALEA_NODE_ID_INVALID if not found
 */
alea_node_id_t alea_surface_node(const alea_system_t* sys, int mcnp_surface_id, bool negative);

#endif // SURFACE_CONV_H