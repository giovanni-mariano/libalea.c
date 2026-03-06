// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_SURFACE_H
#define ALEA_SURFACE_H

/**
 * @file alea_surface.h
 * @brief Surface mapping types and conversion operations
 *
 * Surface-to-primitive mapping, surface vec, and MCNP conversion functions.
 */

#include "alea_types.h"
#include "util/alea_vec.h"


// Forward declarations 
typedef struct mcnp_surface mcnp_surface_t;
typedef struct mcnp_cell mcnp_cell_t;

// ============================================================================
// SURFACE MAPPING
// ============================================================================

/**
 * @brief Maps MCNP surface ID to primitive ID
 */
typedef struct {
    int mc_surface_id;
    uint32_t primitive_id;           // should the name be node_id?
    int transform_id;           // 0 = none, >0 = TRn was applied
    bool inverted;
} surface_mapping_t;

/* alea_surface_entry_t is defined in alea_types.h */

ALEA_VEC_DEFINE(alea_surface_vec, alea_surface_entry_t);

// ============================================================================
// API - CONVERSION (from MCNP)
// ============================================================================

/**
 * @brief Convert MCNP surface to CSG primitive
 *
 * @param sys CSG system
 * @param surf MCNP surface
 * @return 0 on success, -1 on failure
 */
int alea_convert_surface(alea_system_t* sys, const mcnp_surface_t* surf);


// ============================================================================
// API - SURFACE UTILITIES
// ============================================================================

/**
 * @brief Build surface index for all cells
 *
 * Traverses each cell's CSG tree and collects all surface indices
 * that the cell references. This enables faster raycasting by only
 * testing surfaces relevant to each cell instead of all surfaces.
 *
 * Call this after loading/creating geometry and before raycasting.
 *
 * @param sys CSG system
 * @return 0 on success, -1 on failure
 */
int alea_build_cell_surface_index(alea_system_t* sys);


#endif /* ALEA_SURFACE_H */
