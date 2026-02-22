// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef CELL_CONV_H
#define CELL_CONV_H

#include "core/alea_system.h"
#include "mcnp/parser/mcnp_parser.h"
#include <stddef.h>
#include <ctype.h>

/**
 * @file cell_conv.h
 * @brief MCNP cell to CSG tree conversion - Updated for flattened architecture
 */

// ============================================================================
// CELL CONVERSION 
// ============================================================================

typedef struct {
    double imp_n, imp_p, imp_e;
    double vol, tmp;
    int universe_id, fill_universe, fill_transform, lat_type;

    // Material/density overrides (for LIKE BUT)
    int material_id;
    double density;

    // Inline transform support for FILL
    int fill_transform_inline;      // 1 = inline transform, 0 = transform ID reference
    int fill_transform_degrees;     // 1 = *FILL (angles in degrees), 0 = FILL
    int fill_transform_count;       // Number of values parsed: 0, 3, or 12
    double fill_transform_data[12]; // ox oy oz [b1 b2 b3 b4 b5 b6 b7 b8 b9]

    // Lattice FILL support (FILL=i1:i2 j1:j2 k1:k2 followed by universe array)
    int lat_fill_dims[6];           // imin, imax, jmin, jmax, kmin, kmax
    int* lat_fill;                  // Array of universe IDs (dynamically allocated)
    size_t lat_fill_count;          // Number of universes in array
    int lat_fill_is_array;          // 1 = lattice array format

    // TRCL - cell transformation
    int trcl;                       // Transform ID reference
    int trcl_inline;                // 1 = inline transform, 0 = transform ID reference
    int trcl_degrees;               // 1 = *TRCL (angles in degrees)
    int trcl_count;                 // Number of values parsed: 0, 3, or 12
    double trcl_data[12];           // ox oy oz [b1 b2 b3 b4 b5 b6 b7 b8 b9]

    // Flags indicating which parameters were present
    unsigned int has_imp_n : 1;
    unsigned int has_imp_p : 1;
    unsigned int has_imp_e : 1;
    unsigned int has_vol : 1;
    unsigned int has_tmp : 1;
    unsigned int has_universe : 1;
    unsigned int has_fill : 1;
    unsigned int has_lat : 1;
    unsigned int has_trcl : 1;
    unsigned int has_mat : 1;
    unsigned int has_rho : 1;
} alea_cell_params_t;

/**
 * @brief Parse MCNP cell parameters string
 * 
 * Parses the parameters portion of an MCNP cell card (everything after
 * the geometry expression). Extracts importance, universe, fill, lattice,
 * and other parameters.
 * 
 * Example input: "IMP:N=1 U=2 FILL=3"
 * 
 * @param params_str The parameters string to parse (can be NULL)
 * @param out_params Output structure for parsed parameters
 * @return 0 on success, -1 on parse error
 */
int parse_cell_parameters(const char* params_str, alea_cell_params_t* out_params);

/**
 * @brief Convert an MCNP cell to a CSG tree
 * 
 * This function:
 * 1. Parses the MCNP geometry expression
 * 2. Builds a CSG tree from the expression
 * 3. Stores the cell information
 * 
 * @param sys CSG system
 * @param cell MCNP cell to convert
 * @return Root node ID of the CSG tree, or UINT32_MAX on failure
 */
uint32_t alea_convert_cell(alea_system_t* sys, const mcnp_cell_t* cell);

// ============================================================================
// TRCL TRANSFORM APPLICATION
// ============================================================================

/**
 * @brief Apply TRCL transforms to all cells that have them
 *
 * For each cell with has_trcl set, this function:
 * 1. Builds the TRCL transform matrix
 * 2. Clones the cell's CSG tree with the transform applied
 * 3. Updates the cell's root_node_id to the transformed tree
 *
 * This should be called after cell conversion and before universe indexing.
 * After calling this, cells with TRCL will have their geometry in
 * world coordinates - the transform has been "baked in".
 *
 * @param sys CSG system
 * @return Number of cells transformed, or -1 on error
 */
int alea_apply_trcl_transforms(alea_system_t* sys);

/**
 * @brief Resolve LIKE BUT cells by cloning template geometry
 *
 * For each cell with like_cell_id set, this function:
 * 1. Finds the template cell by ID
 * 2. Clones the template's CSG tree
 * 3. Updates the cell's root_node_id with the cloned tree
 *
 * This should be called after all cells are converted and before TRCL
 * transforms are applied. The BUT parameter overrides (MAT, RHO, TRCL,
 * FILL, etc.) are already stored during cell conversion.
 *
 * @param sys CSG system
 * @return Number of cells resolved, or -1 on error
 */
int alea_resolve_like_cells(alea_system_t* sys);

#endif // CELL_CONV_H