// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef CELL_CONV_H
#define CELL_CONV_H

#include "core/alea_system.h"
#include "mcnp/mcnp_model.h"
#include "mcnp/parser/mcnp_parser.h"
#include <stddef.h>
#include <ctype.h>

/**
 * @file cell_conv.h
 * @brief MCNP cell to CSG tree conversion
 */

// ============================================================================
// CELL CONVERSION
// ============================================================================

/**
 * @brief Intermediate parsed cell parameters (from MCNP text)
 *
 * This is a temporary parse-time struct, not the final storage format.
 * After parsing, values are split between alea_cell_entry_t (core fields)
 * and mcnp_cell_params_t (MCNP-specific fields).
 */
typedef struct {
    double imp_n, imp_p, imp_e;

    #define X_FIELD(name, type, kw, prec) type name;
    MCNP_CELL_SIMPLE_PARAMS(X_FIELD)
    #undef X_FIELD

    int universe_id, fill_universe, fill_transform, lat_type;

    // Material/density overrides (for LIKE BUT)
    int material_id;
    double density;

    // Inline transform support for FILL
    int fill_transform_inline;
    int fill_transform_degrees;
    int fill_transform_count;
    double fill_transform_data[12];

    // Lattice FILL support
    int lat_fill_dims[6];
    int* lat_fill;
    size_t lat_fill_count;
    int lat_fill_is_array;

    // TRCL - cell transformation
    int trcl;
    int trcl_inline;
    int trcl_degrees;
    int trcl_count;
    double trcl_data[12];

    // Flags indicating which parameters were present
    unsigned int has_imp_n : 1;
    unsigned int has_imp_p : 1;
    unsigned int has_imp_e : 1;

    #define X_HAS(name, type, kw, prec) unsigned int has_##name : 1;
    MCNP_CELL_SIMPLE_PARAMS(X_HAS)
    #undef X_HAS

    unsigned int has_universe : 1;
    unsigned int has_fill : 1;
    unsigned int has_lat : 1;
    unsigned int has_trcl : 1;
    unsigned int has_mat : 1;
    unsigned int has_rho : 1;
} alea_cell_params_t;

/**
 * @brief Parse MCNP cell parameters string
 */
int parse_cell_parameters(const char* params_str, alea_cell_params_t* out_params,
                          int cell_id);

/**
 * @brief Convert an MCNP cell to a CSG tree
 *
 * Writes core fields to alea_cell_entry_t and MCNP-specific fields
 * to the model's parallel cell_params array.
 *
 * @param sys CSG system
 * @param cell MCNP cell to convert
 * @param model MCNP model (for writing cell params)
 * @return Root node ID of the CSG tree, or UINT32_MAX on failure
 */
uint32_t alea_convert_cell(alea_system_t* sys, const mcnp_cell_t* cell,
                           mcnp_model_t* model);

/**
 * @brief Apply TRCL transforms to all cells that have them
 *
 * Reads TRCL data from model->cell_params.
 *
 * @param sys CSG system
 * @param model MCNP model
 * @return Number of cells transformed, or -1 on error
 */
int alea_apply_trcl_transforms(alea_system_t* sys, mcnp_model_t* model);

/**
 * @brief Resolve LIKE BUT cells by cloning template geometry
 *
 * Reads/writes LIKE data from model->cell_params.
 *
 * @param sys CSG system
 * @param model MCNP model
 * @return Number of cells resolved, or -1 on error
 */
int alea_resolve_like_cells(alea_system_t* sys, mcnp_model_t* model);

#endif // CELL_CONV_H
