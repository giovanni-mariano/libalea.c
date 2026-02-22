// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef GEOM_PARSER_H
#define GEOM_PARSER_H


#include "core/alea_ops.h"
#include "mcnp/conversion/mcnp_conversion.h"
#include "primitives/primitive_create.h"
#include "geom_lexer.h"
/**
 * @file geom_parser.h
 * @brief Recursive descent parser for MCNP geometry expressions
 * 
 * Parses boolean geometry expressions and builds CSG trees.
 * Grammar:
 *   expression   := union
 *   union        := intersection ( ':' intersection )*
 *   intersection := primary ( ' ' primary )*
 *   primary      := NUMBER | '#' primary | '(' expression ')'
 */

// Forward declare types
typedef struct alea_node alea_node_t;
typedef struct alea_node_pool alea_node_pool_t;

/**
 * @brief Set cell context for better error messages
 * 
 * @param lexer Lexer state
 * @param cell_id Cell ID being parsed
 * @param geometry Geometry string being parsed
 */
void lexer_set_cell_context(geometry_lexer_t* lexer, int cell_id, 
                            const char* geometry);

/**
 * @brief Validate geometry syntax before parsing
 * 
 * @param geometry_expr Geometry expression to validate
 * @param error_msg Buffer for error message
 * @param error_size Size of error buffer
 * @return true if syntax is valid, false otherwise
 */
bool validate_geometry_syntax(const char* geometry_expr, 
                              char* error_msg, size_t error_size);

alea_node_id_t parse_geometry_expression(geometry_lexer_t* lexer);

/**
 * Parse MCNP geometry expression and extract LIKE reference if present
 *
 * @param system CSG system
 * @param geometry_expr MCNP geometry expression string
 * @param cell_id ID of cell being parsed
 * @param out_like_cell_id Output: template cell ID if LIKE (0 = not LIKE)
 * @param out_but_clause Output: pointer to BUT clause params (NULL if none)
 * @return Root node of CSG tree, or ALEA_NODE_ID_LIKE_PENDING for LIKE cells
 */
alea_node_id_t parse_mcnp_geometry_with_like(
    alea_system_t* system,
    const char* geometry_expr,
    int cell_id,
    int* out_like_cell_id,
    const char** out_but_clause
);

#endif // GEOM_PARSER_H