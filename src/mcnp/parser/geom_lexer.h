// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef GEOM_LEXER_H
#define GEOM_LEXER_H

#include "mcnp/conversion/mcnp_conversion.h"
#include <stdbool.h>
#include <stddef.h>

/**
 * @file geom_lexer.h
 * @brief Geometry expression lexer for MCNP cell geometry
 * 
 * This module tokenizes MCNP geometry expressions like:
 * "-1 : (2 -3 4)" into tokens (numbers, operators, parentheses)
 */

// Forward declare types

// Token types
typedef enum {
    GEOM_TOKEN_END = 0,           // End of input
    GEOM_TOKEN_SURFACE_REF,       // Surface ID (e.g., 1, -2)
    GEOM_TOKEN_CELL_REF,          // Cell reference (#123)
    GEOM_TOKEN_NUMBER,            // Floating point number
    GEOM_TOKEN_LPAREN,            // (
    GEOM_TOKEN_RPAREN,            // )
    GEOM_TOKEN_UNION,             // : (union operator)
    GEOM_TOKEN_INTERSECTION,      // Space (intersection operator)
    GEOM_TOKEN_COMPLEMENT,        // # (complement operator)
    GEOM_TOKEN_LIKE,              // LIKE keyword
    GEOM_TOKEN_BUT,               // BUT keyword
    GEOM_TOKEN_FILL,              // FILL keyword
    GEOM_TOKEN_UNIVERSE,          // U= keyword
    GEOM_TOKEN_MATERIAL,          // MAT= keyword
    GEOM_TOKEN_DENSITY,           // RHO= keyword
    GEOM_TOKEN_IMPORTANCE,        // IMP: keyword
    GEOM_TOKEN_ERROR              // Parsing error
} geometry_token_type_t;

// Token value union
typedef union {
    int surface_id;      // For GEOM_TOKEN_SURFACE_REF
    int cell_id;         // For GEOM_TOKEN_CELL_REF
    double number;       // For GEOM_TOKEN_NUMBER
} geometry_token_value_t;

// Token structure
typedef struct {
    geometry_token_type_t type;
    geometry_token_value_t value;
    bool negative;           // For negated surface refs
    const char* start;       // Start position in input
    const char* end;         // End position in input
} geometry_token_t;

// Lexer state
typedef struct geometry_lexer {
    const char* input;           // Original input string
    const char* current;         // Current position
    const char* end;             // End of input
    geometry_token_t token;      // Current token
    // Conversion context (unused, kept for reference)
    alea_system_t* system;           //
    alea_node_id_t* surface_nodes;   //
    size_t surface_count;
    int parse_errors;
    char error_msg[256];
    int current_cell_id;
    const char* current_cell_geometry;
    int like_cell_id;            // Template cell ID if LIKE was found (0 = none)
    const char* but_clause;      // Pointer to BUT clause parameters (NULL if none)
} geometry_lexer_t;

/**
 * @brief Initialize the geometry lexer
 * 
 * @param lexer Lexer to initialize
 * @param input Input string to tokenize
 * @param context Context for error reporting (can be NULL)
 */

void geom_lexer_init(geometry_lexer_t* lexer, const char* input,
                     alea_system_t* system);

/**
 * @brief Get the next token from the input
 * 
 * @param lexer Lexer state
 * @return Pointer to the current token (owned by lexer)
 */
const geometry_token_t* geom_lexer_next_token(geometry_lexer_t* lexer);

#endif // GEOM_LEXER_H