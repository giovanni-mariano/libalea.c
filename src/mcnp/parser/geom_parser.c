// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file geom_parser.c
 * @brief Recursive descent parser for MCNP geometry expressions
 * 
 * Implements a recursive descent parser for MCNP boolean geometry.
 */

#include "geom_parser.h"
#include "geom_lexer.h"
#include "../conversion/surface_conv.h"
#include "core/alea_cell_complement.h"
#include "util/alea_log.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// Forward declarations for internal parser functions
// static alea_node_id_t parse_geometry_expression(geometry_lexer_t* lexer);
static alea_node_id_t parse_union_expression(geometry_lexer_t* lexer);
static alea_node_id_t parse_intersection_expression(geometry_lexer_t* lexer);
static alea_node_id_t parse_primary_expression(geometry_lexer_t* lexer);

/**
 * @brief Set cell context for better error reporting
 */
void lexer_set_cell_context(geometry_lexer_t* lexer, int cell_id, const char* geometry) {
    lexer->current_cell_id = cell_id;
    lexer->current_cell_geometry = geometry;
}

/**
 * @brief Top-level expression parser
 */
alea_node_id_t parse_geometry_expression(geometry_lexer_t* lexer) {
    alea_node_id_t result = parse_union_expression(lexer);
    return result;
}

/**
 * @brief Parse union operations (lowest precedence)
 * 
 * Grammar: union := intersection ( ':' intersection )*
 */
static alea_node_id_t parse_union_expression(geometry_lexer_t* lexer) {
    alea_node_id_t left = parse_intersection_expression(lexer);
    if (left == ALEA_NODE_ID_INVALID) return ALEA_NODE_ID_INVALID;
    
    while (lexer->token.type == GEOM_TOKEN_UNION) {
        geom_lexer_next_token(lexer); // consume ':'
        alea_node_id_t right = parse_intersection_expression(lexer);
        
        if (right == ALEA_NODE_ID_INVALID) {
            // Don't fail entirely if right side fails - might be trailing ':'
            break;
        }
        
        alea_node_id_t union_node = alea_create_union(lexer->system, left, right);
        if (union_node == ALEA_NODE_ID_INVALID) {
            return left; // Return what we have so far
        }
        left = union_node;
    }
    
    return left;
}

/**
 * @brief Parse intersection operations (higher precedence)
 * 
 * Grammar: intersection := primary ( ' ' primary )*
 * Note: Intersection is implicit (whitespace)
 */
static alea_node_id_t parse_intersection_expression(geometry_lexer_t* lexer) {
    alea_node_id_t left = parse_primary_expression(lexer);
    if (left == ALEA_NODE_ID_INVALID) return ALEA_NODE_ID_INVALID;
    
    // Look ahead to see if there are more terms without consuming
    while (lexer->token.type == GEOM_TOKEN_SURFACE_REF ||
           lexer->token.type == GEOM_TOKEN_CELL_REF ||
           lexer->token.type == GEOM_TOKEN_COMPLEMENT ||
           lexer->token.type == GEOM_TOKEN_LPAREN) {
        
        alea_node_id_t right = parse_primary_expression(lexer);
        if (right == ALEA_NODE_ID_INVALID) {
            break; // Can't parse right operand, return what we have
        }
        
        alea_node_id_t intersection_node = alea_create_intersection(lexer->system, left, right);
        if (intersection_node == ALEA_NODE_ID_INVALID) {
            return left; // Return what we have
        }
        left = intersection_node;
    }
    
    return left;
}
 
/**
 * @brief Parse primary expressions (highest precedence)
 * 
 * Grammar: primary := '(' expression ')' | '#' primary | surface_ref
 */
static alea_node_id_t parse_primary_expression(geometry_lexer_t* lexer) {
    // Handle parenthesized expressions
    if (lexer->token.type == GEOM_TOKEN_LPAREN) {
        geom_lexer_next_token(lexer); // consume '('
        
        // Parse the expression inside parentheses
        alea_node_id_t result = parse_geometry_expression(lexer);
        
        // Try to find the matching closing parenthesis
        if (lexer->token.type == GEOM_TOKEN_RPAREN) {
            geom_lexer_next_token(lexer); // consume ')'
        } else {
            // Error recovery: try to find matching ')'
            if (lexer->current_cell_id != -1) {
                ALEA_LOG_WARN("Cell %d - Expected ')' but found token type %d",
                       lexer->current_cell_id, lexer->token.type);
            }
            
            // Try to recover by looking for the closing parenthesis
            int paren_depth = 1;
            int max_tokens = 100;
            int tokens_skipped = 0;
            
            while (paren_depth > 0 && 
                   lexer->token.type != GEOM_TOKEN_END && 
                   lexer->token.type != GEOM_TOKEN_ERROR &&
                   tokens_skipped < max_tokens) {
                
                geom_lexer_next_token(lexer);
                tokens_skipped++;
                
                if (lexer->token.type == GEOM_TOKEN_LPAREN) {
                    paren_depth++;
                } else if (lexer->token.type == GEOM_TOKEN_RPAREN) {
                    paren_depth--;
                }
            }
            
            if (paren_depth == 0) {
                geom_lexer_next_token(lexer); // consume the final ')'
            } else {
                ALEA_LOG_ERROR("Could not find matching ')' - giving up");
                lexer->parse_errors++;
                return ALEA_NODE_ID_INVALID;
            }
        }
        
        return result;
    }
    
    // Handle complement operator #
    // #(expr) means NOT(expr) - complement of the expression
    // By De Morgan's law: #(A B) = -A:-B and #(A:B) = -A -B
    if (lexer->token.type == GEOM_TOKEN_COMPLEMENT) {
        geom_lexer_next_token(lexer); // consume '#'
        alea_node_id_t operand = parse_primary_expression(lexer);
        
        if (operand != ALEA_NODE_ID_INVALID) {
            return alea_create_complement(lexer->system, operand);
        }
        return ALEA_NODE_ID_INVALID;
    }
    
    // Handle surface references
    if (lexer->token.type == GEOM_TOKEN_SURFACE_REF) {
        int surface_id = lexer->token.value.surface_id;
        bool negative = lexer->token.negative;
        geom_lexer_next_token(lexer);

        alea_node_id_t node = alea_surface_node(lexer->system, surface_id, negative);
        
        if (node == ALEA_NODE_ID_INVALID) {
            snprintf(lexer->error_msg, sizeof(lexer->error_msg),
                    "Unknown surface %d in cell %d",
                    abs(surface_id), lexer->current_cell_id);
            lexer->parse_errors++;
            return ALEA_NODE_ID_INVALID;
        } 
        

        
        return node;
    }
    
    // Handle cell references (#123)
    if (lexer->token.type == GEOM_TOKEN_CELL_REF) {
        int cell_id = lexer->token.value.cell_id;
        geom_lexer_next_token(lexer);
        
        // Create placeholder and record for later resolution
        // Note: current_cell_id needs to be passed through the parser context
        alea_node_id_t placeholder = alea_record_cell_complement(
            lexer->system, 
            cell_id,
            lexer->current_cell_id  // Add this field to geom_lexer_t
        );
        
        if (placeholder == ALEA_NODE_ID_INVALID) {
            ALEA_LOG_ERROR("Failed to create cell complement for #%d", cell_id);
            return ALEA_NODE_ID_INVALID;
        }
        
        return placeholder;
    }   
    
    // Skip parameter keywords
    if (lexer->token.type == GEOM_TOKEN_FILL ||
        lexer->token.type == GEOM_TOKEN_UNIVERSE ||
        lexer->token.type == GEOM_TOKEN_MATERIAL ||
        lexer->token.type == GEOM_TOKEN_DENSITY ||
        lexer->token.type == GEOM_TOKEN_IMPORTANCE) {
        
        // Skip keyword
        geom_lexer_next_token(lexer);
        
        // Skip value (number or surface reference)
        if (lexer->token.type == GEOM_TOKEN_NUMBER || 
            lexer->token.type == GEOM_TOKEN_SURFACE_REF) {
            geom_lexer_next_token(lexer);
        }
        
        // Continue parsing after skipping parameter
        return parse_primary_expression(lexer);
    }
    
    // Handle LIKE expressions - extract template cell ID for later resolution
    if (lexer->token.type == GEOM_TOKEN_LIKE) {
        geom_lexer_next_token(lexer); // consume 'LIKE'

        // Parse template cell reference
        if (lexer->token.type == GEOM_TOKEN_SURFACE_REF ||
            lexer->token.type == GEOM_TOKEN_NUMBER) {
            int template_cell = (lexer->token.type == GEOM_TOKEN_NUMBER)
                                ? (int)lexer->token.value.number
                                : lexer->token.value.surface_id;
            geom_lexer_next_token(lexer);

            // Store template cell ID for resolution after all cells are parsed
            lexer->like_cell_id = template_cell;

            // Capture BUT clause if present (for parameter overrides)
            if (lexer->token.type == GEOM_TOKEN_BUT) {
                geom_lexer_next_token(lexer); // consume 'BUT'

                // Store pointer to BUT parameters for later parsing
                lexer->but_clause = lexer->current;

                // Skip to end
                while (lexer->token.type != GEOM_TOKEN_END) {
                    geom_lexer_next_token(lexer);
                }
            }

            // Return special marker indicating this is a LIKE cell
            return ALEA_NODE_ID_LIKE_PENDING;
        }

        return ALEA_NODE_ID_INVALID;
    }
    
    return ALEA_NODE_ID_INVALID;
}

/**
 * @brief Validate geometry syntax before parsing
 */
bool validate_geometry_syntax(const char* geometry_expr, char* error_msg, size_t error_size) {
    if (!geometry_expr || !*geometry_expr) {
        snprintf(error_msg, error_size, "Empty geometry expression");
        return false;
    }
    
    // Check for balanced parentheses
    int paren_count = 0;
    bool has_content = false;
    
    for (const char* p = geometry_expr; *p; p++) {
        if (*p == '(') {
            paren_count++;
        } else if (*p == ')') {
            paren_count--;
            if (paren_count < 0) {
                snprintf(error_msg, error_size, 
                        "Unmatched closing parenthesis at position %ld", p - geometry_expr);
                return false;
            }
        } else if (!isspace(*p) && *p != '$') {
            has_content = true;
        }
    }
    
    if (paren_count > 0) {
        snprintf(error_msg, error_size, 
                "Unmatched opening parenthesis (%d unclosed)", paren_count);
        return false;
    }
    
    if (!has_content) {
        snprintf(error_msg, error_size, "No meaningful content in geometry expression");
        return false;
    }
    
    return true;
}

alea_node_id_t parse_mcnp_geometry_with_like(
    alea_system_t* system,
    const char* geometry_expr,
    int cell_id,
    int* out_like_cell_id,
    const char** out_but_clause
) {
    if (out_like_cell_id) *out_like_cell_id = 0;
    if (out_but_clause) *out_but_clause = NULL;

    if (!geometry_expr || !*geometry_expr || !system) {
        return ALEA_NODE_ID_INVALID;
    }

    // Pre-validate geometry expression
    char error_msg[256];
    if (!validate_geometry_syntax(geometry_expr, error_msg, sizeof(error_msg))) {
        ALEA_LOG_ERROR("Geometry syntax error in cell %d: %s", cell_id, error_msg);
        return ALEA_NODE_ID_INVALID;
    }

    // Initialize lexer
    geometry_lexer_t lexer;
    geom_lexer_init(&lexer, geometry_expr, system);

    // Set cell context
    lexer.current_cell_id = cell_id;
    lexer.current_cell_geometry = geometry_expr;

    // Get first token
    geom_lexer_next_token(&lexer);

    if (lexer.token.type == GEOM_TOKEN_END) {
        return ALEA_NODE_ID_INVALID;
    }

    // Parse the geometry expression
    alea_node_id_t result = parse_geometry_expression(&lexer);

    // Check if LIKE was found
    if (lexer.like_cell_id > 0) {
        if (out_like_cell_id) *out_like_cell_id = lexer.like_cell_id;
        if (out_but_clause) *out_but_clause = lexer.but_clause;
        return ALEA_NODE_ID_LIKE_PENDING;
    }

    if (lexer.parse_errors > 0) {
        ALEA_LOG_ERROR("Parse errors in cell %d: %s", cell_id, lexer.error_msg);
        return ALEA_NODE_ID_INVALID;
    }

    return result;
}
