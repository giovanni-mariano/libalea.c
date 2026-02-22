// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file geom_lexer.c
 * @brief Geometry expression lexer for MCNP cell geometry
 * 
 * Tokenizes MCNP geometry expressions like "-1 : (2 -3 4)" into tokens.
 * Extracted from mcnp_geometry_parser.c (lines 18-299)
 */

#include "geom_lexer.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <math.h>

/**
 * @brief Initialize the geometry lexer
 */
void geom_lexer_init(geometry_lexer_t* lexer, const char* input,
                     alea_system_t* system) {
    lexer->input = input;
    lexer->current = input;
    lexer->end = input + strlen(input);
    lexer->system = system;
    lexer->surface_nodes = system->surface_lookup;
    lexer->surface_count = system->surface_lookup_size;
    lexer->parse_errors = 0;
    lexer->error_msg[0] = '\0';
    lexer->current_cell_id = -1;
    lexer->current_cell_geometry = input;
    lexer->like_cell_id = 0;  // 0 = not a LIKE cell
    lexer->but_clause = NULL; // No BUT clause by default
    memset(&lexer->token, 0, sizeof(lexer->token));
}

/**
 * @brief Skip whitespace and comments
 */
static void geom_lexer_skip_whitespace(geometry_lexer_t* lexer) {
    while (lexer->current < lexer->end) {
        char ch = *lexer->current;
        
        // Skip all whitespace characters
        if (isspace(ch)) {
            lexer->current++;
            continue;
        }
        
        // Skip MCNP comments that start with $
        if (ch == '$') {
            // Skip comment to end of line or end of string
            while (lexer->current < lexer->end && *lexer->current != '\n' && 
                   *lexer->current != '\0') {
                lexer->current++;
            }
            continue;
        }
        
        // Skip MCNP comment lines that start with 'c ' or 'C '
        if ((ch == 'c' || ch == 'C') && lexer->current + 1 < lexer->end) {
            char next_ch = lexer->current[1];
            if (isspace(next_ch)) {
                // This is a comment line, skip the rest of the line
                while (lexer->current < lexer->end && *lexer->current != '\n' && 
                       *lexer->current != '\0') {
                    lexer->current++;
                }
                continue;
            }
        }
        
        // If we get here, we've found non-whitespace, non-comment content
        break;
    }
}

/**
 * @brief Parse number (integer or float)
 */
static bool geom_lexer_parse_number(const char* start, const char* end, double* result) {
    if (start >= end) return false;
    
    size_t len = end - start;
    if (len == 0 || len > 31) return false;
    
    char temp[32];
    memcpy(temp, start, len);
    temp[len] = '\0';
    
    char* endptr;
    *result = strtod(temp, &endptr);
    
    // Check if entire string was consumed and result is valid
    return (endptr == temp + len) && !isnan(*result) && !isinf(*result);
}

/**
 * @brief Get next token from the input
 */
const geometry_token_t* geom_lexer_next_token(geometry_lexer_t* lexer) {
    // Skip whitespace aggressively at the start
    geom_lexer_skip_whitespace(lexer);
    
    if (lexer->current >= lexer->end) {
        lexer->token.type = GEOM_TOKEN_END;
        return &lexer->token;
    }
    
    lexer->token.start = lexer->current;
    lexer->token.negative = false;
    lexer->token.end = lexer->current;
    
    char ch = *lexer->current;
    
    // Single character tokens
    switch (ch) {
        case '(':
            lexer->token.type = GEOM_TOKEN_LPAREN;
            lexer->current++;
            lexer->token.end = lexer->current;
            return &lexer->token;
            
        case ')':
            lexer->token.type = GEOM_TOKEN_RPAREN;
            lexer->current++;
            lexer->token.end = lexer->current;
            return &lexer->token;
            
        case ':':
            lexer->token.type = GEOM_TOKEN_UNION;
            lexer->current++;
            lexer->token.end = lexer->current;
            return &lexer->token;
            
        case '#':
            lexer->current++; // Skip #
            lexer->token.type = GEOM_TOKEN_COMPLEMENT;
            
            // Parse cell number after # if present
            geom_lexer_skip_whitespace(lexer);
            if (lexer->current < lexer->end && isdigit(*lexer->current)) {
                const char* num_start = lexer->current;
                while (lexer->current < lexer->end && isdigit(*lexer->current)) {
                    lexer->current++;
                }
                lexer->token.value.cell_id = atoi(num_start);
                lexer->token.type = GEOM_TOKEN_CELL_REF;
            }
            lexer->token.end = lexer->current;
            return &lexer->token;
    }
    
    // Handle signed numbers and surface references
    if (ch == '+' || ch == '-' || isdigit(ch) || ch == '.') {
        bool negative = false;
        const char* number_start = lexer->current;
        
        // Handle sign
        if (ch == '-') {
            negative = true;
            lexer->current++;
        } else if (ch == '+') {
            lexer->current++;
        }
        
        // Skip whitespace after sign but before digits
        if (ch == '-' || ch == '+') {
            geom_lexer_skip_whitespace(lexer);
        }
        
        // Must have at least one digit after sign (and optional whitespace)
        if (lexer->current >= lexer->end || !isdigit(*lexer->current)) {
            // Not a valid number, reset and treat as intersection
            lexer->current = number_start + 1;
            lexer->token.type = GEOM_TOKEN_INTERSECTION;
            lexer->token.end = lexer->current;
            return &lexer->token;
        }
        
        // Parse the numeric part
        const char* num_start = lexer->current;
        bool has_decimal = false;
        bool has_exp = false;
        
        // Parse integer part
        while (lexer->current < lexer->end && isdigit(*lexer->current)) {
            lexer->current++;
        }
        
        // Parse decimal part if present
        if (lexer->current < lexer->end && *lexer->current == '.' && !has_decimal) {
            has_decimal = true;
            lexer->current++;
            while (lexer->current < lexer->end && isdigit(*lexer->current)) {
                lexer->current++;
            }
        }
        
        // Parse exponent if present
        if (lexer->current < lexer->end && 
            (*lexer->current == 'e' || *lexer->current == 'E') && !has_exp) {
            has_exp = true;
            lexer->current++;
            // Handle exponent sign
            if (lexer->current < lexer->end && 
                (*lexer->current == '+' || *lexer->current == '-')) {
                lexer->current++;
            }
            // Parse exponent digits
            while (lexer->current < lexer->end && isdigit(*lexer->current)) {
                lexer->current++;
            }
        }
        
        // Convert to number
        double value;
        if (geom_lexer_parse_number(num_start, lexer->current, &value)) {
            // Determine if it's a floating point number or integer surface reference
            if (has_decimal || has_exp) {
                // Contains decimal point or scientific notation
                lexer->token.type = GEOM_TOKEN_NUMBER;
                lexer->token.value.number = negative ? -value : value;
            } else {
                // Pure integer - treat as surface reference
                lexer->token.type = GEOM_TOKEN_SURFACE_REF;
                lexer->token.value.surface_id = (int)value;
                lexer->token.negative = negative;
            }
            lexer->token.end = lexer->current;
            return &lexer->token;
        } else {
            // Parsing failed
            lexer->token.type = GEOM_TOKEN_ERROR;
            snprintf(lexer->error_msg, sizeof(lexer->error_msg), 
                    "Invalid number at position %ld", number_start - lexer->input);
            lexer->token.end = lexer->current;
            return &lexer->token;
        }
    }
    
    // Handle keywords and identifiers
    if (isalpha(ch)) {
        const char* word_start = lexer->current;
        while (lexer->current < lexer->end && 
               (isalnum(*lexer->current) || *lexer->current == '_' || 
                *lexer->current == '/' || *lexer->current == ':' || 
                *lexer->current == '=')) {
            lexer->current++;
        }
        
        const char* word_end = lexer->current;
        size_t word_len = word_end - word_start;
        
        // Case-insensitive keyword matching
        if (word_len >= 4 && strncasecmp(word_start, "LIKE", 4) == 0) {
            lexer->token.type = GEOM_TOKEN_LIKE;
        } else if (word_len >= 3 && strncasecmp(word_start, "BUT", 3) == 0) {
            lexer->token.type = GEOM_TOKEN_BUT;
        } else if (word_len >= 4 && strncasecmp(word_start, "FILL", 4) == 0) {
            lexer->token.type = GEOM_TOKEN_FILL;
        } else if (word_len >= 5 && strncasecmp(word_start, "*FILL", 5) == 0) {
            lexer->token.type = GEOM_TOKEN_FILL;
        } else if (word_len >= 2 && strncasecmp(word_start, "U=", 2) == 0) {
            lexer->token.type = GEOM_TOKEN_UNIVERSE;
        } else if (word_len >= 4 && strncasecmp(word_start, "MAT=", 4) == 0) {
            lexer->token.type = GEOM_TOKEN_MATERIAL;
        } else if (word_len >= 4 && strncasecmp(word_start, "RHO=", 4) == 0) {
            lexer->token.type = GEOM_TOKEN_DENSITY;
        } else if (word_len >= 4 && strncasecmp(word_start, "IMP:", 4) == 0) {
            lexer->token.type = GEOM_TOKEN_IMPORTANCE;
        } else {
            // Unknown keyword - treat as intersection for now
            lexer->token.type = GEOM_TOKEN_INTERSECTION;
        }
        
        lexer->token.end = lexer->current;
        return &lexer->token;
    }
    
    // Unknown character - skip it and continue
    lexer->current++;
    lexer->token.type = GEOM_TOKEN_INTERSECTION;
    lexer->token.end = lexer->current;
    return &lexer->token;
}