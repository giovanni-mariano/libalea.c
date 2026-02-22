// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef MCNP_LEXER_H
#define MCNP_LEXER_H

#include <stddef.h>

/**
 * @file mcnp_lexer.h
 * @brief MCNP lexer/tokenizer - Extracts tokens from MCNP input
 * 
 * This module handles:
 * - Tokenization of MCNP input lines
 * - Whitespace handling
 * - Token classification (parameters, geometry, etc.)
 */

/**
 * @brief Get the next token from an MCNP line
 * 
 * Advances the cursor past whitespace and returns the next token.
 * 
 * @param cursor [in/out] Pointer to current position (will be advanced)
 * @param token_start [out] Will point to start of found token
 * @return Length of token, or 0 if no token found
 */
size_t mcnp_lexer_get_next_token(const char** cursor, const char** token_start);

/**
 * @brief Check if a token is an MCNP parameter keyword
 * 
 * Parameter keywords include: imp:, u=, vol=, fill=, lat=, tmp=, etc.
 * 
 * @param token Token string
 * @param len Token length
 * @return 1 if parameter keyword, 0 otherwise
 */
int mcnp_lexer_is_parameter_keyword(const char* token, size_t len);

/**
 * @brief Trim trailing whitespace from a line
 * 
 * @param line_end [in/out] Pointer to end of line (will be moved back)
 * @param line_start Start of line (boundary)
 */
void mcnp_lexer_trim_trailing_whitespace(const char** line_end, const char* line_start);

#endif // MCNP_LEXER_H