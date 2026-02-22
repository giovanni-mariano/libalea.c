// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "mcnp_lexer.h"
#include <ctype.h>
#include <string.h>

// Platform-specific case-insensitive comparison
#ifdef _WIN32
#include <strings.h>
#define strncasecmp _strnicmp
#else
#include <strings.h>
#endif

/**
 * @brief Get the next non-whitespace token from a string.
 * 
 * Advances the cursor past the token it found.
 * 
 * @param cursor [in/out] A pointer to a character pointer. It will be advanced.
 * @param token_start [out] Will be set to the start of the found token.
 * @return The length of the token. Returns 0 if no token is found.
 */
size_t mcnp_lexer_get_next_token(const char** cursor, const char** token_start) {
    // 1. Skip leading whitespace
    while (**cursor && isspace((unsigned char)**cursor)) {
        (*cursor)++;
    }

    if (**cursor == '\0') {
        *token_start = NULL;
        return 0; // End of string
    }

    // 2. Mark the start of the token
    *token_start = *cursor;

    // 3. Find the end of the token (the next whitespace)
    while (**cursor && !isspace((unsigned char)**cursor)) {
        (*cursor)++;
    }
    
    return *cursor - *token_start;
}

/**
 * @brief Checks if a token is a known MCNP cell parameter keyword.
 * 
 * This list can be expanded as needed.
 * 
 * @param token Token string
 * @param len Token length  
 * @return 1 if keyword matches, 0 otherwise
 */
int mcnp_lexer_is_parameter_keyword(const char* token, size_t len) {
    (void)len;  /* Token is null-terminated, length available for future bounds checking */
    const char* keywords[] = {
        "imp:", "u=", "vol=", "fill=", "lat=", "tmp=", "pwt=", 
        "ext:", "fcl=", "wwn", "dxc=", "trcl", "mat=", "rho="
        // NOTE: some keywords like 'imp' are followed by a colon, some by equals.
        // The comparison length handles this.
    };
    int num_keywords = sizeof(keywords) / sizeof(keywords[0]);

    for (int i = 0; i < num_keywords; i++) {
        size_t keyword_len = strlen(keywords[i]);
        
        // Use case-insensitive comparison
        if (strncasecmp(token, keywords[i], keyword_len) == 0) {
            return 1; // Match found
        }
    }
    return 0; // No match
}

/**
 * @brief Trim trailing whitespace from a line
 * 
 * @param line_end [in/out] Pointer to end of line (will be moved back)
 * @param line_start Start of line (boundary)
 */
void mcnp_lexer_trim_trailing_whitespace(const char** line_end, const char* line_start) {
    while (*line_end > line_start && isspace(*(*line_end - 1))) {
        (*line_end)--;
    }
}