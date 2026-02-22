// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef MCNP_PARSER_H
#define MCNP_PARSER_H

#include "util/arena.h"
#include <stddef.h>
#include <time.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file mcnp_parser.h
 * @brief MCNP input file parser 
 * 
 * Parses MCNP input files and extracts cell, surface, and material definitions.
 * All data is stored in a memory arena for fast allocation and bulk deallocation.
 * 
 * This is independent of the CSG system - it only handles parsing MCNP text
 * into structured data. Use mcnp_conversion.h to convert parsed data to CSG.
 */

/* ========================================================================== */
/* MCNP DATA STRUCTURES (Allocated from the arena)                            */
/* ========================================================================== */

/**
 * @brief MCNP Cell definition
 */
typedef struct mcnp_cell {
    int cell_id;
    int material_id;
    double density;
    double importance;         // Particle importance (IMP:N)
    char* geometry_definition; // The raw geometry string (e.g., "-1 2 -3")
    char* parameters;          // Other parameters like IMP:N, U, etc.
    char* comments;
} mcnp_cell_t;

/**
 * @brief Surface boundary types
 */
typedef enum {
    SURFACE_TRANSMISSIVE, // Default
    SURFACE_REFLECTIVE,   // From '*' prefix
    SURFACE_WHITE,        // From '+' prefix
    SURFACE_PERIODIC,     // From negative N value
    SURFACE_VACUUM        // Particles escape (detected from graveyard cell)
} mcnp_boundary_type_t;

/**
 * @brief MCNP Surface definition
 */
typedef struct mcnp_surface {
    int surface_id;
    char* mnemonic;            // Surface type (PX, SO, CX, etc.)
    char* coefficients_str;    // Numeric parameters
    char* comments;
    
    // Additional surface properties
    mcnp_boundary_type_t boundary_type;
    int transform_id;        // From positive N value
    int periodic_surface_id; // From negative N value (the paired surface)
} mcnp_surface_t;

/**
 * @brief MCNP Transform definition (TRn card)
 */
typedef struct {
    int transform_id;           // n in TRn
    char* definition;           // Raw parameters (ox oy oz [b1..b9] [M])
    int is_star;                // 1 if *TRn (degrees), 0 if TRn (cosines)
} mcnp_transform_t;

/**
 * @brief MCNP Material definition
 */
typedef struct {
    int material_id;
    char* definition;          // The raw material definition (nuclide pairs)
    char* comments;
} mcnp_material_t;

/* ========================================================================== */
/* MCNP PARSER CONTEXT (Manages the arena and parsed data)                    */
/* ========================================================================== */

/**
 * @brief MCNP parser context
 * 
 * Contains all parsed data from an MCNP file. All data is allocated
 * from a single memory arena for fast allocation and bulk deallocation.
 * 
 * Lifecycle:
 *   1. mcnp_parse_file() - Creates context, allocates from arena
 *   2. Use the data (cells, surfaces, materials)
 *   3. mcnp_context_destroy() - Frees all arena memory at once
 * 
 * Note: This is TEMPORARY data used only during conversion.
 *       After converting to CSG, destroy this context to free memory.
 */
typedef struct mcnp_context {
    arena_t arena;

    // --- File Information ---
    char* source_filename;
    char* title_card;
    time_t parse_time;

    // --- Dynamic arrays of pointers to data within the arena ---
    mcnp_cell_t** cells;
    size_t cell_count;
    size_t cell_capacity;

    mcnp_surface_t** surfaces;
    size_t surface_count;
    size_t surface_capacity;

    mcnp_material_t** materials;
    size_t material_count;
    size_t material_capacity;

    mcnp_transform_t** transforms;
    size_t transform_count;
    size_t transform_capacity;
    
    // Parse statistics
    size_t parse_errors;
    size_t parse_warnings;
} mcnp_context_t;

/* ========================================================================== */
/* PARSER API (Internal - used by conversion layer)                          */
/* ========================================================================== */

/**
 * @brief Creates a new MCNP context with a dedicated memory arena.
 * 
 * Note: This is an internal function. External users should use
 *       mcnp_convert_file() from mcnp_conversion.h instead.
 *
 * @param filename The name of the file being parsed (can be NULL).
 * @return A pointer to the new context, or NULL on failure.
 */
mcnp_context_t* mcnp_context_create(const char* filename);

/**
 * @brief Destroys the MCNP context and frees all associated memory.
 *
 * Thanks to the arena, this is incredibly simple and fast.
 * All cells, surfaces, materials, and strings are freed in one call.
 * 
 * Note: This is an internal function. The conversion layer handles
 *       context cleanup automatically.
 * 
 * @param ctx The context to destroy.
 */
void mcnp_context_destroy(mcnp_context_t* ctx);

/**
 * @brief Parses an MCNP input file and populates a new context.
 *
 * This function orchestrates the memory-mapping and parsing of the file,
 * allocating all resulting data into a new context.
 * 
 * Note: This is an internal function. External users should use
 *       mcnp_convert_file() from mcnp_conversion.h which handles
 *       parsing AND conversion in one call.
 * 
 * Internal usage:
 *   mcnp_context_t* mcnp = NULL;
 *   if (mcnp_parse_file("input.mcnp", &mcnp)) {
 *       // Use mcnp->cells, mcnp->surfaces, etc.
 *       mcnp_context_destroy(mcnp);
 *   }
 *
 * @param filename The path to the MCNP input file.
 * @param out_context A pointer to a context pointer that will be set to the newly created context.
 * @return 1 on success, 0 on failure.
 */
int mcnp_parse_file(const char* filename, mcnp_context_t** out_context);

/**
 * @brief Print summary of parsed data (for debugging)
 * 
 * @param ctx The context to summarize
 */
void mcnp_context_print_summary(const mcnp_context_t* ctx);

#ifdef __cplusplus
}
#endif

#endif // MCNP_PARSER_H