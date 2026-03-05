// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_EXPORT_H
#define ALEA_EXPORT_H

#include "alea_types.h"
#include "util/arena.h"
#include <stdio.h>
#include <stdbool.h>


/**
 * @file alea_export.h
 * @brief CSG geometry export pipeline
 * 
 * Exports CSG systems to various formats (MCNP, OpenMC, Serpent).
 * Key features:
 * - Surface deduplication at export time
 * - Preserves original MCNP surface IDs where possible
 * - Format-specific writers via callbacks
 */

/* Forward declarations */

typedef struct alea_cell_entry alea_cell_entry_t;



/* ============================================================================
 * EXPORT SURFACE ENTRY
 * ============================================================================ */

/**
 * @brief Tracks a surface during export
 */
typedef struct {
    uint32_t primitive_id;      /* Index into sys->primitives */
    int export_id;              /* Surface ID in output file */
    int original_mcnp_id;       /* Original ID if from input, -1 if synthetic */
    bool emitted;               /* Already written to output? */
} export_surface_t;

/* ============================================================================
 * MACROBODY DECOMPOSITION
 * ============================================================================ */

/**
 * @brief Decomposed surface with sense
 * 
 * For expressing macrobody interiors as intersection of half-spaces.
 * sense: +1 means positive side, -1 means negative side
 */
typedef struct {
    int surface_id;    /* MCNP surface ID */
    int8_t sense;      /* +1 or -1 for half-space selection */
} decomposed_surface_t;

/**
 * @brief Macrobody decomposition entry
 * 
 * Maps an original macrobody surface to its constituent primitive surfaces.
 * The interior of the macrobody is the intersection of all surfaces
 * with their specified senses.
 */
typedef struct {
    int original_surface_id;           /* Original macrobody MCNP ID */
    decomposed_surface_t* surfaces;    /* Array of decomposed surfaces */
    size_t surface_count;              /* Number of surfaces (3 for RCC, 6 for BOX) */
} macrobody_decomposition_t;

/* ============================================================================
 * TRCL EXPORT MODE
 * ============================================================================ */

/**
 * @brief How to export cells with TRCL transforms
 */
typedef enum {
    ALEA_TRCL_EXPORT_PRESERVE = 0,  /* Original surfaces + TRCL param (default) */
    ALEA_TRCL_EXPORT_BAKE,          /* Baked: transformed surfaces, no TRCL param */
} alea_trcl_export_mode_t;

/* ============================================================================
 * TRANSFORM EXPORT MODE
 * ============================================================================ */

/**
 * @brief How to export transforms (FILL/TRCL parameters vs TR data cards)
 */
typedef enum {
    ALEA_TR_EXPORT_ORIGINAL = 0, /* Keep original style: inline stays inline, TRn stays card */
    ALEA_TR_EXPORT_INLINE,       /* Force all transforms inline in FILL/TRCL, no TR cards */
    ALEA_TR_EXPORT_CARDS,        /* Force all transforms as TR cards, FILL/TRCL use ID refs */
} alea_transform_export_mode_t;

/* ============================================================================
 * EXPORT CONTEXT
 * ============================================================================ */

/**
 * @brief Export state
 */
typedef struct export_context {
    /* Output */
    FILE* out;
    alea_export_format_t format;
    alea_surface_emit_policy_t surface_policy;  /* How to emit macrobodies */
    bool deduplicate;
    alea_trcl_export_mode_t trcl_export_mode;
    alea_transform_export_mode_t transform_export_mode;

    /* MCNP formatting */
    int mcnp_max_col;
    int mcnp_cont_indent;

    /* Universe depth filtering:
     * -1 = all universes (default)
     *  0 = only base universe (universe 0)
     *  N = include universes up to N levels deep */
    int universe_depth;

    /* FILL expansion (not yet implemented):
     *  0 = no expansion (default)
     *  N = expand N levels
     * -1 = fully flatten */
    int fill_depth;
    
    /* Dedup map: primitive_id -> canonical mcnp_surface_id */
    int* prim_to_surface;
    size_t prim_to_surface_size;

    /* Dedup map: primitive_id -> inverted flag of the canonical surface's pos_node.
     * Used to compute correct cell sense when surface coefficients are
     * un-canonicalized during export (MCNP exporter flips plane coefficients
     * back to original orientation). */
    int8_t* prim_to_surface_inverted;

    /* Macrobody decomposition tracking */
    macrobody_decomposition_t* decompositions;
    size_t decomposition_count;
    size_t decomposition_capacity;
    int next_synthetic_surface_id;  /* For auto-generated surfaces */
    
    /* Statistics */
    size_t surfaces_written;
    size_t cells_written;
    size_t transforms_written;
    size_t dedup_hits;
    size_t macrobodies_decomposed;
    size_t synthetic_surfaces_created;  /* Surfaces created for expanded macrobodies */
    arena_t arena;                     /* Arena for temporary allocations */

    /* Format-specific data (e.g., OpenMC material-density map) */
    void* mat_map;                     /* For OpenMC: pointer to mat_density_map_t */

    /* Module model data (e.g., mcnp_model_t* for MCNP export) */
    const void* module_data;

} export_context_t;

/* ============================================================================
 * EXPORT CONTEXT MANAGEMENT
 * ============================================================================ */

/**
 * @brief Create export context
 *
 * @param format Output format
 * @param surface_policy How to handle macrobodies
 * @param out Output file stream
 * @param deduplicate Enable surface deduplication
 * @param next_synthetic_surface_id Starting ID for synthetic surfaces
 * @param universe_depth Universe depth filter (-1=all, 0=base, N=N levels)
 * @param fill_depth FILL expansion depth (0=none, N=N levels, -1=full)
 */
export_context_t* export_context_create(alea_export_format_t format,
                                        alea_surface_emit_policy_t surface_policy,
                                        FILE* out,
                                        bool deduplicate,
                                        int next_synthetic_surface_id,
                                        int universe_depth,
                                        int fill_depth);

/**
 * @brief Destroy export context
 */
void export_context_destroy(export_context_t* ctx);

/* ============================================================================
 * MACROBODY DECOMPOSITION
 * ============================================================================ */

/**
 * @brief Check if a primitive type is a macrobody
 */
bool is_macrobody(alea_primitive_type_t type);

/**
 * @brief Find decomposition for an original surface ID
 * 
 * @param ctx Export context
 * @param original_surface_id Original macrobody surface ID
 * @return Pointer to decomposition, or NULL if not decomposed
 */
const macrobody_decomposition_t* find_decomposition(const export_context_t* ctx,
                                                     int original_surface_id);

/* ============================================================================
 * MACROBODY TREE EXPANSION
 * ============================================================================ */

/**
 * @brief Replace macrobody nodes in all cells with pre-computed expanded nodes
 *
 * For export formats that don't support macrobodies (OpenMC) or when
 * ALEA_EMIT_SURFACES policy is set. Uses the expanded_pos_node/expanded_neg_node
 * stored in surface entries at parse time.
 *
 * @param sys CSG system (cells are modified)
 * @return Number of cells with expanded macrobodies
 */
int alea_expand_macrobodies_in_cells(alea_system_t* sys);

/**
 * @brief Compute next synthetic surface ID (max existing + 1000)
 */
int alea_next_synthetic_surface_id(const alea_system_t* sys);

/* ============================================================================
 * INTERNAL FUNCTIONS (used by format-specific exporters)
 * ============================================================================ */

/**
 * @brief Assign surface IDs to primitives that don't have them
 * Used for primitives created during macrobody expansion.
 */
void alea_assign_missing_surface_ids(alea_system_t* sys, export_context_t* ctx);

/**
 * @brief Build canonical surface map for deduplication
 */
void alea_build_canonical_surface_map(export_context_t* ctx, const alea_system_t* sys);

/**
 * @brief Expand macrobodies at tree level using pre-computed nodes
 * @return Number of cells with expanded macrobodies, or -1 on error
 */
int alea_expand_macrobodies_tree_level(alea_system_t* sys, export_context_t* ctx);

/**
 * @brief Check if a cell should be exported based on universe depth filter
 */
bool alea_should_export_cell(const alea_system_t* sys, const alea_cell_entry_t* cell,
                            int max_depth, int* depth_cache, size_t cache_size);


#endif /* ALEA_EXPORT_H */