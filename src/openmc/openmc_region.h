// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file openmc_region.h
 * @brief Parser for OpenMC region expressions
 *
 * OpenMC region syntax:
 *   - Space = intersection (AND)
 *   - | = union (OR)
 *   - ~ = complement (NOT)
 *   - -N = inside surface N (negative half-space)
 *   - N = outside surface N (positive half-space)
 *   - (...) = grouping
 *
 * Example: "-1 2 | 3" means (inside surface 1 AND outside surface 2) OR outside surface 3
 */

#ifndef OPENMC_REGION_H
#define OPENMC_REGION_H

#include "core/alea_system.h"
#include "util/arena.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * REGION PARSE CONTEXT
 * ============================================================================ */

/**
 * @brief Pair of pre-created sense nodes for a surface
 */
typedef struct {
    alea_node_id_t pos_node;     /**< +S node (outside), ALEA_NODE_ID_INVALID if unset */
    alea_node_id_t neg_node;     /**< -S node (inside),  ALEA_NODE_ID_INVALID if unset */
} openmc_surface_node_pair_t;

/**
 * @brief Context for region expression parsing
 *
 * Holds reference to CSG system and a surface ID lookup function.
 */
typedef struct {
    alea_system_t* sys;
    arena_t* arena;

    /* Surface lookup: maps OpenMC surface ID -> pre-created sense node pair */
    /* This is populated during conversion */
    openmc_surface_node_pair_t* surface_node_map;
    size_t surface_map_size;

    /* Error state */
    char error_msg[256];
    int error_pos;
} openmc_region_ctx_t;

/**
 * @brief Initialize region parsing context
 */
void openmc_region_ctx_init(openmc_region_ctx_t* ctx, alea_system_t* sys, arena_t* arena);

/**
 * @brief Register a surface ID with its pre-created pos/neg sense nodes
 * @return 0 on success, -1 on error
 */
int openmc_region_register_surface(openmc_region_ctx_t* ctx, int surface_id,
                                   alea_node_id_t pos_node, alea_node_id_t neg_node);

/**
 * @brief Parse a region expression and build CSG tree
 *
 * @param ctx Region context with surface mappings
 * @param region_expr OpenMC region expression string
 * @return Root node ID of CSG tree, or ALEA_NODE_ID_INVALID on error
 */
alea_node_id_t openmc_parse_region(openmc_region_ctx_t* ctx, const char* region_expr);

/**
 * @brief Get error message from failed parse
 */
const char* openmc_region_get_error(const openmc_region_ctx_t* ctx);

#ifdef __cplusplus
}
#endif

#endif /* OPENMC_REGION_H */
