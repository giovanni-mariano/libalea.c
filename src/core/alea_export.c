// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_export.c
 * @brief CSG geometry export pipeline implementation
 *
 * General export infrastructure including:
 * - Export context management
 * - Universe depth computation
 * - Macrobody expansion utilities
 * - Surface ID assignment and deduplication
 *
 * Format-specific export is implemented in:
 * - mcnp/exporter/mcnp_export.c (MCNP format)
 * - openmc/openmc_export.c (OpenMC format)
 */

#include "alea_export.h"
#include "alea_system.h"
#include "core/alea_macrobody.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "util/alea_log.h"

int alea_next_synthetic_surface_id(const alea_system_t* sys) {
    int max_id = 0;
    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        if (sys->surfaces.data[i].mc_surface_id > max_id)
            max_id = sys->surfaces.data[i].mc_surface_id;
    }
    return max_id + 1;
}

/* ============================================================================
 * EXPORT CONTEXT
 * ============================================================================ */

export_context_t* export_context_create(alea_export_format_t format,
                                        alea_surface_emit_policy_t surface_policy,
                                        FILE* out,
                                        bool deduplicate,
                                        int next_synthetic_surface_id,
                                        int universe_depth,
                                        int fill_depth) {
    export_context_t* ctx = calloc(1, sizeof(export_context_t));
    if (!ctx) return NULL;

    ctx->format = format;
    ctx->surface_policy = surface_policy;
    ctx->out = out;
    ctx->deduplicate = deduplicate;
    ctx->universe_depth = universe_depth;
    ctx->fill_depth = fill_depth;
    arena_init(&ctx->arena);

    ctx->next_synthetic_surface_id = next_synthetic_surface_id;

    /* MCNP formatting defaults */
    ctx->mcnp_max_col = 80;
    ctx->mcnp_cont_indent = 5;

    return ctx;
}

void export_context_destroy(export_context_t* ctx) {
    if (!ctx) return;
    free(ctx->prim_to_surface);
    free(ctx->prim_to_surface_inverted);

    if (ctx->decompositions) {
        for (size_t i = 0; i < ctx->decomposition_count; i++) {
            free(ctx->decompositions[i].surfaces);
        }
        free(ctx->decompositions);
    }
    arena_free(&ctx->arena);

    free(ctx);
}

/* ============================================================================
 * UNIVERSE DEPTH COMPUTATION
 * ============================================================================ */

/**
 * @brief Compute the depth of a universe from base universe (0)
 *
 * Depth 0 = base universe (universe 0)
 * Depth 1 = universes directly filled by cells in universe 0
 * Depth N = universes filled by cells in depth N-1 universes
 *
 * Cache values:
 *   -2 = not computed yet
 *   -3 = currently being computed (cycle detection)
 *   -1 = not reachable from base
 *   >= 0 = computed depth
 */
static int compute_universe_depth(const alea_system_t* sys, int universe_id,
                                   int* depth_cache, size_t cache_size) {
    if (universe_id == 0) return 0;

    if (universe_id > 0 && (size_t)universe_id < cache_size) {
        int cached = depth_cache[universe_id];
        if (cached == -3) {
            return -1;  /* Cycle detected */
        }
        if (cached != -2) {
            return cached;
        }
        depth_cache[universe_id] = -3;
    }

    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        const alea_cell_entry_t* cell = &sys->cells.data[i];
        if (cell->fill_universe == universe_id) {
            int parent_depth = compute_universe_depth(sys, cell->universe_id,
                                                       depth_cache, cache_size);
            if (parent_depth >= 0) {
                int depth = parent_depth + 1;
                if (universe_id > 0 && (size_t)universe_id < cache_size) {
                    depth_cache[universe_id] = depth;
                }
                return depth;
            }
        }
    }

    if (universe_id > 0 && (size_t)universe_id < cache_size) {
        depth_cache[universe_id] = -1;
    }
    return -1;
}

bool alea_should_export_cell(const alea_system_t* sys, const alea_cell_entry_t* cell,
                            int max_depth, int* depth_cache, size_t cache_size) {
    if (max_depth < 0) return true;

    int depth = compute_universe_depth(sys, cell->universe_id, depth_cache, cache_size);
    return depth >= 0 && depth <= max_depth;
}

/* ============================================================================
 * MACROBODY HANDLING
 * ============================================================================ */

bool is_macrobody(alea_primitive_type_t type) {
    return alea_is_macrobody(type);
}

/**
 * @brief Find surface entry by pos_node or neg_node
 */
static const alea_surface_entry_t* find_surface_by_node(const alea_system_t* sys,
                                                        alea_node_id_t node_id) {
    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        if (sys->surfaces.data[i].pos_node == node_id ||
            sys->surfaces.data[i].neg_node == node_id) {
            return &sys->surfaces.data[i];
        }
    }
    return NULL;
}

/**
 * @brief Recursively replace macrobody nodes with pre-computed expanded nodes
 */
static alea_node_id_t replace_macrobody_nodes(alea_system_t* sys, alea_node_id_t node_id) {
    if (node_id == ALEA_NODE_ID_INVALID || node_id >= alea_vec_count(&sys->nodes)) {
        return node_id;
    }

    alea_node_t* node = &sys->nodes.data[node_id];
    alea_operation_t op = ALEA_GET_OPERATION(node);

    if (op == ALEA_OP_PRIMITIVE) {
        alea_primitive_type_t ptype = node->primitive.prim_type;

        bool needs_expansion = alea_is_macrobody(ptype);

        /* Check for 1-sheet cone */
        if (!needs_expansion && (ptype == ALEA_PRIMITIVE_CONE_X ||
                                  ptype == ALEA_PRIMITIVE_CONE_Y ||
                                  ptype == ALEA_PRIMITIVE_CONE_Z)) {
            alea_primitive_id_t pid = node->primitive.primitive_id;
            if (pid < alea_vec_count(&sys->primitives)) {
                alea_primitive_data_t pdata;
                if (!alea_primitive_copy_data(sys, pid, &pdata)) {
                    return ALEA_NODE_ID_INVALID;
                }
                switch (ptype) {
                    case ALEA_PRIMITIVE_CONE_X:
                        needs_expansion = (pdata.cone_x.sheet_selection != 0);
                        break;
                    case ALEA_PRIMITIVE_CONE_Y:
                        needs_expansion = (pdata.cone_y.sheet_selection != 0);
                        break;
                    case ALEA_PRIMITIVE_CONE_Z:
                        needs_expansion = (pdata.cone_z.sheet_selection != 0);
                        break;
                    default:
                        break;
                }
            }
        }

        if (!needs_expansion) {
            return node_id;
        }

        const alea_surface_entry_t* surf = find_surface_by_node(sys, node_id);
        if (!surf || surf->expanded_neg_node == ALEA_NODE_ID_INVALID) {
            return node_id;
        }

        int8_t sense = node->primitive.sense;
        if (node->primitive.inverted) {
            sense = -sense;
        }

        return (sense < 0) ? surf->expanded_neg_node : surf->expanded_pos_node;
    }

    if (op == ALEA_OP_COMPLEMENT) {
        alea_node_id_t new_left = replace_macrobody_nodes(sys, node->operation.left);
        if (new_left != node->operation.left) {
            alea_node_id_t new_node = alea_alloc_node(sys);
            if (new_node == ALEA_NODE_ID_INVALID) return ALEA_NODE_ID_INVALID;
            alea_node_t* n = alea_get_node(sys, new_node);
            ALEA_SET_OPERATION(n, ALEA_OP_COMPLEMENT);
            n->operation.left = new_left;
            n->operation.right = ALEA_NODE_ID_INVALID;
            return new_node;
        }
        return node_id;
    }

    if (op == ALEA_OP_UNION || op == ALEA_OP_INTERSECTION || op == ALEA_OP_DIFFERENCE) {
        alea_node_id_t new_left = replace_macrobody_nodes(sys, node->operation.left);
        alea_node_id_t new_right = replace_macrobody_nodes(sys, node->operation.right);

        if (new_left != node->operation.left || new_right != node->operation.right) {
            alea_node_id_t new_node = alea_alloc_node(sys);
            if (new_node == ALEA_NODE_ID_INVALID) return ALEA_NODE_ID_INVALID;
            alea_node_t* n = alea_get_node(sys, new_node);
            ALEA_SET_OPERATION(n, op);
            n->operation.left = new_left;
            n->operation.right = new_right;
            return new_node;
        }
        return node_id;
    }

    return node_id;
}

int alea_expand_macrobodies_tree_level(alea_system_t* sys, export_context_t* ctx) {
    int expanded = 0;

    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        alea_cell_entry_t* cell = &sys->cells.data[i];
        if (cell->root_node_id == ALEA_NODE_ID_INVALID) continue;

        alea_node_id_t new_root = replace_macrobody_nodes(sys, cell->root_node_id);
        if (new_root != cell->root_node_id && new_root != ALEA_NODE_ID_INVALID) {
            cell->root_node_id = new_root;
            expanded++;
        }
    }

    ctx->macrobodies_decomposed = expanded;

    if (expanded > 0) {
        ALEA_LOG_INFO("Tree-level expansion: replaced macrobodies in %d cells using pre-computed nodes", expanded);
    }

    return expanded;
}

int alea_expand_macrobodies_in_cells(alea_system_t* sys) {
    if (!sys) return 0;

    int expanded = 0;

    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        ALEA_CHECK_INTERRUPTED(-1);
        alea_cell_entry_t* cell = &sys->cells.data[i];
        if (cell->root_node_id == ALEA_NODE_ID_INVALID) continue;

        alea_node_id_t new_root = replace_macrobody_nodes(sys, cell->root_node_id);
        if (new_root != cell->root_node_id && new_root != ALEA_NODE_ID_INVALID) {
            cell->root_node_id = new_root;
            expanded++;
        }
    }

    if (expanded > 0) {
        ALEA_LOG_INFO("Expanded macrobodies in %d cells", expanded);
    }

    return expanded;
}

const macrobody_decomposition_t* find_decomposition(const export_context_t* ctx,
                                                     int original_surface_id) {
    if (!ctx || !ctx->decompositions) return NULL;

    for (size_t i = 0; i < ctx->decomposition_count; i++) {
        if (ctx->decompositions[i].original_surface_id == original_surface_id) {
            return &ctx->decompositions[i];
        }
    }
    return NULL;
}

/* ============================================================================
 * SURFACE ID ASSIGNMENT
 * ============================================================================ */

/**
 * @brief Recursively assign surface IDs to primitive nodes that don't have them
 */
static void assign_surface_ids_recursive(alea_system_t* sys, export_context_t* ctx,
                                          alea_node_id_t node_id,
                                          int** prim_surf_map, size_t* map_size) {
    if (node_id == ALEA_NODE_ID_INVALID || node_id >= alea_vec_count(&sys->nodes)) return;

    alea_node_t* node = &sys->nodes.data[node_id];
    alea_operation_t op = ALEA_GET_OPERATION(node);

    if (op == ALEA_OP_PRIMITIVE) {
        if (node->primitive.mc_surface_id == 0) {
            alea_primitive_id_t prim_id = node->primitive.primitive_id;

            if (prim_id >= *map_size) {
                size_t new_size = (prim_id + 1) * 2;
                int* new_map = realloc(*prim_surf_map, new_size * sizeof(int));
                if (!new_map) return;
                for (size_t i = *map_size; i < new_size; i++) {
                    new_map[i] = 0;
                }
                *prim_surf_map = new_map;
                *map_size = new_size;
            }

            if ((*prim_surf_map)[prim_id] != 0) {
                node->primitive.mc_surface_id = (*prim_surf_map)[prim_id];
            } else {
                int new_id = ctx->next_synthetic_surface_id++;
                node->primitive.mc_surface_id = new_id;
                (*prim_surf_map)[prim_id] = new_id;
                ctx->synthetic_surfaces_created++;

                /* Register a surface entry so the canonical surface map
                   and surface emission can find this primitive */
                alea_surface_entry_t* entry = alea_vec_push_uninit(&sys->surfaces, alea_surface_entry_t);
                if (entry) {
                    memset(entry, 0, sizeof(*entry));
                    entry->mc_surface_id = new_id;
                    entry->primitive_id = prim_id;
                    entry->pos_node = node_id;
                    entry->neg_node = ALEA_NODE_ID_INVALID;
                    entry->expanded_pos_node = ALEA_NODE_ID_INVALID;
                    entry->expanded_neg_node = ALEA_NODE_ID_INVALID;
                }

                ALEA_LOG_INFO("Assigned surface ID %d to primitive %u (type %d)",
                       new_id, prim_id, node->primitive.prim_type);
            }
        }
        return;
    }

    if (op == ALEA_OP_COMPLEMENT) {
        assign_surface_ids_recursive(sys, ctx, node->operation.left, prim_surf_map, map_size);
    } else if (op == ALEA_OP_UNION || op == ALEA_OP_INTERSECTION || op == ALEA_OP_DIFFERENCE) {
        assign_surface_ids_recursive(sys, ctx, node->operation.left, prim_surf_map, map_size);
        assign_surface_ids_recursive(sys, ctx, node->operation.right, prim_surf_map, map_size);
    }
}

void alea_assign_missing_surface_ids(alea_system_t* sys, export_context_t* ctx) {
    int* prim_surf_map = NULL;
    size_t map_size = 0;

    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        alea_cell_entry_t* cell = &sys->cells.data[i];
        if (cell->root_node_id == ALEA_NODE_ID_INVALID) continue;

        /* Use original (pre-TRCL) tree when preserving TRCL */
        uint32_t export_root = cell->root_node_id;
        if (cell->original_root_node_id != ALEA_NODE_ID_INVALID
            && ctx->trcl_export_mode == ALEA_TRCL_EXPORT_PRESERVE) {
            export_root = cell->original_root_node_id;
        }

        assign_surface_ids_recursive(sys, ctx, export_root, &prim_surf_map, &map_size);
    }

    free(prim_surf_map);

    if (ctx->synthetic_surfaces_created > 0) {
        ALEA_LOG_INFO("Assigned %zu synthetic surface IDs for expanded macrobodies",
               ctx->synthetic_surfaces_created);
    }
}

/* ============================================================================
 * SURFACE DEDUPLICATION
 * ============================================================================ */

void alea_build_canonical_surface_map(export_context_t* ctx, const alea_system_t* sys) {
    if (!ctx->deduplicate) return;

    uint32_t max_prim_id = 0;
    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        uint32_t prim_id = sys->surfaces.data[i].primitive_id;
        if (prim_id > max_prim_id) {
            max_prim_id = prim_id;
        }
    }

    ctx->prim_to_surface_size = max_prim_id + 1;
    ctx->prim_to_surface = calloc(ctx->prim_to_surface_size, sizeof(int));
    if (!ctx->prim_to_surface) return;

    ctx->prim_to_surface_inverted = calloc(ctx->prim_to_surface_size, sizeof(int8_t));
    if (!ctx->prim_to_surface_inverted) {
        free(ctx->prim_to_surface);
        ctx->prim_to_surface = NULL;
        return;
    }

    for (size_t i = 0; i < ctx->prim_to_surface_size; i++) {
        ctx->prim_to_surface[i] = -1;
    }

    /* Pass 1: pick the lowest MCNP surface ID as canonical for each primitive */
    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        const alea_surface_entry_t* surf = &sys->surfaces.data[i];
        uint32_t prim_id = surf->primitive_id;
        int mcnp_id = surf->mc_surface_id;

        if (ctx->prim_to_surface[prim_id] < 0 ||
            mcnp_id < ctx->prim_to_surface[prim_id]) {
            ctx->prim_to_surface[prim_id] = mcnp_id;

            int8_t inv = 0;
            if (surf->pos_node < alea_vec_count(&sys->nodes)) {
                inv = sys->nodes.data[surf->pos_node].primitive.inverted;
            }
            ctx->prim_to_surface_inverted[prim_id] = inv;
        }
    }

    /* Pass 2: count dedup hits (surfaces that map to a different canonical ID) */
    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        const alea_surface_entry_t* surf = &sys->surfaces.data[i];
        uint32_t prim_id = surf->primitive_id;
        if (ctx->prim_to_surface[prim_id] >= 0 &&
            ctx->prim_to_surface[prim_id] != surf->mc_surface_id) {
            ctx->dedup_hits++;
        }
    }

    ALEA_LOG_INFO("Built canonical surface map: %zu primitives, %zu dedup hits",
           ctx->prim_to_surface_size, ctx->dedup_hits);
}
