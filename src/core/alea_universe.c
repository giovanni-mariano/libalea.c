// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_universe.c
 * @brief Universe instance management and lazy flattening implementation
 */

#include "alea_universe.h"
#include "alea_spatial.h"
#include "core/alea_spatial_hier.h"
#include "core/alea_system.h"
#include "core/alea_ops.h"
#include "core/alea_eval.h"
#include "primitives/bbox.h"
#include "primitives/primitive_desc.h"
#include "core/alea_macrobody.h"
#include "core/alea_transform.h"
#include "raycast/bvh.h"
#include "util/compat.h"
#include "util/alea_bitset.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "util/math.h"
#include "util/alea_log.h"

/* ============================================================================
 * CONSTANTS
 * ============================================================================ */

#define MAX_FLATTEN_DEPTH 100
#define UNIVERSE_POINT_BVH_THRESHOLD 8192
#define UNIVERSE_POINT_BVH_LEAF_SIZE 8

/* Debug trace flag - set via alea_set_debug_point_trace() (thread-local) */
ALEA_THREAD_LOCAL int g_debug_point_trace = 0;

typedef struct {
    uint32_t cell_position;
    alea_bbox_t bbox;
    double centroid[3];
} universe_bvh_item_t;

typedef struct {
    const alea_system_t* sys;
    const alea_universe_t* univ;
    double gx, gy, gz;
    double lx, ly, lz;
    const alea_matrix_t* accumulated_transform;
    int depth;
    alea_size_vec_t* candidate_positions;
    int error;
} universe_point_query_ctx_t;

/* Per-thread accumulator avoids cache-line ping-ponging under OpenMP grid
 * render. The public getter sums across threads (best-effort, no fence). */
static ALEA_THREAD_LOCAL alea_universe_point_bvh_stats_t g_point_bvh_stats;

static size_t universe_point_bvh_threshold(void) {
    const char* env = getenv("ALEA_UNIVERSE_POINT_BVH_THRESHOLD");
    if (env && env[0]) {
        char* end = NULL;
        unsigned long value = strtoul(env, &end, 10);
        if (end != env && value > 0) return (size_t)value;
    }
    return UNIVERSE_POINT_BVH_THRESHOLD;
}

void alea_universe_point_bvh_stats_reset(void) {
    memset(&g_point_bvh_stats, 0, sizeof(g_point_bvh_stats));
}

void alea_universe_point_bvh_stats_get(alea_universe_point_bvh_stats_t* out) {
    if (!out) return;
    *out = g_point_bvh_stats;
}

static int find_all_cells_recursive(const alea_system_t* sys,
                                    double gx, double gy, double gz,
                                    double lx, double ly, double lz,
                                    int universe_id,
                                    const alea_matrix_t* accumulated_transform,
                                    int depth,
                                    alea_cell_hit_t* out_hits,
                                    size_t max_hits,
                                    size_t* hit_count);
static int lattice_rect_lookup(const alea_cell_entry_t* cell,
                               double px, double py, double pz,
                               double* ox, double* oy, double* oz);
int lattice_hex_lookup(const alea_cell_entry_t* cell,
                       double px, double py, double pz,
                       double* ox, double* oy, double* oz);

static alea_bbox_t universe_cell_bbox(alea_system_t* sys, uint32_t cell_index) {
    const alea_cell_entry_t* cell = &sys->cells.data[cell_index];

    if (cell->root_node_id == ALEA_NODE_ID_INVALID) {
        return alea_bbox_empty();
    }

    const alea_node_t* root = &sys->nodes.data[cell->root_node_id];
    if (root->bbox.min_x <= root->bbox.max_x) {
        return root->bbox;
    }

    alea_bbox_t bbox = alea_get_bbox(sys, cell->root_node_id);
    if (bbox.min_x > bbox.max_x ||
        bbox.min_x <= -1e10 || bbox.max_x >= 1e10 ||
        bbox.min_y <= -1e10 || bbox.max_y >= 1e10 ||
        bbox.min_z <= -1e10 || bbox.max_z >= 1e10) {
        bbox.min_x = bbox.min_y = bbox.min_z = -1e6;
        bbox.max_x = bbox.max_y = bbox.max_z = 1e6;
    }

    return bbox;
}

static bool universe_has_lattice_cells(const alea_system_t* sys,
                                       const alea_universe_t* univ) {
    for (size_t i = 0; i < univ->cell_indices.count; i++) {
        size_t cell_idx = univ->cell_indices.data[i];
        const alea_cell_entry_t* cell = &sys->cells.data[cell_idx];
        if (cell->lat_type != 0 && cell->lat_fill) {
            return true;
        }
    }
    return false;
}

static int compare_bvh_item_x(const void* a, const void* b) {
    const universe_bvh_item_t* ia = (const universe_bvh_item_t*)a;
    const universe_bvh_item_t* ib = (const universe_bvh_item_t*)b;
    return (ia->centroid[0] > ib->centroid[0]) - (ia->centroid[0] < ib->centroid[0]);
}

static int compare_bvh_item_y(const void* a, const void* b) {
    const universe_bvh_item_t* ia = (const universe_bvh_item_t*)a;
    const universe_bvh_item_t* ib = (const universe_bvh_item_t*)b;
    return (ia->centroid[1] > ib->centroid[1]) - (ia->centroid[1] < ib->centroid[1]);
}

static int compare_bvh_item_z(const void* a, const void* b) {
    const universe_bvh_item_t* ia = (const universe_bvh_item_t*)a;
    const universe_bvh_item_t* ib = (const universe_bvh_item_t*)b;
    return (ia->centroid[2] > ib->centroid[2]) - (ia->centroid[2] < ib->centroid[2]);
}

static int compare_size_t_ascending(const void* a, const void* b) {
    const size_t* ia = (const size_t*)a;
    const size_t* ib = (const size_t*)b;
    return (*ia > *ib) - (*ia < *ib);
}

static alea_bbox_t universe_bvh_items_bbox(const universe_bvh_item_t* items,
                                           size_t start, size_t end) {
    alea_bbox_t result = alea_bbox_empty();
    for (size_t i = start; i < end; i++) {
        result = alea_bbox_union(&result, &items[i].bbox);
    }
    return result;
}

static int universe_bvh_split_axis(const universe_bvh_item_t* items,
                                   size_t start, size_t end) {
    alea_bbox_t cbox = alea_bbox_empty();
    for (size_t i = start; i < end; i++) {
        double cx = items[i].centroid[0];
        double cy = items[i].centroid[1];
        double cz = items[i].centroid[2];
        if (cx < cbox.min_x) cbox.min_x = cx;
        if (cx > cbox.max_x) cbox.max_x = cx;
        if (cy < cbox.min_y) cbox.min_y = cy;
        if (cy > cbox.max_y) cbox.max_y = cy;
        if (cz < cbox.min_z) cbox.min_z = cz;
        if (cz > cbox.max_z) cbox.max_z = cz;
    }

    double dx = cbox.max_x - cbox.min_x;
    double dy = cbox.max_y - cbox.min_y;
    double dz = cbox.max_z - cbox.min_z;
    if (dx >= dy && dx >= dz) return 0;
    if (dy >= dz) return 1;
    return 2;
}

static int universe_bvh_reserve_node(alea_universe_t* univ) {
    return alea_vec_reserve(&univ->point_bvh_nodes,
                            univ->point_bvh_nodes.count + 1,
                            alea_universe_bvh_node_t);
}

static uint32_t universe_bvh_build_recursive(alea_universe_t* univ,
                                             universe_bvh_item_t* items,
                                             size_t start, size_t end,
                                             int depth) {
    if (universe_bvh_reserve_node(univ) != 0) {
        return UINT32_MAX;
    }

    uint32_t node_idx = (uint32_t)univ->point_bvh_nodes.count++;
    alea_universe_bvh_node_t* node = &univ->point_bvh_nodes.data[node_idx];
    size_t count = end - start;

    if (count <= UNIVERSE_POINT_BVH_LEAF_SIZE || depth > 30) {
        node->bbox = universe_bvh_items_bbox(items, start, end);
        node->left_or_first = (uint32_t)start;
        node->right_child = UINT32_MAX;
        node->count = (uint16_t)count;
        node->axis = 0;
        return node_idx;
    }

    int axis = universe_bvh_split_axis(items, start, end);
    int (*cmp)(const void*, const void*) =
        (axis == 0) ? compare_bvh_item_x :
        (axis == 1) ? compare_bvh_item_y : compare_bvh_item_z;
    qsort(items + start, count, sizeof(universe_bvh_item_t), cmp);

    size_t mid = start + count / 2;
    uint32_t left = universe_bvh_build_recursive(univ, items, start, mid, depth + 1);
    uint32_t right = universe_bvh_build_recursive(univ, items, mid, end, depth + 1);
    if (left == UINT32_MAX || right == UINT32_MAX) {
        return UINT32_MAX;
    }

    node = &univ->point_bvh_nodes.data[node_idx];
    node->bbox = alea_bbox_union(&univ->point_bvh_nodes.data[left].bbox,
                                 &univ->point_bvh_nodes.data[right].bbox);
    node->left_or_first = left;
    node->right_child = right;
    node->count = 0;
    node->axis = (uint8_t)axis;
    return node_idx;
}

/* Forward decl for the eager prebuild. */
static int ensure_universe_point_bvh(alea_system_t* sys, alea_universe_t* univ);

/* Build per-universe point BVHs for every universe sequentially. Called once
 * during query-cache preparation so concurrent point queries (OpenMP grid
 * render, etc.) never trip the lazy build path. */
int alea_prebuild_universe_point_bvhs(alea_system_t* sys) {
    if (!sys) return -1;
    if (!sys->universe_index_built) return 0;
    for (size_t i = 0; i < sys->universes.count; i++) {
        ensure_universe_point_bvh(sys, &sys->universes.data[i]);
        /* Non-fatal: per-universe disabled flag handles failures. */
    }
    return 0;
}

static int ensure_universe_point_bvh(alea_system_t* sys, alea_universe_t* univ) {
    if (!sys || !univ) return -1;
    if (univ->point_bvh_built) return 0;
    if (univ->point_bvh_disabled) return -1;

    /* Serialize the lazy build across threads (OpenMP grid render in the
     * plotter, or any concurrent call). The lock-free fast path above
     * handles the already-built case. Using a plain mutex avoids a
     * link-time dependency on libgomp in non-OpenMP builds. */
    static alea_mutex_t build_mutex = ALEA_MUTEX_INIT;
    alea_mutex_lock(&build_mutex);
    {
    if (univ->point_bvh_built || univ->point_bvh_disabled) goto done;

    size_t cell_count = univ->cell_indices.count;
    if (cell_count <= universe_point_bvh_threshold() ||
        universe_has_lattice_cells(sys, univ)) {
        univ->point_bvh_disabled = true;
        goto done;
    }

    universe_bvh_item_t* items = malloc(cell_count * sizeof(universe_bvh_item_t));
    if (!items) {
        univ->point_bvh_disabled = true;
        goto done;
    }

    for (size_t i = 0; i < cell_count; i++) {
        uint32_t cell_idx = (uint32_t)univ->cell_indices.data[i];
        alea_bbox_t bbox = universe_cell_bbox(sys, cell_idx);
        items[i].cell_position = (uint32_t)i;
        items[i].bbox = bbox;
        items[i].centroid[0] = (bbox.min_x + bbox.max_x) * 0.5;
        items[i].centroid[1] = (bbox.min_y + bbox.max_y) * 0.5;
        items[i].centroid[2] = (bbox.min_z + bbox.max_z) * 0.5;
    }

    alea_vec_init(&univ->point_bvh_nodes);
    uint32_t root = universe_bvh_build_recursive(univ, items, 0, cell_count, 0);
    if (root == UINT32_MAX) {
        free(items);
        alea_vec_free(&univ->point_bvh_nodes);
        univ->point_bvh_disabled = true;
        goto done;
    }

    univ->point_bvh_indices = malloc(cell_count * sizeof(uint32_t));
    if (!univ->point_bvh_indices) {
        free(items);
        alea_vec_free(&univ->point_bvh_nodes);
        univ->point_bvh_disabled = true;
        goto done;
    }

    for (size_t i = 0; i < cell_count; i++) {
        univ->point_bvh_indices[i] = items[i].cell_position;
    }

    free(items);
    univ->point_bvh_built = true;
    g_point_bvh_stats.bvh_builds++;
    done: ;
    } /* end critical section */
    alea_mutex_unlock(&build_mutex);
    return univ->point_bvh_built ? 0 : -1;
}

static int process_cell_for_all_cells_query(const alea_system_t* sys,
                                            size_t cell_idx,
                                            double gx, double gy, double gz,
                                            double lx, double ly, double lz,
                                            int universe_id,
                                            const alea_matrix_t* accumulated_transform,
                                            int depth,
                                            alea_cell_hit_t* out_hits,
                                            size_t max_hits,
                                            size_t* hit_count) {
    const alea_cell_entry_t* cell = &sys->cells.data[cell_idx];

    if (cell->lat_type != 0 && cell->lat_fill) {
        double ox, oy, oz;
        int fill_univ = (cell->lat_type == 2)
            ? lattice_hex_lookup(cell, lx, ly, lz, &ox, &oy, &oz)
            : lattice_rect_lookup(cell, lx, ly, lz, &ox, &oy, &oz);
        if (fill_univ < 0) return 0;

        if (*hit_count < max_hits) {
            alea_cell_hit_t* hit = &out_hits[*hit_count];
            hit->cell_id = cell->mc_cell_id;
            hit->cell_index = (int)cell_idx;
            hit->material_id = cell->material_id;
            hit->universe_id = cell->universe_id;
            hit->fill_universe = fill_univ;
            hit->depth = depth;
            hit->local_x = lx;
            hit->local_y = ly;
            hit->local_z = lz;
            (*hit_count)++;
        }

        double elx = lx - ox, ely = ly - oy, elz = lz - oz;
        return find_all_cells_recursive(sys, gx, gy, gz,
                                        elx, ely, elz,
                                        fill_univ, NULL,
                                        depth + 1,
                                        out_hits, max_hits, hit_count);
    }

    g_point_bvh_stats.exact_cell_tests++;
    if (!alea_contains_point(sys, cell->root_node_id, lx, ly, lz)) {
        return 0;
    }

    if (*hit_count < max_hits) {
        alea_cell_hit_t* hit = &out_hits[*hit_count];
        hit->cell_id = cell->mc_cell_id;
        hit->cell_index = (int)cell_idx;
        hit->material_id = cell->material_id;
        hit->universe_id = cell->universe_id;
        hit->fill_universe = cell->fill_universe;
        hit->depth = depth;
        hit->local_x = lx;
        hit->local_y = ly;
        hit->local_z = lz;
        (*hit_count)++;

        if (g_debug_point_trace) {
            ALEA_LOG_DEBUG("  -> Found cell %d (mat=%d) in universe %d, fill=%d",
                   cell->mc_cell_id, cell->material_id, universe_id, cell->fill_universe);
        }
    }

    if (cell->fill_universe > 0 && *hit_count < max_hits) {
        alea_matrix_t fill_transform;
        if (cell->fill_transform > 0) {
            const alea_transform_t* tr = alea_get_transform(sys, cell->fill_transform);
            if (tr) {
                if (!alea_matrix_from_mcnp(&fill_transform, tr->cosines,
                                           tr->value_count, false)) {
                    return -1;
                }
            } else {
                return -1;
            }
        } else {
            alea_matrix_identity(&fill_transform);
        }

        alea_matrix_t new_accumulated;
        if (accumulated_transform) {
            alea_matrix_multiply(&new_accumulated, accumulated_transform, &fill_transform);
        } else {
            new_accumulated = fill_transform;
        }

        double fill_x = gx, fill_y = gy, fill_z = gz;
        if (!alea_matrix_invert(&new_accumulated)) return -1;
        alea_matrix_transform_point_inverse(&new_accumulated, &fill_x, &fill_y, &fill_z);

        if (g_debug_point_trace) {
            ALEA_LOG_DEBUG("  -> Entering fill %d (transform=%d)",
                   cell->fill_universe, cell->fill_transform);
        }

        return find_all_cells_recursive(sys, gx, gy, gz,
                                        fill_x, fill_y, fill_z,
                                        cell->fill_universe,
                                        &new_accumulated,
                                        depth + 1,
                                        out_hits, max_hits, hit_count);
    }

    return 0;
}

static void traverse_universe_point_bvh_node(universe_point_query_ctx_t* ctx,
                                             uint32_t node_idx) {
    if (ctx->error) return;
    g_point_bvh_stats.bvh_node_visits++;
    const alea_universe_bvh_node_t* node = &ctx->univ->point_bvh_nodes.data[node_idx];

    g_point_bvh_stats.bvh_bbox_tests++;
    if (!alea_bbox_contains_point(&node->bbox, ctx->lx, ctx->ly, ctx->lz)) {
        return;
    }

    if (node->count > 0) {
        for (uint16_t i = 0; i < node->count; i++) {
            uint32_t cell_pos = ctx->univ->point_bvh_indices[node->left_or_first + i];
            uint32_t cell_idx = (uint32_t)ctx->univ->cell_indices.data[cell_pos];
            alea_bbox_t bbox = universe_cell_bbox((alea_system_t*)ctx->sys, cell_idx);
            g_point_bvh_stats.bvh_bbox_tests++;
            if (!alea_bbox_contains_point(&bbox, ctx->lx, ctx->ly, ctx->lz)) {
                continue;
            }
            g_point_bvh_stats.bvh_candidates++;
            if (alea_vec_push(ctx->candidate_positions, (size_t)cell_pos, size_t) != 0) {
                ctx->error = -1;
                return;
            }
        }
        return;
    }

    traverse_universe_point_bvh_node(ctx, node->left_or_first);
    traverse_universe_point_bvh_node(ctx, node->right_child);
}

/**
 * @brief Build primitive->surface mapping after flatten
 * 
 * After flattening, primitives have new IDs due to deduplication and
 * transformation. Original MCNP surface IDs are lost (set to 0 in cloned nodes).
 * 
 * This function:
 * 1. Scans all cells to find which primitives are actually used
 * 2. Assigns NEW sequential surface IDs to each used primitive
 * 3. Updates nodes with the new surface IDs
 * 
 * @param sys Flattened system
 * @return 0 on success, -1 on error
 */
static int build_primitive_to_surface_map_from_nodes(alea_system_t* sys) {
    if (!sys) return -1;

    /* Find max primitive ID */
    uint32_t max_prim_id = 0;
    if (alea_vec_count(&sys->primitives) > 0) {
        max_prim_id = (uint32_t)(alea_vec_count(&sys->primitives) - 1);
    }

    /* Track which primitives are used + preferred surface IDs from existing nodes */
    alea_bitset_t prim_used = alea_bitset_create(max_prim_id + 1);
    int* preferred = calloc(max_prim_id + 1, sizeof(int));
    if (!prim_used.words || !preferred) {
        alea_bitset_destroy(&prim_used);
        free(preferred);
        return -1;
    }

    /* Scan all cell trees to find used primitives and collect preferred surface IDs */
    for (size_t cell_idx = 0; cell_idx < alea_vec_count(&sys->cells); cell_idx++) {
        alea_node_id_t root = sys->cells.data[cell_idx].root_node_id;

        /* Traverse tree (stack-based traversal with dynamic stack) */
        alea_uint32_vec_t stack = ALEA_VEC_INIT;
        alea_vec_push(&stack, root, alea_node_id_t);

        while (alea_vec_count(&stack) > 0) {
            alea_node_id_t node_id = alea_vec_pop(&stack);
            if (node_id == ALEA_NODE_ID_INVALID || node_id >= alea_vec_count(&sys->nodes)) continue;

            const alea_node_t* node = &sys->nodes.data[node_id];
            alea_operation_t op = ALEA_GET_OPERATION(node);

            if (op == ALEA_OP_PRIMITIVE) {
                uint32_t prim_id = node->primitive.primitive_id;
                if (prim_id <= max_prim_id) {
                    alea_bitset_set(&prim_used, prim_id);
                    /* Record preferred surface ID from original mc_surface_id if set */
                    if (node->primitive.mc_surface_id > 0 && preferred[prim_id] == 0) {
                        preferred[prim_id] = node->primitive.mc_surface_id;
                    }
                }
            } else {
                if (op != ALEA_OP_COMPLEMENT) {
                    alea_vec_push(&stack, node->operation.right, alea_node_id_t);
                }
                alea_vec_push(&stack, node->operation.left, alea_node_id_t);
            }
        }
        alea_vec_free(&stack);
    }

    /* Count used primitives */
    size_t used_count = alea_bitset_popcount(&prim_used);

    if (used_count == 0) {
        alea_bitset_destroy(&prim_used);
        free(preferred);
        return 0;
    }

    /* Create prim_id -> new surface_id mapping */
    int* prim_to_surf = malloc((max_prim_id + 1) * sizeof(int));
    if (!prim_to_surf) {
        alea_bitset_destroy(&prim_used);
        free(preferred);
        return -1;
    }
    for (uint32_t i = 0; i <= max_prim_id; i++) {
        prim_to_surf[i] = -1;
    }

    /* Find max possible surface ID for the used-ID bitset */
    int max_possible_id = (int)used_count + 1;
    for (uint32_t i = 0; i <= max_prim_id; i++) {
        if (alea_bitset_test(&prim_used, i) && preferred[i] > max_possible_id)
            max_possible_id = preferred[i];
    }

    alea_bitset_t surf_id_used = alea_bitset_create((size_t)max_possible_id + 1);
    if (!surf_id_used.words) {
        alea_bitset_destroy(&prim_used);
        free(preferred);
        free(prim_to_surf);
        return -1;
    }

    /* First pass: assign preferred surface IDs where available and not colliding */
    for (uint32_t i = 0; i <= max_prim_id; i++) {
        if (alea_bitset_test(&prim_used, i) && preferred[i] > 0 && preferred[i] <= max_possible_id
            && !alea_bitset_test(&surf_id_used, preferred[i])) {
            prim_to_surf[i] = preferred[i];
            alea_bitset_set(&surf_id_used, preferred[i]);
        }
    }

    /* Second pass: auto-assign remaining */
    int next_surf_id = 1;
    for (uint32_t i = 0; i <= max_prim_id; i++) {
        if (alea_bitset_test(&prim_used, i) && prim_to_surf[i] < 0) {
            while (next_surf_id <= max_possible_id && alea_bitset_test(&surf_id_used, next_surf_id))
                next_surf_id++;
            prim_to_surf[i] = next_surf_id;
            if (next_surf_id <= max_possible_id)
                alea_bitset_set(&surf_id_used, next_surf_id);
            next_surf_id++;
        }
    }
    alea_bitset_destroy(&surf_id_used);
    free(preferred);

    /* Find actual max assigned ID for lookup table sizing */
    int max_assigned_id = 0;
    for (uint32_t i = 0; i <= max_prim_id; i++) {
        if (prim_to_surf[i] > max_assigned_id)
            max_assigned_id = prim_to_surf[i];
    }
    
    /* Clear old surfaces */
    alea_vec_clear(&sys->surfaces);

    /* Reserve capacity for surfaces */
    int reserve_res = alea_vec_reserve(&sys->surfaces, used_count, alea_surface_entry_t);
    if (reserve_res != 0) {
        alea_bitset_destroy(&prim_used);
        free(prim_to_surf);
        return -1;
    }
    
    /* First pass: create surface entries and find existing sense nodes */
    /* Track pos/neg nodes per primitive */
    alea_node_id_t* found_pos = calloc(max_prim_id + 1, sizeof(alea_node_id_t));
    alea_node_id_t* found_neg = calloc(max_prim_id + 1, sizeof(alea_node_id_t));
    if (!found_pos || !found_neg) {
        alea_bitset_destroy(&prim_used);
        free(prim_to_surf);
        free(found_pos);
        free(found_neg);
        return -1;
    }
    for (uint32_t i = 0; i <= max_prim_id; i++) {
        found_pos[i] = ALEA_NODE_ID_INVALID;
        found_neg[i] = ALEA_NODE_ID_INVALID;
    }

    /* Track inverted flag per primitive (from existing nodes) */
    int8_t* prim_inverted = calloc(max_prim_id + 1, sizeof(int8_t));
    if (!prim_inverted) {
        alea_bitset_destroy(&prim_used);
        free(prim_to_surf);
        free(found_pos);
        free(found_neg);
        return -1;
    }

    /* Scan existing nodes to find pos/neg sense nodes for each primitive */
    for (size_t i = 0; i < alea_vec_count(&sys->nodes); i++) {
        alea_node_t* node = &sys->nodes.data[i];
        if (ALEA_GET_OPERATION(node) == ALEA_OP_PRIMITIVE) {
            uint32_t prim_id = node->primitive.primitive_id;
            if (prim_id <= max_prim_id && alea_bitset_test(&prim_used, prim_id)) {
                /* Record inverted flag from existing node */
                prim_inverted[prim_id] = node->primitive.inverted;
                if (node->primitive.sense > 0) {
                    if (found_pos[prim_id] == ALEA_NODE_ID_INVALID) {
                        found_pos[prim_id] = (alea_node_id_t)i;
                    }
                } else {
                    if (found_neg[prim_id] == ALEA_NODE_ID_INVALID) {
                        found_neg[prim_id] = (alea_node_id_t)i;
                    }
                }
            }
        }
    }

    /* Create surface entries for each used primitive, creating missing sense nodes */
    for (uint32_t prim_id = 0; prim_id <= max_prim_id; prim_id++) {
        if (!alea_bitset_test(&prim_used, prim_id)) continue;

        int surf_id = prim_to_surf[prim_id];
        int8_t inverted = prim_inverted[prim_id];

        alea_node_id_t pos_node = found_pos[prim_id];
        alea_node_id_t neg_node = found_neg[prim_id];

        /* Create POSITIVE sense node if missing */
        if (pos_node == ALEA_NODE_ID_INVALID) {
            pos_node = alea_add_primitive_node(sys, prim_id, +1, inverted, surf_id);
            if (pos_node == ALEA_NODE_ID_INVALID) {
                alea_bitset_destroy(&prim_used);
                free(prim_to_surf);
                free(found_pos);
                free(found_neg);
                free(prim_inverted);
                return -1;
            }
        }

        /* Create NEGATIVE sense node if missing */
        if (neg_node == ALEA_NODE_ID_INVALID) {
            neg_node = alea_add_primitive_node(sys, prim_id, -1, inverted, surf_id);
            if (neg_node == ALEA_NODE_ID_INVALID) {
                alea_bitset_destroy(&prim_used);
                free(prim_to_surf);
                free(found_pos);
                free(found_neg);
                free(prim_inverted);
                return -1;
            }
        }

        /* Create surface entry using vector API */
        alea_surface_entry_t* surf = alea_vec_push_uninit(&sys->surfaces, alea_surface_entry_t);
        if (!surf) continue;
        memset(surf, 0, sizeof(alea_surface_entry_t));
        surf->mc_surface_id = surf_id;
        surf->primitive_id = prim_id;
        surf->pos_node = pos_node;
        surf->neg_node = neg_node;
        surf->transform_id = 0;
        surf->transform_applied = false;
        surf->boundary_type = ALEA_BOUNDARY_TRANSMISSIVE;
        surf->periodic_surface_id = 0;
    }

    free(found_pos);
    free(found_neg);
    free(prim_inverted);

    /* Update existing nodes with their new surface IDs */
    for (size_t i = 0; i < alea_vec_count(&sys->nodes); i++) {
        alea_node_t* node = &sys->nodes.data[i];
        if (ALEA_GET_OPERATION(node) == ALEA_OP_PRIMITIVE) {
            uint32_t prim_id = node->primitive.primitive_id;
            if (prim_id <= max_prim_id && prim_to_surf[prim_id] > 0) {
                node->primitive.mc_surface_id = prim_to_surf[prim_id];
            }
        }
    }
    
    alea_bitset_destroy(&prim_used);
    free(prim_to_surf);
    
    /* Rebuild surface lookup */
    free(sys->surface_lookup);
    sys->surface_lookup = NULL;

    sys->surface_lookup = malloc((size_t)(max_assigned_id + 1) * sizeof(alea_node_id_t));
    if (!sys->surface_lookup) return -1;

    sys->surface_lookup_size = (size_t)(max_assigned_id + 1);

    for (size_t i = 0; i < sys->surface_lookup_size; i++) {
        sys->surface_lookup[i] = ALEA_NODE_ID_INVALID;
    }

    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        int id = sys->surfaces.data[i].mc_surface_id;
        if (id >= 0 && (size_t)id < sys->surface_lookup_size)
            sys->surface_lookup[id] = (alea_node_id_t)i;
    }

    /* Rebuild macrobody / 1-sheet-cone immediate expansions.
     *
     * For each surface whose primitive is a macrobody (BOX, SPH, TRC, ELL,
     * REC, WED, RHP, ARB) or a 1-sheet cone, recompute the half-space
     * decomposition into expanded_neg_node / expanded_pos_node. The MCNP
     * exporter uses these at write time to swap a macrobody primitive node
     * for its CSG expansion (replace_macrobody_nodes in alea_export.c).
     *
     * Without this step, surface entries created here have
     * expanded_*_node = INVALID, so the exporter falls back to emitting
     * the macrobody surface card directly — which works for RPP and RCC
     * but not for the other macrobody types (the MCNP write_mcnp_surface
     * switch has no case for them) and produces "Unknown primitive type
     * %d" warnings. Building expansions here lets the exporter emit those
     * macrobodies as their component half-spaces. */
    /* Snapshot count up front: alea_expand_macrobody_immediate adds new
     * component primitives (and nodes) via alea_get_or_create_primitive,
     * which can realloc sys->primitives. Iterate using the original surface
     * count so we don't re-expand entries that didn't exist when we started. */
    size_t surface_count_snapshot = alea_vec_count(&sys->surfaces);
    for (size_t i = 0; i < surface_count_snapshot; i++) {
        /* Re-fetch the surface entry pointer through the index each iteration
         * (we don't add surface entries here, but be defensive). Copy the
         * primitive type and data BY VALUE: the expansion below grows
         * primitives, which invalidates any pointer into the primitives
         * vector. */
        uint32_t prim_id = sys->surfaces.data[i].primitive_id;
        if (prim_id >= alea_vec_count(&sys->primitives)) continue;
        alea_primitive_type_t prim_type = sys->primitives.data[prim_id].type;
        alea_primitive_data_t prim_data = sys->primitives.data[prim_id].data;

        alea_node_id_t exp_neg = ALEA_NODE_ID_INVALID;
        alea_node_id_t exp_pos = ALEA_NODE_ID_INVALID;

        if (alea_is_macrobody(prim_type)) {
            if (alea_expand_macrobody_immediate(sys, prim_type, &prim_data,
                                                &exp_neg, &exp_pos) == 0) {
                sys->surfaces.data[i].expanded_neg_node = exp_neg;
                sys->surfaces.data[i].expanded_pos_node = exp_pos;
            }
            continue;
        }

        /* 1-sheet cones: sheet_selection != 0 in the cone data. */
        bool is_1sheet_cone = false;
        switch (prim_type) {
            case ALEA_PRIMITIVE_CONE_X:
                is_1sheet_cone = (prim_data.cone_x.sheet_selection != 0);
                break;
            case ALEA_PRIMITIVE_CONE_Y:
                is_1sheet_cone = (prim_data.cone_y.sheet_selection != 0);
                break;
            case ALEA_PRIMITIVE_CONE_Z:
                is_1sheet_cone = (prim_data.cone_z.sheet_selection != 0);
                break;
            default:
                break;
        }
        if (is_1sheet_cone) {
            if (alea_expand_1sheet_cone_immediate(sys, prim_type, &prim_data,
                                                  &exp_neg, &exp_pos) == 0) {
                sys->surfaces.data[i].expanded_neg_node = exp_neg;
                sys->surfaces.data[i].expanded_pos_node = exp_pos;
            }
        }
    }

    ALEA_LOG_INFO("Built surface mapping for flattened system: %zu surfaces (max ID %d)",
           alea_vec_count(&sys->surfaces), max_assigned_id);

    return 0;
}
/* ============================================================================
 * PRIMITIVE REMAPPING TABLE
 * ============================================================================ */



primitive_remap_t* alea_create_remap_table(size_t src_prim_count) {
    primitive_remap_t* remap = malloc(sizeof(primitive_remap_t));
    if (!remap) return NULL;

    remap->capacity = src_prim_count;
    remap->map = malloc(src_prim_count * sizeof(uint32_t));
    if (!remap->map) {
        free(remap);
        return NULL;
    }
    remap->inverted = calloc(src_prim_count, sizeof(int8_t));
    if (!remap->inverted) {
        free(remap->map);
        free(remap);
        return NULL;
    }

    /* Initialize map to invalid */
    for (size_t i = 0; i < src_prim_count; i++) {
        remap->map[i] = ALEA_PRIMITIVE_ID_INVALID;
    }

    return remap;
}

void alea_destroy_remap_table(primitive_remap_t* remap) {
    if (remap) {
        free(remap->map);
        free(remap->inverted);
        free(remap);
    }
}

static uint32_t get_remapped_primitive(primitive_remap_t* remap,
                                       uint32_t src_id) {
    if (!remap || src_id >= remap->capacity) {
        return ALEA_PRIMITIVE_ID_INVALID;
    }
    return remap->map[src_id];
}


/* ============================================================================
 * PRIMITIVE COPYING
 * ============================================================================ */

/* ============================================================================
 * SURFACE COPYING WITH REMAPPING
 * ============================================================================ */

/**
 * @brief Copy surfaces from src to dst, remapping primitive IDs
 * 
 * Uses the primitive remap table to update primitive_id references
 * in surface entries so they point to the correct primitives in dst.
 */
int alea_copy_surfaces_with_remap(alea_system_t* dst,
                                 const alea_system_t* src,
                                 primitive_remap_t* remap) {
    if (alea_vec_empty(&src->surfaces)) {
        return 0;  /* No surfaces to copy */
    }

    /* Reserve capacity for surfaces */
    int reserve_res = alea_vec_reserve(&dst->surfaces, alea_vec_count(&src->surfaces), alea_surface_entry_t);
    if (reserve_res != 0) {
        return -1;
    }

    /* Copy each surface with remapped primitive ID */
    for (size_t i = 0; i < alea_vec_count(&src->surfaces); i++) {
        const alea_surface_entry_t* src_surf = &src->surfaces.data[i];

        /* CRITICAL: Remap primitive ID */
        uint32_t new_prim_id = get_remapped_primitive(remap, src_surf->primitive_id);

        if (new_prim_id == ALEA_PRIMITIVE_ID_INVALID) {
            /* Primitive wasn't used in flattened geometry - skip this surface */
            ALEA_LOG_WARN("Surface %d references unused primitive %u, skipping",
                    src_surf->mc_surface_id, src_surf->primitive_id);
            continue;
        }

        /* Allocate new surface entry */
        alea_surface_entry_t* dst_surf = alea_vec_push_uninit(&dst->surfaces, alea_surface_entry_t);
        if (!dst_surf) return -1;

        /* Copy surface data */
        *dst_surf = *src_surf;
        dst_surf->primitive_id = new_prim_id;

        /* Note: pos_node and neg_node will need to be recreated if used */
        /* For now, these are invalid in the dst system */
        dst_surf->pos_node = ALEA_NODE_ID_INVALID;
        dst_surf->neg_node = ALEA_NODE_ID_INVALID;
    }

    return 0;
}

/**
 * @brief Rebuild surface lookup table
 * 
 * Creates the MCNP surface ID -> internal index mapping.
 */
/* ============================================================================
 * MATERIAL COPYING UTILITIES
 * ============================================================================ */

/**
 * @brief Deep copy a single nuclide
 */
static int copy_nuclide(alea_nuclide_t* dst, const alea_nuclide_t* src) {
    dst->zaid = src->zaid;
    dst->fraction = src->fraction;
    
    if (src->library) {
        dst->library = alea_strdup(src->library);
        if (!dst->library) return -1;
    } else {
        dst->library = NULL;
    }
    
    return 0;
}

/**
 * @brief Deep copy a thermal law
 */
static int copy_thermal_law(alea_thermal_law_t* dst, const alea_thermal_law_t* src) {
    dst->zaid_match = src->zaid_match;
    
    if (src->identifier) {
        dst->identifier = alea_strdup(src->identifier);
        if (!dst->identifier) return -1;
    } else {
        dst->identifier = NULL;
    }
    
    return 0;
}

/**
 * @brief Deep copy a material
 */
static int copy_material(alea_material_t* dst, const alea_material_t* src) {
    memset(dst, 0, sizeof(*dst));
    
    dst->material_id = src->material_id;
    dst->is_weight_fraction = src->is_weight_fraction;
    dst->avg_atomic_mass = src->avg_atomic_mass;
    dst->properties_valid = src->properties_valid;
    
    /* Copy name */
    if (src->name) {
        dst->name = alea_strdup(src->name);
        if (!dst->name) return -1;
    }
    
    /* Copy comments */
    if (src->comments) {
        dst->comments = alea_strdup(src->comments);
        if (!dst->comments) {
            free(dst->name);
            return -1;
        }
    }
    
    /* Copy nuclides */
    if (alea_vec_count(&src->nuclides) > 0) {
        int r = alea_vec_reserve(&dst->nuclides, alea_vec_count(&src->nuclides), alea_nuclide_t);
        if (r != 0) {
            free(dst->name);
            free(dst->comments);
            return -1;
        }

        for (size_t i = 0; i < alea_vec_count(&src->nuclides); i++) {
            alea_nuclide_t* nuc = alea_vec_push_uninit(&dst->nuclides, alea_nuclide_t);
            if (copy_nuclide(nuc, &src->nuclides.data[i]) != 0) {
                /* Cleanup on error */
                alea_vec_pop_discard(&dst->nuclides);
                for (size_t j = 0; j < alea_vec_count(&dst->nuclides); j++) {
                    free(dst->nuclides.data[j].library);
                }
                alea_vec_free(&dst->nuclides);
                free(dst->name);
                free(dst->comments);
                return -1;
            }
        }
    }

    /* Copy thermal laws */
    if (alea_vec_count(&src->thermal_laws) > 0) {
        int r = alea_vec_reserve(&dst->thermal_laws, alea_vec_count(&src->thermal_laws), alea_thermal_law_t);
        if (r != 0) {
            for (size_t i = 0; i < alea_vec_count(&dst->nuclides); i++) {
                free(dst->nuclides.data[i].library);
            }
            alea_vec_free(&dst->nuclides);
            free(dst->name);
            free(dst->comments);
            return -1;
        }

        for (size_t i = 0; i < alea_vec_count(&src->thermal_laws); i++) {
            alea_thermal_law_t* law = alea_vec_push_uninit(&dst->thermal_laws, alea_thermal_law_t);
            if (copy_thermal_law(law, &src->thermal_laws.data[i]) != 0) {
                alea_vec_pop_discard(&dst->thermal_laws);
                for (size_t j = 0; j < alea_vec_count(&dst->thermal_laws); j++) {
                    free(dst->thermal_laws.data[j].identifier);
                }
                alea_vec_free(&dst->thermal_laws);
                for (size_t j = 0; j < alea_vec_count(&dst->nuclides); j++) {
                    free(dst->nuclides.data[j].library);
                }
                alea_vec_free(&dst->nuclides);
                free(dst->name);
                free(dst->comments);
                return -1;
            }
        }
    }

    return 0;
}

/**
 * @brief Copy all materials from source to destination system
 */
static void free_material_internals(alea_material_t* m) {
    for (size_t k = 0; k < alea_vec_count(&m->nuclides); k++) {
        free(m->nuclides.data[k].library);
    }
    alea_vec_free(&m->nuclides);
    for (size_t k = 0; k < alea_vec_count(&m->elements); k++) {
        free(m->elements.data[k].library);
    }
    alea_vec_free(&m->elements);
    for (size_t k = 0; k < alea_vec_count(&m->thermal_laws); k++) {
        free(m->thermal_laws.data[k].identifier);
    }
    alea_vec_free(&m->thermal_laws);
    free(m->name);
    free(m->comments);
}

static int copy_materials(alea_system_t* dst, const alea_system_t* src) {
    if (alea_vec_empty(&src->materials)) {
        return 0;  /* No materials to copy */
    }

    /* Copy each material */
    for (size_t i = 0; i < alea_vec_count(&src->materials); i++) {
        alea_material_t* mat = alea_vec_push_uninit(&dst->materials, alea_material_t);
        if (!mat) {
            goto cleanup;
        }
        if (copy_material(mat, &src->materials.data[i]) != 0) {
            alea_vec_pop_discard(&dst->materials);
            goto cleanup;
        }
    }

    return 0;

cleanup:
    for (size_t j = 0; j < alea_vec_count(&dst->materials); j++) {
        free_material_internals(&dst->materials.data[j]);
    }
    alea_vec_free(&dst->materials);
    return -1;
}

/* ============================================================================
 * SURFACE COPYING (if needed)
 * 
 * Note: After flattening, surfaces as separate entities might not be needed
 * since we work with transformed primitives directly. However, if MCNP export
 * is needed, we may want to preserve surface information.
 * ============================================================================ */

/* ============================================================================
 * SELECTIVE COPY HELPERS (for extract_universe / extract_region)
 * ============================================================================ */

int alea_copy_referenced_materials(alea_system_t* dst, const alea_system_t* src) {
    if (!dst || !src) return -1;
    if (alea_vec_empty(&src->materials)) return 0;

    /* Collect unique material_ids from dst cells */
    size_t cell_count = alea_vec_count(&dst->cells);
    for (size_t i = 0; i < cell_count; i++) {
        int mat_id = dst->cells.data[i].material_id;
        if (mat_id == 0) continue;  /* void */

        /* Check if already copied */
        bool already_copied = false;
        for (size_t j = 0; j < alea_vec_count(&dst->materials); j++) {
            if (dst->materials.data[j].material_id == mat_id) {
                already_copied = true;
                break;
            }
        }
        if (already_copied) continue;

        /* Find in source */
        for (size_t j = 0; j < alea_vec_count(&src->materials); j++) {
            if (src->materials.data[j].material_id == mat_id) {
                alea_material_t* mat = alea_vec_push_uninit(&dst->materials, alea_material_t);
                if (!mat) return -1;
                if (copy_material(mat, &src->materials.data[j]) != 0) {
                    alea_vec_pop_discard(&dst->materials);
                    return -1;
                }
                break;
            }
        }
    }
    return 0;
}

int alea_copy_referenced_mixtures(alea_system_t* dst, const alea_system_t* src) {
    if (!dst || !src) return -1;
    if (alea_vec_empty(&src->mixtures)) return 0;

    /* Collect material IDs from dst cells */
    size_t cell_count = alea_vec_count(&dst->cells);

    for (size_t i = 0; i < alea_vec_count(&src->mixtures); i++) {
        const alea_mixture_t* mix = &src->mixtures.data[i];

        /* Check if this mixture's mc_material_id matches any cell material_id */
        bool needed = false;
        for (size_t c = 0; c < cell_count; c++) {
            if (dst->cells.data[c].material_id == mix->mc_material_id) {
                needed = true;
                break;
            }
        }
        if (!needed) continue;

        /* Shallow copy (components reference materials by ID) */
        alea_mixture_t* dst_mix = alea_vec_push_uninit(&dst->mixtures, alea_mixture_t);
        if (!dst_mix) return -1;
        *dst_mix = *mix;

        /* Deep copy allocated fields */
        alea_vec_init(&dst_mix->components);
        if (alea_vec_count(&mix->components) > 0) {
            int r = alea_vec_reserve(&dst_mix->components, alea_vec_count(&mix->components), alea_mixture_comp_t);
            if (r != 0) {
                alea_vec_pop_discard(&dst->mixtures);
                return -1;
            }
            memcpy(dst_mix->components.data, mix->components.data,
                   alea_vec_count(&mix->components) * sizeof(alea_mixture_comp_t));
            dst_mix->components.count = alea_vec_count(&mix->components);
        }
        if (mix->name) {
            dst_mix->name = alea_strdup(mix->name);
            if (!dst_mix->name) return -1;
        }
        if (mix->comments) {
            dst_mix->comments = alea_strdup(mix->comments);
            if (!dst_mix->comments) return -1;
        }
    }
    return 0;
}

int alea_copy_referenced_transforms(alea_system_t* dst, const alea_system_t* src) {
    if (!dst || !src) return -1;
    if (alea_vec_empty(&src->transforms)) return 0;

    /* Collect transform IDs needed from dst cells and surfaces */
    /* Use a simple array to track which transform IDs we need */
    size_t max_needed = 256;
    int* needed_ids = malloc(max_needed * sizeof(int));
    if (!needed_ids) return -1;
    size_t needed_count = 0;

    /* Helper to add a unique ID to the needed list */
    #define ADD_NEEDED_ID(id) do { \
        bool dup = false; \
        for (size_t _k = 0; _k < needed_count; _k++) { \
            if (needed_ids[_k] == (id)) { dup = true; break; } \
        } \
        if (!dup) { \
            if (needed_count >= max_needed) { \
                max_needed *= 2; \
                int* tmp = realloc(needed_ids, max_needed * sizeof(int)); \
                if (!tmp) { free(needed_ids); return -1; } \
                needed_ids = tmp; \
            } \
            needed_ids[needed_count++] = (id); \
        } \
    } while(0)

    /* From cells: fill_transform */
    for (size_t i = 0; i < alea_vec_count(&dst->cells); i++) {
        const alea_cell_entry_t* cell = &dst->cells.data[i];
        if (cell->fill_universe > 0 && cell->fill_transform != 0) {
            ADD_NEEDED_ID(cell->fill_transform);
        }
    }

    /* From surfaces: transform_id */
    for (size_t i = 0; i < alea_vec_count(&dst->surfaces); i++) {
        int tr_id = dst->surfaces.data[i].transform_id;
        if (tr_id != 0) {
            ADD_NEEDED_ID(tr_id);
        }
    }
    #undef ADD_NEEDED_ID

    /* Copy each needed transform from src */
    for (size_t i = 0; i < needed_count; i++) {
        int tr_id = needed_ids[i];

        /* Check if already in dst */
        bool exists = false;
        for (size_t j = 0; j < alea_vec_count(&dst->transforms); j++) {
            if (dst->transforms.data[j].transform_id == tr_id) {
                exists = true;
                break;
            }
        }
        if (exists) continue;

        /* Find in src */
        for (size_t j = 0; j < alea_vec_count(&src->transforms); j++) {
            if (src->transforms.data[j].transform_id == tr_id) {
                alea_transform_t* tr = alea_vec_push_uninit(&dst->transforms, alea_transform_t);
                if (!tr) { free(needed_ids); return -1; }
                *tr = src->transforms.data[j];  /* POD copy */
                break;
            }
        }
    }

    free(needed_ids);
    return 0;
}

int alea_copy_referenced_cell_refs(alea_system_t* dst, const alea_system_t* src) {
    if (!dst || !src) return -1;
    if (alea_vec_empty(&src->cell_refs)) return 0;

    size_t dst_node_count = alea_vec_count(&dst->nodes);

    for (size_t i = 0; i < alea_vec_count(&src->cell_refs); i++) {
        const alea_cell_ref_t* ref = &src->cell_refs.data[i];

        /* Check if both referencing and referenced cells exist in dst */
        bool has_referencing = false;
        bool has_referenced = false;
        for (size_t c = 0; c < alea_vec_count(&dst->cells); c++) {
            int cid = dst->cells.data[c].mc_cell_id;
            if (cid == ref->referencing_cell_id) has_referencing = true;
            if (cid == ref->referenced_cell_id) has_referenced = true;
            if (has_referencing && has_referenced) break;
        }
        if (!has_referencing || !has_referenced) continue;

        /* Check that placeholder_node is valid in dst */
        if (ref->placeholder_node >= dst_node_count) continue;

        alea_cell_ref_t* dst_ref = alea_vec_push_uninit(&dst->cell_refs, alea_cell_ref_t);
        if (!dst_ref) return -1;
        *dst_ref = *ref;
    }
    return 0;
}

/* ============================================================================
 * FLATTEN CONFIG DEFAULTS
 * ============================================================================ */


const alea_flatten_config_t ALEA_FLATTEN_DEFAULT = {
    .starting_universe_id = 0,
    .max_depth = 0,
    .starting_cell_id = 1,
    .copy_materials = true,
    .copy_transforms = false,
    .clip_active = false,
    .clip_bbox = {0, 0, 0, 0, 0, 0}
};

/**
 * @brief Parent cell geometry node for tracking boundaries during flatten
 */
typedef struct parent_geom {
    alea_node_id_t node_id;       /* Geometry node in src system */
    alea_matrix_t transform;      /* Accumulated transform at this level */
    struct parent_geom* next;    /* Stack pointer */
} parent_geom_t;

/**
 * @brief Cache entry for parent tree cloning during flatten
 *
 * Caches the result of cloning a parent cell's geometry with a given transform
 * into the destination system. Bounded by (nesting depth * unique parent cells).
 */
typedef struct {
    alea_node_id_t src_node_id;   /* Source geometry node */
    double transform_m[12];      /* Copy of transform matrix */
    alea_node_id_t cloned_root;   /* Result in dst system */
} parent_clone_entry_t;

ALEA_VEC_DEFINE(parent_clone_vec, parent_clone_entry_t);

/**
 * @brief Context for recursive flattening into new system
 */
typedef struct {
    alea_system_t* dst;           /* Destination system */
    alea_system_t* src;           /* Source system */
    const alea_flatten_config_t* config;
    int next_cell_id;
    int error;
    parent_geom_t* parent_stack;
    primitive_remap_t* remap;
    /* Parent clone cache */
    parent_clone_vec_t parent_cache;
    /* Used cell ID tracking for best-effort preservation */
    uint32_t* used_cell_ids;     /* Bitset: bit N set => cell ID N is taken */
    size_t used_cell_ids_size;   /* Number of uint32_t words in bitset */
} flatten_context_t;

/* Bitset helpers for cell ID tracking */
static inline bool cell_id_is_used(const flatten_context_t* ctx, int id) {
    if (id <= 0) return true;
    size_t word = (size_t)id / 32;
    if (word >= ctx->used_cell_ids_size) return false;
    return (ctx->used_cell_ids[word] >> ((uint32_t)id % 32)) & 1;
}

static inline void cell_id_mark_used(flatten_context_t* ctx, int id) {
    if (id <= 0) return;
    size_t word = (size_t)id / 32;
    if (word >= ctx->used_cell_ids_size) {
        /* Grow bitset */
        size_t new_size = word + 1;
        if (new_size < 64) new_size = 64;  /* Min 2K words = ~64K IDs */
        uint32_t* new_bits = realloc(ctx->used_cell_ids, new_size * sizeof(uint32_t));
        if (!new_bits) return;
        memset(new_bits + ctx->used_cell_ids_size, 0,
               (new_size - ctx->used_cell_ids_size) * sizeof(uint32_t));
        ctx->used_cell_ids = new_bits;
        ctx->used_cell_ids_size = new_size;
    }
    ctx->used_cell_ids[word] |= (1u << ((uint32_t)id % 32));
}



/* ============================================================================
 * MATRIX OPERATIONS
 * ============================================================================ */

void alea_matrix_identity(alea_matrix_t* mat) {
    /* Row-major 3x4: identity rotation, zero translation */
    mat->m[0] = 1.0;  mat->m[1] = 0.0;  mat->m[2] = 0.0;  mat->m[3] = 0.0;
    mat->m[4] = 0.0;  mat->m[5] = 1.0;  mat->m[6] = 0.0;  mat->m[7] = 0.0;
    mat->m[8] = 0.0;  mat->m[9] = 0.0;  mat->m[10] = 1.0; mat->m[11] = 0.0;
    
    /* Inverse is also identity */
    memcpy(mat->inv, mat->m, sizeof(mat->m));
    mat->has_inverse = true;
}

bool alea_matrix_from_mcnp(alea_matrix_t* mat, const double* data,
                          int count, bool degrees) {
    if (!mat || !data || count < 1 || count > 12 ||
        (count > 3 && count != 12)) {
        return false;
    }

    alea_matrix_identity(mat);
    
    if (count >= 3) {
        /* Translation: ox, oy, oz */
        mat->m[3] = data[0];   /* Tx */
        mat->m[7] = data[1];   /* Ty */
        mat->m[11] = data[2];  /* Tz */
    }
    
    if (count >= 12) {
        /*
         * MCNP format: ox oy oz b1 b2 b3 b4 b5 b6 b7 b8 b9
         * where b1-b9 are either:
         *   - Direction cosines: b1=cos(main_x,aux_x), b2=cos(main_y,aux_x), ...
         *   - Angles in degrees (*TRn): angles between axes
         *
         * MCNP transform from auxiliary to main coordinates:
         *   x_main = ox + b1*x_aux + b4*y_aux + b7*z_aux
         *   y_main = oy + b2*x_aux + b5*y_aux + b8*z_aux
         *   z_main = oz + b3*x_aux + b6*y_aux + b9*z_aux
         *
         * So the rotation matrix R (aux->main) has COLUMNS [b1,b2,b3], [b4,b5,b6], [b7,b8,b9]:
         *   R = [ b1 b4 b7 ]
         *       [ b2 b5 b8 ]
         *       [ b3 b6 b9 ]
         *   p_main = R * p_aux + T
         */
        if (degrees) {
            /* Convert angles to direction cosines */
            /* MCNP *TRn: angles between primed and unprimed axes */
            double b[9];
            for (int i = 0; i < 9; i++) {
                double angle_deg = data[3 + i];
                b[i] = cos(DEG_TO_RAD(angle_deg));
            }
            if (g_debug_point_trace) {
                ALEA_LOG_DEBUG("alea_matrix_from_mcnp (degrees mode):");
                ALEA_LOG_DEBUG("  Raw angles: %.2f %.2f %.2f  %.2f %.2f %.2f  %.2f %.2f %.2f",
                       data[3], data[4], data[5], data[6], data[7], data[8], data[9], data[10], data[11]);
                ALEA_LOG_DEBUG("  Cosines: %.6f %.6f %.6f  %.6f %.6f %.6f  %.6f %.6f %.6f",
                       b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8]);
            }
            /* MCNP convention: p_main = T + [b1 b4 b7; b2 b5 b8; b3 b6 b9] * p_aux
             * So the rotation matrix is TRANSPOSED from the raw b values */
            mat->m[0] = b[0];  mat->m[1] = b[3];  mat->m[2] = b[6];
            mat->m[4] = b[1];  mat->m[5] = b[4];  mat->m[6] = b[7];
            mat->m[8] = b[2];  mat->m[9] = b[5];  mat->m[10] = b[8];
        } else {
            /* Direct direction cosines - also need to transpose */
            /* MCNP: b1-b3 describe aux x-axis in main coords (column, not row) */
            if (g_debug_point_trace) {
                ALEA_LOG_DEBUG("alea_matrix_from_mcnp (cosines mode):");
                ALEA_LOG_DEBUG("  Raw cosines: %.6f %.6f %.6f  %.6f %.6f %.6f  %.6f %.6f %.6f",
                       data[3], data[4], data[5], data[6], data[7], data[8], data[9], data[10], data[11]);
            }
            mat->m[0] = data[3];  mat->m[1] = data[6];  mat->m[2] = data[9];
            mat->m[4] = data[4];  mat->m[5] = data[7];  mat->m[6] = data[10];
            mat->m[8] = data[5];  mat->m[9] = data[8];  mat->m[10] = data[11];
        }
    }
    
    mat->has_inverse = false;
    return alea_matrix_invert(mat);
}

void alea_matrix_multiply(alea_matrix_t* result, 
                         const alea_matrix_t* a, 
                         const alea_matrix_t* b) {
    /* result = a * b (b applied first, then a) */
    /* 3x4 * 3x4 treating as 4x4 with implicit [0,0,0,1] row */
    
    double r[12];
    
    /* Row 0 */
    r[0] = a->m[0]*b->m[0] + a->m[1]*b->m[4] + a->m[2]*b->m[8];
    r[1] = a->m[0]*b->m[1] + a->m[1]*b->m[5] + a->m[2]*b->m[9];
    r[2] = a->m[0]*b->m[2] + a->m[1]*b->m[6] + a->m[2]*b->m[10];
    r[3] = a->m[0]*b->m[3] + a->m[1]*b->m[7] + a->m[2]*b->m[11] + a->m[3];
    
    /* Row 1 */
    r[4] = a->m[4]*b->m[0] + a->m[5]*b->m[4] + a->m[6]*b->m[8];
    r[5] = a->m[4]*b->m[1] + a->m[5]*b->m[5] + a->m[6]*b->m[9];
    r[6] = a->m[4]*b->m[2] + a->m[5]*b->m[6] + a->m[6]*b->m[10];
    r[7] = a->m[4]*b->m[3] + a->m[5]*b->m[7] + a->m[6]*b->m[11] + a->m[7];
    
    /* Row 2 */
    r[8] = a->m[8]*b->m[0] + a->m[9]*b->m[4] + a->m[10]*b->m[8];
    r[9] = a->m[8]*b->m[1] + a->m[9]*b->m[5] + a->m[10]*b->m[9];
    r[10] = a->m[8]*b->m[2] + a->m[9]*b->m[6] + a->m[10]*b->m[10];
    r[11] = a->m[8]*b->m[3] + a->m[9]*b->m[7] + a->m[10]*b->m[11] + a->m[11];
    
    memcpy(result->m, r, sizeof(r));
    result->has_inverse = false;
}

bool alea_matrix_invert(alea_matrix_t* mat) {
    if (mat->has_inverse) return true;
    
    /* For rigid transforms (rotation + translation): R^-1 = R^T, T^-1 = -R^T * T */
    /* Extract rotation part (3x3) */
    double R[9] = {
        mat->m[0], mat->m[1], mat->m[2],
        mat->m[4], mat->m[5], mat->m[6],
        mat->m[8], mat->m[9], mat->m[10]
    };
    
    /* Check if it's orthogonal (rigid transform) */
    double det = R[0]*(R[4]*R[8] - R[5]*R[7]) 
               - R[1]*(R[3]*R[8] - R[5]*R[6]) 
               + R[2]*(R[3]*R[7] - R[4]*R[6]);
    
    if (fabs(det) < 1e-10) {
        /* Singular matrix */
        return false;
    }
    
    /* For orthogonal matrices, inverse = transpose */
    /* But let's compute general inverse for robustness */
    double inv_det = 1.0 / det;
    
    double Rinv[9];
    Rinv[0] = (R[4]*R[8] - R[5]*R[7]) * inv_det;
    Rinv[1] = (R[2]*R[7] - R[1]*R[8]) * inv_det;
    Rinv[2] = (R[1]*R[5] - R[2]*R[4]) * inv_det;
    Rinv[3] = (R[5]*R[6] - R[3]*R[8]) * inv_det;
    Rinv[4] = (R[0]*R[8] - R[2]*R[6]) * inv_det;
    Rinv[5] = (R[2]*R[3] - R[0]*R[5]) * inv_det;
    Rinv[6] = (R[3]*R[7] - R[4]*R[6]) * inv_det;
    Rinv[7] = (R[1]*R[6] - R[0]*R[7]) * inv_det;
    Rinv[8] = (R[0]*R[4] - R[1]*R[3]) * inv_det;
    
    /* Translation */
    double T[3] = { mat->m[3], mat->m[7], mat->m[11] };
    
    /* Inverse translation: -R^-1 * T */
    double Tinv[3];
    Tinv[0] = -(Rinv[0]*T[0] + Rinv[1]*T[1] + Rinv[2]*T[2]);
    Tinv[1] = -(Rinv[3]*T[0] + Rinv[4]*T[1] + Rinv[5]*T[2]);
    Tinv[2] = -(Rinv[6]*T[0] + Rinv[7]*T[1] + Rinv[8]*T[2]);
    
    /* Store inverse */
    mat->inv[0] = Rinv[0]; mat->inv[1] = Rinv[1]; mat->inv[2] = Rinv[2]; mat->inv[3] = Tinv[0];
    mat->inv[4] = Rinv[3]; mat->inv[5] = Rinv[4]; mat->inv[6] = Rinv[5]; mat->inv[7] = Tinv[1];
    mat->inv[8] = Rinv[6]; mat->inv[9] = Rinv[7]; mat->inv[10] = Rinv[8]; mat->inv[11] = Tinv[2];
    
    mat->has_inverse = true;
    return true;
}

void alea_matrix_transform_point(const alea_matrix_t* mat,
                                double* x, double* y, double* z) {
    double px = *x, py = *y, pz = *z;
    *x = mat->m[0]*px + mat->m[1]*py + mat->m[2]*pz + mat->m[3];
    *y = mat->m[4]*px + mat->m[5]*py + mat->m[6]*pz + mat->m[7];
    *z = mat->m[8]*px + mat->m[9]*py + mat->m[10]*pz + mat->m[11];
}

void alea_matrix_transform_point_inverse(const alea_matrix_t* mat,
                                        double* x, double* y, double* z) {
    if (!mat->has_inverse) {
        ALEA_LOG_WARN("transform_point_inverse called without inverse matrix");
        return;
    }
    double px = *x, py = *y, pz = *z;
    *x = mat->inv[0]*px + mat->inv[1]*py + mat->inv[2]*pz + mat->inv[3];
    *y = mat->inv[4]*px + mat->inv[5]*py + mat->inv[6]*pz + mat->inv[7];
    *z = mat->inv[8]*px + mat->inv[9]*py + mat->inv[10]*pz + mat->inv[11];
}


/* ============================================================================
 * UNIFIED TREE CLONING
 * ============================================================================ */

#define MAX_CLONE_DEPTH 10000

static alea_node_id_t clone_tree_impl(alea_system_t* dst, const alea_system_t* src,
                                      alea_node_id_t root, const alea_matrix_t* mat,
                                      primitive_remap_t* remap, int depth) {
    if (!dst || !src) {
        ALEA_LOG_ERROR("clone_tree_ex: NULL dst or src");
        return ALEA_NODE_ID_INVALID;
    }
    if (root == ALEA_NODE_ID_INVALID || root >= alea_vec_count(&src->nodes)) {
        ALEA_LOG_ERROR("clone_tree_ex: invalid root %u (node_count=%zu)",
                root, alea_vec_count(&src->nodes));
        return ALEA_NODE_ID_INVALID;
    }
    if (depth > MAX_CLONE_DEPTH) {
        ALEA_LOG_ERROR("clone_tree_ex: recursion depth %d exceeds limit (node %u)",
                depth, root);
        return ALEA_NODE_ID_INVALID;
    }

    const alea_node_t* node = &src->nodes.data[root];
    alea_operation_t op = ALEA_GET_OPERATION(node);

    if (op > ALEA_OP_COMPLEMENT) {
        ALEA_LOG_ERROR("clone_tree_ex: invalid op %d at node %u", op, root);
        return ALEA_NODE_ID_INVALID;
    }

    if (op == ALEA_OP_PRIMITIVE) {
        uint32_t src_prim_id = node->primitive.primitive_id;
        if (src_prim_id >= alea_vec_count(&src->primitives)) {
            ALEA_LOG_ERROR("clone_tree_ex: prim_id %u >= primitive_count %zu",
                    src_prim_id, alea_vec_count(&src->primitives));
            return ALEA_NODE_ID_INVALID;
        }

        uint32_t dst_prim_id;
        int8_t transform_inverted = 0;

        if (mat) {
            /* Transform primitive — always re-transform (remap cache is not
             * valid across calls with different transforms) */
            alea_primitive_entry_t* src_prim = &src->primitives.data[src_prim_id];
            alea_primitive_data_t transformed_data;
            alea_primitive_type_t transformed_type;
            bool xform_ok = alea_primitive_transform(src_prim->type, &src_prim->data, mat->m,
                                                     &transformed_type, &transformed_data);
            if (!xform_ok) {
                /* Per-type transform refused to bake (e.g., a rotated torus
                 * whose axis no longer aligns with X/Y/Z — MCNP only
                 * supports axis-aligned tori, and MCNP does NOT accept a
                 * surface-level TR prefix to rescue this case when the TR
                 * itself doesn't preserve axis alignment).
                 *
                 * Producing un-transformed output here would silently
                 * corrupt geometry, so fail hard. The principled fix is
                 * partial-flatten: stop flattening at the fill level whose
                 * accumulated rotation would skew a torus, and emit a fill
                 * cell instead — see docs/PLAN_PARTIAL_FLATTEN_TORUS.md. */
                ALEA_LOG_ERROR("Cannot bake transform into primitive (type %s, src_id=%u, node=%u). "
                               "For a torus this means the accumulated rotation breaks axis alignment, "
                               "which MCNP cannot express. Aborting flatten/extract.",
                               alea_primitive_type_name(src_prim->type), src_prim_id, root);
                return ALEA_NODE_ID_INVALID;
            }

            dst_prim_id = alea_get_or_create_primitive(
                dst, transformed_type, &transformed_data, &transform_inverted);
            if (dst_prim_id == UINT32_MAX) {
                ALEA_LOG_ERROR("clone_tree_ex: transform primitive failed (node %u)", root);
                return ALEA_NODE_ID_INVALID;
            }
        } else if (remap && src_prim_id < remap->capacity &&
                   remap->map[src_prim_id] != ALEA_PRIMITIVE_ID_INVALID) {
            /* No transform — check remap cache (dedup guarantees same result) */
            dst_prim_id = remap->map[src_prim_id];
            transform_inverted = remap->inverted[src_prim_id];
        } else {
            /* No transform, cache miss — copy data directly */
            alea_primitive_entry_t* src_prim = &src->primitives.data[src_prim_id];
            alea_primitive_data_t data_copy = src_prim->data;

            dst_prim_id = alea_get_or_create_primitive(
                dst, src_prim->type, &data_copy, &transform_inverted);
            if (dst_prim_id == UINT32_MAX) {
                ALEA_LOG_ERROR("clone_tree_ex: copy primitive failed (node %u)", root);
                return ALEA_NODE_ID_INVALID;
            }

            if (remap && src_prim_id < remap->capacity) {
                remap->map[src_prim_id] = dst_prim_id;
                remap->inverted[src_prim_id] = transform_inverted;
            }
        }

        /* Preserve mc_surface_id when no transform, clear when transformed */
        int mc_surf_id = mat ? 0 : node->primitive.mc_surface_id;

        int8_t xored_inverted = node->primitive.inverted ^ transform_inverted;
        alea_node_id_t dst_node = alea_add_primitive_node(
            dst, dst_prim_id, node->primitive.sense, xored_inverted, mc_surf_id);
        if (dst_node == ALEA_NODE_ID_INVALID) {
            ALEA_LOG_ERROR("clone_tree_ex: add_primitive_node failed (node %u)", root);
            return ALEA_NODE_ID_INVALID;
        }

        dst->nodes.data[dst_node].material_id = node->material_id;
        return dst_node;

    } else {
        /* Operation node — recurse */
        alea_node_id_t new_left = clone_tree_impl(
            dst, src, node->operation.left, mat, remap, depth + 1);
        if (new_left == ALEA_NODE_ID_INVALID) {
            ALEA_LOG_ERROR("clone_tree_ex: left child failed (node %u op=%d)", root, op);
            return ALEA_NODE_ID_INVALID;
        }

        alea_node_id_t new_right = ALEA_NODE_ID_INVALID;
        if (op != ALEA_OP_COMPLEMENT) {
            new_right = clone_tree_impl(
                dst, src, node->operation.right, mat, remap, depth + 1);
            if (new_right == ALEA_NODE_ID_INVALID) {
                ALEA_LOG_ERROR("clone_tree_ex: right child failed (node %u op=%d)", root, op);
                return ALEA_NODE_ID_INVALID;
            }
        }

        alea_node_id_t new_node;
        switch (op) {
            case ALEA_OP_UNION:
                new_node = alea_create_union(dst, new_left, new_right);
                break;
            case ALEA_OP_INTERSECTION:
                new_node = alea_create_intersection(dst, new_left, new_right);
                break;
            case ALEA_OP_DIFFERENCE:
                new_node = alea_create_difference(dst, new_left, new_right);
                break;
            case ALEA_OP_COMPLEMENT:
                new_node = alea_create_complement(dst, new_left);
                break;
            default:
                return ALEA_NODE_ID_INVALID;
        }
        return new_node;
    }
}

/* Wrapper: same-system clone with transform (raw arrays) */
alea_node_id_t alea_clone_tree_transformed(alea_system_t* sys,
                                         alea_node_id_t root_id,
                                         const alea_matrix_t* mat,
                                         uint32_t* prim_remap,
                                         int8_t* prim_inverted,
                                         size_t remap_size) {
    /* Wrap raw arrays into a stack-local remap struct */
    primitive_remap_t local_remap;
    primitive_remap_t* remap_ptr = NULL;
    if (prim_remap) {
        local_remap.map = prim_remap;
        local_remap.inverted = prim_inverted;
        local_remap.capacity = remap_size;
        remap_ptr = &local_remap;
    }
    return clone_tree_impl(sys, sys, root_id, mat, remap_ptr, 0);
}

/* ============================================================================
 * LATTICE INDEX COMPUTATION
 * ============================================================================ */

/**
 * Look up which universe a point maps to in a rectangular lattice.
 * Returns the universe ID, or -1 if the point is outside the lattice bounds.
 * Writes the element origin into ox, oy, oz for coordinate translation.
 */
static int lattice_rect_lookup(const alea_cell_entry_t* cell,
                               double px, double py, double pz,
                               double* ox, double* oy, double* oz) {
    int ni = cell->lat_fill_dims[1] - cell->lat_fill_dims[0] + 1;
    int nj = cell->lat_fill_dims[3] - cell->lat_fill_dims[2] + 1;
    int nk = cell->lat_fill_dims[5] - cell->lat_fill_dims[4] + 1;

    int i = (int)floor((px - cell->lat_lower_left[0]) / cell->lat_pitch[0]);
    int j = (nj == 1) ? 0 : (int)floor((py - cell->lat_lower_left[1]) / cell->lat_pitch[1]);
    int k = (nk == 1) ? 0 : (int)floor((pz - cell->lat_lower_left[2]) / cell->lat_pitch[2]);

    /* Clamp to bounds */
    if (i < 0 || i >= ni || j < 0 || j >= nj || k < 0 || k >= nk)
        return -1;

    /* Element origin = lower_left + (i+0.5, j+0.5, k+0.5) * pitch */
    *ox = cell->lat_lower_left[0] + (i + 0.5) * cell->lat_pitch[0];
    *oy = cell->lat_lower_left[1] + (j + 0.5) * cell->lat_pitch[1];
    *oz = cell->lat_lower_left[2] + (k + 0.5) * cell->lat_pitch[2];

    size_t idx = (size_t)(i * nj * nk + j * nk + k);
    if (idx >= cell->lat_fill_count) return -1;

    return cell->lat_fill[idx];
}

/**
 * Look up which universe a point maps to in a hexagonal lattice.
 * Flat-top hex orientation, pitch = center-to-center distance.
 * Basis: e1 = (p, 0), e2 = (p/2, p*sqrt(3)/2).
 */
int lattice_hex_lookup(const alea_cell_entry_t* cell,
                              double px, double py, double pz,
                              double* ox, double* oy, double* oz) {
    double p = cell->lat_pitch[0];
    if (p <= 0.0) return -1;

    /* Fractional axial coordinates from inverse of basis matrix.
     * Note: hex lattices use 'center' (not lower_left) as the origin in
     * OpenMC; lat_lower_left stores the bounding box corner and should
     * NOT be subtracted here (unlike rectangular lattices). If a hex
     * center offset is needed, add a lat_center field. */
    double fj = py / (p * M_SQRT3 * 0.5);
    double fi = px / p - fj * 0.5;

    /* Convert to cube coordinates (x + y + z = 0) */
    double cx = fi, cz = fj, cy = -fi - fj;

    /* Round to nearest hex */
    int ri = (int)round(cx);
    int rj = (int)round(cy);
    int rk = (int)round(cz);

    double dx = fabs(ri - cx);
    double dy = fabs(rj - cy);
    double dz = fabs(rk - cz);

    if (dx > dy && dx > dz)
        ri = -rj - rk;
    else if (dy > dz)
        rj = -ri - rk;
    else
        rk = -ri - rj;

    /* Axial result: (i, j) = (ri, rk) */
    int ni = cell->lat_fill_dims[1] - cell->lat_fill_dims[0] + 1;
    int nj = cell->lat_fill_dims[3] - cell->lat_fill_dims[2] + 1;
    int nk = cell->lat_fill_dims[5] - cell->lat_fill_dims[4] + 1;

    int oi = ri - cell->lat_fill_dims[0];
    int oj = rk - cell->lat_fill_dims[2];

    /* Z index (linear, same as rectangular) */
    int ok = (nk == 1) ? 0
           : (int)floor((pz - cell->lat_lower_left[2]) / cell->lat_pitch[2]);

    if (oi < 0 || oi >= ni || oj < 0 || oj >= nj || ok < 0 || ok >= nk)
        return -1;

    /* Element center in Cartesian */
    *ox = ri * p + rk * p * 0.5;
    *oy = rk * p * M_SQRT3 * 0.5;
    *oz = (nk == 1) ? 0.0 : cell->lat_lower_left[2] + (ok + 0.5) * cell->lat_pitch[2];

    size_t idx = (size_t)(oi * nj * nk + oj * nk + ok);
    if (idx >= cell->lat_fill_count) return -1;

    return cell->lat_fill[idx];
}

/* ============================================================================
 * LAZY POINT QUERY
 * ============================================================================ */

static int find_cell_recursive(const alea_system_t* sys,
                            double x, double y, double z,
                            int universe_id,
                            const alea_matrix_t* accumulated_transform,
                            int depth,
                            int max_depth,
                            int* out_cell_index,
                            int* out_cell_id,
                            int* out_material,
                            int* out_depth) {
    if (max_depth > 0 && depth >= max_depth) {
        return -1;  /* Max depth exceeded */
    }
    
    /* Transform point to this universe's coordinate system */
    double lx = x, ly = y, lz = z;
    if (accumulated_transform) {
        if (!accumulated_transform->has_inverse) return -1;
        alea_matrix_transform_point_inverse(accumulated_transform, &lx, &ly, &lz);
    }
    
    /* Find universe */
    const alea_universe_t* univ = alea_get_universe(sys, universe_id);
    if (!univ) return -1;
    
    /* Test each cell in universe */
    for (size_t i = 0; i < univ->cell_indices.count; i++) {
        size_t cell_idx = univ->cell_indices.data[i];
        const alea_cell_entry_t* cell = &sys->cells.data[cell_idx];
        
        /* Lattice cell: bounds check replaces CSG containment test */
        if (cell->lat_type != 0 && cell->lat_fill) {
            double ox, oy, oz;
            int fill_univ = (cell->lat_type == 2)
                ? lattice_hex_lookup(cell, lx, ly, lz, &ox, &oy, &oz)
                : lattice_rect_lookup(cell, lx, ly, lz, &ox, &oy, &oz);
            if (fill_univ < 0) continue;

            double elx = lx - ox, ely = ly - oy, elz = lz - oz;

            int result = find_cell_recursive(
                sys, elx, ely, elz, fill_univ,
                NULL, depth + 1, max_depth,
                out_cell_index, out_cell_id, out_material, out_depth);
            if (result == 0) return 0;
            continue;
        }

        /* Non-lattice: evaluate CSG geometry */
        bool inside = alea_contains_point(sys, cell->root_node_id, lx, ly, lz);

        if (inside) {
            /* Regular FILL: recurse into fill universe */
            if (cell->fill_universe > 0) {
                alea_matrix_t fill_transform;

                if (cell->fill_transform > 0) {
                    const alea_transform_t* tr = alea_get_transform(sys, cell->fill_transform);
                    if (tr) {
                        if (!alea_matrix_from_mcnp(&fill_transform, tr->cosines,
                                                   tr->value_count, false)) {
                            return -1;
                        }
                    } else {
                        return -1;
                    }
                } else {
                    alea_matrix_identity(&fill_transform);
                }

                alea_matrix_t new_accumulated;
                if (accumulated_transform) {
                    alea_matrix_multiply(&new_accumulated, accumulated_transform, &fill_transform);
                } else {
                    new_accumulated = fill_transform;
                }
                if (!alea_matrix_invert(&new_accumulated)) return -1;

                int result = find_cell_recursive(
                    sys, x, y, z, cell->fill_universe,
                    &new_accumulated, depth + 1, max_depth,
                    out_cell_index, out_cell_id, out_material, out_depth);

                if (result == 0) return 0;

            /* Terminal cell */
            } else {
                if (out_cell_index) *out_cell_index = (int)cell_idx;
                if (out_cell_id) *out_cell_id = cell->mc_cell_id;
                if (out_material) *out_material = cell->material_id;
                if (out_depth) *out_depth = depth;
                return 0;
            }
        }
    }
    
    return -1;  /* Not found in any cell (void) */
}

int alea_find_cell_lazy(const alea_system_t* sys,
                       double x, double y, double z,
                       int* out_cell_id,
                       int* out_material,
                       int* out_depth) {
    if (!sys) return -1;
    if (!sys->universe_index_built) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              "universe index is not prepared; call alea_prepare_query_acceleration()");
        return -1;
    }
    
    int cell_id = -1, material = 0;
    int depth = -1;
    int result = find_cell_recursive(
        sys, x, y, z, 0, NULL, 0, MAX_FLATTEN_DEPTH,
        NULL, &cell_id, &material, &depth);
    
    if (out_cell_id) *out_cell_id = cell_id;
    if (out_material) *out_material = material;
    if (out_depth) *out_depth = depth;

    return result;
}

void alea_set_debug_point_trace(int enable) {
    g_debug_point_trace = enable;
}

/* Internal recursive helper for alea_find_all_cells_at_point */
static int find_all_cells_recursive(const alea_system_t* sys,
                                    double gx, double gy, double gz,
                                    double lx, double ly, double lz,
                                    int universe_id,
                                    const alea_matrix_t* accumulated_transform,
                                    int depth,
                                    alea_cell_hit_t* out_hits,
                                    size_t max_hits,
                                    size_t* hit_count) {
    if (*hit_count >= max_hits) return 0;
    if (depth >= MAX_FLATTEN_DEPTH) return 0;

    if (g_debug_point_trace) {
        ALEA_LOG_DEBUG("depth=%d universe=%d global=(%.4f,%.4f,%.4f) local=(%.4f,%.4f,%.4f)",
               depth, universe_id, gx, gy, gz, lx, ly, lz);
    }

    const alea_universe_t* univ = alea_get_universe(sys, universe_id);
    if (!univ) return -1;

    alea_universe_t* mutable_univ = (alea_universe_t*)univ;
    if (ensure_universe_point_bvh((alea_system_t*)sys, mutable_univ) == 0 &&
        mutable_univ->point_bvh_built) {
        alea_size_vec_t candidates = ALEA_VEC_INIT;
        universe_point_query_ctx_t ctx = {
            .sys = sys,
            .univ = univ,
            .gx = gx, .gy = gy, .gz = gz,
            .lx = lx, .ly = ly, .lz = lz,
            .accumulated_transform = accumulated_transform,
            .depth = depth,
            .candidate_positions = &candidates,
            .error = 0
        };
        g_point_bvh_stats.bvh_queries++;
        traverse_universe_point_bvh_node(&ctx, 0);
        if (!ctx.error && candidates.count > 1) {
            qsort(candidates.data, candidates.count, sizeof(size_t),
                  compare_size_t_ascending);
        }
        for (size_t i = 0; !ctx.error && i < candidates.count; i++) {
            size_t cell_idx = univ->cell_indices.data[candidates.data[i]];
            int rc = process_cell_for_all_cells_query(sys, cell_idx,
                                                      gx, gy, gz,
                                                      lx, ly, lz,
                                                      universe_id,
                                                      accumulated_transform,
                                                      depth,
                                                      out_hits,
                                                      max_hits,
                                                      hit_count);
            if (rc < 0) ctx.error = -1;
        }
        alea_vec_free(&candidates);
        return ctx.error ? -1 : 0;
    }

    /* Small/lattice-bearing universes keep the original linear scan. */
    g_point_bvh_stats.linear_universe_scans++;
    g_point_bvh_stats.linear_cell_tests += univ->cell_indices.count;
    for (size_t i = 0; i < univ->cell_indices.count; i++) {
        size_t cell_idx = univ->cell_indices.data[i];
        int rc = process_cell_for_all_cells_query(sys, cell_idx,
                                                  gx, gy, gz,
                                                  lx, ly, lz,
                                                  universe_id,
                                                  accumulated_transform,
                                                  depth,
                                                  out_hits,
                                                  max_hits,
                                                  hit_count);
        if (rc < 0) return -1;
    }

    return 0;
}

static int find_all_cells_at_point_impl(alea_system_t* sys,
                                        double x, double y, double z,
                                        alea_cell_hit_t* out_hits,
                                        size_t max_hits,
                                        bool force_recursive) {
    if (!sys || !out_hits || max_hits == 0) return -1;
    if (!sys->universe_index_built) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              "universe index is not prepared; call alea_prepare_query_acceleration()");
        return -1;
    }

    /* Debug-trace path stays on the recursive walker so it can emit per-step
     * logs. force_recursive (caller-requested cache bypass for overlap
     * detection) is satisfied equally well by the hier path — the hier index
     * has no per-point cache to invalidate, so it always returns full
     * multi-depth hits. Important: the recursive walker lazily builds a
     * per-universe point BVH on first use; that lazy build is NOT thread-safe
     * and segfaults under OpenMP. Routing through hier when available also
     * dodges that race entirely. */
    if (g_debug_point_trace) {
        size_t hit_count = 0;
        int result = find_all_cells_recursive(sys, x, y, z, x, y, z,
                                              0, NULL, 0,
                                              out_hits, max_hits, &hit_count);
        if (result < 0) return -1;
        return (int)hit_count;
    }

    /* Hier-spatial fast path: if the hierarchical index is built, use it
     * instead of the flat spatial cache or the O(N)-per-level recursive
     * fallback. Critical for hier-mode plotter/raycast on large models
     * where the recursive scan would dominate per-pixel cost. */
    if (sys->hier_spatial_index) {
        int n = alea_hier_spatial_find_cells_at_point(sys, x, y, z,
                                                      out_hits, max_hits);
        if (n >= 0) return n;
        /* Fall through to flat/recursive on error. */
    }

    /* No hier index — fall back to recursive (single-threaded plotter path,
     * or when hier index isn't built). The per-universe point BVH inside
     * find_all_cells_recursive remains lazily built; that's fine here because
     * the recursive path is only reached when OpenMP isn't doing parallel
     * grid queries through the hier branch above. */
    if (force_recursive) {
        size_t hit_count = 0;
        int result = find_all_cells_recursive(sys, x, y, z, x, y, z,
                                              0, NULL, 0,
                                              out_hits, max_hits, &hit_count);
        if (result < 0) return -1;
        return (int)hit_count;
    }

    /* Try fast spatial lookup with coherence caching.
     * Will auto-build the spatial index if needed. */
    int result = alea_spatial_find_cells_at_point(sys, x, y, z, out_hits, max_hits);
    if (result >= 0) {
        if (!sys->has_lattice) return result;

        /* Spatial path doesn't handle lattices.  Check if any hit landed on
         * a lattice cell, or if result==0 the point might be in a lattice
         * element outside the base cell geometry. */
        bool need_recursive = (result == 0);
        for (int i = 0; i < result && !need_recursive; i++) {
            const alea_cell_entry_t* c = &sys->cells.data[out_hits[i].cell_index];
            if (c->lat_type != 0 && c->lat_fill)
                need_recursive = true;
        }
        /* Also recurse if the spatial path couldn't fully resolve a fill
         * chain (deepest hit still has fill_universe > 0).  This happens
         * when a fill points to a lattice universe whose synthetic cell
         * has no CSG tree and thus no BVH entry. */
        if (!need_recursive && result > 0) {
            const alea_cell_entry_t* deepest =
                &sys->cells.data[out_hits[result - 1].cell_index];
            if (deepest->fill_universe > 0)
                need_recursive = true;
        }
        if (!need_recursive) return result;
    }

    /* Fallback: recursive approach (O(N) per level) */
    size_t hit_count = 0;

    result = find_all_cells_recursive(sys, x, y, z, x, y, z,
                                      0, NULL, 0,
                                      out_hits, max_hits, &hit_count);

    if (result < 0) return -1;
    return (int)hit_count;
}

int alea_find_all_cells_at_point(alea_system_t* sys,
                                double x, double y, double z,
                                alea_cell_hit_t* out_hits,
                                size_t max_hits) {
    return find_all_cells_at_point_impl(sys, x, y, z, out_hits, max_hits, false);
}

int alea_find_deepest_cell_hit_at_point(alea_system_t* sys,
                                       double x, double y, double z,
                                       alea_cell_hit_t* out_hit) {
    if (!sys || !out_hit) return -1;

    bool has_hierarchy = false;
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        const alea_cell_entry_t* cell = &sys->cells.data[i];
        if (cell->fill_universe > 0 || (cell->lat_type != 0 && cell->lat_fill)) {
            has_hierarchy = true;
            break;
        }
    }

    if (!has_hierarchy) {
        const alea_universe_t* base = alea_get_universe(sys, 0);
        if (!base) return -1;

        for (size_t i = 0; i < base->cell_indices.count; i++) {
            size_t cell_idx = base->cell_indices.data[i];
            const alea_cell_entry_t* cell = &sys->cells.data[cell_idx];
            if (!alea_contains_point(sys, cell->root_node_id, x, y, z)) continue;

            out_hit->cell_id = cell->mc_cell_id;
            out_hit->cell_index = (int)cell_idx;
            out_hit->material_id = cell->material_id;
            out_hit->universe_id = cell->universe_id;
            out_hit->fill_universe = cell->fill_universe;
            out_hit->depth = 0;
            out_hit->local_x = x;
            out_hit->local_y = y;
            out_hit->local_z = z;
            return 0;
        }

        return -1;
    }

    if (!sys->universe_index_built) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              "universe index is not prepared; call alea_prepare_query_acceleration()");
        return -1;
    }

    int cell_index = -1;
    int depth = -1;
    if (find_cell_recursive(sys, x, y, z, 0, NULL, 0, MAX_FLATTEN_DEPTH,
                            &cell_index, &out_hit->cell_id,
                            &out_hit->material_id, &depth) != 0) {
        return -1;
    }

    if (cell_index < 0 || (size_t)cell_index >= alea_vec_count(&sys->cells)) {
        return -1;
    }

    const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
    out_hit->cell_index = cell_index;
    out_hit->universe_id = cell->universe_id;
    out_hit->fill_universe = cell->fill_universe;
    out_hit->depth = depth;
    out_hit->local_x = x;
    out_hit->local_y = y;
    out_hit->local_z = z;
    return 0;
}

/* Debug version that always uses recursive path */
int alea_find_all_cells_at_point_recursive(const alea_system_t* sys,
                                          double x, double y, double z,
                                          alea_cell_hit_t* out_hits,
                                          size_t max_hits) {
    return find_all_cells_at_point_impl((alea_system_t*)sys, x, y, z, out_hits, max_hits, true);
}

/* ============================================================================
 * EXPLICIT FLATTENING
 * ============================================================================ */

/* Forward declarations for static helpers */
static alea_node_id_t clone_tree_to_system_transformed(alea_system_t* dst,
                                                       const alea_system_t* src,
                                                       alea_node_id_t src_root,
                                                       const alea_matrix_t* mat,
                                                       primitive_remap_t* remap);

/**
 * @brief Container bbox of a lattice cell in its own (local) frame.
 *
 * For rectangular lattices this is just lat_lower_left + dims*pitch. For
 * hex lattices it's a conservative axis-aligned cover of all element
 * centers expanded by one pitch (mirrors lattice_container_bbox in
 * alea_spatial_hier.c). Used so the region cull works on the true lattice
 * extent rather than on the stored CSG bbox (which may be a single
 * element's region for MCNP LAT cells).
 */
static alea_bbox_t flatten_lattice_container_bbox(const alea_cell_entry_t* cell) {
    if (cell->lat_type == 1) {
        int ni = cell->lat_fill_dims[1] - cell->lat_fill_dims[0] + 1;
        int nj = cell->lat_fill_dims[3] - cell->lat_fill_dims[2] + 1;
        int nk = cell->lat_fill_dims[5] - cell->lat_fill_dims[4] + 1;
        return (alea_bbox_t){
            cell->lat_lower_left[0],
            cell->lat_lower_left[0] + (double)ni * cell->lat_pitch[0],
            cell->lat_lower_left[1],
            cell->lat_lower_left[1] + (double)nj * cell->lat_pitch[1],
            cell->lat_lower_left[2],
            cell->lat_lower_left[2] + (double)nk * cell->lat_pitch[2]
        };
    }
    /* Hex lattice. */
    double p = cell->lat_pitch[0] > 0.0 ? cell->lat_pitch[0] : 1.0;
    int d0 = cell->lat_fill_dims[0];
    int d1 = cell->lat_fill_dims[1];
    int d2 = cell->lat_fill_dims[2];
    int d3 = cell->lat_fill_dims[3];
    int nk = cell->lat_fill_dims[5] - cell->lat_fill_dims[4] + 1;

    /* Element centers at (ri*p + rk*p/2, rk*p*sqrt(3)/2). Cover all index
     * extremes and pad by one pitch. */
    double cx[4] = {
        d0 * p + d2 * p * 0.5,
        d0 * p + d3 * p * 0.5,
        d1 * p + d2 * p * 0.5,
        d1 * p + d3 * p * 0.5
    };
    double cy[2] = { d2 * p * M_SQRT3 * 0.5, d3 * p * M_SQRT3 * 0.5 };

    double xmin = cx[0], xmax = cx[0];
    for (int i = 1; i < 4; i++) {
        if (cx[i] < xmin) xmin = cx[i];
        if (cx[i] > xmax) xmax = cx[i];
    }
    double ymin = cy[0] < cy[1] ? cy[0] : cy[1];
    double ymax = cy[0] < cy[1] ? cy[1] : cy[0];

    double zmin = (nk == 1) ? -p : cell->lat_lower_left[2];
    double zmax = (nk == 1) ?  p
        : cell->lat_lower_left[2] + (double)nk * cell->lat_pitch[2];

    return (alea_bbox_t){xmin - p, xmax + p, ymin - p, ymax + p, zmin, zmax};
}

/**
 * @brief Project a world-frame bbox conservatively into the lattice cell's
 * local frame by inverting accumulated_transform and AABB-ing the 8 corners.
 * If accumulated_transform is NULL the world frame == local frame.
 */
static alea_bbox_t flatten_project_clip_to_local(const alea_bbox_t* world_clip,
                                                 const alea_matrix_t* accumulated) {
    if (!accumulated) return *world_clip;

    alea_bbox_t result = alea_bbox_empty();
    for (int ix = 0; ix < 2; ix++) {
        for (int iy = 0; iy < 2; iy++) {
            for (int iz = 0; iz < 2; iz++) {
                double x = ix ? world_clip->max_x : world_clip->min_x;
                double y = iy ? world_clip->max_y : world_clip->min_y;
                double z = iz ? world_clip->max_z : world_clip->min_z;
                alea_matrix_transform_point_inverse(accumulated, &x, &y, &z);
                alea_bbox_t p = {x, x, y, y, z, z};
                result = alea_bbox_union(&result, &p);
            }
        }
    }
    return result;
}

/* Forward declaration for the lattice helper that recurses through flatten. */
static void flatten_recursive_to_new(flatten_context_t* ctx,
                                     int universe_id,
                                     const alea_matrix_t* accumulated_transform,
                                     int depth);

/**
 * @brief Expand the elements of a lattice cell that overlap the (optional)
 * clip region, recursing into each element's fill universe.
 *
 * Conventions:
 *   - Element local coords are translated by the element origin so that
 *     each element's fill universe sees coordinates centered on the element.
 *   - The lattice cell's CSG (root_node_id) is pushed onto the parent_stack
 *     so each emitted leaf gets intersected with the lattice shell — same
 *     model as fill cells.
 *   - When ctx->config->clip_active, only elements whose center +/- pitch
 *     could overlap the clip box are visited.
 */
static void flatten_lattice_expand(flatten_context_t* ctx,
                                   const alea_cell_entry_t* cell,
                                   const alea_matrix_t* accumulated_transform,
                                   int depth) {
    int ni = cell->lat_fill_dims[1] - cell->lat_fill_dims[0] + 1;
    int nj = cell->lat_fill_dims[3] - cell->lat_fill_dims[2] + 1;
    int nk = cell->lat_fill_dims[5] - cell->lat_fill_dims[4] + 1;
    if (ni <= 0 || nj <= 0 || nk <= 0) return;
    if (!cell->lat_fill || cell->lat_fill_count == 0) return;

    /* Compute clipped index range. Defaults: full sweep. */
    int oi_lo = 0, oi_hi = ni - 1;
    int oj_lo = 0, oj_hi = nj - 1;
    int ok_lo = 0, ok_hi = nk - 1;

    if (ctx->config->clip_active) {
        const alea_bbox_t world = flatten_lattice_container_bbox(cell);
        const alea_bbox_t world_in_lattice = accumulated_transform
            ? alea_bbox_transform(&world, accumulated_transform)
            : world;
        if (!alea_bbox_intersects(&world_in_lattice, &ctx->config->clip_bbox)) {
            return;
        }

        /* Project the world-frame clip box into this lattice's local frame. */
        const alea_bbox_t local_clip =
            flatten_project_clip_to_local(&ctx->config->clip_bbox,
                                          accumulated_transform);

        if (cell->lat_type == 1) {
            /* Rect: direct division by pitch. Pad by 1 to keep boundary
             * elements when the clip box hugs an element edge. */
            const double *ll = cell->lat_lower_left;
            const double *pp = cell->lat_pitch;
            int i_lo = (int)floor((local_clip.min_x - ll[0]) / pp[0]) - 1;
            int i_hi = (int)floor((local_clip.max_x - ll[0]) / pp[0]) + 1;
            int j_lo = (nj == 1) ? 0
                : (int)floor((local_clip.min_y - ll[1]) / pp[1]) - 1;
            int j_hi = (nj == 1) ? 0
                : (int)floor((local_clip.max_y - ll[1]) / pp[1]) + 1;
            int k_lo = (nk == 1) ? 0
                : (int)floor((local_clip.min_z - ll[2]) / pp[2]) - 1;
            int k_hi = (nk == 1) ? 0
                : (int)floor((local_clip.max_z - ll[2]) / pp[2]) + 1;

            if (i_lo > oi_lo) oi_lo = i_lo;
            if (i_hi < oi_hi) oi_hi = i_hi;
            if (j_lo > oj_lo) oj_lo = j_lo;
            if (j_hi < oj_hi) oj_hi = j_hi;
            if (k_lo > ok_lo) ok_lo = k_lo;
            if (k_hi < ok_hi) ok_hi = k_hi;
        } else if (cell->lat_type == 2) {
            /* Hex: derive a conservative (ri, rk) range from the cartesian
             * corners of the local clip box.
             *   y_local = rk * p * sqrt(3) / 2  =>  rk = 2y / (p*sqrt3)
             *   x_local = ri*p + rk*p/2         =>  ri = (x - rk*p/2) / p
             * Z is rectangular. */
            double p = cell->lat_pitch[0];
            if (p > 0.0) {
                double frk_lo = local_clip.min_y / (p * M_SQRT3 * 0.5);
                double frk_hi = local_clip.max_y / (p * M_SQRT3 * 0.5);
                if (frk_lo > frk_hi) { double t = frk_lo; frk_lo = frk_hi; frk_hi = t; }
                int rk_lo_idx = (int)floor(frk_lo) - 1;
                int rk_hi_idx = (int)floor(frk_hi) + 1;

                /* For ri, sample both rk extremes and both x extremes. */
                double fri_candidates[4] = {
                    (local_clip.min_x - frk_lo * p * 0.5) / p,
                    (local_clip.max_x - frk_lo * p * 0.5) / p,
                    (local_clip.min_x - frk_hi * p * 0.5) / p,
                    (local_clip.max_x - frk_hi * p * 0.5) / p
                };
                double fri_min = fri_candidates[0], fri_max = fri_candidates[0];
                for (int c = 1; c < 4; c++) {
                    if (fri_candidates[c] < fri_min) fri_min = fri_candidates[c];
                    if (fri_candidates[c] > fri_max) fri_max = fri_candidates[c];
                }
                int ri_lo_idx = (int)floor(fri_min) - 1;
                int ri_hi_idx = (int)floor(fri_max) + 1;

                /* Convert ri/rk to oi/oj. */
                int new_oi_lo = ri_lo_idx - cell->lat_fill_dims[0];
                int new_oi_hi = ri_hi_idx - cell->lat_fill_dims[0];
                int new_oj_lo = rk_lo_idx - cell->lat_fill_dims[2];
                int new_oj_hi = rk_hi_idx - cell->lat_fill_dims[2];

                if (new_oi_lo > oi_lo) oi_lo = new_oi_lo;
                if (new_oi_hi < oi_hi) oi_hi = new_oi_hi;
                if (new_oj_lo > oj_lo) oj_lo = new_oj_lo;
                if (new_oj_hi < oj_hi) oj_hi = new_oj_hi;

                if (nk > 1) {
                    int k_lo = (int)floor((local_clip.min_z - cell->lat_lower_left[2])
                                          / cell->lat_pitch[2]) - 1;
                    int k_hi = (int)floor((local_clip.max_z - cell->lat_lower_left[2])
                                          / cell->lat_pitch[2]) + 1;
                    if (k_lo > ok_lo) ok_lo = k_lo;
                    if (k_hi < ok_hi) ok_hi = k_hi;
                }
            }
        }

        /* Clamp to valid range. */
        if (oi_lo < 0) oi_lo = 0;
        if (oj_lo < 0) oj_lo = 0;
        if (ok_lo < 0) ok_lo = 0;
        if (oi_hi >= ni) oi_hi = ni - 1;
        if (oj_hi >= nj) oj_hi = nj - 1;
        if (ok_hi >= nk) ok_hi = nk - 1;
        if (oi_lo > oi_hi || oj_lo > oj_hi || ok_lo > ok_hi) return;
    }

    /* Compose the lattice cell's own fill_transform into accumulated once,
     * if present. Each element then composes a translate on top of that. */
    alea_matrix_t lat_fill_matrix;
    bool have_lat_fill_xform = (cell->fill_transform > 0);
    if (have_lat_fill_xform) {
        const alea_transform_t* tr = alea_get_transform(ctx->src, cell->fill_transform);
        if (!tr || !alea_matrix_from_mcnp(&lat_fill_matrix, tr->cosines,
                                          tr->value_count, false)) {
            ALEA_LOG_ERROR("Invalid matrix for lattice cell %d fill_transform %d",
                    cell->mc_cell_id, cell->fill_transform);
            ctx->error = -1;
            return;
        }
    }

    /* Push the lattice shell onto the parent stack so each element's leaves
     * get intersected with the lattice cell's CSG. The shell's transform is
     * the unmodified accumulated (before per-element translation). */
    parent_geom_t parent_node;
    parent_node.node_id = cell->root_node_id;
    if (accumulated_transform) {
        parent_node.transform = *accumulated_transform;
    } else {
        alea_matrix_identity(&parent_node.transform);
    }
    parent_node.next = ctx->parent_stack;
    parent_geom_t* old_stack = ctx->parent_stack;
    ctx->parent_stack = &parent_node;

    const double p = cell->lat_pitch[0];

    for (int oi = oi_lo; oi <= oi_hi; oi++) {
        for (int oj = oj_lo; oj <= oj_hi; oj++) {
            for (int ok = ok_lo; ok <= ok_hi; ok++) {
                size_t idx = (size_t)oi * (size_t)nj * (size_t)nk
                           + (size_t)oj * (size_t)nk + (size_t)ok;
                if (idx >= cell->lat_fill_count) continue;
                int fill_univ = cell->lat_fill[idx];
                if (fill_univ <= 0) continue;

                /* Element center in lattice-local frame. */
                double ox, oy, oz;
                if (cell->lat_type == 1) {
                    ox = cell->lat_lower_left[0] + (oi + 0.5) * cell->lat_pitch[0];
                    oy = cell->lat_lower_left[1] + (oj + 0.5) * cell->lat_pitch[1];
                    oz = cell->lat_lower_left[2] + (ok + 0.5) * cell->lat_pitch[2];
                } else {
                    int ri = oi + cell->lat_fill_dims[0];
                    int rk = oj + cell->lat_fill_dims[2];
                    ox = ri * p + rk * p * 0.5;
                    oy = rk * p * M_SQRT3 * 0.5;
                    oz = (nk == 1) ? 0.0
                       : cell->lat_lower_left[2] + (ok + 0.5) * cell->lat_pitch[2];
                }

                /* Element translation: element-local -> lattice-local.
                 * alea_matrix_identity sets has_inverse=true with inv = identity;
                 * we then mutate the translation columns, so inv is now stale.
                 * Clear has_inverse so any later alea_matrix_invert recomputes
                 * the real inverse instead of accepting the stale identity. */
                alea_matrix_t translate;
                alea_matrix_identity(&translate);
                translate.m[3]  = ox;
                translate.m[7]  = oy;
                translate.m[11] = oz;
                translate.has_inverse = false;

                /* element-local -> world: accumulated × (lat_fill × translate)
                 * or accumulated × translate when no fill_transform. */
                alea_matrix_t step;
                if (have_lat_fill_xform) {
                    alea_matrix_multiply(&step, &lat_fill_matrix, &translate);
                } else {
                    step = translate;
                }

                alea_matrix_t composed;
                alea_matrix_t* composed_ptr;
                if (accumulated_transform) {
                    alea_matrix_multiply(&composed, accumulated_transform, &step);
                    composed_ptr = &composed;
                } else {
                    composed = step;
                    composed_ptr = &composed;
                }

                if (!alea_matrix_invert(composed_ptr)) {
                    ALEA_LOG_ERROR("Singular lattice element transform "
                                   "(cell %d, oi=%d oj=%d ok=%d)",
                                   cell->mc_cell_id, oi, oj, ok);
                    ctx->error = -1;
                    ctx->parent_stack = old_stack;
                    return;
                }

                flatten_recursive_to_new(ctx, fill_univ, composed_ptr, depth + 1);
                if (ctx->error) {
                    ctx->parent_stack = old_stack;
                    return;
                }
            }
        }
    }

    ctx->parent_stack = old_stack;
}

/**
 * @brief Recursively flatten cells into destination system
 */
static void flatten_recursive_to_new(flatten_context_t* ctx,
                                     int universe_id,
                                     const alea_matrix_t* accumulated_transform,
                                     int depth) {
    if (ctx->error) return;
    if (g_alea_interrupted) { ctx->error = -1; return; }
    if (ctx->config->max_depth > 0 && depth >= ctx->config->max_depth) return;
    
    const alea_universe_t* univ = alea_get_universe(ctx->src, universe_id);
    if (!univ) {
        ctx->error = -1;
        return;
    }
    /* Debug: can enable if needed
    printf("Flattening universe %d at depth %d with %zu cells\n",
           universe_id, depth, univ->cell_indices.count);
    fflush(stdout);
    */

    for (size_t i = 0; i < univ->cell_indices.count; i++) {
        size_t cell_idx = univ->cell_indices.data[i];
        
        /* Bounds check - detect corruption early */
        if (cell_idx >= alea_vec_count(&ctx->src->cells)) {
            ALEA_LOG_ERROR("Invalid cell_idx %zu (cell_count=%zu) at universe %d index %zu",
                    cell_idx, alea_vec_count(&ctx->src->cells), universe_id, i);
            ctx->error = -1;
            return;
        }
        
        const alea_cell_entry_t* cell = &ctx->src->cells.data[cell_idx];
        const bool is_lattice = (cell->lat_type != 0 && cell->lat_fill);

        /* Region-restricted flatten: skip cells whose world-frame bbox does
         * not overlap the configured clip box. For fill cells this prunes
         * the entire sub-universe; for terminal cells it avoids emitting
         * cells that lie outside the region. Conservative under rotation
         * (false positives only — no false negatives).
         *
         * Lattice cells are handled separately because the stored CSG bbox
         * may describe a single element rather than the full lattice extent
         * (MCNP LAT convention). flatten_lattice_expand does its own clip
         * test against the lattice container bbox. */
        if (ctx->config->clip_active && !is_lattice &&
            cell->root_node_id != ALEA_NODE_ID_INVALID) {
            const alea_bbox_t local = ctx->src->nodes.data[cell->root_node_id].bbox;
            const alea_bbox_t world = accumulated_transform
                ? alea_bbox_transform(&local, accumulated_transform)
                : local;
            if (!alea_bbox_intersects(&world, &ctx->config->clip_bbox)) {
                continue;
            }
        }

        if (is_lattice) {
            flatten_lattice_expand(ctx, cell, accumulated_transform, depth);
            if (ctx->error) return;
            continue;
        }

        if (cell->fill_universe > 0) {
            /* Cell has FILL - compose transform and recurse */
            alea_matrix_t fill_matrix;
            
            if (cell->fill_transform > 0) {
                const alea_transform_t* tr = alea_get_transform(ctx->src, cell->fill_transform);
                if (tr) {
                    /* Use tr->cosines which has pre-computed direction cosines */
                    if (!alea_matrix_from_mcnp(&fill_matrix, tr->cosines,
                                               tr->value_count, false)) {
                        ALEA_LOG_ERROR("Invalid matrix for cell %d fill_transform %d",
                                cell->mc_cell_id, cell->fill_transform);
                        ctx->error = -1;
                        return;
                    }
                } else {
                    ALEA_LOG_ERROR("Unknown transform %d for cell %d",
                            cell->fill_transform, cell->mc_cell_id);
                    ctx->error = -1;
                    return;
                }
            } else {
                alea_matrix_identity(&fill_matrix);
            }

            /* Compose: accumulated * fill */
            alea_matrix_t composed;
            alea_matrix_t* composed_ptr = NULL;
            if (accumulated_transform) {
                alea_matrix_multiply(&composed, accumulated_transform, &fill_matrix);
                composed_ptr = &composed;
            } else if (cell->fill_transform > 0) {
                composed = fill_matrix;
                composed_ptr = &composed;
            }
            /* else: both NULL and no fill transform → identity, pass NULL */

            if (composed_ptr && !alea_matrix_invert(composed_ptr)) {
                ALEA_LOG_ERROR("Singular matrix for cell %d fill_transform %d",
                        cell->mc_cell_id, cell->fill_transform);
                ctx->error = -1;
                return;
            }

            /* Push parent cell geometry onto stack (stack-allocated, depth bounded by max_depth) */
            parent_geom_t parent_node;
            parent_node.node_id = cell->root_node_id;

            if (accumulated_transform) {
                parent_node.transform = *accumulated_transform;
            } else {
                alea_matrix_identity(&parent_node.transform);
            }

            parent_node.next = ctx->parent_stack;

            /* Push */
            parent_geom_t* old_stack = ctx->parent_stack;
            ctx->parent_stack = &parent_node;

            /* Recurse */
            flatten_recursive_to_new(ctx, cell->fill_universe, composed_ptr, depth + 1);

            /* Pop */
            ctx->parent_stack = old_stack;
            
            /* Check for errors from recursion - must propagate to stop processing */
            if (ctx->error) {
                return;
            }   
            
        } else {
            /* Terminal cell - clone geometry with transform into dst */
            alea_node_id_t child_root;
            
            if (accumulated_transform) {
                /* Clone tree from src to dst with transform applied */
                child_root = clone_tree_to_system_transformed( 
                    ctx->dst, ctx->src, cell->root_node_id, accumulated_transform, ctx->remap);
            } else {
                /* No transform - just clone */
                child_root = alea_clone_tree_to_system(ctx->dst, ctx->src, cell->root_node_id, ctx->remap);
            }

            if (child_root == ALEA_NODE_ID_INVALID) {
                ctx->error = -1;
                return;
            }

            /* Intersect with all parent boundaries
             * NOTE: bbox optimization disabled - requires proper coordinate transform handling
             * to compare bboxes in the same coordinate system. For now, always intersect.
             */
            alea_node_id_t final_root = child_root;

            for (parent_geom_t* parent = ctx->parent_stack; parent != NULL; parent = parent->next) {
                /* Check parent clone cache first */
                alea_node_id_t parent_root = ALEA_NODE_ID_INVALID;
                for (size_t ci = 0; ci < ctx->parent_cache.count; ci++) {
                    parent_clone_entry_t* e = &ctx->parent_cache.data[ci];
                    if (e->src_node_id == parent->node_id &&
                        memcmp(e->transform_m, parent->transform.m, sizeof(e->transform_m)) == 0) {
                        parent_root = e->cloned_root;
                        break;
                    }
                }

                if (parent_root == ALEA_NODE_ID_INVALID) {
                    /* Cache miss — clone and store */
                    parent_root = clone_tree_to_system_transformed(
                        ctx->dst, ctx->src, parent->node_id, &parent->transform, ctx->remap);

                    if (parent_root == ALEA_NODE_ID_INVALID) {
                        ALEA_LOG_ERROR("Failed to clone parent node %u (src node_count=%zu)",
                                parent->node_id, alea_vec_count(&ctx->src->nodes));
                        ctx->error = -1;
                        return;
                    }

                    /* Cache the clone result */
                    parent_clone_entry_t* entry = alea_vec_push_uninit(&ctx->parent_cache, parent_clone_entry_t);
                    if (entry) {
                        entry->src_node_id = parent->node_id;
                        memcpy(entry->transform_m, parent->transform.m, sizeof(entry->transform_m));
                        entry->cloned_root = parent_root;
                    }
                }

                /* Create intersection: final = final ∩ parent */
                alea_node_id_t new_intersection = alea_create_intersection(ctx->dst, final_root, parent_root);

                if (new_intersection == ALEA_NODE_ID_INVALID) {
                    ALEA_LOG_ERROR("Failed to create intersection");
                    ctx->error = -1;
                    return;
                }

                final_root = new_intersection;
            }

            /* Add cell to destination */
            alea_cell_entry_t* new_cell = alea_vec_push_uninit(&ctx->dst->cells, alea_cell_entry_t);
            if (!new_cell) {
                ctx->error = -1;
                return;
            }
            memset(new_cell, 0, sizeof(*new_cell));

            /* Best-effort: reuse original MCNP cell ID if available and not taken */
            int preferred_id = cell->mc_cell_id;
            if (preferred_id > 0 && !cell_id_is_used(ctx, preferred_id)) {
                new_cell->mc_cell_id = preferred_id;
            } else {
                while (cell_id_is_used(ctx, ctx->next_cell_id))
                    ctx->next_cell_id++;
                new_cell->mc_cell_id = ctx->next_cell_id++;
            }
            cell_id_mark_used(ctx, new_cell->mc_cell_id);
            new_cell->root_node_id = final_root;
            new_cell->material_id = cell->material_id;
            new_cell->density = cell->density;
            new_cell->is_mass_density = cell->is_mass_density;
            new_cell->universe_id = 0;  /* All cells go to universe 0 */
            new_cell->fill_universe = 0;
            new_cell->fill_transform = 0;
            new_cell->temperature = cell->temperature;
            new_cell->has_temperature = cell->has_temperature;
        }
    }
    
}

/* Wrapper: cross-system clone without transform */
alea_node_id_t alea_clone_tree_to_system(alea_system_t* dst,
                                        const alea_system_t* src,
                                        alea_node_id_t src_root,
                                        primitive_remap_t* remap) {
    return clone_tree_impl(dst, src, src_root, NULL, remap, 0);
}

/* Cross-system clone with transform (internal) */
static alea_node_id_t clone_tree_to_system_transformed(alea_system_t* dst,
                                                    const alea_system_t* src,
                                                    alea_node_id_t src_root,
                                                    const alea_matrix_t* mat,
                                                    primitive_remap_t* remap) {
    return clone_tree_impl(dst, src, src_root, mat, remap, 0);
}

alea_system_t* alea_flatten_to_new_system(alea_system_t* src,
                                          const alea_flatten_config_t* config) {
    if (!src) return NULL;
    
    alea_flatten_config_t cfg = config ? *config : ALEA_FLATTEN_DEFAULT;
    
    /* Ensure source has universe index */
    if (!src->universe_index_built) {
        if (alea_build_universe_index(src) < 0) {
            return NULL;
        }
    }
    
    /* Create destination system */
    alea_system_t* dst = alea_system_create();
    if (!dst) return NULL;
    
    dst->config = src->config;
    dst->source = src->source;

    primitive_remap_t* remap = alea_create_remap_table(alea_vec_count(&src->primitives));
    if (!remap) {
        alea_system_destroy(dst);
        return NULL;
    }
  
    
    /* Copy materials if requested */
    if (cfg.copy_materials) {
        if (copy_materials(dst, src) != 0) {
            alea_destroy_remap_table(remap);
            alea_system_destroy(dst);
            return NULL;
        }
    }

    /* Flatten from universe 0 */
    flatten_context_t ctx = {
        .dst = dst,
        .src = src,
        .config = &cfg,
        .next_cell_id = cfg.starting_cell_id,
        .error = 0,
        .parent_stack = NULL,
        .remap = remap,
        .parent_cache = ALEA_VEC_INIT,
        .used_cell_ids = NULL,
        .used_cell_ids_size = 0
    };

    flatten_recursive_to_new(&ctx, cfg.starting_universe_id, NULL, 0);
    alea_vec_free(&ctx.parent_cache);
    free(ctx.used_cell_ids);
    alea_destroy_remap_table(remap);
    if (ctx.error || g_alea_interrupted) {
        alea_system_destroy(dst);
        alea_set_error_detail(ALEA_ERR_INTERRUPTED, "Flatten operation interrupted");
        return NULL;
    }
    /* Build universe index in destination */
    alea_build_universe_index(dst);

    ALEA_LOG_INFO("Flattened %zu cells from %zu universes into new system with %zu cells",
           alea_vec_count(&src->cells), alea_vec_count(&src->universes), alea_vec_count(&dst->cells));
    build_primitive_to_surface_map_from_nodes(dst);
    return dst;
}

int alea_flatten_in_place(alea_system_t* sys,
                         const alea_flatten_config_t* config) {
    if (!sys) return -1;

    /* Flatten to a new system */
    alea_system_t* flat = alea_flatten_to_new_system(sys, config);
    if (!flat) return -1;

    int cell_count = (int)alea_vec_count(&flat->cells);

    /* Destroy old internals, then move flat's internals into sys */
    alea_system_destroy_internals(sys);
    *sys = *flat;

    /* Free flat shell only — internals now owned by sys */
    free(flat);

    return cell_count;
}

