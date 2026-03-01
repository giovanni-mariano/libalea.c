// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_void.c
 * @brief Octree-based void region generation
 *
 * Fixed approach:
 * - Octree provides spatial decomposition (controls complexity)
 * - Void is defined analytically as complement of all material cells
 * - Octree classification is CONSERVATIVE: we only skip boxes that are
 *   DEFINITELY solid, never boxes that MIGHT contain void
 * - Each candidate box gets: regional_void = global_void ∩ region_box
 * - Simplification eliminates surfaces that don't affect each region
 * - If simplification returns INVALID (empty), the region had no void
 *
 * Key insight: the octree probing is an OPTIMIZATION to skip definitely-solid
 * regions, not the source of truth. Let the CSG math decide what's void.
 */

#include "alea_void.h"
#include "alea.h"
#include "core/alea_system.h"
#include "alea_eval.h"
#include "alea_ops.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "util/alea_log.h"
#include "util/alea_bitset.h"

// ============================================================================
// INTERNAL HELPERS
// ============================================================================

#define MERGE_TOL 1e-6

/* Forward declarations */
static bool boxes_can_form_box(const alea_bbox_t* a, const alea_bbox_t* b, alea_bbox_t* merged);

static octree_node_t* octree_node_create(const alea_bbox_t* bbox, int depth) {
    octree_node_t* node = calloc(1, sizeof(octree_node_t));
    if (!node) return NULL;
    
    node->bbox = *bbox;
    node->depth = depth;
    node->classification = OCTREE_UNKNOWN;
    node->cell_index = -1;
    
    return node;
}

void octree_destroy(octree_node_t* node) {
    if (!node) return;
    
    for (int i = 0; i < 8; i++) {
        if (node->children[i]) {
            octree_destroy(node->children[i]);
        }
    }
    free(node);
}

/**
 * Get child bounding box for given octant (0-7)
 */
static alea_bbox_t get_child_bbox(const alea_bbox_t* parent, int octant) {
    double mx = (parent->min_x + parent->max_x) / 2;
    double my = (parent->min_y + parent->max_y) / 2;
    double mz = (parent->min_z + parent->max_z) / 2;
    
    alea_bbox_t child;
    
    // octant bits: x=bit0, y=bit1, z=bit2
    child.min_x = (octant & 1) ? mx : parent->min_x;
    child.max_x = (octant & 1) ? parent->max_x : mx;
    child.min_y = (octant & 2) ? my : parent->min_y;
    child.max_y = (octant & 2) ? parent->max_y : my;
    child.min_z = (octant & 4) ? mz : parent->min_z;
    child.max_z = (octant & 4) ? parent->max_z : mz;
    
    return child;
}

/**
 * Probe points on a regular grid within a node and classify.
 *
 * CONSERVATIVE approach: we only return "definitely solid" if ALL probes
 * are in the SAME material cell. Any void probe means we must consider
 * this box as a void candidate.
 *
 * Returns:
 *   >= 0: all probes in same material cell (cell_index) - SAFE to skip
 *   -1: all probes are void - definitely contains void
 *   -2: mixed (some void, some solid, or multiple cells) - MUST check
 */
static int probe_node(const alea_system_t* sys, const alea_bbox_t* bbox,
                      int probes_per_axis, int* void_count, int* total) {
    double dx = (bbox->max_x - bbox->min_x) / (probes_per_axis - 1);
    double dy = (bbox->max_y - bbox->min_y) / (probes_per_axis - 1);
    double dz = (bbox->max_z - bbox->min_z) / (probes_per_axis - 1);

    if (dx <= 0) dx = 1e-10;
    if (dy <= 0) dy = 1e-10;
    if (dz <= 0) dz = 1e-10;

    int first_cell = -1;
    int voids = 0;
    int count = 0;
    int mixed = 0;

    for (int iz = 0; iz < probes_per_axis && !mixed; iz++) {
        double z = bbox->min_z + iz * dz;
        for (int iy = 0; iy < probes_per_axis && !mixed; iy++) {
            double y = bbox->min_y + iy * dy;
            for (int ix = 0; ix < probes_per_axis && !mixed; ix++) {
                double x = bbox->min_x + ix * dx;
                count++;

                // Find which cell contains this point
                int cell_idx = alea_identify_cell_at_point(sys, x, y, z);

                if (cell_idx < 0) {
                    voids++;
                    if (first_cell >= 0) mixed = 1;  // Was solid, now void
                } else {
                    if (first_cell < 0 && voids == 0) {
                        first_cell = cell_idx;
                    } else if (first_cell != cell_idx) {
                        mixed = 1;  // Different cells
                    } else if (voids > 0) {
                        mixed = 1;  // Was void, now solid
                    }
                }
            }
        }
    }

    *void_count = voids;
    *total = count;

    if (mixed) return -2;
    if (voids == count) return -1;  // All void
    return first_cell;  // All same cell
}

/**
 * Recursive octree building
 *
 * CONSERVATIVE: Only classify as SOLID if all probes land in the same
 * material cell. Any hint of void means we keep it as a candidate.
 */
static void build_octree_recursive(alea_system_t* sys, octree_node_t* node,
                                    const octree_config_t* config,
                                    void_result_t* result) {
    result->total_nodes++;

    // Check minimum size
    double sx = node->bbox.max_x - node->bbox.min_x;
    double sy = node->bbox.max_y - node->bbox.min_y;
    double sz = node->bbox.max_z - node->bbox.min_z;
    double min_dim = sx < sy ? (sx < sz ? sx : sz) : (sy < sz ? sy : sz);

    int probes = config->probes_per_axis;
    if (probes < 2) probes = 2;

    int void_count, total;
    int cell_result = probe_node(sys, &node->bbox, probes, &void_count, &total);

    // Classify based on probing - CONSERVATIVE approach
    if (cell_result == -1) {
        // All probes are void - definitely contains void
        node->classification = OCTREE_VOID;
        node->cell_index = -1;
        result->octree_void_count++;
        return;
    }

    if (cell_result >= 0) {
        // All probes in same material cell.
        // This is the ONLY case where we can confidently skip.
        node->classification = OCTREE_SOLID;
        node->cell_index = cell_result;
        result->solid_nodes++;
        return;
    }

    // Mixed: some void probes, or multiple cells.
    // We CANNOT discard this box - it might contain void.

    bool should_subdivide = (node->depth < config->max_depth) &&
                            (min_dim > config->min_size);

    if (!should_subdivide) {
        // At max depth/size - KEEP AS VOID CANDIDATE regardless of void fraction!
        // The old code discarded boxes with void_fraction < threshold. That was wrong.
        // Let the CSG intersection decide if there's actually void here.
        node->classification = OCTREE_PARTIAL_VOID;
        node->cell_index = -1;
        result->partial_nodes++;
        if (node->depth >= config->max_depth) {
            result->max_depth_reached++;
        }
        return;
    }

    // Subdivide for finer resolution
    node->classification = OCTREE_MIXED;

    for (int i = 0; i < 8; i++) {
        alea_bbox_t child_bbox = get_child_bbox(&node->bbox, i);
        node->children[i] = octree_node_create(&child_bbox, node->depth + 1);
        if (!node->children[i]) return;  // Allocation failure

        build_octree_recursive(sys, node->children[i], config, result);
    }
}

// ============================================================================
// BOX FROM PLANES (for surface-based export)
// ============================================================================

/**
 * Create a box region using 6 plane intersections instead of BOX macrobody.
 * This allows proper surface-based MCNP export.
 *
 * Returns the intersection node, and registers all planes as surfaces.
 */
static alea_node_id_t create_box_from_planes(
    alea_system_t* sys,
    const alea_bbox_t* box,
    int* next_surface_id
) {
    /* Create 6 plane surfaces using public API */
    alea_node_id_t planes[6];

    /* x > min_x: positive sense of plane x - min_x = 0 */
    int px0_idx = alea_plane_surface(sys, (*next_surface_id)++, 1, 0, 0, -box->min_x);
    /* x < max_x: negative sense of plane x - max_x = 0 */
    int px1_idx = alea_plane_surface(sys, (*next_surface_id)++, 1, 0, 0, -box->max_x);
    /* y > min_y */
    int py0_idx = alea_plane_surface(sys, (*next_surface_id)++, 0, 1, 0, -box->min_y);
    /* y < max_y */
    int py1_idx = alea_plane_surface(sys, (*next_surface_id)++, 0, 1, 0, -box->max_y);
    /* z > min_z */
    int pz0_idx = alea_plane_surface(sys, (*next_surface_id)++, 0, 0, 1, -box->min_z);
    /* z < max_z */
    int pz1_idx = alea_plane_surface(sys, (*next_surface_id)++, 0, 0, 1, -box->max_z);

    if (px0_idx < 0 || px1_idx < 0 || py0_idx < 0 ||
        py1_idx < 0 || pz0_idx < 0 || pz1_idx < 0) {
        return ALEA_NODE_ID_INVALID;
    }

    /* Get nodes with appropriate senses for inside box */
    planes[0] = alea_surface_at(sys, px0_idx)->pos_node;  /* x > min_x */
    planes[1] = alea_surface_at(sys, px1_idx)->neg_node;  /* x < max_x */
    planes[2] = alea_surface_at(sys, py0_idx)->pos_node;  /* y > min_y */
    planes[3] = alea_surface_at(sys, py1_idx)->neg_node;  /* y < max_y */
    planes[4] = alea_surface_at(sys, pz0_idx)->pos_node;  /* z > min_z */
    planes[5] = alea_surface_at(sys, pz1_idx)->neg_node;  /* z < max_z */

    /* Create intersection: planes[0] ∩ planes[1] ∩ ... ∩ planes[5] */
    alea_node_id_t result = alea_create_intersection(sys, planes[0], planes[1]);
    for (int i = 2; i < 6; i++) {
        result = alea_create_intersection(sys, result, planes[i]);
    }

    return result;
}

// ============================================================================
// GLOBAL VOID COMPLEMENT
// ============================================================================

/**
 * Compute global void as complement of union of all material cells.
 * void = ¬(Cell1 ∪ Cell2 ∪ ... ∪ CellN)
 *
 * This is the exact analytical definition of void - no approximation.
 */
static alea_node_id_t compute_global_void(alea_system_t* sys, int universe_id) {
    const alea_universe_t* univ = alea_get_universe(sys, universe_id);
    if (!univ || univ->cell_count == 0) {
        return ALEA_NODE_ID_INVALID;
    }

    /* Collect root nodes of all material cells */
    alea_node_id_t* cell_roots = malloc(univ->cell_count * sizeof(alea_node_id_t));
    if (!cell_roots) return ALEA_NODE_ID_INVALID;

    size_t valid_count = 0;
    for (size_t i = 0; i < univ->cell_count; i++) {
        size_t cell_idx = univ->cell_indices[i];
        const alea_cell_entry_t* cell = &sys->cells.data[cell_idx];

        /* Skip void cells (material_id == 0) and invalid cells */
        if (cell->material_id == 0) continue;
        if (cell->root_node_id == ALEA_NODE_ID_INVALID) continue;

        cell_roots[valid_count++] = cell->root_node_id;
    }

    if (valid_count == 0) {
        free(cell_roots);
        return ALEA_NODE_ID_INVALID;
    }

    /* Create union of all material cells */
    alea_node_id_t all_cells_union;
    if (valid_count == 1) {
        all_cells_union = cell_roots[0];
    } else {
        all_cells_union = alea_create_union_many(sys, cell_roots, valid_count);
    }
    free(cell_roots);

    if (all_cells_union == ALEA_NODE_ID_INVALID) {
        return ALEA_NODE_ID_INVALID;
    }

    /* Void = complement of the union */
    alea_node_id_t global_void = alea_create_complement(sys, all_cells_union);

    return global_void;
}

/**
 * Create a regional void cell by intersecting global void with a bounding box.
 * regional_void = global_void ∩ box
 *
 * Uses 6 planes (not BOX macrobody) for proper MCNP surface export.
 * Then simplify to eliminate surfaces that don't affect this region.
 */
static alea_node_id_t create_regional_void(
    alea_system_t* sys,
    alea_node_id_t global_void,
    const alea_bbox_t* box,
    int* next_surface_id
) {
    /* Create box region from 6 planes (registers surfaces automatically) */
    alea_node_id_t box_node = create_box_from_planes(sys, box, next_surface_id);

    if (box_node == ALEA_NODE_ID_INVALID) {
        return ALEA_NODE_ID_INVALID;
    }

    /* Regional void = global_void ∩ box */
    alea_node_id_t regional_void = alea_create_intersection(sys, global_void, box_node);
    return regional_void;
}

/**
 * Collect void boxes with bottom-up collapse.
 *
 * When a MIXED node has all 8 children classified as VOID leaves,
 * emit the parent's bbox instead of 8 smaller boxes. This is a
 * single O(n) pass that replaces the old depth-aware merge.
 *
 * PARTIAL_VOID nodes are always emitted individually (they need CSG verification).
 */
static int collect_and_collapse_void_boxes(
    octree_node_t* node,
    alea_bbox_t** boxes,
    size_t* count,
    size_t* capacity
) {
    if (!node) return 0;

    if (node->classification == OCTREE_VOID ||
        node->classification == OCTREE_PARTIAL_VOID) {
        /* Leaf void node — emit directly */
        if (*count >= *capacity) {
            size_t new_cap = *capacity == 0 ? 256 : (*capacity) * 2;
            alea_bbox_t* new_boxes = realloc(*boxes, new_cap * sizeof(alea_bbox_t));
            if (!new_boxes) return -1;
            *boxes = new_boxes;
            *capacity = new_cap;
        }
        (*boxes)[(*count)++] = node->bbox;
        return 0;
    }

    if (node->classification == OCTREE_MIXED) {
        /* First recurse into all 8 children */
        /* Check if all 8 children are VOID leaves (not PARTIAL_VOID) */
        bool all_void_leaves = true;
        for (int i = 0; i < 8; i++) {
            if (!node->children[i] ||
                node->children[i]->classification != OCTREE_VOID) {
                all_void_leaves = false;
                break;
            }
        }

        if (all_void_leaves) {
            /* Collapse: emit parent bbox instead of 8 children */
            if (*count >= *capacity) {
                size_t new_cap = *capacity == 0 ? 256 : (*capacity) * 2;
                alea_bbox_t* new_boxes = realloc(*boxes, new_cap * sizeof(alea_bbox_t));
                if (!new_boxes) return -1;
                *boxes = new_boxes;
                *capacity = new_cap;
            }
            (*boxes)[(*count)++] = node->bbox;
            return 0;
        }

        /* Not all void — recurse into each child */
        for (int i = 0; i < 8; i++) {
            if (node->children[i]) {
                if (collect_and_collapse_void_boxes(node->children[i], boxes, count, capacity) < 0) {
                    return -1;
                }
            }
        }
    }

    return 0;
}

/**
 * Create void CSG nodes from candidate boxes.
 *
 * For each box:
 * 1. Create regional_void = global_void ∩ box (using 6 planes)
 * 2. If result is INVALID (empty), the box had no actual void - skip it
 * 3. If result is valid, store it as a void region (CSG tree + bbox)
 *
 * This is where the CSG math decides what's actually void, not the octree.
 * Note: These are CSG nodes, not cells. Use alea_void_add_cells() to register
 * them as actual cells in the system.
 */
static int create_void_nodes_from_boxes(
    alea_system_t* sys,
    alea_node_id_t global_void,
    alea_bbox_t* boxes,
    size_t box_count,
    int start_surface_id,
    void_result_t* result
) {
    int next_surface_id = start_surface_id;

    ALEA_LOG_INFO("create_void_nodes_from_boxes: global_void=%u, box_count=%zu",
           global_void, box_count);

    for (size_t i = 0; i < box_count; i++) {
        alea_node_id_t regional_void = create_regional_void(
            sys, global_void, &boxes[i], &next_surface_id
        );

        ALEA_LOG_DEBUG("  Box %zu: [%.1f,%.1f]x[%.1f,%.1f]x[%.1f,%.1f] regional_void=%u",
               i, boxes[i].min_x, boxes[i].max_x,
               boxes[i].min_y, boxes[i].max_y,
               boxes[i].min_z, boxes[i].max_z, regional_void);

        /* If regional void is empty (no actual void in this box), skip it.
         * This is EXPECTED and normal - the octree was conservative, and
         * CSG math determined this box is actually solid. */
        if (regional_void == ALEA_NODE_ID_INVALID) {
            result->empty_regions_skipped++;
            continue;
        }

        /* Store the void region */
        if (result->void_region_count >= result->void_region_capacity) {
            size_t new_cap = result->void_region_capacity == 0 ? 64 :
                             result->void_region_capacity * 2;
            void_region_t* new_regions = realloc(result->void_regions,
                                                 new_cap * sizeof(void_region_t));
            if (!new_regions) return -1;
            result->void_regions = new_regions;
            result->void_region_capacity = new_cap;
        }

        result->void_regions[result->void_region_count].node = regional_void;
        result->void_regions[result->void_region_count].bbox = boxes[i];
        result->void_region_count++;
    }

    result->surfaces_created = next_surface_id - start_surface_id;
    return 0;
}

// ============================================================================
// PUBLIC API
// ============================================================================

void_result_t* alea_generate_void_octree(alea_system_t* sys,
                                         const alea_bbox_t* bounds,
                                         const octree_config_t* config) {
    if (!sys) return NULL;

    // Use defaults if no config
    octree_config_t cfg = config ? *config : OCTREE_DEFAULT_CONFIG;

    // Build universe index if needed
    if (!sys->universe_index_built) {
        if (alea_build_universe_index(sys) < 0) {
            return NULL;
        }
    }

    // Determine bounds
    alea_bbox_t bbox;
    if (bounds) {
        bbox = *bounds;
    } else {
        const alea_universe_t* base = alea_get_universe(sys, 0);
        if (!base) {
            ALEA_LOG_ERROR("No base universe (universe 0)");
            return NULL;
        }
        bbox = base->bbox;
        // Add margin
        double margin = 1.0;
        bbox.min_x -= margin; bbox.max_x += margin;
        bbox.min_y -= margin; bbox.max_y += margin;
        bbox.min_z -= margin; bbox.max_z += margin;
    }

    // Allocate result
    void_result_t* result = calloc(1, sizeof(void_result_t));
    if (!result) return NULL;

    // Create root node
    result->root = octree_node_create(&bbox, 0);
    if (!result->root) {
        free(result);
        return NULL;
    }

    /* Step 1: Build octree using probing for spatial classification */
    build_octree_recursive(sys, result->root, &cfg, result);

    /* Step 2: Compute global void = complement(union(all material cells)) */
    result->global_void = compute_global_void(sys, 0);
    ALEA_LOG_INFO("global_void node = %u", result->global_void);
    if (result->global_void == ALEA_NODE_ID_INVALID) {
        ALEA_LOG_WARN("could not compute global void complement");
        alea_void_print_stats(result);
        return result;
    }

    /* Step 3: Collect candidate boxes from octree with bottom-up collapse.
     * VOID siblings under a MIXED parent are collapsed into the parent box. */
    alea_bbox_t* candidate_boxes = NULL;
    size_t box_count = 0;
    size_t box_capacity = 0;

    if (collect_and_collapse_void_boxes(result->root, &candidate_boxes, &box_count, &box_capacity) < 0) {
        ALEA_LOG_WARN("void box collection incomplete");
        free(candidate_boxes);
        alea_void_print_stats(result);
        return result;
    }

    result->boxes_before_merge = box_count;
    result->boxes_after_merge = box_count;

    ALEA_LOG_INFO("Candidate void boxes from octree: %zu (after collapse)", box_count);

    /* Step 4: Create void regions using CSG intersection.
     * Each box: regional_void = global_void ∩ box (6 planes) */

    /* Find next available surface ID */
    int start_surface_id = 1;
    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        if (sys->surfaces.data[i].mcnp_surface_id >= start_surface_id) {
            start_surface_id = sys->surfaces.data[i].mcnp_surface_id + 1;
        }
    }

    if (create_void_nodes_from_boxes(sys, result->global_void,
                                      candidate_boxes, box_count,
                                      start_surface_id,
                                      result) < 0) {
        ALEA_LOG_WARN("void cell creation incomplete");
    }

    free(candidate_boxes);

    ALEA_LOG_INFO("Void CSG processing:");
    ALEA_LOG_INFO("  Candidate boxes checked: %zu", box_count);
    ALEA_LOG_INFO("  Actual void regions found: %zu", result->void_region_count);
    ALEA_LOG_INFO("  Empty candidates (octree was conservative): %zu", result->empty_regions_skipped);
    ALEA_LOG_INFO("  Plane surfaces created: %zu", result->surfaces_created);

    alea_void_print_stats(result);

    return result;
}

alea_node_id_t alea_void_to_node(alea_system_t* sys, const void_result_t* result) {
    if (!sys || !result || !result->void_regions || result->void_region_count == 0) {
        return ALEA_NODE_ID_INVALID;
    }

    if (result->void_region_count == 1) {
        return result->void_regions[0].node;
    }

    /* Build array of node IDs for union */
    alea_node_id_t* nodes = malloc(result->void_region_count * sizeof(alea_node_id_t));
    if (!nodes) return ALEA_NODE_ID_INVALID;

    for (size_t i = 0; i < result->void_region_count; i++) {
        nodes[i] = result->void_regions[i].node;
    }

    alea_node_id_t result_node = alea_create_union_many(sys, nodes, result->void_region_count);
    free(nodes);

    return result_node;
}

void alea_void_result_destroy(void_result_t* result) {
    if (!result) return;

    octree_destroy(result->root);
    free(result->void_regions);
    free(result);
}

size_t alea_void_box_count(const void_result_t* result) {
    return result ? result->void_region_count : 0;
}

int alea_void_box_get(const void_result_t* result, size_t index, alea_bbox_t* box) {
    if (!result || !box || index >= result->void_region_count) return -1;
    *box = result->void_regions[index].bbox;
    return 0;
}

void alea_void_print_stats(const void_result_t* result) {
    if (!result) return;

    ALEA_LOG_INFO("Octree void generation:");
    ALEA_LOG_INFO("  Total octree nodes: %zu", result->total_nodes);
    ALEA_LOG_INFO("  Void octree nodes: %zu", result->octree_void_count);
    ALEA_LOG_INFO("  Solid nodes: %zu", result->solid_nodes);
    ALEA_LOG_INFO("  Partial void nodes: %zu", result->partial_nodes);
    ALEA_LOG_INFO("  Max depth reached: %zu times", result->max_depth_reached);

    if (result->void_regions && result->void_region_count > 0) {
        ALEA_LOG_INFO("  Boxes before collapse: %zu", result->boxes_before_merge);
        ALEA_LOG_INFO("  Boxes after collapse: %zu (%.1f%% reduction)",
               result->boxes_after_merge,
               result->boxes_before_merge > 0 ?
                   100.0 * (1.0 - (double)result->boxes_after_merge / result->boxes_before_merge) : 0.0);
        ALEA_LOG_INFO("  Void cells created: %zu (using plane surfaces)", result->void_region_count);
        ALEA_LOG_INFO("  Empty regions skipped: %zu", result->empty_regions_skipped);
        ALEA_LOG_INFO("  Plane surfaces created: %zu", result->surfaces_created);
    }

    /* Volume estimation */
    if (result->void_region_count > 0 && result->root) {
        double total_vol = 0;
        for (size_t i = 0; i < result->void_region_count; i++) {
            const alea_bbox_t* b = &result->void_regions[i].bbox;
            total_vol += (b->max_x - b->min_x) *
                         (b->max_y - b->min_y) *
                         (b->max_z - b->min_z);
        }

        const alea_bbox_t* r = &result->root->bbox;
        double total_bbox_vol = (r->max_x - r->min_x) *
                                (r->max_y - r->min_y) *
                                (r->max_z - r->min_z);

        ALEA_LOG_INFO("  Void region volume: %.2f (%.1f%% of bounding box)",
               total_vol, 100.0 * total_vol / total_bbox_vol);
    }
}

int alea_void_add_cells(alea_system_t* sys, void_result_t* result) {
    if (!sys || !result) return -1;
    if (!result->void_regions || result->void_region_count == 0) return 0;

    int next_cell_id = alea_max_cell_id(sys) + 1;

    int cells_added = 0;
    for (size_t i = 0; i < result->void_region_count; i++) {
        ALEA_CHECK_INTERRUPTED(-1);
        if (result->void_regions[i].node != ALEA_NODE_ID_INVALID) {
            int cell_idx = alea_add_cell(sys, next_cell_id++, result->void_regions[i].node, ALEA_MATERIAL_VOID, 0.0, 0);
            if (cell_idx >= 0) {
                cells_added++;
            }
        }
    }

    return cells_added;
}

/* ============================================================================
 * VOID CELL MERGING
 * ============================================================================ */

const alea_void_merge_config_t ALEA_VOID_MERGE_DEFAULT = {
    .cell_weight = 1.0,
    .surface_weight = 0.1,
    .max_surfaces_per_cell = 24,
    .min_cells = 1
};

/* Check if two boxes are adjacent (share a face) */
static bool boxes_adjacent(const alea_bbox_t* a, const alea_bbox_t* b) {
    /* Check X adjacency: a.max_x touches b.min_x or vice versa */
    bool x_touch = (fabs(a->max_x - b->min_x) < MERGE_TOL) ||
                   (fabs(b->max_x - a->min_x) < MERGE_TOL);
    bool x_overlap = (a->min_x < b->max_x - MERGE_TOL) &&
                     (b->min_x < a->max_x - MERGE_TOL);

    bool y_touch = (fabs(a->max_y - b->min_y) < MERGE_TOL) ||
                   (fabs(b->max_y - a->min_y) < MERGE_TOL);
    bool y_overlap = (a->min_y < b->max_y - MERGE_TOL) &&
                     (b->min_y < a->max_y - MERGE_TOL);

    bool z_touch = (fabs(a->max_z - b->min_z) < MERGE_TOL) ||
                   (fabs(b->max_z - a->min_z) < MERGE_TOL);
    bool z_overlap = (a->min_z < b->max_z - MERGE_TOL) &&
                     (b->min_z < a->max_z - MERGE_TOL);

    /* Adjacent = touch on one axis, overlap on other two */
    return (x_touch && y_overlap && z_overlap) ||
           (y_touch && x_overlap && z_overlap) ||
           (z_touch && x_overlap && y_overlap);
}

/* Check if two boxes can merge into a larger box (aligned on 2 axes) */
static bool boxes_can_form_box(const alea_bbox_t* a, const alea_bbox_t* b,
                                alea_bbox_t* merged) {
    /* Try each axis as the merge axis */

    /* Merge along X: same Y and Z extents */
    if (fabs(a->min_y - b->min_y) < MERGE_TOL &&
        fabs(a->max_y - b->max_y) < MERGE_TOL &&
        fabs(a->min_z - b->min_z) < MERGE_TOL &&
        fabs(a->max_z - b->max_z) < MERGE_TOL) {
        /* Check if adjacent along X */
        if (fabs(a->max_x - b->min_x) < MERGE_TOL ||
            fabs(b->max_x - a->min_x) < MERGE_TOL) {
            merged->min_x = fmin(a->min_x, b->min_x);
            merged->max_x = fmax(a->max_x, b->max_x);
            merged->min_y = a->min_y;
            merged->max_y = a->max_y;
            merged->min_z = a->min_z;
            merged->max_z = a->max_z;
            return true;
        }
    }

    /* Merge along Y: same X and Z extents */
    if (fabs(a->min_x - b->min_x) < MERGE_TOL &&
        fabs(a->max_x - b->max_x) < MERGE_TOL &&
        fabs(a->min_z - b->min_z) < MERGE_TOL &&
        fabs(a->max_z - b->max_z) < MERGE_TOL) {
        if (fabs(a->max_y - b->min_y) < MERGE_TOL ||
            fabs(b->max_y - a->min_y) < MERGE_TOL) {
            merged->min_x = a->min_x;
            merged->max_x = a->max_x;
            merged->min_y = fmin(a->min_y, b->min_y);
            merged->max_y = fmax(a->max_y, b->max_y);
            merged->min_z = a->min_z;
            merged->max_z = a->max_z;
            return true;
        }
    }

    /* Merge along Z: same X and Y extents */
    if (fabs(a->min_x - b->min_x) < MERGE_TOL &&
        fabs(a->max_x - b->max_x) < MERGE_TOL &&
        fabs(a->min_y - b->min_y) < MERGE_TOL &&
        fabs(a->max_y - b->max_y) < MERGE_TOL) {
        if (fabs(a->max_z - b->min_z) < MERGE_TOL ||
            fabs(b->max_z - a->min_z) < MERGE_TOL) {
            merged->min_x = a->min_x;
            merged->max_x = a->max_x;
            merged->min_y = a->min_y;
            merged->max_y = a->max_y;
            merged->min_z = fmin(a->min_z, b->min_z);
            merged->max_z = fmax(a->max_z, b->max_z);
            return true;
        }
    }

    return false;
}

/* Count surfaces in a CSG node (approximate) */
static int count_surfaces(const alea_system_t* sys, alea_node_id_t node) {
    if (node == ALEA_NODE_ID_INVALID || node >= alea_vec_count(&sys->nodes)) return 0;

    const alea_node_t* n = &sys->nodes.data[node];
    alea_operation_t op = ALEA_GET_OPERATION(n);

    if (op == ALEA_OP_PRIMITIVE) {
        return 1;
    } else if (op == ALEA_OP_COMPLEMENT) {
        return count_surfaces(sys, n->operation.left);
    } else {
        return count_surfaces(sys, n->operation.left) +
               count_surfaces(sys, n->operation.right);
    }
}

int alea_merge_void_cells(alea_system_t* sys,
                         void_result_t* result,
                         const alea_void_merge_config_t* config) {
    if (!sys || !result) return -1;
    if (result->void_region_count == 0) return 0;

    const alea_void_merge_config_t* cfg = config ? config : &ALEA_VOID_MERGE_DEFAULT;

    size_t n = result->void_region_count;

    /* Track which cells are still alive (not merged into another) */
    alea_bitset_t alive = alea_bitset_create(n);
    /* Track surface count for each cell */
    int* surfaces = calloc(n, sizeof(int));

    if (!alive.words || !surfaces) {
        alea_bitset_destroy(&alive);
        free(surfaces);
        return -1;
    }

    /* Initialize: all alive, count surfaces */
    size_t alive_count = n;
    for (size_t i = 0; i < n; i++) {
        alea_bitset_set(&alive, i);
        surfaces[i] = count_surfaces(sys, result->void_regions[i].node);
    }

    ALEA_LOG_INFO("Void merge: starting with %zu cells", n);

    /* Greedy merge loop */
    bool improved = true;
    while (improved && alive_count > (size_t)cfg->min_cells) {
        improved = false;

        double best_delta = 0;
        size_t best_i = 0, best_j = 0;
        bool best_is_box = false;

        /* Find best merge candidate */
        for (size_t i = 0; i < n; i++) {
            if (!alea_bitset_test(&alive, i)) continue;
            if (g_alea_interrupted) { alea_bitset_destroy(&alive); free(surfaces); ALEA_CHECK_INTERRUPTED(-1); }

            for (size_t j = i + 1; j < n; j++) {
                if (!alea_bitset_test(&alive, j)) continue;

                /* Check adjacency */
                if (!boxes_adjacent(&result->void_regions[i].bbox,
                                    &result->void_regions[j].bbox)) {
                    continue;
                }

                /* Check if aligned merge produces a box (6 surfaces) */
                alea_bbox_t merged_bbox;
                bool is_box = boxes_can_form_box(&result->void_regions[i].bbox,
                                                  &result->void_regions[j].bbox,
                                                  &merged_bbox);
                int new_surfaces = is_box ? 6 : surfaces[i] + surfaces[j];

                /* Check complexity limit */
                if (new_surfaces > cfg->max_surfaces_per_cell) {
                    continue;
                }

                /* Calculate cost delta */
                double before_cost = 2 * cfg->cell_weight +
                                    (surfaces[i] + surfaces[j]) * cfg->surface_weight;
                double after_cost = 1 * cfg->cell_weight +
                                   new_surfaces * cfg->surface_weight;
                double delta = after_cost - before_cost;

                if (delta < best_delta) {
                    best_delta = delta;
                    best_i = i;
                    best_j = j;
                    best_is_box = is_box;
                }
            }
        }

        /* Perform best merge if beneficial */
        if (best_delta < -1e-9) {
            improved = true;

            ALEA_LOG_INFO("  Merging cells %zu + %zu (delta=%.2f)",
                   best_i, best_j, best_delta);

            /* Create merged CSG node as union */
            alea_node_id_t merged_node = alea_create_union(sys,
                result->void_regions[best_i].node,
                result->void_regions[best_j].node);

            /* Update bounding box to encompass both */
            alea_bbox_t* bi = &result->void_regions[best_i].bbox;
            const alea_bbox_t* bj = &result->void_regions[best_j].bbox;
            bi->min_x = fmin(bi->min_x, bj->min_x);
            bi->max_x = fmax(bi->max_x, bj->max_x);
            bi->min_y = fmin(bi->min_y, bj->min_y);
            bi->max_y = fmax(bi->max_y, bj->max_y);
            bi->min_z = fmin(bi->min_z, bj->min_z);
            bi->max_z = fmax(bi->max_z, bj->max_z);

            /* A box merge always has 6 surfaces; union keeps the sum */
            surfaces[best_i] = best_is_box ? 6 : surfaces[best_i] + surfaces[best_j];

            /* Update cell i with merged result */
            result->void_regions[best_i].node = merged_node;

            /* Mark cell j as dead */
            alea_bitset_clear(&alive, best_j);
            alive_count--;
        }
    }

    /* Compact the array to remove dead cells */
    size_t write_idx = 0;
    for (size_t read_idx = 0; read_idx < n; read_idx++) {
        if (alea_bitset_test(&alive, read_idx)) {
            if (write_idx != read_idx) {
                result->void_regions[write_idx] = result->void_regions[read_idx];
            }
            write_idx++;
        }
    }
    result->void_region_count = write_idx;

    ALEA_LOG_INFO("Void merge: reduced to %zu cells (%.1f%% reduction)",
           result->void_region_count,
           100.0 * (1.0 - (double)result->void_region_count / n));

    alea_bitset_destroy(&alive);
    free(surfaces);

    return (int)result->void_region_count;
}