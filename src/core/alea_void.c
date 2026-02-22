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
 * Key insight: the octree sampling is an OPTIMIZATION to skip definitely-solid
 * regions, not the source of truth. Let the CSG math decide what's void.
 */

#include "alea_void.h"
#include "alea.h"
#include "core/alea_system.h"
#include "alea_eval.h"
#include "alea_ops.h"
#include "alea_simplify.h"
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
 * Sample points in a node and classify.
 *
 * CONSERVATIVE approach: we only return "definitely solid" if ALL samples
 * are in the SAME material cell. Any void sample means we must consider
 * this box as a void candidate.
 *
 * Returns:
 *   >= 0: all samples in same material cell (cell_index) - SAFE to skip
 *   -1: all samples are void - definitely contains void
 *   -2: mixed (some void, some solid, or multiple cells) - MUST check
 */
static int sample_node(const alea_system_t* sys, const alea_bbox_t* bbox,
                       int samples_per_axis, int* void_count, int* total) {
    double dx = (bbox->max_x - bbox->min_x) / (samples_per_axis - 1);
    double dy = (bbox->max_y - bbox->min_y) / (samples_per_axis - 1);
    double dz = (bbox->max_z - bbox->min_z) / (samples_per_axis - 1);

    if (dx <= 0) dx = 1e-10;
    if (dy <= 0) dy = 1e-10;
    if (dz <= 0) dz = 1e-10;

    int first_cell = -1;
    int voids = 0;
    int count = 0;
    int mixed = 0;

    for (int iz = 0; iz < samples_per_axis && !mixed; iz++) {
        double z = bbox->min_z + iz * dz;
        for (int iy = 0; iy < samples_per_axis && !mixed; iy++) {
            double y = bbox->min_y + iy * dy;
            for (int ix = 0; ix < samples_per_axis && !mixed; ix++) {
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
 * CONSERVATIVE: Only classify as SOLID if we have HIGH CONFIDENCE that
 * the entire box is inside a single material cell. Any hint of void
 * means we keep it as a candidate.
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

    int samples = (config->samples_per_node <= 27) ? 3 :
                  (int)cbrt(config->samples_per_node);
    if (samples < 2) samples = 2;

    int void_count, total;
    int cell_result = sample_node(sys, &node->bbox, samples, &void_count, &total);

    // Classify based on sampling - CONSERVATIVE approach
    if (cell_result == -1) {
        // All samples are void - definitely contains void
        node->classification = OCTREE_VOID;
        node->cell_index = -1;
        result->octree_void_count++;
        return;
    }

    if (cell_result >= 0) {
        // All samples in same material cell.
        // This is the ONLY case where we can confidently skip.
        node->classification = OCTREE_SOLID;
        node->cell_index = cell_result;
        result->solid_nodes++;
        return;
    }

    // Mixed: some void samples, or multiple cells.
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

/**
 * Collect void boxes from octree leaves
 */
static int collect_void_boxes(octree_node_t* node, void_result_t* result) {
    if (!node) return 0;

    if (node->classification == OCTREE_VOID ||
        node->classification == OCTREE_PARTIAL_VOID) {
        /* Add this box to void collection */
        if (result->void_box_count >= result->void_box_capacity) {
            size_t new_cap = result->void_box_capacity == 0 ? 64 :
                             result->void_box_capacity * 2;
            alea_bbox_t* new_boxes = realloc(result->void_boxes,
                                             new_cap * sizeof(alea_bbox_t));
            if (!new_boxes) {
                /* Allocation failed - keep existing data, signal error */
                return -1;
            }
            result->void_boxes = new_boxes;
            result->void_box_capacity = new_cap;
        }
        result->void_boxes[result->void_box_count++] = node->bbox;
        return 0;
    }

    /* Recurse into children */
    if (node->classification == OCTREE_MIXED) {
        for (int i = 0; i < 8; i++) {
            if (node->children[i]) {
                if (collect_void_boxes(node->children[i], result) < 0) {
                    return -1;  /* Propagate allocation failure */
                }
            }
        }
    }
    return 0;
}

// ============================================================================
// VOID BOX MERGING (DISABLED)
// ============================================================================

/*
 * Box merging has been disabled because it can cause correctness issues:
 *
 * Problem: Merging a VOID box with a PARTIAL_VOID box can create a combined
 * box that covers regions the octree correctly classified as SOLID.
 *
 * The CSG simplification step already handles complexity reduction by
 * pruning surfaces that don't affect each region. Box merging is not
 * needed and only adds risk.
 *
 * If performance requires merging in the future, track classification
 * per-box and only merge boxes with the same classification.
 */

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
    int* next_surface_id,
    alea_simplify_stats_t* stats
) {
    /* Create box region from 6 planes (registers surfaces automatically) */
    alea_node_id_t box_node = create_box_from_planes(sys, box, next_surface_id);

    if (box_node == ALEA_NODE_ID_INVALID) {
        return ALEA_NODE_ID_INVALID;
    }

    /* Regional void = global_void ∩ box */
    alea_node_id_t regional_void = alea_create_intersection(sys, global_void, box_node);
    if (regional_void == ALEA_NODE_ID_INVALID) {
        return ALEA_NODE_ID_INVALID;
    }

    /* Skip simplification - just return the regional void as-is.
     * The simplifier has bugs with void expressions that cause valid void
     * to be incorrectly detected as empty. */
    (void)stats;
    return regional_void;
}

/* Box with depth for merging */
typedef struct {
    alea_bbox_t bbox;
    int depth;
} depth_box_t;

/**
 * Collect void boxes from octree with depth info.
 */
static int collect_octree_void_boxes_with_depth(
    octree_node_t* node,
    depth_box_t** boxes,
    size_t* count,
    size_t* capacity
) {
    if (!node) return 0;

    if (node->classification == OCTREE_VOID ||
        node->classification == OCTREE_PARTIAL_VOID) {

        /* Grow array if needed */
        if (*count >= *capacity) {
            size_t new_cap = *capacity == 0 ? 256 : (*capacity) * 2;
            depth_box_t* new_boxes = realloc(*boxes, new_cap * sizeof(depth_box_t));
            if (!new_boxes) return -1;
            *boxes = new_boxes;
            *capacity = new_cap;
        }

        (*boxes)[*count].bbox = node->bbox;
        (*boxes)[*count].depth = node->depth;
        (*count)++;
        return 0;
    }

    if (node->classification == OCTREE_MIXED) {
        for (int i = 0; i < 8; i++) {
            if (node->children[i]) {
                if (collect_octree_void_boxes_with_depth(node->children[i], boxes, count, capacity) < 0) {
                    return -1;
                }
            }
        }
    }

    return 0;
}

/**
 * Merge boxes at the same depth level.
 * Greedy algorithm: repeatedly find pairs that can form a larger box.
 */
static void merge_boxes_by_depth(depth_box_t* boxes, size_t* count) {
    if (*count <= 1) return;

    /* Find max depth */
    int max_depth = 0;
    for (size_t i = 0; i < *count; i++) {
        if (boxes[i].depth > max_depth) max_depth = boxes[i].depth;
    }

    /* Merge at each depth level, starting from deepest */
    for (int d = max_depth; d >= 0; d--) {
        bool merged_any = true;
        while (merged_any) {
            merged_any = false;
            for (size_t i = 0; i < *count && !merged_any; i++) {
                if (boxes[i].depth != d) continue;
                for (size_t j = i + 1; j < *count && !merged_any; j++) {
                    if (boxes[j].depth != d) continue;

                    alea_bbox_t merged;
                    if (boxes_can_form_box(&boxes[i].bbox, &boxes[j].bbox, &merged)) {
                        /* Replace box i with merged, remove box j */
                        boxes[i].bbox = merged;
                        /* Keep depth the same - merged box is at same level */
                        boxes[j] = boxes[*count - 1];
                        (*count)--;
                        merged_any = true;
                    }
                }
            }
        }
    }

    ALEA_LOG_INFO("After depth-aware merge: %zu boxes", *count);
}

/**
 * Create void CSG nodes from candidate boxes.
 *
 * For each box:
 * 1. Create regional_void = global_void ∩ box (using 6 planes)
 * 2. Simplify using interval arithmetic
 * 3. If result is INVALID (empty), the box had no actual void - skip it
 * 4. If result is valid, store it as a void node (CSG tree)
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
    void_result_t* result,
    alea_simplify_stats_t* total_stats
) {
    int next_surface_id = start_surface_id;

    ALEA_LOG_INFO("create_void_nodes_from_boxes: global_void=%u, box_count=%zu",
           global_void, box_count);

    for (size_t i = 0; i < box_count; i++) {
        alea_simplify_stats_t stats = {0};
        alea_node_id_t regional_void = create_regional_void(
            sys, global_void, &boxes[i], &next_surface_id, &stats
        );

        if (i < 3) {
            ALEA_LOG_INFO("  Box %zu: [%.1f,%.1f]x[%.1f,%.1f]x[%.1f,%.1f] regional_void=%u",
                   i, boxes[i].min_x, boxes[i].max_x,
                   boxes[i].min_y, boxes[i].max_y,
                   boxes[i].min_z, boxes[i].max_z, regional_void);
        }

        /* Accumulate stats */
        if (total_stats) {
            total_stats->absorption_reductions += stats.absorption_reductions;
            total_stats->contradictions_found += stats.contradictions_found;
            total_stats->complements_eliminated += stats.complements_eliminated;
        }

        /* If regional void is empty (no actual void in this box), skip it.
         * This is EXPECTED and normal - the octree was conservative, and
         * CSG math determined this box is actually solid. */
        if (regional_void == ALEA_NODE_ID_INVALID) {
            result->empty_regions_skipped++;
            continue;
        }

        /* Store the void cell and its bounding box */
        if (result->void_node_count >= result->void_node_capacity) {
            size_t new_cap = result->void_node_capacity == 0 ? 64 :
                             result->void_node_capacity * 2;
            alea_node_id_t* new_cells = realloc(result->void_nodes,
                                                new_cap * sizeof(alea_node_id_t));
            if (!new_cells) return -1;
            result->void_nodes = new_cells;

            alea_bbox_t* new_bboxes = realloc(result->void_boxes,
                                              new_cap * sizeof(alea_bbox_t));
            if (!new_bboxes) return -1;
            result->void_boxes = new_bboxes;

            result->void_node_capacity = new_cap;
            result->void_box_capacity = new_cap;  // Keep in sync
        }

        result->void_nodes[result->void_node_count] = regional_void;
        result->void_boxes[result->void_node_count] = boxes[i];
        result->void_node_count++;
        result->void_box_count++;  // Keep in sync with void_node_count
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

    /* Step 1: Build octree using sampling for spatial classification */
    build_octree_recursive(sys, result->root, &cfg, result);

    /* Step 2: Compute global void = complement(union(all material cells)) */
    result->global_void = compute_global_void(sys, 0);
    ALEA_LOG_INFO("global_void node = %u", result->global_void);
    if (result->global_void == ALEA_NODE_ID_INVALID) {
        ALEA_LOG_WARN("could not compute global void complement");
        /* Fall back to old box-based approach */
        if (collect_void_boxes(result->root, result) < 0) {
            ALEA_LOG_WARN("void box collection incomplete");
        }
    } else {
        /* Step 3: Collect candidate boxes from octree leaves with depth info.
         * These include VOID and PARTIAL_VOID nodes - anything that might have void. */
        depth_box_t* depth_boxes = NULL;
        size_t box_count = 0;
        size_t box_capacity = 0;

        if (collect_octree_void_boxes_with_depth(result->root, &depth_boxes, &box_count, &box_capacity) < 0) {
            ALEA_LOG_WARN("void box collection incomplete");
            free(depth_boxes);
        } else {
            result->boxes_before_merge = box_count;

            ALEA_LOG_INFO("Candidate void boxes from octree: %zu", box_count);

            /* Step 3b: Merge boxes at the same depth level.
             * Boxes at the same depth have the same size and can merge into
             * larger boxes if adjacent and aligned. */
            merge_boxes_by_depth(depth_boxes, &box_count);

            result->boxes_after_merge = box_count;

            /* Extract just the bboxes for CSG creation */
            alea_bbox_t* candidate_boxes = malloc(box_count * sizeof(alea_bbox_t));
            if (!candidate_boxes) {
                ALEA_LOG_WARN("malloc failed for candidate_boxes");
                free(depth_boxes);
                return result;
            }
            for (size_t i = 0; i < box_count; i++) {
                candidate_boxes[i] = depth_boxes[i].bbox;
            }
            free(depth_boxes);

            ALEA_LOG_INFO("Boxes after merge: %zu (%.1f%% reduction)",
                   box_count,
                   result->boxes_before_merge > 0 ?
                       100.0 * (1.0 - (double)box_count / result->boxes_before_merge) : 0.0);

            /* Step 4: Create void cells using CSG intersection.
             * Each box: regional_void = global_void ∩ box (6 planes)
             * CSG simplification determines if there's actually void there. */

            /* Find next available surface ID */
            int start_surface_id = 1;
            for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
                if (sys->surfaces.data[i].mcnp_surface_id >= start_surface_id) {
                    start_surface_id = sys->surfaces.data[i].mcnp_surface_id + 1;
                }
            }

            alea_simplify_stats_t simplify_stats = {0};
            if (create_void_nodes_from_boxes(sys, result->global_void,
                                              candidate_boxes, box_count,
                                              start_surface_id,
                                              result, &simplify_stats) < 0) {
                ALEA_LOG_WARN("void cell creation incomplete");
            }

            free(candidate_boxes);

            ALEA_LOG_INFO("Void CSG processing:");
            ALEA_LOG_INFO("  Candidate boxes checked: %zu", box_count);
            ALEA_LOG_INFO("  Actual void regions found: %zu", result->void_node_count);
            ALEA_LOG_INFO("  Empty candidates (octree was conservative): %zu", result->empty_regions_skipped);
            ALEA_LOG_INFO("  Surfaces pruned by simplification: %zu", simplify_stats.absorption_reductions);
            ALEA_LOG_INFO("  Plane surfaces created: %zu", result->surfaces_created);
        }
    }

    alea_void_print_stats(result);

    return result;
}

alea_node_id_t alea_void_to_node(alea_system_t* sys, const void_result_t* result) {
    if (!sys || !result) {
        return ALEA_NODE_ID_INVALID;
    }

    /* Prefer the new void_cells array (proper CSG definitions) */
    if (result->void_nodes && result->void_node_count > 0) {
        if (result->void_node_count == 1) {
            return result->void_nodes[0];
        }
        return alea_create_union_many(sys, result->void_nodes, result->void_node_count);
    }

    /* Fallback: legacy box-based approach */
    if (result->void_box_count == 0) {
        return ALEA_NODE_ID_INVALID;
    }

    alea_node_id_t* box_nodes = malloc(result->void_box_count * sizeof(alea_node_id_t));
    if (!box_nodes) return ALEA_NODE_ID_INVALID;

    /* Start surface IDs after existing surfaces to avoid conflicts */
    int next_surf_id = (int)alea_vec_count(&sys->surfaces) + 90000;

    size_t valid = 0;
    for (size_t i = 0; i < result->void_box_count; i++) {
        const alea_bbox_t* b = &result->void_boxes[i];
        int box_idx = alea_box_surface(sys, next_surf_id++,
            b->min_x, b->max_x, b->min_y, b->max_y, b->min_z, b->max_z);
        if (box_idx >= 0) {
            box_nodes[valid++] = alea_surface_at(sys, box_idx)->neg_node;
        }
    }

    if (valid == 0) {
        free(box_nodes);
        return ALEA_NODE_ID_INVALID;
    }

    alea_node_id_t result_node = alea_create_union_many(sys, box_nodes, valid);
    free(box_nodes);

    return result_node;
}

void alea_void_result_destroy(void_result_t* result) {
    if (!result) return;

    octree_destroy(result->root);
    free(result->void_boxes);
    free(result->void_nodes);
    free(result);
}

size_t alea_void_box_count(const void_result_t* result) {
    return result ? result->void_box_count : 0;
}

int alea_void_box_get(const void_result_t* result, size_t index, alea_bbox_t* box) {
    if (!result || !box || index >= result->void_box_count) return -1;
    *box = result->void_boxes[index];
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

    /* New hybrid approach stats */
    if (result->void_nodes && result->void_node_count > 0) {
        ALEA_LOG_INFO("  Boxes before merge: %zu", result->boxes_before_merge);
        ALEA_LOG_INFO("  Boxes after merge: %zu (%.1f%% reduction)",
               result->boxes_after_merge,
               result->boxes_before_merge > 0 ?
                   100.0 * (1.0 - (double)result->boxes_after_merge / result->boxes_before_merge) : 0.0);
        ALEA_LOG_INFO("  Void cells created: %zu (using plane surfaces)", result->void_node_count);
        ALEA_LOG_INFO("  Empty regions skipped: %zu", result->empty_regions_skipped);
        ALEA_LOG_INFO("  Plane surfaces created: %zu", result->surfaces_created);
    } else {
        ALEA_LOG_INFO("  Void boxes collected: %zu (legacy mode)", result->void_box_count);
    }

    /* Volume estimation */
    size_t box_count = result->void_node_count > 0 ? result->void_node_count : result->void_box_count;
    if (box_count > 0 && result->root) {
        double total_vol = 0;
        for (size_t i = 0; i < box_count; i++) {
            const alea_bbox_t* b = &result->void_boxes[i];
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
    if (!result->void_nodes || result->void_node_count == 0) return 0;

    int next_cell_id = alea_max_cell_id(sys) + 1;

    int cells_added = 0;
    for (size_t i = 0; i < result->void_node_count; i++) {
        ALEA_CHECK_INTERRUPTED(-1);
        if (result->void_nodes[i] != ALEA_NODE_ID_INVALID) {
            int cell_idx = alea_add_cell(sys, next_cell_id++, result->void_nodes[i], 0, 0.0, 0);
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
    if (result->void_node_count == 0) return 0;

    const alea_void_merge_config_t* cfg = config ? config : &ALEA_VOID_MERGE_DEFAULT;

    size_t n = result->void_node_count;

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
        surfaces[i] = count_surfaces(sys, result->void_nodes[i]);
    }

    ALEA_LOG_INFO("Void merge: starting with %zu cells", n);

    /* Greedy merge loop */
    bool improved = true;
    while (improved && alive_count > (size_t)cfg->min_cells) {
        improved = false;

        double best_delta = 0;
        size_t best_i = 0, best_j = 0;

        /* Find best merge candidate */
        for (size_t i = 0; i < n; i++) {
            if (!alea_bitset_test(&alive, i)) continue;
            if (g_alea_interrupted) { alea_bitset_destroy(&alive); free(surfaces); ALEA_CHECK_INTERRUPTED(-1); }

            for (size_t j = i + 1; j < n; j++) {
                if (!alea_bitset_test(&alive, j)) continue;

                /* Check adjacency */
                if (!boxes_adjacent(&result->void_boxes[i], &result->void_boxes[j])) {
                    continue;
                }

                /* Calculate merge benefit - union of two cells */
                int new_surfaces = surfaces[i] + surfaces[j];

                /* Check complexity limit */
                if (new_surfaces > cfg->max_surfaces_per_cell) {
                    continue;
                }

                /* Calculate cost delta */
                /* Before: 2 cells with surfaces[i] + surfaces[j] surfaces */
                /* After: 1 cell with new_surfaces surfaces */
                double before_cost = 2 * cfg->cell_weight +
                                    (surfaces[i] + surfaces[j]) * cfg->surface_weight;
                double after_cost = 1 * cfg->cell_weight +
                                   new_surfaces * cfg->surface_weight;
                double delta = after_cost - before_cost;

                if (delta < best_delta) {
                    best_delta = delta;
                    best_i = i;
                    best_j = j;
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
                result->void_nodes[best_i],
                result->void_nodes[best_j]);

            {
                /* Update bounding box to encompass both */
                result->void_boxes[best_i].min_x = fmin(
                    result->void_boxes[best_i].min_x,
                    result->void_boxes[best_j].min_x);
                result->void_boxes[best_i].max_x = fmax(
                    result->void_boxes[best_i].max_x,
                    result->void_boxes[best_j].max_x);
                result->void_boxes[best_i].min_y = fmin(
                    result->void_boxes[best_i].min_y,
                    result->void_boxes[best_j].min_y);
                result->void_boxes[best_i].max_y = fmax(
                    result->void_boxes[best_i].max_y,
                    result->void_boxes[best_j].max_y);
                result->void_boxes[best_i].min_z = fmin(
                    result->void_boxes[best_i].min_z,
                    result->void_boxes[best_j].min_z);
                result->void_boxes[best_i].max_z = fmax(
                    result->void_boxes[best_i].max_z,
                    result->void_boxes[best_j].max_z);
                surfaces[best_i] = surfaces[best_i] + surfaces[best_j];
            }

            /* Update cell i with merged result */
            result->void_nodes[best_i] = merged_node;

            /* Mark cell j as dead */
            alea_bitset_clear(&alive, best_j);
            alive_count--;
        }
    }

    /* Compact the arrays to remove dead cells */
    size_t write_idx = 0;
    for (size_t read_idx = 0; read_idx < n; read_idx++) {
        if (alea_bitset_test(&alive, read_idx)) {
            if (write_idx != read_idx) {
                result->void_nodes[write_idx] = result->void_nodes[read_idx];
                result->void_boxes[write_idx] = result->void_boxes[read_idx];
            }
            write_idx++;
        }
    }
    result->void_node_count = write_idx;
    result->void_box_count = write_idx;

    ALEA_LOG_INFO("Void merge: reduced to %zu cells (%.1f%% reduction)",
           result->void_node_count,
           100.0 * (1.0 - (double)result->void_node_count / n));

    alea_bitset_destroy(&alive);
    free(surfaces);

    return (int)result->void_node_count;
}