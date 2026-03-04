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
#include "primitives/bbox.h"

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
    const alea_bbox_t* box
) {
    /* Create 6 plane surfaces using public API (auto-assigned IDs) */
    alea_node_id_t planes[6];

    /* x > min_x: positive sense of plane x - min_x = 0 */
    int px0_idx = alea_plane_surface(sys, 0, 1, 0, 0, -box->min_x);
    /* x < max_x: negative sense of plane x - max_x = 0 */
    int px1_idx = alea_plane_surface(sys, 0, 1, 0, 0, -box->max_x);
    /* y > min_y */
    int py0_idx = alea_plane_surface(sys, 0, 0, 1, 0, -box->min_y);
    /* y < max_y */
    int py1_idx = alea_plane_surface(sys, 0, 0, 1, 0, -box->max_y);
    /* z > min_z */
    int pz0_idx = alea_plane_surface(sys, 0, 0, 0, 1, -box->min_z);
    /* z < max_z */
    int pz1_idx = alea_plane_surface(sys, 0, 0, 0, 1, -box->max_z);

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
 * Create a regional void cell by intersecting a void complement with a bounding box.
 * regional_void = void_complement ∩ box
 *
 * Uses 6 planes (not BOX macrobody) for proper MCNP surface export.
 */
static alea_node_id_t create_regional_void(
    alea_system_t* sys,
    alea_node_id_t void_complement,
    const alea_bbox_t* box
) {
    /* Create box region from 6 planes (registers surfaces automatically) */
    alea_node_id_t box_node = create_box_from_planes(sys, box);

    if (box_node == ALEA_NODE_ID_INVALID) {
        return ALEA_NODE_ID_INVALID;
    }

    /* Regional void = void_complement ∩ box */
    alea_node_id_t regional_void = alea_create_intersection(sys, void_complement, box_node);
    return regional_void;
}

// ============================================================================
// PER-REGION COMPLEMENT: cell bbox precomputation and local filtering
// ============================================================================

/**
 * @brief Entry for precomputed cell bounding box
 */
typedef struct {
    alea_node_id_t root;
    alea_bbox_t bbox;
} cell_bbox_entry_t;

/**
 * @brief Collect bboxes of all material cells in a universe
 *
 * @param sys         CSG system
 * @param universe_id Universe to scan
 * @param out_entries Output array (caller frees)
 * @param out_count   Number of entries
 * @return 0 on success, -1 on error
 */
static int collect_cell_bboxes(alea_system_t* sys, int universe_id,
                               cell_bbox_entry_t** out_entries, size_t* out_count) {
    const alea_universe_t* univ = alea_get_universe(sys, universe_id);
    if (!univ || univ->cell_count == 0) {
        *out_entries = NULL;
        *out_count = 0;
        return 0;
    }

    cell_bbox_entry_t* entries = malloc(univ->cell_count * sizeof(cell_bbox_entry_t));
    if (!entries) return -1;

    size_t count = 0;
    for (size_t i = 0; i < univ->cell_count; i++) {
        size_t cell_idx = univ->cell_indices[i];
        const alea_cell_entry_t* cell = &sys->cells.data[cell_idx];

        /* Skip void cells and invalid cells */
        if (cell->material_id == 0) continue;
        if (cell->root_node_id == ALEA_NODE_ID_INVALID) continue;

        entries[count].root = cell->root_node_id;
        entries[count].bbox = alea_get_bbox(sys, cell->root_node_id);
        count++;
    }

    *out_entries = entries;
    *out_count = count;
    return 0;
}

/**
 * @brief Build a local void complement from cells overlapping a box
 *
 * Returns ¬(union of overlapping cells), or ALEA_NODE_ID_INVALID if no cells
 * overlap (meaning the entire box is void — caller should use box_node directly).
 */
static alea_node_id_t compute_local_void(
    alea_system_t* sys,
    const alea_bbox_t* box,
    const cell_bbox_entry_t* cell_bboxes,
    size_t cell_bbox_count
) {
    /* Collect roots of cells whose bbox overlaps this box */
    alea_node_id_t* overlapping = malloc(cell_bbox_count * sizeof(alea_node_id_t));
    if (!overlapping) return ALEA_NODE_ID_INVALID;

    size_t overlap_count = 0;
    for (size_t j = 0; j < cell_bbox_count; j++) {
        if (alea_bbox_intersects(box, &cell_bboxes[j].bbox)) {
            overlapping[overlap_count++] = cell_bboxes[j].root;
        }
    }

    if (overlap_count == 0) {
        free(overlapping);
        return ALEA_NODE_ID_INVALID; /* signal: pure void, no complement needed */
    }

    /* Build local union and complement */
    alea_node_id_t local_union;
    if (overlap_count == 1) {
        local_union = overlapping[0];
    } else {
        local_union = alea_create_union_many(sys, overlapping, overlap_count);
    }
    free(overlapping);

    if (local_union == ALEA_NODE_ID_INVALID) {
        return ALEA_NODE_ID_INVALID;
    }

    return alea_create_complement(sys, local_union);
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
 * Create void CSG nodes from candidate boxes using per-region complement filtering.
 *
 * For each box:
 * 1. Find material cells whose bbox overlaps this box
 * 2. Build local_void = ¬(union of overlapping cells)
 * 3. Create regional_void = local_void ∩ box (using 6 planes)
 * 4. If result is INVALID (empty), the box had no actual void - skip it
 * 5. If result is valid, store it as a void region (CSG tree + bbox)
 *
 * Each void region only references surfaces of nearby cells, not all cells
 * in the system. This naturally bounds per-region complexity.
 */
static int create_void_nodes_filtered(
    alea_system_t* sys,
    const cell_bbox_entry_t* cell_bboxes,
    size_t cell_bbox_count,
    alea_bbox_t* boxes,
    size_t box_count,
    void_result_t* result
) {
    size_t surfaces_before = alea_vec_count(&sys->surfaces);

    ALEA_LOG_INFO("create_void_nodes_filtered: %zu cells, %zu boxes",
           cell_bbox_count, box_count);

    for (size_t i = 0; i < box_count; i++) {
        /* Build per-region complement from overlapping cells only */
        alea_node_id_t local_void = compute_local_void(
            sys, &boxes[i], cell_bboxes, cell_bbox_count
        );

        alea_node_id_t regional_void;
        if (local_void == ALEA_NODE_ID_INVALID) {
            /* No cells overlap this box — entire box is void.
             * The region is just the box itself (no complement needed). */
            alea_node_id_t box_node = create_box_from_planes(
                sys, &boxes[i]
            );
            regional_void = box_node;
        } else {
            regional_void = create_regional_void(
                sys, local_void, &boxes[i]
            );
        }

        ALEA_LOG_DEBUG("  Box %zu: [%.1f,%.1f]x[%.1f,%.1f]x[%.1f,%.1f] regional_void=%u",
               i, boxes[i].min_x, boxes[i].max_x,
               boxes[i].min_y, boxes[i].max_y,
               boxes[i].min_z, boxes[i].max_z, regional_void);

        /* If regional void is empty (no actual void in this box), skip it. */
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

    result->surfaces_created = alea_vec_count(&sys->surfaces) - surfaces_before;
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

    /* Step 2: Collect candidate boxes from octree with bottom-up collapse.
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

    /* Step 3: Precompute cell bboxes for per-region filtering */
    cell_bbox_entry_t* cell_bboxes = NULL;
    size_t cell_bbox_count = 0;
    if (collect_cell_bboxes(sys, 0, &cell_bboxes, &cell_bbox_count) < 0) {
        ALEA_LOG_WARN("cell bbox collection failed");
        free(candidate_boxes);
        alea_void_print_stats(result);
        return result;
    }

    ALEA_LOG_INFO("Material cells for filtering: %zu", cell_bbox_count);

    /* Step 4: Create void regions using per-region complement filtering.
     * Each box gets: local_void = ¬(overlapping cells) ∩ box */

    if (create_void_nodes_filtered(sys, cell_bboxes, cell_bbox_count,
                                   candidate_boxes, box_count,
                                   result) < 0) {
        ALEA_LOG_WARN("void cell creation incomplete");
    }

    free(cell_bboxes);
    free(candidate_boxes);

    /* Compute global_void lazily — only needed for consolidation */
    result->global_void = ALEA_NODE_ID_INVALID;

    ALEA_LOG_INFO("Void CSG processing (per-region filtering):");
    ALEA_LOG_INFO("  Candidate boxes checked: %zu", box_count);
    ALEA_LOG_INFO("  Material cells for filtering: %zu", cell_bbox_count);
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

int alea_void_add_graveyard(alea_system_t* sys, void_result_t* result) {
    if (!sys || !result || !result->root) return -1;

    /* Get outer bounds from the void octree root */
    const alea_bbox_t* bb = &result->root->bbox;

    /* Compute enclosing sphere: center = bbox center, radius = 0.51 * diagonal */
    double cx = (bb->min_x + bb->max_x) * 0.5;
    double cy = (bb->min_y + bb->max_y) * 0.5;
    double cz = (bb->min_z + bb->max_z) * 0.5;

    double dx = bb->max_x - bb->min_x;
    double dy = bb->max_y - bb->min_y;
    double dz = bb->max_z - bb->min_z;
    double diag = sqrt(dx * dx + dy * dy + dz * dz);
    double radius = 0.51 * diag;

    /* Create sphere surface (surface_id=0 for auto-assign) */
    int sphere_idx = alea_sphere_surface(sys, 0, cx, cy, cz, radius);
    if (sphere_idx < 0) return -1;

    /* Create RPP surface for the bounding box */
    int box_idx = alea_box_surface(sys, 0,
                                   bb->min_x, bb->max_x,
                                   bb->min_y, bb->max_y,
                                   bb->min_z, bb->max_z);
    if (box_idx < 0) return -1;

    /* Mark sphere surface as vacuum boundary */
    sys->surfaces.data[sphere_idx].boundary_type = ALEA_BOUNDARY_VACUUM;

    alea_node_id_t inside_sphere = alea_halfspace(sys, sphere_idx, -1);
    alea_node_id_t outside_box   = alea_halfspace(sys, box_idx, +1);
    alea_node_id_t outside_sphere = alea_halfspace(sys, sphere_idx, +1);
    if (inside_sphere == ALEA_NODE_ID_INVALID ||
        outside_box   == ALEA_NODE_ID_INVALID ||
        outside_sphere == ALEA_NODE_ID_INVALID) return -1;

    /* Shell cell: inside sphere but outside bounding box (void, IMP:N=1) */
    alea_node_id_t shell = alea_intersection(sys, inside_sphere, outside_box);
    if (shell == ALEA_NODE_ID_INVALID) return -1;

    int cell_id = alea_max_cell_id(sys) + 1;
    int shell_idx = alea_add_cell(sys, cell_id, shell, ALEA_MATERIAL_VOID, 0.0, 0);
    if (shell_idx < 0) return -1;

    ALEA_LOG_INFO("Shell cell %d: void between bbox and sphere, IMP:N=1", cell_id);

    /* Graveyard cell: outside the sphere (void)
     * Graveyard is identified by vacuum boundary on sphere surface (line 835).
     * MCNP module's on_cell_added callback sets IMP:N=0 on graveyard params. */
    cell_id = alea_max_cell_id(sys) + 1;
    int grav_idx = alea_add_cell(sys, cell_id, outside_sphere, ALEA_MATERIAL_VOID, 0.0, 0);
    if (grav_idx < 0) return -1;

    ALEA_LOG_INFO("Graveyard cell %d: sphere at (%.2f, %.2f, %.2f) R=%.2f",
                  cell_id, cx, cy, cz, radius);

    return 1;
}

/* ============================================================================
 * VOID CELL MERGING
 * ============================================================================ */

const alea_void_merge_config_t ALEA_VOID_MERGE_DEFAULT = {
    .cell_weight = 1.0,
    .surface_weight = 0.1,
    .max_surfaces_per_cell = 24,
    .min_cells = 1,
    .consolidate_max_surfaces = 100
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

/* Compact alive regions and log reduction */
static int compact_and_finish(void_result_t* result, alea_bitset_t* alive,
                              int* surfaces_arr, size_t original_n) {
    size_t write_idx = 0;
    for (size_t read_idx = 0; read_idx < original_n; read_idx++) {
        if (alea_bitset_test(alive, read_idx)) {
            if (write_idx != read_idx) {
                result->void_regions[write_idx] = result->void_regions[read_idx];
            }
            write_idx++;
        }
    }
    result->void_region_count = write_idx;

    ALEA_LOG_INFO("Void merge: reduced to %zu cells (%.1f%% reduction)",
           result->void_region_count,
           100.0 * (1.0 - (double)result->void_region_count / original_n));

    alea_bitset_destroy(alive);
    free(surfaces_arr);

    return (int)result->void_region_count;
}

/* ---------------------------------------------------------------------------
 * Greedy O(n^3) merge — original algorithm
 * ------------------------------------------------------------------------- */

static int merge_void_cells_greedy(alea_system_t* sys,
                                   void_result_t* result,
                                   const alea_void_merge_config_t* cfg) {
    size_t n = result->void_region_count;

    alea_bitset_t alive = alea_bitset_create(n);
    int* surfaces = calloc(n, sizeof(int));

    if (!alive.words || !surfaces) {
        alea_bitset_destroy(&alive);
        free(surfaces);
        return -1;
    }

    size_t alive_count = n;
    for (size_t i = 0; i < n; i++) {
        alea_bitset_set(&alive, i);
        surfaces[i] = count_surfaces(sys, result->void_regions[i].node);
    }

    ALEA_LOG_INFO("Void merge (greedy): starting with %zu cells", n);

    bool improved = true;
    while (improved && alive_count > (size_t)cfg->min_cells) {
        improved = false;

        double best_delta = 0;
        size_t best_i = 0, best_j = 0;
        bool best_is_box = false;

        for (size_t i = 0; i < n; i++) {
            if (!alea_bitset_test(&alive, i)) continue;
            if (g_alea_interrupted) { alea_bitset_destroy(&alive); free(surfaces); ALEA_CHECK_INTERRUPTED(-1); }

            for (size_t j = i + 1; j < n; j++) {
                if (!alea_bitset_test(&alive, j)) continue;

                if (!boxes_adjacent(&result->void_regions[i].bbox,
                                    &result->void_regions[j].bbox)) {
                    continue;
                }

                alea_bbox_t merged_bbox;
                bool is_box = boxes_can_form_box(&result->void_regions[i].bbox,
                                                  &result->void_regions[j].bbox,
                                                  &merged_bbox);
                int new_surfaces = is_box ? 6 : surfaces[i] + surfaces[j];

                if (new_surfaces > cfg->max_surfaces_per_cell) continue;

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

        if (best_delta < -1e-9) {
            improved = true;

            ALEA_LOG_DEBUG("  Merging cells %zu + %zu (delta=%.2f)",
                   best_i, best_j, best_delta);

            alea_node_id_t merged_node = alea_create_union(sys,
                result->void_regions[best_i].node,
                result->void_regions[best_j].node);

            alea_bbox_t* bi = &result->void_regions[best_i].bbox;
            const alea_bbox_t* bj = &result->void_regions[best_j].bbox;
            bi->min_x = fmin(bi->min_x, bj->min_x);
            bi->max_x = fmax(bi->max_x, bj->max_x);
            bi->min_y = fmin(bi->min_y, bj->min_y);
            bi->max_y = fmax(bi->max_y, bj->max_y);
            bi->min_z = fmin(bi->min_z, bj->min_z);
            bi->max_z = fmax(bi->max_z, bj->max_z);

            surfaces[best_i] = best_is_box ? 6 : surfaces[best_i] + surfaces[best_j];
            result->void_regions[best_i].node = merged_node;

            alea_bitset_clear(&alive, best_j);
            alive_count--;
        }
    }

    return compact_and_finish(result, &alive, surfaces, n);
}

/* ---------------------------------------------------------------------------
 * Face-sorted O(n log^2 n) merge
 * ------------------------------------------------------------------------- */

typedef struct {
    uint8_t  axis;        /* 0=X, 1=Y, 2=Z */
    uint8_t  is_max_face; /* 0=min face, 1=max face */
    uint32_t box_index;
    double   coord;       /* face coordinate on this axis */
    double   lo1, hi1;    /* extent on first perpendicular axis */
    double   lo2, hi2;    /* extent on second perpendicular axis */
} face_entry_t;

typedef struct {
    uint32_t i, j;
    bool     is_box;
    alea_bbox_t merged_bbox;
} merge_candidate_t;

/* Emit 6 face entries for one box */
static void emit_faces(const alea_bbox_t* bbox, uint32_t box_index,
                       face_entry_t* faces, size_t* count) {
    /* X-min face (is_max=0): perp axes are Y,Z */
    faces[(*count)++] = (face_entry_t){0, 0, box_index, bbox->min_x,
                                       bbox->min_y, bbox->max_y,
                                       bbox->min_z, bbox->max_z};
    /* X-max face (is_max=1) */
    faces[(*count)++] = (face_entry_t){0, 1, box_index, bbox->max_x,
                                       bbox->min_y, bbox->max_y,
                                       bbox->min_z, bbox->max_z};
    /* Y-min */
    faces[(*count)++] = (face_entry_t){1, 0, box_index, bbox->min_y,
                                       bbox->min_x, bbox->max_x,
                                       bbox->min_z, bbox->max_z};
    /* Y-max */
    faces[(*count)++] = (face_entry_t){1, 1, box_index, bbox->max_y,
                                       bbox->min_x, bbox->max_x,
                                       bbox->min_z, bbox->max_z};
    /* Z-min */
    faces[(*count)++] = (face_entry_t){2, 0, box_index, bbox->min_z,
                                       bbox->min_x, bbox->max_x,
                                       bbox->min_y, bbox->max_y};
    /* Z-max */
    faces[(*count)++] = (face_entry_t){2, 1, box_index, bbox->max_z,
                                       bbox->min_x, bbox->max_x,
                                       bbox->min_y, bbox->max_y};
}

static int face_entry_cmp(const void* va, const void* vb) {
    const face_entry_t* a = (const face_entry_t*)va;
    const face_entry_t* b = (const face_entry_t*)vb;
    if (a->axis != b->axis) return (int)a->axis - (int)b->axis;
    if (a->coord < b->coord - MERGE_TOL) return -1;
    if (a->coord > b->coord + MERGE_TOL) return  1;
    return 0;
}

/* Check if perpendicular extents overlap (shared-face test) */
static bool faces_overlap_2d(const face_entry_t* a, const face_entry_t* b) {
    return (a->lo1 < b->hi1 - MERGE_TOL) && (b->lo1 < a->hi1 - MERGE_TOL) &&
           (a->lo2 < b->hi2 - MERGE_TOL) && (b->lo2 < a->hi2 - MERGE_TOL);
}

/* Check if perpendicular extents are exactly aligned (box merge possible) */
static bool faces_aligned_2d(const face_entry_t* a, const face_entry_t* b) {
    return fabs(a->lo1 - b->lo1) < MERGE_TOL && fabs(a->hi1 - b->hi1) < MERGE_TOL &&
           fabs(a->lo2 - b->lo2) < MERGE_TOL && fabs(a->hi2 - b->hi2) < MERGE_TOL;
}

static int merge_void_cells_face_sorted(alea_system_t* sys,
                                        void_result_t* result,
                                        const alea_void_merge_config_t* cfg) {
    size_t n = result->void_region_count;

    alea_bitset_t alive = alea_bitset_create(n);
    int* surfaces = calloc(n, sizeof(int));

    if (!alive.words || !surfaces) {
        alea_bitset_destroy(&alive);
        free(surfaces);
        return -1;
    }

    size_t alive_count = n;
    for (size_t i = 0; i < n; i++) {
        alea_bitset_set(&alive, i);
        surfaces[i] = count_surfaces(sys, result->void_regions[i].node);
    }

    ALEA_LOG_INFO("Void merge (face-sorted): starting with %zu cells", n);

    /* Allocate face array and candidate array (reused each pass) */
    size_t face_cap = 6 * n;
    face_entry_t* faces = malloc(face_cap * sizeof(face_entry_t));
    size_t cand_cap = n;
    merge_candidate_t* candidates = malloc(cand_cap * sizeof(merge_candidate_t));

    if (!faces || !candidates) {
        free(faces);
        free(candidates);
        alea_bitset_destroy(&alive);
        free(surfaces);
        return -1;
    }

    /* Per-pass bitset: marks boxes already paired this pass */
    alea_bitset_t paired = alea_bitset_create(n);
    if (!paired.words) {
        free(faces);
        free(candidates);
        alea_bitset_destroy(&alive);
        free(surfaces);
        return -1;
    }

    int pass = 0;
    bool merged_any = true;

    while (merged_any && alive_count > (size_t)cfg->min_cells) {
        merged_any = false;
        pass++;

        if (g_alea_interrupted) {
            free(faces); free(candidates);
            alea_bitset_destroy(&paired);
            alea_bitset_destroy(&alive);
            free(surfaces);
            ALEA_CHECK_INTERRUPTED(-1);
        }

        /* Step 1: Emit faces for alive boxes */
        size_t face_count = 0;
        for (size_t i = 0; i < n; i++) {
            if (!alea_bitset_test(&alive, i)) continue;
            /* Grow face array if needed */
            if (face_count + 6 > face_cap) {
                face_cap *= 2;
                face_entry_t* tmp = realloc(faces, face_cap * sizeof(face_entry_t));
                if (!tmp) goto cleanup;
                faces = tmp;
            }
            emit_faces(&result->void_regions[i].bbox, (uint32_t)i,
                       faces, &face_count);
        }

        /* Step 2: Sort faces by (axis, coord) */
        qsort(faces, face_count, sizeof(face_entry_t), face_entry_cmp);

        /* Step 3: Sweep — find merge candidates */
        size_t cand_count = 0;
        alea_bitset_clear_all(&paired);

        size_t grp_start = 0;
        while (grp_start < face_count) {
            /* Find group end: same axis and coord within tolerance */
            size_t grp_end = grp_start + 1;
            while (grp_end < face_count &&
                   faces[grp_end].axis == faces[grp_start].axis &&
                   fabs(faces[grp_end].coord - faces[grp_start].coord) < MERGE_TOL) {
                grp_end++;
            }

            /* Sub-pass A: box merges (max-face paired with min-face, aligned) */
            for (size_t mi = grp_start; mi < grp_end; mi++) {
                if (faces[mi].is_max_face != 1) continue;
                uint32_t bi = faces[mi].box_index;
                if (alea_bitset_test(&paired, bi)) continue;
                if (!alea_bitset_test(&alive, bi)) continue;

                for (size_t mj = grp_start; mj < grp_end; mj++) {
                    if (faces[mj].is_max_face != 0) continue;
                    uint32_t bj = faces[mj].box_index;
                    if (bi == bj) continue;
                    if (alea_bitset_test(&paired, bj)) continue;
                    if (!alea_bitset_test(&alive, bj)) continue;

                    if (!faces_aligned_2d(&faces[mi], &faces[mj])) continue;

                    /* This is a box merge candidate */
                    alea_bbox_t merged_bbox;
                    if (!boxes_can_form_box(&result->void_regions[bi].bbox,
                                           &result->void_regions[bj].bbox,
                                           &merged_bbox)) {
                        continue;
                    }

                    /* Grow candidates array if needed */
                    if (cand_count >= cand_cap) {
                        cand_cap *= 2;
                        merge_candidate_t* tmp = realloc(candidates,
                            cand_cap * sizeof(merge_candidate_t));
                        if (!tmp) goto cleanup;
                        candidates = tmp;
                    }

                    candidates[cand_count++] = (merge_candidate_t){
                        .i = bi, .j = bj, .is_box = true,
                        .merged_bbox = merged_bbox
                    };
                    alea_bitset_set(&paired, bi);
                    alea_bitset_set(&paired, bj);
                    break; /* each max-face pairs with at most one min-face */
                }
            }

            /* Sub-pass B: union merges (remaining unpaired faces) */
            for (size_t mi = grp_start; mi < grp_end; mi++) {
                if (faces[mi].is_max_face != 1) continue;
                uint32_t bi = faces[mi].box_index;
                if (alea_bitset_test(&paired, bi)) continue;
                if (!alea_bitset_test(&alive, bi)) continue;

                for (size_t mj = grp_start; mj < grp_end; mj++) {
                    if (faces[mj].is_max_face != 0) continue;
                    uint32_t bj = faces[mj].box_index;
                    if (bi == bj) continue;
                    if (alea_bitset_test(&paired, bj)) continue;
                    if (!alea_bitset_test(&alive, bj)) continue;

                    if (!faces_overlap_2d(&faces[mi], &faces[mj])) continue;

                    int new_surfaces = surfaces[bi] + surfaces[bj];
                    if (new_surfaces > cfg->max_surfaces_per_cell) continue;

                    /* Cost check: union merges save one cell but add surfaces */
                    double before_cost = 2 * cfg->cell_weight +
                                        (surfaces[bi] + surfaces[bj]) * cfg->surface_weight;
                    double after_cost = 1 * cfg->cell_weight +
                                       new_surfaces * cfg->surface_weight;
                    if (after_cost >= before_cost - 1e-9) continue;

                    /* Grow candidates array if needed */
                    if (cand_count >= cand_cap) {
                        cand_cap *= 2;
                        merge_candidate_t* tmp = realloc(candidates,
                            cand_cap * sizeof(merge_candidate_t));
                        if (!tmp) goto cleanup;
                        candidates = tmp;
                    }

                    alea_bbox_t union_bbox = result->void_regions[bi].bbox;
                    const alea_bbox_t* bbj = &result->void_regions[bj].bbox;
                    union_bbox.min_x = fmin(union_bbox.min_x, bbj->min_x);
                    union_bbox.max_x = fmax(union_bbox.max_x, bbj->max_x);
                    union_bbox.min_y = fmin(union_bbox.min_y, bbj->min_y);
                    union_bbox.max_y = fmax(union_bbox.max_y, bbj->max_y);
                    union_bbox.min_z = fmin(union_bbox.min_z, bbj->min_z);
                    union_bbox.max_z = fmax(union_bbox.max_z, bbj->max_z);

                    candidates[cand_count++] = (merge_candidate_t){
                        .i = bi, .j = bj, .is_box = false,
                        .merged_bbox = union_bbox
                    };
                    alea_bitset_set(&paired, bi);
                    alea_bitset_set(&paired, bj);
                    break;
                }
            }

            grp_start = grp_end;
        }

        /* Step 4: Execute all collected merges */
        for (size_t c = 0; c < cand_count; c++) {
            uint32_t ci = candidates[c].i;
            uint32_t cj = candidates[c].j;

            /* Double-check both still alive (should be, but safety) */
            if (!alea_bitset_test(&alive, ci) || !alea_bitset_test(&alive, cj))
                continue;

            ALEA_LOG_DEBUG("  Pass %d: merging %u + %u (%s)", pass, ci, cj,
                   candidates[c].is_box ? "box" : "union");

            alea_node_id_t merged_node = alea_create_union(sys,
                result->void_regions[ci].node,
                result->void_regions[cj].node);

            result->void_regions[ci].node = merged_node;
            result->void_regions[ci].bbox = candidates[c].merged_bbox;
            surfaces[ci] = candidates[c].is_box ? 6 : surfaces[ci] + surfaces[cj];

            alea_bitset_clear(&alive, cj);
            alive_count--;
            merged_any = true;
        }

        ALEA_LOG_DEBUG("  Pass %d: %zu merges, %zu alive", pass, cand_count, alive_count);
    }

    free(faces);
    free(candidates);
    alea_bitset_destroy(&paired);
    return compact_and_finish(result, &alive, surfaces, n);

cleanup:
    free(faces);
    free(candidates);
    alea_bitset_destroy(&paired);
    alea_bitset_destroy(&alive);
    free(surfaces);
    return -1;
}

/* ---------------------------------------------------------------------------
 * Void consolidation: replace N regions with single global_void ∩ outer_bbox
 * ------------------------------------------------------------------------- */

static void consolidate_void_regions(alea_system_t* sys, void_result_t* result,
                                     int max_surfaces) {
    if (result->void_region_count <= 1) return;

    /* 1. Compute combined bbox = union of all remaining void region bboxes */
    alea_bbox_t combined = result->void_regions[0].bbox;
    for (size_t i = 1; i < result->void_region_count; i++) {
        combined = alea_bbox_union(&combined, &result->void_regions[i].bbox);
    }

    /* 2. Compute global_void lazily (only needed for consolidation) */
    if (result->global_void == ALEA_NODE_ID_INVALID) {
        result->global_void = compute_global_void(sys, 0);
        if (result->global_void == ALEA_NODE_ID_INVALID) {
            ALEA_LOG_INFO("void consolidation: cannot compute global void");
            return;
        }
    }

    /* 3. Create single region: global_void ∩ combined_bbox */
    alea_node_id_t node = create_regional_void(sys, result->global_void,
                                               &combined);
    if (node == ALEA_NODE_ID_INVALID) {
        ALEA_LOG_INFO("void consolidation: node invalid, keeping %zu regions",
                      result->void_region_count);
        return;
    }

    /* 5. Check surface count — reject if too complex */
    int nsurfaces = count_surfaces(sys, node);
    if (nsurfaces > max_surfaces) {
        ALEA_LOG_INFO("void consolidation: %d surfaces > limit %d, keeping %zu regions",
                      nsurfaces, max_surfaces, result->void_region_count);
        return;
    }

    /* 6. Accept: replace all regions with this single one */
    ALEA_LOG_INFO("void consolidation: %zu regions -> 1 (%d surfaces, "
                  "bbox [%.1f,%.1f]x[%.1f,%.1f]x[%.1f,%.1f])",
                  result->void_region_count, nsurfaces,
                  combined.min_x, combined.max_x,
                  combined.min_y, combined.max_y,
                  combined.min_z, combined.max_z);
    result->void_regions[0].node = node;
    result->void_regions[0].bbox = combined;
    result->void_region_count = 1;
}

/* ---------------------------------------------------------------------------
 * Dispatcher
 * ------------------------------------------------------------------------- */

int alea_merge_void_cells(alea_system_t* sys,
                         void_result_t* result,
                         const alea_void_merge_config_t* config) {
    if (!sys || !result) return -1;
    if (result->void_region_count == 0) return 0;

    const alea_void_merge_config_t* cfg = config ? config : &ALEA_VOID_MERGE_DEFAULT;

    int ret;
    if (cfg->use_greedy) {
        ret = merge_void_cells_greedy(sys, result, cfg);
    } else {
        ret = merge_void_cells_face_sorted(sys, result, cfg);
    }

    if (ret >= 0 && cfg->consolidate_max_surfaces > 0 &&
        result->void_region_count > 1) {
        consolidate_void_regions(sys, result, cfg->consolidate_max_surfaces);
        ret = (int)result->void_region_count;
    }

    return ret;
}