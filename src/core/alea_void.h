// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_VOID_H
#define ALEA_VOID_H

#include "alea_types.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file alea_void.h
 * @brief Octree-based void region generation
 *
 * Uses adaptive octree subdivision to efficiently identify void regions.
 *
 * Key design principle: CONSERVATIVE classification.
 * - Octree probing is an OPTIMIZATION to skip definitely-solid regions
 * - It is NOT the source of truth for void identification
 * - Only boxes where ALL probes are in the SAME material cell are skipped
 * - Any box with ANY void probe, or mixed cells, is kept as a candidate
 * - CSG math (complement intersection + simplification) decides what's void
 *
 * This approach may process more boxes than strictly necessary, but it
 * guarantees we never miss void regions due to probing artifacts.
 */

/* Forward declaration */
typedef struct alea_system alea_system_t;

// ============================================================================
// OCTREE NODE CLASSIFICATION
// ============================================================================

typedef enum {
    OCTREE_UNKNOWN = 0,     // Not yet classified
    OCTREE_VOID,            // Entirely void (no cells)
    OCTREE_SOLID,           // Entirely inside one cell
    OCTREE_MIXED,           // Contains boundary, subdivided
    OCTREE_PARTIAL_VOID     // At max depth, contains some void
} octree_class_t;

// ============================================================================
// OCTREE NODE
// ============================================================================

typedef struct octree_node {
    alea_bbox_t bbox;
    octree_class_t classification;
    int cell_index;                 // If SOLID, which cell (-1 for void/mixed)
    int depth;
    struct octree_node* children[8]; // NULL if leaf
} octree_node_t;

// ============================================================================
// OCTREE CONFIGURATION
// ============================================================================

typedef struct octree_config {
    int max_depth;          // Maximum subdivision depth (default: 8)
    int probes_per_axis;    // Probe points per axis for classification (default: 3)
    double min_size;        // Minimum node size, stop subdividing (default: 0.1)
} octree_config_t;

#define OCTREE_DEFAULT_CONFIG ((octree_config_t){ \
    .max_depth = 8, \
    .probes_per_axis = 3, \
    .min_size = 0.1 \
})

// ============================================================================
// VOID REGION
// ============================================================================

/**
 * @brief A void region: CSG node + bounding box
 */
typedef struct {
    alea_node_id_t node;
    alea_bbox_t bbox;
} void_region_t;

// ============================================================================
// OCTREE RESULT
// ============================================================================

/**
 * @brief Result of void generation
 */
typedef struct void_result {
    octree_node_t* root;    // Octree root (caller must free with octree_destroy)

    // Global void complement: ¬(Cell1 ∪ Cell2 ∪ ... ∪ CellN)
    // This is the exact analytical definition of void
    alea_node_id_t global_void;

    // Regional void regions: each is (global_void ∩ region_box).
    // These are CSG tree node IDs, NOT cells. Use alea_void_add_cells() to
    // register them as actual cells in the system.
    void_region_t* void_regions;
    size_t void_region_count;
    size_t void_region_capacity;

    // Statistics
    size_t total_nodes;           // Total octree nodes created
    size_t octree_void_count;     // Octree nodes classified as definitely void
    size_t solid_nodes;           // Nodes classified as definitely solid (skipped)
    size_t partial_nodes;         // Nodes that might contain void (candidates)
    size_t max_depth_reached;     // Times we hit max depth limit
    size_t empty_regions_skipped; // Candidate boxes that CSG determined were solid
    size_t boxes_before_merge;    // Candidate boxes from octree
    size_t boxes_after_merge;     // After any merging
    size_t surfaces_created;      // Plane surfaces created for void boxes
} void_result_t;

// ============================================================================
// API
// ============================================================================

/**
 * @brief Generate void regions using octree
 * 
 * @param sys CSG system (must have cells converted)
 * @param bounds Bounding box to analyze (NULL = auto from universe 0)
 * @param config Configuration (NULL = use defaults)
 * @return Result structure, or NULL on error
 */
void_result_t* alea_generate_void_octree(alea_system_t* sys,
                                         const alea_bbox_t* bounds,
                                         const octree_config_t* config);

/**
 * @brief Create CSG node from void result
 * 
 * Creates union of all void boxes. Use this if you need a CSG tree
 * for boolean operations.
 * 
 * @param sys CSG system
 * @param result Void generation result
 * @return Node ID of void union, or ALEA_NODE_ID_INVALID on error
 */
alea_node_id_t alea_void_to_node(alea_system_t* sys, const void_result_t* result);

/**
 * @brief Add void regions as cells to the system
 *
 * Registers each void node from the result as an actual cell in the system.
 * Cells are numbered sequentially starting from max_cell_id + 1.
 *
 * @param sys CSG system
 * @param result Void generation result
 * @return Number of cells added, or -1 on error
 */
int alea_void_add_cells(alea_system_t* sys, void_result_t* result);

/**
 * @brief Configuration for void cell merging
 */
typedef struct {
    double cell_weight;           /**< Cost per cell (default: 1.0) */
    double surface_weight;        /**< Cost per surface (default: 0.1) */
    int max_surfaces_per_cell;    /**< Hard limit on surfaces per cell (default: 24) */
    int min_cells;                /**< Stop merging below this count (default: 1) */
} alea_void_merge_config_t;

/**
 * @brief Default merge configuration
 */
extern const alea_void_merge_config_t ALEA_VOID_MERGE_DEFAULT;

/**
 * @brief Merge void cells to reduce count while balancing complexity
 *
 * Uses greedy algorithm to merge adjacent void regions when beneficial.
 * Cost function: cost = cell_weight * num_cells + surface_weight * total_surfaces
 *
 * Merging strategy:
 * - Adjacent aligned boxes → larger box (no complexity increase)
 * - Adjacent non-aligned boxes → union (complexity increases)
 * - Stop when no merge reduces cost or limits reached
 *
 * @param sys CSG system
 * @param result Void result to modify in place
 * @param config Merge configuration (NULL for defaults)
 * @return Number of cells after merging, or -1 on error
 */
int alea_merge_void_cells(alea_system_t* sys,
                         void_result_t* result,
                         const alea_void_merge_config_t* config);

/**
 * @brief Free void result
 */
void alea_void_result_destroy(void_result_t* result);

/**
 * @brief Get number of void boxes
 */
size_t alea_void_box_count(const void_result_t* result);

/**
 * @brief Get void box by index
 */
int alea_void_box_get(const void_result_t* result, size_t index, alea_bbox_t* box);

/**
 * @brief Free octree node recursively
 */
void octree_destroy(octree_node_t* node);

/**
 * @brief Print octree statistics
 */
void alea_void_print_stats(const void_result_t* result);

#ifdef __cplusplus
}
#endif

#endif /* ALEA_VOID_H */