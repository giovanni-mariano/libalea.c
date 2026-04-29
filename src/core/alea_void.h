// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_VOID_H
#define ALEA_VOID_H

#include "alea_types.h"
#include "util/alea_vec.h"
#include <stddef.h>


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

ALEA_VEC_DEFINE(void_region_vec, void_region_t);

// ============================================================================
// OCTREE RESULT
// ============================================================================

/**
 * @brief Result of void generation
 *
 * Generation appends surfaces, primitives, and CSG nodes to the owning
 * system. Those mutations are uncommitted until alea_void_add_cells()
 * succeeds. If the result is destroyed without committing,
 * alea_void_result_destroy() rolls the system back to its pre-generation
 * state using the snapshot fields below.
 */
typedef struct void_result {
    octree_node_t* root;    // Octree root (caller must free with octree_destroy)

    // Global void complement: ¬(Cell1 ∪ Cell2 ∪ ... ∪ CellN)
    // This is the exact analytical definition of void
    alea_node_id_t global_void;

    // Regional void regions: each is (global_void ∩ region_box).
    // These are CSG tree node IDs, NOT cells. Use alea_void_add_cells() to
    // register them as actual cells in the system.
    void_region_vec_t void_regions;

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

    // Snapshot for transactional generation. owner_sys is non-NULL while
    // generation owns un-committed mutations on the system. Once
    // alea_void_add_cells() succeeds, committed=true and destroy will not
    // roll back. If owner_sys is NULL (manually-built result, or already
    // rolled back), destroy skips rollback.
    alea_system_t* owner_sys;
    bool   committed;
    size_t snapshot_surfaces;
    size_t snapshot_nodes;
    size_t snapshot_primitives;
    int    snapshot_next_auto_surface_id;
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
 * @brief Add a graveyard cell enclosing the void bounds
 *
 * Creates a sphere that fully encloses the void bounding box, and registers
 * a cell outside the sphere with IMP:N=0, IMP:P=0. This is the standard
 * MCNP "rest of the world" cell that kills escaped particles.
 *
 * Must be called after alea_void_add_cells().
 *
 * @param sys CSG system
 * @param result Void generation result (uses root->bbox for bounds)
 * @return 1 on success, -1 on error
 */
int alea_void_add_graveyard(alea_system_t* sys, void_result_t* result);

/**
 * @brief Configuration for void cell merging
 */
typedef struct {
    double cell_weight;           /**< Cost per cell (default: 1.0) */
    double surface_weight;        /**< Cost per surface (default: 0.1) */
    int max_surfaces_per_cell;    /**< Hard limit on surfaces per cell (default: 24) */
    int min_cells;                /**< Stop merging below this count (default: 1) */
    bool use_greedy;              /**< true = old O(n^3) greedy, false = face-sorted O(n log^2 n) (default: false) */
    int consolidate_max_surfaces; /**< Consolidate void into single cell if it has ≤ this many surfaces (0 = off, default: 100) */
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

/**
 * @brief Install a test-only fault-injection hook for octree allocation.
 *
 * When `should_fail` is non-NULL, every call into `octree_node_create()`
 * invokes it; a true return makes that call simulate an OOM (returns NULL,
 * sets ALEA_ERR_OUT_OF_MEMORY). This is intended exclusively for tests that
 * exercise the partial-allocation rollback path. Pass NULL to disable.
 *
 * @param should_fail Predicate called per allocation, or NULL to disable
 * @param userdata Opaque value passed to should_fail
 */
void alea_void_set_octree_alloc_failure(bool (*should_fail)(void* ud), void* userdata);


#endif /* ALEA_VOID_H */