// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_SPATIAL_H
#define ALEA_SPATIAL_H

/**
 * @file alea_spatial.h
 * @brief Spatial index for cell instances
 *
 * Provides fast spatial queries without flattening the actual geometry.
 * The index stores (cell_index, transform, global_bbox) tuples - just metadata.
 *
 * Key insight: flatten the bboxes, not the geometry.
 *
 * Usage:
 *   alea_spatial_index_build(sys);              // Build once after loading
 *   alea_spatial_query_slice_z(sys, z, ...);    // Fast slice queries
 *   alea_spatial_query_region(sys, bbox, ...);  // Fast region queries
 */

#include "alea_types.h"
#include "alea_universe.h"
#include <stddef.h>
#include <stdbool.h>


/* Forward declarations */
typedef struct alea_system alea_system_t;

/* ============================================================================
 * CELL INSTANCE
 * ============================================================================ */

/**
 * @brief A cell instance in global coordinates
 *
 * Represents a specific occurrence of a cell in the flattened geometry.
 * The same cell can appear multiple times with different transforms
 * (via FILL references).
 */
typedef struct {
    uint32_t cell_index;       /* Index in sys->cells */
    alea_matrix_t transform;    /* Transform from cell-local to global coords */
    alea_bbox_t global_bbox;    /* Bounding box in global coordinates */
    int universe_id;           /* Universe this instance belongs to */
    int depth;                 /* Nesting depth (0 = base universe) */
    bool is_terminal;          /* True if no FILL (actual geometry) */
    uint32_t parent_cell_index; /* Index of parent FILL cell (UINT32_MAX if none) */
} alea_cell_instance_t;

/* ============================================================================
 * SPATIAL INDEX
 * ============================================================================ */

/**
 * @brief BVH node for spatial index (same structure as raycast BVH)
 */
typedef struct {
    alea_bbox_t bbox;           /* Bounding box */
    uint32_t left_or_first;    /* Left child OR first instance index */
    uint32_t right_child;      /* Right child (internal nodes only) */
    uint16_t count;            /* 0 = internal, >0 = leaf with count instances */
    uint8_t axis;              /* Split axis */
    uint8_t pad;
} alea_spatial_node_t;

ALEA_VEC_DEFINE(alea_cell_instance_vec, alea_cell_instance_t);
ALEA_VEC_DEFINE(alea_spatial_node_vec, alea_spatial_node_t);

/**
 * @brief Spatial index over cell instances
 */
typedef struct alea_spatial_index {
    /* Cell instances (flattened bboxes, not geometry) */
    alea_cell_instance_vec_t instances;

    /* BVH over instances */
    alea_spatial_node_vec_t nodes;
    uint32_t* indices;         /* Reordered instance indices */

    /* Global bounds */
    alea_bbox_t bounds;

    bool built;
} alea_spatial_index_t;

/* ============================================================================
 * BUILD API
 * ============================================================================ */

/**
 * @brief Build spatial index by traversing universe hierarchy
 *
 * Traverses all FILL references, computing global bboxes for each cell instance.
 * Does NOT flatten geometry - only stores metadata.
 *
 * @param sys CSG system (must have universe index built)
 * @return 0 on success, -1 on error
 */
int alea_spatial_index_build(alea_system_t* sys);

/**
 * @brief Free spatial index
 */
void alea_spatial_index_free(alea_spatial_index_t* idx);

/**
 * @brief Check if spatial index needs rebuild
 */
bool alea_spatial_index_needs_rebuild(const alea_system_t* sys);

/**
 * @brief Get instance count from spatial index
 * @return Number of instances, or 0 if no index
 */
size_t alea_spatial_get_instance_count(const alea_system_t* sys);

/* ============================================================================
 * QUERY API
 * ============================================================================ */

/**
 * @brief Query result for region/slice queries
 */
typedef struct {
    uint32_t instance_index;   /* Index into spatial_index->instances */
    uint32_t cell_index;       /* Index into sys->cells */
    int cell_id;               /* MCNP cell ID */
    int material_id;           /* Material ID */
    int universe_id;           /* Universe this instance belongs to */
    int depth;                 /* Nesting depth (0 = base universe) */
    bool is_terminal;          /* True if no FILL (actual geometry) */
    alea_matrix_t transform;    /* Transform to apply */
} alea_spatial_hit_t;

/**
 * @brief Find all cell instances intersecting a bounding box
 *
 * @param sys CSG system with spatial index
 * @param query_bbox Region to query
 * @param out_hits Output array
 * @param max_hits Maximum hits to return
 * @return Number of hits, or -1 on error
 */
int alea_spatial_query_region(alea_system_t* sys,
                             const alea_bbox_t* query_bbox,
                             alea_spatial_hit_t* out_hits,
                             size_t max_hits);

/**
 * @brief Find all cell instances intersecting a Z-plane slice
 *
 * Optimized query for axis-aligned slice planes.
 *
 * @param sys CSG system with spatial index
 * @param z Z coordinate of slice plane
 * @param x_min, x_max, y_min, y_max Viewport bounds
 * @param out_hits Output array
 * @param max_hits Maximum hits to return
 * @return Number of hits, or -1 on error
 */
int alea_spatial_query_slice_z(alea_system_t* sys,
                              double z,
                              double x_min, double x_max,
                              double y_min, double y_max,
                              alea_spatial_hit_t* out_hits,
                              size_t max_hits);

/**
 * @brief Find all cell instances whose bboxes contain a point
 *
 * Fast O(log N) query using BVH. Returns candidate cells that might
 * contain the point - caller should verify with actual containment test.
 *
 * @param sys CSG system with spatial index
 * @param x, y, z Point coordinates
 * @param out_hits Output array (candidates, sorted by depth)
 * @param max_hits Maximum hits to return
 * @return Number of candidates, or -1 on error
 */
int alea_spatial_query_point(alea_system_t* sys,
                            double x, double y, double z,
                            alea_spatial_hit_t* out_hits,
                            size_t max_hits);

/**
 * @brief Find all cells containing a point (with coherence caching)
 *
 * High-level function that finds all cells in the containment hierarchy.
 * Uses BVH for candidates + coherence caching for adjacent queries.
 *
 * @param sys CSG system with spatial index
 * @param x, y, z Point coordinates
 * @param out_hits Output array of cell hits
 * @param max_hits Maximum hits to return
 * @return Number of cells found, or -1 on error
 */
int alea_spatial_find_cells_at_point(alea_system_t* sys,
                                    double x, double y, double z,
                                    alea_cell_hit_t* out_hits,
                                    size_t max_hits);

/**
 * @brief Find cells at multiple points (batch query with coherence)
 *
 * Processes N points, writing one result per point into out_hits[i].
 * Exploits spatial coherence: the cache is maintained across consecutive
 * points, so spatially ordered inputs (e.g., scanline pixels) benefit
 * from cache reuse.
 *
 * @param sys CSG system with spatial index
 * @param points Array of 3D coordinates [x0,y0,z0, x1,y1,z1, ...]
 * @param n_points Number of points
 * @param out_hits Output: out_hits[i] receives the result for point i
 * @param max_hits_per_point Maximum hierarchy depth per point
 * @param out_counts Output: out_counts[i] = number of cells found at point i
 * @return 0 on success, -1 on error
 */
int alea_spatial_find_cells_batch(alea_system_t* sys,
                                  const double* points,
                                  size_t n_points,
                                  alea_cell_hit_t* out_hits,
                                  size_t max_hits_per_point,
                                  int* out_counts);

/**
 * @brief Reset coherence cache (call when starting new render)
 */
void alea_spatial_reset_cache(void);

/**
 * @brief Callback for spatial traversal
 */
typedef void (*alea_spatial_callback)(const alea_cell_instance_t* instance,
                                     uint32_t instance_index,
                                     void* userdata);

/**
 * @brief Traverse all instances intersecting a bbox
 *
 * Lower-level than query functions - calls callback for each hit.
 *
 * @param idx Spatial index
 * @param query_bbox Query region
 * @param callback Function to call for each hit
 * @param userdata User data passed to callback
 * @return Number of callbacks made
 */
int alea_spatial_traverse(const alea_spatial_index_t* idx,
                         const alea_bbox_t* query_bbox,
                         alea_spatial_callback callback,
                         void* userdata);


#endif /* ALEA_SPATIAL_H */
