// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_BVH_H
#define ALEA_BVH_H

#include "alea_types.h"
#include "raycast.h"
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file bvh.h
 * @brief Bounding Volume Hierarchy for raycast acceleration
 *
 * Provides O(log N) ray-surface intersection tests instead of O(N).
 * Uses SAH (Surface Area Heuristic) for optimal tree construction.
 */

/* ============================================================================
 * BVH NODE STRUCTURE
 * ============================================================================ */

/**
 * @brief BVH node - 64 bytes, cache-friendly layout
 *
 * For leaf nodes: surface_count > 0, left_or_first is index into surface_indices
 * For internal nodes: surface_count = 0, left_or_first is left child index,
 *                     right_child is right child index
 */
typedef struct {
    alea_bbox_t bbox;           /* Bounding box (48 bytes) */
    uint32_t left_or_first;    /* Left child index OR first surface index */
    uint32_t right_child;      /* Right child index (only for internal nodes) */
    uint16_t surface_count;    /* 0 = internal node, >0 = leaf node */
    uint8_t axis;              /* Split axis (0=X, 1=Y, 2=Z) */
    uint8_t pad;               /* Padding for alignment */
} alea_bvh_node_t;

/**
 * @brief Complete BVH structure
 */
typedef struct alea_bvh {
    alea_bvh_node_t* nodes;      /* Node array */
    uint32_t* surface_indices;  /* Reordered surface indices */
    size_t node_count;          /* Number of nodes used */
    size_t node_capacity;       /* Capacity of node array */
    size_t surface_count;       /* Number of surfaces when built */
} alea_bvh_t;

/* ============================================================================
 * BVH BUILD CONFIGURATION
 * ============================================================================ */

/* Set to 1 to disable BVH and use linear scan (for debugging) */
#define BVH_DISABLED           0  /* Enable BVH acceleration */

#define BVH_LEAF_THRESHOLD     4    /* Max surfaces per leaf node */
#define BVH_SAH_BINS          16    /* Number of SAH split candidates per axis */
#define BVH_TRAVERSAL_COST   1.0    /* Relative cost of traversing a node */
#define BVH_INTERSECT_COST   1.0    /* Relative cost of surface intersection */
#define BVH_MAX_STACK_DEPTH  128    /* Maximum traversal stack depth */

/* ============================================================================
 * BVH API
 * ============================================================================ */

/**
 * @brief Build BVH from system surfaces
 *
 * Constructs an optimized BVH using SAH (Surface Area Heuristic).
 * The BVH is built over all surfaces in the system.
 *
 * @param sys CSG system containing surfaces
 * @return Pointer to new BVH, or NULL on failure
 */
alea_bvh_t* alea_bvh_build(const alea_system_t* sys);

/**
 * @brief Free BVH memory
 *
 * @param bvh BVH to free (can be NULL)
 */
void alea_bvh_free(alea_bvh_t* bvh);

/**
 * @brief Callback type for BVH traversal (per-surface)
 *
 * Called for each surface that the ray potentially intersects.
 * The callback should perform the actual surface intersection test.
 *
 * @param surface_idx Index into sys->surfaces
 * @param userdata User-provided context
 */
typedef void (*alea_bvh_hit_callback)(uint32_t surface_idx, void* userdata);

/**
 * @brief Batch callback type for BVH traversal (per-leaf)
 *
 * Called once per leaf node with all surfaces in the leaf.
 * Reduces function-call overhead for leaf nodes with up to
 * BVH_LEAF_THRESHOLD surfaces.
 *
 * @param surface_indices Pointer into reordered index array
 * @param count Number of surfaces in this leaf (1..BVH_LEAF_THRESHOLD)
 * @param userdata User-provided context
 */
typedef void (*alea_bvh_batch_callback)(const uint32_t* surface_indices,
                                        uint16_t count, void* userdata);

/**
 * @brief Traverse BVH and report candidate surfaces (per-surface callback)
 *
 * @param bvh BVH to traverse
 * @param ray Ray to test
 * @param t_min Minimum ray parameter
 * @param t_max Maximum ray parameter
 * @param callback Function to call for each candidate surface
 * @param userdata User data passed to callback
 * @return Number of surfaces tested (callback invocations)
 */
int alea_bvh_traverse(const alea_bvh_t* bvh,
                     const alea_ray_t* ray,
                     double t_min, double t_max,
                     alea_bvh_hit_callback callback,
                     void* userdata);

/**
 * @brief Traverse BVH with batch callback (per-leaf)
 *
 * Like alea_bvh_traverse but calls the callback once per leaf node
 * with all surfaces in the leaf, reducing function-call overhead.
 *
 * @return Number of surfaces tested
 */
int alea_bvh_traverse_batch(const alea_bvh_t* bvh,
                            const alea_ray_t* ray,
                            double t_min, double t_max,
                            alea_bvh_batch_callback callback,
                            void* userdata);

/* alea_bvh_needs_rebuild removed: always returned 0 due to circular dependency.
 * Rebuild check is handled directly in alea_raycast_ensure_caches(). */

/**
 * @brief Get BVH statistics
 *
 * @param bvh BVH to query
 * @param out_node_count Output: total nodes
 * @param out_leaf_count Output: leaf nodes
 * @param out_max_depth Output: maximum tree depth
 */
void alea_bvh_stats(const alea_bvh_t* bvh,
                   size_t* out_node_count,
                   size_t* out_leaf_count,
                   size_t* out_max_depth);

#ifdef __cplusplus
}
#endif

#endif /* ALEA_BVH_H */
