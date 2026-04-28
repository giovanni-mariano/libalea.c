// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file bvh.c
 * @brief Bounding Volume Hierarchy implementation
 *
 * SAH-based construction and stack-based traversal for ray acceleration.
 */

#include "bvh.h"
#include "ray_epsilon.h"
#include "ray_bbox.h"
#include "core/alea_system.h"
#include "primitives/bbox.h"
#include "util/alea_log.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ============================================================================
 * INTERNAL STRUCTURES
 * ============================================================================ */

/* Surface info for building - stores centroid and bbox */
typedef struct {
    uint32_t index;         /* Original surface index */
    alea_bbox_t bbox;        /* Primitive bounding box */
    double centroid[3];     /* Bbox centroid */
} bvh_surface_info_t;

/* SAH bin for split evaluation */
typedef struct {
    int count;
    alea_bbox_t bbox;
} sah_bin_t;

/* Build context */
typedef struct {
    const alea_system_t* sys;
    bvh_surface_info_t* surfaces;
    size_t surface_count;
    alea_bvh_t* bvh;
} bvh_build_ctx_t;

/* ============================================================================
 * BBOX HELPERS
 * ============================================================================ */

static inline double bbox_centroid_axis(const alea_bbox_t* bbox, int axis) {
    switch (axis) {
        case 0: return (bbox->min_x + bbox->max_x) * 0.5;
        case 1: return (bbox->min_y + bbox->max_y) * 0.5;
        case 2: return (bbox->min_z + bbox->max_z) * 0.5;
        default: return 0;
    }
}

static inline double bbox_min_axis(const alea_bbox_t* bbox, int axis) {
    switch (axis) {
        case 0: return bbox->min_x;
        case 1: return bbox->min_y;
        case 2: return bbox->min_z;
        default: return 0;
    }
}

static inline double bbox_max_axis(const alea_bbox_t* bbox, int axis) {
    switch (axis) {
        case 0: return bbox->max_x;
        case 1: return bbox->max_y;
        case 2: return bbox->max_z;
        default: return 0;
    }
}

static alea_bbox_t compute_surfaces_bbox(const bvh_surface_info_t* surfaces,
                                        size_t start, size_t end) {
    alea_bbox_t result = alea_bbox_empty();
    for (size_t i = start; i < end; i++) {
        result = alea_bbox_union(&result, &surfaces[i].bbox);
    }
    return result;
}

static alea_bbox_t compute_centroids_bbox(const bvh_surface_info_t* surfaces,
                                         size_t start, size_t end) {
    alea_bbox_t result = alea_bbox_empty();
    for (size_t i = start; i < end; i++) {
        double cx = surfaces[i].centroid[0];
        double cy = surfaces[i].centroid[1];
        double cz = surfaces[i].centroid[2];
        if (cx < result.min_x) result.min_x = cx;
        if (cx > result.max_x) result.max_x = cx;
        if (cy < result.min_y) result.min_y = cy;
        if (cy > result.max_y) result.max_y = cy;
        if (cz < result.min_z) result.min_z = cz;
        if (cz > result.max_z) result.max_z = cz;
    }
    return result;
}

/* ============================================================================
 * NODE ALLOCATION
 * ============================================================================ */

static uint32_t alloc_node(alea_bvh_t* bvh) {
    alea_bvh_node_t* node = alea_vec_push_uninit(&bvh->nodes, alea_bvh_node_t);
    if (!node) return UINT32_MAX;
    return (uint32_t)(bvh->nodes.count - 1);
}

/* ============================================================================
 * SAH SPLIT EVALUATION
 * ============================================================================ */

static double evaluate_sah_split(const bvh_surface_info_t* surfaces,
                                 size_t start, size_t end,
                                 int axis,
                                 const alea_bbox_t* centroid_bbox,
                                 const alea_bbox_t* bounds,
                                 double* out_split_pos) {
    size_t count = end - start;
    if (count <= BVH_LEAF_THRESHOLD) {
        *out_split_pos = 0;
        return DBL_MAX;  /* Force leaf */
    }

    /* Compute extent along axis */
    double min_c = bbox_min_axis(centroid_bbox, axis);
    double max_c = bbox_max_axis(centroid_bbox, axis);
    double extent = max_c - min_c;
    if (extent < RAY_EPSILON) {
        *out_split_pos = min_c;  /* Median split */
        return DBL_MAX;
    }

    /* Initialize bins */
    sah_bin_t bins[BVH_SAH_BINS];
    for (int i = 0; i < BVH_SAH_BINS; i++) {
        bins[i].count = 0;
        bins[i].bbox = alea_bbox_empty();
    }

    /* Assign surfaces to bins */
    double scale = BVH_SAH_BINS / extent;

    for (size_t i = start; i < end; i++) {
        double c = surfaces[i].centroid[axis];
        int b = (int)((c - min_c) * scale);
        if (b < 0) b = 0;
        if (b >= BVH_SAH_BINS) b = BVH_SAH_BINS - 1;
        bins[b].count++;
        bins[b].bbox = alea_bbox_union(&bins[b].bbox, &surfaces[i].bbox);
    }

    /* Precompute suffix bounding boxes and counts from right-to-left (O(bins)) */
    alea_bbox_t suffix_bbox[BVH_SAH_BINS];
    int suffix_count[BVH_SAH_BINS];

    suffix_bbox[BVH_SAH_BINS - 1] = bins[BVH_SAH_BINS - 1].bbox;
    suffix_count[BVH_SAH_BINS - 1] = bins[BVH_SAH_BINS - 1].count;
    for (int i = BVH_SAH_BINS - 2; i >= 0; i--) {
        suffix_bbox[i] = alea_bbox_union(&suffix_bbox[i + 1], &bins[i].bbox);
        suffix_count[i] = suffix_count[i + 1] + bins[i].count;
    }

    /* Forward sweep: O(1) per split candidate */
    double costs[BVH_SAH_BINS - 1];
    alea_bbox_t left_bbox = alea_bbox_empty();
    int left_count = 0;
    double parent_area = alea_bbox_surface_area(bounds);

    for (int i = 0; i < BVH_SAH_BINS - 1; i++) {
        left_bbox = alea_bbox_union(&left_bbox, &bins[i].bbox);
        left_count += bins[i].count;

        double left_area = alea_bbox_surface_area(&left_bbox);
        double right_area = alea_bbox_surface_area(&suffix_bbox[i + 1]);
        int right_count = suffix_count[i + 1];

        if (parent_area < RAY_EPSILON) {
            costs[i] = DBL_MAX;
        } else {
            costs[i] = BVH_TRAVERSAL_COST +
                       BVH_INTERSECT_COST * (
                           (left_area / parent_area) * left_count +
                           (right_area / parent_area) * right_count
                       );
        }
    }

    /* Find minimum cost split */
    double min_cost = costs[0];
    int min_bin = 0;
    for (int i = 1; i < BVH_SAH_BINS - 1; i++) {
        if (costs[i] < min_cost) {
            min_cost = costs[i];
            min_bin = i;
        }
    }

    /* Compute actual split position */
    double split_pos = min_c + (min_bin + 1) * (extent / BVH_SAH_BINS);

    *out_split_pos = split_pos;
    return min_cost;
}

/* ============================================================================
 * PARTITION SURFACES
 * ============================================================================ */

static void partition_surfaces(bvh_surface_info_t* surfaces,
                               size_t start, size_t end,
                               int axis, double split_pos,
                               size_t* out_mid) {
    size_t left = start;
    size_t right = end - 1;

    while (left <= right && right < end) {
        while (left < end && surfaces[left].centroid[axis] < split_pos) {
            left++;
        }
        while (right > start && right < end &&
               surfaces[right].centroid[axis] >= split_pos) {
            right--;
        }
        if (left < right) {
            bvh_surface_info_t tmp = surfaces[left];
            surfaces[left] = surfaces[right];
            surfaces[right] = tmp;
            left++;
            right--;
        } else {
            break;
        }
    }

    *out_mid = left;
    if (*out_mid == start) *out_mid = start + 1;
    if (*out_mid == end) *out_mid = end - 1;
}

/* ============================================================================
 * RECURSIVE BUILD
 * ============================================================================ */

static uint32_t build_recursive(bvh_build_ctx_t* ctx,
                                size_t start, size_t end,
                                int depth) {
    uint32_t node_idx = alloc_node(ctx->bvh);
    if (node_idx == UINT32_MAX) return UINT32_MAX;

    alea_bvh_node_t* node = &ctx->bvh->nodes.data[node_idx];
    size_t count = end - start;

    /* Compute bounding box for this node */
    node->bbox = compute_surfaces_bbox(ctx->surfaces, start, end);

    /* Create leaf if few surfaces or max depth */
    if (count <= BVH_LEAF_THRESHOLD || depth > 30) {
        node->left_or_first = (uint32_t)start;
        node->surface_count = (uint16_t)count;
        node->axis = 0;
        return node_idx;
    }

    /* Find best split using SAH */
    alea_bbox_t centroid_bbox = compute_centroids_bbox(ctx->surfaces, start, end);

    double best_cost = DBL_MAX;
    int best_axis = 0;
    double best_split_pos = 0;

    for (int axis = 0; axis < 3; axis++) {
        double split_pos;
        double cost = evaluate_sah_split(ctx->surfaces, start, end,
                                         axis, &centroid_bbox, &node->bbox,
                                         &split_pos);
        if (cost < best_cost) {
            best_cost = cost;
            best_axis = axis;
            best_split_pos = split_pos;
        }
    }

    /* Compare with leaf cost */
    double leaf_cost = count * BVH_INTERSECT_COST;
    if (best_cost >= leaf_cost && count <= BVH_LEAF_THRESHOLD * 2) {
        node->left_or_first = (uint32_t)start;
        node->surface_count = (uint16_t)count;
        node->axis = 0;
        return node_idx;
    }

    /* Partition surfaces using the exact split position from SAH */
    size_t mid;
    partition_surfaces(ctx->surfaces, start, end, best_axis, best_split_pos, &mid);

    /* Ensure valid partition */
    if (mid == start || mid == end) {
        mid = start + count / 2;
    }

    /* Create internal node */
    node->axis = (uint8_t)best_axis;
    node->surface_count = 0;

    /* Build children */
    uint32_t left_idx = build_recursive(ctx, start, mid, depth + 1);
    uint32_t right_idx = build_recursive(ctx, mid, end, depth + 1);

    if (left_idx == UINT32_MAX || right_idx == UINT32_MAX) {
        return UINT32_MAX;
    }

    /* Store both child indices - must reacquire pointer after recursion (may have reallocated) */
    node = &ctx->bvh->nodes.data[node_idx];
    node->left_or_first = left_idx;
    node->right_child = right_idx;

    return node_idx;
}

/* ============================================================================
 * PUBLIC BUILD API
 * ============================================================================ */

alea_bvh_t* alea_bvh_build(const alea_system_t* sys) {
    if (!sys || alea_vec_count(&sys->surfaces) == 0) {
        return NULL;
    }

    alea_bvh_t* bvh = calloc(1, sizeof(alea_bvh_t));
    if (!bvh) return NULL;

    /* Allocate surface info array */
    bvh_surface_info_t* surfaces = malloc(alea_vec_count(&sys->surfaces) *
                                          sizeof(bvh_surface_info_t));
    if (!surfaces) {
        free(bvh);
        return NULL;
    }

    /* Initialize surface info with bboxes and centroids */
    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        const alea_surface_entry_t* surf = &sys->surfaces.data[i];
        const alea_primitive_entry_t* prim = &sys->primitives.data[surf->primitive_id];

        surfaces[i].index = (uint32_t)i;
        surfaces[i].bbox = alea_primitive_bbox(prim->type, &prim->data);

        surfaces[i].centroid[0] = (surfaces[i].bbox.min_x +
                                   surfaces[i].bbox.max_x) * 0.5;
        surfaces[i].centroid[1] = (surfaces[i].bbox.min_y +
                                   surfaces[i].bbox.max_y) * 0.5;
        surfaces[i].centroid[2] = (surfaces[i].bbox.min_z +
                                   surfaces[i].bbox.max_z) * 0.5;
    }

    /* Allocate initial node array */
    size_t initial_capacity = alea_vec_count(&sys->surfaces) * 2;
    alea_vec_init(&bvh->nodes);
    int vres = alea_vec_reserve(&bvh->nodes, initial_capacity, alea_bvh_node_t);
    if (vres != 0) {
        free(surfaces);
        free(bvh);
        return NULL;
    }

    /* Build tree */
    bvh_build_ctx_t ctx = {
        .sys = sys,
        .surfaces = surfaces,
        .surface_count = alea_vec_count(&sys->surfaces),
        .bvh = bvh
    };

    uint32_t root = build_recursive(&ctx, 0, alea_vec_count(&sys->surfaces), 0);
    if (root == UINT32_MAX) {
        free(surfaces);
        alea_bvh_free(bvh);
        return NULL;
    }

    /* Create surface index array (reordered by tree) */
    bvh->surface_indices = malloc(alea_vec_count(&sys->surfaces) * sizeof(uint32_t));
    if (!bvh->surface_indices) {
        free(surfaces);
        alea_bvh_free(bvh);
        return NULL;
    }

    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        bvh->surface_indices[i] = surfaces[i].index;
    }

    bvh->surface_count = alea_vec_count(&sys->surfaces);

    free(surfaces);
    return bvh;
}

void alea_bvh_free(alea_bvh_t* bvh) {
    if (!bvh) return;
    alea_vec_free(&bvh->nodes);
    free(bvh->surface_indices);
    free(bvh);
}

/* ============================================================================
 * RAY-AABB INTERSECTION (robust version for parallel rays)
 * ============================================================================ */

/* ray_bbox_slab from ray_bbox.h replaces the local ray_bbox_intersect */

/* ============================================================================
 * TRAVERSAL
 * ============================================================================ */

int alea_bvh_traverse(const alea_bvh_t* bvh,
                     const alea_ray_t* ray,
                     double t_min, double t_max,
                     alea_bvh_hit_callback callback,
                     void* userdata) {
    if (!bvh || !ray || !callback || bvh->nodes.count == 0) {
        return 0;
    }

    int dir_sign[3];
    dir_sign[0] = ray->dx < 0;
    dir_sign[1] = ray->dy < 0;
    dir_sign[2] = ray->dz < 0;

    /* Stack-based traversal */
    uint32_t stack[BVH_MAX_STACK_DEPTH];
    int sp = 0;
    int surfaces_tested = 0;

    stack[sp++] = 0;  /* Start at root */

    while (sp > 0) {
        uint32_t node_idx = stack[--sp];
        const alea_bvh_node_t* node = &bvh->nodes.data[node_idx];

        /* Test ray against node bbox (uses precomputed inv_d* from ray) */
        if (!ray_bbox_slab(ray, &node->bbox, t_min, t_max)) {
            continue;
        }

        if (node->surface_count > 0) {
            /* Leaf node: report all surfaces */
            for (uint16_t i = 0; i < node->surface_count; i++) {
                uint32_t surf_idx = bvh->surface_indices[node->left_or_first + i];
                callback(surf_idx, userdata);
                surfaces_tested++;
            }
        } else {
            /* Internal node: push children */
            uint32_t left = node->left_or_first;
            uint32_t right = node->right_child;

            /* Push far child first (so near child is processed first) */
            if (sp < BVH_MAX_STACK_DEPTH - 1) {
                if (dir_sign[node->axis]) {
                    stack[sp++] = left;
                    stack[sp++] = right;
                } else {
                    stack[sp++] = right;
                    stack[sp++] = left;
                }
            } else {
                ALEA_LOG_WARN("BVH traversal stack overflow - results may be incomplete");
            }
        }
    }

    return surfaces_tested;
}

int alea_bvh_traverse_batch(const alea_bvh_t* bvh,
                            const alea_ray_t* ray,
                            double t_min, double t_max,
                            alea_bvh_batch_callback callback,
                            void* userdata) {
    if (!bvh || !ray || !callback || bvh->nodes.count == 0) {
        return 0;
    }

    int dir_sign[3];
    dir_sign[0] = ray->dx < 0;
    dir_sign[1] = ray->dy < 0;
    dir_sign[2] = ray->dz < 0;

    uint32_t stack[BVH_MAX_STACK_DEPTH];
    int sp = 0;
    int surfaces_tested = 0;

    stack[sp++] = 0;

    while (sp > 0) {
        uint32_t node_idx = stack[--sp];
        const alea_bvh_node_t* node = &bvh->nodes.data[node_idx];

        if (!ray_bbox_slab(ray, &node->bbox, t_min, t_max)) {
            continue;
        }

        if (node->surface_count > 0) {
            /* Leaf node: report batch of surfaces */
            callback(&bvh->surface_indices[node->left_or_first],
                     node->surface_count, userdata);
            surfaces_tested += node->surface_count;
        } else {
            uint32_t left = node->left_or_first;
            uint32_t right = node->right_child;

            if (sp < BVH_MAX_STACK_DEPTH - 1) {
                if (dir_sign[node->axis]) {
                    stack[sp++] = left;
                    stack[sp++] = right;
                } else {
                    stack[sp++] = right;
                    stack[sp++] = left;
                }
            } else {
                ALEA_LOG_WARN("BVH traversal stack overflow - results may be incomplete");
            }
        }
    }

    return surfaces_tested;
}

/* ============================================================================
 * STATISTICS
 * ============================================================================ */

static void stats_recursive(const alea_bvh_t* bvh, uint32_t node_idx,
                            int depth, size_t* leaf_count, size_t* max_depth) {
    if (node_idx >= bvh->nodes.count) return;

    const alea_bvh_node_t* node = &bvh->nodes.data[node_idx];

    if ((size_t)depth > *max_depth) {
        *max_depth = depth;
    }

    if (node->surface_count > 0) {
        (*leaf_count)++;
    } else {
        stats_recursive(bvh, node->left_or_first, depth + 1, leaf_count, max_depth);
        stats_recursive(bvh, node->right_child, depth + 1, leaf_count, max_depth);
    }
}

void alea_bvh_stats(const alea_bvh_t* bvh,
                   size_t* out_node_count,
                   size_t* out_leaf_count,
                   size_t* out_max_depth) {
    if (!bvh) {
        if (out_node_count) *out_node_count = 0;
        if (out_leaf_count) *out_leaf_count = 0;
        if (out_max_depth) *out_max_depth = 0;
        return;
    }

    if (out_node_count) *out_node_count = bvh->nodes.count;

    if (out_leaf_count || out_max_depth) {
        size_t leaf_count = 0;
        size_t max_depth = 0;
        stats_recursive(bvh, 0, 0, &leaf_count, &max_depth);
        if (out_leaf_count) *out_leaf_count = leaf_count;
        if (out_max_depth) *out_max_depth = max_depth;
    }
}
