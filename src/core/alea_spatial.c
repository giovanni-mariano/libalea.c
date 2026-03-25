// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_spatial.c
 * @brief Spatial index for cell instances
 *
 * The core idea: traverse the universe hierarchy ONCE to collect all cell
 * instances with their global bboxes. Build a BVH over these bboxes.
 * Then queries are fast - just traverse the BVH.
 *
 * This avoids flattening the actual geometry while still enabling
 * fast spatial queries for interactive visualization.
 */

#include "alea_spatial.h"
#include "alea_eval.h"
#include "core/alea_system.h"
#include "core/alea_surface.h"
#include "primitives/bbox.h"
#include "primitives/primitive_desc.h"
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <sys/time.h>
#include "util/alea_log.h"

/* OpenMP removed — centroid init is too trivial to benefit from threading */

#define MAX_DEPTH 32
#define BVH_LEAF_SIZE 4
#define EPSILON 1e-10

/* ============================================================================
 * INTERVAL CULLING (EXPERIMENTAL)
 * ============================================================================
 *
 * Set ALEA_USE_INTERVAL_CULLING=1 to enable interval-based pre-filtering.
 * This uses interval arithmetic on CSG trees to quickly reject candidates
 * that are definitely outside, before doing expensive containment tests.
 *
 * To disable: set to 0 or remove this section entirely.
 */
#ifndef ALEA_USE_INTERVAL_CULLING
#define ALEA_USE_INTERVAL_CULLING 0
#endif

#if ALEA_USE_INTERVAL_CULLING
#include "slice/curve_intersect.h"
static _Thread_local size_t g_interval_culls = 0;
size_t alea_spatial_get_interval_culls(void) { return g_interval_culls; }
#endif

/* ============================================================================
 * BBOX HELPERS
 * ============================================================================ */

/**
 * Transform AABB by decomposing rotation matrix into min/max contributions.
 * Avoids projecting all 8 corners (18 mul + 18 add vs 192 mul-adds).
 * Matrix layout: m[row*4+col], row-major 3x4 with translation in col 3.
 */
static alea_bbox_t bbox_transform(const alea_bbox_t* bbox, const alea_matrix_t* mat) {
    const double bmin[3] = {bbox->min_x, bbox->min_y, bbox->min_z};
    const double bmax[3] = {bbox->max_x, bbox->max_y, bbox->max_z};
    double rmin[3], rmax[3];

    for (int i = 0; i < 3; i++) {
        rmin[i] = rmax[i] = mat->m[i * 4 + 3]; /* translation */
        for (int j = 0; j < 3; j++) {
            double e = mat->m[i * 4 + j] * bmin[j];
            double f = mat->m[i * 4 + j] * bmax[j];
            if (e < f) { rmin[i] += e; rmax[i] += f; }
            else       { rmin[i] += f; rmax[i] += e; }
        }
    }

    return (alea_bbox_t){rmin[0], rmax[0], rmin[1], rmax[1], rmin[2], rmax[2]};
}

static bool bbox_intersects(const alea_bbox_t* a, const alea_bbox_t* b) {
    return a->min_x <= b->max_x && a->max_x >= b->min_x &&
           a->min_y <= b->max_y && a->max_y >= b->min_y &&
           a->min_z <= b->max_z && a->max_z >= b->min_z;
}


/* ============================================================================
 * INSTANCE COLLECTION (UNIVERSE TRAVERSAL)
 * ============================================================================ */

typedef struct {
    alea_system_t* sys;
    alea_spatial_index_t* idx;
    int error;
} collect_ctx_t;

typedef struct {
    alea_system_t* sys;
    size_t total_instances;
    size_t terminal_instances;
    size_t container_instances;
    int max_depth;
    int error;
} count_ctx_t;

static double monotonic_seconds(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;
}

static double bytes_to_mib(size_t bytes) {
    return (double)bytes / (1024.0 * 1024.0);
}

static int ensure_instance_capacity(alea_spatial_index_t* idx, size_t needed) {
    size_t min_cap = idx->instances.count + needed;
    if (min_cap <= idx->instances.capacity) return 0;
    alea_result_t res = alea_vec_reserve(&idx->instances, min_cap, alea_cell_instance_t);
    return ALEA_IS_ERR(res) ? -1 : 0;
}

/**
 * Compute bbox for a cell by examining its CSG tree
 */
static alea_bbox_t compute_cell_bbox(alea_system_t* sys, uint32_t cell_index) {
    const alea_cell_entry_t* cell = &sys->cells.data[cell_index];

    /* If cell has no geometry (void), return empty bbox */
    if (cell->root_node_id == ALEA_NODE_ID_INVALID) {
        return alea_bbox_empty();
    }

    /* Get bbox from the node if it's cached (empty bbox has min > max) */
    const alea_node_t* root = &sys->nodes.data[cell->root_node_id];
    if (root->bbox.min_x <= root->bbox.max_x) {
        return root->bbox;
    }

    /* Compute bbox properly using CSG tree traversal with halfspace awareness.
     * This handles intersections of halfspaces (like plane intersections)
     * correctly, giving finite bboxes for cells bounded by planes. */
    alea_bbox_t bbox = alea_get_bbox(sys, cell->root_node_id);

    /* If still infinite or invalid, use a large default */
    if (bbox.min_x > bbox.max_x ||
        bbox.min_x <= -1e10 || bbox.max_x >= 1e10 ||
        bbox.min_y <= -1e10 || bbox.max_y >= 1e10 ||
        bbox.min_z <= -1e10 || bbox.max_z >= 1e10) {
        bbox.min_x = bbox.min_y = bbox.min_z = -1e6;
        bbox.max_x = bbox.max_y = bbox.max_z = 1e6;
    }

    return bbox;
}

static void count_instances_recursive(count_ctx_t* ctx,
                                      int universe_id,
                                      const alea_matrix_t* accumulated,
                                      int depth,
                                      const alea_bbox_t* parent_global_bbox) {
    if (ctx->error) return;
    if (g_alea_interrupted) { ctx->error = -1; return; }
    if (depth >= MAX_DEPTH) return;

    if (depth > ctx->max_depth) {
        ctx->max_depth = depth;
    }

    const alea_universe_t* univ = alea_get_universe(ctx->sys, universe_id);
    if (!univ) return;

    for (size_t i = 0; i < univ->cell_indices.count; i++) {
        size_t cell_idx = univ->cell_indices.data[i];
        const alea_cell_entry_t* cell = &ctx->sys->cells.data[cell_idx];

        alea_bbox_t local_bbox = compute_cell_bbox(ctx->sys, (uint32_t)cell_idx);
        alea_bbox_t global_bbox;

        if (accumulated) {
            global_bbox = bbox_transform(&local_bbox, accumulated);
        } else {
            global_bbox = local_bbox;
        }

        if (parent_global_bbox) {
            global_bbox = alea_bbox_intersection(&global_bbox, parent_global_bbox);
            if (!alea_bbox_is_valid(&global_bbox)) continue;
        }

        ctx->total_instances++;
        if (cell->fill_universe > 0) ctx->container_instances++;
        else ctx->terminal_instances++;

        if (cell->fill_universe > 0) {
            alea_matrix_t fill_mat;

            if (cell->fill_transform > 0) {
                const alea_transform_t* tr = alea_get_transform(ctx->sys,
                                                              cell->fill_transform);
                if (tr) {
                    alea_matrix_from_mcnp(&fill_mat, tr->cosines,
                                        tr->value_count, false);
                } else {
                    alea_matrix_identity(&fill_mat);
                }
            } else {
                alea_matrix_identity(&fill_mat);
            }

            alea_matrix_t new_accumulated;
            if (accumulated) {
                alea_matrix_multiply(&new_accumulated, accumulated, &fill_mat);
            } else {
                new_accumulated = fill_mat;
            }

            alea_matrix_invert(&new_accumulated);

            count_instances_recursive(ctx, cell->fill_universe,
                                      &new_accumulated, depth + 1,
                                      &global_bbox);
        }
    }
}

/**
 * Recursive traversal to collect all cell instances
 * @param parent_cell_idx Index of parent FILL cell (UINT32_MAX if at depth 0)
 * @param parent_global_bbox Global bbox of the parent FILL cell (NULL at depth 0)
 */
static void collect_instances_recursive(collect_ctx_t* ctx,
                                        int universe_id,
                                        const alea_matrix_t* accumulated,
                                        int depth,
                                        uint32_t parent_cell_idx,
                                        const alea_bbox_t* parent_global_bbox) {
    if (ctx->error) return;
    if (g_alea_interrupted) { ctx->error = -1; return; }
    if (depth >= MAX_DEPTH) return;

    const alea_universe_t* univ = alea_get_universe(ctx->sys, universe_id);
    if (!univ) return;

    for (size_t i = 0; i < univ->cell_indices.count; i++) {
        size_t cell_idx = univ->cell_indices.data[i];
        const alea_cell_entry_t* cell = &ctx->sys->cells.data[cell_idx];

        /* Compute global bbox for this cell */
        alea_bbox_t local_bbox = compute_cell_bbox(ctx->sys, cell_idx);
        alea_bbox_t global_bbox;

        if (accumulated) {
            global_bbox = bbox_transform(&local_bbox, accumulated);
        } else {
            global_bbox = local_bbox;
        }

        /* Clip against parent cell's global bbox (parent constrains children) */
        if (parent_global_bbox) {
            global_bbox = alea_bbox_intersection(&global_bbox, parent_global_bbox);
            if (!alea_bbox_is_valid(&global_bbox)) continue;
        }

        /* Add instance */
        if (ensure_instance_capacity(ctx->idx, 1) != 0) {
            ctx->error = -1;
            return;
        }

        alea_cell_instance_t* inst = &ctx->idx->instances.data[ctx->idx->instances.count++];
        inst->cell_index = (uint32_t)cell_idx;
        inst->global_bbox = global_bbox;
        inst->universe_id = universe_id;
        inst->depth = depth;
        inst->is_terminal = (cell->fill_universe <= 0);
        inst->parent_cell_index = parent_cell_idx;

        if (accumulated) {
            inst->transform = *accumulated;
        } else {
            alea_matrix_identity(&inst->transform);
        }

        /* Update global bounds */
        ctx->idx->bounds = alea_bbox_union(&ctx->idx->bounds, &global_bbox);

        /* Recurse into FILL */
        if (cell->fill_universe > 0) {
            alea_matrix_t fill_mat;

            /* Get fill transform */
            if (cell->fill_transform > 0) {
                const alea_transform_t* tr = alea_get_transform(ctx->sys,
                                                              cell->fill_transform);
                if (tr) {
                    /* Use tr->cosines which has pre-computed direction cosines */
                    alea_matrix_from_mcnp(&fill_mat, tr->cosines,
                                        tr->value_count, false);
                } else {
                    alea_matrix_identity(&fill_mat);
                }
            } else {
                alea_matrix_identity(&fill_mat);
            }

            /* Compose transforms */
            alea_matrix_t new_accumulated;
            if (accumulated) {
                alea_matrix_multiply(&new_accumulated, accumulated, &fill_mat);
            } else {
                new_accumulated = fill_mat;
            }

            alea_matrix_invert(&new_accumulated);

            /* Recurse - pass current cell as parent with its global bbox */
            collect_instances_recursive(ctx, cell->fill_universe,
                                       &new_accumulated, depth + 1,
                                       (uint32_t)cell_idx,
                                       &global_bbox);
        }
    }
}

/* ============================================================================
 * BVH CONSTRUCTION
 * ============================================================================ */

typedef struct {
    uint32_t index;
    double centroid[3];
} bvh_item_t;

static int ensure_node_capacity(alea_spatial_index_t* idx, size_t needed) {
    size_t min_cap = idx->nodes.count + needed;
    if (min_cap <= idx->nodes.capacity) return 0;
    alea_result_t res = alea_vec_reserve(&idx->nodes, min_cap, alea_spatial_node_t);
    return ALEA_IS_ERR(res) ? -1 : 0;
}

static alea_bbox_t compute_items_bbox(const alea_spatial_index_t* idx,
                                      const bvh_item_t* items,
                                      size_t start, size_t end) {
    alea_bbox_t result = alea_bbox_empty();
    for (size_t i = start; i < end; i++) {
        const alea_bbox_t* bbox = &idx->instances.data[items[i].index].global_bbox;
        result = alea_bbox_union(&result, bbox);
    }
    return result;
}

static int find_best_split_axis(const bvh_item_t* items, size_t start, size_t end) {
    alea_bbox_t centroid_bbox = alea_bbox_empty();
    for (size_t i = start; i < end; i++) {
        double cx = items[i].centroid[0];
        double cy = items[i].centroid[1];
        double cz = items[i].centroid[2];
        if (cx < centroid_bbox.min_x) centroid_bbox.min_x = cx;
        if (cx > centroid_bbox.max_x) centroid_bbox.max_x = cx;
        if (cy < centroid_bbox.min_y) centroid_bbox.min_y = cy;
        if (cy > centroid_bbox.max_y) centroid_bbox.max_y = cy;
        if (cz < centroid_bbox.min_z) centroid_bbox.min_z = cz;
        if (cz > centroid_bbox.max_z) centroid_bbox.max_z = cz;
    }

    double dx = centroid_bbox.max_x - centroid_bbox.min_x;
    double dy = centroid_bbox.max_y - centroid_bbox.min_y;
    double dz = centroid_bbox.max_z - centroid_bbox.min_z;

    if (dx >= dy && dx >= dz) return 0;
    if (dy >= dz) return 1;
    return 2;
}

/**
 * Quickselect (Hoare's nth_element) — O(n) amortized median partition.
 * After return, items[mid] is the median and all items in [start, mid)
 * have centroid[axis] <= items[mid], all in [mid, end) have >= items[mid].
 */
static void quickselect(bvh_item_t* items, size_t lo, size_t hi,
                        size_t target, int axis) {
    while (lo + 1 < hi) {
        /* Median-of-three pivot selection */
        size_t mid = lo + (hi - lo) / 2;
        if (items[mid].centroid[axis] < items[lo].centroid[axis]) {
            bvh_item_t tmp = items[lo]; items[lo] = items[mid]; items[mid] = tmp;
        }
        if (items[hi-1].centroid[axis] < items[lo].centroid[axis]) {
            bvh_item_t tmp = items[lo]; items[lo] = items[hi-1]; items[hi-1] = tmp;
        }
        if (items[mid].centroid[axis] < items[hi-1].centroid[axis]) {
            bvh_item_t tmp = items[mid]; items[mid] = items[hi-1]; items[hi-1] = tmp;
        }
        /* Pivot is now at hi-1 */
        double pivot = items[hi-1].centroid[axis];

        size_t i = lo, j = hi - 2;
        for (;;) {
            while (i < hi && items[i].centroid[axis] < pivot) i++;
            while (j > lo && items[j].centroid[axis] > pivot) j--;
            if (i >= j) break;
            bvh_item_t tmp = items[i]; items[i] = items[j]; items[j] = tmp;
            i++; j--;
        }
        /* Place pivot */
        bvh_item_t tmp = items[i]; items[i] = items[hi-1]; items[hi-1] = tmp;

        if (i == target) return;
        if (target < i) hi = i;
        else lo = i + 1;
    }
}

static void partition_items(bvh_item_t* items, size_t start, size_t end,
                           int axis, size_t* out_mid) {
    size_t mid = (start + end) / 2;
    if (end - start > 2) {
        quickselect(items, start, end, mid, axis);
    }
    *out_mid = mid;
}

static size_t estimate_bvh_node_count(size_t item_count) {
    if (item_count == 0) return 0;
    if (item_count <= BVH_LEAF_SIZE) return 1;

    size_t left = item_count / 2;
    size_t right = item_count - left;
    return 1 + estimate_bvh_node_count(left) + estimate_bvh_node_count(right);
}

static uint32_t build_bvh_recursive(alea_spatial_index_t* idx,
                                    bvh_item_t* items,
                                    size_t start, size_t end,
                                    int depth) {
    if (ensure_node_capacity(idx, 1) != 0) {
        return UINT32_MAX;
    }

    uint32_t node_idx = (uint32_t)idx->nodes.count++;
    alea_spatial_node_t* node = &idx->nodes.data[node_idx];
    size_t count = end - start;

    /* Create leaf if small enough */
    if (count <= BVH_LEAF_SIZE || depth > 30) {
        node->bbox = compute_items_bbox(idx, items, start, end);
        node->left_or_first = (uint32_t)start;
        node->count = (uint16_t)count;
        node->axis = 0;
        return node_idx;
    }

    /* Find best split */
    int axis = find_best_split_axis(items, start, end);
    size_t mid;
    partition_items(items, start, end, axis, &mid);

    if (mid == start || mid == end) {
        mid = (start + end) / 2;
    }

    /* Internal node */
    node->axis = (uint8_t)axis;
    node->count = 0;

    uint32_t left = build_bvh_recursive(idx, items, start, mid, depth + 1);
    uint32_t right = build_bvh_recursive(idx, items, mid, end, depth + 1);

    if (left == UINT32_MAX || right == UINT32_MAX) {
        return UINT32_MAX;
    }

    /* Reacquire pointer after potential realloc; compute bbox from children */
    node = &idx->nodes.data[node_idx];
    node->bbox = alea_bbox_union(&idx->nodes.data[left].bbox, &idx->nodes.data[right].bbox);
    node->left_or_first = left;
    node->right_child = right;

    return node_idx;
}

static int build_bvh(alea_spatial_index_t* idx) {
    if (idx->instances.count == 0) {
        return 0;
    }

    /* Prepare centroid refs rather than copying full bboxes */
    bvh_item_t* items = malloc(idx->instances.count * sizeof(bvh_item_t));
    if (!items) return -1;

    for (size_t i = 0; i < idx->instances.count; i++) {
        const alea_bbox_t* bbox = &idx->instances.data[i].global_bbox;
        items[i].index = (uint32_t)i;
        items[i].centroid[0] = (bbox->min_x + bbox->max_x) * 0.5;
        items[i].centroid[1] = (bbox->min_y + bbox->max_y) * 0.5;
        items[i].centroid[2] = (bbox->min_z + bbox->max_z) * 0.5;
    }

    /* Reserve the exact node count produced by this median-split builder. */
    alea_vec_init(&idx->nodes);
    size_t estimated_nodes = estimate_bvh_node_count(idx->instances.count);
    alea_result_t nres = alea_vec_reserve(&idx->nodes, estimated_nodes, alea_spatial_node_t);
    if (ALEA_IS_ERR(nres)) {
        free(items);
        return -1;
    }

    /* Build tree */
    uint32_t root = build_bvh_recursive(idx, items, 0, idx->instances.count, 0);

    /* Store reordered indices */
    idx->indices = malloc(idx->instances.count * sizeof(uint32_t));
    if (!idx->indices) {
        free(items);
        return -1;
    }
    for (size_t i = 0; i < idx->instances.count; i++) {
        idx->indices[i] = items[i].index;
    }

    free(items);
    return (root == UINT32_MAX) ? -1 : 0;
}

/* ============================================================================
 * PUBLIC API
 * ============================================================================ */

/**
 * Internal build — does the actual work.
 */
static int spatial_index_build_impl(alea_system_t* sys) {
    double t_start = monotonic_seconds();

    /* Ensure universe index is built */
    if (!sys->universe_index_built) {
        if (alea_build_universe_index(sys) != 0) {
            return -1;
        }
    }

    /* Ensure cell surface index is built (needed by slice curve generation
     * which accesses cell->surface_indices after spatial queries) */
    if (alea_build_cell_surface_index(sys) != 0) {
        return -1;
    }

    /* Free existing index if any */
    if (sys->spatial_index) {
        alea_spatial_index_free(sys->spatial_index);
        sys->spatial_index = NULL;
    }

    double t_count_start = monotonic_seconds();
    count_ctx_t count_ctx = {
        .sys = sys,
        .max_depth = -1,
        .error = 0
    };
    count_instances_recursive(&count_ctx, 0, NULL, 0, NULL);
    if (count_ctx.error) {
        return -1;
    }
    double t_count_end = monotonic_seconds();

    /* Allocate new index */
    alea_spatial_index_t* idx = calloc(1, sizeof(alea_spatial_index_t));
    if (!idx) return -1;

    alea_vec_init(&idx->instances);
    alea_result_t ires = alea_vec_reserve(&idx->instances, count_ctx.total_instances, alea_cell_instance_t);
    if (ALEA_IS_ERR(ires)) {
        free(idx);
        return -1;
    }

    idx->bounds = alea_bbox_empty();

    /* Collect all instances by traversing universe hierarchy */
    double t_collect_start = monotonic_seconds();
    collect_ctx_t ctx = {
        .sys = sys,
        .idx = idx,
        .error = 0
    };

    collect_instances_recursive(&ctx, 0, NULL, 0, UINT32_MAX, NULL);

    if (ctx.error) {
        alea_spatial_index_free(idx);
        return -1;
    }

    double t_collect_end = monotonic_seconds();

    size_t estimated_nodes = estimate_bvh_node_count(idx->instances.count);
    size_t instance_bytes = idx->instances.capacity * sizeof(alea_cell_instance_t);
    size_t item_bytes = idx->instances.count * sizeof(bvh_item_t);
    size_t node_bytes = estimated_nodes * sizeof(alea_spatial_node_t);
    size_t index_bytes = idx->instances.count * sizeof(uint32_t);

    ALEA_LOG_INFO("Spatial index build: %zu instances (%zu terminal, %zu container), max depth %d",
                  count_ctx.total_instances, count_ctx.terminal_instances,
                  count_ctx.container_instances, count_ctx.max_depth);
    ALEA_LOG_INFO("Spatial index memory estimate: instances=%.1f MiB, bvh_refs=%.1f MiB, nodes=%.1f MiB, indices=%.1f MiB, peak~=%.1f MiB",
                  bytes_to_mib(instance_bytes), bytes_to_mib(item_bytes),
                  bytes_to_mib(node_bytes), bytes_to_mib(index_bytes),
                  bytes_to_mib(instance_bytes + item_bytes + node_bytes + index_bytes));
    ALEA_LOG_DEBUG("Spatial index timing: count=%.3fs collect=%.3fs",
                   t_count_end - t_count_start, t_collect_end - t_collect_start);

    /* Build BVH over instances */
    double t_bvh_start = monotonic_seconds();
    if (build_bvh(idx) != 0) {
        alea_spatial_index_free(idx);
        return -1;
    }
    double t_bvh_end = monotonic_seconds();

    idx->built = true;
    sys->spatial_index = idx;

    ALEA_LOG_INFO("Spatial index ready: %zu nodes, total build %.3fs (bvh %.3fs)",
                  idx->nodes.count, t_bvh_end - t_start, t_bvh_end - t_bvh_start);

    return 0;
}

/**
 * Thread-safe spatial index build.
 *
 * Uses atomic CAS on sys->spatial_build_state to ensure only one thread
 * builds at a time. Other concurrent callers back off until the build
 * is complete.
 *
 * States: 0 = not built, 1 = building, 2 = done
 */
int alea_spatial_index_build(alea_system_t* sys) {
    if (!sys) return -1;

    /* Fast path: already built */
    if (atomic_load(&sys->spatial_build_state) == 2 &&
        sys->spatial_index && sys->spatial_index->built) {
        return 0;
    }

    /* Try to claim the build */
    int expected = 0;
    if (atomic_compare_exchange_strong(&sys->spatial_build_state, &expected, 1)) {
        /* We won the race — do the build */
        int rc = spatial_index_build_impl(sys);
        if (rc == 0) {
            atomic_store(&sys->spatial_build_state, 2);
        } else {
            atomic_store(&sys->spatial_build_state, 0);  /* reset on failure */
        }
        return rc;
    }

    /* Another thread is building — spin until done */
    while (atomic_load(&sys->spatial_build_state) == 1) {
        /* empty — contention is rare (build happens once) */
    }

    /* Check if the builder succeeded */
    return (sys->spatial_index && sys->spatial_index->built) ? 0 : -1;
}

void alea_spatial_index_free(alea_spatial_index_t* idx) {
    if (!idx) return;
    alea_vec_free(&idx->instances);
    alea_vec_free(&idx->nodes);
    free(idx->indices);
    free(idx);
}

size_t alea_spatial_get_instance_count(const alea_system_t* sys) {
    if (!sys || !sys->spatial_index) return 0;
    return sys->spatial_index->instances.count;
}

bool alea_spatial_index_needs_rebuild(const alea_system_t* sys) {
    if (!sys || !sys->spatial_index) return true;
    return !sys->spatial_index->built;
}

/* ============================================================================
 * TRAVERSAL
 * ============================================================================ */

/* BVH build caps at depth 30. Each traversal step pushes at most 2 children
 * and pops 1, so max stack depth = build depth + 2. */
#define TRAVERSE_STACK_SIZE (MAX_DEPTH + 2)

int alea_spatial_traverse(const alea_spatial_index_t* idx,
                         const alea_bbox_t* query_bbox,
                         alea_spatial_callback callback,
                         void* userdata) {
    if (!idx || !query_bbox || !callback || !idx->built) {
        return -1;
    }
    if (idx->nodes.count == 0) return 0;

    uint32_t stack[TRAVERSE_STACK_SIZE];
    int sp = 0;
    int count = 0;
    stack[sp++] = 0; /* root */

    while (sp > 0) {
        uint32_t ni = stack[--sp];
        const alea_spatial_node_t* node = &idx->nodes.data[ni];

        if (!bbox_intersects(&node->bbox, query_bbox))
            continue;

        if (node->count > 0) {
            /* Leaf */
            for (uint16_t i = 0; i < node->count; i++) {
                uint32_t inst_idx = idx->indices[node->left_or_first + i];
                const alea_cell_instance_t* inst = &idx->instances.data[inst_idx];
                if (bbox_intersects(&inst->global_bbox, query_bbox)) {
                    callback(inst, inst_idx, userdata);
                    count++;
                }
            }
        } else {
            /* Internal — push right first so left is popped first */
            assert(sp + 2 <= TRAVERSE_STACK_SIZE &&
                   "BVH deeper than traversal stack");
            stack[sp++] = node->right_child;
            stack[sp++] = node->left_or_first;
        }
    }

    return count;
}

/* ============================================================================
 * QUERY API
 * ============================================================================ */

static alea_spatial_index_t* get_spatial_index(alea_system_t* sys) {
    return sys ? sys->spatial_index : NULL;
}

typedef struct {
    alea_system_t* sys;
    alea_spatial_hit_t* restrict hits;
    size_t max_hits;
    size_t hit_count;
    bool terminal_only;
} query_ctx_t;

static void query_callback(const alea_cell_instance_t* restrict inst,
                           uint32_t inst_idx, void* userdata) {
    query_ctx_t* ctx = (query_ctx_t*)userdata;

    if (ctx->hit_count >= ctx->max_hits) return;
    if (ctx->terminal_only && !inst->is_terminal) return;

    const alea_cell_entry_t* cell = &ctx->sys->cells.data[inst->cell_index];

    alea_spatial_hit_t* restrict hit = &ctx->hits[ctx->hit_count++];
    hit->instance_index = inst_idx;
    hit->cell_index = inst->cell_index;
    hit->cell_id = cell->mc_cell_id;
    hit->material_id = cell->material_id;
    hit->universe_id = inst->universe_id;
    hit->depth = inst->depth;
    hit->is_terminal = inst->is_terminal;
    hit->transform = inst->transform;
}

int alea_spatial_query_region(alea_system_t* sys,
                             const alea_bbox_t* query_bbox,
                             alea_spatial_hit_t* out_hits,
                             size_t max_hits) {
    if (!sys || !query_bbox || !out_hits) return -1;

    alea_spatial_index_t* idx = get_spatial_index(sys);
    if (!idx || !idx->built) {
        if (alea_spatial_index_build(sys) != 0)
            return -1;
        idx = sys->spatial_index;
    }

    query_ctx_t ctx = {
        .sys = sys,
        .hits = out_hits,
        .max_hits = max_hits,
        .hit_count = 0,
        .terminal_only = true
    };

    alea_spatial_traverse(idx, query_bbox, query_callback, &ctx);

    return (int)ctx.hit_count;
}

int alea_spatial_query_slice_z(alea_system_t* sys,
                              double z,
                              double x_min, double x_max,
                              double y_min, double y_max,
                              alea_spatial_hit_t* out_hits,
                              size_t max_hits) {
    /* Convert slice to thin bbox (relative+absolute for large-coordinate robustness) */
    double ez = fabs(z) * EPSILON + EPSILON;
    alea_bbox_t query = {
        .min_x = x_min, .max_x = x_max,
        .min_y = y_min, .max_y = y_max,
        .min_z = z - ez, .max_z = z + ez
    };

    return alea_spatial_query_region(sys, &query, out_hits, max_hits);
}

/* Compare hits by depth, then by instance_index for consistent ordering.
 * Within the same depth, lower instance_index means it was encountered
 * earlier during universe traversal, which matches the recursive approach. */
static int compare_hits_by_depth_and_instance(const void* a, const void* b) {
    const alea_spatial_hit_t* ha = (const alea_spatial_hit_t*)a;
    const alea_spatial_hit_t* hb = (const alea_spatial_hit_t*)b;
    if (ha->depth != hb->depth) {
        return (ha->depth > hb->depth) - (ha->depth < hb->depth);
    }
    /* Same depth: sort by instance_index to match universe traversal order */
    return (ha->instance_index > hb->instance_index) - (ha->instance_index < hb->instance_index);
}

int alea_spatial_query_point(alea_system_t* sys,
                            double x, double y, double z,
                            alea_spatial_hit_t* out_hits,
                            size_t max_hits) {
    if (!sys || !out_hits) return -1;

    alea_spatial_index_t* idx = get_spatial_index(sys);
    if (!idx || !idx->built) {
        if (alea_spatial_index_build(sys) != 0)
            return -1;
        idx = sys->spatial_index;
    }

    /* Query with point bbox (relative+absolute for large-coordinate robustness) */
    double ex = fabs(x) * EPSILON + EPSILON;
    double ey = fabs(y) * EPSILON + EPSILON;
    double ez = fabs(z) * EPSILON + EPSILON;
    alea_bbox_t query = {
        .min_x = x - ex, .max_x = x + ex,
        .min_y = y - ey, .max_y = y + ey,
        .min_z = z - ez, .max_z = z + ez
    };

    query_ctx_t ctx = {
        .sys = sys,
        .hits = out_hits,
        .max_hits = max_hits,
        .hit_count = 0,
        .terminal_only = false
    };

    alea_spatial_traverse(idx, &query, query_callback, &ctx);

    /* Sort by depth so caller can process in hierarchy order */
    if (ctx.hit_count > 1) {
        qsort(out_hits, ctx.hit_count, sizeof(alea_spatial_hit_t), compare_hits_by_depth_and_instance);
    }

    return (int)ctx.hit_count;
}

/* ============================================================================
 * POINT-IN-CELL QUERY WITH COHERENCE CACHING
 * ============================================================================ */

/* Coherence cache: store last found cells to check first */
#define MAX_CACHE_DEPTH 8
typedef struct {
    uint32_t cell_index;
    alea_matrix_t transform;
    int cell_id;
    int material_id;
    int universe_id;
    int fill_universe;
    int depth;
    bool valid;
} cached_cell_t;

static _Thread_local cached_cell_t g_cell_cache[MAX_CACHE_DEPTH];
static _Thread_local int g_cache_count = 0;
static _Thread_local alea_system_t* g_cache_system = NULL;

/* Thread-local candidates buffer for point queries (freed via alea_spatial_reset_cache) */
static _Thread_local alea_spatial_hit_t* g_tls_candidates = NULL;
static _Thread_local size_t g_tls_candidates_cap = 0;

void alea_spatial_reset_cache(void) {
    g_cache_system = NULL;
    g_cache_count = 0;
    for (int i = 0; i < MAX_CACHE_DEPTH; i++) {
        g_cell_cache[i].valid = false;
    }
    free(g_tls_candidates);
    g_tls_candidates = NULL;
    g_tls_candidates_cap = 0;
}

int alea_spatial_find_cells_at_point(alea_system_t* sys,
                                    double x, double y, double z,
                                    alea_cell_hit_t* out_hits,
                                    size_t max_hits) {
    if (!sys || !out_hits || max_hits == 0) return -1;

    /* Check if cache is for a different system - invalidate if so */
    if (g_cache_system != sys) {
        g_cache_count = 0;
        g_cache_system = sys;
    }

    /* First, try the coherence cache - check if point is in same cells as last query */
    size_t hit_count = 0;
    bool cache_valid = true;

    for (int i = 0; i < g_cache_count && cache_valid && hit_count < max_hits; i++) {
        cached_cell_t* cached = &g_cell_cache[i];
        if (!cached->valid) {
            cache_valid = false;
            break;
        }

        /* Transform point to cell-local coordinates */
        double lx = x, ly = y, lz = z;
        alea_matrix_transform_point_inverse(&cached->transform, &lx, &ly, &lz);

        /* Quick containment test */
        const alea_cell_entry_t* cell = &sys->cells.data[cached->cell_index];
        if (!alea_contains_point(sys, cell->root_node_id, lx, ly, lz)) {
            cache_valid = false;
            break;
        }

        /* Cache hit - add to output */
        alea_cell_hit_t* hit = &out_hits[hit_count++];
        hit->cell_id = cached->cell_id;
        hit->cell_index = (int)cached->cell_index;
        hit->material_id = cached->material_id;
        hit->universe_id = cached->universe_id;
        hit->fill_universe = cached->fill_universe;
        hit->depth = cached->depth;
        hit->local_x = lx;
        hit->local_y = ly;
        hit->local_z = lz;
    }

    if (cache_valid && hit_count > 0) {
        return (int)hit_count;
    }

    /* Cache miss - do full BVH query */
    hit_count = 0;

    alea_spatial_index_t* idx = sys->spatial_index;
    size_t max_candidates = idx ? idx->instances.count : 256;
    if (max_candidates < 256) max_candidates = 256;

    if (max_candidates > g_tls_candidates_cap) {
        free(g_tls_candidates);
        g_tls_candidates = malloc(max_candidates * sizeof(alea_spatial_hit_t));
        g_tls_candidates_cap = g_tls_candidates ? max_candidates : 0;
    }
    alea_spatial_hit_t* candidates = g_tls_candidates;
    if (!candidates) return -1;

    int n = alea_spatial_query_point(sys, x, y, z, candidates, max_candidates);
    if (n < 0) {
        return -1;
    }
    if (n == 0) {
        g_cache_count = 0;
        return 0;
    }

    /* Re-read idx: alea_spatial_query_point may have built the index */
    idx = sys->spatial_index;

    /* Verify containment and build hierarchy.
     * Key insight: we must track the hierarchy path. At each depth,
     * we find the containing cell and its fill_universe tells us
     * which universe to search at the next depth.
     *
     * IMPORTANT: When multiple cells fill the same universe with different
     * transforms, we must only consider instances whose parent matches the
     * cell we found at the previous depth.
     *
     * NOTE: correctness depends on cell bboxes being conservative (over-
     * approximating).  If a FILL cell's bbox is too tight, its children
     * will be silently missed because expected_universe never advances. */
    int expected_universe = 0;  /* Start at universe 0 */
    int last_depth = -1;
    uint32_t last_cell_idx = UINT32_MAX;  /* Cell found at previous depth */
    g_cache_count = 0;

    for (int i = 0; i < n && hit_count < max_hits; i++) {
        alea_spatial_hit_t* cand = &candidates[i];

        /* At each new depth, we must find a cell in the expected universe */
        if (cand->depth != last_depth + 1 && cand->depth != 0) {
            /* Skip - we're looking for the next depth level */
            if (cand->depth <= last_depth) continue;
            /* If we jumped depths, we missed a level - stop */
            if (cand->depth > last_depth + 1 && last_depth >= 0) continue;
        }

        /* Skip cells not in the expected universe for this depth */
        if (cand->universe_id != expected_universe) continue;

        /* For depth > 0, check that this instance's parent matches the cell
         * we found at the previous level. This ensures we follow the correct
         * FILL chain when multiple cells fill the same universe. */
        if (cand->depth > 0 && last_cell_idx != UINT32_MAX) {
            alea_cell_instance_t* inst = &idx->instances.data[cand->instance_index];
            if (inst->parent_cell_index != last_cell_idx) {
                continue;  /* Wrong fill chain */
            }
        }

        /* Transform point to cell-local coordinates */
        double lx = x, ly = y, lz = z;
        alea_matrix_transform_point_inverse(&cand->transform, &lx, &ly, &lz);

        const alea_cell_entry_t* cell = &sys->cells.data[cand->cell_index];

#if ALEA_USE_INTERVAL_CULLING
        /* Quick interval-based pre-filter: if the point is definitely outside
         * the cell's CSG tree (based on interval arithmetic), skip expensive
         * containment test. Use a point-sized bbox for the query. */
        {
            alea_bbox_t point_box = { lx, lx, ly, ly, lz, lz };
            alea_box_relation_t rel = alea_tree_box_relation(sys, cell->root_node_id, &point_box);
            if (rel == ALEA_RELATION_POSITIVE) {
                g_interval_culls++;
                continue;  /* Definitely outside - skip */
            }
        }
#endif

        if (!alea_contains_point(sys, cell->root_node_id, lx, ly, lz)) {
            continue;
        }

        /* Found a cell at this depth level */
        alea_cell_hit_t* hit = &out_hits[hit_count++];
        hit->cell_id = cand->cell_id;
        hit->cell_index = (int)cand->cell_index;
        hit->material_id = cand->material_id;
        hit->universe_id = cand->universe_id;
        hit->fill_universe = cell->fill_universe;
        hit->depth = cand->depth;
        hit->local_x = lx;
        hit->local_y = ly;
        hit->local_z = lz;

        /* Update cache */
        if (g_cache_count < MAX_CACHE_DEPTH) {
            cached_cell_t* cached = &g_cell_cache[g_cache_count++];
            cached->cell_index = cand->cell_index;
            cached->transform = cand->transform;
            cached->cell_id = cand->cell_id;
            cached->material_id = cand->material_id;
            cached->universe_id = cand->universe_id;
            cached->fill_universe = cell->fill_universe;
            cached->depth = cand->depth;
            cached->valid = true;
        }

        last_depth = cand->depth;
        last_cell_idx = cand->cell_index;

        /* If this cell has a FILL, the next depth must be in that universe.
         * Don't break on terminal cells - continue looking for overlapping
         * cells at the same depth (like the recursive algorithm does). */
        if (cell->fill_universe > 0) {
            expected_universe = cell->fill_universe;
        }
        /* Continue loop to find other overlapping cells at the same depth */
    }

    return (int)hit_count;
}

int alea_spatial_find_cells_batch(alea_system_t* sys,
                                  const double* points,
                                  size_t n_points,
                                  alea_cell_hit_t* out_hits,
                                  size_t max_hits_per_point,
                                  int* out_counts) {
    if (!sys || !points || !out_hits || !out_counts || max_hits_per_point == 0)
        return -1;

    /* Ensure spatial index is built */
    if (!sys->spatial_index || !sys->spatial_index->built) {
        if (alea_spatial_index_build(sys) != 0)
            return -1;
    }

    for (size_t i = 0; i < n_points; i++) {
        double x = points[i * 3 + 0];
        double y = points[i * 3 + 1];
        double z = points[i * 3 + 2];
        alea_cell_hit_t* dst = &out_hits[i * max_hits_per_point];
        int n = alea_spatial_find_cells_at_point(sys, x, y, z, dst, max_hits_per_point);
        out_counts[i] = (n >= 0) ? n : 0;
    }

    return 0;
}
