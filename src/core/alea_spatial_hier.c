// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "core/alea_spatial_hier.h"
#include "core/alea_system.h"
#include "core/alea_surface.h"
#include "primitives/bbox.h"
#include "util/alea_log.h"
#include <float.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#define HIER_BVH_LEAF_SIZE 4
#define HIER_DEFAULT_BLAS_THRESHOLD 1
#define HIER_MAX_PLACEMENT_DEPTH 64

enum {
    HIER_PLACEMENT_ROOT = 1u << 0,
    HIER_PLACEMENT_FILL = 1u << 1,
    HIER_PLACEMENT_LATTICE = 1u << 2
};

typedef struct {
    alea_bbox_t bbox;
    uint32_t left_or_first;
    uint32_t right_child;
    uint16_t count;
    uint8_t axis;
    uint8_t pad;
} hier_bvh_node_t;

typedef struct {
    uint32_t cell_index;
    alea_bbox_t bbox;
} hier_blas_cell_t;

typedef struct {
    int universe_id;
    size_t cell_count;
    hier_blas_cell_t* cells;
    hier_bvh_node_t* nodes;
    uint32_t* indices;
    size_t node_count;
    size_t node_capacity;
    alea_bbox_t bounds;
    int built;
} hier_universe_blas_t;

typedef struct {
    uint32_t id;
    uint32_t parent_id;
    uint32_t parent_cell_index;
    int universe_id;
    int depth;
    uint32_t transform_index;
    uint32_t flags;
    alea_bbox_t local_bbox;
    alea_bbox_t world_bbox;
} hier_placement_t;

struct alea_hier_spatial_index {
    hier_universe_blas_t* blas;
    size_t blas_count;
    size_t blas_capacity;
    hier_placement_t* placements;
    size_t placement_count;
    size_t placement_capacity;
    alea_matrix_t* transforms;
    size_t transform_count;
    size_t transform_capacity;
    alea_hier_spatial_stats_t stats;
    int built;
};

typedef struct {
    uint32_t index;
    double centroid[3];
} hier_bvh_item_t;

static double monotonic_seconds(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec * 1e-6;
}

static double bytes_to_mib(size_t bytes) {
    return (double)bytes / (1024.0 * 1024.0);
}

static size_t blas_threshold(void) {
    const char* env = getenv("ALEA_HIER_BLAS_THRESHOLD");
    if (env && env[0]) {
        char* end = NULL;
        unsigned long value = strtoul(env, &end, 10);
        if (end != env && value > 0) return (size_t)value;
    }
    return HIER_DEFAULT_BLAS_THRESHOLD;
}

static alea_bbox_t local_cell_bbox(alea_system_t* sys, uint32_t cell_index) {
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

static alea_bbox_t bbox_transform(const alea_bbox_t* bbox, const alea_matrix_t* mat) {
    const double bmin[3] = {bbox->min_x, bbox->min_y, bbox->min_z};
    const double bmax[3] = {bbox->max_x, bbox->max_y, bbox->max_z};
    double rmin[3], rmax[3];

    for (int i = 0; i < 3; i++) {
        rmin[i] = rmax[i] = mat->m[i * 4 + 3];
        for (int j = 0; j < 3; j++) {
            double e = mat->m[i * 4 + j] * bmin[j];
            double f = mat->m[i * 4 + j] * bmax[j];
            if (e < f) {
                rmin[i] += e;
                rmax[i] += f;
            } else {
                rmin[i] += f;
                rmax[i] += e;
            }
        }
    }

    return (alea_bbox_t){rmin[0], rmax[0], rmin[1], rmax[1], rmin[2], rmax[2]};
}

static size_t estimate_bvh_node_count(size_t item_count) {
    if (item_count == 0) return 0;
    if (item_count <= HIER_BVH_LEAF_SIZE) return 1;

    size_t left = item_count / 2;
    size_t right = item_count - left;
    return 1 + estimate_bvh_node_count(left) + estimate_bvh_node_count(right);
}

static int ensure_node_capacity(hier_universe_blas_t* blas, size_t needed) {
    size_t min_cap = blas->node_count + needed;
    if (min_cap <= blas->node_capacity) return 0;

    size_t new_cap = blas->node_capacity ? blas->node_capacity * 2 : 16;
    while (new_cap < min_cap) new_cap *= 2;

    hier_bvh_node_t* nodes = realloc(blas->nodes, new_cap * sizeof(*nodes));
    if (!nodes) return -1;

    blas->nodes = nodes;
    blas->node_capacity = new_cap;
    return 0;
}

static alea_bbox_t item_bbox(const hier_universe_blas_t* blas,
                             const hier_bvh_item_t* items,
                             size_t start,
                             size_t end) {
    alea_bbox_t result = alea_bbox_empty();
    for (size_t i = start; i < end; i++) {
        result = alea_bbox_union(&result, &blas->cells[items[i].index].bbox);
    }
    return result;
}

static int best_split_axis(const hier_bvh_item_t* items,
                           size_t start,
                           size_t end) {
    double minv[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
    double maxv[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};

    for (size_t i = start; i < end; i++) {
        for (int axis = 0; axis < 3; axis++) {
            if (items[i].centroid[axis] < minv[axis]) minv[axis] = items[i].centroid[axis];
            if (items[i].centroid[axis] > maxv[axis]) maxv[axis] = items[i].centroid[axis];
        }
    }

    double dx = maxv[0] - minv[0];
    double dy = maxv[1] - minv[1];
    double dz = maxv[2] - minv[2];

    if (dx >= dy && dx >= dz) return 0;
    if (dy >= dz) return 1;
    return 2;
}

static void quickselect(hier_bvh_item_t* items,
                        size_t lo,
                        size_t hi,
                        size_t target,
                        int axis) {
    while (lo + 1 < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (items[mid].centroid[axis] < items[lo].centroid[axis]) {
            hier_bvh_item_t tmp = items[lo]; items[lo] = items[mid]; items[mid] = tmp;
        }
        if (items[hi - 1].centroid[axis] < items[lo].centroid[axis]) {
            hier_bvh_item_t tmp = items[lo]; items[lo] = items[hi - 1]; items[hi - 1] = tmp;
        }
        if (items[mid].centroid[axis] < items[hi - 1].centroid[axis]) {
            hier_bvh_item_t tmp = items[mid]; items[mid] = items[hi - 1]; items[hi - 1] = tmp;
        }

        double pivot = items[hi - 1].centroid[axis];
        size_t i = lo;
        size_t j = hi - 2;

        for (;;) {
            while (i < hi && items[i].centroid[axis] < pivot) i++;
            while (j > lo && items[j].centroid[axis] > pivot) j--;
            if (i >= j) break;
            hier_bvh_item_t tmp = items[i]; items[i] = items[j]; items[j] = tmp;
            i++;
            j--;
        }

        hier_bvh_item_t tmp = items[i]; items[i] = items[hi - 1]; items[hi - 1] = tmp;

        if (i == target) return;
        if (target < i) hi = i;
        else lo = i + 1;
    }
}

static uint32_t build_bvh_recursive(hier_universe_blas_t* blas,
                                    hier_bvh_item_t* items,
                                    size_t start,
                                    size_t end,
                                    int depth) {
    if (ensure_node_capacity(blas, 1) != 0) return UINT32_MAX;

    uint32_t node_index = (uint32_t)blas->node_count++;
    hier_bvh_node_t* node = &blas->nodes[node_index];
    size_t count = end - start;

    if (count <= HIER_BVH_LEAF_SIZE || depth > 30) {
        node->bbox = item_bbox(blas, items, start, end);
        node->left_or_first = (uint32_t)start;
        node->right_child = 0;
        node->count = (uint16_t)count;
        node->axis = 0;
        node->pad = 0;
        return node_index;
    }

    int axis = best_split_axis(items, start, end);
    size_t mid = (start + end) / 2;
    if (end - start > 2) {
        quickselect(items, start, end, mid, axis);
    }
    if (mid == start || mid == end) mid = (start + end) / 2;

    node->axis = (uint8_t)axis;
    node->count = 0;
    node->pad = 0;

    uint32_t left = build_bvh_recursive(blas, items, start, mid, depth + 1);
    uint32_t right = build_bvh_recursive(blas, items, mid, end, depth + 1);
    if (left == UINT32_MAX || right == UINT32_MAX) return UINT32_MAX;

    node = &blas->nodes[node_index];
    node->bbox = alea_bbox_union(&blas->nodes[left].bbox, &blas->nodes[right].bbox);
    node->left_or_first = left;
    node->right_child = right;

    return node_index;
}

static int build_universe_blas(alea_system_t* sys,
                               const alea_universe_t* univ,
                               hier_universe_blas_t* blas) {
    memset(blas, 0, sizeof(*blas));
    blas->universe_id = univ->universe_id;
    blas->cell_count = univ->cell_indices.count;
    blas->bounds = alea_bbox_empty();

    if (blas->cell_count == 0) {
        blas->built = 1;
        return 0;
    }

    blas->cells = calloc(blas->cell_count, sizeof(*blas->cells));
    hier_bvh_item_t* items = calloc(blas->cell_count, sizeof(*items));
    blas->indices = calloc(blas->cell_count, sizeof(*blas->indices));
    if (!blas->cells || !items || !blas->indices) {
        free(items);
        return -1;
    }

    for (size_t i = 0; i < blas->cell_count; i++) {
        uint32_t cell_index = (uint32_t)univ->cell_indices.data[i];
        alea_bbox_t bbox = local_cell_bbox(sys, cell_index);
        blas->cells[i].cell_index = cell_index;
        blas->cells[i].bbox = bbox;
        blas->bounds = alea_bbox_union(&blas->bounds, &bbox);

        items[i].index = (uint32_t)i;
        items[i].centroid[0] = 0.5 * (bbox.min_x + bbox.max_x);
        items[i].centroid[1] = 0.5 * (bbox.min_y + bbox.max_y);
        items[i].centroid[2] = 0.5 * (bbox.min_z + bbox.max_z);
    }

    size_t estimated_nodes = estimate_bvh_node_count(blas->cell_count);
    blas->nodes = calloc(estimated_nodes ? estimated_nodes : 1, sizeof(*blas->nodes));
    if (!blas->nodes) {
        free(items);
        return -1;
    }
    blas->node_capacity = estimated_nodes;

    uint32_t root = build_bvh_recursive(blas, items, 0, blas->cell_count, 0);
    if (root == UINT32_MAX) {
        free(items);
        return -1;
    }

    for (size_t i = 0; i < blas->cell_count; i++) {
        blas->indices[i] = items[i].index;
    }

    free(items);
    blas->built = 1;
    return 0;
}

static void free_universe_blas(hier_universe_blas_t* blas) {
    if (!blas) return;
    free(blas->cells);
    free(blas->nodes);
    free(blas->indices);
    memset(blas, 0, sizeof(*blas));
}

void alea_hier_spatial_index_free(alea_hier_spatial_index_t* idx) {
    if (!idx) return;
    for (size_t i = 0; i < idx->blas_count; i++) {
        free_universe_blas(&idx->blas[i]);
    }
    free(idx->blas);
    free(idx->placements);
    free(idx->transforms);
    free(idx);
}

const alea_hier_spatial_stats_t*
alea_hier_spatial_index_stats(const alea_hier_spatial_index_t* idx) {
    return idx ? &idx->stats : NULL;
}

static int ensure_blas_capacity(alea_hier_spatial_index_t* idx, size_t needed) {
    if (needed <= idx->blas_capacity) return 0;
    size_t new_cap = idx->blas_capacity ? idx->blas_capacity * 2 : 16;
    while (new_cap < needed) new_cap *= 2;

    hier_universe_blas_t* blas = realloc(idx->blas, new_cap * sizeof(*blas));
    if (!blas) return -1;

    idx->blas = blas;
    idx->blas_capacity = new_cap;
    return 0;
}

static int ensure_placement_capacity(alea_hier_spatial_index_t* idx, size_t needed) {
    if (needed <= idx->placement_capacity) return 0;
    size_t new_cap = idx->placement_capacity ? idx->placement_capacity * 2 : 64;
    while (new_cap < needed) new_cap *= 2;

    hier_placement_t* placements = realloc(idx->placements, new_cap * sizeof(*placements));
    if (!placements) return -1;

    idx->placements = placements;
    idx->placement_capacity = new_cap;
    return 0;
}

static int ensure_transform_capacity(alea_hier_spatial_index_t* idx, size_t needed) {
    if (needed <= idx->transform_capacity) return 0;
    size_t new_cap = idx->transform_capacity ? idx->transform_capacity * 2 : 64;
    while (new_cap < needed) new_cap *= 2;

    alea_matrix_t* transforms = realloc(idx->transforms, new_cap * sizeof(*transforms));
    if (!transforms) return -1;

    idx->transforms = transforms;
    idx->transform_capacity = new_cap;
    return 0;
}

static int append_transform(alea_hier_spatial_index_t* idx,
                            const alea_matrix_t* transform,
                            uint32_t* out_index) {
    if (ensure_transform_capacity(idx, idx->transform_count + 1) != 0) {
        return -1;
    }

    *out_index = (uint32_t)idx->transform_count;
    idx->transforms[idx->transform_count++] = *transform;
    idx->stats.transform_count = idx->transform_count;
    return 0;
}

static int append_placement(alea_hier_spatial_index_t* idx,
                            uint32_t parent_id,
                            uint32_t parent_cell_index,
                            int universe_id,
                            int depth,
                            uint32_t flags,
                            const alea_bbox_t* local_bbox,
                            const alea_bbox_t* world_bbox,
                            const alea_matrix_t* transform,
                            uint32_t* out_id) {
    if (ensure_placement_capacity(idx, idx->placement_count + 1) != 0) {
        return -1;
    }

    uint32_t transform_index = UINT32_MAX;
    if (transform) {
        if (append_transform(idx, transform, &transform_index) != 0) {
            return -1;
        }
    }

    uint32_t id = (uint32_t)idx->placement_count;
    hier_placement_t* placement = &idx->placements[idx->placement_count++];
    placement->id = id;
    placement->parent_id = parent_id;
    placement->parent_cell_index = parent_cell_index;
    placement->universe_id = universe_id;
    placement->depth = depth;
    placement->transform_index = transform_index;
    placement->flags = flags;
    placement->local_bbox = local_bbox ? *local_bbox : alea_bbox_empty();
    placement->world_bbox = world_bbox ? *world_bbox : alea_bbox_empty();

    idx->stats.placement_count = idx->placement_count;
    if (flags & HIER_PLACEMENT_ROOT) idx->stats.root_placement_count++;
    if (flags & HIER_PLACEMENT_FILL) idx->stats.fill_placement_count++;
    if (flags & HIER_PLACEMENT_LATTICE) idx->stats.lattice_placement_count++;
    if (depth > idx->stats.max_placement_depth) idx->stats.max_placement_depth = depth;

    if (out_id) *out_id = id;
    return 0;
}

static int fill_transform_matrix(alea_system_t* sys,
                                 const alea_cell_entry_t* cell,
                                 alea_matrix_t* out) {
    if (cell->fill_transform <= 0) {
        alea_matrix_identity(out);
        return 0;
    }

    const alea_transform_t* tr = alea_get_transform(sys, cell->fill_transform);
    if (!tr) return -1;

    if (!alea_matrix_from_mcnp(out, tr->cosines, tr->value_count, false)) {
        return -1;
    }

    return 0;
}

static int collect_placements_recursive(alea_system_t* sys,
                                        alea_hier_spatial_index_t* idx,
                                        int universe_id,
                                        uint32_t parent_placement_id,
                                        const alea_matrix_t* accumulated,
                                        int depth) {
    if (depth >= HIER_MAX_PLACEMENT_DEPTH) return 0;
    if (g_alea_interrupted) return -1;

    const alea_universe_t* univ = alea_get_universe(sys, universe_id);
    if (!univ) return 0;

    for (size_t i = 0; i < univ->cell_indices.count; i++) {
        uint32_t cell_index = (uint32_t)univ->cell_indices.data[i];
        const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
        alea_bbox_t local_bbox = local_cell_bbox(sys, cell_index);
        alea_bbox_t world_bbox = accumulated ? bbox_transform(&local_bbox, accumulated) : local_bbox;

        if (cell->lat_type != 0 && cell->lat_fill) {
            idx->stats.lattice_cell_count++;
            if (append_placement(idx, parent_placement_id, cell_index,
                                 universe_id, depth + 1,
                                 HIER_PLACEMENT_LATTICE,
                                 &local_bbox, &world_bbox,
                                 accumulated, NULL) != 0) {
                return -1;
            }
            continue;
        }

        if (cell->fill_universe > 0) {
            alea_matrix_t fill_mat;
            if (fill_transform_matrix(sys, cell, &fill_mat) != 0) {
                return -1;
            }

            alea_matrix_t child_transform;
            if (accumulated) {
                alea_matrix_multiply(&child_transform, accumulated, &fill_mat);
            } else {
                child_transform = fill_mat;
            }
            if (!alea_matrix_invert(&child_transform)) {
                return -1;
            }

            uint32_t child_id = UINT32_MAX;
            idx->stats.fill_cell_count++;
            if (append_placement(idx, parent_placement_id, cell_index,
                                 cell->fill_universe, depth + 1,
                                 HIER_PLACEMENT_FILL,
                                 &local_bbox, &world_bbox,
                                 &child_transform, &child_id) != 0) {
                return -1;
            }

            if (collect_placements_recursive(sys, idx, cell->fill_universe,
                                             child_id, &child_transform,
                                             depth + 1) != 0) {
                return -1;
            }
        }
    }
    return 0;
}

static int collect_placements(alea_system_t* sys, alea_hier_spatial_index_t* idx) {
    alea_matrix_t identity;
    alea_matrix_identity(&identity);

    uint32_t root_id = UINT32_MAX;
    alea_bbox_t root_bbox = alea_bbox_empty();
    const alea_universe_t* root = alea_get_universe(sys, 0);
    if (root) root_bbox = root->bbox;

    if (append_placement(idx, UINT32_MAX, UINT32_MAX, 0, 0,
                         HIER_PLACEMENT_ROOT, &root_bbox, &root_bbox,
                         &identity, &root_id) != 0) {
        return -1;
    }

    return collect_placements_recursive(sys, idx, 0, root_id, &identity, 0);
}

int alea_hier_spatial_index_build(alea_system_t* sys) {
    if (!sys) return -1;

    double t_start = monotonic_seconds();
    size_t threshold = blas_threshold();

    if (!sys->universe_index_built) {
        if (alea_build_universe_index(sys) != 0) return -1;
    }

    if (sys->hier_spatial_index) {
        alea_hier_spatial_index_free(sys->hier_spatial_index);
        sys->hier_spatial_index = NULL;
    }

    alea_hier_spatial_index_t* idx = calloc(1, sizeof(*idx));
    if (!idx) return -1;

    idx->stats.universe_count = sys->universes.count;

    if (collect_placements(sys, idx) != 0) {
        alea_hier_spatial_index_free(idx);
        return -1;
    }

    for (size_t i = 0; i < sys->universes.count; i++) {
        const alea_universe_t* univ = &sys->universes.data[i];
        size_t cell_count = univ->cell_indices.count;

        if ((int)cell_count > idx->stats.max_universe_cells) {
            idx->stats.max_universe_cells = (int)cell_count;
            idx->stats.largest_universe_id = univ->universe_id;
        }

        if (cell_count < threshold || cell_count == 0) {
            idx->stats.linear_universe_count++;
            continue;
        }

        if (ensure_blas_capacity(idx, idx->blas_count + 1) != 0) {
            alea_hier_spatial_index_free(idx);
            return -1;
        }

        hier_universe_blas_t* blas = &idx->blas[idx->blas_count];
        if (build_universe_blas(sys, univ, blas) != 0) {
            free_universe_blas(blas);
            alea_hier_spatial_index_free(idx);
            return -1;
        }
        idx->blas_count++;

        idx->stats.blas_count++;
        idx->stats.blas_cell_count += blas->cell_count;
        idx->stats.blas_node_count += blas->node_count;
        idx->stats.memory_bytes += blas->cell_count * sizeof(*blas->cells);
        idx->stats.memory_bytes += blas->cell_count * sizeof(*blas->indices);
        idx->stats.memory_bytes += blas->node_capacity * sizeof(*blas->nodes);
    }

    idx->stats.memory_bytes += idx->blas_capacity * sizeof(*idx->blas);
    idx->stats.memory_bytes += idx->placement_capacity * sizeof(*idx->placements);
    idx->stats.memory_bytes += idx->transform_capacity * sizeof(*idx->transforms);
    idx->built = 1;
    sys->hier_spatial_index = idx;

    double t_end = monotonic_seconds();
    ALEA_LOG_INFO("Hier spatial BLAS build: universes=%zu blas=%zu linear=%zu cells=%zu nodes=%zu memory=%.1f MiB time=%.3fs",
                  idx->stats.universe_count, idx->stats.blas_count,
                  idx->stats.linear_universe_count, idx->stats.blas_cell_count,
                  idx->stats.blas_node_count,
                  bytes_to_mib(idx->stats.memory_bytes), t_end - t_start);
    ALEA_LOG_INFO("Hier spatial placements: total=%zu root=%zu fill=%zu lattice=%zu max_depth=%d transforms=%zu largest_universe=%d/%d cells",
                  idx->stats.placement_count, idx->stats.root_placement_count,
                  idx->stats.fill_placement_count, idx->stats.lattice_placement_count,
                  idx->stats.max_placement_depth, idx->stats.transform_count,
                  idx->stats.largest_universe_id, idx->stats.max_universe_cells);

    return 0;
}
