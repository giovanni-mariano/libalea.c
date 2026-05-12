// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "core/alea_spatial_hier.h"
#include "core/alea_eval.h"
#include "core/alea_system.h"
#include "core/alea_surface.h"
#include "primitives/bbox.h"
#include "util/alea_log.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#define HIER_BVH_LEAF_SIZE 4
#define HIER_DEFAULT_BLAS_THRESHOLD 1
#define HIER_MAX_PLACEMENT_DEPTH 64
#ifndef M_SQRT3
#define M_SQRT3 1.73205080756887729353
#endif

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

typedef struct {
    alea_system_t* sys;
    alea_hier_spatial_index_t* idx;
    const hier_universe_blas_t* blas;
    double x;
    double y;
    double z;
    alea_size_vec_t* candidates;
    int error;
} hier_point_query_ctx_t;

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

static alea_bbox_t lattice_container_bbox(const alea_cell_entry_t* cell) {
    if (cell->lat_type == 1) {
        int ni = cell->lat_fill_dims[1] - cell->lat_fill_dims[0] + 1;
        int nj = cell->lat_fill_dims[3] - cell->lat_fill_dims[2] + 1;
        int nk = cell->lat_fill_dims[5] - cell->lat_fill_dims[4] + 1;
        return (alea_bbox_t){
            cell->lat_lower_left[0],
            cell->lat_lower_left[0] + (double)ni * cell->lat_pitch[0],
            cell->lat_lower_left[1],
            cell->lat_lower_left[1] + (double)nj * cell->lat_pitch[1],
            cell->lat_lower_left[2],
            cell->lat_lower_left[2] + (double)nk * cell->lat_pitch[2]
        };
    }

    /* Hex lattice bounds are intentionally conservative. The exact element
     * lookup still decides whether the point maps to a filled universe. */
    double p = cell->lat_pitch[0] > 0.0 ? cell->lat_pitch[0] : 1.0;
    int ni = cell->lat_fill_dims[1] - cell->lat_fill_dims[0] + 1;
    int nj = cell->lat_fill_dims[3] - cell->lat_fill_dims[2] + 1;
    int nk = cell->lat_fill_dims[5] - cell->lat_fill_dims[4] + 1;
    double extent = (double)(ni + nj + 2) * p;
    double zmax = (nk == 1) ? p : cell->lat_lower_left[2] + (double)nk * cell->lat_pitch[2];
    return (alea_bbox_t){-extent, extent, -extent, extent,
                         cell->lat_lower_left[2], zmax};
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

static alea_bbox_t bbox_transform_inverse_conservative(const alea_bbox_t* bbox,
                                                       const alea_matrix_t* mat) {
    if (!mat || !mat->has_inverse) return alea_bbox_empty();

    alea_bbox_t result = alea_bbox_empty();
    for (int ix = 0; ix < 2; ix++) {
        for (int iy = 0; iy < 2; iy++) {
            for (int iz = 0; iz < 2; iz++) {
                double x = ix ? bbox->max_x : bbox->min_x;
                double y = iy ? bbox->max_y : bbox->min_y;
                double z = iz ? bbox->max_z : bbox->min_z;
                alea_matrix_transform_point_inverse(mat, &x, &y, &z);
                alea_bbox_t p = {x, x, y, y, z, z};
                result = alea_bbox_union(&result, &p);
            }
        }
    }
    return result;
}

static int bbox_intersects_local(const alea_bbox_t* a, const alea_bbox_t* b) {
    return a->min_x <= b->max_x && a->max_x >= b->min_x &&
           a->min_y <= b->max_y && a->max_y >= b->min_y &&
           a->min_z <= b->max_z && a->max_z >= b->min_z;
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

static int compare_size_t_ascending(const void* a, const void* b) {
    const size_t* ia = (const size_t*)a;
    const size_t* ib = (const size_t*)b;
    return (*ia > *ib) - (*ia < *ib);
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
        const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
        alea_bbox_t bbox = (cell->lat_type != 0 && cell->lat_fill)
            ? lattice_container_bbox(cell)
            : local_cell_bbox(sys, cell_index);
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

static int lattice_rect_lookup(const alea_cell_entry_t* cell,
                               double px, double py, double pz,
                               double* ox, double* oy, double* oz) {
    int ni = cell->lat_fill_dims[1] - cell->lat_fill_dims[0] + 1;
    int nj = cell->lat_fill_dims[3] - cell->lat_fill_dims[2] + 1;
    int nk = cell->lat_fill_dims[5] - cell->lat_fill_dims[4] + 1;

    int i = (int)floor((px - cell->lat_lower_left[0]) / cell->lat_pitch[0]);
    int j = (nj == 1) ? 0 : (int)floor((py - cell->lat_lower_left[1]) / cell->lat_pitch[1]);
    int k = (nk == 1) ? 0 : (int)floor((pz - cell->lat_lower_left[2]) / cell->lat_pitch[2]);

    if (i < 0 || i >= ni || j < 0 || j >= nj || k < 0 || k >= nk) {
        return -1;
    }

    *ox = cell->lat_lower_left[0] + (i + 0.5) * cell->lat_pitch[0];
    *oy = cell->lat_lower_left[1] + (j + 0.5) * cell->lat_pitch[1];
    *oz = cell->lat_lower_left[2] + (k + 0.5) * cell->lat_pitch[2];

    size_t idx = (size_t)(i * nj * nk + j * nk + k);
    if (idx >= cell->lat_fill_count) return -1;

    return cell->lat_fill[idx];
}

static int lattice_hex_lookup_local(const alea_cell_entry_t* cell,
                                    double px, double py, double pz,
                                    double* ox, double* oy, double* oz) {
    double p = cell->lat_pitch[0];
    if (p <= 0.0) return -1;

    double fj = py / (p * M_SQRT3 * 0.5);
    double fi = px / p - fj * 0.5;

    double cx = fi;
    double cz = fj;
    double cy = -fi - fj;

    int ri = (int)round(cx);
    int rj = (int)round(cy);
    int rk = (int)round(cz);

    double dx = fabs(ri - cx);
    double dy = fabs(rj - cy);
    double dz = fabs(rk - cz);

    if (dx > dy && dx > dz) {
        ri = -rj - rk;
    } else if (dy > dz) {
        rj = -ri - rk;
    } else {
        rk = -ri - rj;
    }

    int ni = cell->lat_fill_dims[1] - cell->lat_fill_dims[0] + 1;
    int nj = cell->lat_fill_dims[3] - cell->lat_fill_dims[2] + 1;
    int nk = cell->lat_fill_dims[5] - cell->lat_fill_dims[4] + 1;

    int oi = ri - cell->lat_fill_dims[0];
    int oj = rk - cell->lat_fill_dims[2];
    int ok = (nk == 1) ? 0
           : (int)floor((pz - cell->lat_lower_left[2]) / cell->lat_pitch[2]);

    if (oi < 0 || oi >= ni || oj < 0 || oj >= nj || ok < 0 || ok >= nk) {
        return -1;
    }

    *ox = ri * p + rk * p * 0.5;
    *oy = rk * p * M_SQRT3 * 0.5;
    *oz = (nk == 1) ? 0.0 : cell->lat_lower_left[2] + (ok + 0.5) * cell->lat_pitch[2];

    size_t idx = (size_t)(oi * nj * nk + oj * nk + ok);
    if (idx >= cell->lat_fill_count) return -1;

    return cell->lat_fill[idx];
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
        alea_bbox_t local_bbox = (cell->lat_type != 0 && cell->lat_fill)
            ? lattice_container_bbox(cell)
            : local_cell_bbox(sys, cell_index);
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

static const hier_universe_blas_t* find_blas(const alea_hier_spatial_index_t* idx,
                                             int universe_id) {
    if (!idx) return NULL;
    for (size_t i = 0; i < idx->blas_count; i++) {
        if (idx->blas[i].universe_id == universe_id) {
            return &idx->blas[i];
        }
    }
    return NULL;
}

static const alea_matrix_t* placement_transform(const alea_hier_spatial_index_t* idx,
                                                const hier_placement_t* placement) {
    if (!idx || !placement) return NULL;
    if (placement->transform_index >= idx->transform_count) return NULL;
    return &idx->transforms[placement->transform_index];
}

static void query_blas_node(hier_point_query_ctx_t* ctx, uint32_t node_idx) {
    if (ctx->error) return;
    if (node_idx >= ctx->blas->node_count) {
        ctx->error = -1;
        return;
    }

    const hier_bvh_node_t* node = &ctx->blas->nodes[node_idx];
    ctx->idx->stats.point_blas_node_visits++;
    ctx->idx->stats.point_bbox_tests++;
    if (!alea_bbox_contains_point(&node->bbox, ctx->x, ctx->y, ctx->z)) {
        return;
    }

    if (node->count > 0) {
        for (uint16_t i = 0; i < node->count; i++) {
            uint32_t cell_pos = ctx->blas->indices[node->left_or_first + i];
            if (cell_pos >= ctx->blas->cell_count) {
                ctx->error = -1;
                return;
            }
            const hier_blas_cell_t* blas_cell = &ctx->blas->cells[cell_pos];
            ctx->idx->stats.point_bbox_tests++;
            if (!alea_bbox_contains_point(&blas_cell->bbox, ctx->x, ctx->y, ctx->z)) {
                continue;
            }
            ctx->idx->stats.point_candidates++;
            if (alea_vec_push(ctx->candidates, (size_t)cell_pos, size_t) != 0) {
                ctx->error = -1;
                return;
            }
        }
        return;
    }

    query_blas_node(ctx, node->left_or_first);
    query_blas_node(ctx, node->right_child);
}

static int hier_query_universe(alea_system_t* sys,
                               alea_hier_spatial_index_t* idx,
                               int universe_id,
                               double lx,
                               double ly,
                               double lz,
                               int depth,
                               alea_cell_hit_t* out_hits,
                               size_t max_hits,
                               size_t* hit_count);

static int process_point_cell(alea_system_t* sys,
                              alea_hier_spatial_index_t* idx,
                              uint32_t cell_index,
                              double lx,
                              double ly,
                              double lz,
                              int depth,
                              alea_cell_hit_t* out_hits,
                              size_t max_hits,
                              size_t* hit_count) {
    const alea_cell_entry_t* cell = &sys->cells.data[cell_index];

    if (cell->lat_type != 0 && cell->lat_fill) {
        double ox, oy, oz;
        int fill_univ = (cell->lat_type == 2)
            ? lattice_hex_lookup_local(cell, lx, ly, lz, &ox, &oy, &oz)
            : lattice_rect_lookup(cell, lx, ly, lz, &ox, &oy, &oz);
        if (fill_univ < 0) return 0;

        if (*hit_count < max_hits) {
            alea_cell_hit_t* hit = &out_hits[*hit_count];
            hit->cell_id = cell->mc_cell_id;
            hit->cell_index = (int)cell_index;
            hit->material_id = cell->material_id;
            hit->universe_id = cell->universe_id;
            hit->fill_universe = fill_univ;
            hit->depth = depth;
            hit->local_x = lx;
            hit->local_y = ly;
            hit->local_z = lz;
            (*hit_count)++;
        }

        return hier_query_universe(sys, idx, fill_univ,
                                   lx - ox, ly - oy, lz - oz,
                                   depth + 1,
                                   out_hits, max_hits, hit_count);
    }

    idx->stats.point_exact_tests++;
    if (!alea_contains_point(sys, cell->root_node_id, lx, ly, lz)) {
        return 0;
    }

    if (*hit_count < max_hits) {
        alea_cell_hit_t* hit = &out_hits[*hit_count];
        hit->cell_id = cell->mc_cell_id;
        hit->cell_index = (int)cell_index;
        hit->material_id = cell->material_id;
        hit->universe_id = cell->universe_id;
        hit->fill_universe = cell->fill_universe;
        hit->depth = depth;
        hit->local_x = lx;
        hit->local_y = ly;
        hit->local_z = lz;
        (*hit_count)++;
    }

    if (cell->fill_universe > 0 && *hit_count < max_hits) {
        alea_matrix_t fill_transform;
        if (fill_transform_matrix(sys, cell, &fill_transform) != 0) {
            return -1;
        }
        if (!alea_matrix_invert(&fill_transform)) {
            return -1;
        }

        double child_x = lx;
        double child_y = ly;
        double child_z = lz;
        alea_matrix_transform_point_inverse(&fill_transform,
                                            &child_x, &child_y, &child_z);

        return hier_query_universe(sys, idx, cell->fill_universe,
                                   child_x, child_y, child_z,
                                   depth + 1,
                                   out_hits, max_hits, hit_count);
    }

    return 0;
}

static int hier_query_universe(alea_system_t* sys,
                               alea_hier_spatial_index_t* idx,
                               int universe_id,
                               double lx,
                               double ly,
                               double lz,
                               int depth,
                               alea_cell_hit_t* out_hits,
                               size_t max_hits,
                               size_t* hit_count) {
    if (*hit_count >= max_hits) return 0;
    if (depth >= HIER_MAX_PLACEMENT_DEPTH) return 0;

    const alea_universe_t* univ = alea_get_universe(sys, universe_id);
    if (!univ) return -1;

    const hier_universe_blas_t* blas = find_blas(idx, universe_id);
    if (blas && blas->built && blas->node_count > 0) {
        alea_size_vec_t candidates = ALEA_VEC_INIT;
        hier_point_query_ctx_t ctx = {
            .sys = sys,
            .idx = idx,
            .blas = blas,
            .x = lx,
            .y = ly,
            .z = lz,
            .candidates = &candidates,
            .error = 0
        };
        idx->stats.point_blas_queries++;
        query_blas_node(&ctx, 0);
        if (!ctx.error && candidates.count > 1) {
            qsort(candidates.data, candidates.count, sizeof(size_t),
                  compare_size_t_ascending);
        }
        for (size_t i = 0; !ctx.error && i < candidates.count; i++) {
            size_t cell_pos = candidates.data[i];
            if (cell_pos >= blas->cell_count) {
                ctx.error = -1;
                break;
            }
            uint32_t cell_index = blas->cells[cell_pos].cell_index;
            int rc = process_point_cell(sys, idx, cell_index, lx, ly, lz,
                                        depth, out_hits, max_hits, hit_count);
            if (rc < 0) ctx.error = -1;
        }
        alea_vec_free(&candidates);
        return ctx.error ? -1 : 0;
    }

    idx->stats.point_linear_scans++;
    idx->stats.point_linear_cell_tests += univ->cell_indices.count;
    for (size_t i = 0; i < univ->cell_indices.count; i++) {
        uint32_t cell_index = (uint32_t)univ->cell_indices.data[i];
        int rc = process_point_cell(sys, idx, cell_index, lx, ly, lz,
                                    depth, out_hits, max_hits, hit_count);
        if (rc < 0) return -1;
    }

    return 0;
}

int alea_hier_spatial_find_cells_at_point(alea_system_t* sys,
                                          double x,
                                          double y,
                                          double z,
                                          alea_cell_hit_t* out_hits,
                                          size_t max_hits) {
    if (!sys || !out_hits || max_hits == 0) return -1;

    if (!sys->hier_spatial_index || !sys->hier_spatial_index->built) {
        if (alea_hier_spatial_index_build(sys) != 0) {
            return -1;
        }
    }

    alea_hier_spatial_index_t* idx = sys->hier_spatial_index;
    idx->stats.point_queries++;

    size_t hit_count = 0;
    int rc = hier_query_universe(sys, idx, 0, x, y, z, 0,
                                 out_hits, max_hits, &hit_count);
    if (rc < 0) return -1;

    return (int)hit_count;
}

static int append_region_hit(alea_system_t* sys,
                             alea_hier_spatial_index_t* idx,
                             const hier_placement_t* placement,
                             uint32_t cell_index,
                             uint32_t synthetic_index,
                             alea_spatial_hit_t* out_hits,
                             size_t max_hits,
                             size_t* hit_count) {
    if (*hit_count >= max_hits) return 0;

    const alea_matrix_t* transform = placement_transform(idx, placement);
    if (!transform) return -1;

    const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
    alea_spatial_hit_t* hit = &out_hits[*hit_count];
    hit->instance_index = synthetic_index;
    hit->cell_index = cell_index;
    hit->cell_id = cell->mc_cell_id;
    hit->material_id = cell->material_id;
    hit->universe_id = placement->universe_id;
    hit->depth = placement->depth;
    hit->is_terminal = true;
    hit->transform = *transform;
    (*hit_count)++;
    return 0;
}

static int append_region_hit_direct(alea_system_t* sys,
                                    const alea_matrix_t* transform,
                                    int universe_id,
                                    int depth,
                                    uint32_t cell_index,
                                    uint32_t synthetic_index,
                                    alea_spatial_hit_t* out_hits,
                                    size_t max_hits,
                                    size_t* hit_count) {
    if (*hit_count >= max_hits) return 0;
    if (!transform) return -1;

    const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
    alea_spatial_hit_t* hit = &out_hits[*hit_count];
    hit->instance_index = synthetic_index;
    hit->cell_index = cell_index;
    hit->cell_id = cell->mc_cell_id;
    hit->material_id = cell->material_id;
    hit->universe_id = universe_id;
    hit->depth = depth;
    hit->is_terminal = true;
    hit->transform = *transform;
    (*hit_count)++;
    return 0;
}

static void translation_matrix(alea_matrix_t* out, double x, double y, double z) {
    alea_matrix_identity(out);
    out->m[3] = x;
    out->m[7] = y;
    out->m[11] = z;
    out->has_inverse = false;
    alea_matrix_invert(out);
}

static int query_region_blas_node(alea_system_t* sys,
                                  alea_hier_spatial_index_t* idx,
                                  const hier_placement_t* placement,
                                  const hier_universe_blas_t* blas,
                                  uint32_t node_index,
                                  const alea_bbox_t* local_query,
                                  const alea_bbox_t* world_query,
                                  alea_spatial_hit_t* out_hits,
                                  size_t max_hits,
                                  size_t* hit_count) {
    if (node_index >= blas->node_count) return -1;

    const alea_matrix_t* transform = placement_transform(idx, placement);
    if (!transform) return -1;

    const hier_bvh_node_t* node = &blas->nodes[node_index];
    if (!bbox_intersects_local(&node->bbox, local_query)) return 0;

    if (node->count == 0) {
        if (query_region_blas_node(sys, idx, placement, blas, node->left_or_first,
                                   local_query, world_query, out_hits, max_hits,
                                   hit_count) != 0) {
            return -1;
        }
        return query_region_blas_node(sys, idx, placement, blas, node->right_child,
                                      local_query, world_query, out_hits, max_hits,
                                      hit_count);
    }

    for (uint16_t i = 0; i < node->count; i++) {
        uint32_t cell_pos = blas->indices[node->left_or_first + i];
        if (cell_pos >= blas->cell_count) return -1;

        const hier_blas_cell_t* blas_cell = &blas->cells[cell_pos];
        if (!bbox_intersects_local(&blas_cell->bbox, local_query)) continue;

        const alea_cell_entry_t* cell = &sys->cells.data[blas_cell->cell_index];
        if (cell->fill_universe > 0 || (cell->lat_type != 0 && cell->lat_fill)) {
            continue;
        }

        alea_bbox_t world_bbox = bbox_transform(&blas_cell->bbox, transform);
        if (!bbox_intersects_local(&world_bbox, world_query)) continue;

        uint32_t synthetic_index = (uint32_t)(placement->id + cell_pos);
        if (append_region_hit(sys, idx, placement, blas_cell->cell_index,
                              synthetic_index, out_hits, max_hits, hit_count) != 0) {
            return -1;
        }
    }

    return 0;
}

static int query_region_blas(alea_system_t* sys,
                             alea_hier_spatial_index_t* idx,
                             const hier_placement_t* placement,
                             const hier_universe_blas_t* blas,
                             const alea_bbox_t* local_query,
                             const alea_bbox_t* world_query,
                             alea_spatial_hit_t* out_hits,
                             size_t max_hits,
                             size_t* hit_count) {
    if (!blas || blas->node_count == 0) return 0;
    return query_region_blas_node(sys, idx, placement, blas, 0,
                                  local_query, world_query,
                                  out_hits, max_hits, hit_count);
}

static int query_region_linear(alea_system_t* sys,
                               alea_hier_spatial_index_t* idx,
                               const hier_placement_t* placement,
                               const alea_universe_t* univ,
                               const alea_bbox_t* local_query,
                               const alea_bbox_t* world_query,
                               alea_spatial_hit_t* out_hits,
                               size_t max_hits,
                               size_t* hit_count) {
    const alea_matrix_t* transform = placement_transform(idx, placement);
    if (!transform) return -1;

    for (size_t i = 0; i < univ->cell_indices.count; i++) {
        uint32_t cell_index = (uint32_t)univ->cell_indices.data[i];
        const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
        if (cell->fill_universe > 0 || (cell->lat_type != 0 && cell->lat_fill)) {
            continue;
        }

        alea_bbox_t local_bbox = local_cell_bbox(sys, cell_index);
        if (!bbox_intersects_local(&local_bbox, local_query)) continue;

        alea_bbox_t world_bbox = bbox_transform(&local_bbox, transform);
        if (!bbox_intersects_local(&world_bbox, world_query)) continue;

        uint32_t synthetic_index = (uint32_t)(placement->id + i);
        if (append_region_hit(sys, idx, placement, cell_index,
                              synthetic_index, out_hits, max_hits, hit_count) != 0) {
            return -1;
        }
    }

    return 0;
}

static int query_region_universe_direct(alea_system_t* sys,
                                        alea_hier_spatial_index_t* idx,
                                        int universe_id,
                                        const alea_matrix_t* transform,
                                        int depth,
                                        const alea_bbox_t* local_query,
                                        const alea_bbox_t* world_query,
                                        alea_spatial_hit_t* out_hits,
                                        size_t max_hits,
                                        size_t* hit_count);

static int query_lattice_cell_direct(alea_system_t* sys,
                                     alea_hier_spatial_index_t* idx,
                                     const alea_cell_entry_t* cell,
                                     const alea_matrix_t* parent_transform,
                                     int depth,
                                     const alea_bbox_t* local_query,
                                     const alea_bbox_t* world_query,
                                     alea_spatial_hit_t* out_hits,
                                     size_t max_hits,
                                     size_t* hit_count) {
    if (!cell->lat_fill || cell->lat_fill_count == 0) return 0;

    int ni = cell->lat_fill_dims[1] - cell->lat_fill_dims[0] + 1;
    int nj = cell->lat_fill_dims[3] - cell->lat_fill_dims[2] + 1;
    int nk = cell->lat_fill_dims[5] - cell->lat_fill_dims[4] + 1;
    if (ni <= 0 || nj <= 0 || nk <= 0) return 0;

    int i0 = 0, i1 = ni - 1;
    int j0 = 0, j1 = nj - 1;
    int k0 = 0, k1 = nk - 1;

    if (cell->lat_type == 1) {
        if (cell->lat_pitch[0] == 0.0 || cell->lat_pitch[1] == 0.0 ||
            cell->lat_pitch[2] == 0.0) {
            return 0;
        }
        i0 = (int)floor((local_query->min_x - cell->lat_lower_left[0]) / cell->lat_pitch[0]);
        i1 = (int)floor((local_query->max_x - cell->lat_lower_left[0]) / cell->lat_pitch[0]);
        j0 = (nj == 1) ? 0 : (int)floor((local_query->min_y - cell->lat_lower_left[1]) / cell->lat_pitch[1]);
        j1 = (nj == 1) ? 0 : (int)floor((local_query->max_y - cell->lat_lower_left[1]) / cell->lat_pitch[1]);
        k0 = (nk == 1) ? 0 : (int)floor((local_query->min_z - cell->lat_lower_left[2]) / cell->lat_pitch[2]);
        k1 = (nk == 1) ? 0 : (int)floor((local_query->max_z - cell->lat_lower_left[2]) / cell->lat_pitch[2]);

        if (i0 < 0) i0 = 0;
        if (j0 < 0) j0 = 0;
        if (k0 < 0) k0 = 0;
        if (i1 >= ni) i1 = ni - 1;
        if (j1 >= nj) j1 = nj - 1;
        if (k1 >= nk) k1 = nk - 1;
    }

    if (i0 > i1 || j0 > j1 || k0 > k1) return 0;

    for (int i = i0; i <= i1; i++) {
        for (int j = j0; j <= j1; j++) {
            for (int k = k0; k <= k1; k++) {
                size_t fill_index = (size_t)(i * nj * nk + j * nk + k);
                if (fill_index >= cell->lat_fill_count) continue;

                int fill_universe = cell->lat_fill[fill_index];
                if (fill_universe < 0) continue;

                double ox, oy, oz;
                if (cell->lat_type == 2) {
                    int ri = i + cell->lat_fill_dims[0];
                    int rk = j + cell->lat_fill_dims[2];
                    double p = cell->lat_pitch[0] > 0.0 ? cell->lat_pitch[0] : 1.0;
                    ox = ri * p + rk * p * 0.5;
                    oy = rk * p * M_SQRT3 * 0.5;
                    oz = (nk == 1) ? 0.0 : cell->lat_lower_left[2] + (k + 0.5) * cell->lat_pitch[2];
                } else {
                    ox = cell->lat_lower_left[0] + (i + 0.5) * cell->lat_pitch[0];
                    oy = cell->lat_lower_left[1] + (j + 0.5) * cell->lat_pitch[1];
                    oz = cell->lat_lower_left[2] + (k + 0.5) * cell->lat_pitch[2];
                }

                alea_matrix_t element_translation;
                translation_matrix(&element_translation, ox, oy, oz);
                alea_matrix_t element_transform;
                alea_matrix_multiply(&element_transform, parent_transform, &element_translation);
                if (!alea_matrix_invert(&element_transform)) return -1;

                alea_bbox_t element_query = *local_query;
                element_query.min_x -= ox;
                element_query.max_x -= ox;
                element_query.min_y -= oy;
                element_query.max_y -= oy;
                element_query.min_z -= oz;
                element_query.max_z -= oz;

                if (query_region_universe_direct(sys, idx, fill_universe,
                                                 &element_transform, depth + 1,
                                                 &element_query, world_query,
                                                 out_hits, max_hits,
                                                 hit_count) != 0) {
                    return -1;
                }
            }
        }
    }

    return 0;
}

static int query_region_cell_direct(alea_system_t* sys,
                                    alea_hier_spatial_index_t* idx,
                                    uint32_t cell_index,
                                    size_t cell_pos,
                                    int universe_id,
                                    const alea_matrix_t* transform,
                                    int depth,
                                    const alea_bbox_t* cell_bbox,
                                    const alea_bbox_t* local_query,
                                    const alea_bbox_t* world_query,
                                    alea_spatial_hit_t* out_hits,
                                    size_t max_hits,
                                    size_t* hit_count) {
    if (!bbox_intersects_local(cell_bbox, local_query)) return 0;

    const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
    if (cell->lat_type != 0 && cell->lat_fill) {
        return query_lattice_cell_direct(sys, idx, cell, transform, depth,
                                         local_query, world_query,
                                         out_hits, max_hits, hit_count);
    }

    if (cell->fill_universe > 0) {
        alea_matrix_t fill_matrix;
        if (fill_transform_matrix(sys, cell, &fill_matrix) != 0) return -1;
        alea_matrix_t child_transform;
        alea_matrix_multiply(&child_transform, transform, &fill_matrix);
        if (!alea_matrix_invert(&child_transform)) return -1;

        alea_bbox_t child_query = bbox_transform_inverse_conservative(world_query,
                                                                      &child_transform);
        if (!alea_bbox_is_valid(&child_query)) return 0;

        return query_region_universe_direct(sys, idx, cell->fill_universe,
                                            &child_transform, depth + 1,
                                            &child_query, world_query,
                                            out_hits, max_hits, hit_count);
    }

    alea_bbox_t world_bbox = bbox_transform(cell_bbox, transform);
    if (!bbox_intersects_local(&world_bbox, world_query)) return 0;

    uint32_t synthetic_index = (uint32_t)(cell_pos + ((size_t)depth << 24));
    return append_region_hit_direct(sys, transform, universe_id, depth,
                                    cell_index, synthetic_index,
                                    out_hits, max_hits, hit_count);
}

static int query_region_blas_node_direct(alea_system_t* sys,
                                         alea_hier_spatial_index_t* idx,
                                         const hier_universe_blas_t* blas,
                                         uint32_t node_index,
                                         int universe_id,
                                         const alea_matrix_t* transform,
                                         int depth,
                                         const alea_bbox_t* local_query,
                                         const alea_bbox_t* world_query,
                                         alea_spatial_hit_t* out_hits,
                                         size_t max_hits,
                                         size_t* hit_count) {
    if (node_index >= blas->node_count) return -1;

    const hier_bvh_node_t* node = &blas->nodes[node_index];
    if (!bbox_intersects_local(&node->bbox, local_query)) return 0;

    if (node->count == 0) {
        if (query_region_blas_node_direct(sys, idx, blas, node->left_or_first,
                                          universe_id, transform, depth,
                                          local_query, world_query,
                                          out_hits, max_hits, hit_count) != 0) {
            return -1;
        }
        return query_region_blas_node_direct(sys, idx, blas, node->right_child,
                                             universe_id, transform, depth,
                                             local_query, world_query,
                                             out_hits, max_hits, hit_count);
    }

    for (uint16_t i = 0; i < node->count; i++) {
        uint32_t cell_pos = blas->indices[node->left_or_first + i];
        if (cell_pos >= blas->cell_count) return -1;

        const hier_blas_cell_t* blas_cell = &blas->cells[cell_pos];
        if (query_region_cell_direct(sys, idx, blas_cell->cell_index, cell_pos,
                                     universe_id, transform, depth,
                                     &blas_cell->bbox, local_query, world_query,
                                     out_hits, max_hits, hit_count) != 0) {
            return -1;
        }
    }

    return 0;
}

static int query_region_universe_direct(alea_system_t* sys,
                                        alea_hier_spatial_index_t* idx,
                                        int universe_id,
                                        const alea_matrix_t* transform,
                                        int depth,
                                        const alea_bbox_t* local_query,
                                        const alea_bbox_t* world_query,
                                        alea_spatial_hit_t* out_hits,
                                        size_t max_hits,
                                        size_t* hit_count) {
    const hier_universe_blas_t* blas = find_blas(idx, universe_id);
    if (blas && blas->built && blas->node_count > 0) {
        return query_region_blas_node_direct(sys, idx, blas, 0, universe_id,
                                             transform, depth, local_query,
                                             world_query, out_hits, max_hits,
                                             hit_count);
    }

    const alea_universe_t* univ = alea_get_universe(sys, universe_id);
    if (!univ) return 0;
    for (size_t i = 0; i < univ->cell_indices.count; i++) {
        uint32_t cell_index = (uint32_t)univ->cell_indices.data[i];
        const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
        alea_bbox_t bbox = (cell->lat_type != 0 && cell->lat_fill)
            ? lattice_container_bbox(cell)
            : local_cell_bbox(sys, cell_index);
        if (query_region_cell_direct(sys, idx, cell_index, i, universe_id,
                                     transform, depth, &bbox, local_query,
                                     world_query, out_hits, max_hits,
                                     hit_count) != 0) {
            return -1;
        }
    }
    return 0;
}

static int compare_spatial_hits_by_depth_cell(const void* a, const void* b) {
    const alea_spatial_hit_t* ha = (const alea_spatial_hit_t*)a;
    const alea_spatial_hit_t* hb = (const alea_spatial_hit_t*)b;
    if (ha->depth != hb->depth) {
        return (ha->depth > hb->depth) - (ha->depth < hb->depth);
    }
    if (ha->universe_id != hb->universe_id) {
        return (ha->universe_id > hb->universe_id) - (ha->universe_id < hb->universe_id);
    }
    return (ha->cell_index > hb->cell_index) - (ha->cell_index < hb->cell_index);
}

int alea_hier_spatial_query_region(alea_system_t* sys,
                                   const alea_bbox_t* query_bbox,
                                   alea_spatial_hit_t* out_hits,
                                   size_t max_hits) {
    if (!sys || !query_bbox || !out_hits) return -1;

    if (!sys->hier_spatial_index || !sys->hier_spatial_index->built) {
        if (alea_hier_spatial_index_build(sys) != 0) {
            return -1;
        }
    }

    alea_hier_spatial_index_t* idx = sys->hier_spatial_index;
    size_t hit_count = 0;

    for (size_t i = 0; i < idx->placement_count; i++) {
        const hier_placement_t* placement = &idx->placements[i];
        if (alea_bbox_is_valid(&placement->world_bbox) &&
            !bbox_intersects_local(&placement->world_bbox, query_bbox)) {
            continue;
        }

        const alea_matrix_t* transform = placement_transform(idx, placement);
        if (!transform) return -1;
        alea_bbox_t local_query = bbox_transform_inverse_conservative(query_bbox, transform);
        if (!alea_bbox_is_valid(&local_query)) continue;

        if (placement->flags & HIER_PLACEMENT_LATTICE) {
            if (placement->parent_cell_index >= sys->cells.count) return -1;
            const alea_cell_entry_t* cell = &sys->cells.data[placement->parent_cell_index];
            if (query_lattice_cell_direct(sys, idx, cell, transform,
                                          placement->depth, &local_query,
                                          query_bbox, out_hits, max_hits,
                                          &hit_count) != 0) {
                return -1;
            }
            continue;
        }

        const hier_universe_blas_t* blas = find_blas(idx, placement->universe_id);
        if (blas && blas->built) {
            if (query_region_blas(sys, idx, placement, blas, &local_query,
                                  query_bbox, out_hits, max_hits, &hit_count) != 0) {
                return -1;
            }
        } else {
            const alea_universe_t* univ = alea_get_universe(sys, placement->universe_id);
            if (!univ) continue;
            if (query_region_linear(sys, idx, placement, univ, &local_query,
                                    query_bbox, out_hits, max_hits, &hit_count) != 0) {
                return -1;
            }
        }
    }

    if (hit_count > 1) {
        qsort(out_hits, hit_count, sizeof(*out_hits), compare_spatial_hits_by_depth_cell);
    }

    return (int)hit_count;
}

int alea_hier_spatial_query_slice_z(alea_system_t* sys,
                                    double z,
                                    double x_min,
                                    double x_max,
                                    double y_min,
                                    double y_max,
                                    alea_spatial_hit_t* out_hits,
                                    size_t max_hits) {
    double ez = fabs(z) * 1e-10 + 1e-10;
    alea_bbox_t query = {
        .min_x = x_min, .max_x = x_max,
        .min_y = y_min, .max_y = y_max,
        .min_z = z - ez, .max_z = z + ez
    };

    return alea_hier_spatial_query_region(sys, &query, out_hits, max_hits);
}
