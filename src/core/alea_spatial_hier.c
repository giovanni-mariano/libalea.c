// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "core/alea_spatial_hier.h"
#include "core/alea_eval.h"
#include "core/alea_system.h"
#include "core/alea_universe.h"
#include "core/alea_surface.h"
#include "primitives/bbox.h"
#include "alea_types.h"
#include "util/alea_hashmap.h"
#include "util/alea_log.h"
#include "util/compat.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define HIER_BVH_LEAF_SIZE 4
#define HIER_DEFAULT_BLAS_THRESHOLD 1
#define HIER_MAX_PLACEMENT_DEPTH 64
#define HIER_DEFAULT_MEMORY_RESERVE_MB 256
#define HIER_MEMORY_CHECK_INTERVAL 16384

/* ----------------------------------------------------------------------
 * Hier coherence cache: mirror of g_cell_cache in alea_spatial.c.
 * Successive point queries along a ray are typically fractions of a mm
 * apart; without this cache, every step pays the full TLAS+BLAS descent.
 * -------------------------------------------------------------------- */
#define HIER_CACHE_MAX_DEPTH HIER_MAX_PLACEMENT_DEPTH

typedef struct {
    uint32_t cell_index;
    int cell_id;
    int material_id;
    int universe_id;
    int fill_universe;
    int depth;
    /* Forward world->this-depth's-local. Used via transform_point_inverse. */
    alea_matrix_t transform;
    /* Lattice "wrapper" hit (cell->lat_type != 0 && cell->lat_fill). On
     * verify we re-run the lattice lookup and require fill_universe and
     * element origin to match — otherwise the deeper cached transforms
     * (which embed the old element origin) are stale. */
    bool is_lattice;
    int lat_fill_universe;
    int lat_i, lat_j, lat_k;
    double lat_ox, lat_oy, lat_oz;
    bool valid;
} hier_cached_cell_t;

typedef struct {
    hier_cached_cell_t cells[HIER_CACHE_MAX_DEPTH];
    int count;
} hier_path_cache_t;

static ALEA_THREAD_LOCAL hier_cached_cell_t g_hier_cache[HIER_CACHE_MAX_DEPTH];
static ALEA_THREAD_LOCAL int g_hier_cache_count = 0;
static ALEA_THREAD_LOCAL alea_system_t* g_hier_cache_system = NULL;
static ALEA_THREAD_LOCAL uint64_t g_hier_cache_generation = 0;
static ALEA_THREAD_LOCAL uint64_t g_hier_cache_system_id = 0;

static inline void hier_cache_invalidate(void) {
    g_hier_cache_count = 0;
    for (int i = 0; i < HIER_CACHE_MAX_DEPTH; i++) {
        g_hier_cache[i].valid = false;
    }
}

static inline void hier_path_clear(hier_path_cache_t* path) {
    if (!path) return;
    path->count = 0;
    for (int i = 0; i < HIER_CACHE_MAX_DEPTH; i++) {
        path->cells[i].valid = false;
    }
}

void alea_hier_spatial_reset_cache(void) {
    hier_cache_invalidate();
    g_hier_cache_system = NULL;
    g_hier_cache_generation = 0;
    g_hier_cache_system_id = 0;
}

static size_t hier_memory_reserve_bytes(void) {
    const char* env = getenv("ALEA_HIER_MEMORY_RESERVE_MB");
    if (env && env[0]) {
        char* end = NULL;
        unsigned long mb = strtoul(env, &end, 10);
        if (end != env) return (size_t)mb * 1024ULL * 1024ULL;
    }
    return (size_t)HIER_DEFAULT_MEMORY_RESERVE_MB * 1024ULL * 1024ULL;
}

/* Returns 0 if memory is OK or unknown, -1 if low. When low, sets a clear
 * error so the caller can surface it instead of letting the kernel OOM-kill
 * the process. */
static int hier_memory_budget_check(const char* phase) {
    size_t avail = alea_mem_available_bytes();
    if (avail == 0) return 0; /* unknown platform / not enforceable */
    size_t reserve = hier_memory_reserve_bytes();
    if (avail >= reserve) return 0;
    alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
        "hier spatial index aborted in %s: available memory %.1f MiB below "
        "reserve %.1f MiB (set ALEA_HIER_MEMORY_RESERVE_MB to override)",
        phase ? phase : "build",
        (double)avail / (1024.0 * 1024.0),
        (double)reserve / (1024.0 * 1024.0));
    return -1;
}

static inline uint64_t hier_blas_key_hash(int id) {
    return (uint64_t)(unsigned int)id * 0x9E3779B97F4A7C15ULL;
}
static inline bool hier_blas_key_eq(int a, int b) { return a == b; }
ALEA_HASHMAP_DEFINE(hier_blas_map, int, uint32_t,
                    hier_blas_key_hash, hier_blas_key_eq, INT_MIN);
#ifndef M_SQRT3
#define M_SQRT3 1.73205080756887729353
#endif

enum {
    HIER_PLACEMENT_ROOT = 1u << 0,
    HIER_PLACEMENT_FILL = 1u << 1,
    HIER_PLACEMENT_LATTICE = 1u << 2
};

/* Compact bbox used for BVH traversal pruning. Stored as floats with each min
 * rounded down and each max rounded up so the bbox remains conservative (any
 * point/region intersecting the true double-precision bbox also intersects the
 * float version). Exact containment of cells is still verified against the
 * full CSG tree, so float storage here cannot affect correctness, only how
 * many candidates are visited. Halves node and per-cell footprint vs storing
 * 6 doubles. */
typedef struct {
    float min[3];
    float max[3];
} hier_fbbox_t;

static inline hier_fbbox_t hier_fbbox_empty(void) {
    hier_fbbox_t b;
    b.min[0] = b.min[1] = b.min[2] =  FLT_MAX;
    b.max[0] = b.max[1] = b.max[2] = -FLT_MAX;
    return b;
}

static inline hier_fbbox_t hier_fbbox_from_double(const alea_bbox_t* d) {
    if (!alea_bbox_is_valid(d)) return hier_fbbox_empty();
    hier_fbbox_t b;
    b.min[0] = nextafterf((float)d->min_x, -INFINITY);
    b.min[1] = nextafterf((float)d->min_y, -INFINITY);
    b.min[2] = nextafterf((float)d->min_z, -INFINITY);
    b.max[0] = nextafterf((float)d->max_x,  INFINITY);
    b.max[1] = nextafterf((float)d->max_y,  INFINITY);
    b.max[2] = nextafterf((float)d->max_z,  INFINITY);
    return b;
}

static inline hier_fbbox_t hier_fbbox_union(const hier_fbbox_t* a, const hier_fbbox_t* b) {
    hier_fbbox_t r;
    for (int i = 0; i < 3; i++) {
        r.min[i] = a->min[i] < b->min[i] ? a->min[i] : b->min[i];
        r.max[i] = a->max[i] > b->max[i] ? a->max[i] : b->max[i];
    }
    return r;
}

static inline int hier_fbbox_contains_point(const hier_fbbox_t* b,
                                            double x, double y, double z) {
    return x >= (double)b->min[0] && x <= (double)b->max[0] &&
           y >= (double)b->min[1] && y <= (double)b->max[1] &&
           z >= (double)b->min[2] && z <= (double)b->max[2];
}

static inline int hier_fbbox_intersects_dbbox(const hier_fbbox_t* a,
                                              const alea_bbox_t* b) {
    return (double)a->min[0] <= b->max_x && (double)a->max[0] >= b->min_x &&
           (double)a->min[1] <= b->max_y && (double)a->max[1] >= b->min_y &&
           (double)a->min[2] <= b->max_z && (double)a->max[2] >= b->min_z;
}

static inline int hier_fbbox_ray_enter_exit(const hier_fbbox_t* b,
                                            double ox,
                                            double oy,
                                            double oz,
                                            double inv_dx,
                                            double inv_dy,
                                            double inv_dz,
                                            double t_min,
                                            double t_max,
                                            double* out_enter,
                                            double* out_exit) {
    double t1x = ((double)b->min[0] - ox) * inv_dx;
    double t2x = ((double)b->max[0] - ox) * inv_dx;
    double enter = fmin(t1x, t2x);
    double exit = fmax(t1x, t2x);

    double t1y = ((double)b->min[1] - oy) * inv_dy;
    double t2y = ((double)b->max[1] - oy) * inv_dy;
    enter = fmax(enter, fmin(t1y, t2y));
    exit = fmin(exit, fmax(t1y, t2y));

    double t1z = ((double)b->min[2] - oz) * inv_dz;
    double t2z = ((double)b->max[2] - oz) * inv_dz;
    enter = fmax(enter, fmin(t1z, t2z));
    exit = fmin(exit, fmax(t1z, t2z));

    if (enter > t_max || exit < t_min) return 0;
    if (enter < t_min) enter = t_min;
    if (exit > t_max) exit = t_max;
    if (out_enter) *out_enter = enter;
    if (out_exit) *out_exit = exit;
    return 1;
}

static inline alea_bbox_t hier_fbbox_to_double(const hier_fbbox_t* b) {
    return (alea_bbox_t){
        (double)b->min[0], (double)b->max[0],
        (double)b->min[1], (double)b->max[1],
        (double)b->min[2], (double)b->max[2]
    };
}

typedef struct {
    hier_fbbox_t bbox;
    uint32_t left_or_first;
    uint32_t right_child;
    uint16_t count;
    uint8_t axis;
    uint8_t pad;
} hier_bvh_node_t;

typedef struct {
    uint32_t cell_index;
    float bbox_volume;  /* Precomputed for fast candidate ordering. */
    hier_fbbox_t bbox;
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
    uint32_t parent_id;
    uint32_t parent_cell_index;
    int universe_id;
    int depth;
    uint32_t transform_index;
    uint32_t flags;
    alea_bbox_t world_bbox;
} hier_placement_t;

struct alea_hier_spatial_index {
    hier_universe_blas_t* blas;
    size_t blas_count;
    size_t blas_capacity;
    hier_blas_map_t blas_by_universe;
    /* Top-level acceleration structure over placement world bboxes. Reuses the
     * BLAS struct: each "cell" in the TLAS is a placement and cell_index is
     * the placement index. Built once after collect_placements. */
    hier_universe_blas_t tlas;
    int tlas_built;
    hier_placement_t* placements;
    size_t placement_count;
    size_t placement_capacity;
    alea_matrix_t* transforms;
    size_t transform_count;
    size_t transform_capacity;
    /* Per-cell fill matrix cache. cell_fill_matrix_index[cell_index] holds
     * the slot in cell_fill_matrices for cells with fill_universe > 0, or
     * UINT32_MAX if absent. Each cached matrix has both forward and inverse
     * precomputed at build time so query paths avoid repeated alea_matrix_invert. */
    uint32_t* cell_fill_matrix_index;
    size_t cell_fill_matrix_index_size;
    alea_matrix_t* cell_fill_matrices;
    size_t cell_fill_matrix_count;
    alea_hier_spatial_stats_t stats;
    int built;
};

typedef struct {
    uint32_t index;
    double centroid[3];
} hier_bvh_item_t;

#define HIER_CAND_STACK 32

typedef struct {
    uint32_t stack[HIER_CAND_STACK];
    uint32_t* data;
    size_t count;
    size_t cap;
} hier_cand_buf_t;

static inline void hier_cand_init(hier_cand_buf_t* b) {
    b->data = b->stack;
    b->count = 0;
    b->cap = HIER_CAND_STACK;
}

static inline void hier_cand_free(hier_cand_buf_t* b) {
    if (b->data != b->stack) free(b->data);
}

/* Debug counter: candidates pushed across all universes in the current query.
 * Reset/read via alea_hier_debug_candidates_*(). Thread-local to avoid races. */
static ALEA_THREAD_LOCAL size_t g_hier_debug_candidates = 0;
void alea_hier_debug_candidates_reset(void) { g_hier_debug_candidates = 0; }
size_t alea_hier_debug_candidates_get(void) { return g_hier_debug_candidates; }

static int hier_cand_push(hier_cand_buf_t* b, uint32_t v) {
    g_hier_debug_candidates++;
    if (b->count == b->cap) {
        size_t new_cap = b->cap * 2;
        uint32_t* nd;
        if (b->data == b->stack) {
            nd = malloc(new_cap * sizeof(*nd));
            if (!nd) return -1;
            memcpy(nd, b->stack, b->count * sizeof(*nd));
        } else {
            nd = realloc(b->data, new_cap * sizeof(*nd));
            if (!nd) return -1;
        }
        b->data = nd;
        b->cap = new_cap;
    }
    b->data[b->count++] = v;
    return 0;
}

static inline void hier_cand_sort(hier_cand_buf_t* b) {
    /* Insertion sort: candidate counts are almost always small (<= leaf size
     * times a few overlapping leaves). Avoids qsort's function-pointer cost. */
    for (size_t i = 1; i < b->count; i++) {
        uint32_t x = b->data[i];
        size_t j = i;
        while (j > 0 && b->data[j - 1] > x) {
            b->data[j] = b->data[j - 1];
            j--;
        }
        b->data[j] = x;
    }
}

/* qsort comparator: smallest bbox_volume first. The blas pointer is passed
 * via a thread-local so we don't need qsort_r (which is non-portable). */
static ALEA_THREAD_LOCAL const hier_universe_blas_t* g_sort_blas = NULL;
static int hier_cand_cmp_by_volume(const void* a, const void* b) {
    float va = g_sort_blas->cells[*(const uint32_t*)a].bbox_volume;
    float vb = g_sort_blas->cells[*(const uint32_t*)b].bbox_volume;
    return (va > vb) - (va < vb);
}

/* Sort by bbox volume, smallest first. Used by find_cell_in_universe so the
 * most-likely-containing (tightest) cell is tested first and the search can
 * short-circuit at the first match. */
static inline void hier_cand_sort_by_volume(hier_cand_buf_t* b,
                                            const hier_universe_blas_t* blas) {
    g_sort_blas = blas;
    qsort(b->data, b->count, sizeof(*b->data), hier_cand_cmp_by_volume);
    g_sort_blas = NULL;
}

typedef struct {
    alea_system_t* sys;
    alea_hier_spatial_index_t* idx;
    const hier_universe_blas_t* blas;
    double x;
    double y;
    double z;
    hier_cand_buf_t* candidates;
    int error;
} hier_point_query_ctx_t;

static double monotonic_seconds(void) {
    return alea_monotonic_seconds();
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
        return alea_node_bbox_get(&root->bbox);
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
    if (cell->lat_fill_repeating) {
        const double inf = 1e30;
        return (alea_bbox_t){-inf, inf, -inf, inf, -inf, inf};
    }
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

    /* Hex lattice: element centers are at (ri*p + rk*p/2, rk*p*sqrt(3)/2) in
     * local coords (mirrors the canonical lattice mapping). The container bbox is
     * the rectangle that bounds those centers, padded by one pitch on each
     * side to cover the element's own footprint conservatively. Previously
     * this was a giant (ni+nj+2)*p square in both axes, which made hex
     * lattice placements effectively unbounded and defeated TLAS pruning. */
    double p = cell->lat_pitch[0] > 0.0 ? cell->lat_pitch[0] : 1.0;
    int d0 = cell->lat_fill_dims[0];
    int d1 = cell->lat_fill_dims[1];
    int d2 = cell->lat_fill_dims[2];
    int d3 = cell->lat_fill_dims[3];
    int nk = cell->lat_fill_dims[5] - cell->lat_fill_dims[4] + 1;

    double cx[4] = {
        (double)d0 * p + (double)d2 * p * 0.5,
        (double)d0 * p + (double)d3 * p * 0.5,
        (double)d1 * p + (double)d2 * p * 0.5,
        (double)d1 * p + (double)d3 * p * 0.5
    };
    double xmin_c = cx[0], xmax_c = cx[0];
    for (int i = 1; i < 4; i++) {
        if (cx[i] < xmin_c) xmin_c = cx[i];
        if (cx[i] > xmax_c) xmax_c = cx[i];
    }
    double ymin_c = (double)d2 * p * M_SQRT3 * 0.5;
    double ymax_c = (double)d3 * p * M_SQRT3 * 0.5;

    double zmin = (nk == 1) ? -p : cell->lat_lower_left[2];
    double zmax = (nk == 1) ?  p : cell->lat_lower_left[2] + (double)nk * cell->lat_pitch[2];

    return (alea_bbox_t){xmin_c - p, xmax_c + p,
                         ymin_c - p, ymax_c + p,
                         zmin, zmax};
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

static hier_fbbox_t item_bbox(const hier_universe_blas_t* blas,
                              const hier_bvh_item_t* items,
                              size_t start,
                              size_t end) {
    hier_fbbox_t result = hier_fbbox_empty();
    for (size_t i = start; i < end; i++) {
        result = hier_fbbox_union(&result, &blas->cells[items[i].index].bbox);
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

static void translation_matrix(alea_matrix_t* out, double x, double y, double z);

static inline const alea_matrix_t*
cell_fill_matrix_cached(const alea_hier_spatial_index_t* idx, uint32_t cell_index) {
    if (!idx->cell_fill_matrix_index) return NULL;
    if ((size_t)cell_index >= idx->cell_fill_matrix_index_size) return NULL;
    uint32_t slot = idx->cell_fill_matrix_index[cell_index];
    if (slot == UINT32_MAX) return NULL;
    return &idx->cell_fill_matrices[slot];
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
    node->bbox = hier_fbbox_union(&blas->nodes[left].bbox, &blas->nodes[right].bbox);
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
        blas->cells[i].bbox = hier_fbbox_from_double(&bbox);
        {
            float vdx = (float)(bbox.max_x - bbox.min_x);
            float vdy = (float)(bbox.max_y - bbox.min_y);
            float vdz = (float)(bbox.max_z - bbox.min_z);
            blas->cells[i].bbox_volume = vdx * vdy * vdz;
        }
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
    hier_blas_map_destroy(&idx->blas_by_universe);
    free_universe_blas(&idx->tlas);
    free(idx->placements);
    free(idx->transforms);
    free(idx->cell_fill_matrix_index);
    free(idx->cell_fill_matrices);
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
    placement->parent_id = parent_id;
    placement->parent_cell_index = parent_cell_index;
    placement->universe_id = universe_id;
    placement->depth = depth;
    placement->transform_index = transform_index;
    placement->flags = flags;
    placement->world_bbox = world_bbox ? *world_bbox : alea_bbox_empty();
    (void)local_bbox; /* unused: kept in signature for build-time intent */

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

static int build_cell_fill_matrix_cache(alea_system_t* sys,
                                        alea_hier_spatial_index_t* idx) {
    size_t cell_count = sys->cells.count;
    if (cell_count == 0) return 0;

    idx->cell_fill_matrix_index = malloc(cell_count * sizeof(uint32_t));
    if (!idx->cell_fill_matrix_index) return -1;
    idx->cell_fill_matrix_index_size = cell_count;
    for (size_t i = 0; i < cell_count; i++) {
        idx->cell_fill_matrix_index[i] = UINT32_MAX;
    }

    /* Deduplicate by MCNP fill_transform id. All cells with fill_transform == 0
     * (or unset) share one identity slot. Cells with the same MCNP TR id share
     * one slot. This keeps the table size proportional to unique transforms,
     * not to fill-cell count — critical for large models. */
    hier_blas_map_t xform_map = hier_blas_map_create(32);
    if (!xform_map.entries) return -1;

    size_t cap = 16;
    idx->cell_fill_matrices = malloc(cap * sizeof(alea_matrix_t));
    if (!idx->cell_fill_matrices) {
        hier_blas_map_destroy(&xform_map);
        return -1;
    }
    idx->cell_fill_matrix_count = 0;

    int rc = 0;
    for (size_t i = 0; i < cell_count; i++) {
        const alea_cell_entry_t* c = &sys->cells.data[i];
        if (c->fill_universe <= 0) continue;

        int key = c->fill_transform > 0 ? c->fill_transform : 0;
        uint32_t* existing = hier_blas_map_get(&xform_map, key);
        if (existing) {
            idx->cell_fill_matrix_index[i] = *existing;
            continue;
        }

        if (idx->cell_fill_matrix_count == cap) {
            size_t new_cap = cap * 2;
            alea_matrix_t* grown = realloc(idx->cell_fill_matrices,
                                           new_cap * sizeof(alea_matrix_t));
            if (!grown) { rc = -1; break; }
            idx->cell_fill_matrices = grown;
            cap = new_cap;
        }

        uint32_t slot = (uint32_t)idx->cell_fill_matrix_count;
        alea_matrix_t* m = &idx->cell_fill_matrices[slot];
        if (fill_transform_matrix(sys, c, m) != 0) { rc = -1; break; }
        if (!alea_matrix_invert(m)) { rc = -1; break; }
        if (!hier_blas_map_put(&xform_map, key, slot)) { rc = -1; break; }
        idx->cell_fill_matrix_index[i] = slot;
        idx->cell_fill_matrix_count++;
    }

    /* Trim allocation down to actual size; tiny tables stay tiny. */
    if (rc == 0 && idx->cell_fill_matrix_count > 0 &&
        idx->cell_fill_matrix_count < cap) {
        alea_matrix_t* trimmed = realloc(
            idx->cell_fill_matrices,
            idx->cell_fill_matrix_count * sizeof(alea_matrix_t));
        if (trimmed) idx->cell_fill_matrices = trimmed;
    }

    hier_blas_map_destroy(&xform_map);
    return rc;
}

static int collect_placements_recursive(alea_system_t* sys,
                                        alea_hier_spatial_index_t* idx,
                                        int universe_id,
                                        uint32_t parent_placement_id,
                                        const alea_matrix_t* accumulated,
                                        int depth) {
    if (depth >= HIER_MAX_PLACEMENT_DEPTH) return 0;
    if (g_alea_interrupted) return -1;

    /* Periodic memory budget check during the long-growing phase. The
     * placement count grows monotonically, so use it as the polling clock. */
    if ((idx->placement_count & (HIER_MEMORY_CHECK_INTERVAL - 1)) == 0 &&
        idx->placement_count > 0) {
        if (hier_memory_budget_check("collect_placements") != 0) return -1;
    }

    const alea_universe_t* univ = alea_get_universe(sys, universe_id);
    if (!univ) return 0;

    for (size_t i = 0; i < univ->cell_indices.count; i++) {
        uint32_t cell_index = (uint32_t)univ->cell_indices.data[i];
        const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
        alea_bbox_t local_bbox = (cell->lat_type != 0 && cell->lat_fill)
            ? lattice_container_bbox(cell)
            : local_cell_bbox(sys, cell_index);
        alea_bbox_t world_bbox = accumulated ? alea_bbox_transform(&local_bbox, accumulated) : local_bbox;

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
            const alea_matrix_t* fill_mat =
                cell_fill_matrix_cached(idx, cell_index);
            if (!fill_mat) return -1;

            alea_matrix_t child_transform;
            if (accumulated) {
                alea_matrix_multiply(&child_transform, accumulated, fill_mat);
            } else {
                child_transform = *fill_mat;
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

/* Build a top-level BVH over placement world bboxes. Stored as a phantom
 * hier_universe_blas_t where `cell_index` of each entry is the placement id.
 * Region queries traverse this BVH instead of scanning all placements. */
static int build_tlas(alea_hier_spatial_index_t* idx) {
    hier_universe_blas_t* tlas = &idx->tlas;
    memset(tlas, 0, sizeof(*tlas));
    tlas->cell_count = idx->placement_count;
    tlas->bounds = alea_bbox_empty();

    if (idx->placement_count == 0) {
        tlas->built = 1;
        idx->tlas_built = 1;
        return 0;
    }

    tlas->cells = calloc(idx->placement_count, sizeof(*tlas->cells));
    tlas->indices = calloc(idx->placement_count, sizeof(*tlas->indices));
    hier_bvh_item_t* items = calloc(idx->placement_count, sizeof(*items));
    if (!tlas->cells || !tlas->indices || !items) { free(items); return -1; }

    /* Use a finite "infinite" bbox sentinel for placements whose world bbox
     * could not be computed, so they are always traversed (matching the prior
     * flat-scan behavior for those cases) without breaking BVH math. */
    const double inf = 1e30;
    for (size_t i = 0; i < idx->placement_count; i++) {
        alea_bbox_t b = idx->placements[i].world_bbox;
        if (!alea_bbox_is_valid(&b)) {
            b = (alea_bbox_t){-inf, inf, -inf, inf, -inf, inf};
        }
        tlas->cells[i].cell_index = (uint32_t)i;
        tlas->cells[i].bbox = hier_fbbox_from_double(&b);
        tlas->bounds = alea_bbox_union(&tlas->bounds, &b);
        items[i].index = (uint32_t)i;
        items[i].centroid[0] = 0.5 * (b.min_x + b.max_x);
        items[i].centroid[1] = 0.5 * (b.min_y + b.max_y);
        items[i].centroid[2] = 0.5 * (b.min_z + b.max_z);
    }

    size_t est = estimate_bvh_node_count(idx->placement_count);
    tlas->nodes = calloc(est ? est : 1, sizeof(*tlas->nodes));
    if (!tlas->nodes) { free(items); return -1; }
    tlas->node_capacity = est;

    uint32_t root = build_bvh_recursive(tlas, items, 0, idx->placement_count, 0);
    if (root == UINT32_MAX) { free(items); return -1; }

    for (size_t i = 0; i < idx->placement_count; i++) {
        tlas->indices[i] = items[i].index;
    }
    free(items);
    tlas->built = 1;
    idx->tlas_built = 1;
    return 0;
}

int alea_hier_spatial_index_build(alea_system_t* sys) {
    if (!sys) return -1;

    double t_start = monotonic_seconds();
    size_t threshold = blas_threshold();

    if (hier_memory_budget_check("build start") != 0) return -1;

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
    idx->blas_by_universe = hier_blas_map_create(
        sys->universes.count > 0 ? sys->universes.count * 2 : 16);
    if (!idx->blas_by_universe.entries) {
        free(idx);
        return -1;
    }

    if (build_cell_fill_matrix_cache(sys, idx) != 0) {
        alea_hier_spatial_index_free(idx);
        return -1;
    }

    if (collect_placements(sys, idx) != 0) {
        alea_hier_spatial_index_free(idx);
        return -1;
    }

    if (hier_memory_budget_check("tlas build") != 0) {
        alea_hier_spatial_index_free(idx);
        return -1;
    }

    if (build_tlas(idx) != 0) {
        alea_hier_spatial_index_free(idx);
        return -1;
    }

    if (hier_memory_budget_check("blas build") != 0) {
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
        if (!hier_blas_map_put(&idx->blas_by_universe, univ->universe_id,
                               (uint32_t)idx->blas_count)) {
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
    idx->stats.memory_bytes += idx->blas_by_universe.capacity *
                               sizeof(*idx->blas_by_universe.entries);
    idx->stats.memory_bytes += idx->placement_capacity * sizeof(*idx->placements);
    idx->stats.memory_bytes += idx->transform_capacity * sizeof(*idx->transforms);
    idx->stats.memory_bytes += idx->cell_fill_matrix_index_size *
                               sizeof(*idx->cell_fill_matrix_index);
    idx->stats.memory_bytes += idx->cell_fill_matrix_count *
                               sizeof(*idx->cell_fill_matrices);
    idx->stats.memory_bytes += idx->tlas.cell_count * sizeof(*idx->tlas.cells);
    idx->stats.memory_bytes += idx->tlas.cell_count * sizeof(*idx->tlas.indices);
    idx->stats.memory_bytes += idx->tlas.node_capacity * sizeof(*idx->tlas.nodes);
    idx->built = 1;
    sys->hier_spatial_index = idx;
    atomic_fetch_or(&sys->query_cache_state, ALEA_CACHE_HIER_SPATIAL);

    double t_end = monotonic_seconds();
    ALEA_LOG_INFO("Hier spatial cell-fill matrix cache: unique_transforms=%zu cells_indexed=%zu memory=%.2f MiB",
                  idx->cell_fill_matrix_count,
                  idx->cell_fill_matrix_index_size,
                  bytes_to_mib(
                      idx->cell_fill_matrix_count * sizeof(alea_matrix_t) +
                      idx->cell_fill_matrix_index_size * sizeof(uint32_t)));
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
    uint32_t* slot = hier_blas_map_get(&idx->blas_by_universe, universe_id);
    if (!slot) return NULL;
    if (*slot >= idx->blas_count) return NULL;
    return &idx->blas[*slot];
}

static const alea_matrix_t* placement_transform(const alea_hier_spatial_index_t* idx,
                                                const hier_placement_t* placement) {
    if (!idx || !placement) return NULL;
    if (placement->transform_index >= idx->transform_count) return NULL;
    return &idx->transforms[placement->transform_index];
}

int alea_hier_spatial_check_placement_chain(alea_system_t* sys,
                                            uint32_t placement_index,
                                            double x,
                                            double y,
                                            double z) {
    if (!sys || !sys->hier_spatial_index) return -1;
    alea_hier_spatial_index_t* idx = sys->hier_spatial_index;
    if (!idx->built || (size_t)placement_index >= idx->placement_count) {
        return -1;
    }

    uint32_t current = placement_index;
    int guard = 0;
    while ((size_t)current < idx->placement_count &&
           guard++ < HIER_MAX_PLACEMENT_DEPTH) {
        const hier_placement_t* placement = &idx->placements[current];
        if (placement->parent_id == UINT32_MAX ||
            placement->parent_cell_index == UINT32_MAX) {
            return 1;
        }
        if ((size_t)placement->parent_id >= idx->placement_count ||
            (size_t)placement->parent_cell_index >= alea_vec_count(&sys->cells)) {
            return -1;
        }

        const hier_placement_t* parent = &idx->placements[placement->parent_id];
        const alea_matrix_t* parent_transform = placement_transform(idx, parent);
        if (!parent_transform) return -1;

        double lx = x;
        double ly = y;
        double lz = z;
        alea_matrix_transform_point_inverse(parent_transform, &lx, &ly, &lz);

        const alea_cell_entry_t* parent_cell =
            &sys->cells.data[placement->parent_cell_index];
        if (parent_cell->root_node_id == ALEA_NODE_ID_INVALID ||
            !alea_contains_point(sys, parent_cell->root_node_id, lx, ly, lz)) {
            return 0;
        }

        current = placement->parent_id;
    }

    return -1;
}

int alea_hier_spatial_check_lattice_placement_ancestors(
    alea_system_t* sys,
    uint32_t placement_index,
    double x,
    double y,
    double z) {
    if (!sys || !sys->hier_spatial_index) return -1;
    alea_hier_spatial_index_t* idx = sys->hier_spatial_index;
    if (!idx->built || (size_t)placement_index >= idx->placement_count) {
        return -1;
    }

    /* A lattice placement represents the whole finite array.  Its own parent
     * cell is the repeated unit-cell CSG and therefore cannot validate an
     * array-level point.  Begin with the enclosing placement instead. */
    uint32_t current = idx->placements[placement_index].parent_id;
    int guard = 0;
    while ((size_t)current < idx->placement_count &&
           guard++ < HIER_MAX_PLACEMENT_DEPTH) {
        const hier_placement_t* placement = &idx->placements[current];
        if (placement->parent_id == UINT32_MAX ||
            placement->parent_cell_index == UINT32_MAX) {
            return 1;
        }
        if ((size_t)placement->parent_id >= idx->placement_count ||
            (size_t)placement->parent_cell_index >= alea_vec_count(&sys->cells)) {
            return -1;
        }

        const hier_placement_t* parent = &idx->placements[placement->parent_id];
        const alea_matrix_t* parent_transform = placement_transform(idx, parent);
        if (!parent_transform) return -1;

        double lx = x;
        double ly = y;
        double lz = z;
        alea_matrix_transform_point_inverse(parent_transform, &lx, &ly, &lz);

        const alea_cell_entry_t* parent_cell =
            &sys->cells.data[placement->parent_cell_index];
        if (parent_cell->root_node_id == ALEA_NODE_ID_INVALID ||
            !alea_contains_point(sys, parent_cell->root_node_id, lx, ly, lz)) {
            return 0;
        }

        current = placement->parent_id;
    }

    return -1;
}

static void query_blas_node(hier_point_query_ctx_t* ctx, uint32_t node_idx) {
    if (ctx->error) return;
    if (node_idx >= ctx->blas->node_count) {
        ctx->error = -1;
        return;
    }

    const hier_bvh_node_t* node = &ctx->blas->nodes[node_idx];
    /* Per-node/per-bbox stats counters are deliberately absent from this hot
     * path: shared writes to idx->stats cause cache-line ping-pong under
     * OpenMP. Query attribution is collected by worker-local callers. */
    if (!hier_fbbox_contains_point(&node->bbox, ctx->x, ctx->y, ctx->z)) {
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
            if (!hier_fbbox_contains_point(&blas_cell->bbox, ctx->x, ctx->y, ctx->z)) {
                continue;
            }
            if (hier_cand_push(ctx->candidates, cell_pos) != 0) {
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
                               const alea_matrix_t* parent_transform,
                               alea_cell_hit_t* out_hits,
                               size_t max_hits,
                               size_t* hit_count);

static inline void hier_cache_append(uint32_t cell_index,
                                     const alea_cell_hit_t* hit,
                                     const alea_matrix_t* transform,
                                     bool is_lattice,
                                     int lat_fill_universe,
                                     int lat_i,
                                     int lat_j,
                                     int lat_k,
                                     double lat_ox,
                                     double lat_oy,
                                     double lat_oz) {
    if (g_hier_cache_count >= HIER_CACHE_MAX_DEPTH) return;
    hier_cached_cell_t* ent = &g_hier_cache[g_hier_cache_count++];
    ent->cell_index = cell_index;
    ent->cell_id = hit->cell_id;
    ent->material_id = hit->material_id;
    ent->universe_id = hit->universe_id;
    ent->fill_universe = hit->fill_universe;
    ent->depth = hit->depth;
    ent->transform = *transform;
    ent->is_lattice = is_lattice;
    ent->lat_fill_universe = lat_fill_universe;
    ent->lat_i = lat_i;
    ent->lat_j = lat_j;
    ent->lat_k = lat_k;
    ent->lat_ox = lat_ox;
    ent->lat_oy = lat_oy;
    ent->lat_oz = lat_oz;
    ent->valid = true;
}

static inline void hier_path_append(hier_path_cache_t* path,
                                    uint32_t cell_index,
                                    const alea_cell_entry_t* cell,
                                    int fill_universe,
                                    int depth,
                                    const alea_matrix_t* transform,
                                    bool is_lattice,
                                    int lat_fill_universe,
                                    int lat_i,
                                    int lat_j,
                                    int lat_k,
                                    double lat_ox,
                                    double lat_oy,
                                    double lat_oz) {
    if (!path || path->count >= HIER_CACHE_MAX_DEPTH) return;
    hier_cached_cell_t* ent = &path->cells[path->count++];
    ent->cell_index = cell_index;
    ent->cell_id = cell->mc_cell_id;
    ent->material_id = cell->material_id;
    ent->universe_id = cell->universe_id;
    ent->fill_universe = fill_universe;
    ent->depth = depth;
    ent->transform = *transform;
    ent->is_lattice = is_lattice;
    ent->lat_fill_universe = lat_fill_universe;
    ent->lat_i = lat_i;
    ent->lat_j = lat_j;
    ent->lat_k = lat_k;
    ent->lat_ox = lat_ox;
    ent->lat_oy = lat_oy;
    ent->lat_oz = lat_oz;
    ent->valid = true;
}

static inline void hier_path_commit_to_path(hier_path_cache_t* dst,
                                            const hier_path_cache_t* src) {
    if (!dst || !src) return;
    dst->count = src->count;
    for (int i = 0; i < src->count && i < HIER_CACHE_MAX_DEPTH; i++) {
        dst->cells[i] = src->cells[i];
    }
    for (int i = src->count; i < HIER_CACHE_MAX_DEPTH; i++) {
        dst->cells[i].valid = false;
    }
}

static inline void hier_path_commit_to_cache(const hier_path_cache_t* path) {
    hier_cache_invalidate();
    if (!path) return;
    g_hier_cache_count = path->count;
    if (g_hier_cache_count > HIER_CACHE_MAX_DEPTH) {
        g_hier_cache_count = HIER_CACHE_MAX_DEPTH;
    }
    for (int i = 0; i < g_hier_cache_count; i++) {
        g_hier_cache[i] = path->cells[i];
    }
}

static inline void hier_path_export(const hier_path_cache_t* src,
                                    alea_hier_ray_path_t* dst) {
    if (!dst) return;
    dst->count = 0;
    if (!src) return;

    int count = src->count;
    if (count > ALEA_HIER_RAY_PATH_MAX) count = ALEA_HIER_RAY_PATH_MAX;
    dst->count = count;
    for (int i = 0; i < count; i++) {
        const hier_cached_cell_t* in = &src->cells[i];
        alea_hier_ray_path_entry_t* out = &dst->entries[i];
        out->cell_index = in->cell_index;
        out->cell_id = in->cell_id;
        out->material_id = in->material_id;
        out->universe_id = in->universe_id;
        out->fill_universe = in->fill_universe;
        out->depth = in->depth;
        out->is_lattice = in->is_lattice ? 1 : 0;
        out->lat_fill_universe = in->lat_fill_universe;
        out->lat_i = in->lat_i;
        out->lat_j = in->lat_j;
        out->lat_k = in->lat_k;
        out->lat_ox = in->lat_ox;
        out->lat_oy = in->lat_oy;
        out->lat_oz = in->lat_oz;
        out->transform = in->transform;
    }
}

static inline void hier_path_import_prefix(const alea_hier_ray_path_t* src,
                                           int prefix_count,
                                           hier_path_cache_t* dst) {
    hier_path_clear(dst);
    if (!src || !dst || prefix_count <= 0) return;

    int count = prefix_count;
    if (count > src->count) count = src->count;
    if (count > HIER_CACHE_MAX_DEPTH) count = HIER_CACHE_MAX_DEPTH;
    dst->count = count;
    for (int i = 0; i < count; i++) {
        const alea_hier_ray_path_entry_t* in = &src->entries[i];
        hier_cached_cell_t* out = &dst->cells[i];
        out->cell_index = in->cell_index;
        out->cell_id = in->cell_id;
        out->material_id = in->material_id;
        out->universe_id = in->universe_id;
        out->fill_universe = in->fill_universe;
        out->depth = in->depth;
        out->transform = in->transform;
        out->is_lattice = in->is_lattice != 0;
        out->lat_fill_universe = in->lat_fill_universe;
        out->lat_i = in->lat_i;
        out->lat_j = in->lat_j;
        out->lat_k = in->lat_k;
        out->lat_ox = in->lat_ox;
        out->lat_oy = in->lat_oy;
        out->lat_oz = in->lat_oz;
        out->valid = true;
    }
}

static int hier_cache_try(alea_system_t* sys, double x, double y, double z,
                          alea_cell_hit_t* out_hits, size_t max_hits) {
    if (g_hier_cache_count <= 0) return -1;

    /* A cached chain that ends in a fill/lattice container is an
     * undefined-fill fallback: the fill had no content at the cached point.
     * At a different point the same container may resolve deeper, and no
     * containment check on the cached entries can detect that — the missing
     * depth is exactly what is not in the cache. Force a full descent. */
    {
        const hier_cached_cell_t* last = &g_hier_cache[g_hier_cache_count - 1];
        if (!last->valid) return -1;
        if ((size_t)last->cell_index >= alea_vec_count(&sys->cells)) return -1;
        if (!last->is_lattice &&
            alea_cell_entry_is_container(&sys->cells.data[last->cell_index]))
            return -1;
        if (last->is_lattice) return -1;
    }

    size_t hit_count = 0;
    for (int i = 0; i < g_hier_cache_count && hit_count < max_hits; i++) {
        hier_cached_cell_t* ent = &g_hier_cache[i];
        if (!ent->valid) return -1;

        double lx = x, ly = y, lz = z;
        alea_matrix_transform_point_inverse(&ent->transform, &lx, &ly, &lz);

        const alea_cell_entry_t* cell = &sys->cells.data[ent->cell_index];

        if (ent->is_lattice) {
            alea_lattice_location_t location;
            if (alea_lattice_locate_point(sys, cell, lx, ly, lz, &location) != 1 ||
                location.fill_universe != ent->lat_fill_universe) return -1;
            if (location.i != ent->lat_i || location.j != ent->lat_j ||
                location.k != ent->lat_k) return -1;
            /* Lattice lookups derive (ox,oy,oz) from integer indices × pitch,
             * so bitwise equality across calls is the right test. */
            if (location.ox != ent->lat_ox || location.oy != ent->lat_oy ||
                location.oz != ent->lat_oz) {
                return -1;
            }
        } else {
            if (!alea_contains_point(sys, cell->root_node_id, lx, ly, lz)) {
                return -1;
            }
        }

        alea_cell_hit_t* hit = &out_hits[hit_count++];
        hit->cell_id = ent->cell_id;
        hit->cell_index = (int)ent->cell_index;
        hit->material_id = ent->material_id;
        hit->universe_id = ent->universe_id;
        hit->fill_universe = ent->fill_universe;
        hit->depth = ent->depth;
        hit->local_x = lx;
        hit->local_y = ly;
        hit->local_z = lz;
        hit->resolution_flags = alea_cell_entry_is_container(cell)
            ? ALEA_RESOLVE_UNDEFINED_FILL : 0;
    }
    return (int)hit_count;
}

static int process_point_cell(alea_system_t* sys,
                              alea_hier_spatial_index_t* idx,
                              uint32_t cell_index,
                              double lx,
                              double ly,
                              double lz,
                              int depth,
                              const alea_matrix_t* parent_transform,
                              alea_cell_hit_t* out_hits,
                              size_t max_hits,
                              size_t* hit_count) {
    const alea_cell_entry_t* cell = &sys->cells.data[cell_index];

    if (cell->lat_type != 0 && cell->lat_fill) {
        alea_lattice_location_t location;
        if (alea_lattice_locate_point(sys, cell, lx, ly, lz, &location) != 1)
            return 0;

        if (*hit_count < max_hits) {
            alea_cell_hit_t* hit = &out_hits[*hit_count];
            hit->cell_id = cell->mc_cell_id;
            hit->cell_index = (int)cell_index;
            hit->material_id = cell->material_id;
            hit->universe_id = cell->universe_id;
            hit->fill_universe = location.fill_universe;
            hit->depth = depth;
            hit->local_x = lx;
            hit->local_y = ly;
            hit->local_z = lz;
            hit->resolution_flags = ALEA_RESOLVE_UNDEFINED_FILL;
            hier_cache_append(cell_index, hit, parent_transform, true,
                              location.fill_universe,
                              location.i, location.j, location.k,
                              location.ox,
                              location.oy, location.oz);
            (*hit_count)++;
        }

        alea_matrix_t element_translation, child_transform;
        translation_matrix(&element_translation, location.ox, location.oy,
                           location.oz);
        alea_matrix_multiply(&child_transform, parent_transform,
                             &element_translation);
        if (!alea_matrix_invert(&child_transform)) return -1;

        return hier_query_universe(sys, idx, location.fill_universe,
                                   lx - location.ox, ly - location.oy,
                                   lz - location.oz,
                                   depth + 1,
                                   &child_transform,
                                   out_hits, max_hits, hit_count);
    }

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
        hit->resolution_flags = alea_cell_entry_is_container(cell)
            ? ALEA_RESOLVE_UNDEFINED_FILL : 0;
        hier_cache_append(cell_index, hit, parent_transform,
                          false, 0, 0, 0, 0, 0.0, 0.0, 0.0);
        (*hit_count)++;
    }

    if (cell->fill_universe > 0 && *hit_count < max_hits) {
        const alea_matrix_t* fill_transform =
            cell_fill_matrix_cached(idx, cell_index);
        if (!fill_transform) return -1;

        double child_x = lx;
        double child_y = ly;
        double child_z = lz;
        alea_matrix_transform_point_inverse(fill_transform,
                                            &child_x, &child_y, &child_z);

        alea_matrix_t child_transform;
        alea_matrix_multiply(&child_transform, parent_transform, fill_transform);
        if (!alea_matrix_invert(&child_transform)) return -1;

        return hier_query_universe(sys, idx, cell->fill_universe,
                                   child_x, child_y, child_z,
                                   depth + 1,
                                   &child_transform,
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
                               const alea_matrix_t* parent_transform,
                               alea_cell_hit_t* out_hits,
                               size_t max_hits,
                               size_t* hit_count) {
    if (*hit_count >= max_hits) return 0;
    if (depth >= HIER_MAX_PLACEMENT_DEPTH) return 0;

    const alea_universe_t* univ = alea_get_universe(sys, universe_id);
    if (!univ) return -1;

    const hier_universe_blas_t* blas = find_blas(idx, universe_id);
    if (blas && blas->built && blas->node_count > 0) {
        hier_cand_buf_t candidates;
        hier_cand_init(&candidates);
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
        query_blas_node(&ctx, 0);
        if (!ctx.error && candidates.count > 1) {
            hier_cand_sort(&candidates);
        }
        for (size_t i = 0; !ctx.error && i < candidates.count; i++) {
            uint32_t cell_pos = candidates.data[i];
            if (cell_pos >= blas->cell_count) {
                ctx.error = -1;
                break;
            }
            uint32_t cell_index = blas->cells[cell_pos].cell_index;
            int rc = process_point_cell(sys, idx, cell_index, lx, ly, lz,
                                        depth, parent_transform,
                                        out_hits, max_hits, hit_count);
            if (rc < 0) ctx.error = -1;
        }
        hier_cand_free(&candidates);
        return ctx.error ? -1 : 0;
    }

    for (size_t i = 0; i < univ->cell_indices.count; i++) {
        uint32_t cell_index = (uint32_t)univ->cell_indices.data[i];
        int rc = process_point_cell(sys, idx, cell_index, lx, ly, lz,
                                    depth, parent_transform,
                                    out_hits, max_hits, hit_count);
        if (rc < 0) return -1;
    }

    return 0;
}

static void hier_deepest_store_hit(const alea_cell_entry_t* cell,
                                   uint32_t cell_index,
                                   int fill_universe,
                                   double lx,
                                   double ly,
                                   double lz,
                                   int depth,
                                   const alea_matrix_t* transform,
                                   int lattice_cell_index,
                                   const alea_matrix_t* lattice_transform,
                                   alea_hier_cell_hit_t* out_hit,
                                   int* found) {
    out_hit->hit.cell_id = cell->mc_cell_id;
    out_hit->hit.cell_index = (int)cell_index;
    out_hit->hit.material_id = cell->material_id;
    out_hit->hit.universe_id = cell->universe_id;
    out_hit->hit.fill_universe = fill_universe;
    out_hit->hit.depth = depth;
    out_hit->hit.local_x = lx;
    out_hit->hit.local_y = ly;
    out_hit->hit.local_z = lz;
    /* A container survives as the final hit only when its fill descent
     * stores nothing deeper — i.e. an undefined-fill fallback. Deeper
     * stores overwrite this flag along with the rest of the hit. */
    out_hit->hit.resolution_flags = alea_cell_entry_is_container(cell)
        ? ALEA_RESOLVE_UNDEFINED_FILL : 0;
    out_hit->transform = *transform;
    out_hit->lattice_cell_index = lattice_cell_index;
    if (lattice_cell_index >= 0 && lattice_transform) {
        out_hit->lattice_transform = *lattice_transform;
    } else {
        alea_matrix_identity(&out_hit->lattice_transform);
    }
    *found = 1;
}

static int hier_find_deepest_universe(alea_system_t* sys,
                                      alea_hier_spatial_index_t* idx,
                                      int universe_id,
                                      double lx,
                                      double ly,
                                      double lz,
                                      int depth,
                                      const alea_matrix_t* transform,
                                      int lattice_cell_index,
                                      const alea_matrix_t* lattice_transform,
                                      alea_hier_cell_hit_t* out_hit,
                                      int* found,
                                      hier_path_cache_t* path,
                                      hier_path_cache_t* best_path);

/* Returns 1 if the candidate cell contains the point (lattice cells: the
 * lattice lookup resolved an element), 0 if not, -1 on error. Containing
 * candidates store their hit (and descend into fills) before returning. */
static int hier_find_deepest_cell(alea_system_t* sys,
                                  alea_hier_spatial_index_t* idx,
                                  uint32_t cell_index,
                                  double lx,
                                  double ly,
                                  double lz,
                                  int depth,
                                  const alea_matrix_t* transform,
                                  int lattice_cell_index,
                                  const alea_matrix_t* lattice_transform,
                                  alea_hier_cell_hit_t* out_hit,
                                  int* found,
                                  hier_path_cache_t* path,
                                  hier_path_cache_t* best_path) {
    const alea_cell_entry_t* cell = &sys->cells.data[cell_index];

    if (cell->lat_type != 0 && cell->lat_fill) {
        int saved_path_count = path ? path->count : 0;
        alea_lattice_location_t location;
        if (alea_lattice_locate_point(sys, cell, lx, ly, lz, &location) != 1)
            return 0;

        hier_path_append(path, cell_index, cell, location.fill_universe,
                         depth, transform,
                         true, location.fill_universe,
                         location.i, location.j, location.k, location.ox,
                         location.oy, location.oz);

        hier_deepest_store_hit(cell, cell_index, location.fill_universe,
                               lx, ly, lz, depth, transform,
                               lattice_cell_index, lattice_transform,
                               out_hit, found);
        hier_path_commit_to_path(best_path, path);

        alea_matrix_t element_translation;
        translation_matrix(&element_translation, location.ox, location.oy,
                           location.oz);

        alea_matrix_t element_transform;
        alea_matrix_multiply(&element_transform, transform, &element_translation);
        if (!alea_matrix_invert(&element_transform)) {
            if (path) path->count = saved_path_count;
            return -1;
        }

        int rc = hier_find_deepest_universe(sys, idx, location.fill_universe,
                                            lx - location.ox, ly - location.oy,
                                            lz - location.oz,
                                            depth + 1,
                                            &element_transform,
                                            (int)cell_index, transform,
                                            out_hit, found,
                                            path, best_path);
        if (path) path->count = saved_path_count;
        return rc < 0 ? -1 : 1;
    }

    if (!alea_contains_point(sys, cell->root_node_id, lx, ly, lz)) {
        return 0;
    }

    int saved_path_count = path ? path->count : 0;
    hier_path_append(path, cell_index, cell, cell->fill_universe,
                     depth, transform,
                     false, 0, 0, 0, 0, 0.0, 0.0, 0.0);

    hier_deepest_store_hit(cell, cell_index, cell->fill_universe,
                           lx, ly, lz, depth, transform,
                           lattice_cell_index, lattice_transform,
                           out_hit, found);
    hier_path_commit_to_path(best_path, path);

    if (cell->fill_universe > 0) {
        const alea_matrix_t* fill_transform =
            cell_fill_matrix_cached(idx, cell_index);
        if (!fill_transform) {
            if (path) path->count = saved_path_count;
            return -1;
        }

        alea_matrix_t child_transform;
        alea_matrix_multiply(&child_transform, transform, fill_transform);
        if (!alea_matrix_invert(&child_transform)) {
            if (path) path->count = saved_path_count;
            return -1;
        }

        double child_x = lx;
        double child_y = ly;
        double child_z = lz;
        alea_matrix_transform_point_inverse(fill_transform,
                                            &child_x, &child_y, &child_z);

        int rc = hier_find_deepest_universe(sys, idx, cell->fill_universe,
                                            child_x, child_y, child_z,
                                            depth + 1,
                                            &child_transform,
                                            lattice_cell_index,
                                            lattice_transform,
                                            out_hit, found,
                                            path, best_path);
        if (path) path->count = saved_path_count;
        return rc < 0 ? -1 : 1;
    }

    if (path) path->count = saved_path_count;
    return 1;
}

static int hier_find_deepest_universe(alea_system_t* sys,
                                      alea_hier_spatial_index_t* idx,
                                      int universe_id,
                                      double lx,
                                      double ly,
                                      double lz,
                                      int depth,
                                      const alea_matrix_t* transform,
                                      int lattice_cell_index,
                                      const alea_matrix_t* lattice_transform,
                                      alea_hier_cell_hit_t* out_hit,
                                      int* found,
                                      hier_path_cache_t* path,
                                      hier_path_cache_t* best_path) {
    if (depth >= HIER_MAX_PLACEMENT_DEPTH) return 0;

    const alea_universe_t* univ = alea_get_universe(sys, universe_id);
    if (!univ) return -1;

    const hier_universe_blas_t* blas = find_blas(idx, universe_id);
    if (blas && blas->built && blas->node_count > 0) {
        hier_cand_buf_t candidates;
        hier_cand_init(&candidates);
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
        query_blas_node(&ctx, 0);
        if (!ctx.error && candidates.count > 1) {
            hier_cand_sort(&candidates);
        }
        /* cell_pos order is universe definition (deck) order, and the first
         * containing cell wins — the same precedence the canonical resolver
         * (find_cell_recursive) applies when cells overlap. */
        for (size_t i = 0; !ctx.error && i < candidates.count; i++) {
            uint32_t cell_pos = candidates.data[i];
            if (cell_pos >= blas->cell_count) {
                ctx.error = -1;
                break;
            }
            uint32_t candidate_cell = blas->cells[cell_pos].cell_index;
            int rc = hier_find_deepest_cell(sys, idx, candidate_cell,
                                            lx, ly, lz, depth,
                                            transform,
                                            lattice_cell_index,
                                            lattice_transform,
                                            out_hit, found,
                                            path, best_path);
            if (rc < 0) ctx.error = -1;
            if (rc > 0) break;
        }
        hier_cand_free(&candidates);
        return ctx.error ? -1 : 0;
    }

    for (size_t i = 0; i < univ->cell_indices.count; i++) {
        uint32_t candidate_cell = (uint32_t)univ->cell_indices.data[i];
        int rc = hier_find_deepest_cell(sys, idx, candidate_cell,
                                        lx, ly, lz, depth,
                                        transform,
                                        lattice_cell_index,
                                        lattice_transform,
                                        out_hit, found,
                                        path, best_path);
        if (rc < 0) return -1;
        if (rc > 0) break;
    }

    return 0;
}

int alea_hier_spatial_check_cached_containment(alea_system_t* sys,
                                               uint32_t cell_index,
                                               double x, double y, double z) {
    return alea_hier_spatial_get_cached_cell_state(sys, cell_index, x, y, z,
                                                   NULL, NULL, NULL);
}

int alea_hier_spatial_get_cached_cell_state(alea_system_t* sys,
                                            uint32_t cell_index,
                                            double x, double y, double z,
                                            alea_matrix_t* out_transform,
                                            int* out_lattice_cell_index,
                                            alea_matrix_t* out_lattice_transform) {
    if (!sys) return -1;

    uint64_t generation = alea_system_geometry_generation(sys);
    if (g_hier_cache_system != sys ||
        g_hier_cache_generation != generation ||
        g_hier_cache_system_id != sys->system_id) {
        return -1;
    }

    /* Locate the target entry. Return -1 (unknown, caller falls back to a
     * full lookup) if the cell is not on the last query's cached path. */
    int target = -1;
    for (int i = 0; i < g_hier_cache_count; i++) {
        if (!g_hier_cache[i].valid) break;
        if (g_hier_cache[i].cell_index == cell_index) { target = i; break; }
    }
    if (target < 0) return -1;

    /* Validate the FULL ancestor chain. A deep cell's CSG can easily extend
     * outside its parent FILL region; testing the target cell alone gives
     * false positives when the point now sits in a different parent chain.
     * The chain check mirrors hier_cache_try semantics. */
    int lattice_cell_index = -1;
    alea_matrix_t lattice_transform;
    alea_matrix_identity(&lattice_transform);
    for (int i = 0; i <= target; i++) {
        hier_cached_cell_t* ent = &g_hier_cache[i];
        double lx = x, ly = y, lz = z;
        alea_matrix_transform_point_inverse(&ent->transform, &lx, &ly, &lz);

        const alea_cell_entry_t* cell = &sys->cells.data[ent->cell_index];

        if (ent->is_lattice) {
            alea_lattice_location_t location;
            if (alea_lattice_locate_point(sys, cell, lx, ly, lz, &location) != 1 ||
                location.fill_universe != ent->lat_fill_universe) return 0;
            if (location.i != ent->lat_i || location.j != ent->lat_j ||
                location.k != ent->lat_k) return 0;
            if (location.ox != ent->lat_ox || location.oy != ent->lat_oy ||
                location.oz != ent->lat_oz) {
                return 0;
            }
            lattice_cell_index = (int)ent->cell_index;
            lattice_transform = ent->transform;
            continue;
        }

        if (!alea_contains_point(sys, cell->root_node_id, lx, ly, lz))
            return 0;
    }

    if (out_transform) {
        *out_transform = g_hier_cache[target].transform;
    }
    if (out_lattice_cell_index) {
        *out_lattice_cell_index = lattice_cell_index;
    }
    if (out_lattice_transform) {
        *out_lattice_transform = lattice_transform;
    }
    return 1;
}

int alea_hier_spatial_check_path_containment(alea_system_t* sys,
                                             const alea_hier_ray_path_t* path,
                                             int target_entry,
                                             double x,
                                             double y,
                                             double z,
                                             alea_matrix_t* out_transform,
                                             int* out_lattice_cell_index,
                                             alea_matrix_t* out_lattice_transform) {
    if (!sys || !path || target_entry < 0 || target_entry >= path->count)
        return -1;

    int lattice_cell_index = -1;
    alea_matrix_t lattice_transform;
    alea_matrix_identity(&lattice_transform);
    for (int i = 0; i <= target_entry; i++) {
        const alea_hier_ray_path_entry_t* ent = &path->entries[i];
        if ((size_t)ent->cell_index >= alea_vec_count(&sys->cells)) {
            return -1;
        }

        double lx = x;
        double ly = y;
        double lz = z;
        alea_matrix_transform_point_inverse(&ent->transform, &lx, &ly, &lz);

        const alea_cell_entry_t* cell = &sys->cells.data[ent->cell_index];
        if (ent->is_lattice) {
            alea_lattice_location_t location;
            if (alea_lattice_locate_point(sys, cell, lx, ly, lz, &location) != 1 ||
                location.fill_universe != ent->lat_fill_universe) return 0;
            if (location.i != ent->lat_i || location.j != ent->lat_j ||
                location.k != ent->lat_k) return 0;
            if (location.ox != ent->lat_ox || location.oy != ent->lat_oy ||
                location.oz != ent->lat_oz) {
                return 0;
            }
            lattice_cell_index = (int)ent->cell_index;
            lattice_transform = ent->transform;
            continue;
        }

        if (cell->root_node_id == ALEA_NODE_ID_INVALID ||
            !alea_contains_point(sys, cell->root_node_id, lx, ly, lz)) {
            return 0;
        }
    }

    if (out_transform) {
        *out_transform = path->entries[target_entry].transform;
    }
    if (out_lattice_cell_index) {
        *out_lattice_cell_index = lattice_cell_index;
    }
    if (out_lattice_transform) {
        *out_lattice_transform = lattice_transform;
    }
    return 1;
}

static int hier_spatial_find_cells_at_point_impl(alea_system_t* sys,
                                                 double x,
                                                 double y,
                                                 double z,
                                                 alea_cell_hit_t* out_hits,
                                                 size_t max_hits,
                                                 int use_cache) {
    if (!sys || !out_hits || max_hits == 0) return -1;

    if (!sys->hier_spatial_index || !sys->hier_spatial_index->built) {
        if (alea_hier_spatial_index_build(sys) != 0) {
            return -1;
        }
    }

    /* Invalidate cache if the system or its geometry changed. */
    uint64_t generation = alea_system_geometry_generation(sys);
    if (g_hier_cache_system != sys ||
        g_hier_cache_generation != generation ||
        g_hier_cache_system_id != sys->system_id) {
        hier_cache_invalidate();
        g_hier_cache_system = sys;
        g_hier_cache_generation = generation;
        g_hier_cache_system_id = sys->system_id;
    }

    alea_hier_spatial_index_t* idx = sys->hier_spatial_index;

    if (use_cache) {
        int cached = hier_cache_try(sys, x, y, z, out_hits, max_hits);
        if (cached > 0) {
            return cached;
        }
    }

    /* Cache miss — wipe and repopulate during the full descent. */
    hier_cache_invalidate();

    alea_matrix_t identity;
    alea_matrix_identity(&identity);

    size_t hit_count = 0;
    int rc = hier_query_universe(sys, idx, 0, x, y, z, 0,
                                 &identity,
                                 out_hits, max_hits, &hit_count);
    if (rc < 0) {
        hier_cache_invalidate();
        return -1;
    }

    if (!use_cache) {
        hier_cache_invalidate();
    }

    return (int)hit_count;
}

int alea_hier_spatial_find_cells_at_point(alea_system_t* sys,
                                          double x,
                                          double y,
                                          double z,
                                          alea_cell_hit_t* out_hits,
                                          size_t max_hits) {
    return hier_spatial_find_cells_at_point_impl(sys, x, y, z,
                                                out_hits, max_hits, 1);
}

int alea_hier_spatial_find_cells_at_point_uncached(alea_system_t* sys,
                                                   double x,
                                                   double y,
                                                   double z,
                                                   alea_cell_hit_t* out_hits,
                                                   size_t max_hits) {
    return hier_spatial_find_cells_at_point_impl(sys, x, y, z,
                                                out_hits, max_hits, 0);
}

int alea_hier_spatial_find_path_at_point(alea_system_t* sys,
                                         double x,
                                         double y,
                                         double z,
                                         alea_hier_cell_hit_t* out_hit,
                                         alea_hier_ray_path_t* out_path) {
    if (!sys || !out_hit) return -1;
    if (out_path) out_path->count = 0;

    if (!sys->hier_spatial_index || !sys->hier_spatial_index->built) {
        if (alea_hier_spatial_index_build(sys) != 0) {
            return -1;
        }
    }

    alea_hier_spatial_index_t* idx = sys->hier_spatial_index;

    uint64_t generation = alea_system_geometry_generation(sys);
    if (g_hier_cache_system != sys ||
        g_hier_cache_generation != generation ||
        g_hier_cache_system_id != sys->system_id) {
        hier_cache_invalidate();
        g_hier_cache_system = sys;
        g_hier_cache_generation = generation;
        g_hier_cache_system_id = sys->system_id;
    }
    hier_cache_invalidate();

    alea_matrix_t identity;
    alea_matrix_identity(&identity);

    memset(out_hit, 0, sizeof(*out_hit));
    out_hit->lattice_cell_index = -1;
    alea_matrix_identity(&out_hit->lattice_transform);
    int found = 0;
    hier_path_cache_t path;
    hier_path_cache_t best_path;
    hier_path_clear(&path);
    hier_path_clear(&best_path);
    int rc = hier_find_deepest_universe(sys, idx, 0, x, y, z, 0,
                                        &identity, -1, NULL, out_hit, &found,
                                        &path, &best_path);
    if (rc < 0) {
        hier_cache_invalidate();
        return -1;
    }

    if (found) {
        hier_path_commit_to_cache(&best_path);
        hier_path_export(&best_path, out_path);
    } else {
        hier_cache_invalidate();
    }

    return found ? 1 : 0;
}

int alea_hier_spatial_find_deepest_cell_at_point(alea_system_t* sys,
                                                 double x,
                                                 double y,
                                                 double z,
                                                 alea_hier_cell_hit_t* out_hit) {
    return alea_hier_spatial_find_path_at_point(sys, x, y, z, out_hit, NULL);
}

/* Containment plus deck-order selection are distinct operations: callers
 * restarting a cached prefix must prove that its occurrence still owns the
 * point, not merely that its CSG remains true in overlapping geometry. */
static int hier_path_check_canonical_owners(alea_system_t* sys,
                                            const alea_hier_ray_path_t* path,
                                            int last_entry,
                                            double x, double y, double z) {
    if (!sys || !path || last_entry < 0 || last_entry >= path->count) return -1;
    for (int i = 0; i <= last_entry; i++) {
        const alea_hier_ray_path_entry_t* entry = &path->entries[i];
        if ((size_t)entry->cell_index >= alea_vec_count(&sys->cells)) return -1;
        double lx = x, ly = y, lz = z;
        alea_matrix_transform_point_inverse(&entry->transform, &lx, &ly, &lz);
        int owner = alea_hier_spatial_find_ordered_cell_in_universe(
            sys, entry->universe_id, lx, ly, lz, -1);
        if (owner == -2) return -1;
        if (owner != (int)entry->cell_index) return 0;
    }
    return 1;
}

int alea_hier_path_enter_lattice_location(
    alea_system_t* sys,
    const alea_hier_ray_path_t* path,
    int lattice_entry,
    double x, double y, double z,
    const alea_lattice_location_t* location,
    alea_hier_cell_hit_t* out_hit,
    alea_hier_ray_path_t* out_path) {
    if (!sys || !path || !location || !out_hit || lattice_entry < 0 ||
        lattice_entry >= path->count) return -1;
    if (out_path) out_path->count = 0;

    if (!sys->hier_spatial_index || !sys->hier_spatial_index->built) {
        if (alea_hier_spatial_index_build(sys) != 0) return -1;
    }
    const alea_hier_ray_path_entry_t* parent = &path->entries[lattice_entry];
    if ((size_t)parent->cell_index >= alea_vec_count(&sys->cells) ||
        !parent->is_lattice) return -1;

    /* Ancestors retain their exact old placements; the lattice entry itself
     * is intentionally reselected below. */
    if (lattice_entry > 0) {
        int valid = alea_hier_spatial_check_path_containment(
            sys, path, lattice_entry - 1, x, y, z, NULL, NULL, NULL);
        if (valid <= 0) return valid;
        valid = hier_path_check_canonical_owners(sys, path,
                                                  lattice_entry - 1,
                                                  x, y, z);
        if (valid <= 0) return valid;
    }

    double parent_lx = x, parent_ly = y, parent_lz = z;
    alea_matrix_transform_point_inverse(&parent->transform,
                                        &parent_lx, &parent_ly, &parent_lz);
    int owner = alea_hier_spatial_find_ordered_cell_in_universe(
        sys, parent->universe_id, parent_lx, parent_ly, parent_lz, -1);
    if (owner == -2) return -1;
    if (owner != (int)parent->cell_index) return 0;

    const alea_cell_entry_t* cell = &sys->cells.data[parent->cell_index];
    alea_lattice_location_t actual;
    if (alea_lattice_locate_point(sys, cell, parent_lx, parent_ly, parent_lz,
                                  &actual) != 1) return 0;
    if (actual.fill_universe != location->fill_universe ||
        actual.i != location->i || actual.j != location->j ||
        actual.k != location->k || actual.ox != location->ox ||
        actual.oy != location->oy || actual.oz != location->oz) {
        return 0;
    }

    alea_matrix_t element_translation;
    translation_matrix(&element_translation, location->ox, location->oy,
                       location->oz);
    alea_matrix_t transform;
    alea_matrix_multiply(&transform, &parent->transform, &element_translation);
    if (!alea_matrix_invert(&transform)) return -1;

    hier_path_cache_t prefix;
    hier_path_cache_t best_path;
    hier_path_import_prefix(path, lattice_entry + 1, &prefix);
    hier_path_import_prefix(path, lattice_entry + 1, &best_path);
    hier_cached_cell_t* updated = &prefix.cells[lattice_entry];
    updated->fill_universe = location->fill_universe;
    updated->lat_fill_universe = location->fill_universe;
    updated->lat_i = location->i;
    updated->lat_j = location->j;
    updated->lat_k = location->k;
    updated->lat_ox = location->ox;
    updated->lat_oy = location->oy;
    updated->lat_oz = location->oz;
    best_path.cells[lattice_entry] = *updated;

    memset(out_hit, 0, sizeof(*out_hit));
    out_hit->lattice_cell_index = -1;
    alea_matrix_identity(&out_hit->lattice_transform);
    int found = 0;
    int rc = hier_find_deepest_universe(
        sys, sys->hier_spatial_index, location->fill_universe,
        parent_lx - location->ox, parent_ly - location->oy,
        parent_lz - location->oz, parent->depth + 1, &transform,
        (int)parent->cell_index, &parent->transform, out_hit, &found,
        &prefix, &best_path);
    if (rc < 0) return -1;
    if (!found) return 0;

    hier_path_commit_to_cache(&best_path);
    hier_path_export(&best_path, out_path);
    return 1;
}

int alea_hier_spatial_find_path_from_parent(alea_system_t* sys,
                                            const alea_hier_ray_path_t* path,
                                            int parent_entry,
                                            double x,
                                            double y,
                                            double z,
                                            alea_hier_cell_hit_t* out_hit,
                                            alea_hier_ray_path_t* out_path) {
    if (!sys || !path || !out_hit || parent_entry < -1) return -1;
    if (out_path) out_path->count = 0;

    if (!sys->hier_spatial_index || !sys->hier_spatial_index->built) {
        if (alea_hier_spatial_index_build(sys) != 0) {
            return -1;
        }
    }
    alea_hier_spatial_index_t* idx = sys->hier_spatial_index;

    int universe_id = 0;
    int depth = 0;
    double lx = x;
    double ly = y;
    double lz = z;
    alea_matrix_t transform;
    alea_matrix_identity(&transform);
    const alea_hier_ray_path_entry_t* parent = NULL;

    if (parent_entry >= 0) {
        if (parent_entry >= path->count) return -1;
        parent = &path->entries[parent_entry];
        if ((size_t)parent->cell_index >= alea_vec_count(&sys->cells)) {
            return -1;
        }

        if (parent->is_lattice) {
            const alea_cell_entry_t* cell =
                &sys->cells.data[parent->cell_index];
            double parent_lx = x, parent_ly = y, parent_lz = z;
            alea_matrix_transform_point_inverse(&parent->transform,
                                                &parent_lx, &parent_ly,
                                                &parent_lz);
            alea_lattice_location_t location;
            if (alea_lattice_locate_point(sys, cell, parent_lx, parent_ly,
                                          parent_lz, &location) != 1) {
                return 0;
            }
            return alea_hier_path_enter_lattice_location(
                sys, path, parent_entry, x, y, z, &location,
                out_hit, out_path);
        }

        /* Lattice parents returned above through the shared location
         * primitive. Ordinary fills retain and validate the parent itself. */
        int valid = alea_hier_spatial_check_path_containment(
            sys, path, parent_entry, x, y, z, NULL, NULL, NULL);
        if (valid <= 0) return valid;
        valid = hier_path_check_canonical_owners(sys, path, parent_entry,
                                                  x, y, z);
        if (valid <= 0) return valid;

        const alea_cell_entry_t* cell = &sys->cells.data[parent->cell_index];
        double parent_lx = x;
        double parent_ly = y;
        double parent_lz = z;
        alea_matrix_transform_point_inverse(&parent->transform,
                                            &parent_lx, &parent_ly, &parent_lz);

        if (cell->fill_universe <= 0) return 0;
        const alea_matrix_t* fill_transform =
            cell_fill_matrix_cached(idx, parent->cell_index);
        if (!fill_transform) return -1;

        lx = parent_lx;
        ly = parent_ly;
        lz = parent_lz;
        alea_matrix_transform_point_inverse(fill_transform, &lx, &ly, &lz);

        alea_matrix_multiply(&transform, &parent->transform, fill_transform);
        if (!alea_matrix_invert(&transform)) return -1;

        universe_id = cell->fill_universe;
        depth = parent->depth + 1;
    }

    memset(out_hit, 0, sizeof(*out_hit));
    out_hit->lattice_cell_index = -1;
    alea_matrix_identity(&out_hit->lattice_transform);

    hier_path_cache_t prefix;
    hier_path_cache_t best_path;
    hier_path_import_prefix(path, parent_entry + 1, &prefix);
    hier_path_import_prefix(path, parent_entry + 1, &best_path);

    int found = 0;
    int rc = hier_find_deepest_universe(sys, idx, universe_id,
                                        lx, ly, lz, depth,
                                        &transform,
                                        -1, NULL,
                                        out_hit, &found,
                                        &prefix, &best_path);
    if (rc < 0) return -1;
    if (!found) return 0;

    hier_path_commit_to_cache(&best_path);
    hier_path_export(&best_path, out_path);
    return 1;
}

void alea_hier_coherence_state_clear(alea_hier_coherence_state_t* state) {
    if (!state) return;
    memset(state, 0, sizeof(*state));
    state->path.count = 0;
    state->deepest.lattice_cell_index = -1;
    alea_matrix_identity(&state->deepest.lattice_transform);
}

/* Reconstruct only the point-dependent pieces of the deepest result after a
 * complete path reuse.  The concrete placement transform and resolution
 * status stay in the path, while local coordinates must describe this point. */
static int hier_coherence_hit_from_path(alea_system_t* sys,
                                        const alea_hier_ray_path_t* path,
                                        double x, double y, double z,
                                        alea_hier_cell_hit_t* out_hit) {
    if (!sys || !path || !out_hit || path->count <= 0 ||
        path->count > ALEA_HIER_RAY_PATH_MAX) return -1;
    const alea_hier_ray_path_entry_t* leaf = &path->entries[path->count - 1];
    if ((size_t)leaf->cell_index >= alea_vec_count(&sys->cells)) return -1;
    const alea_cell_entry_t* cell = &sys->cells.data[leaf->cell_index];

    memset(out_hit, 0, sizeof(*out_hit));
    out_hit->hit.cell_id = leaf->cell_id;
    out_hit->hit.cell_index = (int)leaf->cell_index;
    out_hit->hit.material_id = leaf->material_id;
    out_hit->hit.universe_id = leaf->universe_id;
    out_hit->hit.fill_universe = leaf->fill_universe;
    out_hit->hit.depth = leaf->depth;
    out_hit->hit.resolution_flags = alea_cell_entry_is_container(cell)
        ? ALEA_RESOLVE_UNDEFINED_FILL : 0;
    out_hit->transform = leaf->transform;
    out_hit->hit.local_x = x;
    out_hit->hit.local_y = y;
    out_hit->hit.local_z = z;
    alea_matrix_transform_point_inverse(&leaf->transform,
                                        &out_hit->hit.local_x,
                                        &out_hit->hit.local_y,
                                        &out_hit->hit.local_z);
    out_hit->lattice_cell_index = -1;
    alea_matrix_identity(&out_hit->lattice_transform);
    for (int i = path->count - 1; i >= 0; i--) {
        const alea_hier_ray_path_entry_t* entry = &path->entries[i];
        if (entry->is_lattice) {
            out_hit->lattice_cell_index = (int)entry->cell_index;
            out_hit->lattice_transform = entry->transform;
            break;
        }
    }
    return 0;
}

static void hier_coherence_store(alea_system_t* sys,
                                 alea_hier_coherence_state_t* state,
                                 const alea_hier_cell_hit_t* hit,
                                 uint64_t generation) {
    state->deepest = *hit;
    state->system_id = sys->system_id;
    state->geometry_generation = generation;
    state->complete = 1;
    state->overflowed = 0;
}

int alea_hier_spatial_resolve_coherent(
    alea_system_t* sys,
    double x, double y, double z,
    const alea_hier_coherence_state_t* previous,
    alea_hier_coherence_state_t* current,
    alea_hier_cell_hit_t* out_hit,
    alea_hier_coherence_kind_t* out_kind) {
    if (!sys || !current || !out_hit) return -1;
    alea_hier_coherence_state_clear(current);
    if (out_kind) *out_kind = ALEA_HIER_COH_ROOT_QUERY;

    const uint64_t generation = alea_system_geometry_generation(sys);
    const bool has_previous = previous && previous->complete &&
        previous->path.count > 0 &&
        previous->path.count <= ALEA_HIER_RAY_PATH_MAX &&
        previous->system_id == sys->system_id &&
        previous->geometry_generation == generation;

    if (has_previous &&
        !(previous->deepest.hit.resolution_flags & ALEA_RESOLVE_UNDEFINED_FILL)) {
        int valid = alea_hier_spatial_check_path_containment(
            sys, &previous->path, previous->path.count - 1, x, y, z,
            NULL, NULL, NULL);
        if (valid > 0)
            valid = hier_path_check_canonical_owners(
                sys, &previous->path, previous->path.count - 1, x, y, z);
        if (valid > 0 && hier_coherence_hit_from_path(
                             sys, &previous->path, x, y, z, out_hit) == 0) {
            current->path.count = previous->path.count;
            memcpy(current->path.entries, previous->path.entries,
                   (size_t)previous->path.count *
                       sizeof(current->path.entries[0]));
            hier_coherence_store(sys, current, out_hit, generation);
            if (out_kind) *out_kind = ALEA_HIER_COH_PATH_REUSED;
            return 1;
        }
    }

    if (has_previous) {
        for (int parent = previous->path.count - 2; parent >= 0; parent--) {
            int rc = alea_hier_spatial_find_path_from_parent(
                sys, &previous->path, parent, x, y, z, out_hit,
                &current->path);
            if (rc < 0) return -1;
            if (rc == 0) continue;
            hier_coherence_store(sys, current, out_hit, generation);
            if (out_kind) {
                const alea_hier_ray_path_entry_t* old_parent =
                    &previous->path.entries[parent];
                const alea_hier_ray_path_entry_t* new_parent =
                    parent < current->path.count
                        ? &current->path.entries[parent] : NULL;
                *out_kind = old_parent->is_lattice && new_parent &&
                    (old_parent->lat_fill_universe != new_parent->lat_fill_universe ||
                     old_parent->lat_i != new_parent->lat_i ||
                     old_parent->lat_j != new_parent->lat_j ||
                     old_parent->lat_k != new_parent->lat_k)
                    ? ALEA_HIER_COH_LATTICE_TRANSITION
                    : ALEA_HIER_COH_PREFIX_RESTART;
            }
            return 1;
        }
    }

    int rc = alea_hier_spatial_find_path_at_point(sys, x, y, z, out_hit,
                                                   &current->path);
    if (rc <= 0) return rc;
    hier_coherence_store(sys, current, out_hit, generation);
    if (out_kind)
        *out_kind = has_previous ? ALEA_HIER_COH_FULL_FALLBACK
                                 : ALEA_HIER_COH_ROOT_QUERY;
    return 1;
}

const alea_matrix_t*
alea_hier_spatial_get_cell_fill_matrix(const alea_system_t* sys,
                                       uint32_t cell_index) {
    if (!sys || !sys->hier_spatial_index) return NULL;
    return cell_fill_matrix_cached(sys->hier_spatial_index, cell_index);
}

static bool hier_cell_contains_local(const alea_system_t* sys,
                                     const alea_cell_entry_t* cell,
                                     double lx, double ly, double lz);

int alea_hier_spatial_find_cell_in_universe(alea_system_t* sys,
                                            int universe_id,
                                            double lx,
                                            double ly,
                                            double lz) {
    if (!sys || !sys->hier_spatial_index) return -2;
    alea_hier_spatial_index_t* idx = sys->hier_spatial_index;
    if (!idx->built) return -2;

    const hier_universe_blas_t* blas = find_blas(idx, universe_id);
    if (blas && blas->built && blas->node_count > 0) {
        hier_cand_buf_t candidates;
        hier_cand_init(&candidates);
        hier_point_query_ctx_t ctx = {
            .sys = sys, .idx = idx, .blas = blas,
            .x = lx, .y = ly, .z = lz,
            .candidates = &candidates, .error = 0
        };
        query_blas_node(&ctx, 0);
        if (ctx.error) { hier_cand_free(&candidates); return -2; }

        /* Sort candidates by bbox volume, smallest first. Cells in a universe
         * should not overlap, so we only need to find the first containing
         * cell — tight bboxes are far more likely to actually contain the
         * point than world-sized "shell" bboxes that just happen to overlap. */
        if (candidates.count > 1) {
            hier_cand_sort_by_volume(&candidates, blas);
        }

        int result = -1;
        for (size_t i = 0; i < candidates.count; i++) {
            uint32_t cell_pos = candidates.data[i];
            if (cell_pos >= blas->cell_count) continue;
            uint32_t cell_index = blas->cells[cell_pos].cell_index;
            const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
            if (hier_cell_contains_local(sys, cell, lx, ly, lz)) {
                result = (int)cell_index;
                break;
            }
        }
        hier_cand_free(&candidates);
        return result;
    }

    /* No BLAS — fall back to linear scan over the universe. */
    const alea_universe_t* univ = alea_get_universe(sys, universe_id);
    if (!univ) return -1;
    for (size_t i = 0; i < univ->cell_indices.count; i++) {
        uint32_t cell_index = (uint32_t)univ->cell_indices.data[i];
        const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
        if (hier_cell_contains_local(sys, cell, lx, ly, lz)) {
            return (int)cell_index;
        }
    }
    return -1;
}

static bool hier_cell_contains_local(const alea_system_t* sys,
                                     const alea_cell_entry_t* cell,
                                     double lx, double ly, double lz) {
    if (cell->lat_type != 0 && cell->lat_fill) {
        alea_lattice_location_t location;
        return alea_lattice_locate_point(sys, cell, lx, ly, lz,
                                         &location) == 1;
    }
    return cell->root_node_id != ALEA_NODE_ID_INVALID &&
        alea_contains_point(sys, cell->root_node_id, lx, ly, lz);
}

/* True if the cell's CSG region references the given MCNP surface id. Used to
 * restrict the expensive containment test to cells sharing a just-crossed
 * boundary. */
static bool cell_references_mc_surface(const alea_system_t* sys,
                                       const alea_cell_entry_t* cell,
                                       int mc_surface_id) {
    if (!cell->surface_indices || mc_surface_id <= 0) return false;
    for (size_t i = 0; i < cell->surface_index_count; i++) {
        uint32_t si = cell->surface_indices[i];
        if (si >= alea_vec_count(&sys->surfaces)) continue;
        if (sys->surfaces.data[si].mc_surface_id == mc_surface_id) return true;
    }
    return false;
}

int alea_hier_spatial_find_ordered_cell_in_universe(alea_system_t* sys,
                                                    int universe_id,
                                                    double lx,
                                                    double ly,
                                                    double lz,
                                                    int crossed_mc_surface_id) {
    if (!sys || !sys->hier_spatial_index) return -2;
    alea_hier_spatial_index_t* idx = sys->hier_spatial_index;
    if (!idx->built) return -2;

    const hier_universe_blas_t* blas = find_blas(idx, universe_id);
    if (blas && blas->built && blas->node_count > 0) {
        hier_cand_buf_t candidates;
        hier_cand_init(&candidates);
        hier_point_query_ctx_t ctx = {
            .sys = sys, .idx = idx, .blas = blas,
            .x = lx, .y = ly, .z = lz,
            .candidates = &candidates, .error = 0
        };
        query_blas_node(&ctx, 0);
        if (ctx.error) {
            hier_cand_free(&candidates);
            return -2;
        }
        if (candidates.count > 1) {
            hier_cand_sort(&candidates);
        }

        int result = -1;
        /* Fast path: the cell across a crossed boundary shares that surface, so
         * test containment only on candidates referencing it (typically tens vs
         * thousands of bbox candidates). */
        if (crossed_mc_surface_id > 0) {
            for (size_t i = 0; i < candidates.count; i++) {
                uint32_t cell_pos = candidates.data[i];
                if (cell_pos >= blas->cell_count) continue;
                uint32_t cell_index = blas->cells[cell_pos].cell_index;
                const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
                if (cell->root_node_id == ALEA_NODE_ID_INVALID &&
                    !(cell->lat_type != 0 && cell->lat_fill)) continue;
                if (!cell_references_mc_surface(sys, cell, crossed_mc_surface_id))
                    continue;
                if (hier_cell_contains_local(sys, cell, lx, ly, lz)) {
                    /* First containing cell in deck order wins, matching the
                     * canonical resolver's precedence for overlapping cells. */
                    result = (int)cell_index;
                    break;
                }
            }
        }
        /* No surface hint, or no referencing candidate contained the point
         * (coincident surfaces, etc.): exhaustive scan. */
        if (result < 0) {
            for (size_t i = 0; i < candidates.count; i++) {
                uint32_t cell_pos = candidates.data[i];
                if (cell_pos >= blas->cell_count) continue;
                uint32_t cell_index = blas->cells[cell_pos].cell_index;
                const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
                if (hier_cell_contains_local(sys, cell, lx, ly, lz)) {
                    result = (int)cell_index;
                    break;
                }
            }
        }
        hier_cand_free(&candidates);
        return result;
    }

    const alea_universe_t* univ = alea_get_universe(sys, universe_id);
    if (!univ) return -1;
    int result = -1;
    for (size_t i = 0; i < univ->cell_indices.count; i++) {
        uint32_t cell_index = (uint32_t)univ->cell_indices.data[i];
        const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
        if (hier_cell_contains_local(sys, cell, lx, ly, lz)) {
            result = (int)cell_index;
            break;
        }
    }
    return result;
}

static int compare_spatial_hits_by_depth_cell(const void* a, const void* b);

static int append_universe_region_hit(alea_system_t* sys,
                                      int universe_id,
                                      int depth,
                                      uint32_t cell_index,
                                      uint32_t synthetic_index,
                                      alea_spatial_hit_t* out_hits,
                                      size_t max_hits,
                                      size_t* hit_count) {
    if (*hit_count >= max_hits) return 0;
    const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
    alea_spatial_hit_t* hit = &out_hits[*hit_count];
    hit->instance_index = synthetic_index;
    hit->cell_index = cell_index;
    hit->cell_id = cell->mc_cell_id;
    hit->material_id = cell->material_id;
    hit->universe_id = universe_id;
    hit->depth = depth;
    hit->is_terminal = (cell->fill_universe <= 0 &&
                        !(cell->lat_type != 0 && cell->lat_fill));
    alea_matrix_identity(&hit->transform);
    (*hit_count)++;
    return 0;
}

static int query_universe_region_blas_node(alea_system_t* sys,
                                           const hier_universe_blas_t* blas,
                                           uint32_t node_index,
                                           int universe_id,
                                           const alea_bbox_t* local_bbox,
                                           alea_spatial_hit_t* out_hits,
                                           size_t max_hits,
                                           size_t* hit_count) {
    if (node_index >= blas->node_count) return -1;
    if (*hit_count >= max_hits) return 0;

    const hier_bvh_node_t* node = &blas->nodes[node_index];
    if (!hier_fbbox_intersects_dbbox(&node->bbox, local_bbox)) return 0;

    if (node->count == 0) {
        if (query_universe_region_blas_node(sys, blas, node->left_or_first,
                                            universe_id, local_bbox,
                                            out_hits, max_hits,
                                            hit_count) != 0) {
            return -1;
        }
        return query_universe_region_blas_node(sys, blas, node->right_child,
                                               universe_id, local_bbox,
                                               out_hits, max_hits, hit_count);
    }

    for (uint16_t i = 0; i < node->count; i++) {
        if (*hit_count >= max_hits) return 0;
        uint32_t cell_pos = blas->indices[node->left_or_first + i];
        if (cell_pos >= blas->cell_count) return -1;
        const hier_blas_cell_t* blas_cell = &blas->cells[cell_pos];
        if (!hier_fbbox_intersects_dbbox(&blas_cell->bbox, local_bbox)) continue;
        if (append_universe_region_hit(sys, universe_id, 0,
                                       blas_cell->cell_index, cell_pos,
                                       out_hits, max_hits, hit_count) != 0) {
            return -1;
        }
    }
    return 0;
}

int alea_hier_spatial_query_universe_region(alea_system_t* sys,
                                            int universe_id,
                                            const alea_bbox_t* local_bbox,
                                            alea_spatial_hit_t* out_hits,
                                            size_t max_hits) {
    if (!sys || !local_bbox || !out_hits || max_hits == 0) return -1;
    if (!sys->hier_spatial_index || !sys->hier_spatial_index->built) {
        if (alea_hier_spatial_index_build(sys) != 0) return -1;
    }

    alea_hier_spatial_index_t* idx = sys->hier_spatial_index;
    size_t hit_count = 0;
    const hier_universe_blas_t* blas = find_blas(idx, universe_id);
    if (blas && blas->built && blas->node_count > 0) {
        int rc = query_universe_region_blas_node(sys, blas, 0,
                                                 universe_id, local_bbox,
                                                 out_hits, max_hits,
                                                 &hit_count);
        if (rc != 0) return -1;
    } else {
        const alea_universe_t* univ = alea_get_universe(sys, universe_id);
        if (!univ) return 0;
        for (size_t i = 0; i < univ->cell_indices.count; i++) {
            if (hit_count >= max_hits) break;
            uint32_t cell_index = (uint32_t)univ->cell_indices.data[i];
            const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
            alea_bbox_t bbox = (cell->lat_type != 0 && cell->lat_fill)
                ? lattice_container_bbox(cell)
                : local_cell_bbox(sys, cell_index);
            if (!bbox_intersects_local(&bbox, local_bbox)) continue;
            if (append_universe_region_hit(sys, universe_id, 0,
                                           cell_index, (uint32_t)i,
                                           out_hits, max_hits,
                                           &hit_count) != 0) {
                return -1;
            }
        }
    }

    if (hit_count > 1) {
        qsort(out_hits, hit_count, sizeof(*out_hits),
              compare_spatial_hits_by_depth_cell);
    }
    return (int)hit_count;
}

static void insert_ray_candidate_sorted(alea_system_t* sys,
                                        uint32_t cell_index,
                                        double t_enter,
                                        double t_exit,
                                        alea_hier_ray_candidate_t* out_hits,
                                        size_t max_hits,
                                        size_t* hit_count) {
    if (max_hits == 0) return;
    size_t pos = *hit_count;
    if (pos >= max_hits) {
        if (t_enter > out_hits[max_hits - 1].t_enter ||
            (t_enter == out_hits[max_hits - 1].t_enter &&
             t_exit >= out_hits[max_hits - 1].t_exit)) {
            return;
        }
        pos = max_hits - 1;
    } else {
        (*hit_count)++;
    }

    while (pos > 0) {
        const alea_hier_ray_candidate_t* prev = &out_hits[pos - 1];
        if (prev->t_enter < t_enter ||
            (prev->t_enter == t_enter && prev->t_exit <= t_exit)) {
            break;
        }
        out_hits[pos] = out_hits[pos - 1];
        pos--;
    }

    const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
    out_hits[pos].cell_index = cell_index;
    out_hits[pos].cell_id = cell->mc_cell_id;
    out_hits[pos].t_enter = t_enter;
    out_hits[pos].t_exit = t_exit;
}

static int query_universe_ray_blas_node(alea_system_t* sys,
                                        const hier_universe_blas_t* blas,
                                        uint32_t node_index,
                                        double ox,
                                        double oy,
                                        double oz,
                                        double inv_dx,
                                        double inv_dy,
                                        double inv_dz,
                                        double t_min,
                                        double t_max,
                                        alea_hier_ray_candidate_t* out_hits,
                                        size_t max_hits,
                                        size_t* hit_count) {
    if (node_index >= blas->node_count) return -1;

    const hier_bvh_node_t* node = &blas->nodes[node_index];
    double node_enter, node_exit;
    if (!hier_fbbox_ray_enter_exit(&node->bbox, ox, oy, oz,
                                   inv_dx, inv_dy, inv_dz,
                                   t_min, t_max,
                                   &node_enter, &node_exit)) {
        return 0;
    }

    if (node->count == 0) {
        if (node->left_or_first >= blas->node_count ||
            node->right_child >= blas->node_count) {
            return -1;
        }
        const hier_bvh_node_t* left = &blas->nodes[node->left_or_first];
        const hier_bvh_node_t* right = &blas->nodes[node->right_child];
        double left_enter = 0.0, left_exit = 0.0;
        double right_enter = 0.0, right_exit = 0.0;
        int hit_left = hier_fbbox_ray_enter_exit(&left->bbox, ox, oy, oz,
                                                 inv_dx, inv_dy, inv_dz,
                                                 t_min, t_max,
                                                 &left_enter, &left_exit);
        int hit_right = hier_fbbox_ray_enter_exit(&right->bbox, ox, oy, oz,
                                                  inv_dx, inv_dy, inv_dz,
                                                  t_min, t_max,
                                                  &right_enter, &right_exit);
        if (hit_left && hit_right && right_enter < left_enter) {
            if (query_universe_ray_blas_node(sys, blas, node->right_child,
                                             ox, oy, oz,
                                             inv_dx, inv_dy, inv_dz,
                                             t_min, t_max,
                                             out_hits, max_hits,
                                             hit_count) != 0) {
                return -1;
            }
            return query_universe_ray_blas_node(sys, blas, node->left_or_first,
                                                ox, oy, oz,
                                                inv_dx, inv_dy, inv_dz,
                                                t_min, t_max,
                                                out_hits, max_hits,
                                                hit_count);
        }
        if (hit_left) {
            if (query_universe_ray_blas_node(sys, blas, node->left_or_first,
                                             ox, oy, oz,
                                             inv_dx, inv_dy, inv_dz,
                                             t_min, t_max,
                                             out_hits, max_hits,
                                             hit_count) != 0) {
                return -1;
            }
        }
        if (hit_right) {
            return query_universe_ray_blas_node(sys, blas, node->right_child,
                                                ox, oy, oz,
                                                inv_dx, inv_dy, inv_dz,
                                                t_min, t_max,
                                                out_hits, max_hits,
                                                hit_count);
        }
        return 0;
    }

    for (uint16_t i = 0; i < node->count; i++) {
        uint32_t cell_pos = blas->indices[node->left_or_first + i];
        if (cell_pos >= blas->cell_count) return -1;
        const hier_blas_cell_t* blas_cell = &blas->cells[cell_pos];
        double t_enter, t_exit;
        if (!hier_fbbox_ray_enter_exit(&blas_cell->bbox, ox, oy, oz,
                                       inv_dx, inv_dy, inv_dz,
                                       t_min, t_max,
                                       &t_enter, &t_exit)) {
            continue;
        }
        insert_ray_candidate_sorted(sys, blas_cell->cell_index,
                                    t_enter, t_exit,
                                    out_hits, max_hits, hit_count);
    }

    return 0;
}

int alea_hier_spatial_query_universe_ray(alea_system_t* sys,
                                         int universe_id,
                                         double ox,
                                         double oy,
                                         double oz,
                                         double dx,
                                         double dy,
                                         double dz,
                                         double inv_dx,
                                         double inv_dy,
                                         double inv_dz,
                                         double t_min,
                                         double t_max,
                                         alea_hier_ray_candidate_t* out_hits,
                                         size_t max_hits) {
    (void)dx;
    (void)dy;
    (void)dz;
    if (!sys || !out_hits || max_hits == 0) return -1;
    if (!sys->hier_spatial_index || !sys->hier_spatial_index->built) {
        if (alea_hier_spatial_index_build(sys) != 0) return -1;
    }

    alea_hier_spatial_index_t* idx = sys->hier_spatial_index;
    size_t hit_count = 0;
    const hier_universe_blas_t* blas = find_blas(idx, universe_id);
    if (blas && blas->built && blas->node_count > 0) {
        int rc = query_universe_ray_blas_node(sys, blas, 0,
                                              ox, oy, oz,
                                              inv_dx, inv_dy, inv_dz,
                                              t_min, t_max,
                                              out_hits, max_hits,
                                              &hit_count);
        if (rc != 0) return -1;
    } else {
        const alea_universe_t* univ = alea_get_universe(sys, universe_id);
        if (!univ) return 0;
        for (size_t i = 0; i < univ->cell_indices.count; i++) {
            uint32_t cell_index = (uint32_t)univ->cell_indices.data[i];
            const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
            alea_bbox_t bbox = (cell->lat_type != 0 && cell->lat_fill)
                ? lattice_container_bbox(cell)
                : local_cell_bbox(sys, cell_index);
            hier_fbbox_t fbbox = hier_fbbox_from_double(&bbox);
            double t_enter, t_exit;
            if (!hier_fbbox_ray_enter_exit(&fbbox, ox, oy, oz,
                                           inv_dx, inv_dy, inv_dz,
                                           t_min, t_max,
                                           &t_enter, &t_exit)) {
                continue;
            }
            insert_ray_candidate_sorted(sys, cell_index, t_enter, t_exit,
                                        out_hits, max_hits, &hit_count);
        }
    }

    return (int)hit_count;
}

static int insert_placement_ray_candidate_sorted(
    const alea_hier_spatial_index_t* idx,
    uint32_t placement_index,
    double t_enter,
    double t_exit,
    alea_hier_placement_ray_candidate_t* out_hits,
    size_t max_hits,
    size_t* hit_count) {
    if (max_hits == 0) return 0;
    if ((size_t)placement_index >= idx->placement_count) return -1;

    size_t pos = *hit_count;
    if (pos >= max_hits) {
        if (t_enter > out_hits[max_hits - 1].t_enter ||
            (t_enter == out_hits[max_hits - 1].t_enter &&
             t_exit >= out_hits[max_hits - 1].t_exit)) {
            return 0;
        }
        pos = max_hits - 1;
    } else {
        (*hit_count)++;
    }

    while (pos > 0) {
        const alea_hier_placement_ray_candidate_t* prev = &out_hits[pos - 1];
        if (prev->t_enter < t_enter ||
            (prev->t_enter == t_enter && prev->t_exit <= t_exit)) {
            break;
        }
        out_hits[pos] = out_hits[pos - 1];
        pos--;
    }

    const hier_placement_t* placement = &idx->placements[placement_index];
    const alea_matrix_t* transform = placement_transform(idx, placement);
    if (!transform) return -1;

    out_hits[pos].placement_index = placement_index;
    out_hits[pos].parent_cell_index = placement->parent_cell_index;
    out_hits[pos].universe_id = placement->universe_id;
    out_hits[pos].depth = placement->depth;
    out_hits[pos].flags = placement->flags;
    out_hits[pos].transform = *transform;
    out_hits[pos].t_enter = t_enter;
    out_hits[pos].t_exit = t_exit;
    return 0;
}

static int query_tlas_ray_node(alea_hier_spatial_index_t* idx,
                               uint32_t node_index,
                               double ox,
                               double oy,
                               double oz,
                               double inv_dx,
                               double inv_dy,
                               double inv_dz,
                               double t_min,
                               double t_max,
                               alea_hier_placement_ray_candidate_t* out_hits,
                               size_t max_hits,
                               size_t* hit_count) {
    const hier_universe_blas_t* tlas = &idx->tlas;
    if (node_index >= tlas->node_count) return -1;

    const hier_bvh_node_t* node = &tlas->nodes[node_index];
    double node_enter, node_exit;
    if (!hier_fbbox_ray_enter_exit(&node->bbox, ox, oy, oz,
                                   inv_dx, inv_dy, inv_dz,
                                   t_min, t_max,
                                   &node_enter, &node_exit)) {
        return 0;
    }

    if (node->count == 0) {
        if (node->left_or_first >= tlas->node_count ||
            node->right_child >= tlas->node_count) {
            return -1;
        }
        const hier_bvh_node_t* left = &tlas->nodes[node->left_or_first];
        const hier_bvh_node_t* right = &tlas->nodes[node->right_child];
        double left_enter = 0.0, left_exit = 0.0;
        double right_enter = 0.0, right_exit = 0.0;
        int hit_left = hier_fbbox_ray_enter_exit(&left->bbox, ox, oy, oz,
                                                 inv_dx, inv_dy, inv_dz,
                                                 t_min, t_max,
                                                 &left_enter, &left_exit);
        int hit_right = hier_fbbox_ray_enter_exit(&right->bbox, ox, oy, oz,
                                                  inv_dx, inv_dy, inv_dz,
                                                  t_min, t_max,
                                                  &right_enter, &right_exit);
        if (hit_left && hit_right && right_enter < left_enter) {
            if (query_tlas_ray_node(idx, node->right_child,
                                    ox, oy, oz,
                                    inv_dx, inv_dy, inv_dz,
                                    t_min, t_max,
                                    out_hits, max_hits,
                                    hit_count) != 0) {
                return -1;
            }
            return query_tlas_ray_node(idx, node->left_or_first,
                                       ox, oy, oz,
                                       inv_dx, inv_dy, inv_dz,
                                       t_min, t_max,
                                       out_hits, max_hits,
                                       hit_count);
        }
        if (hit_left) {
            if (query_tlas_ray_node(idx, node->left_or_first,
                                    ox, oy, oz,
                                    inv_dx, inv_dy, inv_dz,
                                    t_min, t_max,
                                    out_hits, max_hits,
                                    hit_count) != 0) {
                return -1;
            }
        }
        if (hit_right) {
            return query_tlas_ray_node(idx, node->right_child,
                                       ox, oy, oz,
                                       inv_dx, inv_dy, inv_dz,
                                       t_min, t_max,
                                       out_hits, max_hits,
                                       hit_count);
        }
        return 0;
    }

    for (uint16_t i = 0; i < node->count; i++) {
        uint32_t pos = tlas->indices[node->left_or_first + i];
        if (pos >= tlas->cell_count) return -1;
        const hier_blas_cell_t* tlas_cell = &tlas->cells[pos];
        double t_enter, t_exit;
        if (!hier_fbbox_ray_enter_exit(&tlas_cell->bbox, ox, oy, oz,
                                       inv_dx, inv_dy, inv_dz,
                                       t_min, t_max,
                                       &t_enter, &t_exit)) {
            continue;
        }
        if (insert_placement_ray_candidate_sorted(
                idx, tlas_cell->cell_index, t_enter, t_exit,
                out_hits, max_hits, hit_count) != 0) {
            return -1;
        }
    }

    return 0;
}

int alea_hier_spatial_query_placements_ray(
    alea_system_t* sys,
    double ox,
    double oy,
    double oz,
    double dx,
    double dy,
    double dz,
    double inv_dx,
    double inv_dy,
    double inv_dz,
    double t_min,
    double t_max,
    alea_hier_placement_ray_candidate_t* out_hits,
    size_t max_hits) {
    (void)dx;
    (void)dy;
    (void)dz;
    if (!sys || !out_hits || max_hits == 0) return -1;
    if (!sys->hier_spatial_index || !sys->hier_spatial_index->built) {
        if (alea_hier_spatial_index_build(sys) != 0) return -1;
    }

    alea_hier_spatial_index_t* idx = sys->hier_spatial_index;
    if (!idx->tlas_built || idx->tlas.node_count == 0) return 0;

    size_t hit_count = 0;
    int rc = query_tlas_ray_node(idx, 0,
                                 ox, oy, oz,
                                 inv_dx, inv_dy, inv_dz,
                                 t_min, t_max,
                                 out_hits, max_hits,
                                 &hit_count);
    if (rc != 0) return -1;
    return (int)hit_count;
}

static int visit_lattice_placements_ray_node(
    alea_hier_spatial_index_t* idx,
    uint32_t node_index,
    double ox, double oy, double oz,
    double inv_dx, double inv_dy, double inv_dz,
    double t_min, double t_max,
    alea_hier_lattice_placement_ray_visitor_t visitor,
    void* userdata) {
    const hier_universe_blas_t* tlas = &idx->tlas;
    if (node_index >= tlas->node_count) return -1;
    const hier_bvh_node_t* node = &tlas->nodes[node_index];
    double node_enter, node_exit;
    if (!hier_fbbox_ray_enter_exit(&node->bbox, ox, oy, oz,
                                   inv_dx, inv_dy, inv_dz,
                                   t_min, t_max,
                                   &node_enter, &node_exit)) {
        return 0;
    }

    if (node->count == 0) {
        int rc = visit_lattice_placements_ray_node(
            idx, node->left_or_first, ox, oy, oz,
            inv_dx, inv_dy, inv_dz, t_min, t_max, visitor, userdata);
        if (rc != 0) return rc;
        return visit_lattice_placements_ray_node(
            idx, node->right_child, ox, oy, oz,
            inv_dx, inv_dy, inv_dz, t_min, t_max, visitor, userdata);
    }

    for (uint16_t i = 0; i < node->count; i++) {
        uint32_t pos = tlas->indices[node->left_or_first + i];
        if (pos >= tlas->cell_count) return -1;
        const hier_blas_cell_t* tlas_cell = &tlas->cells[pos];
        double t_enter, t_exit;
        if (!hier_fbbox_ray_enter_exit(&tlas_cell->bbox, ox, oy, oz,
                                       inv_dx, inv_dy, inv_dz,
                                       t_min, t_max, &t_enter, &t_exit)) {
            continue;
        }
        uint32_t placement_index = tlas_cell->cell_index;
        if ((size_t)placement_index >= idx->placement_count) return -1;
        const hier_placement_t* placement = &idx->placements[placement_index];
        if (!(placement->flags & HIER_PLACEMENT_LATTICE)) continue;
        const alea_matrix_t* transform = placement_transform(idx, placement);
        if (!transform) return -1;
        int rc = visitor(placement_index, placement->parent_cell_index,
                         transform, t_enter, t_exit, userdata);
        if (rc != 0) return rc;
    }
    return 0;
}

int alea_hier_spatial_visit_lattice_placements_ray(
    alea_system_t* sys,
    double ox, double oy, double oz,
    double inv_dx, double inv_dy, double inv_dz,
    double t_min, double t_max,
    alea_hier_lattice_placement_ray_visitor_t visitor,
    void* userdata) {
    if (!sys || !visitor) return -1;
    if (!sys->hier_spatial_index || !sys->hier_spatial_index->built) {
        if (alea_hier_spatial_index_build(sys) != 0) return -1;
    }
    alea_hier_spatial_index_t* idx = sys->hier_spatial_index;
    if (!idx->tlas_built || idx->tlas.node_count == 0) return 0;
    return visit_lattice_placements_ray_node(
        idx, 0, ox, oy, oz, inv_dx, inv_dy, inv_dz,
        t_min, t_max, visitor, userdata);
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
    /* Pure translation T(x,y,z): inverse is T(-x,-y,-z). Fill both directly
     * instead of going through the general 4x4 invert path. */
    out->m[0] = 1; out->m[1] = 0; out->m[2] = 0; out->m[3] = x;
    out->m[4] = 0; out->m[5] = 1; out->m[6] = 0; out->m[7] = y;
    out->m[8] = 0; out->m[9] = 0; out->m[10] = 1; out->m[11] = z;
    out->inv[0] = 1; out->inv[1] = 0; out->inv[2] = 0; out->inv[3] = -x;
    out->inv[4] = 0; out->inv[5] = 1; out->inv[6] = 0; out->inv[7] = -y;
    out->inv[8] = 0; out->inv[9] = 0; out->inv[10] = 1; out->inv[11] = -z;
    out->has_inverse = true;
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
    if (!hier_fbbox_intersects_dbbox(&node->bbox, local_query)) return 0;

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
        if (!hier_fbbox_intersects_dbbox(&blas_cell->bbox, local_query)) continue;

        const alea_cell_entry_t* cell = &sys->cells.data[blas_cell->cell_index];
        if (cell->fill_universe > 0 || (cell->lat_type != 0 && cell->lat_fill)) {
            continue;
        }

        alea_bbox_t cell_bbox_d = hier_fbbox_to_double(&blas_cell->bbox);
        alea_bbox_t world_bbox = alea_bbox_transform(&cell_bbox_d, transform);
        if (!bbox_intersects_local(&world_bbox, world_query)) continue;

        uint32_t synthetic_index =
            (uint32_t)((placement - idx->placements) + cell_pos);
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

        alea_bbox_t world_bbox = alea_bbox_transform(&local_bbox, transform);
        if (!bbox_intersects_local(&world_bbox, world_query)) continue;

        uint32_t synthetic_index =
            (uint32_t)((placement - idx->placements) + i);
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
        j0 = (cell->lat_fill_repeating || nj > 1) ? (int)floor((local_query->min_y - cell->lat_lower_left[1]) / cell->lat_pitch[1]) : 0;
        j1 = (cell->lat_fill_repeating || nj > 1) ? (int)floor((local_query->max_y - cell->lat_lower_left[1]) / cell->lat_pitch[1]) : 0;
        k0 = (cell->lat_fill_repeating || nk > 1) ? (int)floor((local_query->min_z - cell->lat_lower_left[2]) / cell->lat_pitch[2]) : 0;
        k1 = (cell->lat_fill_repeating || nk > 1) ? (int)floor((local_query->max_z - cell->lat_lower_left[2]) / cell->lat_pitch[2]) : 0;

        if (!cell->lat_fill_repeating) {
            if (i0 < 0) i0 = 0;
            if (j0 < 0) j0 = 0;
            if (k0 < 0) k0 = 0;
            if (i1 >= ni) i1 = ni - 1;
            if (j1 >= nj) j1 = nj - 1;
            if (k1 >= nk) k1 = nk - 1;
        }
    }

    if (i0 > i1 || j0 > j1 || k0 > k1) return 0;

    for (int i = i0; i <= i1; i++) {
        for (int j = j0; j <= j1; j++) {
            for (int k = k0; k <= k1; k++) {
                size_t fill_index = cell->lat_fill_repeating
                    ? 0 : (size_t)(i * nj * nk + j * nk + k);
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
                    oz = cell->lat_fill_zero_element_coords
                       ? (k + cell->lat_fill_dims[4]) * cell->lat_pitch[2]
                       : ((nk == 1) ? 0.0 : cell->lat_lower_left[2] + (k + 0.5) * cell->lat_pitch[2]);
                } else {
                    if (cell->lat_fill_zero_element_coords) {
                        ox = (i + cell->lat_fill_dims[0]) * cell->lat_pitch[0];
                        oy = (j + cell->lat_fill_dims[2]) * cell->lat_pitch[1];
                        oz = (k + cell->lat_fill_dims[4]) * cell->lat_pitch[2];
                    } else {
                        ox = cell->lat_lower_left[0] + (i + 0.5) * cell->lat_pitch[0];
                        oy = cell->lat_lower_left[1] + (j + 0.5) * cell->lat_pitch[1];
                        oz = cell->lat_lower_left[2] + (k + 0.5) * cell->lat_pitch[2];
                    }
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
        const alea_matrix_t* fill_matrix =
            cell_fill_matrix_cached(idx, cell_index);
        if (!fill_matrix) return -1;
        alea_matrix_t child_transform;
        alea_matrix_multiply(&child_transform, transform, fill_matrix);
        if (!alea_matrix_invert(&child_transform)) return -1;

        alea_bbox_t child_query = bbox_transform_inverse_conservative(world_query,
                                                                      &child_transform);
        if (!alea_bbox_is_valid(&child_query)) return 0;

        return query_region_universe_direct(sys, idx, cell->fill_universe,
                                            &child_transform, depth + 1,
                                            &child_query, world_query,
                                            out_hits, max_hits, hit_count);
    }

    alea_bbox_t world_bbox = alea_bbox_transform(cell_bbox, transform);
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
    if (!hier_fbbox_intersects_dbbox(&node->bbox, local_query)) return 0;

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
        alea_bbox_t cell_bbox_d = hier_fbbox_to_double(&blas_cell->bbox);
        if (query_region_cell_direct(sys, idx, blas_cell->cell_index, cell_pos,
                                     universe_id, transform, depth,
                                     &cell_bbox_d, local_query, world_query,
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

typedef struct {
    uint32_t cell_indices[ALEA_HIER_SPATIAL_HIT_CHAIN_MAX];
    alea_matrix_t transforms[ALEA_HIER_SPATIAL_HIT_CHAIN_MAX];
    uint8_t count;
    uint8_t truncated;
} hier_region_chain_t;

static void hier_region_chain_push(hier_region_chain_t* chain,
                                   uint32_t cell_index,
                                   const alea_matrix_t* transform) {
    if (chain->count >= ALEA_HIER_SPATIAL_HIT_CHAIN_MAX) {
        chain->truncated = 1;
        return;
    }
    chain->cell_indices[chain->count] = cell_index;
    chain->transforms[chain->count] = *transform;
    chain->count++;
}

static int append_region_chain_hit(alea_system_t* sys,
                                   const alea_matrix_t* transform,
                                   int universe_id,
                                   int depth,
                                   uint32_t cell_index,
                                   uint32_t synthetic_index,
                                   const hier_region_chain_t* chain,
                                   alea_hier_spatial_chain_hit_t* out_hits,
                                   size_t max_hits,
                                   size_t* hit_count) {
    if (*hit_count >= max_hits) return 0;
    if (!transform || !chain) return -1;

    const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
    alea_hier_spatial_chain_hit_t* out = &out_hits[*hit_count];
    out->hit.instance_index = synthetic_index;
    out->hit.cell_index = cell_index;
    out->hit.cell_id = cell->mc_cell_id;
    out->hit.material_id = cell->material_id;
    out->hit.universe_id = universe_id;
    out->hit.depth = depth;
    out->hit.is_terminal = true;
    out->hit.transform = *transform;
    out->ancestor_count = chain->count;
    out->chain_truncated = chain->truncated;
    for (uint8_t i = 0; i < chain->count; i++) {
        out->ancestor_cell_indices[i] = chain->cell_indices[i];
        out->ancestor_transforms[i] = chain->transforms[i];
    }
    (*hit_count)++;
    return 0;
}

static int query_region_universe_chain(alea_system_t* sys,
                                       alea_hier_spatial_index_t* idx,
                                       int universe_id,
                                       const alea_matrix_t* transform,
                                       int depth,
                                       const alea_bbox_t* local_query,
                                       const alea_bbox_t* world_query,
                                       const hier_region_chain_t* chain,
                                       alea_hier_spatial_chain_hit_t* out_hits,
                                       size_t max_hits,
                                       size_t* hit_count);

static int query_lattice_cell_chain(alea_system_t* sys,
                                    alea_hier_spatial_index_t* idx,
                                    const alea_cell_entry_t* cell,
                                    const alea_matrix_t* parent_transform,
                                    int depth,
                                    const alea_bbox_t* local_query,
                                    const alea_bbox_t* world_query,
                                    const hier_region_chain_t* chain,
                                    alea_hier_spatial_chain_hit_t* out_hits,
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
        j0 = (cell->lat_fill_repeating || nj > 1) ? (int)floor((local_query->min_y - cell->lat_lower_left[1]) / cell->lat_pitch[1]) : 0;
        j1 = (cell->lat_fill_repeating || nj > 1) ? (int)floor((local_query->max_y - cell->lat_lower_left[1]) / cell->lat_pitch[1]) : 0;
        k0 = (cell->lat_fill_repeating || nk > 1) ? (int)floor((local_query->min_z - cell->lat_lower_left[2]) / cell->lat_pitch[2]) : 0;
        k1 = (cell->lat_fill_repeating || nk > 1) ? (int)floor((local_query->max_z - cell->lat_lower_left[2]) / cell->lat_pitch[2]) : 0;

        if (!cell->lat_fill_repeating) {
            if (i0 < 0) i0 = 0;
            if (j0 < 0) j0 = 0;
            if (k0 < 0) k0 = 0;
            if (i1 >= ni) i1 = ni - 1;
            if (j1 >= nj) j1 = nj - 1;
            if (k1 >= nk) k1 = nk - 1;
        }
    }

    if (i0 > i1 || j0 > j1 || k0 > k1) return 0;

    for (int i = i0; i <= i1; i++) {
        for (int j = j0; j <= j1; j++) {
            for (int k = k0; k <= k1; k++) {
                size_t fill_index = cell->lat_fill_repeating
                    ? 0 : (size_t)(i * nj * nk + j * nk + k);
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
                    oz = cell->lat_fill_zero_element_coords
                       ? (k + cell->lat_fill_dims[4]) * cell->lat_pitch[2]
                       : ((nk == 1) ? 0.0 : cell->lat_lower_left[2] + (k + 0.5) * cell->lat_pitch[2]);
                } else {
                    if (cell->lat_fill_zero_element_coords) {
                        ox = (i + cell->lat_fill_dims[0]) * cell->lat_pitch[0];
                        oy = (j + cell->lat_fill_dims[2]) * cell->lat_pitch[1];
                        oz = (k + cell->lat_fill_dims[4]) * cell->lat_pitch[2];
                    } else {
                        ox = cell->lat_lower_left[0] + (i + 0.5) * cell->lat_pitch[0];
                        oy = cell->lat_lower_left[1] + (j + 0.5) * cell->lat_pitch[1];
                        oz = cell->lat_lower_left[2] + (k + 0.5) * cell->lat_pitch[2];
                    }
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

                if (query_region_universe_chain(sys, idx, fill_universe,
                                                &element_transform, depth + 1,
                                                &element_query, world_query,
                                                chain, out_hits, max_hits,
                                                hit_count) != 0) {
                    return -1;
                }
            }
        }
    }

    return 0;
}

static int query_region_cell_chain(alea_system_t* sys,
                                   alea_hier_spatial_index_t* idx,
                                   uint32_t cell_index,
                                   size_t cell_pos,
                                   int universe_id,
                                   const alea_matrix_t* transform,
                                   int depth,
                                   const alea_bbox_t* cell_bbox,
                                   const alea_bbox_t* local_query,
                                   const alea_bbox_t* world_query,
                                   const hier_region_chain_t* chain,
                                   alea_hier_spatial_chain_hit_t* out_hits,
                                   size_t max_hits,
                                   size_t* hit_count) {
    if (!bbox_intersects_local(cell_bbox, local_query)) return 0;

    const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
    if (cell->lat_type != 0 && cell->lat_fill) {
        hier_region_chain_t child_chain = *chain;
        hier_region_chain_push(&child_chain, cell_index, transform);
        return query_lattice_cell_chain(sys, idx, cell, transform, depth,
                                        local_query, world_query,
                                        &child_chain, out_hits, max_hits,
                                        hit_count);
    }

    if (cell->fill_universe > 0) {
        const alea_matrix_t* fill_matrix =
            cell_fill_matrix_cached(idx, cell_index);
        if (!fill_matrix) return -1;
        alea_matrix_t child_transform;
        alea_matrix_multiply(&child_transform, transform, fill_matrix);
        if (!alea_matrix_invert(&child_transform)) return -1;

        alea_bbox_t child_query = bbox_transform_inverse_conservative(world_query,
                                                                      &child_transform);
        if (!alea_bbox_is_valid(&child_query)) return 0;

        hier_region_chain_t child_chain = *chain;
        hier_region_chain_push(&child_chain, cell_index, transform);
        return query_region_universe_chain(sys, idx, cell->fill_universe,
                                           &child_transform, depth + 1,
                                           &child_query, world_query,
                                           &child_chain, out_hits, max_hits,
                                           hit_count);
    }

    alea_bbox_t world_bbox = alea_bbox_transform(cell_bbox, transform);
    if (!bbox_intersects_local(&world_bbox, world_query)) return 0;

    uint32_t synthetic_index = (uint32_t)(cell_pos + ((size_t)depth << 24));
    return append_region_chain_hit(sys, transform, universe_id, depth,
                                   cell_index, synthetic_index, chain,
                                   out_hits, max_hits, hit_count);
}

static int query_region_blas_node_chain(alea_system_t* sys,
                                        alea_hier_spatial_index_t* idx,
                                        const hier_universe_blas_t* blas,
                                        uint32_t node_index,
                                        int universe_id,
                                        const alea_matrix_t* transform,
                                        int depth,
                                        const alea_bbox_t* local_query,
                                        const alea_bbox_t* world_query,
                                        const hier_region_chain_t* chain,
                                        alea_hier_spatial_chain_hit_t* out_hits,
                                        size_t max_hits,
                                        size_t* hit_count) {
    if (node_index >= blas->node_count) return -1;

    const hier_bvh_node_t* node = &blas->nodes[node_index];
    if (!hier_fbbox_intersects_dbbox(&node->bbox, local_query)) return 0;

    if (node->count == 0) {
        if (query_region_blas_node_chain(sys, idx, blas, node->left_or_first,
                                         universe_id, transform, depth,
                                         local_query, world_query, chain,
                                         out_hits, max_hits, hit_count) != 0) {
            return -1;
        }
        return query_region_blas_node_chain(sys, idx, blas, node->right_child,
                                            universe_id, transform, depth,
                                            local_query, world_query, chain,
                                            out_hits, max_hits, hit_count);
    }

    for (uint16_t i = 0; i < node->count; i++) {
        uint32_t cell_pos = blas->indices[node->left_or_first + i];
        if (cell_pos >= blas->cell_count) return -1;

        const hier_blas_cell_t* blas_cell = &blas->cells[cell_pos];
        alea_bbox_t cell_bbox_d = hier_fbbox_to_double(&blas_cell->bbox);
        if (query_region_cell_chain(sys, idx, blas_cell->cell_index, cell_pos,
                                    universe_id, transform, depth,
                                    &cell_bbox_d, local_query, world_query,
                                    chain, out_hits, max_hits,
                                    hit_count) != 0) {
            return -1;
        }
    }

    return 0;
}

static int query_region_universe_chain(alea_system_t* sys,
                                       alea_hier_spatial_index_t* idx,
                                       int universe_id,
                                       const alea_matrix_t* transform,
                                       int depth,
                                       const alea_bbox_t* local_query,
                                       const alea_bbox_t* world_query,
                                       const hier_region_chain_t* chain,
                                       alea_hier_spatial_chain_hit_t* out_hits,
                                       size_t max_hits,
                                       size_t* hit_count) {
    const hier_universe_blas_t* blas = find_blas(idx, universe_id);
    if (blas && blas->built && blas->node_count > 0) {
        return query_region_blas_node_chain(sys, idx, blas, 0, universe_id,
                                            transform, depth, local_query,
                                            world_query, chain,
                                            out_hits, max_hits, hit_count);
    }

    const alea_universe_t* univ = alea_get_universe(sys, universe_id);
    if (!univ) return 0;
    for (size_t i = 0; i < univ->cell_indices.count; i++) {
        uint32_t cell_index = (uint32_t)univ->cell_indices.data[i];
        const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
        alea_bbox_t bbox = (cell->lat_type != 0 && cell->lat_fill)
            ? lattice_container_bbox(cell)
            : local_cell_bbox(sys, cell_index);
        if (query_region_cell_chain(sys, idx, cell_index, i, universe_id,
                                    transform, depth, &bbox, local_query,
                                    world_query, chain, out_hits, max_hits,
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

static int process_region_placement(alea_system_t* sys,
                                    alea_hier_spatial_index_t* idx,
                                    const hier_placement_t* placement,
                                    const alea_bbox_t* query_bbox,
                                    alea_spatial_hit_t* out_hits,
                                    size_t max_hits,
                                    size_t* hit_count) {
    const alea_matrix_t* transform = placement_transform(idx, placement);
    if (!transform) return -1;
    alea_bbox_t local_query = bbox_transform_inverse_conservative(query_bbox, transform);
    if (!alea_bbox_is_valid(&local_query)) return 0;

    if (placement->flags & HIER_PLACEMENT_LATTICE) {
        if (placement->parent_cell_index >= sys->cells.count) return -1;
        const alea_cell_entry_t* cell = &sys->cells.data[placement->parent_cell_index];
        return query_lattice_cell_direct(sys, idx, cell, transform,
                                         placement->depth, &local_query,
                                         query_bbox, out_hits, max_hits,
                                         hit_count);
    }

    const hier_universe_blas_t* blas = find_blas(idx, placement->universe_id);
    if (blas && blas->built) {
        return query_region_blas(sys, idx, placement, blas, &local_query,
                                 query_bbox, out_hits, max_hits, hit_count);
    }
    const alea_universe_t* univ = alea_get_universe(sys, placement->universe_id);
    if (!univ) return 0;
    return query_region_linear(sys, idx, placement, univ, &local_query,
                               query_bbox, out_hits, max_hits, hit_count);
}

static int query_tlas_region_node(alea_system_t* sys,
                                  alea_hier_spatial_index_t* idx,
                                  uint32_t node_idx,
                                  const alea_bbox_t* query_bbox,
                                  alea_spatial_hit_t* out_hits,
                                  size_t max_hits,
                                  size_t* hit_count) {
    if (*hit_count >= max_hits) return 0;
    if (node_idx >= idx->tlas.node_count) return -1;
    const hier_bvh_node_t* node = &idx->tlas.nodes[node_idx];
    if (!hier_fbbox_intersects_dbbox(&node->bbox, query_bbox)) return 0;

    if (node->count > 0) {
        for (uint16_t i = 0; i < node->count; i++) {
            if (*hit_count >= max_hits) return 0;
            uint32_t pos = idx->tlas.indices[node->left_or_first + i];
            if (pos >= idx->placement_count) return -1;
            const hier_placement_t* placement = &idx->placements[pos];
            /* Final per-placement bbox check: the leaf bbox unions multiple
             * placements, so an individual entry may not actually intersect. */
            if (alea_bbox_is_valid(&placement->world_bbox) &&
                !bbox_intersects_local(&placement->world_bbox, query_bbox)) {
                continue;
            }
            if (process_region_placement(sys, idx, placement, query_bbox,
                                         out_hits, max_hits, hit_count) != 0) {
                return -1;
            }
        }
        return 0;
    }

    if (query_tlas_region_node(sys, idx, node->left_or_first, query_bbox,
                               out_hits, max_hits, hit_count) != 0) return -1;
    return query_tlas_region_node(sys, idx, node->right_child, query_bbox,
                                  out_hits, max_hits, hit_count);
}

bool alea_hier_spatial_index_needs_rebuild(const alea_system_t* sys) {
    if (!sys || !sys->hier_spatial_index) return true;
    return !sys->hier_spatial_index->built;
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

    if (idx->tlas_built && idx->tlas.node_count > 0) {
        if (query_tlas_region_node(sys, idx, 0, query_bbox,
                                   out_hits, max_hits, &hit_count) != 0) {
            return -1;
        }
    } else {
        /* Fallback: no TLAS (e.g. zero placements). Walk the flat list. */
        for (size_t i = 0; i < idx->placement_count; i++) {
            if (hit_count >= max_hits) break;
            const hier_placement_t* placement = &idx->placements[i];
            if (alea_bbox_is_valid(&placement->world_bbox) &&
                !bbox_intersects_local(&placement->world_bbox, query_bbox)) {
                continue;
            }
            if (process_region_placement(sys, idx, placement, query_bbox,
                                         out_hits, max_hits, &hit_count) != 0) {
                return -1;
            }
        }
    }

    if (hit_count > 1) {
        qsort(out_hits, hit_count, sizeof(*out_hits), compare_spatial_hits_by_depth_cell);
    }

    return (int)hit_count;
}

int alea_hier_spatial_query_region_direct(alea_system_t* sys,
                                          const alea_bbox_t* query_bbox,
                                          alea_spatial_hit_t* out_hits,
                                          size_t max_hits) {
    if (!sys || !query_bbox || !out_hits) return -1;

    if (!sys->hier_spatial_index || !sys->hier_spatial_index->built) {
        if (alea_hier_spatial_index_build(sys) != 0) {
            return -1;
        }
    }

    alea_matrix_t identity;
    alea_matrix_identity(&identity);

    size_t hit_count = 0;
    int rc = query_region_universe_direct(sys, sys->hier_spatial_index, 0,
                                          &identity, 0,
                                          query_bbox, query_bbox,
                                          out_hits, max_hits, &hit_count);
    if (rc != 0) return -1;

    if (hit_count > 1) {
        qsort(out_hits, hit_count, sizeof(*out_hits),
              compare_spatial_hits_by_depth_cell);
    }

    return (int)hit_count;
}

int alea_hier_spatial_query_region_chain(alea_system_t* sys,
                                         const alea_bbox_t* query_bbox,
                                         alea_hier_spatial_chain_hit_t* out_hits,
                                         size_t max_hits) {
    if (!sys || !query_bbox || !out_hits) return -1;

    if (!sys->hier_spatial_index || !sys->hier_spatial_index->built) {
        if (alea_hier_spatial_index_build(sys) != 0) {
            return -1;
        }
    }

    alea_matrix_t identity;
    alea_matrix_identity(&identity);

    hier_region_chain_t chain;
    memset(&chain, 0, sizeof(chain));

    size_t hit_count = 0;
    int rc = query_region_universe_chain(sys, sys->hier_spatial_index, 0,
                                         &identity, 0,
                                         query_bbox, query_bbox,
                                         &chain, out_hits, max_hits,
                                         &hit_count);
    if (rc != 0) return -1;

    if (hit_count > 1) {
        /* Keep the same stable order as normal region hits. */
        for (size_t i = 1; i < hit_count; i++) {
            alea_hier_spatial_chain_hit_t x = out_hits[i];
            size_t j = i;
            while (j > 0) {
                const alea_spatial_hit_t* a = &out_hits[j - 1].hit;
                const alea_spatial_hit_t* b = &x.hit;
                int cmp = 0;
                if (a->depth != b->depth) {
                    cmp = (a->depth > b->depth) - (a->depth < b->depth);
                } else if (a->universe_id != b->universe_id) {
                    cmp = (a->universe_id > b->universe_id) -
                          (a->universe_id < b->universe_id);
                } else {
                    cmp = (a->cell_index > b->cell_index) -
                          (a->cell_index < b->cell_index);
                }
                if (cmp <= 0) break;
                out_hits[j] = out_hits[j - 1];
                j--;
            }
            out_hits[j] = x;
        }
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
