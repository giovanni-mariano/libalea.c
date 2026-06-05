// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file slice_api.c
 * @brief Public API wrappers for slice/render module
 */

#include "alea.h"
#include "slice/curve_intersect.h"  /* Must come before alea_slice.h for extended enum */
#include "alea_slice.h"
#include "core/alea_system.h"
#include "core/alea_spatial.h"
#include "core/alea_spatial_hier.h"
#include "core/alea_eval.h"
#include <stdio.h>
#include "core/alea_universe.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "util/alea_atomic.h"
#include "util/math.h"
#include "util/alea_log.h"

/* Grid query diagnostic counters — printed when ALEA_GRID_STATS=1.
 * Atomics are cheap vs CSG eval cost; no compile flag needed. */
static _Atomic size_t g_px_total         = 0;
static _Atomic size_t g_px_hint_full     = 0; /* all depths resolved by adjacency */
static _Atomic size_t g_px_hint_partial  = 0; /* adjacency broke mid-hierarchy */
static _Atomic size_t g_px_full_fallback = 0; /* no hint → full alea_find_all_cells_at_point */
static _Atomic size_t g_univ_hint_cell   = 0; /* exact hint cell hit */
static _Atomic size_t g_univ_neighbor    = 0; /* neighbor walk hit */
static _Atomic size_t g_univ_blas        = 0; /* BLAS query used */
static _Atomic size_t g_univ_linear      = 0; /* linear scan fallback */
static _Atomic size_t g_csg_prim_evals   = 0; /* total primitive evaluations */
static _Atomic size_t g_csg_bool_ops     = 0; /* total boolean ops in CSG tree */

static alea_tile_coverage_stats_t g_tile_coverage_stats;
static alea_point_coverage_stats_t g_point_coverage_stats;

void alea_tile_coverage_stats_reset(void) {
    memset(&g_tile_coverage_stats, 0, sizeof(g_tile_coverage_stats));
}

alea_tile_coverage_stats_t alea_tile_coverage_stats_get(void) {
    return g_tile_coverage_stats;
}

void alea_point_coverage_stats_reset(void) {
    memset(&g_point_coverage_stats, 0, sizeof(g_point_coverage_stats));
}

alea_point_coverage_stats_t alea_point_coverage_stats_get(void) {
    return g_point_coverage_stats;
}

/* ============================================================================
 * SLICE CURVES API
 * ============================================================================ */

struct alea_slice_curves {
    alea_curve_collection_t internal;
    double view_u_min, view_u_max;
    double view_v_min, view_v_max;
};

/* ============================================================================
 * SLICE VIEW SETUP
 * ============================================================================ */

static void init_plane_from_axis(alea_slice_plane_t* plane, int axis, double value) {
    alea_slice_plane_init_axis(plane, axis, value);
}

static void init_plane_from_vectors(alea_slice_plane_t* plane,
                                     double ox, double oy, double oz,
                                     double nx, double ny, double nz,
                                     double ux, double uy, double uz) {
    alea_slice_plane_init(plane, ox, oy, oz, nx, ny, nz, ux, uy, uz);
}

void alea_slice_view_axis(alea_slice_view_t* view,
                               int axis, double value,
                               double u_min, double u_max,
                               double v_min, double v_max) {
    if (!view) return;
    init_plane_from_axis(&view->plane, axis, value);
    view->u_min = u_min;
    view->u_max = u_max;
    view->v_min = v_min;
    view->v_max = v_max;
}

void alea_slice_view_init(alea_slice_view_t* view,
                               double ox, double oy, double oz,
                               double nx, double ny, double nz,
                               double ux, double uy, double uz,
                               double u_min, double u_max,
                               double v_min, double v_max) {
    if (!view) return;
    init_plane_from_vectors(&view->plane, ox, oy, oz, nx, ny, nz, ux, uy, uz);
    view->u_min = u_min;
    view->u_max = u_max;
    view->v_min = v_min;
    view->v_max = v_max;
}

/* ============================================================================
 * ANALYTICAL CURVE API
 * ============================================================================ */

alea_slice_curves_t* alea_get_slice_curves(alea_system_t* sys,
                                                    const alea_slice_view_t* view) {
    if (!sys || !view) return NULL;

    alea_slice_curves_t* result = calloc(1, sizeof(alea_slice_curves_t));
    if (!result) return NULL;

    /* Mode-aware: builds flat spatial index or hier index depending on config. */
    if (alea_prepare_query_acceleration(sys) != 0) {
        free(result);
        return NULL;
    }

    int ret = alea_compute_slice_curves_spatial(sys, &view->plane,
                                                view->u_min, view->u_max,
                                                view->v_min, view->v_max,
                                                &result->internal);
    if (ret != 0) {
        free(result);
        return NULL;
    }

    result->view_u_min = view->u_min;
    result->view_u_max = view->u_max;
    result->view_v_min = view->v_min;
    result->view_v_max = view->v_max;
    return result;
}

size_t alea_slice_curves_count(const alea_slice_curves_t* curves) {
    return curves ? curves->internal.curves.count : 0;
}

int alea_slice_curves_get(const alea_slice_curves_t* curves, size_t index, alea_curve_t* out) {
    if (!curves || !out || index >= curves->internal.curves.count) return -1;
    const alea_curve_2d_t* src = &curves->internal.curves.data[index];
    memset(out, 0, sizeof(*out));

    out->surface_id = src->surface_id;
    out->primitive_id = src->primitive_id;
    out->t_min = src->bounds.t_min;
    out->t_max = src->bounds.t_max;

    switch (src->type) {
        case ALEA_CURVE_LINE:
        case ALEA_CURVE_LINE_SEGMENT:
        case ALEA_CURVE_RAY:
            out->type = (src->type == ALEA_CURVE_LINE) ? ALEA_CURVE_LINE : ALEA_CURVE_LINE_SEGMENT;
            out->data.line.point[0] = src->data.line.point[0];
            out->data.line.point[1] = src->data.line.point[1];
            out->data.line.direction[0] = src->data.line.direction[0];
            out->data.line.direction[1] = src->data.line.direction[1];
            break;

        case ALEA_CURVE_CIRCLE:
            out->type = ALEA_CURVE_CIRCLE;
            out->data.circle.center[0] = src->data.circle.center[0];
            out->data.circle.center[1] = src->data.circle.center[1];
            out->data.circle.radius = src->data.circle.radius;
            break;

        case ALEA_CURVE_ARC:
            out->type = ALEA_CURVE_ARC;
            out->data.circle.center[0] = src->data.circle.center[0];
            out->data.circle.center[1] = src->data.circle.center[1];
            out->data.circle.radius = src->data.circle.radius;
            out->t_min = src->bounds.theta_start;
            out->t_max = src->bounds.theta_end;
            break;

        case ALEA_CURVE_ELLIPSE:
        case ALEA_CURVE_ELLIPSE_ARC:
            out->type = (src->type == ALEA_CURVE_ELLIPSE) ? ALEA_CURVE_ELLIPSE : ALEA_CURVE_ELLIPSE_ARC;
            /* Canonical form already computed by finalize_conic_as_ellipse() —
             * just copy data.ellipse directly. */
            out->data.ellipse.center[0] = src->data.ellipse.center[0];
            out->data.ellipse.center[1] = src->data.ellipse.center[1];
            out->data.ellipse.semi_a = src->data.ellipse.semi_a;
            out->data.ellipse.semi_b = src->data.ellipse.semi_b;
            out->data.ellipse.angle = src->data.ellipse.angle;
            break;

        case ALEA_CURVE_POLYGON:
            out->type = ALEA_CURVE_POLYGON;
            out->data.polygon.count = src->data.polygon.vertex_count;
            out->data.polygon.closed = src->data.polygon.closed ? 1 : 0;
            for (int i = 0; i < src->data.polygon.vertex_count && i < 16; i++) {
                out->data.polygon.vertices[i][0] = src->data.polygon.vertices[i][0];
                out->data.polygon.vertices[i][1] = src->data.polygon.vertices[i][1];
            }
            break;

        case ALEA_CURVE_PARALLEL_LINES:
            out->type = ALEA_CURVE_PARALLEL_LINES;
            out->data.parallel_lines.point1[0] = src->data.parallel_lines.point1[0];
            out->data.parallel_lines.point1[1] = src->data.parallel_lines.point1[1];
            out->data.parallel_lines.point2[0] = src->data.parallel_lines.point2[0];
            out->data.parallel_lines.point2[1] = src->data.parallel_lines.point2[1];
            out->data.parallel_lines.direction[0] = src->data.parallel_lines.direction[0];
            out->data.parallel_lines.direction[1] = src->data.parallel_lines.direction[1];
            break;

        default:
            out->type = ALEA_CURVE_NONE;
            break;
    }

    return 0;
}

void alea_slice_curves_bounds(const alea_slice_curves_t* curves,
                              double* u_min, double* u_max,
                              double* v_min, double* v_max) {
    if (!curves) return;
    if (u_min) *u_min = curves->internal.u_min;
    if (u_max) *u_max = curves->internal.u_max;
    if (v_min) *v_min = curves->internal.v_min;
    if (v_max) *v_max = curves->internal.v_max;
}

void alea_slice_curves_free(alea_slice_curves_t* curves) {
    if (curves) {
        alea_curve_collection_free(&curves->internal);
        free(curves);
    }
}

/* ============================================================================
 * GRID-BASED CELL QUERIES
 * ============================================================================ */


/* ============================================================================
 * GRID-BASED CELL QUERIES WITH ERROR DETECTION
 * ============================================================================ */

/* Error codes */
#define GRID_ERR_NONE      0
#define GRID_ERR_OVERLAP   1
#define GRID_ERR_UNDEFINED 2

/* Forward declarations for functions used across sections */
static void get_clipped_param_range(const alea_curve_2d_t* curve,
                                    const alea_slice_view_t* view,
                                    double* t_lo, double* t_hi);
static int dedup_spatial_hits(alea_spatial_hit_t* hits, int hit_count);
typedef struct {
    uint8_t coverage;
    int primary_cell_id;
    int secondary_cell_id;
} point_coverage_t;
static int find_point_coverage_exact(alea_system_t* sys,
                                     double gx, double gy, double gz,
                                     int universe_depth,
                                     point_coverage_t* out);
static bool point_coverage_matches_pair(const point_coverage_t* pc,
                                        int cell_a,
                                        int cell_b);
static int filter_grid_overlap_ambiguities(
    alea_system_t* sys,
    const alea_slice_view_t* view,
    int nu,
    int nv,
    int universe_depth,
    const int* primary_cell_ids,
    int* secondary_cell_ids,
    uint8_t* coverage,
    uint8_t* errors,
    alea_boundary_filter_stats_t* out_stats);

void alea_slice_path_table_free(alea_slice_path_table_t* table) {
    if (!table) return;
    free(table->records);
    table->records = NULL;
    table->count = 0;
    table->capacity = 0;
}

static uint64_t path_hash_mix_u64(uint64_t h, uint64_t v) {
    h ^= v;
    h *= 1099511628211ULL;
    return h;
}

static uint64_t path_hash_matrix_inverse(const alea_matrix_t* m) {
    uint64_t h = 1469598103934665603ULL;
    for (int i = 0; i < 12; i++) {
        union { double d; uint64_t u; } cvt;
        cvt.d = m->inv[i];
        h = path_hash_mix_u64(h, cvt.u);
    }
    return h;
}

static bool path_matrix_inverse_equal_exact(const alea_matrix_t* a,
                                            const double b[12]) {
    for (int i = 0; i < 12; i++) {
        if (a->inv[i] != b[i]) return false;
    }
    return true;
}

static int slice_path_table_intern(alea_slice_path_table_t* table,
                                   int universe_id,
                                   int depth,
                                   const alea_matrix_t* world_to_local,
                                   uint32_t* out_id) {
    if (!table || !world_to_local || !out_id) return -1;
    if (!world_to_local->has_inverse) return -1;

    uint64_t h = 1469598103934665603ULL;
    h = path_hash_mix_u64(h, (uint64_t)(uint32_t)universe_id);
    h = path_hash_mix_u64(h, (uint64_t)(uint32_t)depth);
    h = path_hash_mix_u64(h, path_hash_matrix_inverse(world_to_local));

    for (size_t i = 0; i < table->count; i++) {
        const alea_slice_path_record_t* r = &table->records[i];
        if (r->universe_id == universe_id &&
            r->depth == depth &&
            r->chain_hash == h &&
            path_matrix_inverse_equal_exact(world_to_local,
                                            r->world_to_local)) {
            *out_id = (uint32_t)i;
            return 0;
        }
    }

    if (table->count == table->capacity) {
        size_t new_cap = table->capacity ? table->capacity * 2 : 64;
        alea_slice_path_record_t* grown =
            realloc(table->records, new_cap * sizeof(*grown));
        if (!grown) return -1;
        table->records = grown;
        table->capacity = new_cap;
    }

    if (table->count > UINT32_MAX) return -1;
    uint32_t id = (uint32_t)table->count;
    alea_slice_path_record_t* rec = &table->records[table->count++];
    rec->universe_id = universe_id;
    rec->depth = depth;
    rec->chain_hash = h;
    for (int i = 0; i < 12; i++)
        rec->world_to_local[i] = world_to_local->inv[i];
    *out_id = id;
    return 0;
}

/* ============================================================================
 * PER-UNIVERSE ADJACENCY HINTS
 * ============================================================================ */

/**
 * Maximum nesting depth for multi-level hints.
 * Typical reactor models rarely exceed 4-5 levels.
 */
#define MAX_HINT_DEPTH 8

/**
 * @brief Hint information for a single depth level
 *
 * Stores the cell index, universe, and accumulated transform
 * to allow adjacency walking within a nested universe.
 */
typedef struct {
    int cell_index;              /**< Cell index at this depth (-1 if invalid) */
    int universe_id;             /**< Universe ID at this depth */
    alea_matrix_t transform;      /**< Accumulated transform: world -> local */
} alea_level_hint_t;

/**
 * @brief Multi-level hint structure for per-universe adjacency
 *
 * Tracks the cell path through the hierarchy for coherence-based lookup.
 * When adjacent pixels share the same parent path, we can use adjacency
 * walking at each level instead of full hierarchy traversal.
 */
typedef struct {
    int valid_depth;                      /**< Deepest valid depth (-1 if none) */
    alea_level_hint_t levels[MAX_HINT_DEPTH];  /**< Hints at each depth */
} alea_multilevel_hint_t;

/**
 * @brief Initialize a multi-level hint to empty state
 */
static void multilevel_hint_init(alea_multilevel_hint_t* hint) {
    hint->valid_depth = -1;
    for (int d = 0; d < MAX_HINT_DEPTH; d++) {
        hint->levels[d].cell_index = -1;
        hint->levels[d].universe_id = -1;
    }
}

/**
 * @brief Build multi-level hints from a hierarchy traversal result
 *
 * @param sys CSG system
 * @param hits Array of cell hits from alea_find_all_cells_at_point()
 * @param num_hits Number of hits
 * @param hint Output multi-level hint structure
 */
static void multilevel_hint_from_hits(alea_system_t* sys,
                                       const alea_cell_hit_t* hits,
                                       int num_hits,
                                       alea_multilevel_hint_t* hint) {
    multilevel_hint_init(hint);
    if (num_hits <= 0 || num_hits > MAX_HINT_DEPTH) return;

    /* Build transform chain from root to each level */
    alea_matrix_t accumulated;
    alea_matrix_identity(&accumulated);

    for (int i = 0; i < num_hits && i < MAX_HINT_DEPTH; i++) {
        const alea_cell_hit_t* h = &hits[i];
        alea_level_hint_t* level = &hint->levels[i];

        level->cell_index = h->cell_index;
        level->universe_id = h->universe_id;
        level->transform = accumulated;  /* Copy current accumulated transform */

        /* If this cell has a fill, compose with its transform for next level */
        if (h->fill_universe > 0 && h->cell_index >= 0 &&
            (size_t)h->cell_index < alea_vec_count(&sys->cells)) {
            const alea_cell_entry_t* cell = &sys->cells.data[h->cell_index];
            if (cell->fill_transform > 0) {
                const alea_transform_t* tr = alea_get_transform(sys, cell->fill_transform);
                if (tr) {
                    alea_matrix_t fill_mat;
                    /* Use tr->cosines which has pre-computed direction cosines */
                    if (!alea_matrix_from_mcnp(&fill_mat, tr->cosines,
                                               tr->value_count, false)) {
                        return;
                    }
                    alea_matrix_t new_acc;
                    alea_matrix_multiply(&new_acc, &accumulated, &fill_mat);
                    accumulated = new_acc;
                }
            }
        }

        hint->valid_depth = i;
    }
}

static int slice_path_table_intern_hint(alea_slice_path_table_t* table,
                                        const alea_multilevel_hint_t* hint,
                                        int universe_depth,
                                        uint32_t* out_id) {
    if (!table || !hint || !out_id || hint->valid_depth < 0)
        return -1;

    int depth = universe_depth < 0 ? hint->valid_depth : universe_depth;
    if (depth < 0 || depth > hint->valid_depth || depth >= MAX_HINT_DEPTH)
        return -1;

    const alea_level_hint_t* level = &hint->levels[depth];
    if (level->cell_index < 0)
        return -1;

    alea_matrix_t transform = level->transform;
    if (!transform.has_inverse && !alea_matrix_invert(&transform))
        return -1;

    return slice_path_table_intern(table, level->universe_id, depth,
                                   &transform, out_id);
}

/**
 * @brief Find cell at a specific depth using adjacency walking within that universe
 *
 * Prerequisites:
 * - Parent path (depths 0..depth-1) is already validated
 * - Local coordinates at this depth are known
 *
 * @param sys CSG system
 * @param lx,ly,lz Point in local coordinates at this depth
 * @param universe_id Universe to search in
 * @param hint_cell_idx Previous cell index at this depth (for adjacency walking)
 * @return Cell index found, or -1 if not found
 */
static int find_cell_in_universe_with_hint(alea_system_t* sys,
                                            double lx, double ly, double lz,
                                            int universe_id,
                                            int hint_cell_idx) {
    /* Try hint cell first (adjacency cache).
     * Always verify with alea_contains_point — this is cheap for the common
     * case (point still inside the same cell, single CSG eval with early exit)
     * and immediately detects cell-boundary crossings without losing information. */
    if (hint_cell_idx >= 0 && (size_t)hint_cell_idx < alea_vec_count(&sys->cells)) {
        const alea_cell_entry_t* cell = &sys->cells.data[hint_cell_idx];

        if (cell->universe_id == universe_id &&
            cell->root_node_id != ALEA_NODE_ID_INVALID &&
            alea_contains_point(sys, cell->root_node_id, lx, ly, lz)) {
            atomic_fetch_add_explicit(&g_univ_hint_cell, 1, memory_order_relaxed);
            return hint_cell_idx;
        }

        /* Try neighbors of hint cell */
        if (sys->cell_adjacency_built && cell->neighbors) {
            for (size_t n = 0; n < cell->neighbor_count; n++) {
                int neighbor_idx = (int)cell->neighbors[n].neighbor_index;
                if (neighbor_idx < 0 || (size_t)neighbor_idx >= alea_vec_count(&sys->cells)) continue;

                const alea_cell_entry_t* neighbor = &sys->cells.data[neighbor_idx];

                /* Neighbor must be in same universe */
                if (neighbor->universe_id != universe_id) continue;
                if (neighbor->root_node_id == ALEA_NODE_ID_INVALID) continue;

                /* Quick bbox check */
                const alea_bbox_t* bbox = &sys->nodes.data[neighbor->root_node_id].bbox;
                if (lx < bbox->min_x || lx > bbox->max_x ||
                    ly < bbox->min_y || ly > bbox->max_y ||
                    lz < bbox->min_z || lz > bbox->max_z) {
                    continue;
                }

                if (alea_contains_point(sys, neighbor->root_node_id, lx, ly, lz)) {
                    atomic_fetch_add_explicit(&g_univ_neighbor, 1, memory_order_relaxed);
                    return neighbor_idx;
                }
            }
        }
    }

    /* Hier-spatial fast path: BLAS-pruned lookup before linear scan. On
     * models with very large universes the linear scan would do O(N)
     * alea_contains_point evaluations per pixel after every adjacency miss.
     * The BLAS reduces that to O(log N). */
    if (sys->hier_spatial_index) {
        int found = alea_hier_spatial_find_cell_in_universe(sys, universe_id,
                                                            lx, ly, lz);
        if (found >= 0) { atomic_fetch_add_explicit(&g_univ_blas, 1, memory_order_relaxed); return found; }
        if (found == -1) { atomic_fetch_add_explicit(&g_univ_blas, 1, memory_order_relaxed); return -1; }
        /* found == -2: hier index not usable for this universe, fall through. */
    }

    /* Fall back to linear search within universe */
    const alea_universe_t* univ = alea_get_universe(sys, universe_id);
    if (!univ) return -1;

    atomic_fetch_add_explicit(&g_univ_linear, 1, memory_order_relaxed);
    for (size_t i = 0; i < univ->cell_indices.count; i++) {
        size_t cell_idx = univ->cell_indices.data[i];
        const alea_cell_entry_t* cell = &sys->cells.data[cell_idx];

        if (cell->root_node_id == ALEA_NODE_ID_INVALID) continue;

        if (alea_contains_point(sys, cell->root_node_id, lx, ly, lz)) {
            return (int)cell_idx;
        }
    }

    return -1;  /* Not found in this universe */
}

/**
 * @brief Find cell with multi-level adjacency walking
 *
 * Uses hints at each depth level for coherence-based optimization:
 * 1. At each depth, if parent path matches previous pixel, use adjacency walking
 * 2. If adjacency walk succeeds, continue to next depth with updated transform
 * 3. If it fails at any level, fall back to full search from that level
 *
 * @param sys CSG system
 * @param gx,gy,gz Point in global coordinates
 * @param universe_depth Target depth (-1 for innermost)
 * @param prev_hint Previous pixel's multi-level hint (NULL for no hint)
 * @param out_hint Output hint for current pixel (can be NULL)
 * @param out_cell_id Output cell ID
 * @param out_material_id Output material ID (can be NULL)
 * @param out_error Output error code (can be NULL)
 */
static void find_cell_multilevel(alea_system_t* sys,
                                  double gx, double gy, double gz,
                                  int universe_depth,
                                  const alea_multilevel_hint_t* prev_hint,
                                  alea_multilevel_hint_t* out_hint,
                                  int* out_cell_id,
                                  int* out_material_id,
                                  uint8_t* out_error) {
    *out_cell_id = -1;
    if (out_material_id) *out_material_id = 0;
    if (out_error) *out_error = GRID_ERR_NONE;
    if (out_hint) multilevel_hint_init(out_hint);
    atomic_fetch_add_explicit(&g_px_total, 1, memory_order_relaxed);

    /* Use multi-level walking only when we have valid previous hints */
    bool use_multilevel = (prev_hint != NULL && prev_hint->valid_depth >= 0 &&
                           sys->cell_adjacency_built);

    if (use_multilevel) {
        /* Try to walk through hierarchy using hints at each level */
        alea_matrix_t accumulated;
        alea_matrix_identity(&accumulated);
        int current_universe = 0;  /* Start at root universe */
        double lx = gx, ly = gy, lz = gz;  /* Local coordinates (start = global) */
        int hits_found = 0;
        alea_cell_hit_t local_hits[MAX_HINT_DEPTH];

        for (int depth = 0; depth <= prev_hint->valid_depth && depth < MAX_HINT_DEPTH; depth++) {
            const alea_level_hint_t* level_hint = &prev_hint->levels[depth];

            /* Check if parent path still matches (for depth > 0) */
            bool parent_matches = true;
            if (depth > 0) {
                /* We've already validated depths 0..depth-1 by the loop */
                parent_matches = (level_hint->cell_index >= 0);
            }

            if (!parent_matches) break;

            /* Try adjacency walking at this depth */
            int cell_idx = find_cell_in_universe_with_hint(sys, lx, ly, lz,
                                                           current_universe,
                                                           level_hint->cell_index);

            if (cell_idx < 0) {
                /* Adjacency walk failed - need full search from here */
                break;
            }

            /* Record this hit */
            const alea_cell_entry_t* cell = &sys->cells.data[cell_idx];
            if (hits_found < MAX_HINT_DEPTH) {
                alea_cell_hit_t* hit = &local_hits[hits_found];
                hit->cell_id = cell->mc_cell_id;
                hit->cell_index = cell_idx;
                hit->material_id = cell->material_id;
                hit->universe_id = current_universe;
                hit->fill_universe = cell->fill_universe;
                hit->depth = depth;
                hit->local_x = lx;
                hit->local_y = ly;
                hit->local_z = lz;
                hits_found++;
            }

            /* Check for target depth */
            if (universe_depth >= 0 && depth >= universe_depth) {
                /* Reached target depth — full adjacency success */
                atomic_fetch_add_explicit(&g_px_hint_full, 1, memory_order_relaxed);
                *out_cell_id = cell->mc_cell_id;
                if (out_material_id) *out_material_id = cell->material_id;
                if (out_hint) {
                    multilevel_hint_from_hits(sys, local_hits, hits_found, out_hint);
                }
                return;
            }

            /* If cell has a fill, continue to next level */
            if (cell->fill_universe > 0) {
                /* Prefer the precomputed per-cell fill matrix from the hier
                 * index — it avoids `alea_matrix_from_mcnp` (transform-table
                 * lookup + invert) per pixel per descent. Falls back to the
                 * per-pixel build if the hier cache is unavailable. */
                const alea_matrix_t* cached_fill =
                    alea_hier_spatial_get_cell_fill_matrix(sys, (uint32_t)cell_idx);
                alea_matrix_t fill_transform_local;
                const alea_matrix_t* fill_transform = cached_fill;
                if (!fill_transform) {
                    if (cell->fill_transform > 0) {
                        const alea_transform_t* tr = alea_get_transform(sys, cell->fill_transform);
                        if (tr) {
                            if (!alea_matrix_from_mcnp(&fill_transform_local, tr->cosines,
                                                       tr->value_count, false)) {
                                break;
                            }
                        } else {
                            break;
                        }
                    } else {
                        alea_matrix_identity(&fill_transform_local);
                    }
                    fill_transform = &fill_transform_local;
                }

                /* Compose with accumulated transform */
                alea_matrix_t new_accumulated;
                alea_matrix_multiply(&new_accumulated, &accumulated, fill_transform);
                accumulated = new_accumulated;

                /* Transform point to fill universe coordinates */
                /* Use a copy for inversion to preserve accumulated for next iteration */
                alea_matrix_t inv_accumulated = accumulated;
                if (!alea_matrix_invert(&inv_accumulated)) break;
                lx = gx; ly = gy; lz = gz;
                alea_matrix_transform_point_inverse(&inv_accumulated, &lx, &ly, &lz);

                current_universe = cell->fill_universe;
            } else {
                /* Terminal cell (no fill) - this is the innermost */
                if (universe_depth < 0) {
                    /* Full adjacency success — reached terminal cell */
                    atomic_fetch_add_explicit(&g_px_hint_full, 1, memory_order_relaxed);
                    *out_cell_id = cell->mc_cell_id;
                    if (out_material_id) *out_material_id = cell->material_id;
                    if (out_hint) {
                        multilevel_hint_from_hits(sys, local_hits, hits_found, out_hint);
                    }
                    return;
                }
                break;  /* No more levels to explore */
            }
        }

        /* If we got some hits but didn't complete, we might have partial success */
        if (hits_found > 0 && universe_depth < 0) {
            /* Return deepest hit found */
            alea_cell_hit_t* last = &local_hits[hits_found - 1];
            if (last->fill_universe <= 0) {
                /* This is a terminal cell — partial adjacency success */
                atomic_fetch_add_explicit(&g_px_hint_partial, 1, memory_order_relaxed);
                *out_cell_id = last->cell_id;
                if (out_material_id) *out_material_id = last->material_id;
                if (out_hint) {
                    multilevel_hint_from_hits(sys, local_hits, hits_found, out_hint);
                }
                return;
            }
        }
    }

    /* Fall back to full hierarchy traversal */
    if (use_multilevel)
        atomic_fetch_add_explicit(&g_px_hint_partial, 1, memory_order_relaxed);
    else
        atomic_fetch_add_explicit(&g_px_full_fallback, 1, memory_order_relaxed);

    alea_cell_hit_t hits[32];
    int num_hits = alea_find_all_cells_at_point(sys, gx, gy, gz, hits, 32);

    if (num_hits <= 0) {
        if (out_error) *out_error = GRID_ERR_UNDEFINED;
        return;
    }

    /* Build output hint from hits */
    if (out_hint) {
        multilevel_hint_from_hits(sys, hits, num_hits, out_hint);
    }

    /* Determine target depth */
    int target_idx;
    int target_depth;

    if (universe_depth < 0) {
        target_idx = num_hits - 1;
        target_depth = hits[target_idx].depth;
    } else {
        target_idx = -1;
        target_depth = universe_depth;

        for (int i = 0; i < num_hits; i++) {
            if (hits[i].depth == universe_depth) {
                target_idx = i;
                break;
            }
        }

        if (target_idx < 0) {
            for (int i = num_hits - 1; i >= 0; i--) {
                if (hits[i].depth <= universe_depth) {
                    target_idx = i;
                    target_depth = hits[i].depth;
                    break;
                }
            }
        }

        if (target_idx < 0) {
            if (out_error) *out_error = GRID_ERR_UNDEFINED;
            return;
        }
    }

    /* Check for overlap at target depth */
    if (out_error) {
        int count_at_depth = 0;
        for (int i = 0; i < num_hits; i++) {
            if (hits[i].depth == target_depth) {
                count_at_depth++;
                if (count_at_depth > 1) {
                    *out_error = GRID_ERR_OVERLAP;
                    break;
                }
            }
        }
    }

    *out_cell_id = hits[target_idx].cell_id;
    if (out_material_id) *out_material_id = hits[target_idx].material_id;
}

int alea_find_cells_grid(alea_system_t* sys,
                              const alea_slice_view_t* view,
                              int nu, int nv,
                              int universe_depth,
                              int* out_cell_ids,
                              int* out_material_ids,
                              uint8_t* out_errors) {
    if (!sys || !view || !out_cell_ids || nu <= 0 || nv <= 0) {
        return -1;
    }

    const alea_slice_plane_t* plane = &view->plane;
    double u_min = view->u_min;
    double u_max = view->u_max;
    double v_min = view->v_min;
    double v_max = view->v_max;

    /* Build query acceleration eagerly before the parallel region so that:
     * (a) the internal #pragma omp parallel for gets full thread parallelism
     *     (no nesting penalty)
     * (b) OpenMP workers only read shared cache state */
    if (alea_system_prepare_query_caches(sys, ALEA_CACHE_RAYCAST) != 0)
        return -1;

    double du = (u_max - u_min) / nu;
    double dv = (v_max - v_min) / nv;

    const double* origin = plane->origin;
    const double* u_axis = plane->u_axis;
    const double* v_axis = plane->v_axis;

    /* Initialize error map if provided */
    if (out_errors) {
        memset(out_errors, 0, nu * nv);
    }

    /*
     * Per-universe adjacency optimization:
     * - Track multi-level hints (cell path at each depth)
     * - When parent path matches, use adjacency walking at each level
     * - Works for ALL depths, not just depth=0
     */

#ifdef _OPENMP
    bool stats_en = getenv("ALEA_GRID_STATS") != NULL;
    /* ALEA_GRID_VERIFY_INTERVAL=N: every N-th pixel also runs a full recursive
     * search to detect overlapping geometry. 0 or unset = no periodic scan
     * (relies on boundary-pixel second pass only). Useful for finding geometry
     * errors in interior regions not adjacent to cell boundaries. */
    const char* vi_env = getenv("ALEA_GRID_VERIFY_INTERVAL");
    int overlap_interval = (vi_env && atoi(vi_env) > 0) ? atoi(vi_env) : 0;
    #pragma omp parallel
    {
        if (stats_en) alea_perf_reset();
        #pragma omp for schedule(dynamic, 4)
        for (int j = 0; j < nv; j++) {
            double v = v_min + (j + 0.5) * dv;
            alea_multilevel_hint_t row_hint;
            multilevel_hint_init(&row_hint);

            for (int i = 0; i < nu; i++) {
                double u = u_min + (i + 0.5) * du;
                int idx = j * nu + i;

                double x = origin[0] + u * u_axis[0] + v * v_axis[0];
                double y = origin[1] + u * u_axis[1] + v * v_axis[1];
                double z = origin[2] + u * u_axis[2] + v * v_axis[2];

                int cell_id = -1, material_id = 0;
                uint8_t error = 0;
                alea_multilevel_hint_t curr_hint;

                /* Adjacency-cached cell find: always verifies hint cell with
                 * alea_contains_point, then falls back to BLAS on miss. */
                find_cell_multilevel(sys, x, y, z, universe_depth,
                                     &row_hint, &curr_hint,
                                     &cell_id, &material_id, &error);

                out_cell_ids[idx] = cell_id;
                if (out_material_ids) out_material_ids[idx] = material_id;
                if (out_errors) out_errors[idx] = error;

                /* Periodic full-universe overlap scan. Finds geometry errors in
                 * interior pixels that the boundary-pixel second pass would miss. */
                if (out_errors && overlap_interval > 0 &&
                    (i % overlap_interval == 0)) {
                    alea_cell_hit_t hits[32];
                    int nh = alea_find_all_cells_at_point_recursive(
                                 sys, x, y, z, hits, 32);
                    if (nh > 1) {
                        int td = (universe_depth < 0)
                                     ? hits[nh - 1].depth : universe_depth;
                        int cnt = 0;
                        for (int h = 0; h < nh; h++)
                            if (hits[h].depth == td) cnt++;
                        if (cnt > 1) out_errors[idx] = GRID_ERR_OVERLAP;
                    }
                }

                row_hint = curr_hint;
            }
        }
        /* Collect per-thread CSG counters after the for barrier */
        if (stats_en) {
            alea_perf_counters_t c = alea_perf_get();
            atomic_fetch_add_explicit(&g_csg_prim_evals, c.primitive_evaluations, memory_order_relaxed);
            atomic_fetch_add_explicit(&g_csg_bool_ops,   c.boolean_operations,     memory_order_relaxed);
        }
    }
#else
    /* Sequential version with full coherence (horizontal + vertical hints) */
    alea_multilevel_hint_t* prev_row_hints = calloc(nu, sizeof(alea_multilevel_hint_t));
    alea_multilevel_hint_t* curr_row_hints = calloc(nu, sizeof(alea_multilevel_hint_t));
    if (!prev_row_hints || !curr_row_hints) {
        free(prev_row_hints);
        free(curr_row_hints);
        return -1;
    }
    for (int i = 0; i < nu; i++) {
        multilevel_hint_init(&prev_row_hints[i]);
        multilevel_hint_init(&curr_row_hints[i]);
    }

    for (int j = 0; j < nv; j++) {
        double v = v_min + (j + 0.5) * dv;
        alea_multilevel_hint_t row_hint;
        multilevel_hint_init(&row_hint);

        for (int i = 0; i < nu; i++) {
            double u = u_min + (i + 0.5) * du;
            int idx = j * nu + i;

            double x = origin[0] + u * u_axis[0] + v * v_axis[0];
            double y = origin[1] + u * u_axis[1] + v * v_axis[1];
            double z = origin[2] + u * u_axis[2] + v * v_axis[2];

            /* Choose best hint: prefer horizontal neighbor, then vertical */
            const alea_multilevel_hint_t* hint = NULL;
            if (row_hint.valid_depth >= 0) {
                hint = &row_hint;
            } else if (prev_row_hints[i].valid_depth >= 0) {
                hint = &prev_row_hints[i];
            }

            int cell_id = -1;
            int material_id = 0;
            uint8_t error = 0;

            find_cell_multilevel(sys, x, y, z, universe_depth,
                                 hint, &curr_row_hints[i],
                                 &cell_id, &material_id, &error);

            row_hint = curr_row_hints[i];

            out_cell_ids[idx] = cell_id;
            if (out_material_ids) out_material_ids[idx] = material_id;
            if (out_errors) out_errors[idx] = error;
        }

        /* Swap row buffers */
        alea_multilevel_hint_t* tmp = prev_row_hints;
        prev_row_hints = curr_row_hints;
        curr_row_hints = tmp;
    }

    free(prev_row_hints);
    free(curr_row_hints);
#endif

    /* Second pass: recheck boundary pixels for overlap detection.
     * The adjacency optimization finds one cell per pixel without checking
     * for overlaps. Re-query boundary pixels (where a 4-connected neighbor
     * has a different cell ID) using full hierarchy search.
     *
     * We use the recursive path (bypassing the spatial coherence cache)
     * because the cache only stores one cell path and would mask overlapping
     * cells at the same depth level. */
    if (out_errors) {
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 4)
#endif
        for (int j = 0; j < nv; j++) {
            double v = v_min + (j + 0.5) * dv;
            for (int i = 0; i < nu; i++) {
                int idx = j * nu + i;
                if (out_errors[idx] != GRID_ERR_NONE) continue;

                int cell = out_cell_ids[idx];
                int is_boundary = 0;
                if (i + 1 < nu && out_cell_ids[idx + 1] != cell) is_boundary = 1;
                if (j + 1 < nv && out_cell_ids[(j+1)*nu + i] != cell) is_boundary = 1;
                if (i > 0 && out_cell_ids[idx - 1] != cell) is_boundary = 1;
                if (j > 0 && out_cell_ids[(j-1)*nu + i] != cell) is_boundary = 1;
                if (!is_boundary) continue;

                double u = u_min + (i + 0.5) * du;
                double gx = origin[0] + u * u_axis[0] + v * v_axis[0];
                double gy = origin[1] + u * u_axis[1] + v * v_axis[1];
                double gz = origin[2] + u * u_axis[2] + v * v_axis[2];

                alea_cell_hit_t hits[32];
                int num_hits = alea_find_all_cells_at_point_recursive(sys, gx, gy, gz, hits, 32);
                if (num_hits <= 1) continue;

                int target_depth = (universe_depth < 0)
                    ? hits[num_hits - 1].depth : universe_depth;
                int count = 0;
                for (int h = 0; h < num_hits; h++)
                    if (hits[h].depth == target_depth) count++;
                if (count > 1) out_errors[idx] = GRID_ERR_OVERLAP;
            }
        }
        filter_grid_overlap_ambiguities(sys, view, nu, nv, universe_depth,
                                        out_cell_ids, NULL, NULL, out_errors,
                                        NULL);
    }

    /* Print query stats when ALEA_GRID_STATS=1 */
    if (getenv("ALEA_GRID_STATS")) {
        size_t total    = atomic_load(&g_px_total);
        size_t hint_f   = atomic_load(&g_px_hint_full);
        size_t hint_p   = atomic_load(&g_px_hint_partial);
        size_t fallback = atomic_load(&g_px_full_fallback);
        size_t uh_cell  = atomic_load(&g_univ_hint_cell);
        size_t uh_nbr   = atomic_load(&g_univ_neighbor);
        size_t uh_blas  = atomic_load(&g_univ_blas);
        size_t uh_lin   = atomic_load(&g_univ_linear);
        size_t univ_total = uh_cell + uh_nbr + uh_blas + uh_lin;
        fprintf(stdout, "\n[GRID STATS] pixels=%zu\n", total);
        fprintf(stdout, "  adjacency full:    %6zu (%5.1f%%)\n", hint_f,   total ? 100.0*hint_f/total   : 0.0);
        fprintf(stdout, "  adjacency partial: %6zu (%5.1f%%)\n", hint_p,   total ? 100.0*hint_p/total   : 0.0);
        fprintf(stdout, "  full fallback:     %6zu (%5.1f%%)\n", fallback, total ? 100.0*fallback/total : 0.0);
        fprintf(stdout, "[UNIV LOOKUPS] total=%zu per_pixel=%.1f\n",
                univ_total, total ? (double)univ_total/total : 0.0);
        fprintf(stdout, "  hint cell:  %6zu (%5.1f%%)\n", uh_cell, univ_total ? 100.0*uh_cell/univ_total : 0.0);
        fprintf(stdout, "  neighbor:   %6zu (%5.1f%%)\n", uh_nbr,  univ_total ? 100.0*uh_nbr/univ_total  : 0.0);
        fprintf(stdout, "  BLAS:       %6zu (%5.1f%%)\n", uh_blas, univ_total ? 100.0*uh_blas/univ_total : 0.0);
        fprintf(stdout, "  linear:     %6zu (%5.1f%%)\n", uh_lin,  univ_total ? 100.0*uh_lin/univ_total  : 0.0);
        size_t prim = atomic_load(&g_csg_prim_evals);
        size_t bops = atomic_load(&g_csg_bool_ops);
        fprintf(stdout, "[CSG EVALS] prim=%zu bool_ops=%zu per_pixel=%.0f prim+bool=%.0f\n",
                prim, bops, total ? (double)(prim+bops)/total : 0.0,
                total ? (double)(prim+bops)/total : 0.0);
        fflush(stdout);
        /* Reset for next call */
        atomic_store(&g_px_total, 0); atomic_store(&g_px_hint_full, 0);
        atomic_store(&g_px_hint_partial, 0); atomic_store(&g_px_full_fallback, 0);
        atomic_store(&g_univ_hint_cell, 0); atomic_store(&g_univ_neighbor, 0);
        atomic_store(&g_univ_blas, 0); atomic_store(&g_univ_linear, 0);
        atomic_store(&g_csg_prim_evals, 0); atomic_store(&g_csg_bool_ops, 0);
    }

    return 0;
}

static int alea_find_cells_grid_with_paths(alea_system_t* sys,
                                           const alea_slice_view_t* view,
                                           int nu, int nv,
                                           int universe_depth,
                                           int* out_cell_ids,
                                           int* out_material_ids,
                                           uint8_t* out_errors,
                                           uint32_t* out_path_ids,
                                           alea_slice_path_table_t* out_paths) {
    if (!sys || !view || !out_cell_ids || !out_path_ids || !out_paths ||
        nu <= 0 || nv <= 0) {
        return -1;
    }

    const alea_slice_plane_t* plane = &view->plane;
    double u_min = view->u_min;
    double u_max = view->u_max;
    double v_min = view->v_min;
    double v_max = view->v_max;

    if (alea_system_prepare_query_caches(sys, ALEA_CACHE_RAYCAST) != 0)
        return -1;

    double du = (u_max - u_min) / nu;
    double dv = (v_max - v_min) / nv;
    const double* origin = plane->origin;
    const double* u_axis = plane->u_axis;
    const double* v_axis = plane->v_axis;

    size_t n = (size_t)nu * (size_t)nv;
    if (out_errors) memset(out_errors, 0, n);
    for (size_t i = 0; i < n; i++) out_path_ids[i] = UINT32_MAX;
    alea_slice_path_table_free(out_paths);

    if (!alea_spatial_mode_is_hierarchical(sys)) {
        return alea_find_cells_grid(sys, view, nu, nv, universe_depth,
                                    out_cell_ids, out_material_ids,
                                    out_errors);
    }

    alea_multilevel_hint_t* prev_row_hints =
        calloc((size_t)nu, sizeof(*prev_row_hints));
    alea_multilevel_hint_t* curr_row_hints =
        calloc((size_t)nu, sizeof(*curr_row_hints));
    if (!prev_row_hints || !curr_row_hints) {
        free(prev_row_hints);
        free(curr_row_hints);
        return -1;
    }
    for (int i = 0; i < nu; i++) {
        multilevel_hint_init(&prev_row_hints[i]);
        multilevel_hint_init(&curr_row_hints[i]);
    }

    for (int j = 0; j < nv; j++) {
        double v = v_min + (j + 0.5) * dv;
        alea_multilevel_hint_t row_hint;
        multilevel_hint_init(&row_hint);

        for (int i = 0; i < nu; i++) {
            double u = u_min + (i + 0.5) * du;
            size_t idx = (size_t)j * (size_t)nu + (size_t)i;
            double x = origin[0] + u * u_axis[0] + v * v_axis[0];
            double y = origin[1] + u * u_axis[1] + v * v_axis[1];
            double z = origin[2] + u * u_axis[2] + v * v_axis[2];

            const alea_multilevel_hint_t* hint = NULL;
            if (row_hint.valid_depth >= 0) {
                hint = &row_hint;
            } else if (prev_row_hints[i].valid_depth >= 0) {
                hint = &prev_row_hints[i];
            }

            int cell_id = -1;
            int material_id = 0;
            uint8_t error = 0;
            find_cell_multilevel(sys, x, y, z, universe_depth,
                                 hint, &curr_row_hints[i],
                                 &cell_id, &material_id, &error);

            row_hint = curr_row_hints[i];

            out_cell_ids[idx] = cell_id;
            if (out_material_ids) out_material_ids[idx] = material_id;
            if (out_errors) out_errors[idx] = error;

            if (cell_id >= 0) {
                uint32_t path_id = UINT32_MAX;
                if (slice_path_table_intern_hint(out_paths,
                                                 &curr_row_hints[i],
                                                 universe_depth,
                                                 &path_id) != 0) {
                    free(prev_row_hints);
                    free(curr_row_hints);
                    return -1;
                }
                out_path_ids[idx] = path_id;
            }
        }

        alea_multilevel_hint_t* tmp = prev_row_hints;
        prev_row_hints = curr_row_hints;
        curr_row_hints = tmp;
    }

    free(prev_row_hints);
    free(curr_row_hints);

    if (out_errors) {
        for (int j = 0; j < nv; j++) {
            double v = v_min + (j + 0.5) * dv;
            for (int i = 0; i < nu; i++) {
                size_t idx = (size_t)j * (size_t)nu + (size_t)i;
                if (out_errors[idx] != GRID_ERR_NONE) continue;

                int cell = out_cell_ids[idx];
                int is_boundary = 0;
                if (i + 1 < nu && out_cell_ids[idx + 1] != cell) is_boundary = 1;
                if (j + 1 < nv && out_cell_ids[(j + 1) * nu + i] != cell) is_boundary = 1;
                if (i > 0 && out_cell_ids[idx - 1] != cell) is_boundary = 1;
                if (j > 0 && out_cell_ids[(j - 1) * nu + i] != cell) is_boundary = 1;
                if (!is_boundary) continue;

                double u = u_min + (i + 0.5) * du;
                double gx = origin[0] + u * u_axis[0] + v * v_axis[0];
                double gy = origin[1] + u * u_axis[1] + v * v_axis[1];
                double gz = origin[2] + u * u_axis[2] + v * v_axis[2];

                alea_cell_hit_t hits[32];
                int num_hits = alea_find_all_cells_at_point_recursive(
                    sys, gx, gy, gz, hits, 32);
                if (num_hits <= 1) continue;

                int target_depth = (universe_depth < 0)
                    ? hits[num_hits - 1].depth : universe_depth;
                int count = 0;
                for (int h = 0; h < num_hits; h++)
                    if (hits[h].depth == target_depth) count++;
                if (count > 1) out_errors[idx] = GRID_ERR_OVERLAP;
            }
        }
        filter_grid_overlap_ambiguities(sys, view, nu, nv, universe_depth,
                                        out_cell_ids, NULL, NULL, out_errors,
                                        NULL);
    }

    return 0;
}


/* ============================================================================
 * FULL-GRID OVERLAP CHECK
 * ============================================================================ */

int alea_check_grid_overlaps(alea_system_t* sys,
                                  const alea_slice_view_t* view,
                                  int nu, int nv,
                                  int universe_depth,
                                  const int* cell_ids,
                                  uint8_t* errors) {
    if (!sys || !view || !cell_ids || !errors || nu <= 0 || nv <= 0) {
        return -1;
    }

    const alea_slice_plane_t* plane = &view->plane;
    double u_min = view->u_min;
    double u_max = view->u_max;
    double v_min = view->v_min;
    double v_max = view->v_max;

    double du = (u_max - u_min) / nu;
    double dv = (v_max - v_min) / nv;

    const double* origin = plane->origin;
    const double* u_axis = plane->u_axis;
    const double* v_axis = plane->v_axis;

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 4)
#endif
    for (int j = 0; j < nv; j++) {
        double v = v_min + (j + 0.5) * dv;
        for (int i = 0; i < nu; i++) {
            int idx = j * nu + i;

            /* Skip void pixels and already-flagged pixels */
            if (cell_ids[idx] < 0) continue;
            if (errors[idx] != 0) continue;

            double u = u_min + (i + 0.5) * du;
            double gx = origin[0] + u * u_axis[0] + v * v_axis[0];
            double gy = origin[1] + u * u_axis[1] + v * v_axis[1];
            double gz = origin[2] + u * u_axis[2] + v * v_axis[2];

            alea_cell_hit_t hits[32];
            int num_hits = alea_find_all_cells_at_point_recursive(sys, gx, gy, gz, hits, 32);
            if (num_hits <= 1) continue;

            int target_depth = (universe_depth < 0)
                ? hits[num_hits - 1].depth : universe_depth;
            int count = 0;
            for (int h = 0; h < num_hits; h++)
                if (hits[h].depth == target_depth) count++;
            if (count > 1) errors[idx] = GRID_ERR_OVERLAP;
        }
    }

    return 0;
}

int alea_find_cells_grid_coverage(alea_system_t* sys,
                                  const alea_slice_view_t* view,
                                  int nu, int nv,
                                  int universe_depth,
                                  unsigned flags,
                                  int* out_cell_ids,
                                  int* out_material_ids,
                                  int* out_secondary_cell_ids,
                                  uint8_t* out_coverage,
                                  uint8_t* out_errors) {
    if (!sys || !view || !out_cell_ids || nu <= 0 || nv <= 0) {
        return -1;
    }

    int rc = alea_find_cells_grid(sys, view, nu, nv, universe_depth,
                                  out_cell_ids, out_material_ids, out_errors);
    if (rc != 0) return rc;

    size_t n = (size_t)nu * (size_t)nv;
    if (out_secondary_cell_ids) {
        for (size_t i = 0; i < n; i++) out_secondary_cell_ids[i] = -1;
    }

    if (out_coverage) {
        for (size_t i = 0; i < n; i++) {
            if (out_errors && out_errors[i] == GRID_ERR_OVERLAP) {
                out_coverage[i] = ALEA_COVERAGE_MULTI;
            } else if (out_cell_ids[i] < 0 ||
                       (out_errors && out_errors[i] == GRID_ERR_UNDEFINED)) {
                out_coverage[i] = ALEA_COVERAGE_NONE;
            } else {
                out_coverage[i] = ALEA_COVERAGE_ONE;
            }
        }
    }

    if ((flags & ALEA_GRID_COVERAGE_EXACT) && (out_coverage || out_secondary_cell_ids || out_errors)) {
        alea_point_coverage_stats_reset();
        const alea_slice_plane_t* plane = &view->plane;
        double du = (view->u_max - view->u_min) / nu;
        double dv = (view->v_max - view->v_min) / nv;
        const double* origin = plane->origin;
        const double* u_axis = plane->u_axis;
        const double* v_axis = plane->v_axis;

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 4)
#endif
        for (int j = 0; j < nv; j++) {
            double v = view->v_min + (j + 0.5) * dv;
            for (int i = 0; i < nu; i++) {
                size_t idx = (size_t)j * (size_t)nu + (size_t)i;
                double u = view->u_min + (i + 0.5) * du;
                double gx = origin[0] + u * u_axis[0] + v * v_axis[0];
                double gy = origin[1] + u * u_axis[1] + v * v_axis[1];
                double gz = origin[2] + u * u_axis[2] + v * v_axis[2];

                point_coverage_t pc;
                if (find_point_coverage_exact(sys, gx, gy, gz,
                                              universe_depth, &pc) != 0) {
                    continue;
                }

                if (pc.coverage == ALEA_COVERAGE_NONE) {
                    if (out_coverage) out_coverage[idx] = ALEA_COVERAGE_NONE;
                    if (out_secondary_cell_ids) out_secondary_cell_ids[idx] = -1;
                    if (out_errors) out_errors[idx] = GRID_ERR_UNDEFINED;
                    continue;
                }

                if (pc.coverage == ALEA_COVERAGE_ONE) {
                    if (out_coverage) out_coverage[idx] = ALEA_COVERAGE_ONE;
                    if (out_secondary_cell_ids) out_secondary_cell_ids[idx] = -1;
                    if (out_errors && out_errors[idx] == GRID_ERR_OVERLAP)
                        out_errors[idx] = GRID_ERR_NONE;
                } else {
                    if (out_coverage) out_coverage[idx] = ALEA_COVERAGE_MULTI;
                    if (out_secondary_cell_ids)
                        out_secondary_cell_ids[idx] = pc.secondary_cell_id;
                    if (out_errors) out_errors[idx] = GRID_ERR_OVERLAP;
                }
            }
        }
    }

    return 0;
}

int alea_find_cells_grid_coverage_paths(alea_system_t* sys,
                                        const alea_slice_view_t* view,
                                        int nu, int nv,
                                        int universe_depth,
                                        unsigned flags,
                                        int* out_cell_ids,
                                        int* out_material_ids,
                                        int* out_secondary_cell_ids,
                                        uint8_t* out_coverage,
                                        uint8_t* out_errors,
                                        uint32_t* out_path_ids,
                                        alea_slice_path_table_t* out_paths) {
    if (!sys || !view || !out_cell_ids || nu <= 0 || nv <= 0)
        return -1;

    if ((flags & ALEA_GRID_PATH_IDS) == 0 || !out_path_ids || !out_paths) {
        return alea_find_cells_grid_coverage(
            sys, view, nu, nv, universe_depth, flags,
            out_cell_ids, out_material_ids, out_secondary_cell_ids,
            out_coverage, out_errors);
    }

    int rc = alea_find_cells_grid_with_paths(
        sys, view, nu, nv, universe_depth,
        out_cell_ids, out_material_ids, out_errors,
        out_path_ids, out_paths);
    if (rc != 0) return rc;

    size_t n = (size_t)nu * (size_t)nv;
    if (out_secondary_cell_ids) {
        for (size_t i = 0; i < n; i++) out_secondary_cell_ids[i] = -1;
    }

    if (out_coverage) {
        for (size_t i = 0; i < n; i++) {
            if (out_errors && out_errors[i] == GRID_ERR_OVERLAP) {
                out_coverage[i] = ALEA_COVERAGE_MULTI;
            } else if (out_cell_ids[i] < 0 ||
                       (out_errors && out_errors[i] == GRID_ERR_UNDEFINED)) {
                out_coverage[i] = ALEA_COVERAGE_NONE;
            } else {
                out_coverage[i] = ALEA_COVERAGE_ONE;
            }
        }
    }

    if ((flags & ALEA_GRID_COVERAGE_EXACT) &&
        (out_coverage || out_secondary_cell_ids || out_errors)) {
        alea_point_coverage_stats_reset();
        const alea_slice_plane_t* plane = &view->plane;
        double du = (view->u_max - view->u_min) / nu;
        double dv = (view->v_max - view->v_min) / nv;
        const double* origin = plane->origin;
        const double* u_axis = plane->u_axis;
        const double* v_axis = plane->v_axis;

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 4)
#endif
        for (int j = 0; j < nv; j++) {
            double v = view->v_min + (j + 0.5) * dv;
            for (int i = 0; i < nu; i++) {
                size_t idx = (size_t)j * (size_t)nu + (size_t)i;
                double u = view->u_min + (i + 0.5) * du;
                double gx = origin[0] + u * u_axis[0] + v * v_axis[0];
                double gy = origin[1] + u * u_axis[1] + v * v_axis[1];
                double gz = origin[2] + u * u_axis[2] + v * v_axis[2];

                point_coverage_t pc;
                if (find_point_coverage_exact(sys, gx, gy, gz,
                                              universe_depth, &pc) != 0) {
                    continue;
                }

                if (pc.coverage == ALEA_COVERAGE_NONE) {
                    if (out_coverage) out_coverage[idx] = ALEA_COVERAGE_NONE;
                    if (out_secondary_cell_ids) out_secondary_cell_ids[idx] = -1;
                    if (out_errors) out_errors[idx] = GRID_ERR_UNDEFINED;
                    continue;
                }

                if (pc.coverage == ALEA_COVERAGE_ONE) {
                    if (out_coverage) out_coverage[idx] = ALEA_COVERAGE_ONE;
                    if (out_secondary_cell_ids) out_secondary_cell_ids[idx] = -1;
                    if (out_errors && out_errors[idx] == GRID_ERR_OVERLAP)
                        out_errors[idx] = GRID_ERR_NONE;
                } else {
                    if (out_coverage) out_coverage[idx] = ALEA_COVERAGE_MULTI;
                    if (out_secondary_cell_ids)
                        out_secondary_cell_ids[idx] = pc.secondary_cell_id;
                    if (out_errors) out_errors[idx] = GRID_ERR_OVERLAP;
                }
            }
        }
    }

    return 0;
}

/* ============================================================================
 * CURVE-GUIDED OVERLAP DETECTION
 * ============================================================================ */

/**
 * Check a single pixel for nested overlap via full hierarchy query.
 * Returns 1 if overlap found and error grid updated, 0 otherwise.
 */
static int probe_pixel_for_overlap(alea_system_t* sys,
                                   const alea_slice_plane_t* plane,
                                   double u_min, double v_min,
                                   double du, double dv,
                                   int nu, int nv,
                                   int universe_depth,
                                   const int* cell_ids,
                                   uint8_t* errors,
                                   int pi, int pj) {
    if (pi < 0 || pi >= nu || pj < 0 || pj >= nv) return 0;
    int idx = pj * nu + pi;
    if (cell_ids[idx] < 0) return 0;    /* void */
    if (errors[idx] != GRID_ERR_NONE) return 0;  /* already flagged */

    double u = u_min + (pi + 0.5) * du;
    double v = v_min + (pj + 0.5) * dv;

    double gx = plane->origin[0] + u * plane->u_axis[0] + v * plane->v_axis[0];
    double gy = plane->origin[1] + u * plane->u_axis[1] + v * plane->v_axis[1];
    double gz = plane->origin[2] + u * plane->u_axis[2] + v * plane->v_axis[2];

    alea_cell_hit_t hits[32];
    int num_hits = alea_find_all_cells_at_point_recursive(sys, gx, gy, gz, hits, 32);
    if (num_hits <= 1) return 0;

    int target_depth = (universe_depth < 0)
        ? hits[num_hits - 1].depth : universe_depth;
    int count = 0;
    for (int h = 0; h < num_hits; h++)
        if (hits[h].depth == target_depth) count++;
    if (count > 1) {
        errors[idx] = GRID_ERR_OVERLAP;
        return 1;
    }
    return 0;
}

/**
 * Check if pixel (pi,pj) is a cell-ID boundary in the grid.
 */
static bool pixel_is_boundary(const int* cell_ids, int nu, int nv, int pi, int pj) {
    int idx = pj * nu + pi;
    int cell = cell_ids[idx];
    if (pi + 1 < nu && cell_ids[idx + 1]       != cell) return true;
    if (pj + 1 < nv && cell_ids[(pj+1)*nu + pi] != cell) return true;
    if (pi > 0     && cell_ids[idx - 1]         != cell) return true;
    if (pj > 0     && cell_ids[(pj-1)*nu + pi] != cell) return true;
    return false;
}

static int probe_curve_grid_point_for_overlap(alea_system_t* sys,
                                              const alea_slice_plane_t* plane,
                                              double u_min, double v_min,
                                              double du, double dv,
                                              double inv_du, double inv_dv,
                                              int nu, int nv,
                                              int universe_depth,
                                              const int* cell_ids,
                                              uint8_t* errors,
                                              uint8_t* probed,
                                              double u, double v) {
    int pi = (int)((u - u_min) * inv_du);
    int pj = (int)((v - v_min) * inv_dv);
    if (pi < 0 || pi >= nu || pj < 0 || pj >= nv) return 0;

    /* Skip boundary pixels: alea_find_cells_grid() already rechecks them
     * with a full hierarchy query. The curve refinement is for same-cell
     * pixels where nested overlap can hide. */
    if (pixel_is_boundary(cell_ids, nu, nv, pi, pj)) return 0;

    int new_overlaps = 0;
    for (int dj = -1; dj <= 1; dj++) {
        for (int di = -1; di <= 1; di++) {
            int qi = pi + di, qj = pj + dj;
            if (qi < 0 || qi >= nu || qj < 0 || qj >= nv) continue;

            size_t bit_idx = (size_t)qj * (size_t)nu + (size_t)qi;
            size_t byte_idx = bit_idx / 8;
            uint8_t bit_mask = (uint8_t)(1u << (bit_idx % 8));
            if (probed[byte_idx] & bit_mask) continue;
            probed[byte_idx] |= bit_mask;

            if (pixel_is_boundary(cell_ids, nu, nv, qi, qj)) continue;

            new_overlaps += probe_pixel_for_overlap(
                sys, plane, u_min, v_min, du, dv,
                nu, nv, universe_depth, cell_ids, errors, qi, qj);
        }
    }

    return new_overlaps;
}

static int probe_curve_grid_point_from_coverage(double u_min, double v_min,
                                                double inv_du, double inv_dv,
                                                int nu, int nv,
                                                const uint8_t* coverage,
                                                uint8_t* errors,
                                                uint8_t* probed,
                                                double u, double v) {
    int pi = (int)((u - u_min) * inv_du);
    int pj = (int)((v - v_min) * inv_dv);
    if (pi < 0 || pi >= nu || pj < 0 || pj >= nv) return 0;

    int new_overlaps = 0;
    for (int dj = -1; dj <= 1; dj++) {
        for (int di = -1; di <= 1; di++) {
            int qi = pi + di, qj = pj + dj;
            if (qi < 0 || qi >= nu || qj < 0 || qj >= nv) continue;

            size_t bit_idx = (size_t)qj * (size_t)nu + (size_t)qi;
            size_t byte_idx = bit_idx / 8;
            uint8_t bit_mask = (uint8_t)(1u << (bit_idx % 8));
            if (probed[byte_idx] & bit_mask) continue;
            probed[byte_idx] |= bit_mask;

            size_t idx = (size_t)qj * (size_t)nu + (size_t)qi;
            if (coverage[idx] == ALEA_COVERAGE_MULTI &&
                errors[idx] != GRID_ERR_OVERLAP) {
                errors[idx] = GRID_ERR_OVERLAP;
                new_overlaps++;
            }
        }
    }

    return new_overlaps;
}

static int find_point_coverage_from_hits(const alea_cell_hit_t* hits,
                                         int num_hits,
                                         int universe_depth,
                                         point_coverage_t* out) {
    out->coverage = ALEA_COVERAGE_NONE;
    out->primary_cell_id = -1;
    out->secondary_cell_id = -1;
    if (num_hits <= 0) return 0;

    int target_depth = (universe_depth < 0)
        ? hits[num_hits - 1].depth : universe_depth;
    int count = 0;
    for (int h = 0; h < num_hits; h++) {
        if (hits[h].depth != target_depth) continue;
        if (count == 0) out->primary_cell_id = hits[h].cell_id;
        if (count == 1) out->secondary_cell_id = hits[h].cell_id;
        count++;
        if (count > 1) {
            out->coverage = ALEA_COVERAGE_MULTI;
            return 0;
        }
    }

    out->coverage = (count == 1) ? ALEA_COVERAGE_ONE : ALEA_COVERAGE_NONE;
    return 0;
}

static int find_point_coverage_spatial(alea_system_t* sys,
                                       double gx, double gy, double gz,
                                       int universe_depth,
                                       point_coverage_t* out) {
    out->coverage = ALEA_COVERAGE_NONE;
    out->primary_cell_id = -1;
    out->secondary_cell_id = -1;

    /* Stats counters below are mutated from inside the omp parallel for in
     * alea_find_cells_grid_coverage. Every mutation must be guarded; plain
     * `++` and `+=` race and lose increments under threading. */
    if (sys->has_lattice) {
        #pragma omp atomic
        g_point_coverage_stats.lattice_fallbacks++;
        return -2;
    }

    if (alea_spatial_mode_is_hierarchical(sys)) {
        alea_cell_hit_t hits[32];
        int n = alea_hier_spatial_find_cells_at_point_uncached(sys, gx, gy, gz,
                                                               hits, 32);
        if (n < 0) return -1;
        if (n >= 32) {
            #pragma omp atomic
            g_point_coverage_stats.truncated_fallbacks++;
            return -2;
        }
        #pragma omp atomic
        g_point_coverage_stats.spatial_queries++;
        #pragma omp atomic
        g_point_coverage_stats.candidate_total += (size_t)n;
        #pragma omp critical(pc_stats_candidate_max)
        {
            if ((size_t)n > g_point_coverage_stats.candidate_max)
                g_point_coverage_stats.candidate_max = (size_t)n;
        }
        return find_point_coverage_from_hits(hits, n, universe_depth, out);
    }

    const size_t cap = 4096;
    alea_spatial_hit_t* hits = malloc(cap * sizeof(*hits));
    if (!hits) return -1;

    int n = alea_spatial_query_point(sys, gx, gy, gz, hits, cap);
    if (n < 0 || (size_t)n >= cap) {
        if ((size_t)n >= cap) {
            #pragma omp atomic
            g_point_coverage_stats.truncated_fallbacks++;
        }
        free(hits);
        return -2;
    }

    n = dedup_spatial_hits(hits, n);
    #pragma omp atomic
    g_point_coverage_stats.spatial_queries++;
    #pragma omp atomic
    g_point_coverage_stats.candidate_total += (size_t)n;
    #pragma omp critical(pc_stats_candidate_max)
    {
        if ((size_t)n > g_point_coverage_stats.candidate_max)
            g_point_coverage_stats.candidate_max = (size_t)n;
    }

    if (universe_depth >= 0) {
        int count = 0;
        for (int h = 0; h < n; h++) {
            const alea_spatial_hit_t* hit = &hits[h];
            if (hit->depth != universe_depth) continue;

            double lx = gx, ly = gy, lz = gz;
            alea_matrix_transform_point_inverse(&hit->transform, &lx, &ly, &lz);
            const alea_cell_entry_t* cell = &sys->cells.data[hit->cell_index];
            #pragma omp atomic
            g_point_coverage_stats.contains_tests++;
            if (!alea_contains_point(sys, cell->root_node_id, lx, ly, lz))
                continue;

            if (count == 0) out->primary_cell_id = hit->cell_id;
            if (count == 1) out->secondary_cell_id = hit->cell_id;
            count++;
            if (count > 1) {
                out->coverage = ALEA_COVERAGE_MULTI;
                #pragma omp atomic
                g_point_coverage_stats.spatial_multi_early_exit++;
                free(hits);
                return 0;
            }
        }
        out->coverage = (count == 1) ? ALEA_COVERAGE_ONE : ALEA_COVERAGE_NONE;
        free(hits);
        return 0;
    }

    int max_depth = -1;
    int count_at_max = 0;
    int second_at_max = -1;
    for (int h = n - 1; h >= 0; h--) {
        const alea_spatial_hit_t* hit = &hits[h];
        if (max_depth >= 0 && hit->depth < max_depth)
            break;

        double lx = gx, ly = gy, lz = gz;
        alea_matrix_transform_point_inverse(&hit->transform, &lx, &ly, &lz);
        const alea_cell_entry_t* cell = &sys->cells.data[hit->cell_index];
        #pragma omp atomic
        g_point_coverage_stats.contains_tests++;
        if (!alea_contains_point(sys, cell->root_node_id, lx, ly, lz))
            continue;

        if (hit->depth > max_depth) {
            max_depth = hit->depth;
            count_at_max = 1;
            out->primary_cell_id = hit->cell_id;
            second_at_max = -1;
        } else if (hit->depth == max_depth) {
            if (count_at_max == 1) second_at_max = hit->cell_id;
            count_at_max++;
            if (count_at_max > 1) {
                out->coverage = ALEA_COVERAGE_MULTI;
                out->secondary_cell_id = second_at_max;
                #pragma omp atomic
                g_point_coverage_stats.spatial_multi_early_exit++;
                free(hits);
                return 0;
            }
        }
    }

    out->coverage = (count_at_max == 1) ? ALEA_COVERAGE_ONE : ALEA_COVERAGE_NONE;
    free(hits);
    return 0;
}

static int find_point_coverage_exact(alea_system_t* sys,
                                     double gx, double gy, double gz,
                                     int universe_depth,
                                     point_coverage_t* out) {
    #pragma omp atomic
    g_point_coverage_stats.queries++;
    int rc = find_point_coverage_spatial(sys, gx, gy, gz, universe_depth, out);
    if (rc == 0) return 0;
    if (rc != -2) {
        #pragma omp atomic
        g_point_coverage_stats.query_errors++;
        return rc;
    }

    #pragma omp atomic
    g_point_coverage_stats.recursive_fallbacks++;
    alea_cell_hit_t hits[32];
    int num_hits = alea_find_all_cells_at_point_recursive(sys, gx, gy, gz,
                                                           hits, 32);
    if (num_hits < 0) {
        #pragma omp atomic
        g_point_coverage_stats.query_errors++;
        return -1;
    }
    return find_point_coverage_from_hits(hits, num_hits, universe_depth, out);
}

static int find_point_coverage_exact_uv(alea_system_t* sys,
                                        const alea_slice_plane_t* plane,
                                        double u, double v,
                                        int universe_depth,
                                        point_coverage_t* out) {
    if (!sys || !plane || !out) return -1;
    double gx = plane->origin[0] + u * plane->u_axis[0] + v * plane->v_axis[0];
    double gy = plane->origin[1] + u * plane->u_axis[1] + v * plane->v_axis[1];
    double gz = plane->origin[2] + u * plane->u_axis[2] + v * plane->v_axis[2];
    return find_point_coverage_exact(sys, gx, gy, gz, universe_depth, out);
}

static bool point_coverage_matches_pair(const point_coverage_t* pc,
                                        int cell_a,
                                        int cell_b) {
    if (!pc || pc->coverage != ALEA_COVERAGE_MULTI)
        return false;
    return (pc->primary_cell_id == cell_a &&
            pc->secondary_cell_id == cell_b) ||
           (pc->primary_cell_id == cell_b &&
            pc->secondary_cell_id == cell_a);
}

static int find_point_coverage_exact_world(alea_system_t* sys,
                                           const double p[3],
                                           int universe_depth,
                                           point_coverage_t* out) {
    return find_point_coverage_exact(sys, p[0], p[1], p[2],
                                     universe_depth, out);
}

static int filter_grid_overlap_ambiguities(
    alea_system_t* sys,
    const alea_slice_view_t* view,
    int nu, int nv,
    int universe_depth,
    const int* primary_cell_ids,
    int* secondary_cell_ids,
    uint8_t* coverage,
    uint8_t* errors,
    alea_boundary_filter_stats_t* out_stats) {
    if (!sys || !view || !primary_cell_ids || !errors ||
        nu <= 0 || nv <= 0)
        return -1;

    alea_boundary_filter_stats_t stats;
    memset(&stats, 0, sizeof(stats));

    double du = (view->u_max - view->u_min) / (double)nu;
    double dv = (view->v_max - view->v_min) / (double)nv;
    double eps = 0.05 * fmin(fabs(du), fabs(dv));
    if (eps <= 0.0) return -1;
    double eps_n = eps;
    if (sys->config.abs_tol > 0.0 && eps_n < 10.0 * sys->config.abs_tol)
        eps_n = 10.0 * sys->config.abs_tol;

    const double offsets[4][2] = {
        { 1.0,  0.0 },
        {-1.0,  0.0 },
        { 0.0,  1.0 },
        { 0.0, -1.0 }
    };

    int suppressed = 0;
    size_t n = (size_t)nu * (size_t)nv;
    for (size_t idx = 0; idx < n; idx++) {
        if (errors[idx] != GRID_ERR_OVERLAP)
            continue;

        int primary = primary_cell_ids[idx];
        if (primary < 0) {
            stats.inconclusive++;
            continue;
        }

        int pi = (int)(idx % (size_t)nu);
        int pj = (int)(idx / (size_t)nu);
        double u = view->u_min + (pi + 0.5) * du;
        double v = view->v_min + (pj + 0.5) * dv;
        double center[3] = {
            view->plane.origin[0] + u * view->plane.u_axis[0] + v * view->plane.v_axis[0],
            view->plane.origin[1] + u * view->plane.u_axis[1] + v * view->plane.v_axis[1],
            view->plane.origin[2] + u * view->plane.u_axis[2] + v * view->plane.v_axis[2]
        };

        point_coverage_t center_pc;
        if (find_point_coverage_exact_world(sys, center, universe_depth,
                                            &center_pc) != 0) {
            stats.inconclusive++;
            continue;
        }

        if (center_pc.coverage == ALEA_COVERAGE_NONE) {
            if (coverage) coverage[idx] = ALEA_COVERAGE_NONE;
            if (secondary_cell_ids) secondary_cell_ids[idx] = -1;
            errors[idx] = GRID_ERR_UNDEFINED;
            stats.suppressed++;
            suppressed++;
            continue;
        }

        if (center_pc.coverage == ALEA_COVERAGE_ONE) {
            if (coverage) coverage[idx] = ALEA_COVERAGE_ONE;
            if (secondary_cell_ids) secondary_cell_ids[idx] = -1;
            errors[idx] = GRID_ERR_NONE;
            stats.suppressed++;
            suppressed++;
            continue;
        }

        int pair_primary = center_pc.primary_cell_id >= 0
            ? center_pc.primary_cell_id : primary;
        int secondary = center_pc.secondary_cell_id;
        if (secondary < 0)
            secondary = secondary_cell_ids ? secondary_cell_ids[idx] : -1;
        if (secondary < 0) {
            stats.inconclusive++;
            continue;
        }

        if (coverage) coverage[idx] = ALEA_COVERAGE_MULTI;
        if (secondary_cell_ids) secondary_cell_ids[idx] = secondary;

        stats.checked++;
        bool inconclusive = false;
        int matching_multi = 0;

        point_coverage_t normal_pc[2];
        int normal_rc[2];
        double p_plus[3] = {
            center[0] + eps_n * view->plane.normal[0],
            center[1] + eps_n * view->plane.normal[1],
            center[2] + eps_n * view->plane.normal[2]
        };
        double p_minus[3] = {
            center[0] - eps_n * view->plane.normal[0],
            center[1] - eps_n * view->plane.normal[1],
            center[2] - eps_n * view->plane.normal[2]
        };
        normal_rc[0] = find_point_coverage_exact_world(sys, p_plus,
                                                       universe_depth,
                                                       &normal_pc[0]);
        normal_rc[1] = find_point_coverage_exact_world(sys, p_minus,
                                                       universe_depth,
                                                       &normal_pc[1]);
        if (normal_rc[0] == 0 && normal_rc[1] == 0 &&
            normal_pc[0].coverage != ALEA_COVERAGE_NONE &&
            normal_pc[1].coverage != ALEA_COVERAGE_NONE &&
            !point_coverage_matches_pair(&normal_pc[0], pair_primary, secondary) &&
            !point_coverage_matches_pair(&normal_pc[1], pair_primary, secondary)) {
            if (coverage) coverage[idx] = ALEA_COVERAGE_ONE;
            errors[idx] = GRID_ERR_NONE;
            if (secondary_cell_ids) secondary_cell_ids[idx] = -1;
            stats.suppressed++;
            suppressed++;
            continue;
        }
        if (normal_rc[0] != 0 || normal_rc[1] != 0 ||
            normal_pc[0].coverage == ALEA_COVERAGE_NONE ||
            normal_pc[1].coverage == ALEA_COVERAGE_NONE) {
            inconclusive = true;
        }

        for (int s = 0; s < 4; s++) {
            point_coverage_t pc;
            int rc = find_point_coverage_exact_uv(
                sys, &view->plane,
                u + offsets[s][0] * eps,
                v + offsets[s][1] * eps,
                universe_depth, &pc);
            if (rc != 0 || pc.coverage == ALEA_COVERAGE_NONE) {
                inconclusive = true;
                break;
            }
            if (pc.coverage == ALEA_COVERAGE_MULTI) {
                if (point_coverage_matches_pair(&pc, pair_primary, secondary)) {
                    matching_multi++;
                }
            }
        }

        if (matching_multi >= 3) {
            stats.retained++;
            continue;
        }
        if (inconclusive) {
            stats.inconclusive++;
            continue;
        }

        if (coverage) coverage[idx] = ALEA_COVERAGE_ONE;
        errors[idx] = GRID_ERR_NONE;
        if (secondary_cell_ids) secondary_cell_ids[idx] = -1;
        stats.suppressed++;
        suppressed++;
    }

    if (out_stats) *out_stats = stats;
    return suppressed;
}

int alea_filter_grid_boundary_ambiguities(
    alea_system_t* sys,
    const alea_slice_view_t* view,
    int nu, int nv,
    int universe_depth,
    const int* primary_cell_ids,
    int* secondary_cell_ids,
    uint8_t* coverage,
    uint8_t* errors,
    alea_boundary_filter_stats_t* out_stats) {
    if (!coverage) return -1;
    return filter_grid_overlap_ambiguities(sys, view, nu, nv, universe_depth,
                                           primary_cell_ids, secondary_cell_ids,
                                           coverage, errors, out_stats);
}

static int update_pixel_coverage_exact(alea_system_t* sys,
                                       const alea_slice_plane_t* plane,
                                       double u_min, double v_min,
                                       double du, double dv,
                                       int nu, int nv,
                                       int universe_depth,
                                       int* out_secondary_cell_ids,
                                       uint8_t* coverage,
                                       uint8_t* errors,
                                       int pi, int pj) {
    if (pi < 0 || pi >= nu || pj < 0 || pj >= nv) return 0;

    size_t idx = (size_t)pj * (size_t)nu + (size_t)pi;
    double u = u_min + (pi + 0.5) * du;
    double v = v_min + (pj + 0.5) * dv;
    double gx = plane->origin[0] + u * plane->u_axis[0] + v * plane->v_axis[0];
    double gy = plane->origin[1] + u * plane->u_axis[1] + v * plane->v_axis[1];
    double gz = plane->origin[2] + u * plane->u_axis[2] + v * plane->v_axis[2];

    point_coverage_t pc;
    if (find_point_coverage_exact(sys, gx, gy, gz, universe_depth, &pc) != 0)
        return 0;

    if (pc.coverage == ALEA_COVERAGE_NONE) {
        coverage[idx] = ALEA_COVERAGE_NONE;
        if (out_secondary_cell_ids) out_secondary_cell_ids[idx] = -1;
        errors[idx] = GRID_ERR_UNDEFINED;
        return 0;
    }

    if (pc.coverage == ALEA_COVERAGE_ONE) {
        coverage[idx] = ALEA_COVERAGE_ONE;
        if (out_secondary_cell_ids) out_secondary_cell_ids[idx] = -1;
        if (errors[idx] == GRID_ERR_OVERLAP) errors[idx] = GRID_ERR_NONE;
        return 0;
    }

    int newly_multi = coverage[idx] != ALEA_COVERAGE_MULTI ||
                      errors[idx] != GRID_ERR_OVERLAP;
    coverage[idx] = ALEA_COVERAGE_MULTI;
    if (out_secondary_cell_ids) out_secondary_cell_ids[idx] = pc.secondary_cell_id;
    errors[idx] = GRID_ERR_OVERLAP;
    return newly_multi ? 1 : 0;
}

static bool slice_matrix_equal(const alea_matrix_t* a, const alea_matrix_t* b) {
    const double tol = 1e-10;
    for (int i = 0; i < 12; i++) {
        if (fabs(a->m[i] - b->m[i]) > tol) return false;
    }
    return true;
}

static int dedup_spatial_hits(alea_spatial_hit_t* hits, int hit_count) {
    int unique_count = 0;
    for (int h = 0; h < hit_count; h++) {
        bool duplicate = false;
        for (int k = 0; k < unique_count; k++) {
            if (hits[k].cell_index == hits[h].cell_index &&
                hits[k].depth == hits[h].depth &&
                slice_matrix_equal(&hits[k].transform, &hits[h].transform)) {
                duplicate = true;
                break;
            }
        }
        if (!duplicate) {
            if (unique_count != h) hits[unique_count] = hits[h];
            unique_count++;
        }
    }
    return unique_count;
}

static alea_bbox_t tile_query_bbox(const alea_slice_view_t* view,
                                   double u0, double u1,
                                   double v0, double v1,
                                   double eps) {
    const alea_slice_plane_t* plane = &view->plane;
    double uvals[2] = {u0, u1};
    double vvals[2] = {v0, v1};
    alea_bbox_t bbox = {
        .min_x = DBL_MAX, .max_x = -DBL_MAX,
        .min_y = DBL_MAX, .max_y = -DBL_MAX,
        .min_z = DBL_MAX, .max_z = -DBL_MAX
    };

    for (int ui = 0; ui < 2; ui++) {
        for (int vi = 0; vi < 2; vi++) {
            double x = plane->origin[0] +
                       uvals[ui] * plane->u_axis[0] +
                       vvals[vi] * plane->v_axis[0];
            double y = plane->origin[1] +
                       uvals[ui] * plane->u_axis[1] +
                       vvals[vi] * plane->v_axis[1];
            double z = plane->origin[2] +
                       uvals[ui] * plane->u_axis[2] +
                       vvals[vi] * plane->v_axis[2];
            if (x < bbox.min_x) bbox.min_x = x;
            if (x > bbox.max_x) bbox.max_x = x;
            if (y < bbox.min_y) bbox.min_y = y;
            if (y > bbox.max_y) bbox.max_y = y;
            if (z < bbox.min_z) bbox.min_z = z;
            if (z > bbox.max_z) bbox.max_z = z;
        }
    }

    bbox.min_x -= eps; bbox.max_x += eps;
    bbox.min_y -= eps; bbox.max_y += eps;
    bbox.min_z -= eps; bbox.max_z += eps;
    return bbox;
}

static int update_pixel_coverage_from_candidates(
    alea_system_t* sys,
    const alea_slice_plane_t* plane,
    const alea_spatial_hit_t* hits,
    int hit_count,
    double u_min, double v_min,
    double du, double dv,
    int nu,
    int universe_depth,
    int* out_secondary_cell_ids,
    uint8_t* coverage,
    uint8_t* errors,
    int pi, int pj) {
    size_t idx = (size_t)pj * (size_t)nu + (size_t)pi;
    double u = u_min + (pi + 0.5) * du;
    double v = v_min + (pj + 0.5) * dv;
    double gx = plane->origin[0] + u * plane->u_axis[0] + v * plane->v_axis[0];
    double gy = plane->origin[1] + u * plane->u_axis[1] + v * plane->v_axis[1];
    double gz = plane->origin[2] + u * plane->u_axis[2] + v * plane->v_axis[2];

    int any_inside = 0;
    int count = 0;
    int second_cell = -1;
    int count_at_max_depth = 0;
    int second_cell_at_max_depth = -1;
    int max_depth = -1;

    for (int h = 0; h < hit_count; h++) {
        const alea_spatial_hit_t* hit = &hits[h];
        double lx = gx, ly = gy, lz = gz;
        alea_matrix_transform_point_inverse(&hit->transform, &lx, &ly, &lz);

        const alea_cell_entry_t* cell = &sys->cells.data[hit->cell_index];
        if (!alea_contains_point(sys, cell->root_node_id, lx, ly, lz))
            continue;

        any_inside = 1;
        if (universe_depth >= 0) {
            if (hit->depth != universe_depth) continue;
            if (count == 1) second_cell = hit->cell_id;
            count++;
        } else if (hit->depth > max_depth) {
            max_depth = hit->depth;
            count_at_max_depth = 1;
            second_cell_at_max_depth = -1;
        } else if (hit->depth == max_depth) {
            if (count_at_max_depth == 1)
                second_cell_at_max_depth = hit->cell_id;
            count_at_max_depth++;
        }
    }

    if (!any_inside) {
        coverage[idx] = ALEA_COVERAGE_NONE;
        if (out_secondary_cell_ids) out_secondary_cell_ids[idx] = -1;
        errors[idx] = GRID_ERR_UNDEFINED;
        return 0;
    }

    if (universe_depth < 0) {
        count = count_at_max_depth;
        second_cell = second_cell_at_max_depth;
    }

    if (count <= 0) {
        coverage[idx] = ALEA_COVERAGE_NONE;
        if (out_secondary_cell_ids) out_secondary_cell_ids[idx] = -1;
        errors[idx] = GRID_ERR_UNDEFINED;
        return 0;
    }

    if (count == 1) {
        coverage[idx] = ALEA_COVERAGE_ONE;
        if (out_secondary_cell_ids) out_secondary_cell_ids[idx] = -1;
        if (errors[idx] == GRID_ERR_OVERLAP) errors[idx] = GRID_ERR_NONE;
        return 0;
    }

    int newly_multi = coverage[idx] != ALEA_COVERAGE_MULTI ||
                      errors[idx] != GRID_ERR_OVERLAP;
    coverage[idx] = ALEA_COVERAGE_MULTI;
    if (out_secondary_cell_ids) out_secondary_cell_ids[idx] = second_cell;
    errors[idx] = GRID_ERR_OVERLAP;
    return newly_multi ? 1 : 0;
}

static void path_world_to_local_point(const double m[12],
                                      double* x, double* y, double* z) {
    double px = *x, py = *y, pz = *z;
    *x = m[0] * px + m[1] * py + m[2] * pz + m[3];
    *y = m[4] * px + m[5] * py + m[6] * pz + m[7];
    *z = m[8] * px + m[9] * py + m[10] * pz + m[11];
}

static void slice_grid_world_point(const alea_slice_plane_t* plane,
                                   double u_min, double v_min,
                                   double du, double dv,
                                   int pi, int pj,
                                   double* gx, double* gy, double* gz) {
    double u = u_min + (pi + 0.5) * du;
    double v = v_min + (pj + 0.5) * dv;
    *gx = plane->origin[0] + u * plane->u_axis[0] + v * plane->v_axis[0];
    *gy = plane->origin[1] + u * plane->u_axis[1] + v * plane->v_axis[1];
    *gz = plane->origin[2] + u * plane->u_axis[2] + v * plane->v_axis[2];
}

static int update_path_group_coverage_from_candidates(
    alea_system_t* sys,
    const alea_spatial_hit_t* hits,
    int hit_count,
    size_t group_pixels,
    const size_t* pixel_indices,
    const int* primary_cell_ids,
    const double* local_x,
    const double* local_y,
    const double* local_z,
    int* count,
    int* second_cell,
    int* out_secondary_cell_ids,
    uint8_t* coverage,
    uint8_t* errors,
    size_t* out_bbox_tests,
    size_t* out_bbox_rejects,
    size_t* out_contains_tests,
    size_t* out_early_multi_skips,
    size_t* out_primary_cell_skips) {
    if (!sys || hit_count < 0 || (hit_count > 0 && !hits) ||
        !pixel_indices || !primary_cell_ids ||
        !local_x || !local_y || !local_z || !count || !second_cell ||
        !coverage || !errors)
        return -1;

    for (size_t q = 0; q < group_pixels; q++) {
        count[q] = 0;
        second_cell[q] = -1;
    }

    size_t bbox_tests = 0;
    size_t bbox_rejects = 0;
    size_t contains_tests = 0;
    size_t early_multi_skips = 0;
    size_t primary_cell_skips = 0;
    for (int h = 0; h < hit_count; h++) {
        const alea_spatial_hit_t* hit = &hits[h];
        if ((size_t)hit->cell_index >= sys->cells.count)
            continue;

        const alea_cell_entry_t* cell = &sys->cells.data[hit->cell_index];
        if (cell->root_node_id == ALEA_NODE_ID_INVALID)
            continue;

        const alea_bbox_t* bbox = &sys->nodes.data[cell->root_node_id].bbox;
        for (size_t q = 0; q < group_pixels; q++) {
            if (count[q] > 1) {
                early_multi_skips++;
                continue;
            }
            size_t idx = pixel_indices[q];
            if (hit->cell_id == primary_cell_ids[idx]) {
                primary_cell_skips++;
                continue;
            }

            double lx = local_x[q];
            double ly = local_y[q];
            double lz = local_z[q];
            bbox_tests++;
            if (lx < bbox->min_x || lx > bbox->max_x ||
                ly < bbox->min_y || ly > bbox->max_y ||
                lz < bbox->min_z || lz > bbox->max_z) {
                bbox_rejects++;
                continue;
            }

            contains_tests++;
            if (!alea_contains_point(sys, cell->root_node_id, lx, ly, lz))
                continue;

            if (count[q] == 0) second_cell[q] = hit->cell_id;
            count[q] = 2;
        }
    }

    int updated = 0;
    for (size_t q = 0; q < group_pixels; q++) {
        size_t idx = pixel_indices[q];
        if (count[q] <= 1) {
            coverage[idx] = ALEA_COVERAGE_ONE;
            if (out_secondary_cell_ids) out_secondary_cell_ids[idx] = -1;
            if (errors[idx] == GRID_ERR_OVERLAP) errors[idx] = GRID_ERR_NONE;
        } else {
            if (coverage[idx] != ALEA_COVERAGE_MULTI ||
                errors[idx] != GRID_ERR_OVERLAP)
                updated++;
            coverage[idx] = ALEA_COVERAGE_MULTI;
            if (out_secondary_cell_ids)
                out_secondary_cell_ids[idx] = second_cell[q];
            errors[idx] = GRID_ERR_OVERLAP;
        }
    }

    if (out_bbox_tests) *out_bbox_tests = bbox_tests;
    if (out_bbox_rejects) *out_bbox_rejects = bbox_rejects;
    if (out_contains_tests) *out_contains_tests = contains_tests;
    if (out_early_multi_skips) *out_early_multi_skips = early_multi_skips;
    if (out_primary_cell_skips) *out_primary_cell_skips = primary_cell_skips;
    return updated;
}

typedef struct {
    alea_spatial_hit_t* hits;
    size_t count;
    size_t capacity;
} path_tile_bucket_t;

typedef struct {
    uint32_t path_id;
    int tiles_x;
    int tiles_y;
    path_tile_bucket_t* buckets;
} path_slice_index_t;

static void path_tile_bucket_free(path_tile_bucket_t* bucket) {
    if (!bucket) return;
    free(bucket->hits);
    bucket->hits = NULL;
    bucket->count = 0;
    bucket->capacity = 0;
}

static void path_slice_index_free(path_slice_index_t* index) {
    if (!index || !index->buckets) return;
    size_t n = (size_t)index->tiles_x * (size_t)index->tiles_y;
    for (size_t i = 0; i < n; i++)
        path_tile_bucket_free(&index->buckets[i]);
    free(index->buckets);
    index->buckets = NULL;
}

static int path_tile_bucket_append(path_tile_bucket_t* bucket,
                                   const alea_spatial_hit_t* hit) {
    if (!bucket || !hit) return -1;
    if (bucket->count == bucket->capacity) {
        size_t new_cap = bucket->capacity ? bucket->capacity * 2 : 16;
        alea_spatial_hit_t* grown =
            realloc(bucket->hits, new_cap * sizeof(*grown));
        if (!grown) return -1;
        bucket->hits = grown;
        bucket->capacity = new_cap;
    }
    bucket->hits[bucket->count++] = *hit;
    return 0;
}

static bool path_hit_list_contains_cell(const alea_spatial_hit_t* hits,
                                        size_t count,
                                        uint32_t cell_index) {
    for (size_t i = 0; i < count; i++) {
        if (hits[i].cell_index == cell_index)
            return true;
    }
    return false;
}

#ifdef _OPENMP
static void tile_coverage_stats_add(alea_tile_coverage_stats_t* dst,
                                    const alea_tile_coverage_stats_t* src) {
    if (!dst || !src) return;
    dst->tiles += src->tiles;
    dst->fallback_tiles += src->fallback_tiles;
    dst->query_errors += src->query_errors;
    dst->pixels += src->pixels;
    dst->exact_fallback_pixels += src->exact_fallback_pixels;
    dst->candidate_total += src->candidate_total;
    if (src->candidate_max > dst->candidate_max)
        dst->candidate_max = src->candidate_max;
    dst->dedup_candidate_total += src->dedup_candidate_total;
    if (src->dedup_candidate_max > dst->dedup_candidate_max)
        dst->dedup_candidate_max = src->dedup_candidate_max;
    dst->candidate_pixel_tests += src->candidate_pixel_tests;
    dst->refined_pixels += src->refined_pixels;
    dst->path_groups += src->path_groups;
    if (src->path_group_pixels_max > dst->path_group_pixels_max)
        dst->path_group_pixels_max = src->path_group_pixels_max;
    if (src->path_group_candidates_max > dst->path_group_candidates_max)
        dst->path_group_candidates_max = src->path_group_candidates_max;
    dst->bbox_pixel_tests += src->bbox_pixel_tests;
    dst->bbox_pixel_rejects += src->bbox_pixel_rejects;
    dst->early_multi_skips += src->early_multi_skips;
    dst->primary_cell_skips += src->primary_cell_skips;
    dst->path_2d_verify_queries += src->path_2d_verify_queries;
    dst->path_2d_missing_candidates += src->path_2d_missing_candidates;
    dst->path_2d_missing_tiles += src->path_2d_missing_tiles;
}
#endif

static double path_dot3(double ax, double ay, double az,
                        double bx, double by, double bz) {
    return ax * bx + ay * by + az * bz;
}

static void path_transform_vector(const double m[12],
                                  double* x, double* y, double* z) {
    double px = *x, py = *y, pz = *z;
    *x = m[0] * px + m[1] * py + m[2] * pz;
    *y = m[4] * px + m[5] * py + m[6] * pz;
    *z = m[8] * px + m[9] * py + m[10] * pz;
}

static int path_slice_index_build(alea_system_t* sys,
                                  const alea_slice_view_t* view,
                                  int nu, int nv,
                                  int tile_w, int tile_h,
                                  uint32_t path_id,
                                  const alea_slice_path_record_t* path,
                                  size_t bucket_limit,
                                  int tile_pad,
                                  path_slice_index_t* out_index) {
    if (!sys || !view || !path || !out_index || tile_w <= 0 || tile_h <= 0)
        return -1;
    memset(out_index, 0, sizeof(*out_index));
    out_index->path_id = path_id;
    out_index->tiles_x = (nu + tile_w - 1) / tile_w;
    out_index->tiles_y = (nv + tile_h - 1) / tile_h;
    size_t bucket_count = (size_t)out_index->tiles_x * (size_t)out_index->tiles_y;
    out_index->buckets = calloc(bucket_count, sizeof(*out_index->buckets));
    if (!out_index->buckets) return -1;

    double o[3] = {
        view->plane.origin[0],
        view->plane.origin[1],
        view->plane.origin[2]
    };
    path_world_to_local_point(path->world_to_local, &o[0], &o[1], &o[2]);

    double u_axis[3] = {
        view->plane.u_axis[0],
        view->plane.u_axis[1],
        view->plane.u_axis[2]
    };
    double v_axis[3] = {
        view->plane.v_axis[0],
        view->plane.v_axis[1],
        view->plane.v_axis[2]
    };
    path_transform_vector(path->world_to_local, &u_axis[0], &u_axis[1], &u_axis[2]);
    path_transform_vector(path->world_to_local, &v_axis[0], &v_axis[1], &v_axis[2]);
    double uu = path_dot3(u_axis[0], u_axis[1], u_axis[2],
                          u_axis[0], u_axis[1], u_axis[2]);
    double uv = path_dot3(u_axis[0], u_axis[1], u_axis[2],
                          v_axis[0], v_axis[1], v_axis[2]);
    double vv = path_dot3(v_axis[0], v_axis[1], v_axis[2],
                          v_axis[0], v_axis[1], v_axis[2]);
    double det = uu * vv - uv * uv;
    if (fabs(det) <= 1e-30) {
        path_slice_index_free(out_index);
        return -1;
    }

    const alea_universe_t* univ = alea_get_universe(sys, path->universe_id);
    if (!univ) {
        path_slice_index_free(out_index);
        return -1;
    }

    double du = (view->u_max - view->u_min) / (double)nu;
    double dv = (view->v_max - view->v_min) / (double)nv;
    double eps = 1e-9;

    for (size_t ci = 0; ci < univ->cell_indices.count; ci++) {
        size_t cell_index = univ->cell_indices.data[ci];
        if (cell_index >= sys->cells.count) continue;
        const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
        if (cell->root_node_id == ALEA_NODE_ID_INVALID) continue;

        const alea_bbox_t* bbox = &sys->nodes.data[cell->root_node_id].bbox;
        double u_min_p = DBL_MAX, u_max_p = -DBL_MAX;
        double v_min_p = DBL_MAX, v_max_p = -DBL_MAX;
        for (int ix = 0; ix < 2; ix++) {
            double x = ix ? bbox->max_x : bbox->min_x;
            for (int iy = 0; iy < 2; iy++) {
                double y = iy ? bbox->max_y : bbox->min_y;
                for (int iz = 0; iz < 2; iz++) {
                    double z = iz ? bbox->max_z : bbox->min_z;
                    double rx = x - o[0], ry = y - o[1], rz = z - o[2];
                    double ru = path_dot3(rx, ry, rz,
                                          u_axis[0], u_axis[1], u_axis[2]);
                    double rv = path_dot3(rx, ry, rz,
                                          v_axis[0], v_axis[1], v_axis[2]);
                    double pu = (ru * vv - rv * uv) / det;
                    double pv = (rv * uu - ru * uv) / det;
                    if (pu < u_min_p) u_min_p = pu;
                    if (pu > u_max_p) u_max_p = pu;
                    if (pv < v_min_p) v_min_p = pv;
                    if (pv > v_max_p) v_max_p = pv;
                }
            }
        }
        if (u_max_p < view->u_min || u_min_p > view->u_max ||
            v_max_p < view->v_min || v_min_p > view->v_max)
            continue;

        int i0 = (int)floor((u_min_p - view->u_min - eps) / du);
        int i1 = (int)floor((u_max_p - view->u_min + eps) / du);
        int j0 = (int)floor((v_min_p - view->v_min - eps) / dv);
        int j1 = (int)floor((v_max_p - view->v_min + eps) / dv);
        if (i0 < 0) i0 = 0;
        if (j0 < 0) j0 = 0;
        if (i1 >= nu) i1 = nu - 1;
        if (j1 >= nv) j1 = nv - 1;
        if (i0 > i1 || j0 > j1) continue;

        int ti0 = i0 / tile_w, ti1 = i1 / tile_w;
        int tj0 = j0 / tile_h, tj1 = j1 / tile_h;
        ti0 -= tile_pad; ti1 += tile_pad;
        tj0 -= tile_pad; tj1 += tile_pad;
        if (ti0 < 0) ti0 = 0;
        if (tj0 < 0) tj0 = 0;
        if (ti1 >= out_index->tiles_x) ti1 = out_index->tiles_x - 1;
        if (tj1 >= out_index->tiles_y) tj1 = out_index->tiles_y - 1;
        alea_spatial_hit_t hit;
        memset(&hit, 0, sizeof(hit));
        hit.cell_index = (uint32_t)cell_index;
        hit.cell_id = cell->mc_cell_id;
        hit.material_id = cell->material_id;
        hit.universe_id = path->universe_id;
        hit.depth = 0;
        hit.is_terminal = cell->fill_universe <= 0 && cell->lat_type == 0;
        alea_matrix_identity(&hit.transform);

        for (int tj = tj0; tj <= tj1; tj++) {
            for (int ti = ti0; ti <= ti1; ti++) {
                size_t bidx = (size_t)tj * (size_t)out_index->tiles_x + (size_t)ti;
                path_tile_bucket_t* bucket = &out_index->buckets[bidx];
                if (bucket_limit > 0 && bucket->count >= bucket_limit) {
                    path_slice_index_free(out_index);
                    return -2;
                }
                if (path_tile_bucket_append(bucket, &hit) != 0) {
                    path_slice_index_free(out_index);
                    return -1;
                }
            }
        }
    }

    return 0;
}

#ifdef _OPENMP
static int refine_path_tile_3d_exact(
    alea_system_t* sys,
    const alea_slice_view_t* view,
    int nu, int nv,
    int universe_depth,
    int tile_w, int tile_h,
    int ti, int tj,
    double du, double dv,
    size_t max_hits,
    alea_spatial_hit_t* hits,
    size_t max_group_pixels,
    size_t* group_indices,
    double* group_x,
    double* group_y,
    double* group_z,
    int* group_count,
    int* group_second,
    const int* primary_cell_ids,
    const uint32_t* path_ids,
    const alea_slice_path_table_t* paths,
    int* out_secondary_cell_ids,
    uint8_t* coverage,
    uint8_t* errors,
    alea_tile_coverage_stats_t* stats) {
    if (!sys || !view || !hits || !group_indices || !group_x || !group_y ||
        !group_z || !group_count || !group_second || !primary_cell_ids ||
        !path_ids || !paths || !coverage || !errors || !stats)
        return -1;

    int i_end = ti + tile_w;
    if (i_end > nu) i_end = nu;
    int j_end = tj + tile_h;
    if (j_end > nv) j_end = nv;
    double eps = 1e-9;
    int updated = 0;
    stats->tiles++;

    for (int sj = tj; sj < j_end; sj++) {
        for (int si = ti; si < i_end; si++) {
            size_t seed_idx = (size_t)sj * (size_t)nu + (size_t)si;
            uint32_t path_id = path_ids[seed_idx];
            if (path_id == UINT32_MAX) {
                stats->fallback_tiles++;
                stats->exact_fallback_pixels++;
                updated += update_pixel_coverage_exact(
                    sys, &view->plane, view->u_min, view->v_min,
                    du, dv, nu, nv, universe_depth,
                    out_secondary_cell_ids, coverage, errors, si, sj);
                continue;
            }
            if ((size_t)path_id >= paths->count)
                return -1;

            int already_done = 0;
            for (int pj = tj; pj <= sj && !already_done; pj++) {
                int pi_limit = (pj == sj) ? si : i_end;
                for (int pi = ti; pi < pi_limit; pi++) {
                    size_t pidx = (size_t)pj * (size_t)nu + (size_t)pi;
                    if (path_ids[pidx] == path_id) {
                        already_done = 1;
                        break;
                    }
                }
            }
            if (already_done) continue;

            const alea_slice_path_record_t* path = &paths->records[path_id];
            if (universe_depth >= 0 && path->depth != universe_depth) {
                stats->fallback_tiles++;
                updated += update_pixel_coverage_exact(
                    sys, &view->plane, view->u_min, view->v_min,
                    du, dv, nu, nv, universe_depth,
                    out_secondary_cell_ids, coverage, errors, si, sj);
                continue;
            }

            alea_bbox_t local_bbox = {
                .min_x = DBL_MAX, .max_x = -DBL_MAX,
                .min_y = DBL_MAX, .max_y = -DBL_MAX,
                .min_z = DBL_MAX, .max_z = -DBL_MAX
            };
            size_t group_pixels = 0;
            for (int pj = tj; pj < j_end; pj++) {
                for (int pi = ti; pi < i_end; pi++) {
                    size_t pidx = (size_t)pj * (size_t)nu + (size_t)pi;
                    if (path_ids[pidx] != path_id) continue;
                    if (group_pixels >= max_group_pixels)
                        return -1;

                    double lx, ly, lz;
                    slice_grid_world_point(&view->plane, view->u_min,
                                           view->v_min, du, dv, pi, pj,
                                           &lx, &ly, &lz);
                    path_world_to_local_point(path->world_to_local,
                                              &lx, &ly, &lz);
                    if (lx < local_bbox.min_x) local_bbox.min_x = lx;
                    if (lx > local_bbox.max_x) local_bbox.max_x = lx;
                    if (ly < local_bbox.min_y) local_bbox.min_y = ly;
                    if (ly > local_bbox.max_y) local_bbox.max_y = ly;
                    if (lz < local_bbox.min_z) local_bbox.min_z = lz;
                    if (lz > local_bbox.max_z) local_bbox.max_z = lz;
                    group_indices[group_pixels] = pidx;
                    group_x[group_pixels] = lx;
                    group_y[group_pixels] = ly;
                    group_z[group_pixels] = lz;
                    group_pixels++;
                }
            }
            if (group_pixels == 0) continue;
            local_bbox.min_x -= eps; local_bbox.max_x += eps;
            local_bbox.min_y -= eps; local_bbox.max_y += eps;
            local_bbox.min_z -= eps; local_bbox.max_z += eps;

            int hit_count = alea_hier_spatial_query_universe_region(
                sys, path->universe_id, &local_bbox, hits, max_hits);
            if (hit_count < 0) {
                stats->query_errors++;
                return -1;
            }

            stats->candidate_total += (size_t)hit_count;
            if ((size_t)hit_count > stats->candidate_max)
                stats->candidate_max = (size_t)hit_count;

            if ((size_t)hit_count >= max_hits) {
                stats->fallback_tiles++;
                stats->exact_fallback_pixels += group_pixels;
                for (int pj = tj; pj < j_end; pj++) {
                    for (int pi = ti; pi < i_end; pi++) {
                        size_t pidx = (size_t)pj * (size_t)nu + (size_t)pi;
                        if (path_ids[pidx] != path_id) continue;
                        updated += update_pixel_coverage_exact(
                            sys, &view->plane, view->u_min, view->v_min,
                            du, dv, nu, nv, universe_depth,
                            out_secondary_cell_ids, coverage, errors, pi, pj);
                    }
                }
                continue;
            }

            hit_count = dedup_spatial_hits(hits, hit_count);
            stats->path_groups++;
            if (group_pixels > stats->path_group_pixels_max)
                stats->path_group_pixels_max = group_pixels;
            if ((size_t)hit_count > stats->path_group_candidates_max)
                stats->path_group_candidates_max = (size_t)hit_count;
            stats->dedup_candidate_total += (size_t)hit_count;
            if ((size_t)hit_count > stats->dedup_candidate_max)
                stats->dedup_candidate_max = (size_t)hit_count;
            stats->pixels += group_pixels;

            size_t bbox_tests = 0;
            size_t bbox_rejects = 0;
            size_t contains_tests = 0;
            size_t early_multi_skips = 0;
            size_t primary_cell_skips = 0;
            int rc = update_path_group_coverage_from_candidates(
                sys, hits, hit_count, group_pixels,
                group_indices, primary_cell_ids,
                group_x, group_y, group_z,
                group_count, group_second, out_secondary_cell_ids,
                coverage, errors, &bbox_tests, &bbox_rejects,
                &contains_tests, &early_multi_skips,
                &primary_cell_skips);
            if (rc < 0)
                return -1;

            stats->bbox_pixel_tests += bbox_tests;
            stats->bbox_pixel_rejects += bbox_rejects;
            stats->candidate_pixel_tests += contains_tests;
            stats->early_multi_skips += early_multi_skips;
            stats->primary_cell_skips += primary_cell_skips;
            updated += rc;
        }
    }

    return updated;
}
#endif

static int chain_hit_contains_point(alea_system_t* sys,
                                    const alea_hier_spatial_chain_hit_t* hit,
                                    double gx, double gy, double gz) {
    if (hit->chain_truncated) return -1;

    for (uint8_t a = 0; a < hit->ancestor_count; a++) {
        uint32_t cell_index = hit->ancestor_cell_indices[a];
        if (cell_index >= sys->cells.count) return -1;
        const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
        double lx = gx, ly = gy, lz = gz;
        alea_matrix_transform_point_inverse(&hit->ancestor_transforms[a],
                                            &lx, &ly, &lz);
        if (!alea_contains_point(sys, cell->root_node_id, lx, ly, lz))
            return 0;
    }

    const alea_spatial_hit_t* terminal = &hit->hit;
    if (terminal->cell_index >= sys->cells.count) return -1;
    const alea_cell_entry_t* cell = &sys->cells.data[terminal->cell_index];
    double lx = gx, ly = gy, lz = gz;
    alea_matrix_transform_point_inverse(&terminal->transform, &lx, &ly, &lz);
    return alea_contains_point(sys, cell->root_node_id, lx, ly, lz) ? 1 : 0;
}

static int update_pixel_coverage_from_chain_candidates(
    alea_system_t* sys,
    const alea_slice_plane_t* plane,
    const alea_hier_spatial_chain_hit_t* hits,
    int hit_count,
    double u_min, double v_min,
    double du, double dv,
    int nu,
    int universe_depth,
    int* out_secondary_cell_ids,
    uint8_t* coverage,
    uint8_t* errors,
    int pi, int pj) {
    size_t idx = (size_t)pj * (size_t)nu + (size_t)pi;
    double u = u_min + (pi + 0.5) * du;
    double v = v_min + (pj + 0.5) * dv;
    double gx = plane->origin[0] + u * plane->u_axis[0] + v * plane->v_axis[0];
    double gy = plane->origin[1] + u * plane->u_axis[1] + v * plane->v_axis[1];
    double gz = plane->origin[2] + u * plane->u_axis[2] + v * plane->v_axis[2];

    int count = 0;
    int second_cell = -1;
    int count_at_max_depth = 0;
    int second_cell_at_max_depth = -1;
    int max_depth = -1;

    for (int h = 0; h < hit_count; h++) {
        const alea_hier_spatial_chain_hit_t* hit = &hits[h];
        int inside = chain_hit_contains_point(sys, hit, gx, gy, gz);
        if (inside < 0) return -1;
        if (!inside) continue;

        const alea_spatial_hit_t* sh = &hit->hit;
        if (universe_depth >= 0) {
            if (sh->depth != universe_depth) continue;
            if (count == 1) second_cell = sh->cell_id;
            count++;
            if (count > 1) break;
        } else if (sh->depth > max_depth) {
            max_depth = sh->depth;
            count_at_max_depth = 1;
            second_cell_at_max_depth = -1;
        } else if (sh->depth == max_depth) {
            if (count_at_max_depth == 1)
                second_cell_at_max_depth = sh->cell_id;
            count_at_max_depth++;
            if (count_at_max_depth > 1) break;
        }
    }

    if (universe_depth < 0) {
        count = count_at_max_depth;
        second_cell = second_cell_at_max_depth;
    }

    if (count <= 0) {
        coverage[idx] = ALEA_COVERAGE_NONE;
        if (out_secondary_cell_ids) out_secondary_cell_ids[idx] = -1;
        errors[idx] = GRID_ERR_UNDEFINED;
        return 0;
    }

    if (count == 1) {
        coverage[idx] = ALEA_COVERAGE_ONE;
        if (out_secondary_cell_ids) out_secondary_cell_ids[idx] = -1;
        if (errors[idx] == GRID_ERR_OVERLAP) errors[idx] = GRID_ERR_NONE;
        return 0;
    }

    int newly_multi = coverage[idx] != ALEA_COVERAGE_MULTI ||
                      errors[idx] != GRID_ERR_OVERLAP;
    coverage[idx] = ALEA_COVERAGE_MULTI;
    if (out_secondary_cell_ids) out_secondary_cell_ids[idx] = second_cell;
    errors[idx] = GRID_ERR_OVERLAP;
    return newly_multi ? 1 : 0;
}

static int update_tile_coverage_from_chain_candidates_bitset(
    alea_system_t* sys,
    const alea_slice_plane_t* plane,
    const alea_hier_spatial_chain_hit_t* hits,
    int hit_count,
    double u_min, double v_min,
    double du, double dv,
    int nu,
    int universe_depth,
    int* out_secondary_cell_ids,
    uint8_t* coverage,
    uint8_t* errors,
    int ti, int tj,
    int i_end, int j_end) {
    int tw = i_end - ti;
    int th = j_end - tj;
    int tile_pixels = tw * th;
    if (tile_pixels <= 0 || tile_pixels > 64) return -2;

    double gx[64], gy[64], gz[64];
    int pi_arr[64], pj_arr[64];
    int count[64] = {0};
    int second_cell[64];
    int max_depth[64];
    for (int p = 0; p < tile_pixels; p++) {
        second_cell[p] = -1;
        max_depth[p] = -1;
    }

    int p = 0;
    for (int j = tj; j < j_end; j++) {
        double v = v_min + (j + 0.5) * dv;
        for (int i = ti; i < i_end; i++, p++) {
            double u = u_min + (i + 0.5) * du;
            gx[p] = plane->origin[0] + u * plane->u_axis[0] + v * plane->v_axis[0];
            gy[p] = plane->origin[1] + u * plane->u_axis[1] + v * plane->v_axis[1];
            gz[p] = plane->origin[2] + u * plane->u_axis[2] + v * plane->v_axis[2];
            pi_arr[p] = i;
            pj_arr[p] = j;
        }
    }

    uint64_t all_mask = tile_pixels == 64 ? UINT64_MAX : ((1ULL << tile_pixels) - 1ULL);

    for (int h = 0; h < hit_count; h++) {
        const alea_hier_spatial_chain_hit_t* hit = &hits[h];
        const alea_spatial_hit_t* sh = &hit->hit;
        uint64_t active = all_mask;

        if (universe_depth >= 0) {
            if (sh->depth != universe_depth) continue;
            for (int q = 0; q < tile_pixels; q++) {
                if (count[q] > 1) active &= ~(1ULL << q);
            }
        } else {
            for (int q = 0; q < tile_pixels; q++) {
                if (sh->depth < max_depth[q] ||
                    (sh->depth == max_depth[q] && count[q] > 1)) {
                    active &= ~(1ULL << q);
                }
            }
        }
        if (!active) continue;

        uint64_t inside = 0;
        for (int q = 0; q < tile_pixels; q++) {
            if ((active & (1ULL << q)) == 0) continue;
            int rc = chain_hit_contains_point(sys, hit, gx[q], gy[q], gz[q]);
            if (rc < 0) return -1;
            if (rc) inside |= (1ULL << q);
        }
        if (!inside) continue;

        for (int q = 0; q < tile_pixels; q++) {
            if ((inside & (1ULL << q)) == 0) continue;
            if (universe_depth >= 0) {
                if (count[q] == 1) second_cell[q] = sh->cell_id;
                count[q]++;
            } else if (sh->depth > max_depth[q]) {
                max_depth[q] = sh->depth;
                count[q] = 1;
                second_cell[q] = -1;
            } else if (sh->depth == max_depth[q]) {
                if (count[q] == 1) second_cell[q] = sh->cell_id;
                count[q]++;
            }
        }
    }

    int updated = 0;
    for (int q = 0; q < tile_pixels; q++) {
        size_t idx = (size_t)pj_arr[q] * (size_t)nu + (size_t)pi_arr[q];
        if (count[q] <= 0) {
            coverage[idx] = ALEA_COVERAGE_NONE;
            if (out_secondary_cell_ids) out_secondary_cell_ids[idx] = -1;
            errors[idx] = GRID_ERR_UNDEFINED;
        } else if (count[q] == 1) {
            coverage[idx] = ALEA_COVERAGE_ONE;
            if (out_secondary_cell_ids) out_secondary_cell_ids[idx] = -1;
            if (errors[idx] == GRID_ERR_OVERLAP) errors[idx] = GRID_ERR_NONE;
        } else {
            if (coverage[idx] != ALEA_COVERAGE_MULTI ||
                errors[idx] != GRID_ERR_OVERLAP) {
                updated++;
            }
            coverage[idx] = ALEA_COVERAGE_MULTI;
            if (out_secondary_cell_ids) out_secondary_cell_ids[idx] = second_cell[q];
            errors[idx] = GRID_ERR_OVERLAP;
        }
    }

    return updated;
}

static int probe_curve_grid_point_for_exact_coverage(
    alea_system_t* sys,
    const alea_slice_plane_t* plane,
    double u_min, double v_min,
    double du, double dv,
    double inv_du, double inv_dv,
    int nu, int nv,
    int universe_depth,
    int* out_secondary_cell_ids,
    uint8_t* coverage,
    uint8_t* errors,
    uint8_t* probed,
    double u, double v) {
    int pi = (int)((u - u_min) * inv_du);
    int pj = (int)((v - v_min) * inv_dv);
    if (pi < 0 || pi >= nu || pj < 0 || pj >= nv) return 0;

    int updated = 0;
    for (int dj = -1; dj <= 1; dj++) {
        for (int di = -1; di <= 1; di++) {
            int qi = pi + di, qj = pj + dj;
            if (qi < 0 || qi >= nu || qj < 0 || qj >= nv) continue;

            size_t bit_idx = (size_t)qj * (size_t)nu + (size_t)qi;
            size_t byte_idx = bit_idx / 8;
            uint8_t bit_mask = (uint8_t)(1u << (bit_idx % 8));
            if (probed[byte_idx] & bit_mask) continue;
            probed[byte_idx] |= bit_mask;

            updated += update_pixel_coverage_exact(
                sys, plane, u_min, v_min, du, dv, nu, nv, universe_depth,
                out_secondary_cell_ids, coverage, errors, qi, qj);
        }
    }

    return updated;
}

static int refine_quartic_curve_overlaps(alea_system_t* sys,
                                         const alea_slice_view_t* view,
                                         const alea_curve_2d_t* curve,
                                         double sample_spacing,
                                         double u_min, double v_min,
                                         double du, double dv,
                                         double inv_du, double inv_dv,
                                         int nu, int nv,
                                         int universe_depth,
                                         const int* cell_ids,
                                         uint8_t* errors,
                                         uint8_t* probed) {
    double bbox_umin, bbox_umax, bbox_vmin, bbox_vmax;
    alea_curve_bbox(curve, &bbox_umin, &bbox_umax, &bbox_vmin, &bbox_vmax);

    double v_lo = fmax(bbox_vmin, view->v_min);
    double v_hi = fmin(bbox_vmax, view->v_max);
    if (v_lo >= v_hi) return 0;

    int n_scanlines = (int)((v_hi - v_lo) / sample_spacing);
    if (n_scanlines < 2) n_scanlines = 2;
    if (n_scanlines > 10000) n_scanlines = 10000;

    int new_overlaps = 0;
    for (int i = 0; i <= n_scanlines; i++) {
        double v_scan = v_lo + (v_hi - v_lo) * (double)i / (double)n_scanlines;
        double u_vals[16];
        int n_isect = alea_curve_scanline_intersect(curve, v_scan, u_vals, 16);

        for (int j = 0; j < n_isect; j++) {
            double u = u_vals[j];
            if (u < view->u_min || u > view->u_max) continue;
            new_overlaps += probe_curve_grid_point_for_overlap(
                sys, &view->plane, u_min, v_min, du, dv, inv_du, inv_dv,
                nu, nv, universe_depth, cell_ids, errors, probed, u, v_scan);
        }
    }

    return new_overlaps;
}

static int refine_polygon_curve_overlaps(alea_system_t* sys,
                                         const alea_slice_view_t* view,
                                         const alea_curve_2d_t* curve,
                                         double sample_spacing,
                                         double u_min, double v_min,
                                         double du, double dv,
                                         double inv_du, double inv_dv,
                                         int nu, int nv,
                                         int universe_depth,
                                         const int* cell_ids,
                                         uint8_t* errors,
                                         uint8_t* probed) {
    const alea_polygon_2d_t* poly = &curve->data.polygon;
    int nverts = poly->vertex_count;
    if (nverts < 2) return 0;

    int n_edges = poly->closed ? nverts : (nverts - 1);
    int new_overlaps = 0;

    for (int e = 0; e < n_edges; e++) {
        int i0 = e;
        int i1 = (e + 1) % nverts;
        double x0 = poly->vertices[i0][0], y0 = poly->vertices[i0][1];
        double x1 = poly->vertices[i1][0], y1 = poly->vertices[i1][1];
        double dx = x1 - x0;
        double dy = y1 - y0;
        double edge_len = sqrt(dx * dx + dy * dy);
        if (edge_len < 1e-15) continue;

        int n_samples = (int)(edge_len / sample_spacing);
        if (n_samples < 2) n_samples = 2;
        if (n_samples > 10000) n_samples = 10000;

        for (int s = 0; s <= n_samples; s++) {
            double frac = (double)s / (double)n_samples;
            double u = x0 + frac * dx;
            double v = y0 + frac * dy;
            if (u < view->u_min || u > view->u_max ||
                v < view->v_min || v > view->v_max) {
                continue;
            }
            new_overlaps += probe_curve_grid_point_for_overlap(
                sys, &view->plane, u_min, v_min, du, dv, inv_du, inv_dv,
                nu, nv, universe_depth, cell_ids, errors, probed, u, v);
        }
    }

    return new_overlaps;
}

int alea_check_grid_overlaps_curves(alea_system_t* sys,
                                    const alea_slice_view_t* view,
                                    const alea_slice_curves_t* curves,
                                    int nu, int nv,
                                    int universe_depth,
                                    const int* cell_ids,
                                    uint8_t* errors) {
    if (!sys || !view || !curves || !cell_ids || !errors || nu <= 0 || nv <= 0)
        return -1;

    const alea_curve_collection_t* coll = &curves->internal;
    if (coll->curves.count == 0) return 0;

    const alea_slice_plane_t* plane = &view->plane;
    double u_min = view->u_min, u_max = view->u_max;
    double v_min = view->v_min, v_max = view->v_max;
    double du = (u_max - u_min) / nu;
    double dv = (v_max - v_min) / nv;
    double inv_du = 1.0 / du;
    double inv_dv = 1.0 / dv;

    /* Bitmap to avoid re-probing the same pixel from multiple curves */
    size_t bitmap_size = ((size_t)nu * nv + 7) / 8;
    uint8_t* probed = calloc(bitmap_size, 1);
    if (!probed) return -1;

    double sample_spacing = fmin(du, dv) * 0.7;  /* sub-pixel stepping */

    int new_overlaps = 0;

    for (size_t ci = 0; ci < coll->curves.count; ci++) {
        const alea_curve_2d_t* curve = &coll->curves.data[ci];

        if (curve->type == ALEA_CURVE_NONE || curve->type == ALEA_CURVE_POINT)
            continue;

        if (curve->type == ALEA_CURVE_QUARTIC) {
            new_overlaps += refine_quartic_curve_overlaps(
                sys, view, curve, sample_spacing,
                u_min, v_min, du, dv, inv_du, inv_dv,
                nu, nv, universe_depth, cell_ids, errors, probed);
            continue;
        }

        if (curve->type == ALEA_CURVE_POLYGON) {
            new_overlaps += refine_polygon_curve_overlaps(
                sys, view, curve, sample_spacing,
                u_min, v_min, du, dv, inv_du, inv_dv,
                nu, nv, universe_depth, cell_ids, errors, probed);
            continue;
        }

        double t_lo, t_hi;
        get_clipped_param_range(curve, view, &t_lo, &t_hi);
        if (t_lo >= t_hi) continue;

        double arc_approx = t_hi - t_lo;
        if (curve->type == ALEA_CURVE_CIRCLE || curve->type == ALEA_CURVE_ARC) {
            arc_approx = (t_hi - t_lo) * curve->data.circle.radius;
        } else if (curve->type == ALEA_CURVE_ELLIPSE || curve->type == ALEA_CURVE_ELLIPSE_ARC) {
            double avg_r = (curve->data.ellipse.semi_a + curve->data.ellipse.semi_b) * 0.5;
            arc_approx = (t_hi - t_lo) * avg_r;
        }

        int n_samples = (int)(arc_approx / sample_spacing);
        if (n_samples < 2) n_samples = 2;
        if (n_samples > 10000) n_samples = 10000;

        for (int s = 0; s <= n_samples; s++) {
            double t = t_lo + (t_hi - t_lo) * s / n_samples;

            double u, v;
            if (!alea_curve_eval(curve, t, &u, &v)) continue;
            if (u < u_min || u > u_max || v < v_min || v > v_max) continue;

            new_overlaps += probe_curve_grid_point_for_overlap(
                sys, plane, u_min, v_min, du, dv, inv_du, inv_dv,
                nu, nv, universe_depth, cell_ids, errors, probed, u, v);
        }
    }

    free(probed);
    return new_overlaps;
}

int alea_check_grid_overlaps_curves_coverage(
    const alea_slice_view_t* view,
    const alea_slice_curves_t* curves,
    int nu, int nv,
    const uint8_t* coverage,
    uint8_t* errors) {
    if (!view || !curves || !coverage || !errors || nu <= 0 || nv <= 0)
        return -1;

    const alea_curve_collection_t* coll = &curves->internal;
    if (coll->curves.count == 0) return 0;

    double u_min = view->u_min, u_max = view->u_max;
    double v_min = view->v_min, v_max = view->v_max;
    double du = (u_max - u_min) / nu;
    double dv = (v_max - v_min) / nv;
    double inv_du = 1.0 / du;
    double inv_dv = 1.0 / dv;

    size_t bitmap_size = ((size_t)nu * nv + 7) / 8;
    uint8_t* probed = calloc(bitmap_size, 1);
    if (!probed) return -1;

    double sample_spacing = fmin(du, dv) * 0.7;
    int new_overlaps = 0;

    for (size_t ci = 0; ci < coll->curves.count; ci++) {
        const alea_curve_2d_t* curve = &coll->curves.data[ci];

        if (curve->type == ALEA_CURVE_NONE || curve->type == ALEA_CURVE_POINT)
            continue;

        if (curve->type == ALEA_CURVE_QUARTIC) {
            double bbox_umin, bbox_umax, bbox_vmin, bbox_vmax;
            alea_curve_bbox(curve, &bbox_umin, &bbox_umax, &bbox_vmin, &bbox_vmax);
            double v_lo = fmax(bbox_vmin, view->v_min);
            double v_hi = fmin(bbox_vmax, view->v_max);
            if (v_lo >= v_hi) continue;

            int n_scanlines = (int)((v_hi - v_lo) / sample_spacing);
            if (n_scanlines < 2) n_scanlines = 2;
            if (n_scanlines > 10000) n_scanlines = 10000;

            for (int i = 0; i <= n_scanlines; i++) {
                double v_scan = v_lo + (v_hi - v_lo) * (double)i / (double)n_scanlines;
                double u_vals[16];
                int n_isect = alea_curve_scanline_intersect(curve, v_scan, u_vals, 16);
                for (int j = 0; j < n_isect; j++) {
                    double u = u_vals[j];
                    if (u < view->u_min || u > view->u_max) continue;
                    new_overlaps += probe_curve_grid_point_from_coverage(
                        u_min, v_min, inv_du, inv_dv,
                        nu, nv, coverage, errors, probed, u, v_scan);
                }
            }
            continue;
        }

        if (curve->type == ALEA_CURVE_POLYGON) {
            const alea_polygon_2d_t* poly = &curve->data.polygon;
            int nverts = poly->vertex_count;
            if (nverts < 2) continue;
            int n_edges = poly->closed ? nverts : (nverts - 1);

            for (int e = 0; e < n_edges; e++) {
                int i0 = e;
                int i1 = (e + 1) % nverts;
                double x0 = poly->vertices[i0][0], y0 = poly->vertices[i0][1];
                double x1 = poly->vertices[i1][0], y1 = poly->vertices[i1][1];
                double dx = x1 - x0;
                double dy = y1 - y0;
                double edge_len = sqrt(dx * dx + dy * dy);
                if (edge_len < 1e-15) continue;

                int n_samples = (int)(edge_len / sample_spacing);
                if (n_samples < 2) n_samples = 2;
                if (n_samples > 10000) n_samples = 10000;

                for (int s = 0; s <= n_samples; s++) {
                    double frac = (double)s / (double)n_samples;
                    double u = x0 + frac * dx;
                    double v = y0 + frac * dy;
                    if (u < view->u_min || u > view->u_max ||
                        v < view->v_min || v > view->v_max)
                        continue;
                    new_overlaps += probe_curve_grid_point_from_coverage(
                        u_min, v_min, inv_du, inv_dv,
                        nu, nv, coverage, errors, probed, u, v);
                }
            }
            continue;
        }

        double t_lo, t_hi;
        get_clipped_param_range(curve, view, &t_lo, &t_hi);
        if (t_lo >= t_hi) continue;

        double arc_approx = t_hi - t_lo;
        if (curve->type == ALEA_CURVE_CIRCLE || curve->type == ALEA_CURVE_ARC) {
            arc_approx = (t_hi - t_lo) * curve->data.circle.radius;
        } else if (curve->type == ALEA_CURVE_ELLIPSE ||
                   curve->type == ALEA_CURVE_ELLIPSE_ARC) {
            double avg_r = (curve->data.ellipse.semi_a + curve->data.ellipse.semi_b) * 0.5;
            arc_approx = (t_hi - t_lo) * avg_r;
        }

        int n_samples = (int)(arc_approx / sample_spacing);
        if (n_samples < 2) n_samples = 2;
        if (n_samples > 10000) n_samples = 10000;

        for (int s = 0; s <= n_samples; s++) {
            double t = t_lo + (t_hi - t_lo) * s / n_samples;
            double u, v;
            if (!alea_curve_eval(curve, t, &u, &v)) continue;
            if (u < u_min || u > u_max || v < v_min || v > v_max) continue;

            new_overlaps += probe_curve_grid_point_from_coverage(
                u_min, v_min, inv_du, inv_dv,
                nu, nv, coverage, errors, probed, u, v);
        }
    }

    free(probed);
    return new_overlaps;
}

int alea_refine_grid_coverage_curves_exact(
    alea_system_t* sys,
    const alea_slice_view_t* view,
    const alea_slice_curves_t* curves,
    int nu, int nv,
    int universe_depth,
    int* out_secondary_cell_ids,
    uint8_t* coverage,
    uint8_t* errors) {
    if (!sys || !view || !curves || !coverage || !errors || nu <= 0 || nv <= 0)
        return -1;
    alea_point_coverage_stats_reset();

    const alea_curve_collection_t* coll = &curves->internal;
    if (coll->curves.count == 0) return 0;

    const alea_slice_plane_t* plane = &view->plane;
    double u_min = view->u_min, u_max = view->u_max;
    double v_min = view->v_min, v_max = view->v_max;
    double du = (u_max - u_min) / nu;
    double dv = (v_max - v_min) / nv;
    double inv_du = 1.0 / du;
    double inv_dv = 1.0 / dv;

    size_t bitmap_size = ((size_t)nu * nv + 7) / 8;
    uint8_t* probed = calloc(bitmap_size, 1);
    if (!probed) return -1;

    double sample_spacing = fmin(du, dv) * 0.7;
    int updated = 0;

    for (size_t ci = 0; ci < coll->curves.count; ci++) {
        const alea_curve_2d_t* curve = &coll->curves.data[ci];
        if (curve->type == ALEA_CURVE_NONE || curve->type == ALEA_CURVE_POINT)
            continue;

        if (curve->type == ALEA_CURVE_QUARTIC) {
            double bbox_umin, bbox_umax, bbox_vmin, bbox_vmax;
            alea_curve_bbox(curve, &bbox_umin, &bbox_umax, &bbox_vmin, &bbox_vmax);
            double v_lo = fmax(bbox_vmin, view->v_min);
            double v_hi = fmin(bbox_vmax, view->v_max);
            if (v_lo >= v_hi) continue;

            int n_scanlines = (int)((v_hi - v_lo) / sample_spacing);
            if (n_scanlines < 2) n_scanlines = 2;
            if (n_scanlines > 10000) n_scanlines = 10000;

            for (int i = 0; i <= n_scanlines; i++) {
                double v_scan = v_lo + (v_hi - v_lo) * (double)i / (double)n_scanlines;
                double u_vals[16];
                int n_isect = alea_curve_scanline_intersect(curve, v_scan, u_vals, 16);
                for (int j = 0; j < n_isect; j++) {
                    double u = u_vals[j];
                    if (u < view->u_min || u > view->u_max) continue;
                    updated += probe_curve_grid_point_for_exact_coverage(
                        sys, plane, u_min, v_min, du, dv, inv_du, inv_dv,
                        nu, nv, universe_depth, out_secondary_cell_ids,
                        coverage, errors, probed, u, v_scan);
                }
            }
            continue;
        }

        if (curve->type == ALEA_CURVE_POLYGON) {
            const alea_polygon_2d_t* poly = &curve->data.polygon;
            int nverts = poly->vertex_count;
            if (nverts < 2) continue;
            int n_edges = poly->closed ? nverts : (nverts - 1);
            for (int e = 0; e < n_edges; e++) {
                int i0 = e;
                int i1 = (e + 1) % nverts;
                double x0 = poly->vertices[i0][0], y0 = poly->vertices[i0][1];
                double x1 = poly->vertices[i1][0], y1 = poly->vertices[i1][1];
                double dx = x1 - x0;
                double dy = y1 - y0;
                double edge_len = sqrt(dx * dx + dy * dy);
                if (edge_len < 1e-15) continue;

                int n_samples = (int)(edge_len / sample_spacing);
                if (n_samples < 2) n_samples = 2;
                if (n_samples > 10000) n_samples = 10000;

                for (int s = 0; s <= n_samples; s++) {
                    double frac = (double)s / (double)n_samples;
                    double u = x0 + frac * dx;
                    double v = y0 + frac * dy;
                    if (u < view->u_min || u > view->u_max ||
                        v < view->v_min || v > view->v_max)
                        continue;
                    updated += probe_curve_grid_point_for_exact_coverage(
                        sys, plane, u_min, v_min, du, dv, inv_du, inv_dv,
                        nu, nv, universe_depth, out_secondary_cell_ids,
                        coverage, errors, probed, u, v);
                }
            }
            continue;
        }

        double t_lo, t_hi;
        get_clipped_param_range(curve, view, &t_lo, &t_hi);
        if (t_lo >= t_hi) continue;

        double arc_approx = t_hi - t_lo;
        if (curve->type == ALEA_CURVE_CIRCLE || curve->type == ALEA_CURVE_ARC) {
            arc_approx = (t_hi - t_lo) * curve->data.circle.radius;
        } else if (curve->type == ALEA_CURVE_ELLIPSE ||
                   curve->type == ALEA_CURVE_ELLIPSE_ARC) {
            double avg_r = (curve->data.ellipse.semi_a + curve->data.ellipse.semi_b) * 0.5;
            arc_approx = (t_hi - t_lo) * avg_r;
        }

        int n_samples = (int)(arc_approx / sample_spacing);
        if (n_samples < 2) n_samples = 2;
        if (n_samples > 10000) n_samples = 10000;

        for (int s = 0; s <= n_samples; s++) {
            double t = t_lo + (t_hi - t_lo) * s / n_samples;
            double u, v;
            if (!alea_curve_eval(curve, t, &u, &v)) continue;
            if (u < u_min || u > u_max || v < v_min || v > v_max) continue;

            updated += probe_curve_grid_point_for_exact_coverage(
                sys, plane, u_min, v_min, du, dv, inv_du, inv_dv,
                nu, nv, universe_depth, out_secondary_cell_ids,
                coverage, errors, probed, u, v);
        }
    }

    free(probed);
    return updated;
}

int alea_refine_grid_coverage_tiles_exact(
    alea_system_t* sys,
    const alea_slice_view_t* view,
    int nu, int nv,
    int universe_depth,
    int tile_w, int tile_h,
    int* out_secondary_cell_ids,
    uint8_t* coverage,
    uint8_t* errors) {
    if (!sys || !view || !coverage || !errors || nu <= 0 || nv <= 0)
        return -1;
    if (tile_w <= 0) tile_w = 16;
    if (tile_h <= 0) tile_h = 16;
    alea_tile_coverage_stats_reset();
    alea_point_coverage_stats_reset();

    if (alea_prepare_query_acceleration(sys) != 0)
        return -1;

    size_t max_hits = 4096;
    const char* max_hits_env = getenv("ALEA_TILE_MAX_CANDIDATES");
    if (max_hits_env && max_hits_env[0]) {
        unsigned long value = strtoul(max_hits_env, NULL, 10);
        if (value > 0) max_hits = (size_t)value;
    }
    if (max_hits < 4096) max_hits = 4096;
    alea_spatial_hit_t* hits = malloc(max_hits * sizeof(*hits));
    if (!hits) return -1;
    alea_hier_spatial_chain_hit_t* chain_hits = NULL;

    double u_range = view->u_max - view->u_min;
    double v_range = view->v_max - view->v_min;
    double du = u_range / nu;
    double dv = v_range / nv;
    double eps = fmax(fabs(du), fabs(dv)) * 2.0 + 1e-10;
    bool use_hier = alea_spatial_mode_is_hierarchical(sys);
    bool use_hier_direct_region =
        use_hier && getenv("ALEA_TILE_DIRECT_REGION") != NULL;
    bool use_hier_chain_candidates =
        use_hier && getenv("ALEA_TILE_CHAIN_CANDIDATES") != NULL;
    bool use_chain_bitset =
        use_hier_chain_candidates && getenv("ALEA_TILE_CHAIN_BITSET") != NULL;
    if (use_hier_chain_candidates) {
        chain_hits = malloc(max_hits * sizeof(*chain_hits));
        if (!chain_hits) {
            free(hits);
            return -1;
        }
    }
    size_t cap_storm_threshold = 16;
    const char* cap_storm_env = getenv("ALEA_TILE_CAP_STORM_THRESHOLD");
    if (cap_storm_env && cap_storm_env[0]) {
        cap_storm_threshold = (size_t)strtoul(cap_storm_env, NULL, 10);
    }
    size_t eval_candidate_limit = 4096;
    const char* eval_limit_env = getenv("ALEA_TILE_EVAL_MAX_CANDIDATES");
    if (eval_limit_env && eval_limit_env[0]) {
        unsigned long value = strtoul(eval_limit_env, NULL, 10);
        if (value > 0) eval_candidate_limit = (size_t)value;
    }
    size_t capped_tile_queries = 0;
    size_t evaluable_tile_queries = 0;
    bool exact_remaining_tiles = false;
    int updated = 0;

    for (int tj = 0; tj < nv; tj += tile_h) {
        int j_end = tj + tile_h;
        if (j_end > nv) j_end = nv;
        double v0 = view->v_min + tj * dv;
        double v1 = view->v_min + j_end * dv;

        for (int ti = 0; ti < nu; ti += tile_w) {
            int i_end = ti + tile_w;
            if (i_end > nu) i_end = nu;
            double u0 = view->u_min + ti * du;
            double u1 = view->u_min + i_end * du;
            size_t tile_pixels = (size_t)(i_end - ti) * (size_t)(j_end - tj);
            g_tile_coverage_stats.tiles++;

            if (exact_remaining_tiles) {
                g_tile_coverage_stats.fallback_tiles++;
                g_tile_coverage_stats.exact_fallback_pixels += tile_pixels;
                for (int j = tj; j < j_end; j++) {
                    for (int i = ti; i < i_end; i++) {
                        updated += update_pixel_coverage_exact(
                            sys, &view->plane, view->u_min, view->v_min,
                            du, dv, nu, nv, universe_depth,
                            out_secondary_cell_ids, coverage, errors, i, j);
                    }
                }
                continue;
            }

            alea_bbox_t query = tile_query_bbox(view, u0, u1, v0, v1, eps);
            int hit_count;
            if (use_hier_chain_candidates) {
                hit_count = alea_hier_spatial_query_region_chain(
                    sys, &query, chain_hits, max_hits);
            } else if (use_hier_direct_region) {
                hit_count = alea_hier_spatial_query_region_direct(
                    sys, &query, hits, max_hits);
            } else if (use_hier) {
                hit_count = alea_hier_spatial_query_region(
                    sys, &query, hits, max_hits);
            } else {
                hit_count = alea_spatial_query_region(
                    sys, &query, hits, max_hits);
            }
            if (hit_count < 0) {
                g_tile_coverage_stats.query_errors++;
                free(chain_hits);
                free(hits);
                return -1;
            }

            g_tile_coverage_stats.candidate_total += (size_t)hit_count;
            if ((size_t)hit_count > g_tile_coverage_stats.candidate_max)
                g_tile_coverage_stats.candidate_max = (size_t)hit_count;

            if ((size_t)hit_count >= max_hits) {
                capped_tile_queries++;
                g_tile_coverage_stats.fallback_tiles++;
                g_tile_coverage_stats.exact_fallback_pixels += tile_pixels;
                for (int j = tj; j < j_end; j++) {
                    for (int i = ti; i < i_end; i++) {
                        updated += update_pixel_coverage_exact(
                            sys, &view->plane, view->u_min, view->v_min,
                            du, dv, nu, nv, universe_depth,
                            out_secondary_cell_ids, coverage, errors, i, j);
                    }
                }
                if (cap_storm_threshold > 0 &&
                    evaluable_tile_queries == 0 &&
                    capped_tile_queries >= cap_storm_threshold) {
                    exact_remaining_tiles = true;
                }
                continue;
            }

            if (use_hier_chain_candidates &&
                eval_candidate_limit > 0 &&
                (size_t)hit_count > eval_candidate_limit) {
                g_tile_coverage_stats.fallback_tiles++;
                g_tile_coverage_stats.exact_fallback_pixels += tile_pixels;
                for (int j = tj; j < j_end; j++) {
                    for (int i = ti; i < i_end; i++) {
                        updated += update_pixel_coverage_exact(
                            sys, &view->plane, view->u_min, view->v_min,
                            du, dv, nu, nv, universe_depth,
                            out_secondary_cell_ids, coverage, errors, i, j);
                    }
                }
                continue;
            }

            evaluable_tile_queries++;
            if (!use_hier_chain_candidates)
                hit_count = dedup_spatial_hits(hits, hit_count);
            g_tile_coverage_stats.dedup_candidate_total += (size_t)hit_count;
            if ((size_t)hit_count > g_tile_coverage_stats.dedup_candidate_max)
                g_tile_coverage_stats.dedup_candidate_max = (size_t)hit_count;
            g_tile_coverage_stats.pixels += tile_pixels;
            g_tile_coverage_stats.candidate_pixel_tests +=
                (size_t)hit_count * tile_pixels;

            if (use_chain_bitset && tile_pixels <= 64) {
                int rc = update_tile_coverage_from_chain_candidates_bitset(
                    sys, &view->plane, chain_hits, hit_count,
                    view->u_min, view->v_min, du, dv, nu,
                    universe_depth, out_secondary_cell_ids,
                    coverage, errors, ti, tj, i_end, j_end);
                if (rc >= 0) {
                    updated += rc;
                    continue;
                }
                g_tile_coverage_stats.fallback_tiles++;
                g_tile_coverage_stats.exact_fallback_pixels += tile_pixels;
                for (int j = tj; j < j_end; j++) {
                    for (int i = ti; i < i_end; i++) {
                        updated += update_pixel_coverage_exact(
                            sys, &view->plane, view->u_min, view->v_min,
                            du, dv, nu, nv, universe_depth,
                            out_secondary_cell_ids, coverage, errors, i, j);
                    }
                }
                continue;
            }

            for (int j = tj; j < j_end; j++) {
                for (int i = ti; i < i_end; i++) {
                    if (use_hier_chain_candidates) {
                        int rc = update_pixel_coverage_from_chain_candidates(
                            sys, &view->plane, chain_hits, hit_count,
                            view->u_min, view->v_min, du, dv, nu,
                            universe_depth, out_secondary_cell_ids,
                            coverage, errors, i, j);
                        if (rc < 0) {
                            g_tile_coverage_stats.fallback_tiles++;
                            g_tile_coverage_stats.exact_fallback_pixels++;
                            rc = update_pixel_coverage_exact(
                                sys, &view->plane, view->u_min, view->v_min,
                                du, dv, nu, nv, universe_depth,
                                out_secondary_cell_ids, coverage, errors, i, j);
                        }
                        updated += rc;
                    } else {
                        updated += update_pixel_coverage_from_candidates(
                            sys, &view->plane, hits, hit_count,
                            view->u_min, view->v_min, du, dv, nu,
                            universe_depth, out_secondary_cell_ids,
                            coverage, errors, i, j);
                    }
                }
            }
        }
    }

    free(chain_hits);
    free(hits);
    g_tile_coverage_stats.refined_pixels = (size_t)updated;
    return updated;
}

int alea_refine_grid_coverage_paths_exact(
    alea_system_t* sys,
    const alea_slice_view_t* view,
    int nu, int nv,
    int universe_depth,
    int tile_w, int tile_h,
    const int* primary_cell_ids,
    const uint32_t* path_ids,
    const alea_slice_path_table_t* paths,
    int* out_secondary_cell_ids,
    uint8_t* coverage,
    uint8_t* errors) {
    if (!sys || !view || !primary_cell_ids || !path_ids || !paths ||
        !coverage || !errors || nu <= 0 || nv <= 0)
        return -1;
    if (tile_w <= 0) tile_w = 16;
    if (tile_h <= 0) tile_h = 16;

    alea_tile_coverage_stats_reset();
    alea_point_coverage_stats_reset();

    if (!alea_spatial_mode_is_hierarchical(sys))
        return -1;
    if (alea_prepare_query_acceleration(sys) != 0)
        return -1;

    size_t max_hits = 4096;
    const char* max_hits_env = getenv("ALEA_PATH_MAX_CANDIDATES");
    if (!max_hits_env || !max_hits_env[0])
        max_hits_env = getenv("ALEA_TILE_MAX_CANDIDATES");
    if (max_hits_env && max_hits_env[0]) {
        unsigned long value = strtoul(max_hits_env, NULL, 10);
        if (value > 0) max_hits = (size_t)value;
    }
    if (max_hits < 128) max_hits = 128;

    alea_spatial_hit_t* hits = malloc(max_hits * sizeof(*hits));
    if (!hits) return -1;
    size_t max_group_pixels = (size_t)tile_w * (size_t)tile_h;
    size_t* group_indices = malloc(max_group_pixels * sizeof(*group_indices));
    double* group_x = malloc(max_group_pixels * sizeof(*group_x));
    double* group_y = malloc(max_group_pixels * sizeof(*group_y));
    double* group_z = malloc(max_group_pixels * sizeof(*group_z));
    int* group_count = malloc(max_group_pixels * sizeof(*group_count));
    int* group_second = malloc(max_group_pixels * sizeof(*group_second));
    if (!group_indices || !group_x || !group_y || !group_z ||
        !group_count || !group_second) {
        free(group_second);
        free(group_count);
        free(group_z);
        free(group_y);
        free(group_x);
        free(group_indices);
        free(hits);
        return -1;
    }
    bool use_path_2d_index = getenv("ALEA_PATH_2D_INDEX") != NULL;
    bool verify_path_2d_index =
        use_path_2d_index && getenv("ALEA_PATH_2D_VERIFY") != NULL;
    size_t verify_2d_log_limit = 20;
    const char* verify_limit_env = getenv("ALEA_PATH_2D_VERIFY_LIMIT");
    if (verify_limit_env && verify_limit_env[0]) {
        unsigned long value = strtoul(verify_limit_env, NULL, 10);
        verify_2d_log_limit = (size_t)value;
    }
    size_t verify_2d_log_count = 0;
    path_tile_bucket_t verify_bucket = {0};
    path_slice_index_t* path_indexes = NULL;
    bool* path_index_ready = NULL;
    bool* path_index_disabled = NULL;
    size_t path_bucket_limit = max_hits;
    const char* bucket_limit_env = getenv("ALEA_PATH_2D_BUCKET_LIMIT");
    if (bucket_limit_env && bucket_limit_env[0]) {
        unsigned long value = strtoul(bucket_limit_env, NULL, 10);
        if (value > 0) path_bucket_limit = (size_t)value;
    }
    int path_tile_pad = 1;
    const char* tile_pad_env = getenv("ALEA_PATH_2D_TILE_PAD");
    if (tile_pad_env && tile_pad_env[0]) {
        int value = atoi(tile_pad_env);
        if (value >= 0) path_tile_pad = value;
    }
    if (use_path_2d_index && paths->count > 0) {
        path_indexes = calloc(paths->count, sizeof(*path_indexes));
        path_index_ready = calloc(paths->count, sizeof(*path_index_ready));
        path_index_disabled = calloc(paths->count, sizeof(*path_index_disabled));
        if (!path_indexes || !path_index_ready || !path_index_disabled) {
            free(path_index_disabled);
            free(path_index_ready);
            free(path_indexes);
            free(group_second);
            free(group_count);
            free(group_z);
            free(group_y);
            free(group_x);
            free(group_indices);
            free(hits);
            return -1;
        }
    }

    double u_range = view->u_max - view->u_min;
    double v_range = view->v_max - view->v_min;
    double du = u_range / nu;
    double dv = v_range / nv;
    double eps = 1e-9;
    int updated = 0;

#ifdef _OPENMP
    if (!use_path_2d_index && getenv("ALEA_PATH_SERIAL") == NULL) {
        _Atomic int abort_flag = 0;
        int parallel_updated = 0;
        alea_tile_coverage_stats_t parallel_stats;
        memset(&parallel_stats, 0, sizeof(parallel_stats));

        #pragma omp parallel
        {
            alea_spatial_hit_t* local_hits =
                malloc(max_hits * sizeof(*local_hits));
            size_t* local_indices =
                malloc(max_group_pixels * sizeof(*local_indices));
            double* local_x = malloc(max_group_pixels * sizeof(*local_x));
            double* local_y = malloc(max_group_pixels * sizeof(*local_y));
            double* local_z = malloc(max_group_pixels * sizeof(*local_z));
            int* local_count = malloc(max_group_pixels * sizeof(*local_count));
            int* local_second = malloc(max_group_pixels * sizeof(*local_second));
            alea_tile_coverage_stats_t local_stats;
            memset(&local_stats, 0, sizeof(local_stats));
            int local_updated = 0;

            if (!local_hits || !local_indices || !local_x || !local_y ||
                !local_z || !local_count || !local_second) {
                atomic_store(&abort_flag, 1);
            }

            #pragma omp for collapse(2) schedule(dynamic, 1)
            for (int tj = 0; tj < nv; tj += tile_h) {
                for (int ti = 0; ti < nu; ti += tile_w) {
                    if (atomic_load(&abort_flag))
                        continue;
                    int rc = refine_path_tile_3d_exact(
                        sys, view, nu, nv, universe_depth,
                        tile_w, tile_h, ti, tj, du, dv, max_hits,
                        local_hits, max_group_pixels, local_indices,
                        local_x, local_y, local_z, local_count,
                        local_second, primary_cell_ids, path_ids, paths,
                        out_secondary_cell_ids, coverage, errors,
                        &local_stats);
                    if (rc < 0) {
                        atomic_store(&abort_flag, 1);
                    } else {
                        local_updated += rc;
                    }
                }
            }

            #pragma omp critical(path_refine_stats_merge)
            {
                parallel_updated += local_updated;
                tile_coverage_stats_add(&parallel_stats, &local_stats);
            }

            free(local_second);
            free(local_count);
            free(local_z);
            free(local_y);
            free(local_x);
            free(local_indices);
            free(local_hits);
        }

        free(path_index_disabled);
        free(path_index_ready);
        free(path_indexes);
        free(group_second);
        free(group_count);
        free(group_z);
        free(group_y);
        free(group_x);
        free(group_indices);
        path_tile_bucket_free(&verify_bucket);
        free(hits);

        if (atomic_load(&abort_flag))
            return -1;

        g_tile_coverage_stats = parallel_stats;
        g_tile_coverage_stats.refined_pixels = (size_t)parallel_updated;
        return parallel_updated;
    }
#endif

    for (int tj = 0; tj < nv; tj += tile_h) {
        int j_end = tj + tile_h;
        if (j_end > nv) j_end = nv;

        for (int ti = 0; ti < nu; ti += tile_w) {
            int i_end = ti + tile_w;
            if (i_end > nu) i_end = nu;
            g_tile_coverage_stats.tiles++;

            for (int sj = tj; sj < j_end; sj++) {
                for (int si = ti; si < i_end; si++) {
                    size_t seed_idx = (size_t)sj * (size_t)nu + (size_t)si;
                    uint32_t path_id = path_ids[seed_idx];
                    if (path_id == UINT32_MAX) {
                        g_tile_coverage_stats.fallback_tiles++;
                        g_tile_coverage_stats.exact_fallback_pixels++;
                        updated += update_pixel_coverage_exact(
                            sys, &view->plane, view->u_min, view->v_min,
                            du, dv, nu, nv, universe_depth,
                            out_secondary_cell_ids, coverage, errors, si, sj);
                        continue;
                    }
                    if ((size_t)path_id >= paths->count) {
                        path_tile_bucket_free(&verify_bucket);
                        free(hits);
                        return -1;
                    }

                    int already_done = 0;
                    for (int pj = tj; pj <= sj && !already_done; pj++) {
                        int pi_limit = (pj == sj) ? si : i_end;
                        for (int pi = ti; pi < pi_limit; pi++) {
                            size_t pidx = (size_t)pj * (size_t)nu + (size_t)pi;
                            if (path_ids[pidx] == path_id) {
                                already_done = 1;
                                break;
                            }
                        }
                    }
                    if (already_done) continue;

                    const alea_slice_path_record_t* path =
                        &paths->records[path_id];
                    if (universe_depth >= 0 && path->depth != universe_depth) {
                        g_tile_coverage_stats.fallback_tiles++;
                        updated += update_pixel_coverage_exact(
                            sys, &view->plane, view->u_min, view->v_min,
                            du, dv, nu, nv, universe_depth,
                            out_secondary_cell_ids, coverage, errors, si, sj);
                        continue;
                    }

                    alea_bbox_t local_bbox = {
                        .min_x = DBL_MAX, .max_x = -DBL_MAX,
                        .min_y = DBL_MAX, .max_y = -DBL_MAX,
                        .min_z = DBL_MAX, .max_z = -DBL_MAX
                    };
                    size_t group_pixels = 0;
                    for (int pj = tj; pj < j_end; pj++) {
                        for (int pi = ti; pi < i_end; pi++) {
                            size_t pidx = (size_t)pj * (size_t)nu + (size_t)pi;
                            if (path_ids[pidx] != path_id) continue;
                            if (group_pixels >= max_group_pixels) {
                                path_tile_bucket_free(&verify_bucket);
                                free(group_second);
                                free(group_count);
                                free(group_z);
                                free(group_y);
                                free(group_x);
                                free(group_indices);
                                free(hits);
                                return -1;
                            }
                            double lx, ly, lz;
                            slice_grid_world_point(&view->plane,
                                                   view->u_min, view->v_min,
                                                   du, dv, pi, pj,
                                                   &lx, &ly, &lz);
                            path_world_to_local_point(path->world_to_local,
                                                      &lx, &ly, &lz);
                            if (lx < local_bbox.min_x) local_bbox.min_x = lx;
                            if (lx > local_bbox.max_x) local_bbox.max_x = lx;
                            if (ly < local_bbox.min_y) local_bbox.min_y = ly;
                            if (ly > local_bbox.max_y) local_bbox.max_y = ly;
                            if (lz < local_bbox.min_z) local_bbox.min_z = lz;
                            if (lz > local_bbox.max_z) local_bbox.max_z = lz;
                            group_indices[group_pixels] = pidx;
                            group_x[group_pixels] = lx;
                            group_y[group_pixels] = ly;
                            group_z[group_pixels] = lz;
                            group_pixels++;
                        }
                    }
                    if (group_pixels == 0) continue;
                    local_bbox.min_x -= eps; local_bbox.max_x += eps;
                    local_bbox.min_y -= eps; local_bbox.max_y += eps;
                    local_bbox.min_z -= eps; local_bbox.max_z += eps;

                    alea_spatial_hit_t* eval_hits = hits;
                    int hit_count = -1;
                    bool using_2d_bucket = false;
                    bool candidate_query_saturated = false;
                    if (use_path_2d_index && path_indexes &&
                        !path_index_disabled[path_id]) {
                        if (!path_index_ready[path_id]) {
                            int brc = path_slice_index_build(
                                sys, view, nu, nv, tile_w, tile_h,
                                path_id, path, path_bucket_limit,
                                path_tile_pad,
                                &path_indexes[path_id]);
                            if (brc == 0) {
                                path_index_ready[path_id] = true;
                            } else {
                                path_index_disabled[path_id] = true;
                            }
                        }
                        if (path_index_ready[path_id]) {
                            path_slice_index_t* pidx = &path_indexes[path_id];
                            int tile_i = ti / tile_w;
                            int tile_j = tj / tile_h;
                            size_t bidx = (size_t)tile_j * (size_t)pidx->tiles_x +
                                          (size_t)tile_i;
                            path_tile_bucket_t* bucket = &pidx->buckets[bidx];
                            if (bucket->count <= (size_t)INT_MAX) {
                                eval_hits = bucket->hits;
                                hit_count = (int)bucket->count;
                                using_2d_bucket = true;
                            } else {
                                path_index_disabled[path_id] = true;
                            }
                        }
                    }
                    if (using_2d_bucket && verify_path_2d_index) {
                        verify_bucket.count = 0;
                        for (int h = 0; h < hit_count; h++) {
                            if (path_tile_bucket_append(&verify_bucket,
                                                        &eval_hits[h]) != 0) {
                                path_tile_bucket_free(&verify_bucket);
                                if (path_indexes) {
                                    for (size_t pi = 0; pi < paths->count; pi++)
                                        path_slice_index_free(&path_indexes[pi]);
                                }
                                free(path_index_disabled);
                                free(path_index_ready);
                                free(path_indexes);
                                free(group_second);
                                free(group_count);
                                free(group_z);
                                free(group_y);
                                free(group_x);
                                free(group_indices);
                                free(hits);
                                return -1;
                            }
                        }

                        int exact_hit_count =
                            alea_hier_spatial_query_universe_region(
                                sys, path->universe_id, &local_bbox,
                                hits, max_hits);
                        g_tile_coverage_stats.path_2d_verify_queries++;
                        if (exact_hit_count < 0) {
                            g_tile_coverage_stats.query_errors++;
                            path_tile_bucket_free(&verify_bucket);
                            if (path_indexes) {
                                for (size_t pi = 0; pi < paths->count; pi++)
                                    path_slice_index_free(&path_indexes[pi]);
                            }
                            free(path_index_disabled);
                            free(path_index_ready);
                            free(path_indexes);
                            free(group_second);
                            free(group_count);
                            free(group_z);
                            free(group_y);
                            free(group_x);
                            free(group_indices);
                            free(hits);
                            return -1;
                        }
                        bool group_missing = false;
                        for (int h = 0; h < exact_hit_count; h++) {
                            if (path_hit_list_contains_cell(
                                    verify_bucket.hits, verify_bucket.count,
                                    hits[h].cell_index))
                                continue;
                            if (path_tile_bucket_append(&verify_bucket,
                                                        &hits[h]) != 0 ||
                                verify_bucket.count > (size_t)INT_MAX) {
                                path_tile_bucket_free(&verify_bucket);
                                if (path_indexes) {
                                    for (size_t pi = 0; pi < paths->count; pi++)
                                        path_slice_index_free(&path_indexes[pi]);
                                }
                                free(path_index_disabled);
                                free(path_index_ready);
                                free(path_indexes);
                                free(group_second);
                                free(group_count);
                                free(group_z);
                                free(group_y);
                                free(group_x);
                                free(group_indices);
                                free(hits);
                                return -1;
                            }
                            if (!group_missing) {
                                g_tile_coverage_stats.path_2d_missing_tiles++;
                                group_missing = true;
                            }
                            g_tile_coverage_stats.path_2d_missing_candidates++;
                            if (verify_2d_log_count < verify_2d_log_limit) {
                                fprintf(stderr,
                                        "[ALEA_PATH_2D_VERIFY] path=%u universe=%d tile=(%d,%d) missing cell_id=%d cell_index=%u bucket_hits=%d exact_hits=%d bbox=[%.9g %.9g]x[%.9g %.9g]x[%.9g %.9g]\n",
                                        path_id, path->universe_id,
                                        ti / tile_w, tj / tile_h,
                                        hits[h].cell_id,
                                        hits[h].cell_index,
                                        hit_count, exact_hit_count,
                                        local_bbox.min_x, local_bbox.max_x,
                                        local_bbox.min_y, local_bbox.max_y,
                                        local_bbox.min_z, local_bbox.max_z);
                                verify_2d_log_count++;
                            }
                        }
                        eval_hits = verify_bucket.hits;
                        hit_count = (int)verify_bucket.count;
                        if ((size_t)exact_hit_count >= max_hits)
                            candidate_query_saturated = true;
                    }
                    if (!using_2d_bucket) {
                        hit_count = alea_hier_spatial_query_universe_region(
                            sys, path->universe_id, &local_bbox, hits, max_hits);
                        if (hit_count < 0) {
                            g_tile_coverage_stats.query_errors++;
                            path_tile_bucket_free(&verify_bucket);
                            if (path_indexes) {
                                for (size_t pi = 0; pi < paths->count; pi++)
                                    path_slice_index_free(&path_indexes[pi]);
                            }
                            free(path_index_disabled);
                            free(path_index_ready);
                            free(path_indexes);
                            free(group_second);
                            free(group_count);
                            free(group_z);
                            free(group_y);
                            free(group_x);
                            free(group_indices);
                            free(hits);
                            return -1;
                        }
                        if ((size_t)hit_count >= max_hits)
                            candidate_query_saturated = true;
                    }

                    g_tile_coverage_stats.candidate_total += (size_t)hit_count;
                    if ((size_t)hit_count > g_tile_coverage_stats.candidate_max)
                        g_tile_coverage_stats.candidate_max = (size_t)hit_count;

                    if (candidate_query_saturated) {
                        g_tile_coverage_stats.fallback_tiles++;
                        g_tile_coverage_stats.exact_fallback_pixels += group_pixels;
                        for (int pj = tj; pj < j_end; pj++) {
                            for (int pi = ti; pi < i_end; pi++) {
                                size_t pidx = (size_t)pj * (size_t)nu + (size_t)pi;
                                if (path_ids[pidx] != path_id) continue;
                                updated += update_pixel_coverage_exact(
                                    sys, &view->plane, view->u_min,
                                    view->v_min, du, dv, nu, nv,
                                    universe_depth, out_secondary_cell_ids,
                                    coverage, errors, pi, pj);
                            }
                        }
                        continue;
                    }

                    if (!using_2d_bucket)
                        hit_count = dedup_spatial_hits(hits, hit_count);
                    g_tile_coverage_stats.path_groups++;
                    if (group_pixels >
                        g_tile_coverage_stats.path_group_pixels_max)
                        g_tile_coverage_stats.path_group_pixels_max =
                            group_pixels;
                    if ((size_t)hit_count >
                        g_tile_coverage_stats.path_group_candidates_max)
                        g_tile_coverage_stats.path_group_candidates_max =
                            (size_t)hit_count;
                    g_tile_coverage_stats.dedup_candidate_total +=
                        (size_t)hit_count;
                    if ((size_t)hit_count >
                        g_tile_coverage_stats.dedup_candidate_max)
                        g_tile_coverage_stats.dedup_candidate_max =
                            (size_t)hit_count;
                    g_tile_coverage_stats.pixels += group_pixels;
                    size_t bbox_tests = 0;
                    size_t bbox_rejects = 0;
                    size_t contains_tests = 0;
                    size_t early_multi_skips = 0;
                    size_t primary_cell_skips = 0;
                    int rc = update_path_group_coverage_from_candidates(
                        sys, eval_hits, hit_count, group_pixels,
                        group_indices, primary_cell_ids,
                        group_x, group_y, group_z,
                        group_count, group_second, out_secondary_cell_ids,
                        coverage, errors, &bbox_tests, &bbox_rejects,
                        &contains_tests, &early_multi_skips,
                        &primary_cell_skips);
                    if (rc < 0) {
                        path_tile_bucket_free(&verify_bucket);
                        if (path_indexes) {
                            for (size_t pi = 0; pi < paths->count; pi++)
                                path_slice_index_free(&path_indexes[pi]);
                        }
                        free(path_index_disabled);
                        free(path_index_ready);
                        free(path_indexes);
                        free(group_second);
                        free(group_count);
                        free(group_z);
                        free(group_y);
                        free(group_x);
                        free(group_indices);
                        free(hits);
                        return -1;
                    }
                    g_tile_coverage_stats.bbox_pixel_tests += bbox_tests;
                    g_tile_coverage_stats.bbox_pixel_rejects += bbox_rejects;
                    g_tile_coverage_stats.candidate_pixel_tests +=
                        contains_tests;
                    g_tile_coverage_stats.early_multi_skips +=
                        early_multi_skips;
                    g_tile_coverage_stats.primary_cell_skips +=
                        primary_cell_skips;
                    updated += rc;
                }
            }
        }
    }

    if (path_indexes) {
        for (size_t pi = 0; pi < paths->count; pi++)
            path_slice_index_free(&path_indexes[pi]);
    }
    free(path_index_disabled);
    free(path_index_ready);
    free(path_indexes);
    free(group_second);
    free(group_count);
    free(group_z);
    free(group_y);
    free(group_x);
    free(group_indices);
    path_tile_bucket_free(&verify_bucket);
    free(hits);
    g_tile_coverage_stats.refined_pixels = (size_t)updated;
    return updated;
}

/* ============================================================================
 * LABEL POSITION COMPUTATION
 *
 * Uses connected component labeling to handle:
 * 1. Disconnected regions with same cell ID (e.g., torus -> two circles)
 * 2. Interior-only sampling to avoid placing labels on borders
 * ============================================================================ */

#define LABEL_MAX_SAMPLES 64   /* Sample points per region for label placement */
#define LABEL_MAX_REGIONS 10000

typedef struct {
    int id;                   /* Cell/material ID */
    double sum_x, sum_y;      /* For centroid approximation */
    int pixel_count;          /* Total pixels in component */
    int interior_count;       /* Interior pixels (for sampling probability) */
    /* Sample points that are interior (not on any border) */
    int sample_x[LABEL_MAX_SAMPLES];
    int sample_y[LABEL_MAX_SAMPLES];
    int sample_count;
    uint32_t rng_state;       /* Per-region LCG state for reservoir sampling */
} label_region_t;

/** @brief Fast LCG — returns value in [0, bound) */
static inline int label_rng_bounded(uint32_t* state, int bound) {
    *state = *state * 1664525u + 1013904223u;
    return (int)(((uint64_t)*state * (uint64_t)bound) >> 32);
}

/**
 * @brief Check if pixel is interior (all 8 neighbors have same cell ID)
 */
static int is_interior_pixel(const int* ids, int width, int height, int x, int y) {
    int id = ids[y * width + x];
    for (int dy = -1; dy <= 1; dy++) {
        for (int dx = -1; dx <= 1; dx++) {
            if (dx == 0 && dy == 0) continue;
            int nx = x + dx, ny = y + dy;
            if (nx < 0 || nx >= width || ny < 0 || ny >= height) return 0;
            if (ids[ny * width + nx] != id) return 0;
        }
    }
    return 1;
}

/**
 * @brief Add an interior sample point using reservoir sampling
 */
static void label_region_add_interior_sample(label_region_t* r, int x, int y) {
    r->interior_count++;
    if (r->sample_count < LABEL_MAX_SAMPLES) {
        r->sample_x[r->sample_count] = x;
        r->sample_y[r->sample_count] = y;
        r->sample_count++;
    } else {
        /* Reservoir sampling: replace with probability LABEL_MAX_SAMPLES/interior_count */
        int idx = label_rng_bounded(&r->rng_state, r->interior_count);
        if (idx < LABEL_MAX_SAMPLES) {
            r->sample_x[idx] = x;
            r->sample_y[idx] = y;
        }
    }
}

/**
 * @brief Find sample point closest to centroid
 */
static void label_find_best_position(const label_region_t* r, int* out_x, int* out_y) {
    if (r->pixel_count == 0) {
        *out_x = 0;
        *out_y = 0;
        return;
    }

    double centroid_x = r->sum_x / r->pixel_count;
    double centroid_y = r->sum_y / r->pixel_count;

    if (r->sample_count == 0) {
        /* No interior samples - fallback to centroid (shouldn't happen for large regions) */
        *out_x = (int)centroid_x;
        *out_y = (int)centroid_y;
        return;
    }

    /* Find interior sample point closest to centroid */
    double best_dist = 1e30;
    int best_idx = 0;

    for (int i = 0; i < r->sample_count; i++) {
        double dx = r->sample_x[i] - centroid_x;
        double dy = r->sample_y[i] - centroid_y;
        double dist = dx*dx + dy*dy;
        if (dist < best_dist) {
            best_dist = dist;
            best_idx = i;
        }
    }

    *out_x = r->sample_x[best_idx];
    *out_y = r->sample_y[best_idx];
}

/**
 * @brief Iterative flood-fill to label a connected component
 *
 * Uses a queue-based approach to avoid stack overflow on large regions.
 */
static void flood_fill_component(const int* ids, int* component_map,
                                  int width, int height,
                                  int start_x, int start_y,
                                  int component_id, int target_id,
                                  int* queue, int queue_cap) {
    /* Uses caller-provided queue [x0,y0, x1,y1, ...] to avoid
     * per-component allocation. queue_cap must be >= width*height*2. */
    int queue_size = 0;

    queue[queue_size++] = start_x;
    queue[queue_size++] = start_y;
    component_map[start_y * width + start_x] = component_id;

    int head = 0;
    while (head < queue_size) {
        int x = queue[head++];
        int y = queue[head++];

        const int dx[] = {-1, 1, 0, 0};
        const int dy[] = {0, 0, -1, 1};

        for (int d = 0; d < 4; d++) {
            int nx = x + dx[d];
            int ny = y + dy[d];

            if (nx < 0 || nx >= width || ny < 0 || ny >= height) continue;

            int nidx = ny * width + nx;
            if (ids[nidx] == target_id && component_map[nidx] == 0) {
                component_map[nidx] = component_id;
                if (queue_size + 2 <= queue_cap) {
                    queue[queue_size++] = nx;
                    queue[queue_size++] = ny;
                }
            }
        }
    }
}

int alea_find_label_positions(
    const int* ids,
    int width, int height,
    int min_pixels,
    alea_label_position_t** out_labels,
    int* out_count)
{
    if (!ids || !out_labels || !out_count || width <= 0 || height <= 0) {
        return -1;
    }

    *out_labels = NULL;
    *out_count = 0;

    int num_pixels = width * height;

    /* Step 1: Connected component labeling */
    int* component_map = calloc(num_pixels, sizeof(int));
    if (!component_map) return -1;

    /* Pre-allocate flood fill queue once (worst case: all pixels in one component) */
    int queue_cap = num_pixels * 2;
    int* queue = malloc(queue_cap * sizeof(int));
    if (!queue) {
        free(component_map);
        return -1;
    }

    int next_component = 1;  /* Component IDs start at 1; 0 = unlabeled */

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int idx = y * width + x;
            int id = ids[idx];

            /* Skip void and already-labeled pixels */
            if (id < 0 || component_map[idx] != 0) continue;

            /* Start a new component */
            flood_fill_component(ids, component_map, width, height,
                                x, y, next_component, id,
                                queue, queue_cap);
            next_component++;

            if (next_component >= LABEL_MAX_REGIONS) break;
        }
        if (next_component >= LABEL_MAX_REGIONS) break;
    }

    free(queue);

    int num_components = next_component - 1;
    if (num_components == 0) {
        free(component_map);
        return 0;
    }

    /* Step 2: Collect statistics for each component */
    label_region_t* regions = calloc(num_components, sizeof(label_region_t));
    if (!regions) {
        free(component_map);
        return -1;
    }

    /* Initialize regions */
    for (int i = 0; i < num_components; i++) {
        regions[i].id = -1;
        regions[i].sum_x = regions[i].sum_y = 0;
        regions[i].pixel_count = 0;
        regions[i].interior_count = 0;
        regions[i].sample_count = 0;
    }

    /* Collect per-component data */
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int idx = y * width + x;
            int comp = component_map[idx];
            if (comp == 0) continue;

            label_region_t* r = &regions[comp - 1];

            /* Store cell ID (first pixel sets it) */
            if (r->id < 0) {
                r->id = ids[idx];
                r->rng_state = (uint32_t)(comp * 2654435761u);  /* seed from component ID */
            }

            /* Accumulate for centroid */
            r->sum_x += x;
            r->sum_y += y;
            r->pixel_count++;

            /* Only sample interior pixels */
            if (is_interior_pixel(ids, width, height, x, y)) {
                label_region_add_interior_sample(r, x, y);
            }
        }
    }

    free(component_map);

    /* Step 3: Count valid components (meeting minimum size AND having interior samples) */
    int valid_count = 0;
    for (int i = 0; i < num_components; i++) {
        if (regions[i].pixel_count >= min_pixels && regions[i].sample_count > 0) {
            valid_count++;
        }
    }

    if (valid_count == 0) {
        free(regions);
        return 0;
    }

    /* Step 4: Allocate and fill output array */
    alea_label_position_t* labels = malloc(valid_count * sizeof(alea_label_position_t));
    if (!labels) {
        free(regions);
        return -1;
    }

    int out_idx = 0;
    for (int i = 0; i < num_components; i++) {
        if (regions[i].pixel_count >= min_pixels && regions[i].sample_count > 0) {
            labels[out_idx].id = regions[i].id;
            labels[out_idx].pixel_count = regions[i].pixel_count;
            label_find_best_position(&regions[i], &labels[out_idx].px, &labels[out_idx].py);
            out_idx++;
        }
    }

    free(regions);
    *out_labels = labels;
    *out_count = valid_count;
    return 0;
}

/**
 * @brief Evaluate curve position at parameter t
 */
static int eval_curve_at_t(const alea_curve_t* curve, double t,
                            double* out_x, double* out_y) {
    switch (curve->type) {
        case ALEA_CURVE_LINE:
        case ALEA_CURVE_LINE_SEGMENT:
            *out_x = curve->data.line.point[0] + t * curve->data.line.direction[0];
            *out_y = curve->data.line.point[1] + t * curve->data.line.direction[1];
            return 0;

        case ALEA_CURVE_CIRCLE:
        case ALEA_CURVE_ARC:
            *out_x = curve->data.circle.center[0] + curve->data.circle.radius * cos(t);
            *out_y = curve->data.circle.center[1] + curve->data.circle.radius * sin(t);
            return 0;

        case ALEA_CURVE_ELLIPSE:
        case ALEA_CURVE_ELLIPSE_ARC: {
            double c = cos(curve->data.ellipse.angle);
            double s = sin(curve->data.ellipse.angle);
            double lx = curve->data.ellipse.semi_a * cos(t);
            double ly = curve->data.ellipse.semi_b * sin(t);
            *out_x = curve->data.ellipse.center[0] + c * lx - s * ly;
            *out_y = curve->data.ellipse.center[1] + s * lx + c * ly;
            return 0;
        }

        default:
            return -1;
    }
}

/**
 * @brief Compute curve length estimate (for minimum length filtering)
 */
static double estimate_curve_length(const alea_curve_t* curve,
                                     double x_min, double x_max,
                                     double y_min, double y_max,
                                     int width, int height) {
    double dx = (x_max - x_min) / width;

    /* For circles: full circle if t_min == t_max, otherwise use arc angle */
    if (curve->type == ALEA_CURVE_CIRCLE) {
        /* Full circle - circumference in pixels */
        double r_pixels = curve->data.circle.radius / dx;
        return 2.0 * M_PI * r_pixels;
    }

    if (curve->type == ALEA_CURVE_ARC) {
        double r_pixels = curve->data.circle.radius / dx;
        double arc_angle = curve->t_max - curve->t_min;
        if (arc_angle <= 0) arc_angle += 2 * M_PI;
        return r_pixels * arc_angle;
    }

    /* For ellipses: approximate using average radius */
    if (curve->type == ALEA_CURVE_ELLIPSE) {
        double avg_r = (curve->data.ellipse.semi_a + curve->data.ellipse.semi_b) / 2.0;
        double r_pixels = avg_r / dx;
        return 2.0 * M_PI * r_pixels;
    }

    if (curve->type == ALEA_CURVE_ELLIPSE_ARC) {
        double avg_r = (curve->data.ellipse.semi_a + curve->data.ellipse.semi_b) / 2.0;
        double r_pixels = avg_r / dx;
        double arc_angle = curve->t_max - curve->t_min;
        if (arc_angle <= 0) arc_angle += 2 * M_PI;
        return r_pixels * arc_angle;
    }

    /* For lines: use endpoint distance */
    if (curve->type == ALEA_CURVE_LINE ||
        curve->type == ALEA_CURVE_LINE_SEGMENT) {
        double x1, y1, x2, y2;
        if (eval_curve_at_t(curve, curve->t_min, &x1, &y1) != 0) return 0;
        if (eval_curve_at_t(curve, curve->t_max, &x2, &y2) != 0) return 0;

        double dy = (y_max - y_min) / height;
        double px1 = (x1 - x_min) / dx;
        double py1 = (y1 - y_min) / dy;
        double px2 = (x2 - x_min) / dx;
        double py2 = (y2 - y_min) / dy;

        double pdx = px2 - px1;
        double pdy = py2 - py1;
        double len = sqrt(pdx*pdx + pdy*pdy);

        /* For unbounded lines, return large value */
        if (curve->type == ALEA_CURVE_LINE) return 1e6;
        return len;
    }

    /* Default: assume large enough */
    return 1e6;
}

static int contour_drawn_at_grid_pixel(const int* boundary_ids,
                                       int width, int height,
                                       int ix, int iy) {
    if (!boundary_ids || ix < 0 || ix >= width || iy < 0 || iy >= height) {
        return 0;
    }

    int idx = iy * width + ix;
    int id = boundary_ids[idx];

    if (ix + 1 < width && boundary_ids[idx + 1] != id) return 1;
    if (iy > 0 && boundary_ids[(iy - 1) * width + ix] != id) return 1;
    return 0;
}

static int point_has_drawn_contour_nearby(const int* boundary_ids,
                                          int width, int height,
                                          int ix, int iy) {
    if (!boundary_ids) return 1;

    /* Match draw_contours_ex(), with a one-pixel tolerance because a curve
     * point may quantize to either side of the grid edge where the contour is
     * actually rasterized. */
    for (int dy = -1; dy <= 1; dy++) {
        for (int dx = -1; dx <= 1; dx++) {
            if (contour_drawn_at_grid_pixel(boundary_ids, width, height,
                                            ix + dx, iy + dy)) {
                return 1;
            }
        }
    }
    return 0;
}

static void surface_label_param_range(const alea_curve_t* curve,
                                      double* t_start,
                                      double* t_end,
                                      double* t_preferred) {
    if (curve->type == ALEA_CURVE_CIRCLE || curve->type == ALEA_CURVE_ELLIPSE) {
        *t_start = 0.0;
        *t_end = 2.0 * M_PI;
        *t_preferred = M_PI / 4.0;
    } else if (curve->t_min == curve->t_max) {
        *t_start = 0.0;
        *t_end = 0.0;
        *t_preferred = 0.0;
    } else {
        *t_start = curve->t_min;
        *t_end = curve->t_max;
        *t_preferred = (curve->t_min + curve->t_max) / 2.0;
    }
}

static int surface_label_candidate(const alea_curve_t* curve,
                                   const int* boundary_ids,
                                   double x_min, double y_min,
                                   double dx, double dy,
                                   int width, int height,
                                   int margin,
                                   int* out_ix, int* out_iy) {
    double t_start, t_end, t_preferred;
    surface_label_param_range(curve, &t_start, &t_end, &t_preferred);

    int samples = boundary_ids ? 32 : 0;
    double best_score = 1e30;
    int found = 0;

    for (int s = 0; s <= samples; s++) {
        double t;
        if (samples == 0) {
            t = t_preferred;
        } else {
            t = t_start + (t_end - t_start) * (double)s / (double)samples;
        }

        double px, py;
        if (eval_curve_at_t(curve, t, &px, &py) != 0) continue;

        int ix = (int)((px - x_min) / dx);
        int iy = (int)((py - y_min) / dy);

        if (ix < margin || ix >= width - margin ||
            iy < margin || iy >= height - margin) {
            continue;
        }

        if (!point_has_drawn_contour_nearby(boundary_ids, width, height, ix, iy)) {
            continue;
        }

        double score = fabs(t - t_preferred);
        if (curve->type == ALEA_CURVE_CIRCLE || curve->type == ALEA_CURVE_ELLIPSE) {
            double period = 2.0 * M_PI;
            score = fmod(score, period);
            if (score > period * 0.5) score = period - score;
        }

        if (!found || score < best_score) {
            found = 1;
            best_score = score;
            *out_ix = ix;
            *out_iy = iy;
        }
    }

    return found;
}

int alea_find_surface_label_positions_on_boundaries(
    const alea_slice_curves_t* curves,
    const int* boundary_ids,
    double x_min, double x_max,
    double y_min, double y_max,
    int width, int height,
    int margin,
    alea_label_position_t** out_labels,
    int* out_count)
{
    if (!curves || !out_labels || !out_count || width <= 0 || height <= 0) {
        return -1;
    }

    *out_labels = NULL;
    *out_count = 0;

    size_t num_curves = alea_slice_curves_count(curves);
    if (num_curves == 0) return 0;

    double dx = (x_max - x_min) / width;
    double dy = (y_max - y_min) / height;

    /* Minimum curve length in pixels to get a label (avoid tiny fragments) */
    const double MIN_CURVE_LENGTH = 30.0;

    /* Minimum distance between labels of same surface (in pixels) */
    const double MIN_LABEL_SPACING = 50.0;

    /* Temporary storage for valid label positions */
    typedef struct {
        int surf_id;
        int px, py;
    } temp_label_t;

    temp_label_t* temp_labels = malloc(num_curves * sizeof(temp_label_t));
    if (!temp_labels) return -1;
    int temp_count = 0;

    /* First pass: collect all valid label positions */
    for (size_t i = 0; i < num_curves; i++) {
        alea_curve_t curve;
        if (alea_slice_curves_get(curves, i, &curve) != 0) continue;

        int surf_id = curve.surface_id;
        if (surf_id <= 0) continue;

        /* Skip curves that are too short */
        double curve_len = estimate_curve_length(&curve, x_min, x_max, y_min, y_max,
                                                  width, height);
        if (curve_len < MIN_CURVE_LENGTH) continue;

        int ix, iy;
        if (!surface_label_candidate(&curve, boundary_ids,
                                     x_min, y_min, dx, dy,
                                     width, height, margin,
                                     &ix, &iy)) {
            continue;
        }

        /* Check spacing from existing labels of the same surface */
        int too_close = 0;
        for (int j = 0; j < temp_count; j++) {
            if (temp_labels[j].surf_id == surf_id) {
                double ddx = ix - temp_labels[j].px;
                double ddy = iy - temp_labels[j].py;
                if (ddx*ddx + ddy*ddy < MIN_LABEL_SPACING * MIN_LABEL_SPACING) {
                    too_close = 1;
                    break;
                }
            }
        }
        if (too_close) continue;

        temp_labels[temp_count].surf_id = surf_id;
        temp_labels[temp_count].px = ix;
        temp_labels[temp_count].py = iy;
        temp_count++;
    }

    if (temp_count == 0) {
        free(temp_labels);
        return 0;
    }

    /* Allocate and fill output */
    alea_label_position_t* labels = malloc(temp_count * sizeof(alea_label_position_t));
    if (!labels) {
        free(temp_labels);
        return -1;
    }

    for (int i = 0; i < temp_count; i++) {
        labels[i].id = temp_labels[i].surf_id;
        labels[i].px = temp_labels[i].px;
        labels[i].py = temp_labels[i].py;
        labels[i].pixel_count = 0;
    }

    free(temp_labels);
    *out_labels = labels;
    *out_count = temp_count;
    return 0;
}

int alea_find_surface_label_positions(
    const alea_slice_curves_t* curves,
    double x_min, double x_max,
    double y_min, double y_max,
    int width, int height,
    int margin,
    alea_label_position_t** out_labels,
    int* out_count)
{
    return alea_find_surface_label_positions_on_boundaries(
        curves, NULL, x_min, x_max, y_min, y_max,
        width, height, margin, out_labels, out_count);
}

/* ============================================================================
 * ANALYTICAL ERROR LINE CHECKING
 * ============================================================================ */

/**
 * Count cells at a specific depth in the hit array.
 * depth == -1 means innermost (max depth hit).
 * Returns the number of distinct cells at that depth.
 */
/**
 * Classify a sample point on a curve by offsetting to both sides and
 * counting cells. Returns 0 = OK, ALEA_SLICE_ERR_OVERLAP, or ALEA_SLICE_ERR_GAP.
 */
/**
 * Count how many cells in a given universe contain a 3D point.
 * Returns the count (0 = gap, 1 = OK, >1 = overlap).
 *
 * Uses the spatial index (BVH + coherence cache) for fast lookup
 * instead of brute-forcing all cells in the universe.
 */
static int count_cells_in_universe_at_point(alea_system_t* sys,
                                             int universe_id,
                                             double x, double y, double z) {
    enum { MAX_HITS = 64 };
    alea_cell_hit_t hits[MAX_HITS];
    int n;

    if (alea_spatial_mode_is_hierarchical(sys)) {
        n = alea_hier_spatial_find_cells_at_point_uncached(sys, x, y, z,
                                                           hits, MAX_HITS);
    } else {
        n = alea_spatial_find_cells_at_point(sys, x, y, z, hits, MAX_HITS);
    }
    if (n <= 0) return 0;

    int count = 0;
    for (int i = 0; i < n; i++) {
        if (hits[i].universe_id == universe_id) {
            count++;
            if (count > 1) return count;  /* Early exit on overlap */
        }
    }
    if (n >= MAX_HITS && count > 0) return 2;  /* Conservative on truncation. */
    return count;
}

static int classify_sample(alea_system_t* sys,
                           const alea_slice_view_t* view,
                           double u, double v,
                           double nu, double nv,
                           double eps,
                           int universe_id) {
    const alea_slice_plane_t* pl = &view->plane;

    /* Offset points in 2D */
    double u_plus  = u + eps * nu;
    double v_plus  = v + eps * nv;
    double u_minus = u - eps * nu;
    double v_minus = v - eps * nv;

    /* Convert to 3D */
    double x_plus  = pl->origin[0] + u_plus * pl->u_axis[0] + v_plus * pl->v_axis[0];
    double y_plus  = pl->origin[1] + u_plus * pl->u_axis[1] + v_plus * pl->v_axis[1];
    double z_plus  = pl->origin[2] + u_plus * pl->u_axis[2] + v_plus * pl->v_axis[2];

    double x_minus = pl->origin[0] + u_minus * pl->u_axis[0] + v_minus * pl->v_axis[0];
    double y_minus = pl->origin[1] + u_minus * pl->u_axis[1] + v_minus * pl->v_axis[1];
    double z_minus = pl->origin[2] + u_minus * pl->u_axis[2] + v_minus * pl->v_axis[2];

    /* Count cells in this universe on both sides */
    int cells_plus  = count_cells_in_universe_at_point(sys, universe_id,
                                                        x_plus, y_plus, z_plus);
    int cells_minus = count_cells_in_universe_at_point(sys, universe_id,
                                                        x_minus, y_minus, z_minus);

    if (cells_plus > 1 || cells_minus > 1)
        return ALEA_SLICE_ERR_OVERLAP;
    if (cells_plus == 0 || cells_minus == 0)
        return ALEA_SLICE_ERR_GAP;
    return 0; /* OK */
}

/**
 * Compute parameter range for a curve clipped to the viewport.
 * For bounded curves, returns the stored bounds.
 * For infinite curves (LINE, PARABOLA, HYPERBOLA), clips to viewport.
 */
static void get_clipped_param_range(const alea_curve_2d_t* curve,
                                    const alea_slice_view_t* view,
                                    double* t_lo, double* t_hi) {
    switch (curve->type) {
        case ALEA_CURVE_LINE: {
            /* For infinite lines, find parameter range that covers viewport */
            const alea_line_2d_t* ln = &curve->data.line;
            double du = ln->direction[0];
            double dv = ln->direction[1];

            /* Find t range where point + t*dir intersects viewport */
            double t_min_val = -1e10, t_max_val = 1e10;

            if (fabs(du) > 1e-15) {
                double t1 = (view->u_min - ln->point[0]) / du;
                double t2 = (view->u_max - ln->point[0]) / du;
                if (t1 > t2) { double tmp = t1; t1 = t2; t2 = tmp; }
                if (t1 > t_min_val) t_min_val = t1;
                if (t2 < t_max_val) t_max_val = t2;
            }
            if (fabs(dv) > 1e-15) {
                double t1 = (view->v_min - ln->point[1]) / dv;
                double t2 = (view->v_max - ln->point[1]) / dv;
                if (t1 > t2) { double tmp = t1; t1 = t2; t2 = tmp; }
                if (t1 > t_min_val) t_min_val = t1;
                if (t2 < t_max_val) t_max_val = t2;
            }

            *t_lo = t_min_val;
            *t_hi = t_max_val;
            return;
        }

        case ALEA_CURVE_PARALLEL_LINES:
        case ALEA_CURVE_RAY: {
            /* For rays and parallel lines, clip similar to lines */
            const alea_line_2d_t* ln = &curve->data.line;
            double du = ln->direction[0];
            double dv = ln->direction[1];
            double t_min_val = -1e10, t_max_val = 1e10;

            if (fabs(du) > 1e-15) {
                double t1 = (view->u_min - ln->point[0]) / du;
                double t2 = (view->u_max - ln->point[0]) / du;
                if (t1 > t2) { double tmp = t1; t1 = t2; t2 = tmp; }
                if (t1 > t_min_val) t_min_val = t1;
                if (t2 < t_max_val) t_max_val = t2;
            }
            if (fabs(dv) > 1e-15) {
                double t1 = (view->v_min - ln->point[1]) / dv;
                double t2 = (view->v_max - ln->point[1]) / dv;
                if (t1 > t2) { double tmp = t1; t1 = t2; t2 = tmp; }
                if (t1 > t_min_val) t_min_val = t1;
                if (t2 < t_max_val) t_max_val = t2;
            }

            if (curve->type == ALEA_CURVE_RAY && t_min_val < curve->bounds.t_min)
                t_min_val = curve->bounds.t_min;

            *t_lo = t_min_val;
            *t_hi = t_max_val;
            return;
        }

        case ALEA_CURVE_CIRCLE:
        case ALEA_CURVE_ELLIPSE:
            /* Full closed curves: parameter range is [0, 2*PI] */
            *t_lo = 0;
            *t_hi = 2.0 * M_PI;
            return;

        default:
            /* Bounded curves: use stored parameter range */
            *t_lo = curve->bounds.t_min;
            *t_hi = curve->bounds.t_max;
            return;
    }
}

/* Public accessors so external consumers (geometry validator) can sample
 * boundary curves without reaching into slice internals. */
int alea_slice_curve_eval(const alea_slice_curves_t* curves, size_t index,
                          double t, double* u, double* v) {
    if (!curves || !u || !v || index >= curves->internal.curves.count)
        return -1;
    const alea_curve_2d_t* c = &curves->internal.curves.data[index];
    return alea_curve_eval(c, t, u, v) ? 0 : -1;
}

void alea_slice_curve_param_range(const alea_slice_curves_t* curves, size_t index,
                                  const alea_slice_view_t* view,
                                  double* t_lo, double* t_hi) {
    if (t_lo) *t_lo = 0.0;
    if (t_hi) *t_hi = 0.0;
    if (!curves || !view || !t_lo || !t_hi ||
        index >= curves->internal.curves.count)
        return;
    get_clipped_param_range(&curves->internal.curves.data[index], view,
                            t_lo, t_hi);
}

/**
 * Check if a 2D point is inside the viewport (with small margin).
 */
static bool point_in_viewport(double u, double v, const alea_slice_view_t* view) {
    return u >= view->u_min && u <= view->u_max &&
           v >= view->v_min && v <= view->v_max;
}

/* Dynamic array for error segments */
ALEA_VEC_DEFINE(error_vec, alea_slice_error_t);
ALEA_VEC_DEFINE(plot_error_component_vec, alea_plot_error_component_t);

static void error_vec_push(error_vec_t* vec, const alea_slice_error_t* err) {
    alea_vec_push(vec, *err, alea_slice_error_t);
}

static void plot_error_component_vec_push(plot_error_component_vec_t* vec,
                                          const alea_plot_error_component_t* comp) {
    alea_vec_push(vec, *comp, alea_plot_error_component_t);
}

/**
 * Process a single parametric curve: sample, classify, merge error segments.
 */
static void check_parametric_curve(alea_system_t* sys,
                                   const alea_slice_view_t* view,
                                   const alea_curve_2d_t* curve,
                                   size_t curve_index,
                                   double sample_spacing,
                                   double eps,
                                   int universe_id,
                                   error_vec_t* errors) {
    double t_lo, t_hi;
    get_clipped_param_range(curve, view, &t_lo, &t_hi);

    if (t_lo >= t_hi) return;

    double arc_approx = t_hi - t_lo;
    /* For angular curves (circle/ellipse), scale by radius */
    if (curve->type == ALEA_CURVE_CIRCLE || curve->type == ALEA_CURVE_ARC) {
        arc_approx = (t_hi - t_lo) * curve->data.circle.radius;
    } else if (curve->type == ALEA_CURVE_ELLIPSE || curve->type == ALEA_CURVE_ELLIPSE_ARC) {
        double avg_r = (curve->data.ellipse.semi_a + curve->data.ellipse.semi_b) * 0.5;
        arc_approx = (t_hi - t_lo) * avg_r;
    }

    int n_samples = (int)(arc_approx / sample_spacing);
    if (n_samples < 2) n_samples = 2;
    if (n_samples > 2000) n_samples = 2000;

    double dt_finite = (t_hi - t_lo) * 1e-6;
    if (dt_finite < 1e-15) dt_finite = 1e-15;


    /* Track current error segment */
    int cur_type = 0;      /* 0 = OK or not started */
    double seg_t_start = 0;

    for (int i = 0; i <= n_samples; i++) {
        double t = t_lo + (t_hi - t_lo) * i / n_samples;

        double u, v;
        if (!alea_curve_eval(curve, t, &u, &v))
            continue;

        if (!point_in_viewport(u, v, view))
            continue;

        /* Compute tangent via finite difference */
        double u1, v1, u2, v2;
        double t_back = t - dt_finite;
        double t_fwd  = t + dt_finite;
        if (t_back < t_lo) t_back = t_lo;
        if (t_fwd > t_hi) t_fwd = t_hi;

        if (!alea_curve_eval(curve, t_back, &u1, &v1) ||
            !alea_curve_eval(curve, t_fwd,  &u2, &v2))
            continue;

        double tang_u = u2 - u1;
        double tang_v = v2 - v1;
        double tang_len = sqrt(tang_u * tang_u + tang_v * tang_v);
        if (tang_len < 1e-20) continue;

        /* Normal = (-tang_v, tang_u) / len */
        double nu = -tang_v / tang_len;
        double nv =  tang_u / tang_len;

        int classification = classify_sample(sys, view, u, v, nu, nv,
                                             eps, universe_id);

        if (classification != cur_type) {
            /* Close previous segment if it was an error */
            if (cur_type != 0) {
                alea_slice_error_t err = {
                    .curve_index = curve_index,
                    .surface_id = curve->surface_id,
                    .type = (alea_slice_error_type_t)cur_type,
                    .t_start = seg_t_start,
                    .t_end = t
                };
                error_vec_push(errors, &err);
            }
            cur_type = classification;
            seg_t_start = t;
        }
    }

    /* Close final segment */
    if (cur_type != 0) {
        alea_slice_error_t err = {
            .curve_index = curve_index,
            .surface_id = curve->surface_id,
            .type = (alea_slice_error_type_t)cur_type,
            .t_start = seg_t_start,
            .t_end = t_hi
        };
        error_vec_push(errors, &err);
    }
}

/**
 * Process a quartic (torus) curve via scanline sampling.
 */
static void check_quartic_curve(alea_system_t* sys,
                                const alea_slice_view_t* view,
                                const alea_curve_2d_t* curve,
                                size_t curve_index,
                                double sample_spacing,
                                double eps,
                                int universe_id,
                                error_vec_t* errors) {
    /* Get bounding box for this curve */
    double bbox_umin, bbox_umax, bbox_vmin, bbox_vmax;
    alea_curve_bbox(curve, &bbox_umin, &bbox_umax, &bbox_vmin, &bbox_vmax);

    /* Clip to viewport */
    double v_lo = fmax(bbox_vmin, view->v_min);
    double v_hi = fmin(bbox_vmax, view->v_max);
    if (v_lo >= v_hi) return;

    int n_scanlines = (int)((v_hi - v_lo) / sample_spacing);
    if (n_scanlines < 2) n_scanlines = 2;
    if (n_scanlines > 2000) n_scanlines = 2000;

    for (int i = 0; i <= n_scanlines; i++) {
        double v_scan = v_lo + (v_hi - v_lo) * i / n_scanlines;

        double u_vals[16];
        int n_isect = alea_curve_scanline_intersect(curve, v_scan, u_vals, 16);

        for (int j = 0; j < n_isect; j++) {
            double u = u_vals[j];
            if (u < view->u_min || u > view->u_max) continue;

            /* For quartic, estimate normal from scanline direction (vertical) */
            /* Use a small offset in v to approximate tangent */
            double u_vals2[16];
            double dv = sample_spacing * 0.01;
            int n2 = alea_curve_scanline_intersect(curve, v_scan + dv, u_vals2, 16);

            /* Find closest u in the next scanline */
            double tang_u = 0, tang_v = dv;
            double best_dist = 1e20;
            for (int k = 0; k < n2; k++) {
                double dist = fabs(u_vals2[k] - u);
                if (dist < best_dist) {
                    best_dist = dist;
                    tang_u = u_vals2[k] - u;
                }
            }

            double tang_len = sqrt(tang_u * tang_u + tang_v * tang_v);
            if (tang_len < 1e-20) {
                /* Fallback: use horizontal normal */
                tang_u = 0;
                tang_v = 1;
                tang_len = 1;
            }

            double nu = -tang_v / tang_len;
            double nv =  tang_u / tang_len;

            int classification = classify_sample(sys, view, u, v_scan,
                                                 nu, nv, eps, universe_id);

            if (classification != 0) {
                /* For quartic, use v as parameter proxy */
                alea_slice_error_t err = {
                    .curve_index = curve_index,
                    .surface_id = curve->surface_id,
                    .type = (alea_slice_error_type_t)classification,
                    .t_start = v_scan,
                    .t_end = v_scan
                };
                error_vec_push(errors, &err);
            }
        }
    }
}

/**
 * Process a polygon curve: check each edge.
 */
static void check_polygon_curve(alea_system_t* sys,
                                const alea_slice_view_t* view,
                                const alea_curve_2d_t* curve,
                                size_t curve_index,
                                double sample_spacing,
                                double eps,
                                int universe_id,
                                error_vec_t* errors) {
    const alea_polygon_2d_t* poly = &curve->data.polygon;
    int nv = poly->vertex_count;
    if (nv < 2) return;

    int n_edges = poly->closed ? nv : (nv - 1);

    for (int e = 0; e < n_edges; e++) {
        int i0 = e;
        int i1 = (e + 1) % nv;
        double x0 = poly->vertices[i0][0], y0 = poly->vertices[i0][1];
        double x1 = poly->vertices[i1][0], y1 = poly->vertices[i1][1];

        double edge_len = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
        if (edge_len < 1e-15) continue;

        int n_samples = (int)(edge_len / sample_spacing);
        if (n_samples < 2) n_samples = 2;
        if (n_samples > 500) n_samples = 500;

        /* Edge normal */
        double dx = x1 - x0, dy = y1 - y0;
        double nu = -dy / edge_len;
        double nv_dir = dx / edge_len;

        int cur_type = 0;
        double seg_t_start = 0;

        for (int s = 0; s <= n_samples; s++) {
            double frac = (double)s / n_samples;
            double u = x0 + frac * dx;
            double v = y0 + frac * dy;
            double t_param = i0 + frac;

            if (!point_in_viewport(u, v, view))
                continue;

            int classification = classify_sample(sys, view, u, v,
                                                 nu, nv_dir, eps,
                                                 universe_id);

            if (classification != cur_type) {
                if (cur_type != 0) {
                    alea_slice_error_t err = {
                        .curve_index = curve_index,
                        .surface_id = curve->surface_id,
                        .type = (alea_slice_error_type_t)cur_type,
                        .t_start = seg_t_start,
                        .t_end = t_param
                    };
                    error_vec_push(errors, &err);
                }
                cur_type = classification;
                seg_t_start = t_param;
            }
        }

        if (cur_type != 0) {
            alea_slice_error_t err = {
                .curve_index = curve_index,
                .surface_id = curve->surface_id,
                .type = (alea_slice_error_type_t)cur_type,
                .t_start = seg_t_start,
                .t_end = (double)(i0 + 1)
            };
            error_vec_push(errors, &err);
        }
    }
}

alea_slice_error_result_t* alea_check_slice_errors(
    alea_system_t* sys,
    const alea_slice_view_t* view,
    const alea_slice_curves_t* curves,
    int universe_depth)
{
    (void)universe_depth;  /* universe_id is now stored per-curve */
    if (!sys || !view || !curves) return NULL;

    const alea_curve_collection_t* coll = &curves->internal;
    if (coll->curves.count == 0) return NULL;

    /* Reset spatial coherence cache so it starts fresh for curve walking */
    alea_spatial_reset_cache();

    /* Compute viewport-derived parameters */
    double vp_width  = view->u_max - view->u_min;
    double vp_height = view->v_max - view->v_min;
    double vp_diag = sqrt(vp_width * vp_width + vp_height * vp_height);

    double sample_spacing = vp_diag / 200.0;
    double eps = vp_diag * 1e-6;

    error_vec_t errvec = ALEA_VEC_INIT;

    for (size_t i = 0; i < coll->curves.count; i++) {
        const alea_curve_2d_t* curve = &coll->curves.data[i];

        switch (curve->type) {
            case ALEA_CURVE_NONE:
            case ALEA_CURVE_POINT:
                /* Skip: no meaningful boundary to check */
                continue;

            case ALEA_CURVE_QUARTIC:
                check_quartic_curve(sys, view, curve, i, sample_spacing,
                                    eps, curve->universe_id, &errvec);
                continue;

            case ALEA_CURVE_POLYGON:
                check_polygon_curve(sys, view, curve, i, sample_spacing,
                                    eps, curve->universe_id, &errvec);
                continue;

            default:
                check_parametric_curve(sys, view, curve, i, sample_spacing,
                                       eps, curve->universe_id, &errvec);
                continue;
        }
    }

    /* Build result */
    alea_slice_error_result_t* result = malloc(sizeof(alea_slice_error_result_t));
    if (!result) {
        free(errvec.data);
        return NULL;
    }

    result->errors = errvec.data;
    result->error_count = errvec.count;
    return result;
}

/* ============================================================================
 * GRID-BASED ERROR CHECKING (fast O(1) per sample)
 * ============================================================================ */

/** Grid context passed to grid-based classify */
typedef struct {
    const int* cell_ids;
    const uint8_t* coverage;
    const uint8_t* grid_errors;
    int nu, nv;
    double u_min, v_min;
    double inv_du, inv_dv;  /* 1/du, 1/dv for fast coord→pixel */
    /* For CSG fallback on same-cell cases (nested overlaps) */
    alea_system_t* sys;
    const alea_slice_plane_t* plane;
} grid_ctx_t;

/** Map 2D coordinate to pixel index, return -1 if out of bounds */
static inline int grid_pixel_at(const grid_ctx_t* g, double u, double v) {
    int i = (int)((u - g->u_min) * g->inv_du);
    int j = (int)((v - g->v_min) * g->inv_dv);
    if (i < 0 || i >= g->nu || j < 0 || j >= g->nv) return -1;
    return j * g->nu + i;
}

/**
 * Grid-based classification of a sample point on a curve.
 * Returns 0=OK/skip, ALEA_SLICE_ERR_OVERLAP, or ALEA_SLICE_ERR_GAP.
 *
 * For the common case (different cells on each side), this is pure O(1)
 * grid lookup. When both sides show the same cell (possible nested overlap),
 * falls back to a single CSG query to check for multiple cells.
 */
static int classify_sample_grid(const grid_ctx_t* g,
                                double u, double v,
                                double nu, double nv,
                                double eps,
                                int universe_id) {
    double u_plus  = u + eps * nu;
    double v_plus  = v + eps * nv;
    double u_minus = u - eps * nu;
    double v_minus = v - eps * nv;

    int idx_plus  = grid_pixel_at(g, u_plus, v_plus);
    int idx_minus = grid_pixel_at(g, u_minus, v_minus);

    /* Out of grid bounds → can't classify */
    if (idx_plus < 0 || idx_minus < 0) return 0;

    /* Check grid-level overlap flags first */
    if (g->grid_errors) {
        if (g->grid_errors[idx_plus] == GRID_ERR_OVERLAP ||
            g->grid_errors[idx_minus] == GRID_ERR_OVERLAP)
            return ALEA_SLICE_ERR_OVERLAP;
    }

    if (g->coverage) {
        if (g->coverage[idx_plus] == ALEA_COVERAGE_MULTI ||
            g->coverage[idx_minus] == ALEA_COVERAGE_MULTI)
            return ALEA_SLICE_ERR_OVERLAP;
        if (g->coverage[idx_plus] == ALEA_COVERAGE_NONE ||
            g->coverage[idx_minus] == ALEA_COVERAGE_NONE)
            return ALEA_SLICE_ERR_GAP;
    }

    int cell_plus  = g->cell_ids[idx_plus];
    int cell_minus = g->cell_ids[idx_minus];

    /* Different cells on each side → proper boundary, OK */
    if (cell_plus != cell_minus && cell_plus >= 0 && cell_minus >= 0)
        return 0;

    /* One or both sides void → gap */
    if (cell_plus < 0 || cell_minus < 0)
        return ALEA_SLICE_ERR_GAP;

    /* Same winning cell on both sides. In coverage-driven mode, a ONE/ONE
     * result is authoritative for this plot grid, so avoid the old CSG
     * fallback. */
    if (cell_plus == cell_minus && g->coverage)
        return 0;

    /* Same cell on both sides — this curve might be inside a nested overlap
     * (e.g., inner sphere fully inside outer sphere). The grid only stores
     * the "winning" cell so it can't see the overlap. Do a single CSG query
     * on one side to check if multiple cells claim this region. */
    if (g->sys && g->plane) {
        const alea_slice_plane_t* pl = g->plane;
        double x = pl->origin[0] + u_plus * pl->u_axis[0] + v_plus * pl->v_axis[0];
        double y = pl->origin[1] + u_plus * pl->u_axis[1] + v_plus * pl->v_axis[1];
        double z = pl->origin[2] + u_plus * pl->u_axis[2] + v_plus * pl->v_axis[2];

        int count = count_cells_in_universe_at_point(g->sys, universe_id, x, y, z);
        if (count > 1) return ALEA_SLICE_ERR_OVERLAP;
    }

    return 0;
}

/** Grid-based version of check_parametric_curve */
static void check_parametric_curve_grid(const grid_ctx_t* g,
                                        const alea_slice_view_t* view,
                                        const alea_curve_2d_t* curve,
                                        size_t curve_index,
                                        double sample_spacing,
                                        double eps,
                                        error_vec_t* errors) {
    double t_lo, t_hi;
    get_clipped_param_range(curve, view, &t_lo, &t_hi);
    if (t_lo >= t_hi) return;

    double arc_approx = t_hi - t_lo;
    if (curve->type == ALEA_CURVE_CIRCLE || curve->type == ALEA_CURVE_ARC) {
        arc_approx = (t_hi - t_lo) * curve->data.circle.radius;
    } else if (curve->type == ALEA_CURVE_ELLIPSE || curve->type == ALEA_CURVE_ELLIPSE_ARC) {
        double avg_r = (curve->data.ellipse.semi_a + curve->data.ellipse.semi_b) * 0.5;
        arc_approx = (t_hi - t_lo) * avg_r;
    }

    int n_samples = (int)(arc_approx / sample_spacing);
    if (n_samples < 2) n_samples = 2;
    if (n_samples > 2000) n_samples = 2000;

    double dt_finite = (t_hi - t_lo) * 1e-6;
    if (dt_finite < 1e-15) dt_finite = 1e-15;

    int cur_type = 0;
    double seg_t_start = 0;

    for (int i = 0; i <= n_samples; i++) {
        double t = t_lo + (t_hi - t_lo) * i / n_samples;

        double u, v;
        if (!alea_curve_eval(curve, t, &u, &v)) continue;
        if (!point_in_viewport(u, v, view)) continue;

        double u1, v1, u2, v2;
        double t_back = t - dt_finite;
        double t_fwd  = t + dt_finite;
        if (t_back < t_lo) t_back = t_lo;
        if (t_fwd > t_hi) t_fwd = t_hi;

        if (!alea_curve_eval(curve, t_back, &u1, &v1) ||
            !alea_curve_eval(curve, t_fwd,  &u2, &v2))
            continue;

        double tang_u = u2 - u1;
        double tang_v = v2 - v1;
        double tang_len = sqrt(tang_u * tang_u + tang_v * tang_v);
        if (tang_len < 1e-20) continue;

        double nu = -tang_v / tang_len;
        double nv =  tang_u / tang_len;

        int classification = classify_sample_grid(g, u, v, nu, nv, eps, curve->universe_id);

        if (classification != cur_type) {
            if (cur_type != 0) {
                alea_slice_error_t err = {
                    .curve_index = curve_index,
                    .surface_id = curve->surface_id,
                    .type = (alea_slice_error_type_t)cur_type,
                    .t_start = seg_t_start,
                    .t_end = t
                };
                error_vec_push(errors, &err);
            }
            cur_type = classification;
            seg_t_start = t;
        }
    }

    if (cur_type != 0) {
        alea_slice_error_t err = {
            .curve_index = curve_index,
            .surface_id = curve->surface_id,
            .type = (alea_slice_error_type_t)cur_type,
            .t_start = seg_t_start,
            .t_end = t_hi
        };
        error_vec_push(errors, &err);
    }
}

/** Grid-based version of check_quartic_curve */
static void check_quartic_curve_grid(const grid_ctx_t* g,
                                     const alea_slice_view_t* view,
                                     const alea_curve_2d_t* curve,
                                     size_t curve_index,
                                     double sample_spacing,
                                     double eps,
                                     error_vec_t* errors) {
    double bbox_umin, bbox_umax, bbox_vmin, bbox_vmax;
    alea_curve_bbox(curve, &bbox_umin, &bbox_umax, &bbox_vmin, &bbox_vmax);

    double v_lo = fmax(bbox_vmin, view->v_min);
    double v_hi = fmin(bbox_vmax, view->v_max);
    if (v_lo >= v_hi) return;

    int n_scanlines = (int)((v_hi - v_lo) / sample_spacing);
    if (n_scanlines < 2) n_scanlines = 2;
    if (n_scanlines > 2000) n_scanlines = 2000;

    for (int i = 0; i <= n_scanlines; i++) {
        double v_scan = v_lo + (v_hi - v_lo) * i / n_scanlines;

        double u_vals[16];
        int n_isect = alea_curve_scanline_intersect(curve, v_scan, u_vals, 16);

        for (int j = 0; j < n_isect; j++) {
            double u = u_vals[j];
            if (u < view->u_min || u > view->u_max) continue;

            double u_vals2[16];
            double dv = sample_spacing * 0.01;
            int n2 = alea_curve_scanline_intersect(curve, v_scan + dv, u_vals2, 16);

            double tang_u = 0, tang_v = dv;
            double best_dist = 1e20;
            for (int k = 0; k < n2; k++) {
                double dist = fabs(u_vals2[k] - u);
                if (dist < best_dist) {
                    best_dist = dist;
                    tang_u = u_vals2[k] - u;
                }
            }

            double tang_len = sqrt(tang_u * tang_u + tang_v * tang_v);
            if (tang_len < 1e-20) { tang_u = 0; tang_v = 1; tang_len = 1; }

            double nu = -tang_v / tang_len;
            double nv =  tang_u / tang_len;

            int classification = classify_sample_grid(g, u, v_scan, nu, nv, eps, curve->universe_id);

            if (classification != 0) {
                alea_slice_error_t err = {
                    .curve_index = curve_index,
                    .surface_id = curve->surface_id,
                    .type = (alea_slice_error_type_t)classification,
                    .t_start = v_scan,
                    .t_end = v_scan
                };
                error_vec_push(errors, &err);
            }
        }
    }
}

/** Grid-based version of check_polygon_curve */
static void check_polygon_curve_grid(const grid_ctx_t* g,
                                     const alea_slice_view_t* view,
                                     const alea_curve_2d_t* curve,
                                     size_t curve_index,
                                     double sample_spacing,
                                     double eps,
                                     error_vec_t* errors) {
    const alea_polygon_2d_t* poly = &curve->data.polygon;
    int nverts = poly->vertex_count;
    if (nverts < 2) return;

    int n_edges = poly->closed ? nverts : (nverts - 1);

    for (int e = 0; e < n_edges; e++) {
        int i0 = e;
        int i1 = (e + 1) % nverts;
        double x0 = poly->vertices[i0][0], y0 = poly->vertices[i0][1];
        double x1 = poly->vertices[i1][0], y1 = poly->vertices[i1][1];

        double edge_len = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
        if (edge_len < 1e-15) continue;

        int n_samples = (int)(edge_len / sample_spacing);
        if (n_samples < 2) n_samples = 2;
        if (n_samples > 500) n_samples = 500;

        double dx = x1 - x0, dy = y1 - y0;
        double nu = -dy / edge_len;
        double nv_dir = dx / edge_len;

        int cur_type = 0;
        double seg_t_start = 0;

        for (int s = 0; s <= n_samples; s++) {
            double frac = (double)s / n_samples;
            double u = x0 + frac * dx;
            double v = y0 + frac * dy;
            double t_param = i0 + frac;

            if (!point_in_viewport(u, v, view)) continue;

            int classification = classify_sample_grid(g, u, v,
                                                      nu, nv_dir, eps,
                                                      curve->universe_id);

            if (classification != cur_type) {
                if (cur_type != 0) {
                    alea_slice_error_t err = {
                        .curve_index = curve_index,
                        .surface_id = curve->surface_id,
                        .type = (alea_slice_error_type_t)cur_type,
                        .t_start = seg_t_start,
                        .t_end = t_param
                    };
                    error_vec_push(errors, &err);
                }
                cur_type = classification;
                seg_t_start = t_param;
            }
        }

        if (cur_type != 0) {
            alea_slice_error_t err = {
                .curve_index = curve_index,
                .surface_id = curve->surface_id,
                .type = (alea_slice_error_type_t)cur_type,
                .t_start = seg_t_start,
                .t_end = (double)(i0 + 1)
            };
            error_vec_push(errors, &err);
        }
    }
}

static alea_slice_error_result_t* check_slice_errors_grid_impl(
    alea_system_t* sys,
    const alea_slice_view_t* view,
    const alea_slice_curves_t* curves,
    const int* cell_ids,
    const uint8_t* coverage,
    const uint8_t* grid_errors,
    int nu, int nv)
{
    if (!view || !curves || !cell_ids || nu <= 0 || nv <= 0) return NULL;

    const alea_curve_collection_t* coll = &curves->internal;
    if (coll->curves.count == 0) return NULL;

    /* Set up grid context */
    double u_range = view->u_max - view->u_min;
    double v_range = view->v_max - view->v_min;

    grid_ctx_t g = {
        .cell_ids    = cell_ids,
        .coverage    = coverage,
        .grid_errors = grid_errors,
        .nu = nu,
        .nv = nv,
        .u_min = view->u_min,
        .v_min = view->v_min,
        .inv_du = nu / u_range,
        .inv_dv = nv / v_range,
        .sys   = sys,
        .plane = &view->plane,
    };

    /* Offset must be at least 1.5 pixels so the two lookups land in different
     * pixels even when the normal is diagonal */
    double pixel_size = fmin(u_range / nu, v_range / nv);
    double eps = pixel_size * 1.5;

    double vp_diag = sqrt(u_range * u_range + v_range * v_range);
    double sample_spacing = vp_diag / 200.0;

    error_vec_t errvec = ALEA_VEC_INIT;

    for (size_t i = 0; i < coll->curves.count; i++) {
        const alea_curve_2d_t* curve = &coll->curves.data[i];

        switch (curve->type) {
            case ALEA_CURVE_NONE:
            case ALEA_CURVE_POINT:
                continue;
            case ALEA_CURVE_QUARTIC:
                check_quartic_curve_grid(&g, view, curve, i,
                                         sample_spacing, eps, &errvec);
                continue;
            case ALEA_CURVE_POLYGON:
                check_polygon_curve_grid(&g, view, curve, i,
                                         sample_spacing, eps, &errvec);
                continue;
            default:
                check_parametric_curve_grid(&g, view, curve, i,
                                            sample_spacing, eps, &errvec);
                continue;
        }
    }

    alea_slice_error_result_t* result = malloc(sizeof(alea_slice_error_result_t));
    if (!result) {
        free(errvec.data);
        return NULL;
    }

    result->errors = errvec.data;
    result->error_count = errvec.count;
    return result;
}

alea_slice_error_result_t* alea_check_slice_errors_grid(
    alea_system_t* sys,
    const alea_slice_view_t* view,
    const alea_slice_curves_t* curves,
    const int* cell_ids,
    const uint8_t* grid_errors,
    int nu, int nv)
{
    return check_slice_errors_grid_impl(sys, view, curves, cell_ids, NULL,
                                        grid_errors, nu, nv);
}

alea_slice_error_result_t* alea_check_slice_errors_grid_ex(
    const alea_slice_view_t* view,
    const alea_slice_curves_t* curves,
    const int* cell_ids,
    const uint8_t* coverage,
    const uint8_t* grid_errors,
    int nu, int nv)
{
    return check_slice_errors_grid_impl(NULL, view, curves, cell_ids, coverage,
                                        grid_errors, nu, nv);
}

void alea_slice_errors_free(alea_slice_error_result_t* result) {
    if (!result) return;
    free(result->errors);
    free(result);
}

static int component_pixel_matches(const uint8_t* coverage,
                                   uint8_t target,
                                   int idx) {
    return coverage[idx] == target;
}

alea_plot_error_component_result_t* alea_classify_plot_error_components(
    const int* cell_ids,
    const int* secondary_cell_ids,
    const uint8_t* coverage,
    int nu, int nv) {
    if (!cell_ids || !coverage || nu <= 0 || nv <= 0) return NULL;

    size_t n = (size_t)nu * (size_t)nv;
    uint8_t* visited = calloc(n, 1);
    int* queue = malloc(n * sizeof(*queue));
    if (!visited || !queue) {
        free(visited);
        free(queue);
        return NULL;
    }

    plot_error_component_vec_t comps = ALEA_VEC_INIT;
    const int di[4] = { 1, -1, 0, 0 };
    const int dj[4] = { 0, 0, 1, -1 };

    for (int j0 = 0; j0 < nv; j0++) {
        for (int i0 = 0; i0 < nu; i0++) {
            int start = j0 * nu + i0;
            uint8_t cov = coverage[start];
            if (visited[start]) continue;
            if (cov != ALEA_COVERAGE_NONE && cov != ALEA_COVERAGE_MULTI)
                continue;

            int head = 0, tail = 0;
            visited[start] = 1;
            queue[tail++] = start;

            alea_plot_error_component_t comp;
            memset(&comp, 0, sizeof(comp));
            comp.kind = (cov == ALEA_COVERAGE_NONE)
                ? ALEA_PLOT_ERR_UNDEFINED_REGION
                : ALEA_PLOT_ERR_TOTAL_OVERLAP;
            comp.primary_cell_id = -1;
            comp.secondary_cell_id = -1;
            comp.min_i = comp.max_i = i0;
            comp.min_j = comp.max_j = j0;

            bool primary_visible = false;
            bool secondary_visible = false;
            bool other_visible = false;

            while (head < tail) {
                int idx = queue[head++];
                int i = idx % nu;
                int j = idx / nu;
                comp.pixel_count++;
                if (i < comp.min_i) comp.min_i = i;
                if (i > comp.max_i) comp.max_i = i;
                if (j < comp.min_j) comp.min_j = j;
                if (j > comp.max_j) comp.max_j = j;

                if (cov == ALEA_COVERAGE_MULTI) {
                    if (comp.primary_cell_id < 0 && cell_ids[idx] >= 0)
                        comp.primary_cell_id = cell_ids[idx];
                    if (secondary_cell_ids &&
                        comp.secondary_cell_id < 0 &&
                        secondary_cell_ids[idx] >= 0) {
                        comp.secondary_cell_id = secondary_cell_ids[idx];
                    }
                }

                for (int k = 0; k < 4; k++) {
                    int ni = i + di[k];
                    int nj = j + dj[k];
                    if (ni < 0 || ni >= nu || nj < 0 || nj >= nv)
                        continue;
                    int nidx = nj * nu + ni;
                    if (component_pixel_matches(coverage, cov, nidx)) {
                        if (!visited[nidx]) {
                            visited[nidx] = 1;
                            queue[tail++] = nidx;
                        }
                        continue;
                    }

                    if (cov != ALEA_COVERAGE_MULTI)
                        continue;

                    int neighbor_cell = cell_ids[nidx];
                    if (neighbor_cell < 0 ||
                        coverage[nidx] == ALEA_COVERAGE_MULTI)
                        continue;
                    if (neighbor_cell == comp.primary_cell_id) {
                        primary_visible = true;
                    } else if (comp.secondary_cell_id >= 0 &&
                               neighbor_cell == comp.secondary_cell_id) {
                        secondary_visible = true;
                    } else {
                        other_visible = true;
                    }
                }
            }

            if (cov == ALEA_COVERAGE_MULTI) {
                primary_visible = false;
                secondary_visible = false;
                other_visible = false;
                for (int q = 0; q < tail; q++) {
                    int idx = queue[q];
                    int i = idx % nu;
                    int j = idx / nu;
                    for (int k = 0; k < 4; k++) {
                        int ni = i + di[k];
                        int nj = j + dj[k];
                        if (ni < 0 || ni >= nu || nj < 0 || nj >= nv)
                            continue;
                        int nidx = nj * nu + ni;
                        if (coverage[nidx] == ALEA_COVERAGE_MULTI)
                            continue;
                        int neighbor_cell = cell_ids[nidx];
                        if (neighbor_cell < 0) continue;
                        if (neighbor_cell == comp.primary_cell_id) {
                            primary_visible = true;
                        } else if (comp.secondary_cell_id >= 0 &&
                                   neighbor_cell == comp.secondary_cell_id) {
                            secondary_visible = true;
                        } else {
                            other_visible = true;
                        }
                    }
                }

                if (secondary_visible ||
                    (comp.secondary_cell_id < 0 && other_visible) ||
                    (primary_visible && other_visible)) {
                    comp.kind = ALEA_PLOT_ERR_PARTIAL_OVERLAP;
                } else {
                    comp.kind = ALEA_PLOT_ERR_TOTAL_OVERLAP;
                }
            }

            plot_error_component_vec_push(&comps, &comp);
        }
    }

    free(visited);
    free(queue);

    alea_plot_error_component_result_t* result =
        malloc(sizeof(*result));
    if (!result) {
        free(comps.data);
        return NULL;
    }
    result->components = comps.data;
    result->component_count = comps.count;
    return result;
}

void alea_plot_error_components_free(
    alea_plot_error_component_result_t* result) {
    if (!result) return;
    free(result->components);
    free(result);
}
