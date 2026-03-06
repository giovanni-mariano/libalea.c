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
#include "core/alea_eval.h"
#include "core/alea_universe.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util/math.h"
#include "util/alea_log.h"

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

alea_slice_curves_t* alea_get_slice_curves(const alea_system_t* sys,
                                                    const alea_slice_view_t* view) {
    if (!sys || !view) return NULL;

    alea_slice_curves_t* result = calloc(1, sizeof(alea_slice_curves_t));
    if (!result) return NULL;

    if (!sys->spatial_index) {
        if (alea_build_spatial_index((alea_system_t*)sys) != 0) {
            free(result);
            return NULL;
        }
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
    return curves ? curves->internal.curve_count : 0;
}

int alea_slice_curves_get(const alea_slice_curves_t* curves, size_t index, alea_curve_t* out) {
    if (!curves || !out || index >= curves->internal.curve_count) return -1;
    const alea_curve_2d_t* src = &curves->internal.curves[index];
    memset(out, 0, sizeof(*out));

    out->surface_id = src->surface_id;
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
static void multilevel_hint_from_hits(const alea_system_t* sys,
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
                    alea_matrix_from_mcnp(&fill_mat, tr->cosines, tr->value_count, false);
                    alea_matrix_t new_acc;
                    alea_matrix_multiply(&new_acc, &accumulated, &fill_mat);
                    accumulated = new_acc;
                }
            }
        }

        hint->valid_depth = i;
    }
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
static int find_cell_in_universe_with_hint(const alea_system_t* sys,
                                            double lx, double ly, double lz,
                                            int universe_id,
                                            int hint_cell_idx) {
    /* Try hint cell first (coherence) */
    if (hint_cell_idx >= 0 && (size_t)hint_cell_idx < alea_vec_count(&sys->cells)) {
        const alea_cell_entry_t* cell = &sys->cells.data[hint_cell_idx];

        /* Verify hint cell is in the correct universe */
        if (cell->universe_id == universe_id &&
            cell->root_node_id != ALEA_NODE_ID_INVALID &&
            alea_contains_point(sys, cell->root_node_id, lx, ly, lz)) {
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
                    return neighbor_idx;
                }
            }
        }
    }

    /* Fall back to linear search within universe */
    const alea_universe_t* univ = alea_get_universe(sys, universe_id);
    if (!univ) return -1;

    for (size_t i = 0; i < univ->cell_count; i++) {
        size_t cell_idx = univ->cell_indices[i];
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
static void find_cell_multilevel(const alea_system_t* sys,
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
                /* Reached target depth */
                *out_cell_id = cell->mc_cell_id;
                if (out_material_id) *out_material_id = cell->material_id;
                if (out_hint) {
                    multilevel_hint_from_hits(sys, local_hits, hits_found, out_hint);
                }
                return;
            }

            /* If cell has a fill, continue to next level */
            if (cell->fill_universe > 0) {
                /* Build transform for the fill */
                alea_matrix_t fill_transform;
                if (cell->fill_transform > 0) {
                    const alea_transform_t* tr = alea_get_transform(sys, cell->fill_transform);
                    if (tr) {
                        /* Use tr->cosines which has pre-computed direction cosines */
                        alea_matrix_from_mcnp(&fill_transform, tr->cosines,
                                            tr->value_count, false);
                    } else {
                        alea_matrix_identity(&fill_transform);
                    }
                } else {
                    alea_matrix_identity(&fill_transform);
                }

                /* Compose with accumulated transform */
                alea_matrix_t new_accumulated;
                alea_matrix_multiply(&new_accumulated, &accumulated, &fill_transform);
                accumulated = new_accumulated;

                /* Transform point to fill universe coordinates */
                /* Use a copy for inversion to preserve accumulated for next iteration */
                alea_matrix_t inv_accumulated = accumulated;
                alea_matrix_invert(&inv_accumulated);
                lx = gx; ly = gy; lz = gz;
                alea_matrix_transform_point_inverse(&inv_accumulated, &lx, &ly, &lz);

                current_universe = cell->fill_universe;
            } else {
                /* Terminal cell (no fill) - this is the innermost */
                if (universe_depth < 0) {
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
                /* This is a terminal cell */
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

int alea_find_cells_grid(const alea_system_t* sys,
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

    /* Build cell adjacency if not already built - needed for ALL depths now */
    if (!sys->cell_adjacency_built) {
        alea_build_cell_adjacency((alea_system_t*)sys);
    }

    /* Build spatial index eagerly before the parallel region so that:
     * (a) the internal #pragma omp parallel for gets full thread parallelism
     *     (no nesting penalty)
     * (b) OpenMP workers don't race into the lazy-build path */
    alea_spatial_index_build((alea_system_t*)sys);

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
    #pragma omp parallel for schedule(dynamic, 4)
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

            int cell_id = -1;
            int material_id = 0;
            uint8_t error = 0;
            alea_multilevel_hint_t curr_hint;

            find_cell_multilevel(sys, x, y, z, universe_depth,
                                 &row_hint, &curr_hint,
                                 &cell_id, &material_id, &error);

            row_hint = curr_hint;  /* Use current as hint for next pixel */

            out_cell_ids[idx] = cell_id;
            if (out_material_ids) out_material_ids[idx] = material_id;
            if (out_errors) out_errors[idx] = error;
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
    }

    return 0;
}


/* ============================================================================
 * FULL-GRID OVERLAP CHECK
 * ============================================================================ */

int alea_check_grid_overlaps(const alea_system_t* sys,
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

/* ============================================================================
 * CURVE-GUIDED OVERLAP DETECTION
 * ============================================================================ */

/**
 * Check a single pixel for nested overlap via full hierarchy query.
 * Returns 1 if overlap found and error grid updated, 0 otherwise.
 */
static int probe_pixel_for_overlap(const alea_system_t* sys,
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

int alea_check_grid_overlaps_curves(const alea_system_t* sys,
                                    const alea_slice_view_t* view,
                                    const alea_slice_curves_t* curves,
                                    int nu, int nv,
                                    int universe_depth,
                                    const int* cell_ids,
                                    uint8_t* errors) {
    if (!sys || !view || !curves || !cell_ids || !errors || nu <= 0 || nv <= 0)
        return -1;

    const alea_curve_collection_t* coll = &curves->internal;
    if (coll->curve_count == 0) return 0;

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

    for (size_t ci = 0; ci < coll->curve_count; ci++) {
        const alea_curve_2d_t* curve = &coll->curves[ci];

        if (curve->type == ALEA_CURVE_NONE || curve->type == ALEA_CURVE_POINT)
            continue;

        /* For parametric curves: walk the parameter range and rasterize */
        if (curve->type == ALEA_CURVE_QUARTIC || curve->type == ALEA_CURVE_POLYGON)
            continue;  /* TODO: handle these if needed */

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

            int pi = (int)((u - u_min) * inv_du);
            int pj = (int)((v - v_min) * inv_dv);
            if (pi < 0 || pi >= nu || pj < 0 || pj >= nv) continue;

            /* Skip boundary pixels — already checked by pass 2 */
            if (pixel_is_boundary(cell_ids, nu, nv, pi, pj)) continue;

            /* Probe this pixel and its immediate neighbors (to catch
             * near-boundary cases where the curve is between pixels) */
            for (int dj = -1; dj <= 1; dj++) {
                for (int di = -1; di <= 1; di++) {
                    int qi = pi + di, qj = pj + dj;
                    if (qi < 0 || qi >= nu || qj < 0 || qj >= nv) continue;

                    size_t bit_idx = (size_t)qj * nu + qi;
                    size_t byte_idx = bit_idx / 8;
                    uint8_t bit_mask = 1u << (bit_idx % 8);
                    if (probed[byte_idx] & bit_mask) continue;
                    probed[byte_idx] |= bit_mask;

                    /* Skip if already a boundary pixel */
                    if (pixel_is_boundary(cell_ids, nu, nv, qi, qj)) continue;

                    new_overlaps += probe_pixel_for_overlap(
                        sys, plane, u_min, v_min, du, dv,
                        nu, nv, universe_depth, cell_ids, errors, qi, qj);
                }
            }
        }
    }

    free(probed);
    return new_overlaps;
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

int alea_find_surface_label_positions(
    const alea_slice_curves_t* curves,
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

        /* Compute position at curve midpoint */
        double t_mid;
        if (curve.type == ALEA_CURVE_CIRCLE) {
            /* Full circle: pick angle π/4 (45°) for label position */
            t_mid = M_PI / 4.0;
        } else if (curve.type == ALEA_CURVE_ELLIPSE) {
            /* Full ellipse: pick angle π/4 */
            t_mid = M_PI / 4.0;
        } else if (curve.t_min == curve.t_max) {
            /* Degenerate bounds - use 0 */
            t_mid = 0;
        } else {
            t_mid = (curve.t_min + curve.t_max) / 2.0;
        }

        double px, py;
        if (eval_curve_at_t(&curve, t_mid, &px, &py) != 0) continue;

        /* Convert to pixel coordinates */
        int ix = (int)((px - x_min) / dx);
        int iy = (int)((py - y_min) / dy);

        /* Check margin */
        if (ix < margin || ix >= width - margin ||
            iy < margin || iy >= height - margin) {
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
static int count_cells_in_universe_at_point(const alea_system_t* sys,
                                             int universe_id,
                                             double x, double y, double z) {
    alea_cell_hit_t hits[16];
    int n = alea_spatial_find_cells_at_point(sys, x, y, z, hits, 16);
    if (n <= 0) return 0;

    int count = 0;
    for (int i = 0; i < n; i++) {
        if (hits[i].universe_id == universe_id) {
            count++;
            if (count > 1) return count;  /* Early exit on overlap */
        }
    }
    return count;
}

static int classify_sample(const alea_system_t* sys,
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

/**
 * Check if a 2D point is inside the viewport (with small margin).
 */
static bool point_in_viewport(double u, double v, const alea_slice_view_t* view) {
    return u >= view->u_min && u <= view->u_max &&
           v >= view->v_min && v <= view->v_max;
}

/* Dynamic array for error segments */
typedef struct {
    alea_slice_error_t* data;
    size_t count;
    size_t capacity;
} error_vec_t;

static void error_vec_push(error_vec_t* vec, const alea_slice_error_t* err) {
    if (vec->count >= vec->capacity) {
        vec->capacity = vec->capacity ? vec->capacity * 2 : 32;
        vec->data = realloc(vec->data, vec->capacity * sizeof(alea_slice_error_t));
    }
    vec->data[vec->count++] = *err;
}

/**
 * Process a single parametric curve: sample, classify, merge error segments.
 */
static void check_parametric_curve(const alea_system_t* sys,
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
static void check_quartic_curve(const alea_system_t* sys,
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
static void check_polygon_curve(const alea_system_t* sys,
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
    const alea_system_t* sys,
    const alea_slice_view_t* view,
    const alea_slice_curves_t* curves,
    int universe_depth)
{
    (void)universe_depth;  /* universe_id is now stored per-curve */
    if (!sys || !view || !curves) return NULL;

    const alea_curve_collection_t* coll = &curves->internal;
    if (coll->curve_count == 0) return NULL;

    /* Reset spatial coherence cache so it starts fresh for curve walking */
    alea_spatial_reset_cache();

    /* Compute viewport-derived parameters */
    double vp_width  = view->u_max - view->u_min;
    double vp_height = view->v_max - view->v_min;
    double vp_diag = sqrt(vp_width * vp_width + vp_height * vp_height);

    double sample_spacing = vp_diag / 200.0;
    double eps = vp_diag * 1e-6;

    error_vec_t errvec = { NULL, 0, 0 };

    for (size_t i = 0; i < coll->curve_count; i++) {
        const alea_curve_2d_t* curve = &coll->curves[i];

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
    const uint8_t* grid_errors;
    int nu, nv;
    double u_min, v_min;
    double inv_du, inv_dv;  /* 1/du, 1/dv for fast coord→pixel */
    /* For CSG fallback on same-cell cases (nested overlaps) */
    const alea_system_t* sys;
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

    int cell_plus  = g->cell_ids[idx_plus];
    int cell_minus = g->cell_ids[idx_minus];

    /* Different cells on each side → proper boundary, OK */
    if (cell_plus != cell_minus && cell_plus >= 0 && cell_minus >= 0)
        return 0;

    /* One or both sides void → gap */
    if (cell_plus < 0 || cell_minus < 0)
        return ALEA_SLICE_ERR_GAP;

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

alea_slice_error_result_t* alea_check_slice_errors_grid(
    const alea_system_t* sys,
    const alea_slice_view_t* view,
    const alea_slice_curves_t* curves,
    const int* cell_ids,
    const uint8_t* grid_errors,
    int nu, int nv)
{
    if (!view || !curves || !cell_ids || nu <= 0 || nv <= 0) return NULL;

    const alea_curve_collection_t* coll = &curves->internal;
    if (coll->curve_count == 0) return NULL;

    /* Set up grid context */
    double u_range = view->u_max - view->u_min;
    double v_range = view->v_max - view->v_min;

    grid_ctx_t g = {
        .cell_ids    = cell_ids,
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

    error_vec_t errvec = { NULL, 0, 0 };

    for (size_t i = 0; i < coll->curve_count; i++) {
        const alea_curve_2d_t* curve = &coll->curves[i];

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

void alea_slice_errors_free(alea_slice_error_result_t* result) {
    if (!result) return;
    free(result->errors);
    free(result);
}

