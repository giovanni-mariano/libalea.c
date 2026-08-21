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
#include "raycast/raycast.h"
#include "raycast/ray_epsilon.h"
#include "core/alea_system.h"
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

/* Grid query diagnostics are opt-in. Each OpenMP worker updates a private
 * instance and merges it once after the grid loop, avoiding shared cache-line
 * traffic in the per-pixel hot path. A NULL pointer disables collection. */
typedef struct {
    size_t px_total;
    size_t rich_root_queries;
    size_t rich_path_reused;
    size_t rich_lattice_transitions;
    size_t rich_prefix_restarts;
    size_t rich_full_fallbacks;
    size_t csg_prim_evals;    /* total primitive evaluations */
    size_t csg_bool_ops;      /* total boolean ops in CSG tree */
} grid_query_stats_t;

#define GRID_STAT_INC(stats, field) \
    do { if ((stats) != NULL) (stats)->field++; } while (0)

#ifdef _OPENMP
static void grid_query_stats_merge(grid_query_stats_t* dst,
                                   const grid_query_stats_t* src) {
    dst->px_total += src->px_total;
    dst->rich_root_queries += src->rich_root_queries;
    dst->rich_path_reused += src->rich_path_reused;
    dst->rich_lattice_transitions += src->rich_lattice_transitions;
    dst->rich_prefix_restarts += src->rich_prefix_restarts;
    dst->rich_full_fallbacks += src->rich_full_fallbacks;
    dst->csg_prim_evals += src->csg_prim_evals;
    dst->csg_bool_ops += src->csg_bool_ops;
}
#endif

static alea_tile_coverage_stats_t g_tile_coverage_stats;
static alea_point_coverage_stats_t g_point_coverage_stats;
static alea_sparse_surface_label_stats_t g_sparse_surface_label_stats;

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

alea_sparse_surface_label_stats_t alea_sparse_surface_label_stats_get(void) {
    return g_sparse_surface_label_stats;
}

/* ============================================================================
 * SLICE CURVES API
 * ============================================================================ */

struct alea_slice_curves {
    alea_curve_collection_t internal;
    double view_u_min, view_u_max;
    double view_v_min, view_v_max;
};

struct alea_slice_surface_boundary_map {
    int width;
    int height;
    alea_slice_boundary_status_t* status; /* [2 * width * height] */
    size_t* surface_offsets;              /* [2 * width * height + 1] */
    int* surface_ids;                     /* CSR participant IDs */
    size_t surface_count;
    size_t surface_capacity;
    size_t* group_offsets;                /* [2 * width * height + 1] */
    double* group_fractions;              /* CSR group positions */
    size_t* group_surface_offsets;        /* [group_count + 1] */
    int* group_surface_ids;               /* CSR participant IDs per group */
    size_t group_count;
    size_t group_surface_count;
    size_t group_capacity;
    size_t group_surface_capacity;
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

        /* The public structure intentionally does not expose the internal
         * conic/torus coefficients, but consumers still need the curve family
         * for sampled rendering and diagnostics. */
        case ALEA_CURVE_PARABOLA:
        case ALEA_CURVE_HYPERBOLA:
        case ALEA_CURVE_QUARTIC:
            out->type = src->type;
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
#define GRID_ERR_NONE           0
#define GRID_ERR_OVERLAP        1
#define GRID_ERR_UNDEFINED      2
/* Deepest containing cell is a fill/lattice container whose filling universe
 * has no cell at the point (MCNP undefined region): the pixel shows the
 * container, flagged. */
#define GRID_ERR_UNDEFINED_FILL 3

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
static int find_point_coverage_exact_uv(alea_system_t* sys,
                                        const alea_slice_plane_t* plane,
                                        double u, double v,
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

/* Select the same target depth policy as the grid's primary cell output,
 * but take its concrete placement transform from the canonical rich path.
 * In particular this keeps occurrence identity for repeated cells and lattice
 * elements; a slice-local hint cannot reconstruct that information safely. */
static int slice_project_hier_path(const alea_hier_ray_path_t* path,
                                   int universe_depth,
                                   const alea_hier_ray_path_entry_t** out_entry) {
    if (!path || !out_entry || path->count <= 0) return -1;

    int target = -1;
    if (universe_depth < 0) {
        target = path->count - 1;
    } else {
        for (int i = 0; i < path->count; i++) {
            if (path->entries[i].depth == universe_depth) {
                target = i;
                break;
            }
        }
        if (target < 0) {
            for (int i = path->count - 1; i >= 0; i--) {
                if (path->entries[i].depth <= universe_depth) {
                    target = i;
                    break;
                }
            }
        }
    }

    if (target < 0) return -1;
    *out_entry = &path->entries[target];
    return 0;
}

/* Primary grid projection from a canonical rich path. Overlap reporting stays
 * separate because a first-owner path intentionally contains one branch only. */
static int slice_resolve_grid_rich(alea_system_t* sys,
                                   const alea_hier_coherence_state_t* previous,
                                   double x, double y, double z,
                                   int universe_depth,
                                   alea_hier_coherence_state_t* current,
                                   int* out_cell_id,
                                   int* out_material_id,
                                   uint8_t* out_error,
                                   grid_query_stats_t* stats) {
    if (!out_cell_id || !out_material_id || !out_error) return -1;
    *out_cell_id = -1;
    *out_material_id = 0;
    *out_error = GRID_ERR_NONE;

    alea_hier_cell_hit_t deepest;
    alea_hier_coherence_kind_t kind;
    /* Coherent ownership: re-deriving the deck-first owner per point would only
     * change the cell reported inside an overlap, and those pixels are already
     * flagged by the boundary/periodic passes below, which query overlaps
     * directly instead of through the cached path. */
    int rc = alea_hier_spatial_resolve_coherent(
        sys, x, y, z, previous, ALEA_HIER_COH_OWNERSHIP_COHERENT,
        current, &deepest, &kind);
    GRID_STAT_INC(stats, px_total);
    if (rc > 0) {
        switch (kind) {
        case ALEA_HIER_COH_ROOT_QUERY:
            GRID_STAT_INC(stats, rich_root_queries);
            break;
        case ALEA_HIER_COH_PATH_REUSED:
            GRID_STAT_INC(stats, rich_path_reused);
            break;
        case ALEA_HIER_COH_LATTICE_TRANSITION:
            GRID_STAT_INC(stats, rich_lattice_transitions);
            break;
        case ALEA_HIER_COH_PREFIX_RESTART:
            GRID_STAT_INC(stats, rich_prefix_restarts);
            break;
        case ALEA_HIER_COH_FULL_FALLBACK:
            GRID_STAT_INC(stats, rich_full_fallbacks);
            break;
        }
    }
    if (rc < 0) return -1;
    if (rc == 0) {
        *out_error = GRID_ERR_UNDEFINED;
        return 0;
    }

    const alea_hier_ray_path_entry_t* target = NULL;
    if (slice_project_hier_path(&current->path, universe_depth, &target) != 0)
        return -1;
    *out_cell_id = target->cell_id;
    *out_material_id = target->material_id;
    if (universe_depth < 0 &&
        (deepest.hit.resolution_flags & ALEA_RESOLVE_UNDEFINED_FILL)) {
        *out_error = GRID_ERR_UNDEFINED_FILL;
    }
    return 1;
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
    if (alea_interrupted()) return -1;

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

    const bool stats_en = getenv("ALEA_GRID_STATS") != NULL;
    grid_query_stats_t grid_stats = {0};

#ifdef _OPENMP
    /* ALEA_GRID_VERIFY_INTERVAL=N: every N-th pixel also runs a full recursive
     * search to detect overlapping geometry. 0 or unset = no periodic scan
     * (relies on boundary-pixel second pass only). Useful for finding geometry
     * errors in interior regions not adjacent to cell boundaries. */
    const char* vi_env = getenv("ALEA_GRID_VERIFY_INTERVAL");
    int overlap_interval = (vi_env && atoi(vi_env) > 0) ? atoi(vi_env) : 0;
    #pragma omp parallel
    {
        grid_query_stats_t thread_stats = {0};
        grid_query_stats_t* active_stats = stats_en ? &thread_stats : NULL;
        if (stats_en) alea_perf_reset();
        #pragma omp for schedule(dynamic, 4)
        for (int j = 0; j < nv; j++) {
            if (alea_interrupted()) continue;
            double v = v_min + (j + 0.5) * dv;
            alea_hier_coherence_state_t state_a;
            alea_hier_coherence_state_t state_b;
            alea_hier_coherence_state_t* previous_state = &state_a;
            alea_hier_coherence_state_t* current_state = &state_b;
            alea_hier_coherence_state_clear(previous_state);
            alea_hier_coherence_state_clear(current_state);

            for (int i = 0; i < nu; i++) {
                if (alea_interrupted()) break;
                double u = u_min + (i + 0.5) * du;
                int idx = j * nu + i;

                double x = origin[0] + u * u_axis[0] + v * v_axis[0];
                double y = origin[1] + u * u_axis[1] + v * v_axis[1];
                double z = origin[2] + u * u_axis[2] + v * v_axis[2];

                int cell_id = -1, material_id = 0;
                uint8_t error = GRID_ERR_NONE;
                int rich_rc = slice_resolve_grid_rich(
                    sys, previous_state, x, y, z, universe_depth, current_state,
                    &cell_id, &material_id, &error, active_stats);
                if (rich_rc < 0) {
                    cell_id = -1;
                    material_id = 0;
                    error = GRID_ERR_UNDEFINED;
                    alea_hier_coherence_state_clear(current_state);
                }
                out_cell_ids[idx] = cell_id;
                if (out_material_ids) out_material_ids[idx] = material_id;
                if (out_errors) out_errors[idx] = error;

                alea_hier_coherence_state_t* state_tmp = previous_state;
                previous_state = current_state;
                current_state = state_tmp;

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
            }
        }
        /* Collect per-thread CSG counters after the for barrier */
        if (stats_en) {
            alea_perf_counters_t c = alea_perf_get();
            thread_stats.csg_prim_evals = c.primitive_evaluations;
            thread_stats.csg_bool_ops = c.boolean_operations;
            #pragma omp critical(alea_grid_stats_merge)
            grid_query_stats_merge(&grid_stats, &thread_stats);
        }
    }
#else
    /* Sequential version with row-local rich-path coherence. */
    if (stats_en) alea_perf_reset();

    for (int j = 0; j < nv; j++) {
        if (alea_interrupted()) break;
        double v = v_min + (j + 0.5) * dv;
        alea_hier_coherence_state_t state_a;
        alea_hier_coherence_state_t state_b;
        alea_hier_coherence_state_t* previous_state = &state_a;
        alea_hier_coherence_state_t* current_state = &state_b;
        alea_hier_coherence_state_clear(previous_state);
        alea_hier_coherence_state_clear(current_state);

        for (int i = 0; i < nu; i++) {
            if (alea_interrupted()) break;
            double u = u_min + (i + 0.5) * du;
            int idx = j * nu + i;

            double x = origin[0] + u * u_axis[0] + v * v_axis[0];
            double y = origin[1] + u * u_axis[1] + v * v_axis[1];
            double z = origin[2] + u * u_axis[2] + v * v_axis[2];

            int cell_id = -1;
            int material_id = 0;
            uint8_t error = GRID_ERR_NONE;

            int rich_rc = slice_resolve_grid_rich(
                sys, previous_state, x, y, z, universe_depth, current_state,
                &cell_id, &material_id, &error,
                stats_en ? &grid_stats : NULL);
            if (rich_rc < 0) {
                cell_id = -1;
                material_id = 0;
                error = GRID_ERR_UNDEFINED;
                alea_hier_coherence_state_clear(current_state);
            }
            out_cell_ids[idx] = cell_id;
            if (out_material_ids) out_material_ids[idx] = material_id;
            if (out_errors) out_errors[idx] = error;

            alea_hier_coherence_state_t* state_tmp = previous_state;
            previous_state = current_state;
            current_state = state_tmp;
        }

    }

    if (stats_en) {
        alea_perf_counters_t c = alea_perf_get();
        grid_stats.csg_prim_evals = c.primitive_evaluations;
        grid_stats.csg_bool_ops = c.boolean_operations;
    }

#endif

    if (alea_interrupted()) return -1;

    /* Second pass: recheck boundary pixels for overlap detection.
     * The first-owner coherent resolver finds one cell per pixel without checking
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
            if (alea_interrupted()) continue;
            double v = v_min + (j + 0.5) * dv;
            for (int i = 0; i < nu; i++) {
                if (alea_interrupted()) break;
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
        if (alea_interrupted()) return -1;
        filter_grid_overlap_ambiguities(sys, view, nu, nv, universe_depth,
                                        out_cell_ids, NULL, NULL, out_errors,
                                        NULL);
    }

    /* Print query stats when ALEA_GRID_STATS=1 */
    if (stats_en) {
        const size_t total = grid_stats.px_total;
        fprintf(stdout, "\n[GRID STATS] pixels=%zu\n", total);
        fprintf(stdout, "  root query:         %6zu (%5.1f%%)\n",
                grid_stats.rich_root_queries,
                total ? 100.0 * grid_stats.rich_root_queries / total : 0.0);
        fprintf(stdout, "  complete reuse:     %6zu (%5.1f%%)\n",
                grid_stats.rich_path_reused,
                total ? 100.0 * grid_stats.rich_path_reused / total : 0.0);
        fprintf(stdout, "  lattice transition: %6zu (%5.1f%%)\n",
                grid_stats.rich_lattice_transitions,
                total ? 100.0 * grid_stats.rich_lattice_transitions / total : 0.0);
        fprintf(stdout, "  prefix restart:     %6zu (%5.1f%%)\n",
                grid_stats.rich_prefix_restarts,
                total ? 100.0 * grid_stats.rich_prefix_restarts / total : 0.0);
        fprintf(stdout, "  full fallback:      %6zu (%5.1f%%)\n",
                grid_stats.rich_full_fallbacks,
                total ? 100.0 * grid_stats.rich_full_fallbacks / total : 0.0);
        const size_t prim = grid_stats.csg_prim_evals;
        const size_t bops = grid_stats.csg_bool_ops;
        fprintf(stdout, "[CSG EVALS] prim=%zu bool_ops=%zu per_pixel=%.0f prim+bool=%.0f\n",
                prim, bops, total ? (double)(prim+bops)/total : 0.0,
                total ? (double)(prim+bops)/total : 0.0);
        fflush(stdout);
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

    for (int j = 0; j < nv; j++) {
        double v = v_min + (j + 0.5) * dv;
        alea_hier_coherence_state_t state_a;
        alea_hier_coherence_state_t state_b;
        alea_hier_coherence_state_t* previous_state = &state_a;
        alea_hier_coherence_state_t* current_state = &state_b;
        alea_hier_coherence_state_clear(previous_state);
        alea_hier_coherence_state_clear(current_state);

        for (int i = 0; i < nu; i++) {
            double u = u_min + (i + 0.5) * du;
            size_t idx = (size_t)j * (size_t)nu + (size_t)i;
            double x = origin[0] + u * u_axis[0] + v * v_axis[0];
            double y = origin[1] + u * u_axis[1] + v * v_axis[1];
            double z = origin[2] + u * u_axis[2] + v * v_axis[2];

            int cell_id = -1;
            int material_id = 0;
            uint8_t error = GRID_ERR_NONE;

            int rich_rc = slice_resolve_grid_rich(
                sys, previous_state, x, y, z, universe_depth, current_state,
                &cell_id, &material_id, &error, NULL);
            if (rich_rc < 0) {
                return -1;
            }

            out_cell_ids[idx] = cell_id;
            if (out_material_ids) out_material_ids[idx] = material_id;
            if (out_errors) out_errors[idx] = error;

            if (rich_rc > 0) {
                const alea_hier_ray_path_entry_t* target = NULL;
                uint32_t path_id = UINT32_MAX;
                if (slice_project_hier_path(&current_state->path, universe_depth,
                                            &target) != 0 ||
                    slice_path_table_intern(out_paths, target->universe_id,
                                            target->depth, &target->transform,
                                            &path_id) != 0) {
                    return -1;
                }
                out_path_ids[idx] = path_id;
            }

            alea_hier_coherence_state_t* state_tmp = previous_state;
            previous_state = current_state;
            current_state = state_tmp;
        }

    }

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
            if (alea_interrupted()) continue;
            double v = view->v_min + (j + 0.5) * dv;
            for (int i = 0; i < nu; i++) {
                if (alea_interrupted()) break;
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
        if (alea_interrupted()) return -1;
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
            if (alea_interrupted()) continue;
            double v = view->v_min + (j + 0.5) * dv;
            for (int i = 0; i < nu; i++) {
                if (alea_interrupted()) break;
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
        if (alea_interrupted()) return -1;
    }

    return 0;
}

int alea_find_cells_grid_paths_selected(
    alea_system_t* sys, const alea_slice_view_t* view,
    int nu, int nv, int universe_depth, int tile_w, int tile_h,
    const uint8_t* tile_mask, size_t tile_mask_count,
    int* inout_cell_ids, uint32_t* out_path_ids,
    alea_slice_path_table_t* out_paths) {
    if (!sys || !view || !inout_cell_ids || !out_path_ids || !out_paths ||
        nu <= 0 || nv <= 0 || tile_w <= 0 || tile_h <= 0)
        return -1;

    const size_t tile_cols =
        ((size_t)nu + (size_t)tile_w - 1) / (size_t)tile_w;
    const size_t tile_rows =
        ((size_t)nv + (size_t)tile_h - 1) / (size_t)tile_h;
    if (tile_mask &&
        (tile_cols > SIZE_MAX / tile_rows ||
         tile_mask_count != tile_cols * tile_rows))
        return -1;
    if (alea_system_prepare_query_caches(sys, ALEA_CACHE_RAYCAST) != 0)
        return -1;

    const size_t pixel_count = (size_t)nu * (size_t)nv;
    for (size_t i = 0; i < pixel_count; i++)
        out_path_ids[i] = UINT32_MAX;
    alea_slice_path_table_free(out_paths);

    const alea_slice_plane_t* plane = &view->plane;
    const double du = (view->u_max - view->u_min) / (double)nu;
    const double dv = (view->v_max - view->v_min) / (double)nv;
    for (int j = 0; j < nv; j++) {
        alea_hier_coherence_state_t state_a;
        alea_hier_coherence_state_t state_b;
        alea_hier_coherence_state_t* previous = &state_a;
        alea_hier_coherence_state_t* current = &state_b;
        alea_hier_coherence_state_clear(previous);
        alea_hier_coherence_state_clear(current);
        bool have_previous = false;
        const double v = view->v_min + (j + 0.5) * dv;

        for (int i = 0; i < nu; i++) {
            const size_t tile_index = (size_t)(j / tile_h) * tile_cols +
                                      (size_t)(i / tile_w);
            if (tile_mask && !tile_mask[tile_index]) {
                have_previous = false;
                continue;
            }
            if (!have_previous)
                alea_hier_coherence_state_clear(previous);

            const double u = view->u_min + (i + 0.5) * du;
            const double x = plane->origin[0] + u * plane->u_axis[0] +
                             v * plane->v_axis[0];
            const double y = plane->origin[1] + u * plane->u_axis[1] +
                             v * plane->v_axis[1];
            const double z = plane->origin[2] + u * plane->u_axis[2] +
                             v * plane->v_axis[2];
            int cell_id = -1;
            int material_id = 0;
            uint8_t error = GRID_ERR_NONE;
            int rc = slice_resolve_grid_rich(
                sys, have_previous ? previous : NULL, x, y, z,
                universe_depth, current, &cell_id, &material_id, &error, NULL);
            if (rc < 0) {
                alea_slice_path_table_free(out_paths);
                return -1;
            }
            const size_t pixel_index = (size_t)j * (size_t)nu + (size_t)i;
            inout_cell_ids[pixel_index] = cell_id;
            if (rc > 0) {
                const alea_hier_ray_path_entry_t* target = NULL;
                uint32_t path_id = UINT32_MAX;
                if (slice_project_hier_path(&current->path, universe_depth,
                                            &target) != 0 ||
                    slice_path_table_intern(out_paths, target->universe_id,
                                            target->depth, &target->transform,
                                            &path_id) != 0) {
                    alea_slice_path_table_free(out_paths);
                    return -1;
                }
                out_path_ids[pixel_index] = path_id;
            }
            alea_hier_coherence_state_t* tmp = previous;
            previous = current;
            current = tmp;
            have_previous = true;
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
            #pragma omp atomic
            g_point_coverage_stats.spatial_multi_early_exit++;
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

static int find_point_coverage_exact(alea_system_t* sys,
                                     double gx, double gy, double gz,
                                     int universe_depth,
                                     point_coverage_t* out) {
    #pragma omp atomic
    g_point_coverage_stats.queries++;
    int rc = find_point_coverage_spatial(sys, gx, gy, gz, universe_depth, out);
    if (rc == 0) return 0;
    /* A hierarchy candidate query can be unavailable for a freshly-built
     * programmatic model even though the recursive ownership query is valid.
     * Exact local diagnostics must retain correctness in that case; record
     * the spatial failure but use the same conservative fallback as a capped
     * candidate list. */
    if (rc != -2) {
        #pragma omp atomic
        g_point_coverage_stats.query_errors++;
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

        const alea_bbox_t bbox_v = alea_node_bbox_get(&sys->nodes.data[cell->root_node_id].bbox);
        const alea_bbox_t* bbox = &bbox_v;
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
    dst->skipped_tiles += src->skipped_tiles;
    dst->query_errors += src->query_errors;
    dst->pixels += src->pixels;
    dst->exact_fallback_pixels += src->exact_fallback_pixels;
    dst->skipped_pixels += src->skipped_pixels;
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
    dst->incomplete = dst->incomplete || src->incomplete;
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

        const alea_bbox_t bbox_v = alea_node_bbox_get(&sys->nodes.data[cell->root_node_id].bbox);
        const alea_bbox_t* bbox = &bbox_v;
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
    bool tile_used_fallback = false;
    stats->tiles++;

    size_t missing_path_pixels = 0;
    for (int pj = tj; pj < j_end; pj++) {
        for (int pi = ti; pi < i_end; pi++) {
            size_t pidx = (size_t)pj * (size_t)nu + (size_t)pi;
            if (path_ids[pidx] == UINT32_MAX)
                group_indices[missing_path_pixels++] = pidx;
        }
    }
    if (missing_path_pixels > 0) {
        tile_used_fallback = true;
        stats->exact_fallback_pixels += missing_path_pixels;
        for (size_t p = 0; p < missing_path_pixels; p++) {
            int pi = (int)(group_indices[p] % (size_t)nu);
            int pj = (int)(group_indices[p] / (size_t)nu);
            updated += update_pixel_coverage_exact(
                sys, &view->plane, view->u_min, view->v_min,
                du, dv, nu, nv, universe_depth,
                out_secondary_cell_ids, coverage, errors, pi, pj);
        }
    }

    for (int sj = tj; sj < j_end; sj++) {
        for (int si = ti; si < i_end; si++) {
            size_t seed_idx = (size_t)sj * (size_t)nu + (size_t)si;
            uint32_t path_id = path_ids[seed_idx];
            if (path_id == UINT32_MAX) continue;
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
            if (universe_depth >= 0 && path->depth != universe_depth) {
                tile_used_fallback = true;
                stats->exact_fallback_pixels += group_pixels;
                for (size_t p = 0; p < group_pixels; p++) {
                    int pi = (int)(group_indices[p] % (size_t)nu);
                    int pj = (int)(group_indices[p] / (size_t)nu);
                    updated += update_pixel_coverage_exact(
                        sys, &view->plane, view->u_min, view->v_min,
                        du, dv, nu, nv, universe_depth,
                        out_secondary_cell_ids, coverage, errors, pi, pj);
                }
                continue;
            }
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
                tile_used_fallback = true;
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

    if (tile_used_fallback)
        stats->fallback_tiles++;

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
        if (hit->ancestor_is_lattice[a]) {
            alea_lattice_location_t location;
            if (alea_lattice_locate_point(sys, cell, lx, ly, lz,
                                          &location) != 1 ||
                location.fill_universe !=
                    hit->ancestor_lattice_fill_universes[a] ||
                location.i != hit->ancestor_lattice_i[a] ||
                location.j != hit->ancestor_lattice_j[a] ||
                location.k != hit->ancestor_lattice_k[a] ||
                location.ox != hit->ancestor_lattice_ox[a] ||
                location.oy != hit->ancestor_lattice_oy[a] ||
                location.oz != hit->ancestor_lattice_oz[a]) {
                return 0;
            }
            continue;
        }
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

void alea_tile_refinement_options_init(
    alea_tile_refinement_options_t* options) {
    if (!options) return;
    memset(options, 0, sizeof(*options));
    options->max_candidates = 4096;
    options->max_evaluated_candidates = 4096;
    options->cap_storm_threshold = 16;
    options->path_2d_verify_limit = 20;
    options->path_2d_tile_pad = 1;
    options->parallel = true;
}

/* Exact fallback is an oracle, not an entitlement to silently turn a bounded
 * diagnostic pass into a full raster scan.  Leave over-budget tiles
 * provisional and make that fact observable through tile coverage stats. */
static int refine_tile_exact_with_budget(
    alea_system_t* sys, const alea_slice_view_t* view,
    int nu, int nv, int universe_depth,
    double du, double dv, int ti, int tj, int i_end, int j_end,
    size_t max_exact_fallback_pixels, int* out_secondary_cell_ids,
    uint8_t* coverage, uint8_t* errors, uint8_t* provisional_mask) {
    size_t tile_pixels = (size_t)(i_end - ti) * (size_t)(j_end - tj);
    g_tile_coverage_stats.fallback_tiles++;
    if (max_exact_fallback_pixels > 0 &&
        (g_tile_coverage_stats.exact_fallback_pixels >
             max_exact_fallback_pixels ||
         tile_pixels > max_exact_fallback_pixels -
             g_tile_coverage_stats.exact_fallback_pixels)) {
        g_tile_coverage_stats.skipped_tiles++;
        g_tile_coverage_stats.skipped_pixels += tile_pixels;
        g_tile_coverage_stats.incomplete = true;
        if (provisional_mask) {
            for (int j = tj; j < j_end; j++)
                for (int i = ti; i < i_end; i++)
                    provisional_mask[(size_t)j * (size_t)nu + (size_t)i] = 1;
        }
        return 0;
    }
    g_tile_coverage_stats.exact_fallback_pixels += tile_pixels;
    int updated = 0;
    for (int j = tj; j < j_end; j++) {
        for (int i = ti; i < i_end; i++) {
            updated += update_pixel_coverage_exact(
                sys, &view->plane, view->u_min, view->v_min, du, dv,
                nu, nv, universe_depth, out_secondary_cell_ids, coverage,
                errors, i, j);
        }
    }
    return updated;
}

static int refine_grid_indices_exact_with_budget(
    alea_system_t* sys, const alea_slice_view_t* view,
    int nu, int nv, int universe_depth, double du, double dv,
    const size_t* indices, size_t count, size_t max_exact_fallback_pixels,
    int* out_secondary_cell_ids, uint8_t* coverage, uint8_t* errors,
    uint8_t* provisional_mask, bool* out_skipped) {
    if (out_skipped) *out_skipped = false;
    size_t exact_count = count;
    if (max_exact_fallback_pixels > 0) {
        size_t available =
            g_tile_coverage_stats.exact_fallback_pixels <
                    max_exact_fallback_pixels
                ? max_exact_fallback_pixels -
                      g_tile_coverage_stats.exact_fallback_pixels
                : 0;
        if (exact_count > available)
            exact_count = available;
    }
    if (exact_count < count) {
        size_t skipped_count = count - exact_count;
        g_tile_coverage_stats.skipped_pixels += skipped_count;
        g_tile_coverage_stats.incomplete = true;
        if (out_skipped) *out_skipped = true;
        if (provisional_mask)
            for (size_t i = exact_count; i < count; i++)
                provisional_mask[indices[i]] = 1;
    }
    g_tile_coverage_stats.exact_fallback_pixels += exact_count;
    int updated = 0;
    for (size_t i = 0; i < exact_count; i++) {
        int px = (int)(indices[i] % (size_t)nu);
        int py = (int)(indices[i] / (size_t)nu);
        updated += update_pixel_coverage_exact(
            sys, &view->plane, view->u_min, view->v_min, du, dv, nu, nv,
            universe_depth, out_secondary_cell_ids, coverage, errors, px, py);
    }
    return updated;
}

int alea_refine_grid_coverage_tiles_exact_ex(
    alea_system_t* sys,
    const alea_slice_view_t* view,
    int nu, int nv,
    int universe_depth,
    int tile_w, int tile_h,
    const alea_tile_refinement_options_t* options,
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

    alea_tile_refinement_options_t default_options;
    if (!options) {
        alea_tile_refinement_options_init(&default_options);
        options = &default_options;
    }
    size_t max_hits = options->max_candidates ? options->max_candidates : 4096;
    if (max_hits < 4096) max_hits = 4096;
    alea_spatial_hit_t* hits = malloc(max_hits * sizeof(*hits));
    if (!hits) return -1;
    alea_hier_spatial_chain_hit_t* chain_hits = NULL;

    double u_range = view->u_max - view->u_min;
    double v_range = view->v_max - view->v_min;
    double du = u_range / nu;
    double dv = v_range / nv;
    double eps = fmax(fabs(du), fabs(dv)) * 2.0 + 1e-10;
    bool use_hier_direct_region = options->use_hier_direct_region;
    bool use_hier_chain_candidates = options->use_hier_chain_candidates;
    bool use_chain_bitset = use_hier_chain_candidates && options->use_chain_bitset;
    if (use_hier_chain_candidates) {
        chain_hits = malloc(max_hits * sizeof(*chain_hits));
        if (!chain_hits) {
            free(hits);
            return -1;
        }
    }
    size_t cap_storm_threshold = options->cap_storm_threshold;
    size_t eval_candidate_limit = options->max_evaluated_candidates;
    size_t max_exact_fallback_pixels = options->max_exact_fallback_pixels;
    size_t capped_tile_queries = 0;
    size_t evaluable_tile_queries = 0;
    bool exact_remaining_tiles = false;
    size_t tile_cols = ((size_t)nu + (size_t)tile_w - 1) / (size_t)tile_w;
    size_t tile_rows = ((size_t)nv + (size_t)tile_h - 1) / (size_t)tile_h;
    if (options->tile_mask &&
        (tile_cols > SIZE_MAX / tile_rows ||
         options->tile_mask_count != tile_cols * tile_rows)) {
        free(chain_hits);
        free(hits);
        return -1;
    }
    size_t pixel_count = (size_t)nu * (size_t)nv;
    if (options->provisional_mask &&
        options->provisional_mask_count != pixel_count) {
        free(chain_hits);
        free(hits);
        return -1;
    }
    if (options->provisional_mask)
        memset(options->provisional_mask, 0, pixel_count);
    uint8_t* provisional_mask = options->provisional_mask;
    int updated = 0;

    for (int tj = 0; tj < nv && !alea_interrupted(); tj += tile_h) {
        int j_end = tj + tile_h;
        if (j_end > nv) j_end = nv;
        double v0 = view->v_min + tj * dv;
        double v1 = view->v_min + j_end * dv;

        for (int ti = 0; ti < nu && !alea_interrupted(); ti += tile_w) {
            size_t tile_index = (size_t)(tj / tile_h) * tile_cols +
                                (size_t)(ti / tile_w);
            if (options->tile_mask && !options->tile_mask[tile_index])
                continue;
            int i_end = ti + tile_w;
            if (i_end > nu) i_end = nu;
            double u0 = view->u_min + ti * du;
            double u1 = view->u_min + i_end * du;
            size_t tile_pixels = (size_t)(i_end - ti) * (size_t)(j_end - tj);
            g_tile_coverage_stats.tiles++;

            if (exact_remaining_tiles) {
                updated += refine_tile_exact_with_budget(
                    sys, view, nu, nv, universe_depth, du, dv, ti, tj,
                    i_end, j_end, max_exact_fallback_pixels,
                    out_secondary_cell_ids, coverage, errors, provisional_mask);
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
            } else {
                hit_count = alea_hier_spatial_query_region(
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
                updated += refine_tile_exact_with_budget(
                    sys, view, nu, nv, universe_depth, du, dv, ti, tj,
                    i_end, j_end, max_exact_fallback_pixels,
                    out_secondary_cell_ids, coverage, errors, provisional_mask);
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
                updated += refine_tile_exact_with_budget(
                    sys, view, nu, nv, universe_depth, du, dv, ti, tj,
                    i_end, j_end, max_exact_fallback_pixels,
                    out_secondary_cell_ids, coverage, errors, provisional_mask);
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
                updated += refine_tile_exact_with_budget(
                    sys, view, nu, nv, universe_depth, du, dv, ti, tj,
                    i_end, j_end, max_exact_fallback_pixels,
                    out_secondary_cell_ids, coverage, errors, provisional_mask);
                continue;
            }

            for (int j = tj; j < j_end && !alea_interrupted(); j++) {
                for (int i = ti; i < i_end && !alea_interrupted(); i++) {
                    if (use_hier_chain_candidates) {
                        int rc = update_pixel_coverage_from_chain_candidates(
                            sys, &view->plane, chain_hits, hit_count,
                            view->u_min, view->v_min, du, dv, nu,
                            universe_depth, out_secondary_cell_ids,
                            coverage, errors, i, j);
                        if (rc < 0) {
                            g_tile_coverage_stats.fallback_tiles++;
                            if (max_exact_fallback_pixels > 0 &&
                                g_tile_coverage_stats.exact_fallback_pixels >=
                                    max_exact_fallback_pixels) {
                                g_tile_coverage_stats.skipped_pixels++;
                                g_tile_coverage_stats.incomplete = true;
                                if (provisional_mask)
                                    provisional_mask[(size_t)j * (size_t)nu +
                                                     (size_t)i] = 1;
                                rc = 0;
                            } else {
                                g_tile_coverage_stats.exact_fallback_pixels++;
                                rc = update_pixel_coverage_exact(
                                    sys, &view->plane, view->u_min,
                                    view->v_min, du, dv, nu, nv,
                                    universe_depth, out_secondary_cell_ids,
                                    coverage, errors, i, j);
                            }
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
    if (alea_interrupted())
        return -1;
    g_tile_coverage_stats.refined_pixels = (size_t)updated;
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
    alea_tile_refinement_options_t options;
    alea_tile_refinement_options_init(&options);
    return alea_refine_grid_coverage_tiles_exact_ex(
        sys, view, nu, nv, universe_depth, tile_w, tile_h, &options,
        out_secondary_cell_ids, coverage, errors);
}

int alea_refine_grid_coverage_paths_exact_ex(
    alea_system_t* sys,
    const alea_slice_view_t* view,
    int nu, int nv,
    int universe_depth,
    int tile_w, int tile_h,
    const int* primary_cell_ids,
    const uint32_t* path_ids,
    const alea_slice_path_table_t* paths,
    const alea_tile_refinement_options_t* options,
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

    if (alea_prepare_query_acceleration(sys) != 0)
        return -1;

    alea_tile_refinement_options_t default_options;
    if (!options) {
        alea_tile_refinement_options_init(&default_options);
        options = &default_options;
    }
    size_t max_hits = options->max_candidates ? options->max_candidates : 4096;
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
    bool use_path_2d_index = options->use_path_2d_index;
    bool verify_path_2d_index =
        use_path_2d_index && options->verify_path_2d_index;
    size_t verify_2d_log_limit = options->path_2d_verify_limit;
    size_t verify_2d_log_count = 0;
    path_tile_bucket_t verify_bucket = {0};
    path_slice_index_t* path_indexes = NULL;
    bool* path_index_ready = NULL;
    bool* path_index_disabled = NULL;
    size_t path_bucket_limit = options->path_2d_bucket_limit
        ? options->path_2d_bucket_limit : max_hits;
    int path_tile_pad = options->path_2d_tile_pad;
    if (path_tile_pad < 0) path_tile_pad = 0;
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
    size_t max_exact_fallback_pixels = options->max_exact_fallback_pixels;
    size_t pixel_count = (size_t)nu * (size_t)nv;
    size_t tile_cols = ((size_t)nu + (size_t)tile_w - 1) / (size_t)tile_w;
    size_t tile_rows = ((size_t)nv + (size_t)tile_h - 1) / (size_t)tile_h;
    if (options->tile_mask &&
        (tile_cols > SIZE_MAX / tile_rows ||
         options->tile_mask_count != tile_cols * tile_rows)) {
        free(path_index_disabled); free(path_index_ready); free(path_indexes);
        free(group_second); free(group_count); free(group_z); free(group_y);
        free(group_x); free(group_indices); free(hits);
        return -1;
    }
    if (options->provisional_mask &&
        options->provisional_mask_count != pixel_count) {
        free(path_index_disabled); free(path_index_ready); free(path_indexes);
        free(group_second); free(group_count); free(group_z); free(group_y);
        free(group_x); free(group_indices); free(hits);
        return -1;
    }
    if (options->provisional_mask)
        memset(options->provisional_mask, 0, pixel_count);
    uint8_t* provisional_mask = options->provisional_mask;
    int updated = 0;

#ifdef _OPENMP
    if (!use_path_2d_index && options->parallel &&
        max_exact_fallback_pixels == 0) {
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
                    if (atomic_load(&abort_flag) || alea_interrupted()) {
                        atomic_store(&abort_flag, 1);
                        continue;
                    }
                    size_t tile_index = (size_t)(tj / tile_h) * tile_cols +
                                        (size_t)(ti / tile_w);
                    if (options->tile_mask &&
                        !options->tile_mask[tile_index])
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

        if (atomic_load(&abort_flag) || alea_interrupted())
            return -1;

        g_tile_coverage_stats = parallel_stats;
        g_tile_coverage_stats.refined_pixels = (size_t)parallel_updated;
        return parallel_updated;
    }
#endif

    for (int tj = 0; tj < nv && !alea_interrupted(); tj += tile_h) {
        int j_end = tj + tile_h;
        if (j_end > nv) j_end = nv;

        for (int ti = 0; ti < nu && !alea_interrupted(); ti += tile_w) {
            size_t tile_index = (size_t)(tj / tile_h) * tile_cols +
                                (size_t)(ti / tile_w);
            if (options->tile_mask && !options->tile_mask[tile_index])
                continue;
            int i_end = ti + tile_w;
            if (i_end > nu) i_end = nu;
            g_tile_coverage_stats.tiles++;
            bool tile_used_fallback = false;
            bool tile_skipped_fallback = false;

            size_t missing_path_pixels = 0;
            for (int pj = tj; pj < j_end; pj++) {
                for (int pi = ti; pi < i_end; pi++) {
                    size_t pidx = (size_t)pj * (size_t)nu + (size_t)pi;
                    if (path_ids[pidx] == UINT32_MAX)
                        group_indices[missing_path_pixels++] = pidx;
                }
            }
            if (missing_path_pixels > 0) {
                bool skipped = false;
                tile_used_fallback = true;
                updated += refine_grid_indices_exact_with_budget(
                    sys, view, nu, nv, universe_depth, du, dv,
                    group_indices, missing_path_pixels,
                    max_exact_fallback_pixels, out_secondary_cell_ids,
                    coverage, errors, provisional_mask, &skipped);
                tile_skipped_fallback = tile_skipped_fallback || skipped;
            }

            for (int sj = tj; sj < j_end; sj++) {
                for (int si = ti; si < i_end; si++) {
                    size_t seed_idx = (size_t)sj * (size_t)nu + (size_t)si;
                    uint32_t path_id = path_ids[seed_idx];
                    if (path_id == UINT32_MAX) continue;
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
                    if (universe_depth >= 0 &&
                        path->depth != universe_depth) {
                        bool skipped = false;
                        tile_used_fallback = true;
                        updated += refine_grid_indices_exact_with_budget(
                            sys, view, nu, nv, universe_depth, du, dv,
                            group_indices, group_pixels,
                            max_exact_fallback_pixels,
                            out_secondary_cell_ids, coverage, errors,
                            provisional_mask, &skipped);
                        tile_skipped_fallback =
                            tile_skipped_fallback || skipped;
                        continue;
                    }
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
                        bool skipped = false;
                        tile_used_fallback = true;
                        updated += refine_grid_indices_exact_with_budget(
                            sys, view, nu, nv, universe_depth, du, dv,
                            group_indices, group_pixels,
                            max_exact_fallback_pixels, out_secondary_cell_ids,
                            coverage, errors, provisional_mask, &skipped);
                        tile_skipped_fallback =
                            tile_skipped_fallback || skipped;
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
            if (tile_used_fallback)
                g_tile_coverage_stats.fallback_tiles++;
            if (tile_skipped_fallback)
                g_tile_coverage_stats.skipped_tiles++;
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
    if (alea_interrupted())
        return -1;
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
    alea_tile_refinement_options_t options;
    alea_tile_refinement_options_init(&options);
    return alea_refine_grid_coverage_paths_exact_ex(
        sys, view, nu, nv, universe_depth, tile_w, tile_h,
        primary_cell_ids, path_ids, paths, &options, out_secondary_cell_ids,
        coverage, errors);
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

/** Compute a pixel-space length estimate for parametrically evaluable curves. */
static double estimate_parametric_curve_length(const alea_curve_2d_t* curve,
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
        double arc_angle = curve->bounds.theta_end - curve->bounds.theta_start;
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
        double arc_angle = curve->bounds.theta_end - curve->bounds.theta_start;
        if (arc_angle <= 0) arc_angle += 2 * M_PI;
        return r_pixels * arc_angle;
    }

    /* For lines: use endpoint distance */
    if (curve->type == ALEA_CURVE_LINE ||
        curve->type == ALEA_CURVE_LINE_SEGMENT) {
        double x1, y1, x2, y2;
        if (!alea_curve_eval(curve, curve->bounds.t_min, &x1, &y1)) return 0;
        if (!alea_curve_eval(curve, curve->bounds.t_max, &x2, &y2)) return 0;

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

static void surface_label_param_range(const alea_curve_2d_t* curve,
                                      double x_min, double x_max,
                                      double y_min, double y_max,
                                      double* t_start, double* t_end,
                                      double* t_preferred) {
    if (curve->type == ALEA_CURVE_CIRCLE || curve->type == ALEA_CURVE_ELLIPSE) {
        *t_start = 0.0;
        *t_end = 2.0 * M_PI;
        /* Full closed curves have no natural start.  Use one fixed geometric
         * anchor so labels are predictable; the candidate search below moves
         * only colliding labels to the nearest available point. */
        *t_preferred = M_PI / 4.0;
    } else if (curve->type == ALEA_CURVE_LINE ||
               curve->type == ALEA_CURVE_RAY ||
               curve->type == ALEA_CURVE_PARALLEL_LINES) {
        /* Infinite lines intentionally have no stored curve bounds. Clip the
         * parameter to the viewport before selecting or verifying labels;
         * otherwise every candidate collapses to t=0. */
        const alea_line_2d_t* line = &curve->data.line;
        double lo = -DBL_MAX, hi = DBL_MAX;
        if (fabs(line->direction[0]) > 1e-15) {
            double a = (x_min - line->point[0]) / line->direction[0];
            double b = (x_max - line->point[0]) / line->direction[0];
            if (a > b) { double swap = a; a = b; b = swap; }
            if (a > lo) lo = a;
            if (b < hi) hi = b;
        }
        if (fabs(line->direction[1]) > 1e-15) {
            double a = (y_min - line->point[1]) / line->direction[1];
            double b = (y_max - line->point[1]) / line->direction[1];
            if (a > b) { double swap = a; a = b; b = swap; }
            if (a > lo) lo = a;
            if (b < hi) hi = b;
        }
        if (curve->type == ALEA_CURVE_RAY && lo < curve->bounds.t_min)
            lo = curve->bounds.t_min;
        if (!(lo <= hi) || !isfinite(lo) || !isfinite(hi)) lo = hi = 0.0;
        *t_start = lo;
        *t_end = hi;
        *t_preferred = 0.5 * (lo + hi);
    } else if (curve->bounds.t_min == curve->bounds.t_max) {
        *t_start = 0.0;
        *t_end = 0.0;
        *t_preferred = 0.0;
    } else {
        *t_start = curve->bounds.t_min;
        *t_end = curve->bounds.t_max;
        *t_preferred = (curve->bounds.t_min + curve->bounds.t_max) / 2.0;
    }
}

typedef struct {
    int surf_id;
    int px, py;
    /* Retained internally so provenance need not reconstruct the analytical
     * point from an integer label anchor. Both candidate paths retain an
     * exact changed-edge descriptor when one can be resolved. */
    double u, v;
    int edge_x, edge_y;
    alea_slice_edge_orientation_t edge_orientation;
    double edge_fraction;
    int provenance_group;
} surface_label_candidate_t;

typedef struct surface_label_provenance_cache surface_label_provenance_cache_t;

/* Find a rendered finite right edge intersected by a parametric curve near
 * `target_t`.  The sampling only locates a bracket; the reported fraction is
 * refined against the analytical curve, rather than taken from a pixel guess. */
static double parametric_candidate_find_right_edge(
    const alea_curve_2d_t* curve, double t_start, double t_end, double target_t,
    const int* ids, int width, int height, double x_min, double y_min,
    double dx, double dy, surface_label_candidate_t* out) {
    if (!ids || t_end <= t_start) return DBL_MAX;
    double best = DBL_MAX, pu, pv;
    if (!alea_curve_eval(curve, t_start, &pu, &pv)) return DBL_MAX;
    for (int s = 1; s <= 128; s++) {
        double t = t_start + (t_end - t_start) * (double)s / 128.0, u, v;
        if (!alea_curve_eval(curve, t, &u, &v)) continue;
        int row0 = (int)ceil((fmin(pv, v) - y_min) / dy - 0.5);
        int row1 = (int)floor((fmax(pv, v) - y_min) / dy - 0.5);
        for (int row = row0; row <= row1; row++) {
            if (row < 0 || row >= height || (pv - (y_min + (row + .5) * dy)) *
                (v - (y_min + (row + .5) * dy)) > 0.0) continue;
            double lo = t_start + (t_end - t_start) * (double)(s - 1) / 128.0;
            double hi = t, line = y_min + (row + .5) * dy;
            for (int it = 0; it < 24; it++) {
                double mid = .5 * (lo + hi), mu, mv;
                if (!alea_curve_eval(curve, mid, &mu, &mv)) break;
                double lu, lv;
                alea_curve_eval(curve, lo, &lu, &lv);
                if ((lv - line) * (mv - line) <= 0.0) hi = mid; else lo = mid;
            }
            double hit_t = .5 * (lo + hi), hu, hv;
            if (!alea_curve_eval(curve, hit_t, &hu, &hv)) continue;
            double pixel_u = (hu - x_min) / dx - .5;
            int edge_x = (int)floor(pixel_u);
            if (edge_x < 0 || edge_x + 1 >= width ||
                ids[row * width + edge_x] == ids[row * width + edge_x + 1]) continue;
            double score = fabs(hit_t - target_t);
            if (score < best) {
                best = score; out->edge_x = edge_x; out->edge_y = row;
                out->edge_orientation = ALEA_SLICE_EDGE_RIGHT;
                out->edge_fraction = fmax(0.0, fmin(1.0,
                                                     pixel_u - edge_x));
            }
        }
        pu = u; pv = v;
    }
    return best;
}

/* As above, but resolve the other raster edge family.  A DOWN edge starts at
 * (column, row) and ends at (column, row - 1), hence its fraction decreases
 * with the pixel-space y coordinate. */
static void parametric_candidate_find_down_edge(
    const alea_curve_2d_t* curve, double t_start, double t_end, double target_t,
    const int* ids, int width, int height, double x_min, double y_min,
    double dx, double dy, double best, surface_label_candidate_t* out) {
    if (!ids || t_end <= t_start) return;
    double pu, pv;
    if (!alea_curve_eval(curve, t_start, &pu, &pv)) return;
    for (int s = 1; s <= 128; s++) {
        double t = t_start + (t_end - t_start) * (double)s / 128.0, u, v;
        if (!alea_curve_eval(curve, t, &u, &v)) continue;
        int col0 = (int)ceil((fmin(pu, u) - x_min) / dx - 0.5);
        int col1 = (int)floor((fmax(pu, u) - x_min) / dx - 0.5);
        for (int col = col0; col <= col1; col++) {
            double line = x_min + (col + .5) * dx;
            if (col < 0 || col >= width || (pu - line) * (u - line) > 0.0)
                continue;
            double lo = t_start + (t_end - t_start) * (double)(s - 1) / 128.0;
            double hi = t;
            for (int it = 0; it < 24; it++) {
                double mid = .5 * (lo + hi), mu, mv, lu, lv;
                if (!alea_curve_eval(curve, mid, &mu, &mv) ||
                    !alea_curve_eval(curve, lo, &lu, &lv)) break;
                if ((lu - line) * (mu - line) <= 0.0) hi = mid; else lo = mid;
            }
            double hit_t = .5 * (lo + hi), hu, hv;
            if (!alea_curve_eval(curve, hit_t, &hu, &hv)) continue;
            double pixel_v = (hv - y_min) / dy - .5;
            int edge_y = (int)ceil(pixel_v);
            if (edge_y <= 0 || edge_y >= height ||
                ids[edge_y * width + col] == ids[(edge_y - 1) * width + col])
                continue;
            double score = fabs(hit_t - target_t);
            if (score < best) {
                best = score; out->edge_x = col; out->edge_y = edge_y;
                out->edge_orientation = ALEA_SLICE_EDGE_DOWN;
                out->edge_fraction = fmax(0.0, fmin(1.0,
                    edge_y - pixel_v));
            }
        }
        pu = u; pv = v;
    }
}

static int parametric_surface_label_candidate(const alea_curve_2d_t* curve,
                                              const int* boundary_ids,
                                              double x_min, double y_min,
                                              double dx, double dy,
                                              int width, int height,
                                              int margin,
                                              const surface_label_candidate_t* existing_labels,
                                              int existing_count,
                                              double min_label_spacing,
                                              double* out_visible_fraction,
                                              surface_label_candidate_t* out_candidate) {
    double t_start, t_end, t_preferred;
    surface_label_param_range(curve, x_min, x_min + width * dx,
                              y_min, y_min + height * dy,
                              &t_start, &t_end, &t_preferred);

    int samples = boundary_ids ? 32 : 0;
    typedef struct { double t, u, v; int ix, iy, visible; }
        parametric_label_sample_t;
    parametric_label_sample_t points[33] = {{0}};
    int visible = 0, run_start = -1, run_length = 0;
    int best_run_start = -1, best_run_length = 0;
    for (int s = 0; s <= samples; s++) {
        double t = samples ? t_start + (t_end - t_start) * (double)s / samples
                           : t_preferred;
        points[s].t = t;
        if (!alea_curve_eval(curve, t, &points[s].u, &points[s].v)) continue;
        points[s].ix = (int)((points[s].u - x_min) / dx);
        points[s].iy = (int)((points[s].v - y_min) / dy);
        /* Visibility deliberately precedes the label margin.  A clipped arc
         * gets one run identifier even when its endpoint is too close to the
         * viewport for text. */
        points[s].visible = point_has_drawn_contour_nearby(
            boundary_ids, width, height, points[s].ix, points[s].iy);
        if (!points[s].visible) {
            if (run_length > best_run_length) {
                best_run_start = run_start;
                best_run_length = run_length;
            }
            run_start = -1;
            run_length = 0;
            continue;
        }
        visible++;
        if (run_start < 0) run_start = s;
        run_length++;
    }
    if (run_length > best_run_length) {
        best_run_start = run_start;
        best_run_length = run_length;
    }

    double best_score = DBL_MAX;
    double best_t = t_preferred;
    int found = 0;
    int full_closed = (curve->type == ALEA_CURVE_CIRCLE ||
                       curve->type == ALEA_CURVE_ELLIPSE) &&
                      visible == samples + 1;
    double run_preferred = full_closed ? t_preferred :
        (best_run_length > 0
            ? points[best_run_start + (best_run_length - 1) / 2].t
            : t_preferred);
    for (int s = best_run_start; s >= 0 && s < best_run_start + best_run_length;
         s++) {
        int ix = points[s].ix, iy = points[s].iy;
        if (ix < margin || ix >= width - margin ||
            iy < margin || iy >= height - margin)
            continue;
        double score = fabs(points[s].t - run_preferred);
        if (full_closed) {
            double period = 2.0 * M_PI;
            score = fmod(score, period);
            if (score > period * 0.5) score = period - score;
        }

        /* Prefer another point on a closed curve over putting its text on an
         * already placed surface label.  The large penalty retains the
         * geometry-derived preferred point whenever it is clear, yet gives
         * concentric or nearly coincident curves a usable alternative. */
        for (int i = 0; i < existing_count; i++) {
            double ddx = ix - existing_labels[i].px;
            double ddy = iy - existing_labels[i].py;
            double distance_sq = ddx * ddx + ddy * ddy;
            double spacing_sq = min_label_spacing * min_label_spacing;
            if (distance_sq < spacing_sq)
                score += 10000.0 + spacing_sq - distance_sq;
        }

        if (!found || score < best_score) {
            found = 1;
            best_score = score;
            out_candidate->px = ix;
            out_candidate->py = iy;
            out_candidate->u = points[s].u;
            out_candidate->v = points[s].v;
            out_candidate->edge_x = -1;
            out_candidate->edge_y = -1;
            out_candidate->provenance_group = -1;
            best_t = points[s].t;
        }
    }

    if (found) {
        double edge_score = parametric_candidate_find_right_edge(
            curve, t_start, t_end, best_t, boundary_ids, width, height,
            x_min, y_min, dx, dy, out_candidate);
        parametric_candidate_find_down_edge(
            curve, t_start, t_end, best_t, boundary_ids, width, height,
            x_min, y_min, dx, dy, edge_score, out_candidate);
    }

    *out_visible_fraction = (double)visible / (double)(samples + 1);
    return found;
}

static void sort_surface_label_roots(double* roots, int count) {
    for (int i = 1; i < count; i++) {
        double value = roots[i];
        int j = i - 1;
        while (j >= 0 && roots[j] > value) {
            roots[j + 1] = roots[j];
            j--;
        }
        roots[j + 1] = value;
    }
}

/* General conics and torus quartics have no stable single parameter in the
 * slice API. Walk horizontal viewport scanlines using the same intersection
 * solver as rasterization, then select a central, non-tangent visible point. */
static int scanline_surface_label_candidate(const alea_curve_2d_t* curve,
                                            const int* boundary_ids,
                                            double x_min, double x_max,
                                            double y_min, double y_max,
                                            double dx, double dy,
                                            int width, int height,
                                            int margin,
                                            const surface_label_candidate_t* existing_labels,
                                            int existing_count,
                                            double min_label_spacing,
                                            double* out_length,
    surface_label_candidate_t* out_candidate) {
    int visible_samples = 0;
    (void)y_max;
    int run_start[8], run_length[8], last_row[8];
    for (int r = 0; r < 8; r++) {
        run_start[r] = -1;
        run_length[r] = 0;
        last_row[r] = -2;
    }
    int best_branch = -1, best_start = -1, best_length = 0;

    /* First pass identifies visible connected runs. Root ordinal is stable
     * between adjacent scanlines except at a genuine branch merge/split,
     * which intentionally terminates a run rather than joining unrelated
     * pieces of a quartic or hyperbola. */
    for (int iy = margin; iy < height - margin; iy++) {
        double v = y_min + (iy + 0.5) * dy;
        double roots[8];
        int count = alea_curve_scanline_intersect(curve, v, roots, 8);
        if (count <= 0) continue;
        sort_surface_label_roots(roots, count);

        for (int r = 0; r < count; r++) {
            double u = roots[r];
            if (!isfinite(u) || u < x_min || u >= x_max) continue;
            int ix = (int)((u - x_min) / dx);
            if (ix < margin || ix >= width - margin) continue;
            if (!point_has_drawn_contour_nearby(boundary_ids, width, height,
                                                ix, iy)) {
                continue;
            }
            visible_samples++;
            if (last_row[r] == iy - 1) {
                run_length[r]++;
            } else {
                run_start[r] = iy;
                run_length[r] = 1;
            }
            last_row[r] = iy;
            if (run_length[r] > best_length) {
                best_branch = r;
                best_start = run_start[r];
                best_length = run_length[r];
            }
        }
    }

    double best_score = DBL_MAX;
    int found = 0;
    double preferred_y = best_start + 0.5 * (best_length - 1);
    for (int iy = best_start; best_branch >= 0 && iy < best_start + best_length;
         iy++) {
        double v = y_min + (iy + 0.5) * dy;
        double roots[8];
        int count = alea_curve_scanline_intersect(curve, v, roots, 8);
        if (count <= best_branch) continue;
        sort_surface_label_roots(roots, count);
        double u = roots[best_branch];
        if (!isfinite(u) || u < x_min || u >= x_max) continue;
        int ix = (int)((u - x_min) / dx);
        if (ix < margin || ix >= width - margin ||
            !point_has_drawn_contour_nearby(boundary_ids, width, height, ix, iy))
            continue;
        double score = fabs(iy - preferred_y);
        for (int i = 0; i < existing_count; i++) {
            double ddx = ix - existing_labels[i].px;
            double ddy = iy - existing_labels[i].py;
            double distance_sq = ddx * ddx + ddy * ddy;
            double spacing_sq = min_label_spacing * min_label_spacing;
            if (distance_sq < spacing_sq)
                score += 10000.0 + spacing_sq - distance_sq;
        }
        double nearest = DBL_MAX;
        if (best_branch > 0)
            nearest = fmin(nearest,
                           fabs(roots[best_branch] - roots[best_branch - 1]) / dx);
        if (best_branch + 1 < count)
            nearest = fmin(nearest,
                           fabs(roots[best_branch + 1] - roots[best_branch]) / dx);
        if (nearest < 2.0) score += (double)width * height * 4.0;
        if (!found || score < best_score) {
            found = 1;
            best_score = score;
            out_candidate->px = ix;
            out_candidate->py = iy;
            out_candidate->u = u;
            out_candidate->v = v;
            /* A scanline root lies exactly on this finite horizontal
             * pixel-centre edge when the visible grid transition is the
             * right edge below. Retain it for the exact-edge verifier. */
            double pixel_u = (u - x_min) / dx - 0.5;
            int edge_x = (int)floor(pixel_u);
            if (boundary_ids && edge_x >= 0 && edge_x + 1 < width &&
                boundary_ids[iy * width + edge_x] !=
                    boundary_ids[iy * width + edge_x + 1]) {
                out_candidate->edge_x = edge_x;
                out_candidate->edge_y = iy;
                out_candidate->edge_orientation = ALEA_SLICE_EDGE_RIGHT;
                out_candidate->edge_fraction = pixel_u - edge_x;
            } else {
                out_candidate->edge_x = -1;
                out_candidate->edge_y = -1;
                out_candidate->provenance_group = -1;
            }
        }
    }

    /* Each retained crossing accounts for at least one pixel of vertical arc
     * length. This conservative lower bound naturally handles multiple conic
     * branches and multiple quartic loops. */
    *out_length = (double)visible_samples;
    return found;
}

static int find_surface_label_candidates_on_boundaries(
    const alea_slice_curves_t* curves,
    const int* boundary_ids,
    double x_min, double x_max,
    double y_min, double y_max,
    int width, int height,
    int margin,
    surface_label_candidate_t** out_labels,
    int* out_count)
{
    if (!curves || !out_labels || !out_count || width <= 0 || height <= 0) {
        return -1;
    }

    *out_labels = NULL;
    *out_count = 0;

    const alea_curve_collection_t* collection = &curves->internal;
    size_t num_curves = collection->curves.count;
    if (num_curves == 0) return 0;

    double dx = (x_max - x_min) / width;
    double dy = (y_max - y_min) / height;

    /* Minimum length in pixels of the part of a curve that is actually on
     * display, before it earns a label. Both candidate searches report the
     * visible length, not the analytical one: a curve the CSG leaves as a
     * sliver must not be labelled as if the whole conic were drawn. */
    const double MIN_CURVE_LENGTH = 30.0;

    /* Minimum distance between labels of same surface (in pixels) */
    const double MIN_LABEL_SPACING = 50.0;

    /* Temporary storage for valid label positions */
    surface_label_candidate_t* temp_labels =
        malloc(num_curves * sizeof(*temp_labels));
    if (!temp_labels) return -1;
    int temp_count = 0;

    /* First pass: collect all valid label positions */
    for (size_t i = 0; i < num_curves; i++) {
        const alea_curve_2d_t* curve = &collection->curves.data[i];
        int surf_id = curve->surface_id;
        if (surf_id <= 0) continue;

        surface_label_candidate_t candidate = {
            .surf_id = surf_id, .edge_x = -1, .edge_y = -1,
            .provenance_group = -1
        };
        double curve_len;
        int found;
        if (curve->type == ALEA_CURVE_PARABOLA ||
            curve->type == ALEA_CURVE_HYPERBOLA ||
            curve->type == ALEA_CURVE_QUARTIC) {
            found = scanline_surface_label_candidate(
                curve, boundary_ids, x_min, x_max, y_min, y_max,
                dx, dy, width, height, margin, temp_labels, temp_count,
                MIN_LABEL_SPACING, &curve_len, &candidate);
        } else {
            /* Uniform parameter samples approximate arc length well for
             * circles and acceptably for the eccentric ellipses a slice
             * produces; the estimate only has to resolve slivers. */
            double visible_fraction = 0.0;
            found = parametric_surface_label_candidate(
                curve, boundary_ids, x_min, y_min, dx, dy,
                width, height, margin, temp_labels, temp_count,
                MIN_LABEL_SPACING, &visible_fraction, &candidate);
            curve_len = estimate_parametric_curve_length(
                curve, x_min, x_max, y_min, y_max, width, height) * visible_fraction;
        }
        if (!found || curve_len < MIN_CURVE_LENGTH) {
            continue;
        }

        /* A surface may generate several analytical curve pieces.  Keep the
         * historical one-label-per-surface spacing rule for those pieces;
         * different surfaces are instead allowed to search for separate
         * points above, so a valid label is never silently lost. */
        int too_close = 0;
        for (int j = 0; j < temp_count; j++) {
            if (temp_labels[j].surf_id != surf_id) continue;
            double ddx = candidate.px - temp_labels[j].px;
            double ddy = candidate.py - temp_labels[j].py;
            if (ddx*ddx + ddy*ddy < MIN_LABEL_SPACING * MIN_LABEL_SPACING) {
                too_close = 1;
                break;
            }
        }
        if (too_close) continue;

        temp_labels[temp_count] = candidate;
        temp_count++;
    }

    if (temp_count == 0) {
        free(temp_labels);
        return 0;
    }

    *out_labels = temp_labels;
    *out_count = temp_count;
    return 0;
}

int alea_find_surface_label_positions_on_boundaries(
    const alea_slice_curves_t* curves, const int* boundary_ids,
    double x_min, double x_max, double y_min, double y_max,
    int width, int height, int margin,
    alea_label_position_t** out_labels, int* out_count) {
    surface_label_candidate_t* candidates = NULL;
    int count = 0;
    int rc = find_surface_label_candidates_on_boundaries(
        curves, boundary_ids, x_min, x_max, y_min, y_max, width, height,
        margin, &candidates, &count);
    if (rc != 0 || count == 0) return rc;
    alea_label_position_t* labels = calloc((size_t)count, sizeof(*labels));
    if (!labels) { free(candidates); return -1; }
    for (int i = 0; i < count; i++) {
        labels[i].id = candidates[i].surf_id;
        labels[i].px = candidates[i].px;
        labels[i].py = candidates[i].py;
    }
    free(candidates);
    *out_labels = labels;
    *out_count = count;
    return 0;
}

static void slice_world_point(const alea_slice_view_t* view, double u, double v,
                              double out[3]);

/* Defined beside the boundary-map canonical trace below.  Keeping the fast
 * label verifier on that contract prevents a raw event at a coincident or
 * cell-only transition from being accepted merely because it has the right
 * surface ID. */
static int label_edge_has_surface_canonical(
    alea_system_t* sys, const double start[3], const double end[3],
    double fraction, double fraction_tolerance, int surface_id,
    int material_boundary, surface_label_candidate_t* candidate,
    surface_label_provenance_cache_t* cache);
static surface_label_provenance_cache_t*
surface_label_provenance_cache_create(void);
static void surface_label_provenance_cache_destroy(
    surface_label_provenance_cache_t* cache);

static int label_point_has_surface_provenance(
    alea_system_t* sys, const alea_slice_view_t* view,
    const int* boundary_ids, int width, int height,
    double x_min, double y_min, double du, double dv, double u, double v,
    int surface_id, int material_boundary,
    surface_label_candidate_t* candidate,
    surface_label_provenance_cache_t* cache) {
    if (!boundary_ids) return 0;
    double px = (u - x_min) / du - 0.5;
    double py = (v - y_min) / dv - 0.5;
    int cx = (int)floor(px), cy = (int)floor(py);
    int best_x = -1, best_y = -1, best_orient = -1;
    double best_fraction = 0.0, best_distance = DBL_MAX;
    if (candidate && candidate->edge_x >= 0) {
        best_x = candidate->edge_x;
        best_y = candidate->edge_y;
        best_orient = candidate->edge_orientation;
        best_fraction = candidate->edge_fraction;
        best_distance = 0.0;
    }
    for (int y = cy - 1; best_orient < 0 && y <= cy + 1; y++)
    for (int x = cx - 1; x <= cx + 1; x++) {
        if (x < 0 || x >= width || y < 0 || y >= height) continue;
        if (x + 1 < width && boundary_ids[y * width + x] != boundary_ids[y * width + x + 1]) {
            double fraction = fmax(0.0, fmin(1.0, px - x));
            double ex = x + fraction, ey = y;
            double distance = hypot(px - ex, py - ey);
            if (distance < best_distance) {
                best_distance = distance; best_fraction = fraction;
                best_x = x; best_y = y; best_orient = ALEA_SLICE_EDGE_RIGHT;
            }
        }
        if (y > 0 && boundary_ids[y * width + x] != boundary_ids[(y - 1) * width + x]) {
            double fraction = fmax(0.0, fmin(1.0, y - py));
            double ex = x, ey = y - fraction;
            double distance = hypot(px - ex, py - ey);
            if (distance < best_distance) {
                best_distance = distance; best_fraction = fraction;
                best_x = x; best_y = y; best_orient = ALEA_SLICE_EDGE_DOWN;
            }
        }
    }
    if (best_orient < 0 || best_distance > 0.75) return 0;
    int nx = best_x + (best_orient == ALEA_SLICE_EDGE_RIGHT);
    int ny = best_y - (best_orient == ALEA_SLICE_EDGE_DOWN);
    double su = x_min + (best_x + 0.5) * du;
    double sv = y_min + (best_y + 0.5) * dv;
    double eu = x_min + (nx + 0.5) * du;
    double ev = y_min + (ny + 0.5) * dv;
    double origin[3], end[3];
    slice_world_point(view, su, sv, origin);
    slice_world_point(view, eu, ev, end);
    double length = 0.0;
    for (int c = 0; c < 3; c++) {
        double delta = end[c] - origin[c];
        length += delta * delta;
    }
    length = sqrt(length);
    if (!(length > 0.0)) return 0;
    return label_edge_has_surface_canonical(
        sys, origin, end, best_fraction,
        candidate && candidate->edge_x >= 0 ? 1e-4 : 0.30,
        surface_id, material_boundary, candidate, cache);
}

int alea_find_surface_label_positions_with_provenance(
    alea_system_t* sys, const alea_slice_view_t* view,
    const alea_slice_curves_t* curves, const int* boundary_ids,
    double x_min, double x_max, double y_min, double y_max,
    int width, int height, int margin, int material_boundary,
    alea_label_position_t** out_labels, int* out_count) {
    surface_label_candidate_t* candidates = NULL;
    int candidate_count = 0;
    int rc = find_surface_label_candidates_on_boundaries(
        curves, boundary_ids, x_min, x_max, y_min, y_max,
        width, height, margin, &candidates, &candidate_count);
    if (rc != 0 || candidate_count == 0) return rc;
    *out_labels = calloc((size_t)candidate_count, sizeof(**out_labels));
    if (!*out_labels) { free(candidates); return -1; }
    *out_count = candidate_count;
    for (int i = 0; i < candidate_count; i++) {
        (*out_labels)[i].id = candidates[i].surf_id;
        (*out_labels)[i].px = candidates[i].px;
        (*out_labels)[i].py = candidates[i].py;
        (*out_labels)[i].provenance_edge_x = -1;
        (*out_labels)[i].provenance_edge_y = -1;
        (*out_labels)[i].provenance_orientation = -1;
        (*out_labels)[i].provenance_group = -1;
    }
    double du = (x_max - x_min) / width, dv = (y_max - y_min) / height;
    surface_label_provenance_cache_t* provenance_cache =
        surface_label_provenance_cache_create();
    int write = 0;
    for (int i = 0; i < *out_count; i++) {
        int verified = label_point_has_surface_provenance(
            sys, view, boundary_ids, width, height, x_min, y_min, du, dv,
            candidates[i].u, candidates[i].v, candidates[i].surf_id,
            material_boundary, &candidates[i], provenance_cache);
        for (size_t c = 0; c < curves->internal.curves.count && !verified; c++) {
            const alea_curve_2d_t* curve = &curves->internal.curves.data[c];
            if (curve->surface_id != (*out_labels)[i].id) continue;
            if (curve->type == ALEA_CURVE_PARABOLA ||
                curve->type == ALEA_CURVE_HYPERBOLA ||
                curve->type == ALEA_CURVE_QUARTIC) {
                /* These curves have no stable public parameter. Keep a
                 * deterministic shortlist of visible scanline roots rather
                 * than requiring the original x coordinate. The original
                 * root can be a cell-only crossing while another branch is
                 * the displayed boundary for this surface. */
                typedef struct { surface_label_candidate_t point; double score; }
                    scanline_label_candidate_t;
                scanline_label_candidate_t candidates[64];
                int candidate_count = 0;
                for (int delta = 0; delta < height && candidate_count < 64;
                     delta++) {
                    int rows[2] = {(*out_labels)[i].py + delta,
                                   (*out_labels)[i].py - delta};
                    int row_count = delta == 0 ? 1 : 2;
                    for (int row = 0; row < row_count && candidate_count < 64;
                         row++) {
                    int iy = rows[row];
                    if (iy < 0 || iy >= height) continue;
                    double v = y_min + (iy + 0.5) * dv;
                    double roots[8];
                    int root_count = alea_curve_scanline_intersect(
                        curve, v, roots, (int)(sizeof(roots) / sizeof(roots[0])));
                    for (int r = 0; r < root_count && candidate_count < 64;
                         r++) {
                        double u = roots[r];
                        int px = (int)((u - x_min) / du);
                        if (!isfinite(u) || px < 0 || px >= width ||
                            !point_has_drawn_contour_nearby(
                                boundary_ids, width, height, px, iy))
                            continue;
                        double pixel_u = (u - x_min) / du - 0.5;
                        int edge_x = (int)floor(pixel_u);
                        surface_label_candidate_t point = {
                            .surf_id = curve->surface_id, .px = px, .py = iy,
                            .u = u, .v = v, .edge_x = -1, .edge_y = -1,
                            .provenance_group = -1};
                        if (edge_x >= 0 && edge_x + 1 < width &&
                            boundary_ids[iy * width + edge_x] !=
                                boundary_ids[iy * width + edge_x + 1]) {
                            point.edge_x = edge_x;
                            point.edge_y = iy;
                            point.edge_orientation = ALEA_SLICE_EDGE_RIGHT;
                            point.edge_fraction = pixel_u - edge_x;
                        }
                        candidates[candidate_count++] = (scanline_label_candidate_t){
                            .point = point,
                            .score = hypot(px - (*out_labels)[i].px,
                                           iy - (*out_labels)[i].py)};
                    }
                }
                }
                for (int a = 1; a < candidate_count; a++) {
                    scanline_label_candidate_t value = candidates[a];
                    int b = a - 1;
                    while (b >= 0 && candidates[b].score > value.score) {
                        candidates[b + 1] = candidates[b]; b--;
                    }
                    candidates[b + 1] = value;
                }
                for (int c = 0; c < candidate_count; c++) {
                    verified = label_point_has_surface_provenance(
                        sys, view, boundary_ids, width, height, x_min, y_min,
                        du, dv, candidates[c].point.u, candidates[c].point.v,
                        curve->surface_id, material_boundary,
                        &candidates[c].point, provenance_cache);
                    if (verified) {
                        (*out_labels)[i].px = candidates[c].point.px;
                        (*out_labels)[i].py = candidates[c].point.py;
                        break;
                    }
                }
                continue;
            }
            double lo, hi, preferred;
            surface_label_param_range(curve, x_min, x_max, y_min, y_max,
                                      &lo, &hi, &preferred);
            typedef struct {
                surface_label_candidate_t point;
                double score;
            } param_label_candidate_t;
            param_label_candidate_t candidates[33];
            int candidate_count = 0;
            for (int s = 0; s <= 32; s++) {
                double t = lo + (hi - lo) * (double)s / 32.0, u, v;
                if (!alea_curve_eval(curve, t, &u, &v)) continue;
                int px = (int)((u - x_min) / du), py = (int)((v - y_min) / dv);
                if (px < 0 || px >= width || py < 0 || py >= height ||
                    !point_has_drawn_contour_nearby(boundary_ids, width, height,
                                                     px, py))
                    continue;
                surface_label_candidate_t point = {
                    .surf_id = curve->surface_id, .px = px, .py = py,
                    .u = u, .v = v, .edge_x = -1, .edge_y = -1,
                    .provenance_group = -1
                };
                double edge_score = parametric_candidate_find_right_edge(
                    curve, lo, hi, t, boundary_ids, width, height,
                    x_min, y_min, du, dv, &point);
                parametric_candidate_find_down_edge(
                    curve, lo, hi, t, boundary_ids, width, height,
                    x_min, y_min, du, dv, edge_score, &point);
                if (point.edge_x < 0) continue;
                candidates[candidate_count++] = (param_label_candidate_t){
                    .point = point,
                    .score = hypot(px - (*out_labels)[i].px,
                                   py - (*out_labels)[i].py)};
            }
            for (int a = 1; a < candidate_count; a++) {
                param_label_candidate_t value = candidates[a];
                int b = a - 1;
                while (b >= 0 && candidates[b].score > value.score) {
                    candidates[b + 1] = candidates[b]; b--;
                }
                candidates[b + 1] = value;
            }
            for (int s = 0; s < candidate_count; s++) {
                verified = label_point_has_surface_provenance(
                    sys, view, boundary_ids, width, height, x_min, y_min, du, dv,
                    candidates[s].point.u, candidates[s].point.v,
                    curve->surface_id, material_boundary,
                    &candidates[s].point, provenance_cache);
                if (verified) {
                    (*out_labels)[i].px = candidates[s].point.px;
                    (*out_labels)[i].py = candidates[s].point.py;
                    break;
                }
            }
        }
        if (verified) {
            if (candidates[i].provenance_group >= 0) {
                (*out_labels)[i].provenance_edge_x = candidates[i].edge_x;
                (*out_labels)[i].provenance_edge_y = candidates[i].edge_y;
                (*out_labels)[i].provenance_orientation =
                    candidates[i].edge_orientation;
                (*out_labels)[i].provenance_group =
                    candidates[i].provenance_group;
            }
            (*out_labels)[write++] = (*out_labels)[i];
        }
    }
    free(candidates);
    surface_label_provenance_cache_destroy(provenance_cache);
    *out_count = write;
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

    n = alea_hier_spatial_find_cells_at_point_uncached(sys, x, y, z,
                                                       hits, MAX_HITS);
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
    alea_hier_spatial_reset_cache();

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
            comp.representative_i = i0;
            comp.representative_j = j0;

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

int alea_find_local_coverage_components(
    alea_system_t* sys, const alea_slice_view_t* view,
    int nu, int nv, int universe_depth,
    size_t max_pixels, size_t max_scratch_bytes, int max_workers,
    alea_plot_error_component_result_t** out_components,
    alea_local_coverage_stats_t* out_stats) {
    if (out_components) *out_components = NULL;
    if (out_stats) memset(out_stats, 0, sizeof(*out_stats));
    if (!sys || !view || !out_components || nu <= 0 || nv <= 0 || max_workers < 0)
        return -1;

    const size_t width = (size_t)nu;
    const size_t height = (size_t)nv;
    if (height > SIZE_MAX / width) return -1;
    const size_t pixels = width * height;
    /* primary + secondary IDs, coverage + error bytes, then classifier's
     * visited bitmap + int queue. Keep this exact so a caller can enforce a
     * predictable hard cap before any allocation. */
    if (pixels > SIZE_MAX / (3u * sizeof(int) + 3u)) return -1;
    const size_t scratch_bytes = pixels * (3u * sizeof(int) + 3u);
    if ((max_pixels && pixels > max_pixels) ||
        (max_scratch_bytes && scratch_bytes > max_scratch_bytes))
        return ALEA_LOCAL_COVERAGE_BUDGET_EXCEEDED;
    if (out_stats) {
        out_stats->pixels = pixels;
        out_stats->scratch_bytes = scratch_bytes;
        out_stats->worker_limit = max_workers > 0 ? max_workers : 1;
    }

    int* cell_ids = malloc(pixels * sizeof(*cell_ids));
    int* secondary_ids = malloc(pixels * sizeof(*secondary_ids));
    uint8_t* coverage = malloc(pixels * sizeof(*coverage));
    uint8_t* errors = malloc(pixels * sizeof(*errors));
    if (!cell_ids || !secondary_ids || !coverage || !errors) {
        free(cell_ids); free(secondary_ids); free(coverage); free(errors);
        return -1;
    }

    /* Do not first construct a coherent fast grid and then replace every
     * pixel with exact coverage: an exact local scan has no use for that
     * intermediate raster. Each point uses the hierarchical spatial candidate
     * path directly, with recursive traversal reserved for explicit fallback
     * cases (lattices or saturated candidate sets). */
    if (alea_system_prepare_query_caches(sys, ALEA_CACHE_RAYCAST) != 0) {
        free(cell_ids); free(secondary_ids); free(coverage); free(errors);
        return -1;
    }
    alea_point_coverage_stats_reset();
    const alea_slice_plane_t* plane = &view->plane;
    const double du = (view->u_max - view->u_min) / (double)nu;
    const double dv = (view->v_max - view->v_min) / (double)nv;
    int query_failed = 0;
#ifdef _OPENMP
    const int worker_limit = max_workers > 0 ? max_workers : 1;
    #pragma omp parallel for schedule(dynamic, 4) num_threads(worker_limit)
#endif
    for (int j = 0; j < nv; j++) {
        if (alea_interrupted()) continue;
        const double v = view->v_min + (j + 0.5) * dv;
        for (int i = 0; i < nu; i++) {
            if (alea_interrupted()) break;
            const size_t idx = (size_t)j * width + (size_t)i;
            const double u = view->u_min + (i + 0.5) * du;
            point_coverage_t pc;
            if (find_point_coverage_exact_uv(
                    sys, plane, u, v, universe_depth, &pc) != 0) {
                #pragma omp atomic write
                query_failed = 1;
                break;
            }
            cell_ids[idx] = pc.primary_cell_id;
            secondary_ids[idx] = pc.secondary_cell_id;
            coverage[idx] = pc.coverage;
            errors[idx] = (pc.coverage == ALEA_COVERAGE_MULTI)
                ? GRID_ERR_OVERLAP
                : (pc.coverage == ALEA_COVERAGE_NONE)
                    ? GRID_ERR_UNDEFINED : GRID_ERR_NONE;
        }
    }
    if (query_failed || alea_interrupted()) {
        free(cell_ids); free(secondary_ids); free(coverage); free(errors);
        return -1;
    }
    /* The filter is conservative: an individual ambiguous probe can decline
     * refinement while leaving the original exact classification intact.
     * Match the full-grid contract and reserve failure for cancellation. */
    (void)alea_filter_grid_boundary_ambiguities(
        sys, view, nu, nv, universe_depth,
        cell_ids, secondary_ids, coverage, errors, NULL);
    if (alea_interrupted()) {
        free(cell_ids); free(secondary_ids); free(coverage); free(errors);
        return -1;
    }

    alea_plot_error_component_result_t* components =
        alea_classify_plot_error_components(
            cell_ids, secondary_ids, coverage, nu, nv);
    free(cell_ids); free(secondary_ids); free(coverage); free(errors);
    if (!components) return -1;
    if (out_stats)
        out_stats->point_coverage = alea_point_coverage_stats_get();
    *out_components = components;
    return 0;
}

/* ============================================================================
 * SURFACE BOUNDARY PROVENANCE
 * ============================================================================ */

typedef struct {
    double fraction;
    int ids[64];
    size_t count;
} boundary_trace_group_t;

typedef struct {
    int ids[64];
    size_t count;
    size_t group_count;
    boundary_trace_group_t groups[64];
    int saw_synthetic;
    int saw_gap;
    int saw_overlap;
    int saw_unresolved;
} boundary_trace_t;

typedef struct {
    int edge_x, edge_y;
    alea_slice_edge_orientation_t orientation;
    double fraction;
    int surface_id;
    int material_boundary;
    int verified;
    int group;
} surface_label_provenance_cache_entry_t;

struct surface_label_provenance_cache {
    surface_label_provenance_cache_entry_t* entries;
    size_t count, capacity;
};

static surface_label_provenance_cache_t*
surface_label_provenance_cache_create(void) {
    return calloc(1, sizeof(surface_label_provenance_cache_t));
}

static void surface_label_provenance_cache_destroy(
    surface_label_provenance_cache_t* cache) {
    if (!cache) return;
    free(cache->entries);
    free(cache);
}

static surface_label_provenance_cache_entry_t*
surface_label_provenance_cache_find(
    surface_label_provenance_cache_t* cache,
    const surface_label_candidate_t* candidate, int surface_id,
    int material_boundary) {
    if (!cache || !candidate || candidate->edge_x < 0) return NULL;
    for (size_t i = 0; i < cache->count; i++) {
        surface_label_provenance_cache_entry_t* entry = &cache->entries[i];
        if (entry->edge_x == candidate->edge_x &&
            entry->edge_y == candidate->edge_y &&
            entry->orientation == candidate->edge_orientation &&
            fabs(entry->fraction - candidate->edge_fraction) <= 1e-4 &&
            entry->surface_id == surface_id &&
            entry->material_boundary == material_boundary)
            return entry;
    }
    if (cache->count == cache->capacity) {
        size_t capacity = cache->capacity ? cache->capacity * 2 : 64;
        surface_label_provenance_cache_entry_t* entries = realloc(
            cache->entries, capacity * sizeof(*entries));
        if (!entries) return NULL;
        cache->entries = entries;
        cache->capacity = capacity;
    }
    surface_label_provenance_cache_entry_t* entry = &cache->entries[cache->count++];
    *entry = (surface_label_provenance_cache_entry_t){
        .edge_x = candidate->edge_x, .edge_y = candidate->edge_y,
        .orientation = candidate->edge_orientation,
        .fraction = candidate->edge_fraction, .surface_id = surface_id,
        .group = -1,
        .material_boundary = material_boundary
    };
    return entry;
}

static void boundary_trace_add_id(boundary_trace_t* trace, int surface_id) {
    if (surface_id <= 0) return;
    for (size_t i = 0; i < trace->count; i++)
        if (trace->ids[i] == surface_id) return;
    if (trace->count == sizeof(trace->ids) / sizeof(trace->ids[0])) {
        trace->saw_unresolved = 1;
        return;
    }
    trace->ids[trace->count++] = surface_id;
}

static void boundary_trace_group_add_id(boundary_trace_group_t* group,
                                        int surface_id) {
    if (surface_id <= 0) return;
    for (size_t i = 0; i < group->count; i++)
        if (group->ids[i] == surface_id) return;
    if (group->count < sizeof(group->ids) / sizeof(group->ids[0]))
        group->ids[group->count++] = surface_id;
}

static void slice_world_point(const alea_slice_view_t* view, double u, double v,
                              double out[3]) {
    for (int i = 0; i < 3; i++)
        out[i] = view->plane.origin[i] + u * view->plane.u_axis[i] +
                 v * view->plane.v_axis[i];
}

/* Reduce ordered ownership events to the physical crossing groups which
 * actually change the caller's displayed identity.  The event stream may
 * contain several coincident physical IDs at one parameter, or several
 * distinct transitions inside one coarse pixel edge.  Keeping the group count
 * separate from the flattened compatibility ID list prevents coincident
 * surfaces from being mistaken for spatially separate crossings. */
static int trace_boundary_event_groups(
    alea_system_t* sys, const alea_ray_boundary_event_t* events,
    size_t event_count, double t_offset, const double start[3],
    const double end[3], alea_slice_classify_point_fn classify,
    void* userdata, boundary_trace_t* out) {
    double dx = end[0] - start[0], dy = end[1] - start[1], dz = end[2] - start[2];
    double length = sqrt(dx * dx + dy * dy + dz * dz);
    if (!(length > 0.0)) return -1;
    dx /= length; dy /= length; dz /= length;

    const double endpoint_eps = length * 1e-8;
    for (size_t i = 0; i < event_count;) {
        const double t = events[i].t - t_offset;
        size_t group_end = i + 1;
        while (group_end < event_count &&
               fabs(events[group_end].t - events[i].t) <= RAY_EPSILON)
            group_end++;

        if (t <= endpoint_eps || t >= length - endpoint_eps) {
            i = group_end;
            continue;
        }

        int has_physical = 0;
        for (size_t k = i; k < group_end; k++) {
            if (events[k].kind == ALEA_RAY_BOUNDARY_EVENT_SYNTHETIC_LATTICE)
                out->saw_synthetic = 1;
            else if (events[k].kind == ALEA_RAY_BOUNDARY_EVENT_UNRESOLVED)
                out->saw_unresolved = 1;
            else if (events[k].surface_id > 0)
                has_physical = 1;
        }
        if (!has_physical) {
            i = group_end;
            continue;
        }

        /* Classify only inside the open interval adjacent to this group. A
         * fixed percentage of the edge can jump across a neighbouring thin
         * region and attribute its surface to this group. */
        double previous = 0.0, following = length;
        for (size_t k = i; k > 0;) {
            k--;
            double candidate = events[k].t - t_offset;
            if (candidate < t - RAY_EPSILON) {
                if (candidate > 0.0) previous = candidate;
                break;
            }
        }
        for (size_t k = group_end; k < event_count; k++) {
            double candidate = events[k].t - t_offset;
            if (candidate > t + RAY_EPSILON) {
                if (candidate < length) following = candidate;
                break;
            }
        }
        double eps = fmin(length * 5e-2,
                          fmin(t - previous, following - t) * 0.5);
        if (eps <= endpoint_eps) {
            out->saw_unresolved = 1;
            i = group_end;
            continue;
        }
        double before[3], after[3];
        for (int c = 0; c < 3; c++) {
            double direction = c == 0 ? dx : c == 1 ? dy : dz;
            before[c] = start[c] + direction * (t - eps);
            after[c] = start[c] + direction * (t + eps);
        }
        alea_slice_classification_t a = {0}, b = {0};
        if (classify(sys, before[0], before[1], before[2], userdata, &a) != 0 ||
            classify(sys, after[0], after[1], after[2], userdata, &b) != 0) {
            out->saw_unresolved = 1;
            i = group_end;
            continue;
        }
        if (a.status == ALEA_SLICE_SAMPLE_UNRESOLVED ||
            b.status == ALEA_SLICE_SAMPLE_UNRESOLVED) {
            out->saw_unresolved = 1;
            i = group_end;
            continue;
        }

        /* A material contour intentionally ignores a cell transition whose
         * two projected material IDs are equal.  The same comparison also
         * makes custom display classifiers authoritative. */
        /* The exterior is a real displayed identity (zero), not an
         * indeterminate gap.  Treating it as such preserves outer CSG
         * boundaries while still discarding a 0-to-0 traversal. */
        int a_identity = a.status == ALEA_SLICE_SAMPLE_GAP ? 0 : a.identity;
        int b_identity = b.status == ALEA_SLICE_SAMPLE_GAP ? 0 : b.identity;
        if ((a.status == ALEA_SLICE_SAMPLE_SINGLE ||
             a.status == ALEA_SLICE_SAMPLE_GAP) &&
            (b.status == ALEA_SLICE_SAMPLE_SINGLE ||
             b.status == ALEA_SLICE_SAMPLE_GAP) && a_identity == b_identity) {
            i = group_end;
            continue;
        }
        if (a.status == ALEA_SLICE_SAMPLE_OVERLAP ||
            b.status == ALEA_SLICE_SAMPLE_OVERLAP)
            out->saw_overlap = 1;

        if (out->group_count == sizeof(out->groups) / sizeof(out->groups[0])) {
            out->saw_unresolved = 1;
            i = group_end;
            continue;
        }
        boundary_trace_group_t* group = &out->groups[out->group_count++];
        group->fraction = t / length;
        for (size_t k = i; k < group_end; k++)
            if (events[k].kind == ALEA_RAY_BOUNDARY_EVENT_PHYSICAL) {
                boundary_trace_add_id(out, events[k].surface_id);
                boundary_trace_group_add_id(group, events[k].surface_id);
            }
        i = group_end;
    }
    return 0;
}

static int trace_boundary_from_cached_events(
    alea_system_t* sys, const alea_slice_directional_event_cache_t* cache,
    alea_slice_edge_orientation_t orientation, int reverse, int line,
    double t_first, const double start[3], const double end[3],
    alea_slice_classify_point_fn classify, void* userdata,
    boundary_trace_t* out) {
    const alea_ray_boundary_event_t* events;
    size_t event_count;
    if (line < 0 || alea_slice_directional_event_cache_line_events(
            cache, orientation, reverse, (size_t)line, &events, &event_count) != 0)
        return -1;
    return trace_boundary_event_groups(sys, events, event_count, t_first,
                                       start, end, classify, userdata, out);
}

/* The long directional streams are an acceleration only.  This is the
 * canonical short-edge contract used when a stream is incomplete or the two
 * directions disagree: collect physical/synthetic events on this exact edge
 * and classify immediately around each crossing. */
static int trace_boundary_short_canonical(
    alea_system_t* sys, const double start[3], const double end[3],
    alea_slice_classify_point_fn classify, void* userdata,
    boundary_trace_t* out) {
    double dx = end[0] - start[0], dy = end[1] - start[1], dz = end[2] - start[2];
    double length = sqrt(dx * dx + dy * dy + dz * dz);
    alea_raycast_result_t trace;
    alea_ray_boundary_event_result_t event_result;
    alea_ray_boundary_event_options_internal_t options = {
        .include_all_coincident_physical = true
    };
    if (!(length > 0.0)) return -1;
    dx /= length; dy /= length; dz /= length;
    alea_ray_t ray;
    alea_ray_init_normalized(&ray, start[0], start[1], start[2], dx, dy, dz);
    alea_raycast_result_init(&trace);
    alea_ray_boundary_event_result_init(&event_result);
    if (alea_raycast_boundary_events_with_options(sys, &ray, length, &options,
                                                  &trace, &event_result) != 0) {
        alea_ray_boundary_event_result_free(&event_result);
        alea_raycast_result_free(&trace);
        return -1;
    }
    int rc = trace_boundary_event_groups(sys, event_result.events.data,
                                         event_result.events.count, 0.0,
                                         start, end, classify, userdata, out);
    alea_ray_boundary_event_result_free(&event_result);
    alea_raycast_result_free(&trace);
    return rc;
}

/* Convert one compact boundary-event batch row back to the canonical event
 * representation consumed by the shared ownership-group reducer. The sparse
 * caller bounds the batch event count before this allocation is reachable. */
static int trace_boundary_from_batch_row(
    alea_system_t* sys, const alea_ray_boundary_event_batch_result_t* batch,
    size_t row, const double start[3], const double end[3],
    alea_slice_classify_point_fn classify, void* userdata,
    boundary_trace_t* out) {
    if (!batch || row >= batch->ray_count || !batch->ray_offsets) return -1;
    uint64_t begin64 = batch->ray_offsets[row];
    uint64_t end64 = batch->ray_offsets[row + 1];
    if (end64 < begin64 || end64 > SIZE_MAX) return -1;
    size_t count = (size_t)(end64 - begin64);
    alea_ray_boundary_event_t* events = NULL;
    if (count) {
        events = malloc(count * sizeof(*events));
        if (!events) return -1;
        size_t begin = (size_t)begin64;
        for (size_t i = 0; i < count; i++) {
            const size_t source = begin + i;
            events[i] = (alea_ray_boundary_event_t){
                .t = batch->t[source],
                .kind = (alea_ray_boundary_event_kind_t)batch->kinds[source],
                .surface_id = batch->surface_ids[source],
                .primitive_id = batch->primitive_ids
                    ? batch->primitive_ids[source] : UINT32_MAX,
                .cell_before = batch->cell_before[source],
                .cell_after = batch->cell_after[source],
                .material_before = batch->material_before[source],
                .material_after = batch->material_after[source],
                .resolution_flags = batch->resolution_flags[source]
            };
        }
    }
    int rc = trace_boundary_event_groups(sys, events, count, 0.0, start, end,
                                         classify, userdata, out);
    free(events);
    return rc;
}

static int boundary_trace_same_ids(const boundary_trace_t* a,
                                   const boundary_trace_t* b);

static int trace_boundary_shared_root_surface(
    alea_system_t* sys, int cell_a, int cell_b,
    const double start[3], const double end[3], boundary_trace_t* forward,
    boundary_trace_t* reverse) {
    alea_hier_cell_hit_t hit_a, hit_b;
    alea_hier_ray_path_t path_a, path_b;
    alea_ray_t ray;
    double dx = end[0] - start[0], dy = end[1] - start[1], dz = end[2] - start[2];
    const double length = sqrt(dx * dx + dy * dy + dz * dz);
    int surface_id = -1;
    double fraction = 0.0;
    if (!(length > 0.0) ||
        alea_hier_spatial_find_path_at_point(
            sys, start[0], start[1], start[2], &hit_a, &path_a) <= 0 ||
        alea_hier_spatial_find_path_at_point(
            sys, end[0], end[1], end[2], &hit_b, &path_b) <= 0 ||
        path_a.count <= 0 || path_b.count <= 0)
        return 0;
    g_sparse_surface_label_stats.local_path_pairs_resolved++;
    const alea_hier_ray_path_entry_t* terminal_a = &path_a.entries[path_a.count - 1];
    const alea_hier_ray_path_entry_t* terminal_b = &path_b.entries[path_b.count - 1];
    if (terminal_a->cell_id != cell_a || terminal_b->cell_id != cell_b ||
        terminal_a->universe_id != terminal_b->universe_id)
        return 0;
    g_sparse_surface_label_stats.local_path_pairs_same_universe++;
    if (terminal_a->depth != terminal_b->depth ||
        memcmp(terminal_a->transform.inv, terminal_b->transform.inv,
               sizeof(terminal_a->transform.inv)) != 0)
        return 0;
    g_sparse_surface_label_stats.local_path_pairs_same_transform++;
    double local_start[3] = {start[0], start[1], start[2]};
    alea_matrix_transform_point_inverse(&terminal_a->transform,
                                        &local_start[0], &local_start[1],
                                        &local_start[2]);
    const double local_dx = terminal_a->transform.inv[0] * dx +
        terminal_a->transform.inv[1] * dy + terminal_a->transform.inv[2] * dz;
    const double local_dy = terminal_a->transform.inv[4] * dx +
        terminal_a->transform.inv[5] * dy + terminal_a->transform.inv[6] * dz;
    const double local_dz = terminal_a->transform.inv[8] * dx +
        terminal_a->transform.inv[9] * dy + terminal_a->transform.inv[10] * dz;
    const double local_length = sqrt(local_dx * local_dx + local_dy * local_dy +
                                     local_dz * local_dz);
    if (!(local_length > 0.0)) return 0;
    if (alea_ray_init(&ray, local_start[0], local_start[1], local_start[2],
                      local_dx, local_dy, local_dz) != 0 ||
        alea_raycast_shared_terminal_surface_nocache(
            sys, cell_a, cell_b, &ray, local_length, &surface_id, &fraction) != 1)
        return 0;
    *forward = (boundary_trace_t){ .ids = {surface_id}, .count = 1,
        .group_count = 1, .groups = {{.fraction = fraction,
                                      .ids = {surface_id}, .count = 1}} };
    *reverse = (boundary_trace_t){ .ids = {surface_id}, .count = 1,
        .group_count = 1, .groups = {{.fraction = 1.0 - fraction,
                                      .ids = {surface_id}, .count = 1}} };
    return 1;
}

static int trace_boundary_one_sided_surface(
    alea_system_t* sys, int cell_a, int cell_b, const double start[3],
    const double end[3], boundary_trace_t* forward, boundary_trace_t* reverse) {
    alea_hier_cell_hit_t hit;
    alea_hier_ray_path_t path_a, path_b;
    const double dx = end[0] - start[0], dy = end[1] - start[1], dz = end[2] - start[2];
    const double world_length = sqrt(dx * dx + dy * dy + dz * dz);
    int selected_cell = -1; const double* selected_point = NULL;
    int first = alea_hier_spatial_find_path_at_point(sys, start[0], start[1], start[2], &hit, &path_a);
    if (first > 0 && path_a.count > 0) { selected_cell = cell_a; selected_point = start; }
    int second = alea_hier_spatial_find_path_at_point(sys, end[0], end[1], end[2], &hit, &path_b);
    if ((first > 0) == (second > 0) || !(world_length > 0.0) || selected_cell < 0)
        return 0;
    const alea_hier_ray_path_t* path = &path_a;
    if (second > 0) { selected_cell = cell_b; selected_point = end; path = &path_b; }
    if (path->count <= 0 || path->entries[path->count - 1].cell_id != selected_cell) return 0;
    const alea_hier_ray_path_entry_t* terminal = &path->entries[path->count - 1];
    double local[3] = {selected_point[0], selected_point[1], selected_point[2]};
    alea_matrix_transform_point_inverse(&terminal->transform, &local[0], &local[1], &local[2]);
    double ldx = terminal->transform.inv[0]*dx + terminal->transform.inv[1]*dy + terminal->transform.inv[2]*dz;
    double ldy = terminal->transform.inv[4]*dx + terminal->transform.inv[5]*dy + terminal->transform.inv[6]*dz;
    double ldz = terminal->transform.inv[8]*dx + terminal->transform.inv[9]*dy + terminal->transform.inv[10]*dz;
    double local_length = sqrt(ldx*ldx + ldy*ldy + ldz*ldz), fraction; int surface; alea_ray_t ray;
    if (!(local_length > 0.0) || alea_ray_init(&ray, local[0], local[1], local[2], ldx, ldy, ldz) != 0 ||
        alea_raycast_terminal_surface_nocache(sys, selected_cell, &ray, local_length, &surface, &fraction) != 1)
        return 0;
    if (selected_point == end) fraction = 1.0 - fraction;
    *forward = (boundary_trace_t){.ids={surface},.count=1,.group_count=1,.groups={{.fraction=fraction,.ids={surface},.count=1}}};
    *reverse = (boundary_trace_t){.ids={surface},.count=1,.group_count=1,.groups={{.fraction=1.0-fraction,.ids={surface},.count=1}}};
    return 1;
}

/* The sparse curve path uses this exact-edge check after it has selected a
 * candidate.  It deliberately shares the full boundary map's physical-group
 * reduction and adjacent-point identity test: surface membership alone is not
 * enough for a visible material boundary, nor for a coincident ambiguous
 * crossing. */
static int label_edge_has_surface_canonical(
    alea_system_t* sys, const double start[3], const double end[3],
    double fraction, double fraction_tolerance, int surface_id,
    int material_boundary, surface_label_candidate_t* candidate,
    surface_label_provenance_cache_t* cache) {
    surface_label_provenance_cache_entry_t* cached =
        surface_label_provenance_cache_find(cache, candidate, surface_id,
                                            material_boundary);
    if (cached && cached->verified) {
        if (cached->verified > 0 && candidate)
            candidate->provenance_group = cached->group;
        return cached->verified > 0;
    }
    boundary_trace_t forward = {0}, reverse = {0};
    alea_slice_classify_point_fn classify = material_boundary
        ? alea_slice_classify_material : alea_slice_classify_cell;
    int valid = trace_boundary_short_canonical(sys, start, end, classify,
                                                   NULL, &forward) == 0 &&
        trace_boundary_short_canonical(sys, end, start, classify, NULL,
                                           &reverse) == 0 &&
        !forward.saw_synthetic && !reverse.saw_synthetic &&
        !forward.saw_gap && !reverse.saw_gap &&
        !forward.saw_unresolved && !reverse.saw_unresolved &&
        boundary_trace_same_ids(&forward, &reverse);
    if (!valid) {
        if (cached) cached->verified = -1;
        return 0;
    }
    for (size_t g = 0; g < forward.group_count; g++) {
        if (fabs(forward.groups[g].fraction - fraction) > fraction_tolerance)
            continue;
        for (size_t i = 0; i < forward.groups[g].count; i++)
            if (forward.groups[g].ids[i] == surface_id) {
                if (candidate) candidate->provenance_group = (int)g;
                if (cached) {
                    cached->verified = 1;
                    cached->group = (int)g;
                }
                return 1;
            }
    }
    if (cached) cached->verified = -1;
    return 0;
}

static int boundary_trace_same_ids(const boundary_trace_t* a,
                                   const boundary_trace_t* b) {
    if (a->count != b->count || a->group_count != b->group_count) return 0;
    for (size_t i = 0; i < a->count; i++) {
        int found = 0;
        for (size_t j = 0; j < b->count; j++)
            if (a->ids[i] == b->ids[j]) { found = 1; break; }
        if (!found) return 0;
    }
    /* Reverse tracing reports fractions from the opposite endpoint and visits
     * groups in reverse order. Agreement is a causal-group contract, not just
     * an equality check on the flattened surface-ID summary. */
    for (size_t i = 0; i < a->group_count; i++) {
        const boundary_trace_group_t* forward = &a->groups[i];
        const boundary_trace_group_t* reverse =
            &b->groups[b->group_count - 1 - i];
        if (forward->count != reverse->count ||
            fabs(forward->fraction - (1.0 - reverse->fraction)) > RAY_EPSILON)
            return 0;
        for (size_t s = 0; s < forward->count; s++) {
            int found = 0;
            for (size_t r = 0; r < reverse->count; r++)
                if (forward->ids[s] == reverse->ids[r]) {
                    found = 1;
                    break;
                }
            if (!found) return 0;
        }
    }
    return 1;
}

/* One changed rendered edge selected by the bounded sparse-grid sampler. */
typedef struct {
    int x, y;
    int orientation;
    int identity_lo, identity_hi;
    size_t support;
    double centre_distance;
} sparse_grid_edge_t;

typedef struct {
    alea_label_position_t label;
    size_t support;
} sparse_grid_observation_t;

typedef struct {
    alea_label_position_t label;
    size_t observations;
    size_t support;
    double centre_distance;
} sparse_grid_ranked_label_t;

static int compare_sparse_grid_edge_pair(const void* lhs, const void* rhs) {
    const sparse_grid_edge_t* a = lhs;
    const sparse_grid_edge_t* b = rhs;
    if (a->identity_lo != b->identity_lo)
        return a->identity_lo < b->identity_lo ? -1 : 1;
    if (a->identity_hi != b->identity_hi)
        return a->identity_hi < b->identity_hi ? -1 : 1;
    if (a->centre_distance != b->centre_distance)
        return a->centre_distance < b->centre_distance ? -1 : 1;
    if (a->y != b->y) return a->y < b->y ? -1 : 1;
    if (a->x != b->x) return a->x < b->x ? -1 : 1;
    return a->orientation - b->orientation;
}

static int sparse_grid_edge_better(const sparse_grid_edge_t* a,
                                   const sparse_grid_edge_t* b) {
    if (a->support != b->support) return a->support > b->support;
    if (a->centre_distance != b->centre_distance)
        return a->centre_distance < b->centre_distance;
    if (a->identity_lo != b->identity_lo)
        return a->identity_lo < b->identity_lo;
    return a->identity_hi < b->identity_hi;
}

static int compare_sparse_grid_observation(const void* lhs, const void* rhs) {
    const sparse_grid_observation_t* a = lhs;
    const sparse_grid_observation_t* b = rhs;
    if (a->label.id != b->label.id)
        return a->label.id < b->label.id ? -1 : 1;
    if (a->label.py != b->label.py)
        return a->label.py < b->label.py ? -1 : 1;
    if (a->label.px != b->label.px)
        return a->label.px < b->label.px ? -1 : 1;
    if (a->label.provenance_orientation != b->label.provenance_orientation)
        return a->label.provenance_orientation - b->label.provenance_orientation;
    return a->label.provenance_group - b->label.provenance_group;
}

static int compare_sparse_grid_ranked_label(const void* lhs, const void* rhs) {
    const sparse_grid_ranked_label_t* a = lhs;
    const sparse_grid_ranked_label_t* b = rhs;
    if (a->observations != b->observations)
        return a->observations > b->observations ? -1 : 1;
    if (a->support != b->support)
        return a->support > b->support ? -1 : 1;
    if (a->centre_distance != b->centre_distance)
        return a->centre_distance < b->centre_distance ? -1 : 1;
    if (a->label.id == b->label.id) return 0;
    return a->label.id < b->label.id ? -1 : 1;
}

static int sparse_grid_append_observation(
    sparse_grid_observation_t** observations, size_t* count, size_t* capacity,
    int surface_id, const sparse_grid_edge_t* edge, int group) {
    if (*count == *capacity) {
        size_t next = *capacity ? *capacity * 2 : 256;
        sparse_grid_observation_t* grown = realloc(
            *observations, next * sizeof(**observations));
        if (!grown) return -1;
        *observations = grown;
        *capacity = next;
    }
    (*observations)[(*count)++] = (sparse_grid_observation_t){
        .label = {
            .id = surface_id,
            .px = edge->x,
            .py = edge->y,
            .pixel_count = edge->support > INT_MAX ? INT_MAX : (int)edge->support,
            .provenance_edge_x = edge->x,
            .provenance_edge_y = edge->y,
            .provenance_orientation = edge->orientation,
            .provenance_group = group
        },
        .support = edge->support
    };
    return 0;
}

int alea_find_surface_labels_sparse_on_grid(
    alea_system_t* sys, const alea_slice_view_t* view,
    int width, int height, const int* grid_ids,
    alea_slice_classify_point_fn classify, void* classify_userdata,
    int margin, size_t max_queries, size_t max_labels,
    alea_label_position_t** out_labels, int* out_count) {
    if (!sys || !view || !grid_ids || !classify || !out_labels || !out_count ||
        width <= 0 || height <= 0 || margin < 0 ||
        max_queries == 0 || max_labels == 0)
        return -1;
    memset(&g_sparse_surface_label_stats, 0,
           sizeof(g_sparse_surface_label_stats));
    *out_labels = NULL;
    *out_count = 0;
    /* Build mutable query caches before independent edge traces enter the
     * parallel region; lazy first-use preparation from several workers would
     * violate the public query-cache lifecycle. */
    if (alea_prepare_query_acceleration(sys) != 0) return -1;

    /* Use up to two dominant displayed-identity pairs per tile.  The tile size
     * is increased until the number of canonical edge queries is bounded by
     * max_queries, independent of grid resolution and model surface count. */
    const size_t slots_per_tile = max_queries >= 2 ? 2 : 1;
    const size_t target_tiles = max_queries / slots_per_tile;
    int tile_size = 1;
    for (;;) {
        size_t tile_cols = ((size_t)width + tile_size - 1) / tile_size;
        size_t tile_rows = ((size_t)height + tile_size - 1) / tile_size;
        if (tile_cols <= target_tiles / tile_rows) break;
        tile_size++;
    }
    size_t tile_cols = ((size_t)width + tile_size - 1) / tile_size;
    size_t tile_rows = ((size_t)height + tile_size - 1) / tile_size;
    size_t candidate_capacity = tile_cols * tile_rows * slots_per_tile;
    if (candidate_capacity > max_queries) candidate_capacity = max_queries;
    sparse_grid_edge_t* candidates = malloc(
        candidate_capacity * sizeof(*candidates));
    size_t local_capacity = (size_t)tile_size * tile_size * 2;
    sparse_grid_edge_t* local = malloc(local_capacity * sizeof(*local));
    if (!candidates || !local) {
        free(candidates); free(local);
        return -1;
    }

    size_t candidate_count = 0;
    for (size_t tile_y = 0; tile_y < tile_rows; tile_y++) {
        int y0 = (int)(tile_y * tile_size);
        int y1 = y0 + tile_size;
        if (y1 > height) y1 = height;
        for (size_t tile_x = 0; tile_x < tile_cols; tile_x++) {
            g_sparse_surface_label_stats.tiles_examined++;
            int x0 = (int)(tile_x * tile_size);
            int x1 = x0 + tile_size;
            if (x1 > width) x1 = width;
            double centre_x = 0.5 * (x0 + x1 - 1);
            double centre_y = 0.5 * (y0 + y1 - 1);
            size_t local_count = 0;
            for (int y = y0; y < y1; y++) {
                for (int x = x0; x < x1; x++) {
                    if (x < margin || x >= width - margin ||
                        y < margin || y >= height - margin)
                        continue;
                    for (int orient = ALEA_SLICE_EDGE_RIGHT;
                         orient <= ALEA_SLICE_EDGE_DOWN; orient++) {
                        int nx = x + (orient == ALEA_SLICE_EDGE_RIGHT);
                        int ny = y - (orient == ALEA_SLICE_EDGE_DOWN);
                        if (nx < 0 || nx >= width || ny < 0 || ny >= height)
                            continue;
                        int a = grid_ids[(size_t)y * width + x];
                        int b = grid_ids[(size_t)ny * width + nx];
                        if (a == b) continue;
                        g_sparse_surface_label_stats.changed_edges++;
                        if (local_count == local_capacity) {
                            free(local); free(candidates);
                            return -1;
                        }
                        double dx = x - centre_x, dy = y - centre_y;
                        local[local_count++] = (sparse_grid_edge_t){
                            .x = x, .y = y, .orientation = orient,
                            .identity_lo = a < b ? a : b,
                            .identity_hi = a < b ? b : a,
                            .centre_distance = dx * dx + dy * dy
                        };
                    }
                }
            }
            if (local_count == 0) continue;
            qsort(local, local_count, sizeof(*local),
                  compare_sparse_grid_edge_pair);

            sparse_grid_edge_t best[2];
            size_t best_count = 0;
            for (size_t first = 0; first < local_count;) {
                size_t last = first + 1;
                while (last < local_count &&
                       local[last].identity_lo == local[first].identity_lo &&
                       local[last].identity_hi == local[first].identity_hi)
                    last++;
                sparse_grid_edge_t group = local[first];
                group.support = last - first;
                size_t pos = best_count;
                if (pos < slots_per_tile) {
                    best[best_count++] = group;
                } else {
                    pos = slots_per_tile - 1;
                    if (sparse_grid_edge_better(&group, &best[pos]))
                        best[pos] = group;
                    else {
                        first = last;
                        continue;
                    }
                }
                while (pos > 0 && sparse_grid_edge_better(&best[pos], &best[pos - 1])) {
                    sparse_grid_edge_t swap = best[pos - 1];
                    best[pos - 1] = best[pos];
                    best[pos] = swap;
                    pos--;
                }
                first = last;
            }
            for (size_t i = 0;
                 i < best_count && candidate_count < candidate_capacity; i++)
                candidates[candidate_count++] = best[i];
        }
    }
    free(local);
    g_sparse_surface_label_stats.candidate_edges = candidate_count;

    sparse_grid_observation_t* observations = NULL;
    size_t observation_count = 0, observation_capacity = 0;
    int observation_failed = 0;
    double du = (view->u_max - view->u_min) / width;
    double dv = (view->v_max - view->v_min) / height;
    /* Batch both directions for each selected edge. The batch owns one
     * reusable trace/event pair per capped worker, avoiding per-edge rich
     * result churn while retaining the exact scalar event contract. Its fixed
     * event cap falls back to the scalar path rather than publishing a partial
     * provenance answer. */
    boundary_trace_t* batch_forward = NULL;
    boundary_trace_t* batch_reverse = NULL;
    int batch_ready = 0;
    if (candidate_count && candidate_count <= SIZE_MAX / 6 &&
        candidate_count <= SIZE_MAX / (2 * sizeof(*batch_forward))) {
        g_sparse_surface_label_stats.batch_attempts++;
        const size_t ray_count = candidate_count * 2;
        double* origins = calloc(ray_count * 3, sizeof(*origins));
        double* directions = calloc(ray_count * 3, sizeof(*directions));
        double* t_maxs = calloc(ray_count, sizeof(*t_maxs));
        batch_forward = calloc(candidate_count, sizeof(*batch_forward));
        batch_reverse = calloc(candidate_count, sizeof(*batch_reverse));
        if (origins && directions && t_maxs && batch_forward && batch_reverse) {
            for (size_t i = 0; i < candidate_count; i++) {
                const sparse_grid_edge_t* edge = &candidates[i];
                int nx = edge->x + (edge->orientation == ALEA_SLICE_EDGE_RIGHT);
                int ny = edge->y - (edge->orientation == ALEA_SLICE_EDGE_DOWN);
                double start[3], end[3];
                slice_world_point(view, view->u_min + (edge->x + 0.5) * du,
                                  view->v_min + (edge->y + 0.5) * dv, start);
                slice_world_point(view, view->u_min + (nx + 0.5) * du,
                                  view->v_min + (ny + 0.5) * dv, end);
                double dx = end[0] - start[0], dy = end[1] - start[1];
                double dz = end[2] - start[2];
                double length = sqrt(dx * dx + dy * dy + dz * dz);
                if (!(length > 0.0)) continue;
                const size_t forward_row = i * 2, reverse_row = forward_row + 1;
                for (int c = 0; c < 3; c++) {
                    const double direction = (end[c] - start[c]) / length;
                    origins[forward_row * 3 + c] = start[c];
                    directions[forward_row * 3 + c] = direction;
                    origins[reverse_row * 3 + c] = end[c];
                    directions[reverse_row * 3 + c] = -direction;
                }
                t_maxs[forward_row] = length;
                t_maxs[reverse_row] = length;
            }
            const alea_ray_batch_query_t query = {
                .kind = ALEA_RAY_QUERY_BOUNDARY_EVENTS,
                .material_filter = -1,
                .t_maxs = t_maxs,
                .max_events = ray_count <= UINT64_MAX / 256 ? ray_count * 256 : 0,
                .max_output_bytes = 8u * 1024u * 1024u,
                .max_workers = 1,
                .include_all_coincident_physical = true
            };
            alea_ray_boundary_event_batch_result_t batch;
            alea_ray_boundary_event_batch_result_init(&batch);
            if (alea_raycast_boundary_events_batch_nocache(
                    sys, origins, directions, ray_count, &query, &batch) == 0) {
                batch_ready = 1;
                g_sparse_surface_label_stats.breakpoint_hits =
                    batch.breakpoint_hits;
                g_sparse_surface_label_stats.selected_segments =
                    batch.selected_segments;
                for (size_t i = 0; i < candidate_count && batch_ready; i++) {
                    const sparse_grid_edge_t* edge = &candidates[i];
                    int nx = edge->x + (edge->orientation == ALEA_SLICE_EDGE_RIGHT);
                    int ny = edge->y - (edge->orientation == ALEA_SLICE_EDGE_DOWN);
                    double start[3], end[3];
                    slice_world_point(view, view->u_min + (edge->x + 0.5) * du,
                                      view->v_min + (edge->y + 0.5) * dv, start);
                    slice_world_point(view, view->u_min + (nx + 0.5) * du,
                                      view->v_min + (ny + 0.5) * dv, end);
                    if (trace_boundary_from_batch_row(
                            sys, &batch, i * 2, start, end, classify,
                            classify_userdata, &batch_forward[i]) != 0 ||
                        trace_boundary_from_batch_row(
                            sys, &batch, i * 2 + 1, end, start, classify,
                            classify_userdata, &batch_reverse[i]) != 0)
                        batch_ready = 0;
                }
            }
            alea_ray_boundary_event_batch_result_free(&batch);
        }
        free(origins); free(directions); free(t_maxs);
        if (!batch_ready) {
            free(batch_forward); free(batch_reverse);
            batch_forward = NULL; batch_reverse = NULL;
        } else {
            g_sparse_surface_label_stats.batch_traces_used = ray_count;
        }
    }
    for (long long candidate_index = 0;
         candidate_index < (long long)candidate_count; candidate_index++) {
        size_t i = (size_t)candidate_index;
        const sparse_grid_edge_t* edge = &candidates[i];
        int nx = edge->x + (edge->orientation == ALEA_SLICE_EDGE_RIGHT);
        int ny = edge->y - (edge->orientation == ALEA_SLICE_EDGE_DOWN);
        double u0 = view->u_min + (edge->x + 0.5) * du;
        double v0 = view->v_min + (edge->y + 0.5) * dv;
        double u1 = view->u_min + (nx + 0.5) * du;
        double v1 = view->v_min + (ny + 0.5) * dv;
        double start[3], end[3];
        slice_world_point(view, u0, v0, start);
        slice_world_point(view, u1, v1, end);
        boundary_trace_t forward = {0}, reverse = {0};
        int forward_rc = 0, reverse_rc = 0;
        const int cell_a = grid_ids[(size_t)edge->y * (size_t)width + edge->x];
        const int cell_b = grid_ids[(size_t)ny * (size_t)width + nx];
        /* Endpoint identities may bracket several selected intervals.  They
         * are therefore never sufficient provenance for a coarse grid edge;
         * attribution begins with the complete canonical trace until the
         * selected-walker group certificate is available. */
        (void)cell_a;
        (void)cell_b;
        if (batch_ready) {
            forward = batch_forward[i];
            reverse = batch_reverse[i];
            g_sparse_surface_label_stats.forward_trace_calls++;
            g_sparse_surface_label_stats.reverse_trace_calls++;
        } else {
            g_sparse_surface_label_stats.forward_trace_calls++;
            forward_rc = trace_boundary_short_canonical(
                sys, start, end, classify, classify_userdata, &forward);
            reverse_rc = -1;
            if (forward_rc == 0) {
                g_sparse_surface_label_stats.reverse_trace_calls++;
                reverse_rc = trace_boundary_short_canonical(
                    sys, end, start, classify, classify_userdata, &reverse);
            }
        }
        int valid = forward_rc == 0 && reverse_rc == 0 &&
            !forward.saw_synthetic && !reverse.saw_synthetic &&
            !forward.saw_gap && !reverse.saw_gap &&
            !forward.saw_overlap && !reverse.saw_overlap &&
            !forward.saw_unresolved && !reverse.saw_unresolved &&
            forward.group_count != 0 &&
            boundary_trace_same_ids(&forward, &reverse);
        if (!valid) continue;
        g_sparse_surface_label_stats.accepted_edges++;
        /* Ray attribution is independent per edge.  Only the compact accepted
         * observations share storage, so keep that append in a short critical
         * section while the expensive canonical traces run in parallel. */
        #pragma omp critical(sparse_surface_label_observations)
        {
            if (!observation_failed) {
                for (size_t group = 0; group < forward.group_count; group++) {
                    for (size_t participant = 0;
                         participant < forward.groups[group].count; participant++) {
                        if (sparse_grid_append_observation(
                                &observations, &observation_count,
                                &observation_capacity,
                                forward.groups[group].ids[participant], edge,
                                (int)group) != 0) {
                            observation_failed = 1;
                            break;
                        }
                    }
                    if (observation_failed) break;
                }
            }
        }
    }
    free(candidates);
    free(batch_forward);
    free(batch_reverse);
    if (observation_failed) {
        free(observations);
        return -1;
    }
    if (observation_count == 0) {
        free(observations);
        return 0;
    }
    g_sparse_surface_label_stats.observations = observation_count;

    qsort(observations, observation_count, sizeof(*observations),
          compare_sparse_grid_observation);
    sparse_grid_ranked_label_t* ranked = malloc(
        observation_count * sizeof(*ranked));
    if (!ranked) { free(observations); return -1; }
    size_t ranked_count = 0;
    double grid_cx = 0.5 * (width - 1), grid_cy = 0.5 * (height - 1);
    for (size_t first = 0; first < observation_count;) {
        size_t last = first + 1;
        while (last < observation_count &&
               observations[last].label.id == observations[first].label.id)
            last++;
        size_t best = first, total_support = 0;
        double best_distance = DBL_MAX;
        for (size_t i = first; i < last; i++) {
            double dx = observations[i].label.px - grid_cx;
            double dy = observations[i].label.py - grid_cy;
            double distance = dx * dx + dy * dy;
            if (observations[i].support > observations[best].support ||
                (observations[i].support == observations[best].support &&
                 distance < best_distance)) {
                best = i;
                best_distance = distance;
            }
            if (SIZE_MAX - total_support < observations[i].support)
                total_support = SIZE_MAX;
            else
                total_support += observations[i].support;
        }
        ranked[ranked_count++] = (sparse_grid_ranked_label_t){
            .label = observations[best].label,
            .observations = last - first,
            .support = total_support,
            .centre_distance = best_distance
        };
        first = last;
    }
    free(observations);
    qsort(ranked, ranked_count, sizeof(*ranked),
          compare_sparse_grid_ranked_label);
    if (ranked_count > max_labels) ranked_count = max_labels;
    alea_label_position_t* labels = malloc(ranked_count * sizeof(*labels));
    if (!labels) { free(ranked); return -1; }
    for (size_t i = 0; i < ranked_count; i++) {
        labels[i] = ranked[i].label;
        labels[i].pixel_count = ranked[i].support > INT_MAX
            ? INT_MAX : (int)ranked[i].support;
    }
    free(ranked);
    *out_labels = labels;
    *out_count = (int)ranked_count;
    g_sparse_surface_label_stats.labels = ranked_count;
    return 0;
}

static int boundary_map_append_ids(alea_slice_surface_boundary_map_t* map,
                                   const boundary_trace_t* a,
                                   const boundary_trace_t* b,
                                   size_t edge_start) {
    for (size_t pass = 0; pass < 2; pass++) {
        const boundary_trace_t* trace = pass == 0 ? a : b;
        for (size_t i = 0; i < trace->count; i++) {
            int duplicate = 0;
            for (size_t j = edge_start; j < map->surface_count; j++) {
                if (map->surface_ids[j] == trace->ids[i]) { duplicate = 1; break; }
            }
            if (duplicate) continue;
            if (map->surface_count == map->surface_capacity) {
                size_t capacity = map->surface_capacity ?
                    map->surface_capacity * 2 : 64;
                int* next = realloc(map->surface_ids,
                                    capacity * sizeof(*next));
                if (!next) return -1;
                map->surface_ids = next;
                map->surface_capacity = capacity;
            }
            map->surface_ids[map->surface_count++] = trace->ids[i];
        }
    }
    return 0;
}

static int boundary_map_append_groups(alea_slice_surface_boundary_map_t* map,
                                      const boundary_trace_t* trace) {
    for (size_t g = 0; g < trace->group_count; g++) {
        if (map->group_count == map->group_capacity) {
            size_t capacity = map->group_capacity ? map->group_capacity * 2 : 64;
            double* fractions = realloc(map->group_fractions,
                                        capacity * sizeof(*fractions));
            if (!fractions) return -1;
            size_t* offsets = realloc(map->group_surface_offsets,
                                      (capacity + 1) * sizeof(*offsets));
            if (!offsets) {
                map->group_fractions = fractions;
                return -1;
            }
            map->group_fractions = fractions;
            map->group_surface_offsets = offsets;
            map->group_capacity = capacity;
        }
        map->group_surface_offsets[map->group_count] = map->group_surface_count;
        map->group_fractions[map->group_count] = trace->groups[g].fraction;
        for (size_t i = 0; i < trace->groups[g].count; i++) {
            if (map->group_surface_count == map->group_surface_capacity) {
                size_t capacity = map->group_surface_capacity ?
                    map->group_surface_capacity * 2 : 64;
                int* next_ids = realloc(map->group_surface_ids,
                                        capacity * sizeof(*next_ids));
                if (!next_ids) return -1;
                map->group_surface_ids = next_ids;
                map->group_surface_capacity = capacity;
            }
            map->group_surface_ids[map->group_surface_count++] =
                trace->groups[g].ids[i];
        }
        map->group_count++;
        map->group_surface_offsets[map->group_count] = map->group_surface_count;
    }
    return 0;
}

static int slice_classify_builtin(alea_system_t* sys,
                                  double x, double y, double z,
                                  int material_identity, int universe_depth,
                                  alea_slice_classification_t* out) {
    if (!sys || !out) return -1;
    out->status = ALEA_SLICE_SAMPLE_UNRESOLVED;
    out->identity = -1;
    out->secondary_identity = -1;
    alea_cell_hit_t hits[64];
    int count = alea_find_all_cells_at_point(sys, x, y, z, hits, 64);
    if (count < 0 || count >= (int)(sizeof(hits) / sizeof(hits[0]))) return -1;
    if (count == 0) {
        out->status = ALEA_SLICE_SAMPLE_GAP;
        return 0;
    }
    int target_depth = -1;
    if (universe_depth < 0) {
        target_depth = hits[count - 1].depth;
    } else {
        /* Match alea_find_cells_grid: select the requested hierarchy depth
         * when present, otherwise the deepest available owner no deeper than
         * that request. */
        for (int i = 0; i < count; i++)
            if (hits[i].depth <= universe_depth &&
                (target_depth < 0 || hits[i].depth > target_depth))
                target_depth = hits[i].depth;
    }
    if (target_depth < 0) {
        out->status = ALEA_SLICE_SAMPLE_GAP;
        return 0;
    }
    int at_depth = 0;
    for (int i = 0; i < count; i++) {
        if (hits[i].depth != target_depth) continue;
        int identity = material_identity ? hits[i].material_id : hits[i].cell_id;
        if (at_depth == 0) out->identity = identity;
        else if (at_depth == 1) out->secondary_identity = identity;
        at_depth++;
    }
    if (at_depth == 0) {
        out->status = ALEA_SLICE_SAMPLE_GAP;
    } else if (at_depth == 1) {
        out->status = ALEA_SLICE_SAMPLE_SINGLE;
    } else {
        out->status = ALEA_SLICE_SAMPLE_OVERLAP;
    }
    return 0;
}

int alea_slice_classify_cell(alea_system_t* sys, double x, double y, double z,
                             void* userdata, alea_slice_classification_t* out) {
    (void)userdata;
    return slice_classify_builtin(sys, x, y, z, 0, -1, out);
}

int alea_slice_classify_material(alea_system_t* sys, double x, double y, double z,
                                 void* userdata, alea_slice_classification_t* out) {
    (void)userdata;
    return slice_classify_builtin(sys, x, y, z, 1, -1, out);
}

int alea_slice_classify_cell_at_depth(
    alea_system_t* sys, double x, double y, double z, void* userdata,
    alea_slice_classification_t* out) {
    int universe_depth = userdata ? *(const int*)userdata : -1;
    return slice_classify_builtin(sys, x, y, z, 0, universe_depth, out);
}

int alea_slice_classify_material_at_depth(
    alea_system_t* sys, double x, double y, double z, void* userdata,
    alea_slice_classification_t* out) {
    int universe_depth = userdata ? *(const int*)userdata : -1;
    return slice_classify_builtin(sys, x, y, z, 1, universe_depth, out);
}

static int slice_surface_boundary_map_create_impl(
    alea_system_t* sys, const alea_slice_view_t* view,
    int width, int height, const int* grid_ids,
    alea_slice_classify_point_fn classify, void* classify_userdata,
    const alea_slice_directional_event_cache_t* supplied_cache,
    alea_slice_surface_boundary_map_t** out_map) {
    if (!sys || !view || !grid_ids || !classify || !out_map ||
        width <= 0 || height <= 0)
        return -1;
    *out_map = NULL;
    size_t pixels = (size_t)width * (size_t)height;
    if (pixels > (SIZE_MAX / 2) - 1) return -1;
    alea_slice_surface_boundary_map_t* map = calloc(1, sizeof(*map));
    if (!map) return -1;
    map->width = width;
    map->height = height;
    map->status = calloc(pixels * 2, sizeof(*map->status));
    map->surface_offsets = calloc(pixels * 2 + 1, sizeof(*map->surface_offsets));
    map->group_offsets = calloc(pixels * 2 + 1, sizeof(*map->group_offsets));
    if (!map->status || !map->surface_offsets || !map->group_offsets) {
        alea_slice_surface_boundary_map_free(map);
        return -1;
    }

    double du = (view->u_max - view->u_min) / width;
    double dv = (view->v_max - view->v_min) / height;
    alea_slice_directional_event_cache_t* owned_cache = NULL;
    const alea_slice_directional_event_cache_t* cache = supplied_cache;
    if (cache) {
        if (!alea_slice_directional_event_cache_matches(cache, sys, view, width, height)) {
            alea_slice_surface_boundary_map_free(map);
            return -1;
        }
    } else {
        owned_cache = alea_slice_directional_event_cache_create(sys, view, width, height);
        if (!owned_cache) { alea_slice_surface_boundary_map_free(map); return -1; }
        cache = owned_cache;
    }
    for (int orient = ALEA_SLICE_EDGE_RIGHT;
         orient <= ALEA_SLICE_EDGE_DOWN; orient++) {
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                size_t edge = (size_t)orient * pixels + (size_t)y * width + x;
                map->surface_offsets[edge] = map->surface_count;
                map->group_offsets[edge] = map->group_count;
                int nx = x + (orient == ALEA_SLICE_EDGE_RIGHT);
                int ny = y - (orient == ALEA_SLICE_EDGE_DOWN);
                if (nx < 0 || nx >= width || ny < 0 || ny >= height ||
                    grid_ids[y * width + x] == grid_ids[ny * width + nx])
                    continue;
                double u0 = view->u_min + (x + 0.5) * du;
                double v0 = view->v_min + (y + 0.5) * dv;
                double u1 = view->u_min + (nx + 0.5) * du;
                double v1 = view->v_min + (ny + 0.5) * dv;
                double start[3], end[3];
                slice_world_point(view, u0, v0, start);
                slice_world_point(view, u1, v1, end);
                boundary_trace_t forward = {0}, reverse = {0};
                int line = orient == ALEA_SLICE_EDGE_RIGHT ? y : x;
                double forward_t = orient == ALEA_SLICE_EDGE_RIGHT ? x * du
                    : (height - 1 - y) * dv;
                double reverse_t = orient == ALEA_SLICE_EDGE_RIGHT
                    ? (width - 2 - x) * du : (y - 1) * dv;
                int forward_rc = trace_boundary_from_cached_events(
                    sys, cache, orient, 0, line, forward_t, start, end,
                    classify, classify_userdata, &forward);
                int reverse_rc = trace_boundary_from_cached_events(
                    sys, cache, orient, 1, line, reverse_t, end, start,
                    classify, classify_userdata, &reverse);
                if (forward_rc != 0 || reverse_rc != 0) {
                    forward.saw_unresolved = 1;
                }
                if (forward_rc != 0 || reverse_rc != 0 ||
                    forward.saw_unresolved || reverse.saw_unresolved ||
                    !boundary_trace_same_ids(&forward, &reverse)) {
                    /* Discard suspect long-line evidence; the precise short
                     * edge is the reporting authority for this map entry. */
                    forward = (boundary_trace_t){0};
                    reverse = (boundary_trace_t){0};
                    if (trace_boundary_short_canonical(sys, start, end, classify,
                                                       classify_userdata, &forward) != 0 ||
                        trace_boundary_short_canonical(sys, end, start, classify,
                                                       classify_userdata, &reverse) != 0)
                        forward.saw_unresolved = 1;
                }
                if (boundary_map_append_ids(map, &forward, &reverse,
                                            map->surface_offsets[edge]) != 0) {
                    alea_slice_directional_event_cache_destroy(owned_cache);
                    alea_slice_surface_boundary_map_free(map);
                    return -1;
                }
                if (boundary_map_append_groups(map, &forward) != 0) {
                    alea_slice_directional_event_cache_destroy(owned_cache);
                    alea_slice_surface_boundary_map_free(map);
                    return -1;
                }
                if (forward.saw_gap || reverse.saw_gap)
                    map->status[edge] = ALEA_SLICE_BOUNDARY_GAP;
                else if (forward.saw_overlap || reverse.saw_overlap)
                    map->status[edge] = ALEA_SLICE_BOUNDARY_OVERLAP;
                else if (forward.saw_unresolved || reverse.saw_unresolved)
                    map->status[edge] = ALEA_SLICE_BOUNDARY_UNRESOLVED;
                else if (!boundary_trace_same_ids(&forward, &reverse))
                    map->status[edge] = ALEA_SLICE_BOUNDARY_AMBIGUOUS;
                else if (forward.group_count > 1)
                    map->status[edge] = ALEA_SLICE_BOUNDARY_MULTI_HIT;
                else if (forward.group_count == 1)
                    map->status[edge] = ALEA_SLICE_BOUNDARY_VALID;
                else if (forward.saw_synthetic || reverse.saw_synthetic)
                    map->status[edge] = ALEA_SLICE_BOUNDARY_SYNTHETIC;
                else
                    map->status[edge] = ALEA_SLICE_BOUNDARY_UNRESOLVED;
        }
        }
    }
    map->surface_offsets[pixels * 2] = map->surface_count;
    map->group_offsets[pixels * 2] = map->group_count;
    alea_slice_directional_event_cache_destroy(owned_cache);
    *out_map = map;
    return 0;
}

int alea_slice_surface_boundary_map_create(
    alea_system_t* sys, const alea_slice_view_t* view,
    int width, int height, const int* grid_ids,
    alea_slice_classify_point_fn classify, void* classify_userdata,
    alea_slice_surface_boundary_map_t** out_map) {
    return slice_surface_boundary_map_create_impl(
        sys, view, width, height, grid_ids, classify, classify_userdata,
        NULL, out_map);
}

int alea_slice_surface_boundary_map_create_with_event_cache(
    alea_system_t* sys, const alea_slice_view_t* view, int width, int height,
    const int* grid_ids, alea_slice_classify_point_fn classify,
    void* classify_userdata, const alea_slice_directional_event_cache_t* cache,
    alea_slice_surface_boundary_map_t** out_map) {
    return slice_surface_boundary_map_create_impl(
        sys, view, width, height, grid_ids, classify, classify_userdata,
        cache, out_map);
}

int alea_slice_surface_boundary_map_create_with_directional_cache(
    alea_system_t* sys, const alea_slice_view_t* view, int width, int height,
    const int* grid_ids, alea_slice_classify_point_fn classify,
    void* classify_userdata, const alea_slice_directional_trace_cache_t* cache,
    alea_slice_surface_boundary_map_t** out_map) {
    return alea_slice_surface_boundary_map_create_with_event_cache(
        sys, view, width, height, grid_ids, classify, classify_userdata,
        cache, out_map);
}

void alea_slice_surface_boundary_map_free(alea_slice_surface_boundary_map_t* map) {
    if (!map) return;
    free(map->status);
    free(map->surface_offsets);
    free(map->surface_ids);
    free(map->group_offsets);
    free(map->group_fractions);
    free(map->group_surface_offsets);
    free(map->group_surface_ids);
    free(map);
}

static int boundary_map_edge_index(const alea_slice_surface_boundary_map_t* map,
                                   int x, int y,
                                   alea_slice_edge_orientation_t orientation,
                                   size_t* out) {
    if (!map || x < 0 || x >= map->width || y < 0 || y >= map->height ||
        (orientation != ALEA_SLICE_EDGE_RIGHT &&
         orientation != ALEA_SLICE_EDGE_DOWN))
        return -1;
    *out = (size_t)orientation * (size_t)map->width * map->height +
           (size_t)y * map->width + x;
    return 0;
}

alea_slice_boundary_status_t alea_slice_surface_boundary_status(
    const alea_slice_surface_boundary_map_t* map, int x, int y,
    alea_slice_edge_orientation_t orientation) {
    size_t edge;
    return boundary_map_edge_index(map, x, y, orientation, &edge) == 0
        ? map->status[edge] : ALEA_SLICE_BOUNDARY_UNRESOLVED;
}

size_t alea_slice_surface_boundary_surface_count(
    const alea_slice_surface_boundary_map_t* map, int x, int y,
    alea_slice_edge_orientation_t orientation) {
    size_t edge;
    if (boundary_map_edge_index(map, x, y, orientation, &edge) != 0) return 0;
    return map->surface_offsets[edge + 1] - map->surface_offsets[edge];
}

int alea_slice_surface_boundary_surface_id(
    const alea_slice_surface_boundary_map_t* map, int x, int y,
    alea_slice_edge_orientation_t orientation, size_t index) {
    size_t edge;
    if (boundary_map_edge_index(map, x, y, orientation, &edge) != 0 ||
        index >= map->surface_offsets[edge + 1] - map->surface_offsets[edge])
        return -1;
    return map->surface_ids[map->surface_offsets[edge] + index];
}

size_t alea_slice_surface_boundary_group_count(
    const alea_slice_surface_boundary_map_t* map, int x, int y,
    alea_slice_edge_orientation_t orientation) {
    size_t edge;
    if (boundary_map_edge_index(map, x, y, orientation, &edge) != 0) return 0;
    return map->group_offsets[edge + 1] - map->group_offsets[edge];
}

static int boundary_map_group_index(const alea_slice_surface_boundary_map_t* map,
                                    int x, int y,
                                    alea_slice_edge_orientation_t orientation,
                                    size_t group_index, size_t* out) {
    size_t edge;
    if (boundary_map_edge_index(map, x, y, orientation, &edge) != 0 ||
        group_index >= map->group_offsets[edge + 1] - map->group_offsets[edge])
        return -1;
    *out = map->group_offsets[edge] + group_index;
    return 0;
}

double alea_slice_surface_boundary_group_fraction(
    const alea_slice_surface_boundary_map_t* map, int x, int y,
    alea_slice_edge_orientation_t orientation, size_t group_index) {
    size_t group;
    return boundary_map_group_index(map, x, y, orientation, group_index,
                                    &group) == 0
        ? map->group_fractions[group] : -1.0;
}

size_t alea_slice_surface_boundary_group_surface_count(
    const alea_slice_surface_boundary_map_t* map, int x, int y,
    alea_slice_edge_orientation_t orientation, size_t group_index) {
    size_t group;
    if (boundary_map_group_index(map, x, y, orientation, group_index,
                                 &group) != 0)
        return 0;
    return map->group_surface_offsets[group + 1] -
           map->group_surface_offsets[group];
}

int alea_slice_surface_boundary_group_surface_id(
    const alea_slice_surface_boundary_map_t* map, int x, int y,
    alea_slice_edge_orientation_t orientation, size_t group_index,
    size_t surface_index) {
    size_t group;
    if (boundary_map_group_index(map, x, y, orientation, group_index,
                                 &group) != 0 ||
        surface_index >= map->group_surface_offsets[group + 1] -
                         map->group_surface_offsets[group])
        return -1;
    return map->group_surface_ids[map->group_surface_offsets[group] +
                                  surface_index];
}

static int boundary_status_is_labelable(alea_slice_boundary_status_t status) {
    return status == ALEA_SLICE_BOUNDARY_VALID ||
           status == ALEA_SLICE_BOUNDARY_MULTI_HIT;
}

/* One participation of a physical surface in a labelable boundary edge,
 * reduced to the pixel that carries it. */
typedef struct {
    int id;
    int x, y;
} boundary_label_point_t;

static int compare_boundary_label_points(const void* a, const void* b) {
    const boundary_label_point_t* p = (const boundary_label_point_t*)a;
    const boundary_label_point_t* q = (const boundary_label_point_t*)b;
    if (p->id != q->id) return p->id < q->id ? -1 : 1;
    if (p->y != q->y) return p->y < q->y ? -1 : 1;
    if (p->x != q->x) return p->x < q->x ? -1 : 1;
    return 0;
}

/* Locate (x,y) inside the already sorted [lo,hi) range of a single surface. */
static int boundary_label_find_point(const boundary_label_point_t* points,
                                     size_t lo, size_t hi, int x, int y,
                                     size_t* out_index) {
    size_t end = hi;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (points[mid].y < y || (points[mid].y == y && points[mid].x < x)) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    if (lo >= end || points[lo].x != x || points[lo].y != y) return 0;
    *out_index = lo;
    return 1;
}

static size_t boundary_label_root(size_t* parent, size_t i) {
    while (parent[i] != i) {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    return i;
}

static void boundary_label_join(size_t* parent, size_t a, size_t b) {
    a = boundary_label_root(parent, a);
    b = boundary_label_root(parent, b);
    if (a != b) parent[b] = a;
}

static int boundary_label_collect_points(
    const alea_slice_surface_boundary_map_t* map,
    boundary_label_point_t** out_points, size_t* out_count) {
    size_t pixels = (size_t)map->width * (size_t)map->height;
    size_t edges = pixels * 2;
    boundary_label_point_t* points = NULL;
    size_t count = 0, capacity = 0;

    for (size_t edge = 0; edge < edges; edge++) {
        if (!boundary_status_is_labelable(map->status[edge])) continue;
        int x = (int)((edge % pixels) % (size_t)map->width);
        int y = (int)((edge % pixels) / (size_t)map->width);
        for (size_t si = map->surface_offsets[edge];
             si < map->surface_offsets[edge + 1]; si++) {
            if (count == capacity) {
                size_t next_capacity = capacity ? capacity * 2 : 256;
                boundary_label_point_t* next =
                    realloc(points, next_capacity * sizeof(*next));
                if (!next) { free(points); return -1; }
                points = next;
                capacity = next_capacity;
            }
            points[count].id = map->surface_ids[si];
            points[count].x = x;
            points[count].y = y;
            count++;
        }
    }

    *out_points = points;
    *out_count = count;
    return 0;
}

int alea_find_surface_labels_on_boundary_map(
    const alea_slice_surface_boundary_map_t* map, int margin,
    alea_label_position_t** out_labels, int* out_count) {
    if (!map || !out_labels || !out_count || margin < 0) return -1;
    *out_labels = NULL;
    *out_count = 0;

    /* Minimum edges in one connected arc before it earns a label. Applied per
     * arc rather than per surface so that a surface which appears as several
     * disjoint pieces is labelled on each of them, and so that an incidental
     * few-pixel sliver is labelled on none. */
    const size_t MIN_ARC_EDGE_COUNT = 30;

    boundary_label_point_t* points = NULL;
    size_t count = 0;
    if (boundary_label_collect_points(map, &points, &count) != 0)
        return -1;
    if (count == 0) { free(points); return 0; }

    qsort(points, count, sizeof(*points), compare_boundary_label_points);

    /* A surface owning both edges of one pixel contributes it only once. */
    size_t unique = 1;
    for (size_t i = 1; i < count; i++) {
        if (compare_boundary_label_points(&points[i], &points[unique - 1]) != 0)
            points[unique++] = points[i];
    }
    count = unique;

    size_t* parent = malloc(count * sizeof(*parent));
    double* sum_x = calloc(count, sizeof(*sum_x));
    double* sum_y = calloc(count, sizeof(*sum_y));
    size_t* members = calloc(count, sizeof(*members));
    double* best_distance = malloc(count * sizeof(*best_distance));
    size_t* best_point = malloc(count * sizeof(*best_point));
    if (!parent || !sum_x || !sum_y || !members || !best_distance || !best_point) {
        free(points); free(parent); free(sum_x); free(sum_y);
        free(members); free(best_distance); free(best_point);
        return -1;
    }
    for (size_t i = 0; i < count; i++) parent[i] = i;

    /* Group the points of each surface into connected arcs. Only neighbours
     * that precede a point in sort order need to be joined; the reverse
     * direction is covered when that neighbour is visited. */
    static const int NEIGHBOR_DX[4] = { -1, -1, 0, 1 };
    static const int NEIGHBOR_DY[4] = { 0, -1, -1, -1 };
    size_t group_lo = 0;
    while (group_lo < count) {
        size_t group_hi = group_lo;
        int id = points[group_lo].id;
        while (group_hi < count && points[group_hi].id == id) group_hi++;
        for (size_t i = group_lo; i < group_hi; i++) {
            for (int d = 0; d < 4; d++) {
                size_t j;
                if (boundary_label_find_point(points, group_lo, group_hi,
                                              points[i].x + NEIGHBOR_DX[d],
                                              points[i].y + NEIGHBOR_DY[d], &j))
                    boundary_label_join(parent, i, j);
            }
        }
        group_lo = group_hi;
    }

    for (size_t i = 0; i < count; i++) {
        size_t root = boundary_label_root(parent, i);
        sum_x[root] += points[i].x;
        sum_y[root] += points[i].y;
        members[root]++;
    }

    /* Place each arc's label on the arc point nearest its centroid. The margin
     * constrains where the text may sit, but arcs are formed and measured
     * without it: clipping first would split one arc crossing the margin into
     * two, and label it twice. */
    for (size_t i = 0; i < count; i++) best_distance[i] = DBL_MAX;
    for (size_t i = 0; i < count; i++) {
        size_t root = boundary_label_root(parent, i);
        if (members[root] < MIN_ARC_EDGE_COUNT) continue;
        if (points[i].x < margin || points[i].x >= map->width - margin ||
            points[i].y < margin || points[i].y >= map->height - margin)
            continue;
        double dx = points[i].x - sum_x[root] / (double)members[root];
        double dy = points[i].y - sum_y[root] / (double)members[root];
        double distance = dx * dx + dy * dy;
        if (distance < best_distance[root]) {
            best_distance[root] = distance;
            best_point[root] = i;
        }
    }

    /* An arc entirely inside the margin has no place to put its label. */
    size_t label_count = 0;
    for (size_t i = 0; i < count; i++)
        if (boundary_label_root(parent, i) == i &&
            members[i] >= MIN_ARC_EDGE_COUNT && best_distance[i] < DBL_MAX)
            label_count++;
    if (label_count == 0) {
        free(points); free(parent); free(sum_x); free(sum_y);
        free(members); free(best_distance); free(best_point);
        return 0;
    }

    alea_label_position_t* labels = calloc(label_count, sizeof(*labels));
    if (!labels) {
        free(points); free(parent); free(sum_x); free(sum_y);
        free(members); free(best_distance); free(best_point);
        return -1;
    }

    size_t emitted = 0;
    for (size_t i = 0; i < count; i++) {
        if (boundary_label_root(parent, i) != i) continue;
        if (members[i] < MIN_ARC_EDGE_COUNT || best_distance[i] == DBL_MAX) continue;
        size_t chosen = best_point[i];
        labels[emitted++] = (alea_label_position_t){
            .id = points[chosen].id,
            .px = points[chosen].x,
            .py = points[chosen].y,
            .pixel_count = (int)members[i]
        };
    }

    free(points); free(parent); free(sum_x); free(sum_y);
    free(members); free(best_distance); free(best_point);
    *out_labels = labels;
    *out_count = (int)emitted;
    return 0;
}
