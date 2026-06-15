// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file raycast.c
 * @brief Ray casting through CSG geometry
 *
 * Main raycast logic: find surface intersections, determine cell segments.
 */

#include "raycast.h"
#include "alea.h"
#include "ray_intersect.h"
#include "ray_epsilon.h"
#include "ray_bbox.h"
#include "bvh.h"
#include "core/alea_system.h"
#include "core/alea_universe.h"
#include "util/compat.h"
#include "core/alea_spatial_hier.h"
#include "core/alea_eval.h"
#include "util/alea_log.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "util/math.h"
#include <stdio.h>

#define INITIAL_CAPACITY 32
#define MAX_FILL_RAYCAST_DEPTH 32

/* ============================================================================
 * CACHE PRE-BUILD (thread safety)
 *
 * All shared caches (surface BVH, spatial index, cell adjacency) must be
 * built before any concurrent access.  Call once from each public entry
 * point so that OpenMP threads never race on lazy initialisation.
 * ============================================================================ */

int alea_raycast_ensure_caches(alea_system_t* sys) {
    return alea_system_prepare_query_caches(sys, ALEA_CACHE_RAYCAST);
}

int alea_raycast_ensure_hier_caches(alea_system_t* sys) {
    return alea_system_prepare_query_caches(sys, ALEA_CACHE_RAYCAST_HIER);
}

static bool raycast_prefers_hier_mode(const alea_system_t* sys) {
    (void)sys;
    return true;  /* hierarchical spatial index is the only backend */
}

static int raycast_cell_aware_impl(alea_system_t* sys,
                                   const alea_ray_t* ray,
                                   double effective_t_max,
                                   bool use_hier_lookup,
                                   bool emit_hits,
                                   alea_raycast_result_t* result);

/* ============================================================================
 * RAY UTILITIES
 * ============================================================================ */

int alea_ray_init(alea_ray_t* ray,
                  double ox, double oy, double oz,
                  double dx, double dy, double dz) {
    ray->ox = ox;
    ray->oy = oy;
    ray->oz = oz;

    /* Normalize direction */
    double len = sqrt(dx * dx + dy * dy + dz * dz);
    if (len > RAY_EPSILON) {
        ray->dx = dx / len;
        ray->dy = dy / len;
        ray->dz = dz / len;
    } else {
        /* Zero-length direction: set fallback so struct is valid */
        ray->dx = 0;
        ray->dy = 0;
        ray->dz = 1;
    }

    /* IEEE 754: 1.0/0.0 = +inf, 1.0/(-0.0) = -inf — correct for slab tests */
    ray->inv_dx = 1.0 / ray->dx;
    ray->inv_dy = 1.0 / ray->dy;
    ray->inv_dz = 1.0 / ray->dz;

    return (len > RAY_EPSILON) ? 0 : -1;
}

void alea_ray_init_normalized(alea_ray_t* ray,
                              double ox, double oy, double oz,
                              double dx, double dy, double dz) {
    ray->ox = ox; ray->oy = oy; ray->oz = oz;
    ray->dx = dx; ray->dy = dy; ray->dz = dz;
    ray->inv_dx = 1.0 / dx;
    ray->inv_dy = 1.0 / dy;
    ray->inv_dz = 1.0 / dz;
}

/* ============================================================================
 * RESULT MANAGEMENT
 * ============================================================================ */

void alea_raycast_result_init(alea_raycast_result_t* result) {
    memset(result, 0, sizeof(*result));
}

/* Zero every per-ray counter (everything in the result except the ray and the
 * hit/segment vectors). Shared by free and clear so new counters only need to
 * be added in one place. */
static void raycast_result_reset_counters(alea_raycast_result_t* result) {
    result->surfaces_tested = 0;
    result->bbox_culled = 0;
    result->point_lookups = 0;
    result->step_iterations = 0;
    result->blas_placement_candidates = 0;
    result->blas_placements_pruned = 0;
    result->blas_universe_queries = 0;
    result->blas_cell_candidates = 0;
    result->blas_cells_tested = 0;
    result->blas_hits_before_dedup = 0;
    result->terminal_surfaces_tested = 0;
    result->lattice_surfaces_tested = 0;
    result->ancestor_surfaces_tested = 0;
    result->crossed_cell_count = 0;
    result->max_cell_surface_count = 0;
    result->sum_cell_surface_count = 0;
    memset(result->prim_type_tests, 0, sizeof(result->prim_type_tests));
}

void alea_raycast_result_free(alea_raycast_result_t* result) {
    alea_vec_free(&result->hits);
    alea_vec_free(&result->segments);
    raycast_result_reset_counters(result);
}

void alea_raycast_result_clear(alea_raycast_result_t* result) {
    alea_vec_clear(&result->hits);
    alea_vec_clear(&result->segments);
    raycast_result_reset_counters(result);
}

void alea_raycast_result_reserve(alea_raycast_result_t* result,
                                size_t hit_cap, size_t seg_cap) {
    alea_vec_reserve(&result->hits, hit_cap, alea_ray_hit_t);
    alea_vec_reserve(&result->segments, seg_cap, alea_ray_segment_t);
}

static int add_hit(alea_raycast_result_t* result, const alea_ray_hit_t* hit) {
    int res = alea_vec_push(&result->hits, *hit, alea_ray_hit_t);
    return res != 0 ? -1 : 0;
}

static bool raycast_primitive_copy_payload(const alea_system_t* sys,
                                           uint32_t primitive_id,
                                           alea_primitive_type_t type,
                                           alea_primitive_data_t* out) {
    if (!out) return false;
    const void* payload = alea_primitive_payload_const(sys, primitive_id);
    if (!payload) return false;

    switch (type) {
        case ALEA_PRIMITIVE_PLANE:       out->plane = *(const alea_plane_data_t*)payload; return true;
        case ALEA_PRIMITIVE_SPHERE:      out->sphere = *(const alea_sphere_data_t*)payload; return true;
        case ALEA_PRIMITIVE_CYLINDER_X:  out->cyl_x = *(const alea_cylinder_x_data_t*)payload; return true;
        case ALEA_PRIMITIVE_CYLINDER_Y:  out->cyl_y = *(const alea_cylinder_y_data_t*)payload; return true;
        case ALEA_PRIMITIVE_CYLINDER_Z:  out->cyl_z = *(const alea_cylinder_z_data_t*)payload; return true;
        case ALEA_PRIMITIVE_CONE_X:      out->cone_x = *(const alea_cone_x_data_t*)payload; return true;
        case ALEA_PRIMITIVE_CONE_Y:      out->cone_y = *(const alea_cone_y_data_t*)payload; return true;
        case ALEA_PRIMITIVE_CONE_Z:      out->cone_z = *(const alea_cone_z_data_t*)payload; return true;
        case ALEA_PRIMITIVE_RPP:         out->box = *(const alea_box_data_t*)payload; return true;
        case ALEA_PRIMITIVE_QUADRIC:     out->quadric = *(const alea_quadric_data_t*)payload; return true;
        case ALEA_PRIMITIVE_TORUS_X:
        case ALEA_PRIMITIVE_TORUS_Y:
        case ALEA_PRIMITIVE_TORUS_Z:     out->torus = *(const alea_torus_data_t*)payload; return true;
        case ALEA_PRIMITIVE_RCC:         out->rcc = *(const alea_rcc_data_t*)payload; return true;
        case ALEA_PRIMITIVE_BOX:         out->box_general = *(const alea_box_general_data_t*)payload; return true;
        case ALEA_PRIMITIVE_SPH:         out->sph = *(const alea_sph_data_t*)payload; return true;
        case ALEA_PRIMITIVE_TRC:         out->trc = *(const alea_trc_data_t*)payload; return true;
        case ALEA_PRIMITIVE_ELL:         out->ell = *(const alea_ell_data_t*)payload; return true;
        case ALEA_PRIMITIVE_REC:         out->rec = *(const alea_rec_data_t*)payload; return true;
        case ALEA_PRIMITIVE_WED:         out->wed = *(const alea_wed_data_t*)payload; return true;
        case ALEA_PRIMITIVE_RHP:         out->rhp = *(const alea_rhp_data_t*)payload; return true;
        case ALEA_PRIMITIVE_ARB:         out->arb = *(const alea_arb_data_t*)payload; return true;
        default: return false;
    }
}

static int add_segment(alea_raycast_result_t* result, const alea_ray_segment_t* seg) {
    int res = alea_vec_push(&result->segments, *seg, alea_ray_segment_t);
    return res != 0 ? -1 : 0;
}

/* ============================================================================
 * COMPARISON FOR SORTING
 * ============================================================================ */

static int compare_hits(const void* a, const void* b) {
    const alea_ray_hit_t* ha = (const alea_ray_hit_t*)a;
    const alea_ray_hit_t* hb = (const alea_ray_hit_t*)b;
    if (ha->t < hb->t) return -1;
    if (ha->t > hb->t) return 1;
    /* When t values are equal, put DDA boundary hits (surface_id=0) LAST.
     * This ensures the segment builder correctly identifies intervals after
     * coincident surface+DDA hits as crossing a DDA boundary, not a surface.
     * Critical for correct lattice element transitions on all platforms. */
    if (ha->surface_id == 0 && hb->surface_id != 0) return 1;
    if (ha->surface_id != 0 && hb->surface_id == 0) return -1;
    return 0;
}

/* Inline comparison for insertion sort (avoids function-call overhead) */
static inline int hit_less(const alea_ray_hit_t* a, const alea_ray_hit_t* b) {
    if (a->t != b->t) return a->t < b->t;
    /* DDA boundary hits (surface_id=0) sort after real surface hits */
    return (a->surface_id != 0 && b->surface_id == 0);
}

/**
 * Sort hits by distance. Uses insertion sort for small arrays (<=64),
 * qsort for larger arrays. Insertion sort has lower overhead and is
 * cache-friendly for the small arrays typical of raycast results.
 */
static void sort_hits(alea_ray_hit_t* hits, size_t count) {
    if (count <= 1) return;

    if (count <= 64) {
        /* Insertion sort — O(n^2) but faster than qsort for small n */
        for (size_t i = 1; i < count; i++) {
            alea_ray_hit_t key = hits[i];
            size_t j = i;
            while (j > 0 && hit_less(&key, &hits[j - 1])) {
                hits[j] = hits[j - 1];
                j--;
            }
            hits[j] = key;
        }
    } else {
        qsort(hits, count, sizeof(alea_ray_hit_t), compare_hits);
    }
}

/* ray_bbox_slab and ray_bbox_slab_enter_exit are in ray_bbox.h */

/* ============================================================================
 * MAIN RAYCAST - SURFACE INTERSECTIONS
 * ============================================================================ */

/* Context for BVH traversal callback */
typedef struct {
    alea_system_t* sys;
    const alea_ray_t* ray;
    double t_min;
    double t_max;
    alea_raycast_result_t* result;
} bvh_raycast_ctx_t;

/* BVH batch callback: test surface intersections for a leaf node */
static void bvh_surface_batch_callback(const uint32_t* surface_indices,
                                        uint16_t count, void* userdata) {
    bvh_raycast_ctx_t* ctx = (bvh_raycast_ctx_t*)userdata;

    for (uint16_t si = 0; si < count; si++) {
        uint32_t surface_idx = surface_indices[si];
        const alea_surface_entry_t* surf = &ctx->sys->surfaces.data[surface_idx];
        const alea_primitive_entry_t* prim = &ctx->sys->primitives.data[surf->primitive_id];
        alea_primitive_data_t prim_data;
        if (!raycast_primitive_copy_payload(ctx->sys, surf->primitive_id,
                                            prim->type, &prim_data)) continue;

        ctx->result->surfaces_tested++;

        /* Find intersections */
        double t[4];
        int n = ray_intersect_primitive(ctx->ray, prim->type, &prim_data, t);

        /* Add valid hits */
        for (int j = 0; j < n; j++) {
            if (t[j] >= ctx->t_min && t[j] <= ctx->t_max) {
                alea_ray_hit_t hit;
                hit.t = t[j];
                hit.surface_id = surf->mc_surface_id;
                hit.primitive_id = surf->primitive_id;

                /* Compute normal at hit point */
                double px, py, pz;
                alea_ray_point_at(ctx->ray, t[j], &px, &py, &pz);
                primitive_normal_at(prim->type, &prim_data, px, py, pz,
                                   &hit.nx, &hit.ny, &hit.nz);

                if (add_hit(ctx->result, &hit) != 0) {
                    ALEA_LOG_WARN("add_hit failed (out of memory) - raycast results may be incomplete");
                    return;
                }
            }
        }
    }
}

/* Linear scan fallback (used when BVH not available) */
static int raycast_surfaces_linear(alea_system_t* sys,
                                   const alea_ray_t* ray,
                                   double t_min, double t_max,
                                   alea_raycast_result_t* result) {
    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        
        const alea_surface_entry_t* surf = &sys->surfaces.data[i];
        const alea_primitive_entry_t* prim = &sys->primitives.data[surf->primitive_id];
        alea_primitive_data_t prim_data;
        if (!raycast_primitive_copy_payload(sys, surf->primitive_id,
                                            prim->type, &prim_data)) continue;

        result->surfaces_tested++;

        double t[4];
        int count = ray_intersect_primitive(ray, prim->type, &prim_data, t);

        for (int j = 0; j < count; j++) {
            if (t[j] >= t_min && t[j] <= t_max) {
                alea_ray_hit_t hit;
                hit.t = t[j];
                hit.surface_id = surf->mc_surface_id;
                hit.primitive_id = surf->primitive_id;

                double px, py, pz;
                alea_ray_point_at(ray, t[j], &px, &py, &pz);
                primitive_normal_at(prim->type, &prim_data, px, py, pz,
                                   &hit.nx, &hit.ny, &hit.nz);

                if (add_hit(result, &hit) != 0) {
                    ALEA_LOG_WARN("add_hit failed (out of memory) - raycast results may be incomplete");
                    return -1;
                }
            }
        }
    }
    return 0;
}

/* Shared implementation: intersect, sort, dedup */
static int raycast_surfaces_impl(alea_system_t* sys,
                                  const alea_ray_t* ray,
                                  double t_min, double t_max,
                                  alea_raycast_result_t* result) {
    alea_raycast_result_clear(result);
    result->ray = *ray;

    /* Use BVH if available, otherwise fall back to linear scan */
    if (!BVH_DISABLED && sys->surface_bvh) {
        bvh_raycast_ctx_t ctx = {
            .sys = sys,
            .ray = ray,
            .t_min = t_min,
            .t_max = t_max,
            .result = result
        };
        alea_bvh_traverse_batch(sys->surface_bvh, ray, t_min, t_max,
                               bvh_surface_batch_callback, &ctx);
    } else {
        raycast_surfaces_linear(sys, ray, t_min, t_max, result);
    }

    /* Sort hits by distance */
    sort_hits(result->hits.data, result->hits.count);

    /* Remove duplicate hits (same t AND same surface_id within epsilon) */
    if (result->hits.count > 1) {
        size_t write = 1;
        for (size_t read = 1; read < result->hits.count; read++) {
            int same_t = fabs(result->hits.data[read].t - result->hits.data[write - 1].t) <= DEDUP_EPSILON;
            int same_surf = result->hits.data[read].surface_id == result->hits.data[write - 1].surface_id;
            if (!(same_t && same_surf)) {
                result->hits.data[write++] = result->hits.data[read];
            }
        }
        result->hits.count = write;
    }

    return 0;
}

int alea_raycast_surfaces(alea_system_t* sys,
                         const alea_ray_t* ray,
                         double t_min, double t_max,
                         alea_raycast_result_t* result) {
    if (!sys || !ray || !result) return -1;

    if (!alea_system_query_cache_ready(sys, ALEA_CACHE_RAYCAST)) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              "raycast caches are not prepared; call alea_prepare_query_acceleration()");
        return -1;
    }

    return raycast_surfaces_impl(sys, ray, t_min, t_max, result);
}

int alea_raycast_surfaces_nocache(alea_system_t* sys,
                                  const alea_ray_t* ray,
                                  double t_min, double t_max,
                                  alea_raycast_result_t* result) {
    if (!sys || !ray || !result) return -1;
    return raycast_surfaces_impl(sys, ray, t_min, t_max, result);
}

/* ============================================================================
 * CONVERT HITS TO CELL SEGMENTS (with cell adjacency optimization)
 * ============================================================================ */

/**
 * @brief Find cell via neighbor lookup (O(1) when adjacency is built)
 * @return 1 if found, 0 if full lookup needed
 */
static int raycast_find_neighbor(alea_system_t* sys,
                                 int current_cell_idx,
                                 int surface_id,
                                 int* out_cell_id, int* out_cell_idx,
                                 int* out_material_id, double* out_density) {
    if (!sys->cell_adjacency_built || current_cell_idx < 0) {
        return 0;
    }

    int neighbor_idx = alea_find_neighbor_cell(sys, (uint32_t)current_cell_idx, surface_id);
    if (neighbor_idx < 0) {
        return 0;  /* No neighbor (exterior/void boundary) */
    }

    const alea_cell_entry_t* neighbor = &sys->cells.data[neighbor_idx];
    *out_cell_id = neighbor->mc_cell_id;
    *out_cell_idx = neighbor_idx;
    *out_material_id = neighbor->material_id;
    *out_density = neighbor->density;
    return 1;
}

static int cell_references_surface_id(const alea_system_t* sys,
                                      const alea_cell_entry_t* cell,
                                      int surface_id) {
    if (!sys || !cell || surface_id <= 0 ||
        !cell->surface_indices || cell->surface_index_count == 0) {
        return 0;
    }

    for (size_t i = 0; i < cell->surface_index_count; i++) {
        uint32_t surf_idx = cell->surface_indices[i];
        if (surf_idx >= alea_vec_count(&sys->surfaces))
            continue;
        if (sys->surfaces.data[surf_idx].mc_surface_id == surface_id)
            return 1;
    }
    return 0;
}

/**
 * @brief Full cell lookup via point-in-cell search
 */
static void raycast_find_cell_full(alea_system_t* sys,
                                   double px, double py, double pz,
                                   int* out_cell_id, int* out_cell_idx,
                                   int* out_material_id, double* out_density) {
    *out_cell_id = -1;
    *out_cell_idx = -1;
    *out_material_id = 0;
    *out_density = 0;

    alea_cell_hit_t hits[32];
    int num_hits = alea_find_all_cells_at_point(sys, px, py, pz, hits, 32);

    if (num_hits > 0) {
        const alea_cell_hit_t* hit = &hits[num_hits - 1];
        *out_cell_id = hit->cell_id;
        *out_cell_idx = hit->cell_index;
        *out_material_id = hit->material_id;

        if (hit->cell_index >= 0 && (size_t)hit->cell_index < alea_vec_count(&sys->cells)) {
            *out_density = sys->cells.data[hit->cell_index].density;
        }
    }
}

static void raycast_find_cell_full_hier(alea_system_t* sys,
                                        double px, double py, double pz,
                                        int* out_cell_id, int* out_cell_idx,
                                        int* out_material_id, double* out_density) {
    *out_cell_id = -1;
    *out_cell_idx = -1;
    *out_material_id = 0;
    *out_density = 0;

    alea_cell_hit_t hits[32];
    int num_hits = alea_hier_spatial_find_cells_at_point(sys, px, py, pz,
                                                         hits, 32);

    if (num_hits > 0) {
        const alea_cell_hit_t* hit = &hits[num_hits - 1];
        *out_cell_id = hit->cell_id;
        *out_cell_idx = hit->cell_index;
        *out_material_id = hit->material_id;

        if (hit->cell_index >= 0 && (size_t)hit->cell_index < alea_vec_count(&sys->cells)) {
            *out_density = sys->cells.data[hit->cell_index].density;
        }
    }
}

/**
 * 3-tier cell lookup after crossing a surface:
 * 1) Neighbor lookup (O(1) via adjacency)
 * 2) Coherence check (O(tree_depth) — is point still in previous cell?)
 * 3) Full point-in-cell search (O(n_candidates) via spatial index)
 */
static void find_cell_after_crossing(alea_system_t* sys,
                                     const alea_ray_t* ray,
                                     double t_prev, double t_curr,
                                     int prev_cell_idx,
                                     int crossed_surface_id,
                                     bool use_hier_lookup,
                                     int* out_cell_id, int* out_cell_idx,
                                     int* out_material_id, double* out_density) {
    *out_cell_id = -1;
    *out_cell_idx = -1;
    *out_material_id = 0;
    *out_density = 0;
    int found = 0;

    /* Tier 1: neighbor lookup */
    if (prev_cell_idx >= 0 && crossed_surface_id > 0) {
        const alea_cell_entry_t* prev_cell = &sys->cells.data[prev_cell_idx];
        if (prev_cell->universe_id == 0 &&
            cell_references_surface_id(sys, prev_cell, crossed_surface_id)) {
            found = raycast_find_neighbor(sys, prev_cell_idx, crossed_surface_id,
                                          out_cell_id, out_cell_idx,
                                          out_material_id, out_density);
        }
    }

    int crossed_dda = (crossed_surface_id == 0);

    /* Tier 2: coherence check (skip for DDA boundaries) */
    if (!found && prev_cell_idx >= 0 && !crossed_dda) {
        const alea_cell_entry_t* prev_cell = &sys->cells.data[prev_cell_idx];

        /* Nested-cell case: in flat mode Tier 3 is already fast via the
         * flat point-query coherence cache, so skip Tier 2 there. In hier
         * mode we consult the hier cache for the cell's world→local
         * transform. */
        if (prev_cell->universe_id != 0 && !use_hier_lookup)
            goto full_lookup;

        double t_sample = t_prev + fmin(0.5 * (t_curr - t_prev), SURFACE_SAMPLE_OFFSET);
        double px, py, pz;
        alea_ray_point_at(ray, t_sample, &px, &py, &pz);

        int in_cell;
        if (prev_cell->universe_id == 0) {
            in_cell = alea_contains_point(sys, prev_cell->root_node_id,
                                          px, py, pz) ? 1 : 0;
        } else {
            in_cell = alea_hier_spatial_check_cached_containment(
                sys, (uint32_t)prev_cell_idx, px, py, pz);
            if (in_cell < 0) goto full_lookup;
        }

        if (in_cell) {
            *out_cell_id = prev_cell->mc_cell_id;
            *out_cell_idx = prev_cell_idx;
            *out_material_id = prev_cell->material_id;
            *out_density = prev_cell->density;
            found = 1;
        }
    }

    /* Tier 3: full lookup */
full_lookup:
    if (!found) {
        double max_offset = crossed_dda ? DDA_SAMPLE_OFFSET : SURFACE_SAMPLE_OFFSET;
        double t_sample = t_prev + fmin(0.5 * (t_curr - t_prev), max_offset);
        double px, py, pz;
        alea_ray_point_at(ray, t_sample, &px, &py, &pz);
        if (use_hier_lookup) {
            raycast_find_cell_full_hier(sys, px, py, pz,
                                        out_cell_id, out_cell_idx,
                                        out_material_id, out_density);
        } else {
            raycast_find_cell_full(sys, px, py, pz,
                                   out_cell_id, out_cell_idx,
                                   out_material_id, out_density);
        }
    }
}

static int raycast_to_segments_impl(alea_system_t* sys,
                                    double t_max,
                                    alea_raycast_result_t* result,
                                    bool use_hier_lookup) {
    if (!sys || !result) return -1;

    result->segments.count = 0;

    const alea_ray_t* ray = &result->ray;
    double t_prev = 0;
    int prev_cell_id = -2;  /* Use -2 as "no previous" since -1 is valid (void) */
    int prev_cell_idx = -1; /* Track cell index for neighbor lookup */
    double t_terminus = (t_max > 0) ? t_max : DBL_MAX;

    /* Process intervals between hits */
    for (size_t i = 0; i <= result->hits.count; i++) {
        double t_curr = (i < result->hits.count) ? result->hits.data[i].t : t_terminus;
        if (t_curr > t_terminus) t_curr = t_terminus;

        /* Only process if there's a real interval */
        if (t_curr > t_prev + RAY_EPSILON) {
            int cell_id, cell_idx, material_id;
            double density;
            int crossed_surface = (i > 0) ? result->hits.data[i - 1].surface_id : -1;
            find_cell_after_crossing(sys, ray, t_prev, t_curr,
                                     prev_cell_idx, crossed_surface,
                                     use_hier_lookup,
                                     &cell_id, &cell_idx, &material_id, &density);

            /* Update tracking for next iteration */
            prev_cell_idx = cell_idx;

            /* Extend previous segment or start new one */
            if (cell_id == prev_cell_id && result->segments.count > 0) {
                /* Extend */
                alea_ray_segment_t* prev_seg =
                    &result->segments.data[result->segments.count - 1];
                prev_seg->t_exit = t_curr;
                prev_seg->exit_surface_id =
                    (i < result->hits.count &&
                     result->hits.data[i].t <= t_terminus + RAY_EPSILON)
                        ? result->hits.data[i].surface_id
                        : -1;
            } else {
                /* New segment */
                alea_ray_segment_t seg;
                seg.t_enter = t_prev;
                seg.t_exit = t_curr;
                seg.cell_id = cell_id;
                seg.material_id = material_id;
                seg.density = density;
                seg.enter_surface_id = (i > 0) ? result->hits.data[i - 1].surface_id : -1;
                seg.exit_surface_id =
                    (i < result->hits.count &&
                     result->hits.data[i].t <= t_terminus + RAY_EPSILON)
                        ? result->hits.data[i].surface_id
                        : -1;
                seg.enter_hit_index = (i > 0) ? (int)(i - 1) : -1;

                add_segment(result, &seg);
                prev_cell_id = cell_id;
            }
        }

        t_prev = t_curr;
    }

    return 0;
}

int alea_raycast_to_segments(alea_system_t* sys,
                            double t_max,
                            alea_raycast_result_t* result) {
    return raycast_to_segments_impl(sys, t_max, result, false);
}

/* ============================================================================
 * LATTICE DDA RAYCAST
 * ============================================================================ */

/**
 * Walk a CSG tree and test ray intersection with every primitive leaf.
 * Adds hits to result. Duplicates from shared surfaces are handled by
 * the global dedup pass after all hits are collected.
 */
static void raycast_tree_primitives(alea_system_t* sys,
                                    const alea_ray_t* ray,
                                    alea_node_id_t node_id,
                                    double t_min, double t_max,
                                    alea_raycast_result_t* result) {
    if (node_id >= alea_vec_count(&sys->nodes)) return;
    const alea_node_t* node = &sys->nodes.data[node_id];
    alea_operation_t op = ALEA_GET_OPERATION(node);

    if (op == ALEA_OP_PRIMITIVE) {
        const alea_primitive_entry_t* prim =
            &sys->primitives.data[node->primitive.primitive_id];
        alea_primitive_data_t prim_data;
        if (!raycast_primitive_copy_payload(sys, node->primitive.primitive_id,
                                            prim->type, &prim_data)) return;
        double t[4];
        result->surfaces_tested++;
        int count = ray_intersect_primitive(ray, prim->type, &prim_data, t);
        for (int j = 0; j < count; j++) {
            if (t[j] >= t_min && t[j] <= t_max) {
                alea_ray_hit_t hit;
                hit.t = t[j];
                hit.surface_id = node->primitive.mc_surface_id;
                hit.primitive_id = node->primitive.primitive_id;
                double px, py, pz;
                alea_ray_point_at(ray, t[j], &px, &py, &pz);
                primitive_normal_at(prim->type, &prim_data, px, py, pz,
                                   &hit.nx, &hit.ny, &hit.nz);
                if (add_hit(result, &hit) != 0) {
                    ALEA_LOG_WARN("add_hit failed (out of memory) - raycast results may be incomplete");
                    return;
                }
            }
        }
        return;
    }

    if (op == ALEA_OP_COMPLEMENT) {
        raycast_tree_primitives(sys, ray, node->operation.left,
                                t_min, t_max, result);
    } else {
        raycast_tree_primitives(sys, ray, node->operation.left,
                                t_min, t_max, result);
        raycast_tree_primitives(sys, ray, node->operation.right,
                                t_min, t_max, result);
    }
}

static void raycast_cell_indexed_surface_hits(alea_system_t* sys,
                                              const alea_ray_t* ray,
                                              const alea_cell_entry_t* cell,
                                              double t_min,
                                              double t_max,
                                              alea_raycast_result_t* result) {
    if (!cell->surface_indices || cell->surface_index_count == 0) {
        raycast_tree_primitives(sys, ray, cell->root_node_id,
                                t_min, t_max, result);
        return;
    }

    for (size_t i = 0; i < cell->surface_index_count; i++) {
        uint32_t surf_idx = cell->surface_indices[i];
        if (surf_idx >= alea_vec_count(&sys->surfaces)) continue;

        const alea_surface_entry_t* surf = &sys->surfaces.data[surf_idx];
        if (surf->primitive_id >= alea_vec_count(&sys->primitives)) continue;

        const alea_primitive_entry_t* prim =
            &sys->primitives.data[surf->primitive_id];
        alea_primitive_data_t prim_data;
        if (!raycast_primitive_copy_payload(sys, surf->primitive_id,
                                            prim->type, &prim_data)) {
            continue;
        }

        result->surfaces_tested++;
        double t[4];
        int count = ray_intersect_primitive(ray, prim->type, &prim_data, t);
        for (int j = 0; j < count; j++) {
            if (t[j] >= t_min && t[j] <= t_max) {
                alea_ray_hit_t hit;
                hit.t = t[j];
                hit.surface_id = surf->mc_surface_id;
                hit.primitive_id = surf->primitive_id;
                double px, py, pz;
                alea_ray_point_at(ray, t[j], &px, &py, &pz);
                primitive_normal_at(prim->type, &prim_data, px, py, pz,
                                    &hit.nx, &hit.ny, &hit.nz);
                if (add_hit(result, &hit) != 0) {
                    ALEA_LOG_WARN("add_hit failed (out of memory) - raycast results may be incomplete");
                    return;
                }
            }
        }
    }
}

/**
 * Raycast all surfaces in a universe using a translated ray.
 * The ray is in element-local coordinates (origin shifted by -element_center).
 * Since translation preserves t, hits can be added directly to the result.
 */
static void raycast_universe_surfaces(alea_system_t* sys,
                                      const alea_ray_t* local_ray,
                                      int universe_id,
                                      double t_min, double t_max,
                                      alea_raycast_result_t* result) {
    const alea_universe_t* univ = alea_get_universe(sys, universe_id);
    if (!univ) return;

    for (size_t c = 0; c < univ->cell_indices.count; c++) {
        size_t cell_idx = univ->cell_indices.data[c];
        const alea_cell_entry_t* cell = &sys->cells.data[cell_idx];
        raycast_tree_primitives(sys, local_ray, cell->root_node_id,
                                t_min, t_max, result);
    }
}

static bool system_has_lattice_cells(const alea_system_t* sys);
static bool system_has_fill_cells(const alea_system_t* sys);
static int dedup_sorted_hits(alea_raycast_result_t* result);
static void transform_ray_inverse(const alea_matrix_t* mat,
                                  const alea_ray_t* ray,
                                  alea_ray_t* local_ray);
static void raycast_lattice_rect(alea_system_t* sys,
                                 const alea_ray_t* ray,
                                 const alea_cell_entry_t* lat_cell,
                                 double t_min, double t_max,
                                 alea_raycast_result_t* result);
static void raycast_lattice_hex(alea_system_t* sys,
                                const alea_ray_t* ray,
                                const alea_cell_entry_t* lat_cell,
                                double t_min, double t_max,
                                alea_raycast_result_t* result);

static int raycast_root_blas_surfaces(alea_system_t* sys,
                                      const alea_ray_t* ray,
                                      double t_min,
                                      double t_max,
                                      alea_raycast_result_t* result) {
    alea_raycast_result_clear(result);
    result->ray = *ray;

    const alea_hier_spatial_stats_t* stats =
        alea_hier_spatial_index_stats(sys->hier_spatial_index);
    if (!stats || stats->placement_count == 0 ||
        stats->max_universe_cells <= 0) {
        return 0;
    }

    alea_hier_placement_ray_candidate_t* placements =
        malloc(stats->placement_count * sizeof(*placements));
    alea_hier_ray_candidate_t* candidates =
        malloc((size_t)stats->max_universe_cells * sizeof(*candidates));
    if (!placements || !candidates) {
        free(placements);
        free(candidates);
        return -1;
    }

    int placement_count = alea_hier_spatial_query_placements_ray(
        sys,
        ray->ox, ray->oy, ray->oz,
        ray->dx, ray->dy, ray->dz,
        ray->inv_dx, ray->inv_dy, ray->inv_dz,
        t_min, t_max,
        placements, stats->placement_count);
    if (placement_count < 0) {
        free(placements);
        free(candidates);
        return -1;
    }
    result->blas_placement_candidates = (size_t)placement_count;

    for (int p = 0; p < placement_count; p++) {
        alea_hier_placement_ray_candidate_t* placement = &placements[p];
        const alea_universe_t* univ =
            alea_get_universe(sys, placement->universe_id);
        if (!univ || univ->cell_indices.count == 0) continue;

        alea_ray_t local_ray;
        transform_ray_inverse(&placement->transform, ray, &local_ray);

        size_t max_candidates = univ->cell_indices.count;
        if (max_candidates > (size_t)stats->max_universe_cells) {
            max_candidates = (size_t)stats->max_universe_cells;
        }

        double placement_min = fmax(t_min, placement->t_enter - RAY_EPSILON);
        double placement_max = fmin(t_max, placement->t_exit + RAY_EPSILON);
        bool lattice_placement = false;
        if ((size_t)placement->parent_cell_index < alea_vec_count(&sys->cells)) {
            const alea_cell_entry_t* parent_cell =
                &sys->cells.data[placement->parent_cell_index];
            lattice_placement = parent_cell->lat_type != 0 && parent_cell->lat_fill;
        }
        if (placement->placement_index != 0 && !lattice_placement) {
            double samples[3] = {
                placement_min + RAY_EPSILON,
                0.5 * (placement_min + placement_max),
                placement_max - RAY_EPSILON
            };
            int chain_valid = 0;
            for (int s = 0; s < 3; s++) {
                if (samples[s] < placement_min) samples[s] = placement_min;
                if (samples[s] > placement_max) samples[s] = placement_max;
                double sx, sy, sz;
                alea_ray_point_at(ray, samples[s], &sx, &sy, &sz);
                int valid = alea_hier_spatial_check_placement_chain(
                    sys, placement->placement_index, sx, sy, sz);
                if (valid < 0) {
                    free(placements);
                    free(candidates);
                    return -1;
                }
                if (valid) {
                    chain_valid = 1;
                    break;
                }
            }
            if (!chain_valid) {
                result->blas_placements_pruned++;
                continue;
            }
        }

        result->blas_universe_queries++;
        int count = alea_hier_spatial_query_universe_ray(
            sys, placement->universe_id,
            local_ray.ox, local_ray.oy, local_ray.oz,
            local_ray.dx, local_ray.dy, local_ray.dz,
            local_ray.inv_dx, local_ray.inv_dy, local_ray.inv_dz,
            placement_min, placement_max,
            candidates, max_candidates);
        if (count < 0) {
            free(placements);
            free(candidates);
            return -1;
        }
        result->blas_cell_candidates += (size_t)count;

        for (int i = 0; i < count; i++) {
            uint32_t cell_index = candidates[i].cell_index;
            if ((size_t)cell_index >= alea_vec_count(&sys->cells)) continue;
            const alea_cell_entry_t* cell = &sys->cells.data[cell_index];
            if (cell->root_node_id == ALEA_NODE_ID_INVALID) continue;
            result->blas_cells_tested++;

            double local_min = fmax(placement_min,
                                    candidates[i].t_enter - RAY_EPSILON);
            double local_max = fmin(placement_max,
                                    candidates[i].t_exit + RAY_EPSILON);
            raycast_cell_indexed_surface_hits(sys, &local_ray, cell,
                                              local_min, local_max, result);
            if (cell->lat_type == 1 && cell->lat_fill) {
                raycast_lattice_rect(sys, &local_ray, cell,
                                     local_min, local_max, result);
            } else if (cell->lat_type == 2 && cell->lat_fill) {
                raycast_lattice_hex(sys, &local_ray, cell,
                                    local_min, local_max, result);
            }
        }
    }

    free(placements);
    free(candidates);
    result->blas_hits_before_dedup = result->hits.count;
    dedup_sorted_hits(result);
    return 0;
}

static void transform_ray_inverse(const alea_matrix_t* mat,
                                  const alea_ray_t* ray,
                                  alea_ray_t* local_ray) {
    double ox = ray->ox;
    double oy = ray->oy;
    double oz = ray->oz;
    alea_matrix_transform_point_inverse(mat, &ox, &oy, &oz);

    double dx = mat->inv[0] * ray->dx + mat->inv[1] * ray->dy + mat->inv[2] * ray->dz;
    double dy = mat->inv[4] * ray->dx + mat->inv[5] * ray->dy + mat->inv[6] * ray->dz;
    double dz = mat->inv[8] * ray->dx + mat->inv[9] * ray->dy + mat->inv[10] * ray->dz;

    if (alea_ray_init(local_ray, ox, oy, oz, dx, dy, dz) != 0) {
        alea_ray_init_normalized(local_ray, ox, oy, oz, ray->dx, ray->dy, ray->dz);
    }
}

static void raycast_fill_universe_hits_recursive(alea_system_t* sys,
                                                 const alea_ray_t* global_ray,
                                                 int universe_id,
                                                 const alea_matrix_t* accumulated,
                                                 double t_min, double t_max,
                                                 int depth,
                                                 alea_raycast_result_t* result) {
    if (!sys || !global_ray || !accumulated || !result)
        return;
    if (depth >= MAX_FILL_RAYCAST_DEPTH)
        return;

    const alea_universe_t* univ = alea_get_universe(sys, universe_id);
    if (!univ)
        return;

    alea_ray_t parent_local_ray;
    transform_ray_inverse(accumulated, global_ray, &parent_local_ray);

    for (size_t i = 0; i < univ->cell_indices.count; i++) {
        size_t cell_idx = univ->cell_indices.data[i];
        if (cell_idx >= alea_vec_count(&sys->cells))
            continue;

        const alea_cell_entry_t* cell = &sys->cells.data[cell_idx];
        if (cell->fill_universe <= 0)
            continue;

        if (cell->root_node_id != ALEA_NODE_ID_INVALID &&
            cell->root_node_id < alea_vec_count(&sys->nodes)) {
            const alea_bbox_t bbox_v = alea_node_bbox_get(&sys->nodes.data[cell->root_node_id].bbox);
            if (!ray_bbox_slab(&parent_local_ray, &bbox_v, t_min, t_max))
                continue;
        }

        alea_matrix_t fill_transform;
        if (cell->fill_transform > 0) {
            const alea_transform_t* tr = alea_get_transform(sys, cell->fill_transform);
            if (!tr ||
                !alea_matrix_from_mcnp(&fill_transform, tr->cosines,
                                       tr->value_count, false)) {
                continue;
            }
        } else {
            alea_matrix_identity(&fill_transform);
        }

        alea_matrix_t fill_accumulated;
        alea_matrix_multiply(&fill_accumulated, accumulated, &fill_transform);
        if (!fill_accumulated.has_inverse && !alea_matrix_invert(&fill_accumulated))
            continue;

        alea_ray_t fill_local_ray;
        transform_ray_inverse(&fill_accumulated, global_ray, &fill_local_ray);

        raycast_universe_surfaces(sys, &fill_local_ray, cell->fill_universe,
                                  t_min, t_max, result);
        raycast_fill_universe_hits_recursive(sys, global_ray, cell->fill_universe,
                                             &fill_accumulated, t_min, t_max,
                                             depth + 1, result);
    }
}

static void raycast_add_fill_hits(alea_system_t* sys,
                                  const alea_ray_t* ray,
                                  double t_min, double t_max,
                                  alea_raycast_result_t* result) {
    alea_matrix_t identity;
    alea_matrix_identity(&identity);
    raycast_fill_universe_hits_recursive(sys, ray, 0, &identity,
                                         t_min, t_max, 0, result);
}

/* ray_bbox_slab_enter_exit in ray_bbox.h replaces ray_bbox_slab_enter_exit */

/**
 * DDA raycast through a rectangular lattice.
 * Steps through elements along the ray, raycasting base universe surfaces
 * in each element's local coordinate system.
 */
static void raycast_lattice_rect(alea_system_t* sys,
                                 const alea_ray_t* ray,
                                 const alea_cell_entry_t* lat_cell,
                                 double t_min, double t_max,
                                 alea_raycast_result_t* result) {
    int ni = lat_cell->lat_fill_dims[1] - lat_cell->lat_fill_dims[0] + 1;
    int nj = lat_cell->lat_fill_dims[3] - lat_cell->lat_fill_dims[2] + 1;
    int nk = lat_cell->lat_fill_dims[5] - lat_cell->lat_fill_dims[4] + 1;

    double px = lat_cell->lat_pitch[0];
    double py = lat_cell->lat_pitch[1];
    double pz = lat_cell->lat_pitch[2];
    const double* ll = lat_cell->lat_lower_left;

    /* Lattice bounding box */
    alea_bbox_t lat_bbox = {
        .min_x = ll[0], .max_x = ll[0] + ni * px,
        .min_y = ll[1], .max_y = ll[1] + nj * py,
        .min_z = ll[2], .max_z = ll[2] + nk * pz,
    };

    double t_enter, t_exit;
    if (!ray_bbox_slab_enter_exit(ray, &lat_bbox, t_min, t_max, &t_enter, &t_exit))
        return;

    /* Starting point (nudge slightly inside) */
    double start_t = t_enter + RAY_EPSILON;
    double sx = ray->ox + start_t * ray->dx;
    double sy = ray->oy + start_t * ray->dy;
    double sz = ray->oz + start_t * ray->dz;

    int i = (int)floor((sx - ll[0]) / px);
    int j = (nj == 1) ? 0 : (int)floor((sy - ll[1]) / py);
    int k = (nk == 1) ? 0 : (int)floor((sz - ll[2]) / pz);

    /* Clamp to valid range (handles edge cases) */
    if (i < 0) i = 0;
    if (i >= ni) i = ni - 1;
    if (j < 0) j = 0;
    if (j >= nj) j = nj - 1;
    if (k < 0) k = 0;
    if (k >= nk) k = nk - 1;

    /* DDA step direction and delta-t per cell crossing */
    int step_i = (ray->dx > 0) ? 1 : -1;
    int step_j = (ray->dy > 0) ? 1 : -1;
    int step_k = (ray->dz > 0) ? 1 : -1;

    /* Distance to next grid boundary in each dimension */
    double t_next_i, t_next_j, t_next_k;
    if (fabs(ray->dx) > RAY_EPSILON) {
        double boundary = ll[0] + ((ray->dx > 0) ? (i + 1) : i) * px;
        t_next_i = (boundary - ray->ox) / ray->dx;
    } else {
        t_next_i = DBL_MAX;
    }
    if (nj > 1 && fabs(ray->dy) > RAY_EPSILON) {
        double boundary = ll[1] + ((ray->dy > 0) ? (j + 1) : j) * py;
        t_next_j = (boundary - ray->oy) / ray->dy;
    } else {
        t_next_j = DBL_MAX;
    }
    if (nk > 1 && fabs(ray->dz) > RAY_EPSILON) {
        double boundary = ll[2] + ((ray->dz > 0) ? (k + 1) : k) * pz;
        t_next_k = (boundary - ray->oz) / ray->dz;
    } else {
        t_next_k = DBL_MAX;
    }

    /* Walk through lattice elements */
    double t_cur = t_enter;
    int max_steps = 2 * (ni + nj + nk) + 4;  /* Safety limit: diagonal ray crosses ~3N boundaries */

    for (int step = 0; step < max_steps; step++) {
        if (i < 0 || i >= ni || j < 0 || j >= nj || k < 0 || k >= nk)
            break;

        /* Emit synthetic hit at internal element boundaries so the
         * segment builder knows to re-query the cell at transitions
         * between lattice elements.  surface_id=0 causes the segment
         * builder to skip the (universe-unaware) neighbor lookup and
         * fall through to a full point query. */
        if (step > 0 && t_cur > t_enter + RAY_EPSILON) {
            alea_ray_hit_t bnd = { .t = t_cur, .surface_id = 0,
                                   .primitive_id = ALEA_PRIMITIVE_ID_INVALID };
            if (add_hit(result, &bnd) != 0) {
                ALEA_LOG_WARN("add_hit failed (out of memory) - lattice raycast incomplete");
                return;
            }
        }

        /* t range for this element */
        double t_min_next = t_next_i;
        if (t_next_j < t_min_next) t_min_next = t_next_j;
        if (t_next_k < t_min_next) t_min_next = t_next_k;
        if (t_min_next > t_exit) t_min_next = t_exit;

        /* Get universe for this element */
        size_t idx = (size_t)(i * nj * nk + j * nk + k);
        if (idx < lat_cell->lat_fill_count) {
            int univ_id = lat_cell->lat_fill[idx];

            /* Translate ray to element-local coordinates */
            double cx = ll[0] + (i + 0.5) * px;
            double cy = ll[1] + (j + 0.5) * py;
            double cz = ll[2] + (k + 0.5) * pz;

            alea_ray_t local_ray = *ray;
            local_ray.ox -= cx;
            local_ray.oy -= cy;
            local_ray.oz -= cz;

            raycast_universe_surfaces(sys, &local_ray, univ_id,
                                      t_cur, t_min_next, result);
        }

        /* Step to next element */
        t_cur = t_min_next;
        if (t_cur >= t_exit) break;

        /* Step indices and recalculate boundary positions from scratch
         * to avoid accumulating floating-point errors. */
        if (t_next_i <= t_next_j && t_next_i <= t_next_k) {
            i += step_i;
            if (fabs(ray->dx) > RAY_EPSILON) {
                double boundary = ll[0] + ((ray->dx > 0) ? (i + 1) : i) * px;
                t_next_i = (boundary - ray->ox) / ray->dx;
            }
        } else if (t_next_j <= t_next_k) {
            j += step_j;
            if (fabs(ray->dy) > RAY_EPSILON) {
                double boundary = ll[1] + ((ray->dy > 0) ? (j + 1) : j) * py;
                t_next_j = (boundary - ray->oy) / ray->dy;
            }
        } else {
            k += step_k;
            if (fabs(ray->dz) > RAY_EPSILON) {
                double boundary = ll[2] + ((ray->dz > 0) ? (k + 1) : k) * pz;
                t_next_k = (boundary - ray->oz) / ray->dz;
            }
        }
    }
}

/**
 * Hex DDA: walk a ray through a hex lattice, raycasting fill‐universe
 * surfaces in each element.  Uses 3-axis DDA on the oblique hex
 * coordinate system (fi, fj, fi+fj) to find element transitions.
 */
static void raycast_lattice_hex(alea_system_t* sys,
                                const alea_ray_t* ray,
                                const alea_cell_entry_t* lat_cell,
                                double t_min, double t_max,
                                alea_raycast_result_t* result) {
    double p = lat_cell->lat_pitch[0];
    if (p <= 0.0) return;

    /* M_M_SQRT3 from util/math.h */

    int imin = lat_cell->lat_fill_dims[0], imax = lat_cell->lat_fill_dims[1];
    int jmin = lat_cell->lat_fill_dims[2], jmax = lat_cell->lat_fill_dims[3];
    int nk   = lat_cell->lat_fill_dims[5] - lat_cell->lat_fill_dims[4] + 1;
    int ni   = imax - imin + 1;
    int nj   = jmax - jmin + 1;

    /* AABB of all hex element centers, expanded by circumradius */
    double bb_min_x =  DBL_MAX, bb_max_x = -DBL_MAX;
    double bb_min_y =  DBL_MAX, bb_max_y = -DBL_MAX;
    double R = p / M_SQRT3;  /* circumradius */

    for (int ci = imin; ci <= imax; ci++) {
        for (int cj = jmin; cj <= jmax; cj++) {
            double cx = ci * p + cj * p * 0.5;
            double cy = cj * p * M_SQRT3 * 0.5;
            if (cx - R < bb_min_x) bb_min_x = cx - R;
            if (cx + R > bb_max_x) bb_max_x = cx + R;
            if (cy - p * 0.5 < bb_min_y) bb_min_y = cy - p * 0.5;
            if (cy + p * 0.5 > bb_max_y) bb_max_y = cy + p * 0.5;
        }
    }

    alea_bbox_t lat_bbox = {
        .min_x = bb_min_x, .max_x = bb_max_x,
        .min_y = bb_min_y, .max_y = bb_max_y,
        .min_z = (nk > 1) ? lat_cell->lat_lower_left[2] : -1e30,
        .max_z = (nk > 1) ? lat_cell->lat_lower_left[2]
                           + nk * lat_cell->lat_pitch[2] : 1e30,
    };

    double t_enter, t_exit;
    if (!ray_bbox_slab_enter_exit(ray, &lat_bbox, t_min, t_max, &t_enter, &t_exit))
        return;

    /* ---- oblique coordinate rates along the ray ---- */
    double inv_p = 1.0 / p;
    double inv_ps = inv_p / M_SQRT3;          /* 1 / (p * sqrt3) */
    double dq = ray->dx * inv_p - ray->dy * inv_ps;          /* dfi/dt */
    double dr = 2.0 * ray->dy * inv_ps;                      /* dfj/dt */
    double ds = -dq - dr;                                     /* d(-fi-fj)/dt */

    /* Starting oblique coords at t_enter */
    double sx = ray->ox + t_enter * ray->dx;
    double sy = ray->oy + t_enter * ray->dy;
    double q0 = sx * inv_p - sy * inv_ps;
    double r0 = 2.0 * sy * inv_ps;
    double s0 = -q0 - r0;

    /* Helper: compute first half-integer boundary crossing after val in
     * the direction of step.  Boundaries are at n + 0.5. */
    #define NEXT_HALF(val, rate, t_out) do {                       \
        if (fabs(rate) > RAY_EPSILON) {                            \
            double bnd;                                            \
            if ((rate) > 0) {                                      \
                bnd = floor((val) + 0.5) + 0.5;                   \
                if (bnd <= (val) + RAY_EPSILON) bnd += 1.0;        \
            } else {                                               \
                bnd = ceil((val) - 0.5) - 0.5;                    \
                if (bnd >= (val) - RAY_EPSILON) bnd -= 1.0;        \
            }                                                      \
            (t_out) = t_enter + (bnd - (val)) / (rate);            \
        } else {                                                   \
            (t_out) = DBL_MAX;                                     \
        }                                                          \
    } while(0)

    double t_next_q;  NEXT_HALF(q0, dq, t_next_q);
    double t_next_r;  NEXT_HALF(r0, dr, t_next_r);
    double t_next_s;  NEXT_HALF(s0, ds, t_next_s);
    #undef NEXT_HALF

    /* ---- Walk through hex elements ---- */
    double t_cur = t_enter;
    double prev_cx = DBL_MAX;   /* track previous element center */
    double prev_cy = DBL_MAX;
    int max_steps = 3 * (ni + nj + nk) + 10;

    for (int step = 0; step < max_steps; step++) {
        if (t_cur >= t_exit - RAY_EPSILON) break;

        /* Nearest boundary crossing among the 3 axes */
        double t_next = t_next_q;
        if (t_next_r < t_next) t_next = t_next_r;
        if (t_next_s < t_next) t_next = t_next_s;
        if (t_next > t_exit) t_next = t_exit;

        /* Identify hex element at midpoint of this interval */
        double t_mid = 0.5 * (t_cur + t_next);
        double mx = ray->ox + t_mid * ray->dx;
        double my = ray->oy + t_mid * ray->dy;
        double mz = ray->oz + t_mid * ray->dz;

        double ox, oy, oz;
        int univ = lattice_hex_lookup(lat_cell, mx, my, mz, &ox, &oy, &oz);

        if (univ >= 0) {
            /* Emit boundary hit when we enter a new hex element */
            int new_elem = (fabs(ox - prev_cx) > RAY_EPSILON ||
                            fabs(oy - prev_cy) > RAY_EPSILON);
            if (new_elem && prev_cx < DBL_MAX && t_cur > t_enter + RAY_EPSILON) {
                alea_ray_hit_t bnd = { .t = t_cur, .surface_id = 0,
                                       .primitive_id = ALEA_PRIMITIVE_ID_INVALID };
                if (add_hit(result, &bnd) != 0) {
                    ALEA_LOG_WARN("add_hit failed (out of memory) - hex lattice raycast incomplete");
                    return;
                }
            }
            if (new_elem) { prev_cx = ox; prev_cy = oy; }

            /* Raycast fill universe in element-local coordinates */
            alea_ray_t local_ray = *ray;
            local_ray.ox -= ox;
            local_ray.oy -= oy;
            local_ray.oz -= oz;
            raycast_universe_surfaces(sys, &local_ray, univ,
                                      t_cur, t_next, result);
        }

        /* Advance DDA — recompute next boundaries from current position
         * to avoid accumulating floating-point drift from repeated += dt */
        t_cur = t_next;
        if (fabs(t_next - t_next_q) < RAY_EPSILON) {
            double cur_q = (ray->ox + t_cur * ray->dx) * inv_p
                         - (ray->oy + t_cur * ray->dy) * inv_ps;
            double bnd = (dq > 0) ? floor(cur_q + 0.5) + 0.5
                                  : ceil(cur_q - 0.5) - 0.5;
            if ((dq > 0 && bnd <= cur_q + RAY_EPSILON) ||
                (dq < 0 && bnd >= cur_q - RAY_EPSILON))
                bnd += (dq > 0) ? 1.0 : -1.0;
            t_next_q = t_cur + (bnd - cur_q) / dq;
        }
        if (fabs(t_next - t_next_r) < RAY_EPSILON) {
            double cur_r = 2.0 * (ray->oy + t_cur * ray->dy) * inv_ps;
            double bnd = (dr > 0) ? floor(cur_r + 0.5) + 0.5
                                  : ceil(cur_r - 0.5) - 0.5;
            if ((dr > 0 && bnd <= cur_r + RAY_EPSILON) ||
                (dr < 0 && bnd >= cur_r - RAY_EPSILON))
                bnd += (dr > 0) ? 1.0 : -1.0;
            t_next_r = t_cur + (bnd - cur_r) / dr;
        }
        if (fabs(t_next - t_next_s) < RAY_EPSILON) {
            double cur_q2 = (ray->ox + t_cur * ray->dx) * inv_p
                          - (ray->oy + t_cur * ray->dy) * inv_ps;
            double cur_r2 = 2.0 * (ray->oy + t_cur * ray->dy) * inv_ps;
            double cur_s = -cur_q2 - cur_r2;
            double bnd = (ds > 0) ? floor(cur_s + 0.5) + 0.5
                                  : ceil(cur_s - 0.5) - 0.5;
            if ((ds > 0 && bnd <= cur_s + RAY_EPSILON) ||
                (ds < 0 && bnd >= cur_s - RAY_EPSILON))
                bnd += (ds > 0) ? 1.0 : -1.0;
            t_next_s = t_cur + (bnd - cur_s) / ds;
        }
    }
}

/**
 * Add lattice surface hits for all lattice cells the ray may cross.
 * Called after alea_raycast_surfaces() and before sorting/dedup.
 */
static void raycast_add_lattice_hits(alea_system_t* sys,
                                     const alea_ray_t* ray,
                                     double t_min, double t_max,
                                     alea_raycast_result_t* result) {
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        const alea_cell_entry_t* cell = &sys->cells.data[i];
        if (cell->lat_type == 0 || !cell->lat_fill) continue;

        if (cell->lat_type == 1) {
            raycast_lattice_rect(sys, ray, cell, t_min, t_max, result);
        } else if (cell->lat_type == 2) {
            raycast_lattice_hex(sys, ray, cell, t_min, t_max, result);
        }
    }
}

static bool system_has_lattice_cells(const alea_system_t* sys) {
    if (!sys) return false;
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        const alea_cell_entry_t* cell = &sys->cells.data[i];
        if (cell->lat_type != 0 && cell->lat_fill && cell->lat_fill_count > 0)
            return true;
    }
    return false;
}

static bool system_has_fill_cells(const alea_system_t* sys) {
    if (!sys) return false;
    for (size_t i = 0; i < alea_vec_count(&sys->cells); i++) {
        if (sys->cells.data[i].fill_universe > 0)
            return true;
    }
    return false;
}

static int dedup_sorted_hits(alea_raycast_result_t* result) {
    if (!result) return -1;
    sort_hits(result->hits.data, result->hits.count);

    if (result->hits.count > 1) {
        size_t write = 1;
        for (size_t read = 1; read < result->hits.count; read++) {
            int same_t = fabs(result->hits.data[read].t -
                              result->hits.data[write - 1].t) <= DEDUP_EPSILON;
            int same_surf = result->hits.data[read].surface_id ==
                            result->hits.data[write - 1].surface_id;
            if (!(same_t && same_surf)) {
                result->hits.data[write++] = result->hits.data[read];
            }
        }
        result->hits.count = write;
    }

    return 0;
}

static int raycast_global_pipeline(alea_system_t* sys,
                                   const alea_ray_t* ray,
                                   double t_min, double t_max,
                                   bool include_lattice_hits,
                                   bool use_hier_lookup,
                                   alea_raycast_result_t* result) {
    if (use_hier_lookup) {
        raycast_surfaces_impl(sys, ray, t_min, t_max, result);
    } else {
        int rc = alea_raycast_surfaces(sys, ray, t_min, t_max, result);
        if (rc != 0) return rc;
    }

    if (system_has_fill_cells(sys))
        raycast_add_fill_hits(sys, ray, t_min, t_max, result);

    if (include_lattice_hits)
        raycast_add_lattice_hits(sys, ray, t_min, t_max, result);

    dedup_sorted_hits(result);
    return raycast_to_segments_impl(sys, t_max, result, use_hier_lookup);
}

/* ============================================================================
 * CONVENIENCE FUNCTIONS
 * ============================================================================ */

int alea_raycast(alea_system_t* sys,
                double ox, double oy, double oz,
                double dx, double dy, double dz,
                double t_max,
                alea_raycast_result_t* result) {
    if (!result) return -1;

    /* alea_raycast() is the hit-producing surface+segment pipeline (uses the
     * surface BVH). For segment-only fast tracing use alea_raycast_hier(); for
     * hit-producing hierarchical tracing use alea_raycast_hier_with_hits(). */

    /* Free any prior allocations then reinitialize (safe for reuse) */
    alea_raycast_result_free(result);

    alea_ray_t ray;
    if (alea_ray_init(&ray, ox, oy, oz, dx, dy, dz) != 0) {
        return -1;  /* Zero-length direction */
    }

    /* t_max=0 means infinite */
    double effective_t_max = (t_max <= 0) ? DBL_MAX : t_max;
    return raycast_global_pipeline(sys, &ray, 0, effective_t_max,
                                   system_has_lattice_cells(sys), false,
                                   result);
}

int alea_raycast_hier(alea_system_t* sys,
                      double ox, double oy, double oz,
                      double dx, double dy, double dz,
                      double t_max,
                      alea_raycast_result_t* result) {
    if (!sys || !result) return -1;

    if (alea_system_prepare_query_caches(sys,
            ALEA_CACHE_HIER_SPATIAL | ALEA_CACHE_CELL_SURFACES) != 0)
        return -1;

    alea_raycast_result_free(result);

    alea_ray_t ray;
    if (alea_ray_init(&ray, ox, oy, oz, dx, dy, dz) != 0) {
        return -1;
    }

    double effective_t_max = (t_max <= 0) ? DBL_MAX : t_max;
    result->ray = ray;
    return raycast_cell_aware_impl(sys, &ray, effective_t_max, true, false,
                                   result);
}

int alea_raycast_hier_blas_experimental(alea_system_t* sys,
                                        double ox, double oy, double oz,
                                        double dx, double dy, double dz,
                                        double t_max,
                                        alea_raycast_result_t* result) {
    if (!sys || !result) return -1;
    if (alea_raycast_ensure_hier_caches(sys) != 0) return -1;

    alea_raycast_result_free(result);

    alea_ray_t ray;
    if (alea_ray_init(&ray, ox, oy, oz, dx, dy, dz) != 0) {
        return -1;
    }
    result->ray = ray;

    double effective_t_max = (t_max <= 0) ? DBL_MAX : t_max;
    int rc = raycast_root_blas_surfaces(sys, &ray, 0.0, effective_t_max,
                                        result);
    if (rc == 1) {
        return alea_raycast_hier(sys, ox, oy, oz, dx, dy, dz, t_max, result);
    }
    if (rc != 0) return rc;

    return raycast_to_segments_impl(sys, effective_t_max, result, true);
}

int alea_ray_first_cell(alea_system_t* sys,
                        double ox, double oy, double oz,
                       double dx, double dy, double dz,
                       double t_max,
                       double* out_t) {
    alea_raycast_result_t result;
    alea_raycast_result_init(&result);

    int rc = alea_raycast(sys, ox, oy, oz, dx, dy, dz, t_max, &result);
    if (rc != 0) {
        alea_raycast_result_free(&result);
        return -1;
    }

    int first_cell = -1;
    if (result.segments.count > 0) {
        /* Find first non-void segment */
        for (size_t i = 0; i < result.segments.count; i++) {
            if (result.segments.data[i].cell_id >= 0) {
                first_cell = result.segments.data[i].cell_id;
                if (out_t) *out_t = result.segments.data[i].t_enter;
                break;
            }
        }
    }

    alea_raycast_result_free(&result);
    return first_cell;
}

int alea_ray_is_occluded(alea_system_t* sys,
                        double ox, double oy, double oz,
                        double dx, double dy, double dz,
                        double t_max) {
    if (!sys) return 0;

    if (raycast_prefers_hier_mode(sys)) {
        alea_raycast_result_t result;
        alea_raycast_result_init(&result);
        int rc = alea_raycast_hier(sys, ox, oy, oz, dx, dy, dz, t_max, &result);
        if (rc != 0) {
            alea_raycast_result_free(&result);
            return 0;
        }
        int occluded = 0;
        for (size_t i = 0; i < result.segments.count; i++) {
            if (result.segments.data[i].cell_id >= 0 &&
                result.segments.data[i].material_id != 0) {
                occluded = 1;
                break;
            }
        }
        alea_raycast_result_free(&result);
        return occluded;
    }

    /* Reuse a thread-local result to avoid malloc per shadow ray */
    static ALEA_THREAD_LOCAL alea_raycast_result_t tls_result;
    static ALEA_THREAD_LOCAL int tls_init = 0;
    if (!tls_init) {
        alea_raycast_result_init(&tls_result);
        alea_raycast_result_reserve(&tls_result, 64, 32);
        tls_init = 1;
    }
    alea_raycast_result_clear(&tls_result);

    alea_ray_t ray;
    alea_ray_init(&ray, ox, oy, oz, dx, dy, dz);

    double effective_t_max = (t_max <= 0) ? DBL_MAX : t_max;

    /* Get surface intersections */
    if (alea_raycast_surfaces(sys, &ray, 0, effective_t_max, &tls_result) != 0)
        return 0;

    /* Sort hits by distance */
    sort_hits(tls_result.hits.data, tls_result.hits.count);

    /* Walk intervals and early-out on first non-void cell */
    double t_prev = 0;
    int prev_cell_idx = -1;

    for (size_t i = 0; i <= tls_result.hits.count; i++) {
        double t_curr = (i < tls_result.hits.count) ? tls_result.hits.data[i].t : effective_t_max;

        if (t_curr > t_prev + RAY_EPSILON) {
            int cell_id, cell_idx, material_id;
            double density;
            int crossed_surface = (i > 0) ? tls_result.hits.data[i - 1].surface_id : -1;
            find_cell_after_crossing(sys, &ray, t_prev, t_curr,
                                     prev_cell_idx, crossed_surface,
                                     false,
                                     &cell_id, &cell_idx, &material_id, &density);

            prev_cell_idx = cell_idx;

            /* Early out: found a non-void cell with non-void material */
            if (cell_id >= 0 && material_id != 0) {
                return 1;
            }
        }

        t_prev = t_curr;
    }

    return 0;
}

double alea_raycast_path_length(const alea_raycast_result_t* result,
                               int material_id) {
    if (!result) return 0;

    double total = 0;
    for (size_t i = 0; i < result->segments.count; i++) {
        const alea_ray_segment_t* seg = &result->segments.data[i];

        /* Skip infinite segments */
        if (seg->t_exit >= DBL_MAX - 1) continue;

        int match = (material_id < 0) ||  /* -1 = all materials */
                    (seg->material_id == material_id);

        if (match) {
            total += seg->t_exit - seg->t_enter;
        }
    }

    return total;
}

/* ============================================================================
 * CELL-AWARE RAYCAST
 * ============================================================================ */

/**
 * One boundary crossing produced by a single stepper iteration.
 *
 * The segment-only path consumes only `t` / `surface_id`; the hit-producing
 * path (alea_raycast_hier_with_hits) also needs to materialize a world-space
 * surface normal. Rather than carry the primitive payload and local ray by
 * value through the hot loop, we keep the winning primitive id plus the frame
 * transform (primitive-local -> world). At emit time the payload is fetched
 * from `primitive_id` and the local hit point is recovered from the world hit
 * point via the transform, so the segment-only path pays nothing for this.
 */
typedef struct {
    double t;                /* distance of the boundary along the world ray */
    int surface_id;          /* MCNP id: -1 none, 0 synthetic lattice, >0 physical */
    uint32_t primitive_id;   /* winning primitive, ALEA_PRIMITIVE_ID_INVALID if none */
    alea_matrix_t transform; /* primitive-local -> world frame of the winning surface */
    bool has_physical_surface;
    bool is_synthetic_lattice_boundary;
} alea_raycast_boundary_event_t;

/* Map a primitive-local normal back to world space.
 *
 * Normals transform by the inverse-transpose of the local->world linear part.
 * `mat->inv` is the world->local linear map, so its transpose is exactly that
 * inverse-transpose. For an identity transform this is a no-op, which keeps
 * non-lattice/non-fill normals bit-identical to the flat path. */
static void boundary_event_world_normal(const alea_matrix_t* mat,
                                        double nlx, double nly, double nlz,
                                        double* nx, double* ny, double* nz) {
    double wx = mat->inv[0] * nlx + mat->inv[4] * nly + mat->inv[8] * nlz;
    double wy = mat->inv[1] * nlx + mat->inv[5] * nly + mat->inv[9] * nlz;
    double wz = mat->inv[2] * nlx + mat->inv[6] * nly + mat->inv[10] * nlz;
    double len = sqrt(wx * wx + wy * wy + wz * wz);
    if (len > RAY_EPSILON) {
        *nx = wx / len;
        *ny = wy / len;
        *nz = wz / len;
    } else {
        *nx = nlx;
        *ny = nly;
        *nz = nlz;
    }
}

/* Emit one hit for the winning boundary event of a stepper iteration.
 * Returns the index of the pushed hit in result->hits, or -1 if no physical
 * surface hit was produced. */
static int boundary_event_emit_hit(alea_system_t* sys,
                                    const alea_ray_t* ray,
                                    const alea_raycast_boundary_event_t* ev,
                                    alea_raycast_result_t* result) {
    if (!ev->has_physical_surface ||
        ev->primitive_id == ALEA_PRIMITIVE_ID_INVALID ||
        ev->primitive_id >= alea_vec_count(&sys->primitives)) {
        return -1;
    }

    const alea_primitive_entry_t* prim = &sys->primitives.data[ev->primitive_id];
    alea_primitive_data_t prim_data;
    if (!raycast_primitive_copy_payload(sys, ev->primitive_id,
                                        prim->type, &prim_data)) {
        return -1;
    }

    double wx, wy, wz;
    alea_ray_point_at(ray, ev->t, &wx, &wy, &wz);

    /* Recover the primitive-local hit point and evaluate the local normal. */
    double lx = wx, ly = wy, lz = wz;
    alea_matrix_transform_point_inverse(&ev->transform, &lx, &ly, &lz);
    double nlx, nly, nlz;
    primitive_normal_at(prim->type, &prim_data, lx, ly, lz, &nlx, &nly, &nlz);

    alea_ray_hit_t hit;
    hit.t = ev->t;
    hit.surface_id = ev->surface_id;
    hit.primitive_id = ev->primitive_id;
    boundary_event_world_normal(&ev->transform, nlx, nly, nlz,
                                &hit.nx, &hit.ny, &hit.nz);

    int hit_index = (int)alea_vec_count(&result->hits);
    if (add_hit(result, &hit) != 0) {
        ALEA_LOG_WARN("add_hit failed (out of memory) - hier hits may be incomplete");
        return -1;
    }
    return hit_index;
}

/**
 * Test ray against a specific cell's surfaces only.
 * Returns closest hit distance, or DBL_MAX if no hit.
 *
 * When `out_primitive_id` is non-NULL it receives the canonical primitive id
 * of the closest surface (ALEA_PRIMITIVE_ID_INVALID if none), so the caller
 * can build a boundary event without re-scanning the cell.
 */
static double raycast_cell_surfaces(alea_system_t* sys,
                                    const alea_ray_t* ray,
                                    const alea_cell_entry_t* cell,
                                    double t_min, double t_max,
                                    alea_raycast_result_t* result,
                                    int* out_surface_id,
                                    uint32_t* out_primitive_id) {
    double closest_t = DBL_MAX;
    *out_surface_id = -1;
    if (out_primitive_id) *out_primitive_id = ALEA_PRIMITIVE_ID_INVALID;

    /* Safety check: surface index must be built */
    if (!cell->surface_indices || cell->surface_index_count == 0) {
        return DBL_MAX;  /* No surfaces indexed for this cell */
    }

    /* Test each surface belonging to this cell */
    for (size_t i = 0; i < cell->surface_index_count; i++) {
        uint32_t surf_idx = cell->surface_indices[i];
        if (surf_idx >= alea_vec_count(&sys->surfaces)) continue;

        const alea_surface_entry_t* surf = &sys->surfaces.data[surf_idx];
        const alea_primitive_entry_t* prim = &sys->primitives.data[surf->primitive_id];
        alea_primitive_data_t prim_data;
        if (!raycast_primitive_copy_payload(sys, surf->primitive_id,
                                            prim->type, &prim_data)) continue;
        if (result) {
            result->surfaces_tested++;
            if ((unsigned)prim->type < ALEA_RAYCAST_PRIM_TYPE_BINS)
                result->prim_type_tests[prim->type]++;
        }

        double t[4];
        int count = ray_intersect_primitive(ray, prim->type, &prim_data, t);

        for (int j = 0; j < count; j++) {
            if (t[j] > t_min && t[j] < t_max && t[j] < closest_t) {
                closest_t = t[j];
                *out_surface_id = surf->mc_surface_id;
                if (out_primitive_id) *out_primitive_id = surf->primitive_id;
            }
        }
    }

    return closest_t;
}

static double raycast_hier_path_ancestor_surfaces(alea_system_t* sys,
                                                  const alea_ray_t* ray,
                                                  const alea_hier_ray_path_t* path,
                                                  uint32_t terminal_cell_index,
                                                  int already_tested_lattice_cell,
                                                  double t_min,
                                                  double t_max,
                                                  alea_raycast_result_t* result,
                                                  int* out_surface_id,
                                                  uint32_t* out_primitive_id,
                                                  alea_matrix_t* out_transform) {
    double closest_t = DBL_MAX;
    *out_surface_id = -1;
    if (out_primitive_id) *out_primitive_id = ALEA_PRIMITIVE_ID_INVALID;

    if (!path || path->count <= 1) return closest_t;

    for (int i = 0; i < path->count - 1; i++) {
        const alea_hier_ray_path_entry_t* entry = &path->entries[i];
        if (entry->cell_index == terminal_cell_index) continue;
        if ((int)entry->cell_index == already_tested_lattice_cell) continue;
        if ((size_t)entry->cell_index >= alea_vec_count(&sys->cells)) continue;

        const alea_cell_entry_t* cell = &sys->cells.data[entry->cell_index];
        if (!cell->surface_indices || cell->surface_index_count == 0) continue;

        alea_ray_t local_ray;
        transform_ray_inverse(&entry->transform, ray, &local_ray);

        int surface_id = -1;
        uint32_t prim_id = ALEA_PRIMITIVE_ID_INVALID;
        double t = raycast_cell_surfaces(sys, &local_ray, cell, t_min, t_max,
                                         result, &surface_id, &prim_id);
        if (t < closest_t) {
            closest_t = t;
            *out_surface_id = surface_id;
            if (out_primitive_id) *out_primitive_id = prim_id;
            if (out_transform) *out_transform = entry->transform;
        }
    }

    return closest_t;
}

static double lattice_rect_next_boundary(const alea_ray_t* ray,
                                         const alea_cell_entry_t* lat_cell,
                                         double t_min,
                                         double t_max) {
    int ni = lat_cell->lat_fill_dims[1] - lat_cell->lat_fill_dims[0] + 1;
    int nj = lat_cell->lat_fill_dims[3] - lat_cell->lat_fill_dims[2] + 1;
    int nk = lat_cell->lat_fill_dims[5] - lat_cell->lat_fill_dims[4] + 1;
    double px = lat_cell->lat_pitch[0];
    double py = lat_cell->lat_pitch[1];
    double pz = lat_cell->lat_pitch[2];
    const double* ll = lat_cell->lat_lower_left;

    if (px <= 0.0 || py <= 0.0 || pz <= 0.0) return t_max;

    alea_bbox_t lat_bbox = {
        .min_x = ll[0], .max_x = ll[0] + ni * px,
        .min_y = ll[1], .max_y = ll[1] + nj * py,
        .min_z = ll[2], .max_z = ll[2] + nk * pz,
    };

    double t_enter, t_exit;
    if (!ray_bbox_slab_enter_exit(ray, &lat_bbox, 0.0, t_max,
                                  &t_enter, &t_exit)) {
        return t_max;
    }

    double sample_t = t_min + RAY_EPSILON;
    double sx = ray->ox + sample_t * ray->dx;
    double sy = ray->oy + sample_t * ray->dy;
    double sz = ray->oz + sample_t * ray->dz;

    int i = (int)floor((sx - ll[0]) / px);
    int j = (nj == 1) ? 0 : (int)floor((sy - ll[1]) / py);
    int k = (nk == 1) ? 0 : (int)floor((sz - ll[2]) / pz);

    if (i < 0) i = 0;
    if (i >= ni) i = ni - 1;
    if (j < 0) j = 0;
    if (j >= nj) j = nj - 1;
    if (k < 0) k = 0;
    if (k >= nk) k = nk - 1;

    (void)t_enter;
    (void)t_exit;
    double t_next = t_max;
    if (fabs(ray->dx) > RAY_EPSILON) {
        int crosses_outer_edge = (ray->dx > 0.0 && i + 1 >= ni) ||
                                 (ray->dx < 0.0 && i <= 0);
        if (!crosses_outer_edge) {
            double boundary = ll[0] + ((ray->dx > 0.0) ? (i + 1) : i) * px;
            double t = (boundary - ray->ox) / ray->dx;
            if (t > t_min + RAY_EPSILON && t < t_next) t_next = t;
        }
    }
    if (nj > 1 && fabs(ray->dy) > RAY_EPSILON) {
        int crosses_outer_edge = (ray->dy > 0.0 && j + 1 >= nj) ||
                                 (ray->dy < 0.0 && j <= 0);
        if (!crosses_outer_edge) {
            double boundary = ll[1] + ((ray->dy > 0.0) ? (j + 1) : j) * py;
            double t = (boundary - ray->oy) / ray->dy;
            if (t > t_min + RAY_EPSILON && t < t_next) t_next = t;
        }
    }
    if (nk > 1 && fabs(ray->dz) > RAY_EPSILON) {
        int crosses_outer_edge = (ray->dz > 0.0 && k + 1 >= nk) ||
                                 (ray->dz < 0.0 && k <= 0);
        if (!crosses_outer_edge) {
            double boundary = ll[2] + ((ray->dz > 0.0) ? (k + 1) : k) * pz;
            double t = (boundary - ray->oz) / ray->dz;
            if (t > t_min + RAY_EPSILON && t < t_next) t_next = t;
        }
    }

    return (t_next < t_max) ? t_next : t_max;
}

static double lattice_hex_next_boundary(const alea_ray_t* ray,
                                        const alea_cell_entry_t* lat_cell,
                                        double t_min,
                                        double t_max) {
    double p = lat_cell->lat_pitch[0];
    if (p <= 0.0) return t_max;

    int imin = lat_cell->lat_fill_dims[0], imax = lat_cell->lat_fill_dims[1];
    int jmin = lat_cell->lat_fill_dims[2], jmax = lat_cell->lat_fill_dims[3];
    int nk = lat_cell->lat_fill_dims[5] - lat_cell->lat_fill_dims[4] + 1;

    double bb_min_x = DBL_MAX, bb_max_x = -DBL_MAX;
    double bb_min_y = DBL_MAX, bb_max_y = -DBL_MAX;
    double r = p / M_SQRT3;

    for (int ci = imin; ci <= imax; ci++) {
        for (int cj = jmin; cj <= jmax; cj++) {
            double cx = ci * p + cj * p * 0.5;
            double cy = cj * p * M_SQRT3 * 0.5;
            if (cx - r < bb_min_x) bb_min_x = cx - r;
            if (cx + r > bb_max_x) bb_max_x = cx + r;
            if (cy - p * 0.5 < bb_min_y) bb_min_y = cy - p * 0.5;
            if (cy + p * 0.5 > bb_max_y) bb_max_y = cy + p * 0.5;
        }
    }

    alea_bbox_t lat_bbox = {
        .min_x = bb_min_x, .max_x = bb_max_x,
        .min_y = bb_min_y, .max_y = bb_max_y,
        .min_z = (nk > 1) ? lat_cell->lat_lower_left[2] : -1e30,
        .max_z = (nk > 1) ? lat_cell->lat_lower_left[2]
                           + nk * lat_cell->lat_pitch[2] : 1e30,
    };

    double t_enter, t_exit;
    if (!ray_bbox_slab_enter_exit(ray, &lat_bbox, 0.0, t_max,
                                  &t_enter, &t_exit)) {
        return t_max;
    }

    double inv_p = 1.0 / p;
    double inv_ps = inv_p / M_SQRT3;
    double dq = ray->dx * inv_p - ray->dy * inv_ps;
    double dr = 2.0 * ray->dy * inv_ps;
    double ds = -dq - dr;

    double sample_t = t_min + RAY_EPSILON;
    double sx = ray->ox + sample_t * ray->dx;
    double sy = ray->oy + sample_t * ray->dy;
    double q = sx * inv_p - sy * inv_ps;
    double rr = 2.0 * sy * inv_ps;
    double s = -q - rr;

    (void)t_enter;
    (void)t_exit;
    double t_next = t_max;
#define CONSIDER_NEXT_HALF(val, rate) do {                         \
        if (fabs(rate) > RAY_EPSILON) {                            \
            double bnd;                                            \
            if ((rate) > 0.0) {                                    \
                bnd = floor((val) + 0.5) + 0.5;                   \
                if (bnd <= (val) + RAY_EPSILON) bnd += 1.0;        \
            } else {                                               \
                bnd = ceil((val) - 0.5) - 0.5;                    \
                if (bnd >= (val) - RAY_EPSILON) bnd -= 1.0;        \
            }                                                      \
            double t = sample_t + (bnd - (val)) / (rate);          \
            if (t > t_min + RAY_EPSILON && t < t_next) t_next = t; \
        }                                                          \
    } while (0)

    CONSIDER_NEXT_HALF(q, dq);
    CONSIDER_NEXT_HALF(rr, dr);
    CONSIDER_NEXT_HALF(s, ds);
#undef CONSIDER_NEXT_HALF

    if (t_next < t_max) {
        double probe_t = t_next + RAY_EPSILON;
        double px = ray->ox + probe_t * ray->dx;
        double py = ray->oy + probe_t * ray->dy;
        double pz = ray->oz + probe_t * ray->dz;
        double ox, oy, oz;
        if (lattice_hex_lookup(lat_cell, px, py, pz, &ox, &oy, &oz) < 0) {
            return t_max;
        }
    }

    return (t_next < t_max) ? t_next : t_max;
}

static double lattice_next_boundary(const alea_ray_t* ray,
                                    const alea_cell_entry_t* lat_cell,
                                    double t_min,
                                    double t_max) {
    if (lat_cell->lat_type == 1)
        return lattice_rect_next_boundary(ray, lat_cell, t_min, t_max);
    if (lat_cell->lat_type == 2)
        return lattice_hex_next_boundary(ray, lat_cell, t_min, t_max);
    return t_max;
}

/**
 * Find the cell at a point, returning cell index or -1 for void.
 */
static int find_cell_at_point(alea_system_t* sys,
                              double px, double py, double pz,
                              int* out_material_id,
                              double* out_density) {
    alea_cell_hit_t hits[32];
    int num_hits = alea_find_all_cells_at_point(sys, px, py, pz, hits, 32);

    if (num_hits > 0) {
        /* Use deepest cell (innermost in hierarchy) */
        const alea_cell_hit_t* hit = &hits[num_hits - 1];
        if (out_material_id) *out_material_id = hit->material_id;
        if (out_density && hit->cell_index >= 0 &&
            (size_t)hit->cell_index < alea_vec_count(&sys->cells)) {
            *out_density = sys->cells.data[hit->cell_index].density;
        }
        return hit->cell_index;
    }

    if (out_material_id) *out_material_id = 0;
    if (out_density) *out_density = 0;
    return -1;
}

static int find_cell_at_point_hier_path(alea_system_t* sys,
                                        double px, double py, double pz,
                                        int* out_material_id,
                                        double* out_density,
                                        alea_matrix_t* out_transform,
                                        int* out_lattice_cell_index,
                                        alea_matrix_t* out_lattice_transform,
                                        alea_hier_ray_path_t* out_path);

static int find_cell_at_point_hier_path(alea_system_t* sys,
                                        double px, double py, double pz,
                                        int* out_material_id,
                                        double* out_density,
                                        alea_matrix_t* out_transform,
                                        int* out_lattice_cell_index,
                                        alea_matrix_t* out_lattice_transform,
                                        alea_hier_ray_path_t* out_path) {
    alea_hier_cell_hit_t hit_with_transform;
    int found = alea_hier_spatial_find_path_at_point(sys, px, py, pz,
                                                     &hit_with_transform,
                                                     out_path);

    if (found > 0) {
        const alea_cell_hit_t* hit = &hit_with_transform.hit;
        if (out_material_id) *out_material_id = hit->material_id;
        if (out_density && hit->cell_index >= 0 &&
            (size_t)hit->cell_index < alea_vec_count(&sys->cells)) {
            *out_density = sys->cells.data[hit->cell_index].density;
        }
        if (out_transform) *out_transform = hit_with_transform.transform;
        if (out_lattice_cell_index)
            *out_lattice_cell_index = hit_with_transform.lattice_cell_index;
        if (out_lattice_transform)
            *out_lattice_transform = hit_with_transform.lattice_transform;
        return hit->cell_index;
    }

    if (out_material_id) *out_material_id = 0;
    if (out_density) *out_density = 0;
    if (out_transform) alea_matrix_identity(out_transform);
    if (out_lattice_cell_index) *out_lattice_cell_index = -1;
    if (out_lattice_transform) alea_matrix_identity(out_lattice_transform);
    return -1;
}

static int find_cell_from_existing_hier_path(alea_system_t* sys,
                                             const alea_hier_ray_path_t* current_path,
                                             double px, double py, double pz,
                                             int* out_material_id,
                                             double* out_density,
                                             alea_matrix_t* out_transform,
                                             int* out_lattice_cell_index,
                                             alea_matrix_t* out_lattice_transform,
                                             alea_hier_ray_path_t* out_path) {
    if (!sys || !current_path || current_path->count <= 1) return -1;

    for (int parent = current_path->count - 2; parent >= 0; parent--) {
        alea_hier_cell_hit_t hit_with_transform;
        alea_hier_ray_path_t candidate_path;
        int found = alea_hier_spatial_find_path_from_parent(
            sys, current_path, parent, px, py, pz,
            &hit_with_transform, &candidate_path);
        if (found < 0) return -2;
        if (found == 0) continue;

        const alea_cell_hit_t* hit = &hit_with_transform.hit;
        if (out_material_id) *out_material_id = hit->material_id;
        if (out_density && hit->cell_index >= 0 &&
            (size_t)hit->cell_index < alea_vec_count(&sys->cells)) {
            *out_density = sys->cells.data[hit->cell_index].density;
        }
        if (out_transform) *out_transform = hit_with_transform.transform;
        if (out_lattice_cell_index) {
            *out_lattice_cell_index = hit_with_transform.lattice_cell_index;
        }
        if (out_lattice_transform) {
            *out_lattice_transform = hit_with_transform.lattice_transform;
        }
        if (out_path) *out_path = candidate_path;
        return hit->cell_index;
    }

    return -1;
}

static int find_cell_from_root_universe(alea_system_t* sys,
                                        double px, double py, double pz,
                                        int crossed_surface_id,
                                        int* out_material_id,
                                        double* out_density,
                                        alea_matrix_t* out_transform,
                                        int* out_lattice_cell_index,
                                        alea_matrix_t* out_lattice_transform,
                                        alea_hier_ray_path_t* out_path) {
    if (!sys) return -2;

    int root_cell = alea_hier_spatial_find_ordered_cell_in_universe(
        sys, 0, px, py, pz, crossed_surface_id);
    if (root_cell < 0) return root_cell == -2 ? -2 : -1;
    if ((size_t)root_cell >= alea_vec_count(&sys->cells)) return -2;

    const alea_cell_entry_t* cell = &sys->cells.data[root_cell];
    if (cell->lat_type != 0 && cell->lat_fill) {
        return -1;
    }

    alea_hier_ray_path_t root_path;
    root_path.count = 1;
    alea_hier_ray_path_entry_t* ent = &root_path.entries[0];
    ent->cell_index = (uint32_t)root_cell;
    ent->cell_id = cell->mc_cell_id;
    ent->material_id = cell->material_id;
    ent->universe_id = cell->universe_id;
    ent->fill_universe = cell->fill_universe;
    ent->depth = 0;
    ent->is_lattice = 0;
    ent->lat_fill_universe = 0;
    ent->lat_ox = 0.0;
    ent->lat_oy = 0.0;
    ent->lat_oz = 0.0;
    alea_matrix_identity(&ent->transform);

    if (cell->fill_universe > 0) {
        alea_hier_cell_hit_t hit_with_transform;
        alea_hier_ray_path_t candidate_path;
        int found = alea_hier_spatial_find_path_from_parent(
            sys, &root_path, 0, px, py, pz,
            &hit_with_transform, &candidate_path);
        if (found < 0) return -2;
        if (found > 0) {
            const alea_cell_hit_t* hit = &hit_with_transform.hit;
            if (out_material_id) *out_material_id = hit->material_id;
            if (out_density && hit->cell_index >= 0 &&
                (size_t)hit->cell_index < alea_vec_count(&sys->cells)) {
                *out_density = sys->cells.data[hit->cell_index].density;
            }
            if (out_transform) *out_transform = hit_with_transform.transform;
            if (out_lattice_cell_index) {
                *out_lattice_cell_index = hit_with_transform.lattice_cell_index;
            }
            if (out_lattice_transform) {
                *out_lattice_transform = hit_with_transform.lattice_transform;
            }
            if (out_path) *out_path = candidate_path;
            return hit->cell_index;
        }
    }

    if (out_material_id) *out_material_id = cell->material_id;
    if (out_density) *out_density = cell->density;
    if (out_transform) alea_matrix_identity(out_transform);
    if (out_lattice_cell_index) *out_lattice_cell_index = -1;
    if (out_lattice_transform) alea_matrix_identity(out_lattice_transform);
    if (out_path) *out_path = root_path;
    return root_cell;
}

/**
 * Context for finding closest intersection during BVH traversal.
 */
typedef struct {
    alea_system_t* sys;
    const alea_ray_t* ray;
    double t_min;
    double closest_t;
    int closest_surface_id;
    uint32_t closest_primitive_id;
} closest_hit_ctx_t;

/**
 * BVH callback that finds the closest intersection.
 */
static void find_closest_callback(uint32_t surface_idx, void* userdata) {
    closest_hit_ctx_t* ctx = (closest_hit_ctx_t*)userdata;
    const alea_surface_entry_t* surf = &ctx->sys->surfaces.data[surface_idx];
    const alea_primitive_entry_t* prim = &ctx->sys->primitives.data[surf->primitive_id];
    alea_primitive_data_t prim_data;
    if (!raycast_primitive_copy_payload(ctx->sys, surf->primitive_id,
                                        prim->type, &prim_data)) return;

    double t[4];
    int count = ray_intersect_primitive(ctx->ray, prim->type, &prim_data, t);

    for (int j = 0; j < count; j++) {
        if (t[j] > ctx->t_min && t[j] < ctx->closest_t) {
            ctx->closest_t = t[j];
            ctx->closest_surface_id = surf->mc_surface_id;
            ctx->closest_primitive_id = surf->primitive_id;
        }
    }
}

/**
 * Find the closest surface intersection using BVH or linear scan.
 */
static double find_closest_intersection(alea_system_t* sys,
                                        const alea_ray_t* ray,
                                        double t_min, double t_max,
                                        int* out_surface_id,
                                        uint32_t* out_primitive_id) {
    *out_surface_id = -1;
    if (out_primitive_id) *out_primitive_id = ALEA_PRIMITIVE_ID_INVALID;

    /* Use BVH if available */
#if !BVH_DISABLED
    if (sys->surface_bvh) {
        closest_hit_ctx_t ctx = {
            .sys = sys,
            .ray = ray,
            .t_min = t_min,
            .closest_t = t_max,
            .closest_surface_id = -1,
            .closest_primitive_id = ALEA_PRIMITIVE_ID_INVALID
        };
        alea_bvh_traverse(sys->surface_bvh, ray, t_min, t_max,
                         find_closest_callback, &ctx);
        *out_surface_id = ctx.closest_surface_id;
        if (out_primitive_id) *out_primitive_id = ctx.closest_primitive_id;
        return ctx.closest_t;
    }
#endif

    /* Fallback: linear scan */
    double closest_t = t_max;
    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        const alea_surface_entry_t* surf = &sys->surfaces.data[i];
        const alea_primitive_entry_t* prim = &sys->primitives.data[surf->primitive_id];
        alea_primitive_data_t prim_data;
        if (!raycast_primitive_copy_payload(sys, surf->primitive_id,
                                            prim->type, &prim_data)) continue;

        double t[4];
        int count = ray_intersect_primitive(ray, prim->type, &prim_data, t);

        for (int j = 0; j < count; j++) {
            if (t[j] > t_min && t[j] < closest_t) {
                closest_t = t[j];
                *out_surface_id = surf->mc_surface_id;
                if (out_primitive_id) *out_primitive_id = surf->primitive_id;
            }
        }
    }
    return closest_t;
}

static int raycast_cell_aware_impl(alea_system_t* sys,
                                   const alea_ray_t* ray,
                                   double effective_t_max,
                                   bool use_hier_lookup,
                                   bool emit_hits,
                                   alea_raycast_result_t* result) {
    /* Current position along ray */
    double t_current = 0;
    int prev_cell_idx = -2;  /* -2 = no previous */
    int prev_surface_id = -1;  /* Surface ID we're about to cross */
    /* Hit index of the surface crossed at t_current (the next segment's
     * enter boundary), or -1. Only maintained when emit_hits is set. */
    int pending_enter_hit_index = -1;
    alea_hier_ray_path_t current_path;
    current_path.count = 0;

    /* Safety limit to prevent infinite loops */
    int max_iterations = 10000;

    while (t_current < effective_t_max && max_iterations-- > 0) {
        result->step_iterations++;
        int cell_idx = -1;
        int cell_id = -1;
        int material_id = 0;
        double density = 0;
        int found_via_neighbor = 0;
        alea_matrix_t cell_transform;
        alea_matrix_t lattice_transform;
        int lattice_cell_index = -1;
        alea_matrix_identity(&cell_transform);
        alea_matrix_identity(&lattice_transform);

        /* Try neighbor lookup first if we have previous cell and surface info */
        if (prev_cell_idx >= 0 && prev_surface_id > 0) {
            const alea_cell_entry_t* prev_cell = &sys->cells.data[prev_cell_idx];
            if (!use_hier_lookup || prev_cell->universe_id == 0) {
                found_via_neighbor = raycast_find_neighbor(sys, prev_cell_idx,
                                                           prev_surface_id,
                                                           &cell_id, &cell_idx,
                                                           &material_id,
                                                           &density);
            }
        }

        /* In hierarchical mode, avoid a full root-to-deepest lookup when the
         * sampled point remains inside the cached placement path. */
        if (!found_via_neighbor && use_hier_lookup && prev_cell_idx >= 0 &&
            (size_t)prev_cell_idx < alea_vec_count(&sys->cells)) {
            const alea_cell_entry_t* prev_cell = &sys->cells.data[prev_cell_idx];
            double px, py, pz;
            alea_ray_point_at(ray, t_current + RAY_EPSILON, &px, &py, &pz);
            int in_cell = 0;
            if (prev_cell->universe_id == 0) {
                in_cell = prev_cell->root_node_id != ALEA_NODE_ID_INVALID &&
                    alea_contains_point(sys, prev_cell->root_node_id,
                                        px, py, pz);
                alea_matrix_identity(&cell_transform);
                lattice_cell_index = -1;
                alea_matrix_identity(&lattice_transform);
            } else {
                in_cell = alea_hier_spatial_check_path_containment(
                    sys, &current_path, (uint32_t)prev_cell_idx, px, py, pz,
                    &cell_transform, &lattice_cell_index, &lattice_transform);
                if (in_cell < 0) in_cell = 0;
            }
            if (in_cell) {
                cell_idx = prev_cell_idx;
                cell_id = prev_cell->mc_cell_id;
                material_id = prev_cell->material_id;
                density = prev_cell->density;
                found_via_neighbor = 1;
            }
        }

        /* Fall back to full lookup if neighbor lookup failed or not available */
        if (!found_via_neighbor) {
            double px, py, pz;
            alea_ray_point_at(ray, t_current + RAY_EPSILON, &px, &py, &pz);
            if (use_hier_lookup) {
                cell_idx = find_cell_from_existing_hier_path(
                    sys, &current_path, px, py, pz,
                    &material_id, &density,
                    &cell_transform, &lattice_cell_index,
                    &lattice_transform, &current_path);
                if (cell_idx == -2) return -1;
            }
            if (use_hier_lookup && cell_idx < 0) {
                cell_idx = find_cell_from_root_universe(
                    sys, px, py, pz, prev_surface_id,
                    &material_id, &density,
                    &cell_transform, &lattice_cell_index,
                    &lattice_transform, &current_path);
                if (cell_idx == -2) return -1;
            }
            if (cell_idx < 0) {
                result->point_lookups++;
                cell_idx = use_hier_lookup
                    ? find_cell_at_point_hier_path(
                        sys, px, py, pz, &material_id, &density,
                        &cell_transform, &lattice_cell_index,
                        &lattice_transform, &current_path)
                    : find_cell_at_point(sys, px, py, pz,
                                         &material_id, &density);
                if (use_hier_lookup && cell_idx < 0) current_path.count = 0;
            }
            if (cell_idx >= 0 && (size_t)cell_idx < alea_vec_count(&sys->cells)) {
                cell_id = sys->cells.data[cell_idx].mc_cell_id;
            }
        }

        /* Find next surface crossing:
         * Use per-cell surface index when inside a cell (the whole point of cell-aware),
         * fall back to global search only for void regions. */
        int hit_surface_id = -1;
        int next_enter_surface_id = -1;
        double t_next;
        double t_lattice_next = DBL_MAX;
        /* Winning boundary event for this step. Only populated/consumed when
         * emit_hits is set; the segment-only path never touches it, so it adds
         * no work to the hot loop. */
        alea_raycast_boundary_event_t bevent;
        if (emit_hits) {
            bevent.t = 0;
            bevent.surface_id = -1;
            bevent.primitive_id = ALEA_PRIMITIVE_ID_INVALID;
            bevent.has_physical_surface = false;
            bevent.is_synthetic_lattice_boundary = false;
            alea_matrix_identity(&bevent.transform);
        }
        if (cell_idx >= 0 && (size_t)cell_idx < alea_vec_count(&sys->cells) &&
            sys->cells.data[cell_idx].surface_indices) {
            alea_ray_t local_ray;
            const alea_ray_t* surface_ray = ray;
            if (use_hier_lookup) {
                transform_ray_inverse(&cell_transform, ray, &local_ray);
                surface_ray = &local_ray;
            }
            /* Instrumentation: attribute this cell's surface tests and record
             * its size. Snapshotting surfaces_tested keeps raycast_cell_surfaces
             * free of a category parameter. */
            uint32_t terminal_surf_count =
                sys->cells.data[cell_idx].surface_index_count;
            result->crossed_cell_count++;
            result->sum_cell_surface_count += terminal_surf_count;
            if (terminal_surf_count > result->max_cell_surface_count)
                result->max_cell_surface_count = terminal_surf_count;
            int surfaces_tested_before = result->surfaces_tested;

            uint32_t terminal_prim_id = ALEA_PRIMITIVE_ID_INVALID;
            t_next = raycast_cell_surfaces(sys, surface_ray,
                                           &sys->cells.data[cell_idx],
                                           t_current + RAY_EPSILON, effective_t_max,
                                           result,
                                           &hit_surface_id,
                                           &terminal_prim_id);
            result->terminal_surfaces_tested +=
                result->surfaces_tested - surfaces_tested_before;
            if (emit_hits) {
                bevent.t = t_next;
                bevent.surface_id = hit_surface_id;
                bevent.primitive_id = terminal_prim_id;
                bevent.has_physical_surface = (hit_surface_id > 0);
                bevent.is_synthetic_lattice_boundary = false;
                bevent.transform = cell_transform;
            }
            if (use_hier_lookup &&
                lattice_cell_index >= 0 &&
                (size_t)lattice_cell_index < alea_vec_count(&sys->cells)) {
                const alea_cell_entry_t* lattice_cell =
                    &sys->cells.data[lattice_cell_index];
                alea_ray_t lattice_ray;
                transform_ray_inverse(&lattice_transform, ray, &lattice_ray);
                int lattice_surface_id = -1;
                uint32_t lattice_prim_id = ALEA_PRIMITIVE_ID_INVALID;
                int lattice_tested_before = result->surfaces_tested;
                double t_lattice_surface =
                    raycast_cell_surfaces(sys, &lattice_ray, lattice_cell,
                                          t_current + RAY_EPSILON,
                                          effective_t_max,
                                          result,
                                          &lattice_surface_id,
                                          &lattice_prim_id);
                result->lattice_surfaces_tested +=
                    result->surfaces_tested - lattice_tested_before;
                t_lattice_next = lattice_next_boundary(&lattice_ray, lattice_cell,
                                                       t_current + RAY_EPSILON,
                                                       effective_t_max);
                if (t_lattice_surface < t_next - RAY_EPSILON) {
                    t_next = t_lattice_surface;
                    hit_surface_id = lattice_surface_id;
                    if (emit_hits) {
                        bevent.t = t_lattice_surface;
                        bevent.surface_id = lattice_surface_id;
                        bevent.primitive_id = lattice_prim_id;
                        bevent.has_physical_surface = (lattice_surface_id > 0);
                        bevent.is_synthetic_lattice_boundary = false;
                        bevent.transform = lattice_transform;
                    }
                }
                if (t_lattice_next < effective_t_max - RAY_EPSILON &&
                    t_lattice_next < t_next - RAY_EPSILON) {
                    t_next = t_lattice_next;
                    if (fabs(t_lattice_surface - t_lattice_next) <= RAY_EPSILON) {
                        hit_surface_id = lattice_surface_id;
                        next_enter_surface_id = 0;
                        if (emit_hits) {
                            bevent.t = t_lattice_next;
                            bevent.surface_id = lattice_surface_id;
                            bevent.primitive_id = lattice_prim_id;
                            bevent.has_physical_surface = (lattice_surface_id > 0);
                            bevent.is_synthetic_lattice_boundary = true;
                            bevent.transform = lattice_transform;
                        }
                    } else {
                        hit_surface_id = 0;
                        if (emit_hits) {
                            /* Pure synthetic DDA boundary: no physical surface. */
                            bevent.t = t_lattice_next;
                            bevent.surface_id = 0;
                            bevent.primitive_id = ALEA_PRIMITIVE_ID_INVALID;
                            bevent.has_physical_surface = false;
                            bevent.is_synthetic_lattice_boundary = true;
                        }
                    }
                } else if (fabs(t_lattice_surface - t_lattice_next) <= RAY_EPSILON &&
                           fabs(t_lattice_next - t_next) <= RAY_EPSILON &&
                           lattice_surface_id >= 0) {
                    hit_surface_id = lattice_surface_id;
                    next_enter_surface_id = 0;
                    if (emit_hits) {
                        bevent.surface_id = lattice_surface_id;
                        bevent.primitive_id = lattice_prim_id;
                        bevent.has_physical_surface = (lattice_surface_id > 0);
                        bevent.is_synthetic_lattice_boundary = true;
                        bevent.transform = lattice_transform;
                    }
                }
            }
            if (use_hier_lookup && current_path.count > 1) {
                int ancestor_surface_id = -1;
                uint32_t ancestor_prim_id = ALEA_PRIMITIVE_ID_INVALID;
                alea_matrix_t ancestor_transform;
                alea_matrix_identity(&ancestor_transform);
                int ancestor_tested_before = result->surfaces_tested;
                double t_ancestor = raycast_hier_path_ancestor_surfaces(
                    sys, ray, &current_path, (uint32_t)cell_idx,
                    lattice_cell_index, t_current + RAY_EPSILON,
                    effective_t_max, result, &ancestor_surface_id,
                    &ancestor_prim_id, &ancestor_transform);
                result->ancestor_surfaces_tested +=
                    result->surfaces_tested - ancestor_tested_before;
                if (t_ancestor < t_next - RAY_EPSILON) {
                    t_next = t_ancestor;
                    hit_surface_id = ancestor_surface_id;
                    next_enter_surface_id = ancestor_surface_id;
                    if (emit_hits) {
                        bevent.t = t_ancestor;
                        bevent.surface_id = ancestor_surface_id;
                        bevent.primitive_id = ancestor_prim_id;
                        bevent.has_physical_surface = (ancestor_surface_id > 0);
                        bevent.is_synthetic_lattice_boundary = false;
                        bevent.transform = ancestor_transform;
                    }
                }
            }
        } else {
            uint32_t void_prim_id = ALEA_PRIMITIVE_ID_INVALID;
            t_next = find_closest_intersection(sys, ray,
                                               t_current + RAY_EPSILON, effective_t_max,
                                               &hit_surface_id, &void_prim_id);
            if (emit_hits) {
                /* Void-region global search runs in the world frame. */
                bevent.t = t_next;
                bevent.surface_id = hit_surface_id;
                bevent.primitive_id = void_prim_id;
                bevent.has_physical_surface = (hit_surface_id > 0);
                bevent.is_synthetic_lattice_boundary = false;
                alea_matrix_identity(&bevent.transform);
            }
        }

        /* Ensure we make progress at any scale:
         * absolute 1e-10 dominates near origin; relative 1e-6 at large t */
        if (t_next <= t_current + RAY_EPSILON) {
            t_next = t_current * (1.0 + 1e-6) + RAY_EPSILON;
        }

        /* raycast_cell_surfaces returns DBL_MAX when no hit lies before t_max;
         * clamp so the terminal segment stops at the user's max_distance. */
        if (t_next > effective_t_max) t_next = effective_t_max;

        /* Emit a boundary hit for this step's physical surface crossing (if any)
         * before recording the segment, so the segment can reference it. The
         * progress-guard and t_max clamps only fire when no physical surface was
         * found (has_physical_surface == false), so they never spawn a hit. */
        int exit_hit_index = -1;
        if (emit_hits) {
            bevent.t = t_next;
            exit_hit_index = boundary_event_emit_hit(sys, ray, &bevent, result);
        }

        /* Add or extend segment */
        if (cell_idx == prev_cell_idx && result->segments.count > 0) {
            /* Extend previous segment */
            alea_ray_segment_t* prev_seg =
                &result->segments.data[result->segments.count - 1];
            prev_seg->t_exit = t_next;
            prev_seg->exit_surface_id = hit_surface_id;
        } else {
            /* Create new segment */
            alea_ray_segment_t seg;
            seg.t_enter = t_current;
            seg.t_exit = t_next;
            seg.cell_id = cell_id;
            seg.material_id = material_id;
            seg.density = density;
            seg.enter_surface_id = prev_surface_id;
            seg.exit_surface_id = hit_surface_id;
            seg.enter_hit_index = pending_enter_hit_index;
            add_segment(result, &seg);
            prev_cell_idx = cell_idx;
        }

        /* The hit at t_next is the enter boundary of the next segment. */
        if (emit_hits) pending_enter_hit_index = exit_hit_index;

        if (next_enter_surface_id < 0)
            next_enter_surface_id = hit_surface_id;
        prev_surface_id = next_enter_surface_id;

        /* Move past the intersection */
        t_current = t_next;

        /* If we hit nothing, we're done */
        if (t_next >= effective_t_max - RAY_EPSILON) {
            break;
        }
    }

    return 0;
}

int alea_raycast_cell_aware(alea_system_t* sys,
                           double ox, double oy, double oz,
                           double dx, double dy, double dz,
                           double t_max,
                           alea_raycast_result_t* result) {
    if (!sys || !result) return -1;

    if (system_has_lattice_cells(sys)) {
        /* Lattice transport requires synthetic DDA boundary hits and
         * element-local universe raycasts. Use the canonical lattice-aware
         * pipeline so this public entry point cannot diverge semantically. */
        return alea_raycast(sys, ox, oy, oz, dx, dy, dz, t_max, result);
    }

    alea_raycast_result_free(result);

    alea_ray_t ray;
    if (alea_ray_init(&ray, ox, oy, oz, dx, dy, dz) != 0) {
        return -1;
    }
    result->ray = ray;

    double effective_t_max = (t_max <= 0) ? DBL_MAX : t_max;

    if (!alea_system_query_cache_ready(sys, ALEA_CACHE_RAYCAST)) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              "raycast caches are not prepared; call alea_prepare_query_acceleration()");
        return -1;
    }

    return raycast_cell_aware_impl(sys, &ray, effective_t_max, false, false,
                                   result);
}

int alea_raycast_hier_cell_aware(alea_system_t* sys,
                                 double ox, double oy, double oz,
                                 double dx, double dy, double dz,
                                 double t_max,
                                 alea_raycast_result_t* result) {
    return alea_raycast_hier_fast_segments(sys, ox, oy, oz,
                                           dx, dy, dz, t_max,
                                           result);
}

int alea_raycast_hier_fast_segments(alea_system_t* sys,
                                    double ox, double oy, double oz,
                                    double dx, double dy, double dz,
                                    double t_max,
                                    alea_raycast_result_t* result) {
    if (!sys || !result) return -1;

    if (alea_system_prepare_query_caches(sys,
            ALEA_CACHE_HIER_SPATIAL | ALEA_CACHE_CELL_SURFACES) != 0)
        return -1;

    alea_raycast_result_free(result);

    alea_ray_t ray;
    if (alea_ray_init(&ray, ox, oy, oz, dx, dy, dz) != 0) {
        return -1;
    }
    result->ray = ray;

    double effective_t_max = (t_max <= 0) ? DBL_MAX : t_max;
    return raycast_cell_aware_impl(sys, &ray, effective_t_max, true, false,
                                   result);
}

int alea_raycast_hier_with_hits(alea_system_t* sys,
                                double ox, double oy, double oz,
                                double dx, double dy, double dz,
                                double t_max,
                                alea_raycast_result_t* result) {
    if (!sys || !result) return -1;

    if (alea_system_prepare_query_caches(sys,
            ALEA_CACHE_HIER_SPATIAL | ALEA_CACHE_CELL_SURFACES) != 0)
        return -1;

    alea_raycast_result_free(result);

    alea_ray_t ray;
    if (alea_ray_init(&ray, ox, oy, oz, dx, dy, dz) != 0) {
        return -1;
    }
    result->ray = ray;

    double effective_t_max = (t_max <= 0) ? DBL_MAX : t_max;
    return raycast_cell_aware_impl(sys, &ray, effective_t_max, true, true,
                                   result);
}

/* Buffer-reuse hierarchical variants: take a pre-normalized ray, assume query
 * caches are already built (ALEA_CACHE_RAYCAST), and do NOT free the result
 * buffers. The caller must clear the result (count reset, capacity retained)
 * between rays. These are the hot-loop entry points for renderers that trace
 * one ray per pixel and must avoid per-pixel malloc/free churn. */
int alea_raycast_hier_with_hits_nocache(alea_system_t* sys,
                                        const alea_ray_t* ray,
                                        double t_max,
                                        alea_raycast_result_t* result) {
    if (!sys || !ray || !result) return -1;
    double effective_t_max = (t_max <= 0) ? DBL_MAX : t_max;
    result->ray = *ray;
    return raycast_cell_aware_impl(sys, ray, effective_t_max, true, true,
                                   result);
}

int alea_raycast_hier_segments_nocache(alea_system_t* sys,
                                       const alea_ray_t* ray,
                                       double t_max,
                                       alea_raycast_result_t* result) {
    if (!sys || !ray || !result) return -1;
    double effective_t_max = (t_max <= 0) ? DBL_MAX : t_max;
    result->ray = *ray;
    return raycast_cell_aware_impl(sys, ray, effective_t_max, true, false,
                                   result);
}
