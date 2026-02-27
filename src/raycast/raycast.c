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
#include "ray_intersect.h"
#include "ray_epsilon.h"
#include "ray_bbox.h"
#include "bvh.h"
#include "core/alea_system.h"
#include "core/alea_universe.h"
#include "core/alea_spatial.h"
#include "core/alea_eval.h"
#include "util/alea_log.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "util/math.h"
#include <stdio.h>

#define INITIAL_CAPACITY 32

/* ============================================================================
 * CACHE PRE-BUILD (thread safety)
 *
 * All shared caches (surface BVH, spatial index, cell adjacency) must be
 * built before any concurrent access.  Call once from each public entry
 * point so that OpenMP threads never race on lazy initialisation.
 * ============================================================================ */

int alea_raycast_ensure_caches(const alea_system_t* sys) {
    alea_system_t* m = (alea_system_t*)sys;
    int did_build = 0;

#if !BVH_DISABLED
    if (m->bvh_dirty ||
        (m->surface_bvh &&
         m->surface_bvh->surface_count != alea_vec_count(&sys->surfaces))) {
        if (m->surface_bvh) {
            alea_bvh_free(m->surface_bvh);
            m->surface_bvh = NULL;
        }
        if (alea_vec_count(&sys->surfaces) > 0) {
            m->surface_bvh = alea_bvh_build(sys);
        }
        m->bvh_dirty = false;
        did_build = 1;
    }
#endif

    if (!m->spatial_index || !m->spatial_index->built) {
        alea_spatial_index_build(m);
        did_build = 1;
    }

    if (!m->cell_adjacency_built) {
        did_build = 1;
    }
    alea_build_cell_adjacency(m);

    return did_build;
}

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

void alea_raycast_result_free(alea_raycast_result_t* result) {
    free(result->hits);
    free(result->segments);
    memset(result, 0, sizeof(*result));
}

void alea_raycast_result_clear(alea_raycast_result_t* result) {
    result->hit_count = 0;
    result->segment_count = 0;
    result->surfaces_tested = 0;
    result->bbox_culled = 0;
}

void alea_raycast_result_reserve(alea_raycast_result_t* result,
                                size_t hit_cap, size_t seg_cap) {
    if (hit_cap > result->hit_capacity) {
        alea_ray_hit_t* new_hits = realloc(result->hits, hit_cap * sizeof(alea_ray_hit_t));
        if (new_hits) {
            result->hits = new_hits;
            result->hit_capacity = hit_cap;
        }
    }
    if (seg_cap > result->segment_capacity) {
        alea_ray_segment_t* new_segs = realloc(result->segments, seg_cap * sizeof(alea_ray_segment_t));
        if (new_segs) {
            result->segments = new_segs;
            result->segment_capacity = seg_cap;
        }
    }
}

static int add_hit(alea_raycast_result_t* result, const alea_ray_hit_t* hit) {
    if (result->hit_count >= result->hit_capacity) {
        size_t new_cap = result->hit_capacity ? result->hit_capacity * 2 : INITIAL_CAPACITY;
        alea_ray_hit_t* new_hits = realloc(result->hits, new_cap * sizeof(alea_ray_hit_t));
        if (!new_hits) return -1;
        result->hits = new_hits;
        result->hit_capacity = new_cap;
    }
    result->hits[result->hit_count++] = *hit;
    return 0;
}

static int add_segment(alea_raycast_result_t* result, const alea_ray_segment_t* seg) {
    if (result->segment_count >= result->segment_capacity) {
        size_t new_cap = result->segment_capacity ? result->segment_capacity * 2 : INITIAL_CAPACITY;
        alea_ray_segment_t* new_segs = realloc(result->segments, new_cap * sizeof(alea_ray_segment_t));
        if (!new_segs) return -1;
        result->segments = new_segs;
        result->segment_capacity = new_cap;
    }
    result->segments[result->segment_count++] = *seg;
    return 0;
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
    const alea_system_t* sys;
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

        ctx->result->surfaces_tested++;

        /* Find intersections */
        double t[4];
        int n = ray_intersect_primitive(ctx->ray, prim->type, &prim->data, t);

        /* Add valid hits */
        for (int j = 0; j < n; j++) {
            if (t[j] >= ctx->t_min && t[j] <= ctx->t_max) {
                alea_ray_hit_t hit;
                hit.t = t[j];
                hit.surface_id = surf->mcnp_surface_id;

                /* Compute normal at hit point */
                double px, py, pz;
                alea_ray_point_at(ctx->ray, t[j], &px, &py, &pz);
                primitive_normal_at(prim->type, &prim->data, px, py, pz,
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
static int raycast_surfaces_linear(const alea_system_t* sys,
                                   const alea_ray_t* ray,
                                   double t_min, double t_max,
                                   alea_raycast_result_t* result) {
    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        
        const alea_surface_entry_t* surf = &sys->surfaces.data[i];
        const alea_primitive_entry_t* prim = &sys->primitives.data[surf->primitive_id];

        result->surfaces_tested++;

        double t[4];
        int count = ray_intersect_primitive(ray, prim->type, &prim->data, t);

        for (int j = 0; j < count; j++) {
            if (t[j] >= t_min && t[j] <= t_max) {
                alea_ray_hit_t hit;
                hit.t = t[j];
                hit.surface_id = surf->mcnp_surface_id;

                double px, py, pz;
                alea_ray_point_at(ray, t[j], &px, &py, &pz);
                primitive_normal_at(prim->type, &prim->data, px, py, pz,
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
static int raycast_surfaces_impl(const alea_system_t* sys,
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
    sort_hits(result->hits, result->hit_count);

    /* Remove duplicate hits (same t AND same surface_id within epsilon) */
    if (result->hit_count > 1) {
        size_t write = 1;
        for (size_t read = 1; read < result->hit_count; read++) {
            int same_t = fabs(result->hits[read].t - result->hits[write - 1].t) <= DEDUP_EPSILON;
            int same_surf = result->hits[read].surface_id == result->hits[write - 1].surface_id;
            if (!(same_t && same_surf)) {
                result->hits[write++] = result->hits[read];
            }
        }
        result->hit_count = write;
    }

    return 0;
}

int alea_raycast_surfaces(const alea_system_t* sys,
                         const alea_ray_t* ray,
                         double t_min, double t_max,
                         alea_raycast_result_t* result) {
    if (!sys || !ray || !result) return -1;

    if (alea_raycast_ensure_caches(sys)) {
        ALEA_LOG_WARN("Lazy-building raycast caches. Call alea_build_spatial_index() "
                     "before concurrent raycast calls to avoid data races.");
    }

    return raycast_surfaces_impl(sys, ray, t_min, t_max, result);
}

int alea_raycast_surfaces_nocache(const alea_system_t* sys,
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
static int raycast_find_neighbor(const alea_system_t* sys,
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
    *out_cell_id = neighbor->mcnp_cell_id;
    *out_cell_idx = neighbor_idx;
    *out_material_id = neighbor->material_id;
    *out_density = neighbor->density;
    return 1;
}

/**
 * @brief Full cell lookup via point-in-cell search
 */
static void raycast_find_cell_full(const alea_system_t* sys,
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

/**
 * 3-tier cell lookup after crossing a surface:
 * 1) Neighbor lookup (O(1) via adjacency)
 * 2) Coherence check (O(tree_depth) — is point still in previous cell?)
 * 3) Full point-in-cell search (O(n_candidates) via spatial index)
 */
static void find_cell_after_crossing(const alea_system_t* sys,
                                     const alea_ray_t* ray,
                                     double t_prev, double t_curr,
                                     int prev_cell_idx,
                                     int crossed_surface_id,
                                     int* out_cell_id, int* out_cell_idx,
                                     int* out_material_id, double* out_density) {
    *out_cell_id = -1;
    *out_cell_idx = -1;
    *out_material_id = 0;
    *out_density = 0;
    int found = 0;

    /* Tier 1: neighbor lookup */
    if (prev_cell_idx >= 0 && crossed_surface_id > 0) {
        found = raycast_find_neighbor(sys, prev_cell_idx, crossed_surface_id,
                                      out_cell_id, out_cell_idx,
                                      out_material_id, out_density);
    }

    int crossed_dda = (crossed_surface_id == 0);

    /* Tier 2: coherence check (skip for DDA boundaries) */
    if (!found && prev_cell_idx >= 0 && !crossed_dda) {
        double t_sample = t_prev + fmin(0.5 * (t_curr - t_prev), SURFACE_SAMPLE_OFFSET);
        double px, py, pz;
        alea_ray_point_at(ray, t_sample, &px, &py, &pz);
        const alea_cell_entry_t* prev_cell = &sys->cells.data[prev_cell_idx];
        if (alea_contains_point(sys, prev_cell->root_node_id, px, py, pz)) {
            *out_cell_id = prev_cell->mcnp_cell_id;
            *out_cell_idx = prev_cell_idx;
            *out_material_id = prev_cell->material_id;
            *out_density = prev_cell->density;
            found = 1;
        }
    }

    /* Tier 3: full lookup */
    if (!found) {
        double max_offset = crossed_dda ? DDA_SAMPLE_OFFSET : SURFACE_SAMPLE_OFFSET;
        double t_sample = t_prev + fmin(0.5 * (t_curr - t_prev), max_offset);
        double px, py, pz;
        alea_ray_point_at(ray, t_sample, &px, &py, &pz);
        raycast_find_cell_full(sys, px, py, pz,
                               out_cell_id, out_cell_idx, out_material_id, out_density);
    }
}

int alea_raycast_to_segments(const alea_system_t* sys,
                            alea_raycast_result_t* result) {
    if (!sys || !result) return -1;

    result->segment_count = 0;

    const alea_ray_t* ray = &result->ray;
    double t_prev = 0;
    int prev_cell_id = -2;  /* Use -2 as "no previous" since -1 is valid (void) */
    int prev_cell_idx = -1; /* Track cell index for neighbor lookup */

    /* Process intervals between hits */
    for (size_t i = 0; i <= result->hit_count; i++) {
        double t_curr = (i < result->hit_count) ? result->hits[i].t : DBL_MAX;

        /* Only process if there's a real interval */
        if (t_curr > t_prev + RAY_EPSILON) {
            int cell_id, cell_idx, material_id;
            double density;
            int crossed_surface = (i > 0) ? result->hits[i - 1].surface_id : -1;
            find_cell_after_crossing(sys, ray, t_prev, t_curr,
                                     prev_cell_idx, crossed_surface,
                                     &cell_id, &cell_idx, &material_id, &density);

            /* Update tracking for next iteration */
            prev_cell_idx = cell_idx;

            /* Extend previous segment or start new one */
            if (cell_id == prev_cell_id && result->segment_count > 0) {
                /* Extend */
                result->segments[result->segment_count - 1].t_exit = t_curr;
            } else {
                /* New segment */
                alea_ray_segment_t seg;
                seg.t_enter = t_prev;
                seg.t_exit = t_curr;
                seg.cell_id = cell_id;
                seg.material_id = material_id;
                seg.density = density;
                seg.enter_hit_index = (i > 0) ? (int)(i - 1) : -1;

                add_segment(result, &seg);
                prev_cell_id = cell_id;
            }
        }

        t_prev = t_curr;
    }

    return 0;
}

/* ============================================================================
 * LATTICE DDA RAYCAST
 * ============================================================================ */

/**
 * Walk a CSG tree and test ray intersection with every primitive leaf.
 * Adds hits to result. Duplicates from shared surfaces are handled by
 * the global dedup pass after all hits are collected.
 */
static void raycast_tree_primitives(const alea_system_t* sys,
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
        double t[4];
        int count = ray_intersect_primitive(ray, prim->type, &prim->data, t);
        for (int j = 0; j < count; j++) {
            if (t[j] >= t_min && t[j] <= t_max) {
                alea_ray_hit_t hit;
                hit.t = t[j];
                hit.surface_id = node->primitive.mcnp_surface_id;
                double px, py, pz;
                alea_ray_point_at(ray, t[j], &px, &py, &pz);
                primitive_normal_at(prim->type, &prim->data, px, py, pz,
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

/**
 * Raycast all surfaces in a universe using a translated ray.
 * The ray is in element-local coordinates (origin shifted by -element_center).
 * Since translation preserves t, hits can be added directly to the result.
 */
static void raycast_universe_surfaces(const alea_system_t* sys,
                                      const alea_ray_t* local_ray,
                                      int universe_id,
                                      double t_min, double t_max,
                                      alea_raycast_result_t* result) {
    const alea_universe_t* univ = alea_get_universe(sys, universe_id);
    if (!univ) return;

    for (size_t c = 0; c < univ->cell_count; c++) {
        size_t cell_idx = univ->cell_indices[c];
        const alea_cell_entry_t* cell = &sys->cells.data[cell_idx];
        raycast_tree_primitives(sys, local_ray, cell->root_node_id,
                                t_min, t_max, result);
    }
}

/* ray_bbox_slab_enter_exit in ray_bbox.h replaces ray_bbox_slab_enter_exit */

/**
 * DDA raycast through a rectangular lattice.
 * Steps through elements along the ray, raycasting base universe surfaces
 * in each element's local coordinate system.
 */
static void raycast_lattice_rect(const alea_system_t* sys,
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
            alea_ray_hit_t bnd = { .t = t_cur, .surface_id = 0 };
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
static void raycast_lattice_hex(const alea_system_t* sys,
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
                alea_ray_hit_t bnd = { .t = t_cur, .surface_id = 0 };
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
static void raycast_add_lattice_hits(const alea_system_t* sys,
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

/* ============================================================================
 * CONVENIENCE FUNCTIONS
 * ============================================================================ */

int alea_raycast(const alea_system_t* sys,
                double ox, double oy, double oz,
                double dx, double dy, double dz,
                double t_max,
                alea_raycast_result_t* result) {
    if (!result) return -1;

    /* Free any prior allocations then reinitialize (safe for reuse) */
    alea_raycast_result_free(result);

    alea_ray_t ray;
    if (alea_ray_init(&ray, ox, oy, oz, dx, dy, dz) != 0) {
        return -1;  /* Zero-length direction */
    }

    /* t_max=0 means infinite */
    double effective_t_max = (t_max <= 0) ? DBL_MAX : t_max;
    int rc = alea_raycast_surfaces(sys, &ray, 0, effective_t_max, result);
    if (rc != 0) return rc;

    /* Add surface hits from lattice elements (DDA) */
    raycast_add_lattice_hits(sys, &ray, 0, effective_t_max, result);

    /* Re-sort and dedup after adding lattice hits */
    if (result->hit_count > 1) {
        sort_hits(result->hits, result->hit_count);
        size_t write = 1;
        for (size_t read = 1; read < result->hit_count; read++) {
            int same_t = fabs(result->hits[read].t - result->hits[write - 1].t) <= DEDUP_EPSILON;
            int same_surf = result->hits[read].surface_id == result->hits[write - 1].surface_id;
            if (!(same_t && same_surf)) {
                result->hits[write++] = result->hits[read];
            }
        }
        result->hit_count = write;
    }

    return alea_raycast_to_segments(sys, result);
}

int alea_ray_first_cell(const alea_system_t* sys,
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
    if (result.segment_count > 0) {
        /* Find first non-void segment */
        for (size_t i = 0; i < result.segment_count; i++) {
            if (result.segments[i].cell_id >= 0) {
                first_cell = result.segments[i].cell_id;
                if (out_t) *out_t = result.segments[i].t_enter;
                break;
            }
        }
    }

    alea_raycast_result_free(&result);
    return first_cell;
}

int alea_ray_is_occluded(const alea_system_t* sys,
                        double ox, double oy, double oz,
                        double dx, double dy, double dz,
                        double t_max) {
    if (!sys) return 0;

    /* Reuse a thread-local result to avoid malloc per shadow ray */
    static _Thread_local alea_raycast_result_t tls_result;
    static _Thread_local int tls_init = 0;
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
    sort_hits(tls_result.hits, tls_result.hit_count);

    /* Walk intervals and early-out on first non-void cell */
    double t_prev = 0;
    int prev_cell_idx = -1;

    for (size_t i = 0; i <= tls_result.hit_count; i++) {
        double t_curr = (i < tls_result.hit_count) ? tls_result.hits[i].t : effective_t_max;

        if (t_curr > t_prev + RAY_EPSILON) {
            int cell_id, cell_idx, material_id;
            double density;
            int crossed_surface = (i > 0) ? tls_result.hits[i - 1].surface_id : -1;
            find_cell_after_crossing(sys, &ray, t_prev, t_curr,
                                     prev_cell_idx, crossed_surface,
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
    for (size_t i = 0; i < result->segment_count; i++) {
        const alea_ray_segment_t* seg = &result->segments[i];

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
 * Test ray against a specific cell's surfaces only.
 * Returns closest hit distance, or DBL_MAX if no hit.
 */
static double raycast_cell_surfaces(const alea_system_t* sys,
                                    const alea_ray_t* ray,
                                    const alea_cell_entry_t* cell,
                                    double t_min, double t_max,
                                    int* out_surface_id) {
    double closest_t = DBL_MAX;
    *out_surface_id = -1;

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

        double t[4];
        int count = ray_intersect_primitive(ray, prim->type, &prim->data, t);

        for (int j = 0; j < count; j++) {
            if (t[j] > t_min && t[j] < t_max && t[j] < closest_t) {
                closest_t = t[j];
                *out_surface_id = surf->mcnp_surface_id;
            }
        }
    }

    return closest_t;
}

/**
 * Find the cell at a point, returning cell index or -1 for void.
 */
static int find_cell_at_point(const alea_system_t* sys,
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

/**
 * Context for finding closest intersection during BVH traversal.
 */
typedef struct {
    const alea_system_t* sys;
    const alea_ray_t* ray;
    double t_min;
    double closest_t;
    int closest_surface_id;
} closest_hit_ctx_t;

/**
 * BVH callback that finds the closest intersection.
 */
static void find_closest_callback(uint32_t surface_idx, void* userdata) {
    closest_hit_ctx_t* ctx = (closest_hit_ctx_t*)userdata;
    const alea_surface_entry_t* surf = &ctx->sys->surfaces.data[surface_idx];
    const alea_primitive_entry_t* prim = &ctx->sys->primitives.data[surf->primitive_id];

    double t[4];
    int count = ray_intersect_primitive(ctx->ray, prim->type, &prim->data, t);

    for (int j = 0; j < count; j++) {
        if (t[j] > ctx->t_min && t[j] < ctx->closest_t) {
            ctx->closest_t = t[j];
            ctx->closest_surface_id = surf->mcnp_surface_id;
        }
    }
}

/**
 * Find the closest surface intersection using BVH or linear scan.
 */
static double find_closest_intersection(const alea_system_t* sys,
                                        const alea_ray_t* ray,
                                        double t_min, double t_max,
                                        int* out_surface_id) {
    *out_surface_id = -1;

    /* Use BVH if available */
#if !BVH_DISABLED
    if (sys->surface_bvh) {
        closest_hit_ctx_t ctx = {
            .sys = sys,
            .ray = ray,
            .t_min = t_min,
            .closest_t = t_max,
            .closest_surface_id = -1
        };
        alea_bvh_traverse(sys->surface_bvh, ray, t_min, t_max,
                         find_closest_callback, &ctx);
        *out_surface_id = ctx.closest_surface_id;
        return ctx.closest_t;
    }
#endif

    /* Fallback: linear scan */
    double closest_t = t_max;
    for (size_t i = 0; i < alea_vec_count(&sys->surfaces); i++) {
        const alea_surface_entry_t* surf = &sys->surfaces.data[i];
        const alea_primitive_entry_t* prim = &sys->primitives.data[surf->primitive_id];

        double t[4];
        int count = ray_intersect_primitive(ray, prim->type, &prim->data, t);

        for (int j = 0; j < count; j++) {
            if (t[j] > t_min && t[j] < closest_t) {
                closest_t = t[j];
                *out_surface_id = surf->mcnp_surface_id;
            }
        }
    }
    return closest_t;
}

int alea_raycast_cell_aware(const alea_system_t* sys,
                           double ox, double oy, double oz,
                           double dx, double dy, double dz,
                           double t_max,
                           alea_raycast_result_t* result) {
    if (!sys || !result) return -1;

    /* Free any prior allocations then reinitialize (safe for reuse) */
    alea_raycast_result_free(result);

    alea_ray_t ray;
    if (alea_ray_init(&ray, ox, oy, oz, dx, dy, dz) != 0) {
        return -1;  /* Zero-length direction */
    }
    result->ray = ray;

    double effective_t_max = (t_max <= 0) ? DBL_MAX : t_max;

    if (alea_raycast_ensure_caches(sys)) {
        ALEA_LOG_WARN("Lazy-building raycast caches. Call alea_build_spatial_index() "
                     "before concurrent raycast calls to avoid data races.");
    }

    /* Current position along ray */
    double t_current = 0;
    int prev_cell_idx = -2;  /* -2 = no previous */
    int prev_surface_id = -1;  /* Surface ID we're about to cross */

    /* Safety limit to prevent infinite loops */
    int max_iterations = 10000;

    while (t_current < effective_t_max && max_iterations-- > 0) {
        int cell_idx = -1;
        int cell_id = -1;
        int material_id = 0;
        double density = 0;
        int found_via_neighbor = 0;

        /* Try neighbor lookup first if we have previous cell and surface info */
        if (prev_cell_idx >= 0 && prev_surface_id >= 0) {
            found_via_neighbor = raycast_find_neighbor(sys, prev_cell_idx,
                                                       prev_surface_id,
                                                       &cell_id, &cell_idx,
                                                       &material_id, &density);
        }

        /* Fall back to full lookup if neighbor lookup failed or not available */
        if (!found_via_neighbor) {
            double px, py, pz;
            alea_ray_point_at(&ray, t_current + RAY_EPSILON, &px, &py, &pz);
            cell_idx = find_cell_at_point(sys, px, py, pz, &material_id, &density);
            if (cell_idx >= 0 && (size_t)cell_idx < alea_vec_count(&sys->cells)) {
                cell_id = sys->cells.data[cell_idx].mcnp_cell_id;
            }
        }

        /* Find next surface crossing:
         * Use per-cell surface index when inside a cell (the whole point of cell-aware),
         * fall back to global search only for void regions. */
        int hit_surface_id = -1;
        double t_next;
        if (cell_idx >= 0 && (size_t)cell_idx < alea_vec_count(&sys->cells) &&
            sys->cells.data[cell_idx].surface_indices) {
            t_next = raycast_cell_surfaces(sys, &ray,
                                           &sys->cells.data[cell_idx],
                                           t_current + RAY_EPSILON, effective_t_max,
                                           &hit_surface_id);
        } else {
            t_next = find_closest_intersection(sys, &ray,
                                               t_current + RAY_EPSILON, effective_t_max,
                                               &hit_surface_id);
        }

        /* Ensure we make progress at any scale:
         * absolute 1e-10 dominates near origin; relative 1e-6 at large t */
        if (t_next <= t_current + RAY_EPSILON) {
            t_next = t_current * (1.0 + 1e-6) + RAY_EPSILON;
        }

        /* Add or extend segment */
        if (cell_idx == prev_cell_idx && result->segment_count > 0) {
            /* Extend previous segment */
            result->segments[result->segment_count - 1].t_exit = t_next;
        } else {
            /* Create new segment */
            alea_ray_segment_t seg;
            seg.t_enter = t_current;
            seg.t_exit = t_next;
            seg.cell_id = cell_id;
            seg.material_id = material_id;
            seg.density = density;
            seg.enter_hit_index = -1;  /* cell-aware path doesn't track hit indices */
            add_segment(result, &seg);
            prev_cell_idx = cell_idx;
        }

        /* Remember surface we're crossing for next iteration's neighbor lookup */
        prev_surface_id = hit_surface_id;

        /* Move past the intersection */
        t_current = t_next;

        /* If we hit nothing, we're done */
        if (t_next >= effective_t_max - RAY_EPSILON) {
            break;
        }
    }

    return 0;
}
