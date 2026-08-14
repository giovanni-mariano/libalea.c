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
#include <limits.h>
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
                                   alea_ray_first_visible_result_t* first_visible,
                                   int* first_cell_id,
                                   double* first_cell_t,
                                   double first_cell_t_min,
                                   int first_cell_material_filter,
                                   double visible_t_min,
                                   int visible_material_filter,
                                   bool visible_wants_normal,
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
    result->ancestor_surface_queries = 0;
    result->ancestor_path_entries_examined = 0;
    result->ancestor_max_path_depth = 0;
    result->ancestor_unattributed_surface_tests = 0;
    result->ancestor_unattributed_queries = 0;
    result->ancestor_unattributed_winning_events = 0;
    memset(result->ancestor_hot_cells, 0, sizeof(result->ancestor_hot_cells));
    result->lattice_entry_calls = 0;
    result->lattice_entry_candidates = 0;
    result->lattice_entry_dda_steps = 0;
    result->lattice_entry_no_entry_results = 0;
    result->lattice_entry_future_entry_results = 0;
    result->lattice_entry_already_inside_results = 0;
    result->lattice_entry_ancestor_surface_tests = 0;
    result->lattice_entry_ancestor_events = 0;
    result->lattice_entry_canonical_rejections = 0;
}

void alea_raycast_result_free(alea_raycast_result_t* result) {
    alea_vec_free(&result->hits);
    alea_vec_free(&result->segments);
    alea_vec_free(&result->paths);
    alea_vec_free(&result->path_entries);
    raycast_result_reset_counters(result);
}

void alea_raycast_result_clear(alea_raycast_result_t* result) {
    alea_vec_clear(&result->hits);
    alea_vec_clear(&result->segments);
    alea_vec_clear(&result->paths);
    alea_vec_clear(&result->path_entries);
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

static int reserve_batch_budget(atomic_uint_fast64_t* counter, uint64_t limit,
                                uint64_t amount) {
    if (!counter || limit == 0 || amount == 0) return 0;
    uint_fast64_t used = atomic_load_explicit(counter, memory_order_relaxed);
    for (;;) {
        if (used > limit || amount > limit - used) return -1;
        if (atomic_compare_exchange_strong(counter, &used, used + amount)) {
            return 0;
        }
    }
}

static int add_segment(alea_raycast_result_t* result, const alea_ray_segment_t* seg) {
    if (reserve_batch_budget(result->segment_counter, result->segment_limit, 1) != 0) {
        result->segment_limit_exceeded = 1;
        return -1;
    }
    if (result->path_entry_counter && result->path_entry_limit != 0 &&
        seg->path_index != UINT32_MAX && seg->path_index < result->paths.count) {
        uint64_t path_count = result->paths.data[seg->path_index].count;
        if (reserve_batch_budget(result->path_entry_counter,
                                 result->path_entry_limit, path_count) != 0) {
            result->path_entry_limit_exceeded = 1;
            return -1;
        }
    }
    int res = alea_vec_push(&result->segments, *seg, alea_ray_segment_t);
    return res != 0 ? -1 : 0;
}

static uint64_t path_hash_bytes(uint64_t hash, const void* data, size_t size) {
    const unsigned char* bytes = (const unsigned char*)data;
    for (size_t i = 0; i < size; ++i) {
        hash ^= bytes[i];
        hash *= UINT64_C(1099511628211);
    }
    return hash;
}

static uint64_t path_entry_occurrence_key(uint64_t parent_key,
                                           const alea_hier_ray_path_entry_t* entry) {
    uint64_t hash = parent_key ? parent_key : UINT64_C(1469598103934665603);
    hash = path_hash_bytes(hash, &entry->cell_index, sizeof(entry->cell_index));
    hash = path_hash_bytes(hash, &entry->cell_id, sizeof(entry->cell_id));
    hash = path_hash_bytes(hash, &entry->universe_id, sizeof(entry->universe_id));
    hash = path_hash_bytes(hash, &entry->fill_universe, sizeof(entry->fill_universe));
    hash = path_hash_bytes(hash, &entry->depth, sizeof(entry->depth));
    hash = path_hash_bytes(hash, &entry->is_lattice, sizeof(entry->is_lattice));
    hash = path_hash_bytes(hash, &entry->lat_fill_universe, sizeof(entry->lat_fill_universe));
    hash = path_hash_bytes(hash, &entry->lat_i, sizeof(entry->lat_i));
    hash = path_hash_bytes(hash, &entry->lat_j, sizeof(entry->lat_j));
    hash = path_hash_bytes(hash, &entry->lat_k, sizeof(entry->lat_k));
    hash = path_hash_bytes(hash, &entry->lat_ox, sizeof(entry->lat_ox));
    hash = path_hash_bytes(hash, &entry->lat_oy, sizeof(entry->lat_oy));
    hash = path_hash_bytes(hash, &entry->lat_oz, sizeof(entry->lat_oz));
    hash = path_hash_bytes(hash, entry->transform.m, sizeof(entry->transform.m));
    return hash;
}

static int ray_path_entry_equal(const alea_ray_path_entry_t* a,
                                const alea_ray_path_entry_t* b) {
    return a->cell_index == b->cell_index &&
           a->cell_id == b->cell_id &&
           a->material_id == b->material_id &&
           a->universe_id == b->universe_id &&
           a->fill_universe == b->fill_universe &&
           a->depth == b->depth &&
           a->is_lattice == b->is_lattice &&
           a->lattice_origin[0] == b->lattice_origin[0] &&
           a->lattice_origin[1] == b->lattice_origin[1] &&
           a->lattice_origin[2] == b->lattice_origin[2] &&
           a->occurrence_key == b->occurrence_key;
}

static int capture_hier_path(alea_raycast_result_t* result,
                             const alea_hier_ray_path_t* path,
                             uint32_t* out_path_index) {
    *out_path_index = UINT32_MAX;
    if (!result->capture_paths || !path || path->count <= 0) return 0;

    size_t count = (size_t)path->count;
    alea_ray_path_entry_t entries[ALEA_HIER_RAY_PATH_MAX];
    uint64_t parent_key = 0;
    for (size_t i = 0; i < count; ++i) {
        const alea_hier_ray_path_entry_t* src = &path->entries[i];
        alea_ray_path_entry_t* dst = &entries[i];
        dst->cell_index = src->cell_index;
        dst->cell_id = src->cell_id;
        dst->material_id = src->material_id;
        dst->universe_id = src->universe_id;
        dst->fill_universe = src->fill_universe;
        dst->depth = src->depth;
        dst->is_lattice = src->is_lattice;
        dst->lattice_origin[0] = src->lat_ox;
        dst->lattice_origin[1] = src->lat_oy;
        dst->lattice_origin[2] = src->lat_oz;
        parent_key = path_entry_occurrence_key(parent_key, src);
        dst->occurrence_key = parent_key;
    }

    for (size_t i = 0; i < result->paths.count; ++i) {
        const alea_ray_path_t* candidate = &result->paths.data[i];
        if (candidate->count != count) continue;
        int same = 1;
        for (size_t j = 0; j < count; ++j) {
            if (!ray_path_entry_equal(&result->path_entries.data[candidate->offset + j],
                                      &entries[j])) {
                same = 0;
                break;
            }
        }
        if (same) {
            *out_path_index = (uint32_t)i;
            return 0;
        }
    }

    alea_ray_path_t record = {
        .offset = (uint32_t)result->path_entries.count,
        .count = (uint16_t)count,
    };
    for (size_t i = 0; i < count; ++i) {
        if (alea_vec_push(&result->path_entries, entries[i], alea_ray_path_entry_t) != 0)
            return -1;
    }
    if (alea_vec_push(&result->paths, record, alea_ray_path_t) != 0)
        return -1;
    *out_path_index = (uint32_t)(result->paths.count - 1);
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
                                 double px, double py, double pz,
                                 int* out_cell_id, int* out_cell_idx,
                                 int* out_material_id, double* out_density) {
    if (!sys->cell_adjacency_built || current_cell_idx < 0) {
        return 0;
    }

    /* A surface may bound several cells, and crossing a surface referenced
     * by a non-convex cell does not always leave that cell, so the first
     * adjacency entry for surface_id is only a candidate: accept a neighbor
     * only if it contains the sample point past the crossing. Callers gate
     * this on universe-0 cells, whose local frame is the world frame. */
    const alea_cell_entry_t* cell = &sys->cells.data[current_cell_idx];
    for (size_t i = 0; i < cell->neighbor_count; i++) {
        if (cell->neighbors[i].surface_id != surface_id) continue;
        uint32_t nb_idx = cell->neighbors[i].neighbor_index;
        if (nb_idx >= alea_vec_count(&sys->cells)) continue;
        const alea_cell_entry_t* neighbor = &sys->cells.data[nb_idx];
        if (neighbor->root_node_id == ALEA_NODE_ID_INVALID) continue;
        if (!alea_contains_point(sys, neighbor->root_node_id, px, py, pz))
            continue;
        *out_cell_id = neighbor->mc_cell_id;
        *out_cell_idx = (int)nb_idx;
        *out_material_id = neighbor->material_id;
        *out_density = neighbor->density;
        return 1;
    }
    return 0;  /* No containing neighbor (re-entrant crossing/exterior) */
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

/* Hits from the all-cells queries are in DFS preorder. Follow the first
 * strictly-deepening chain and take its deepest hit, so overlapping
 * same-depth cells resolve to the first containing cell in deck order —
 * matching the canonical resolver and the hier descent. */
static int deepest_first_chain_index(const alea_cell_hit_t* hits, int n) {
    int idx = 0;
    while (idx + 1 < n && hits[idx + 1].depth == hits[idx].depth + 1) idx++;
    return idx;
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
        const alea_cell_hit_t* hit =
            &hits[deepest_first_chain_index(hits, num_hits)];
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
        const alea_cell_hit_t* hit =
            &hits[deepest_first_chain_index(hits, num_hits)];
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
            double t_sample = t_prev +
                fmin(0.5 * (t_curr - t_prev), SURFACE_SAMPLE_OFFSET);
            double px, py, pz;
            alea_ray_point_at(ray, t_sample, &px, &py, &pz);
            found = raycast_find_neighbor(sys, prev_cell_idx, crossed_surface_id,
                                          px, py, pz,
                                          out_cell_id, out_cell_idx,
                                          out_material_id, out_density);
            /* In hier mode a filled container is not a terminal answer:
             * reject it so the full lookup descends into the fill. */
            if (found && use_hier_lookup) {
                const alea_cell_entry_t* nb = &sys->cells.data[*out_cell_idx];
                if (nb->fill_universe > 0 ||
                    (nb->lat_type != 0 && nb->lat_fill)) {
                    found = 0;
                    *out_cell_id = -1;
                    *out_cell_idx = -1;
                    *out_material_id = 0;
                    *out_density = 0;
                }
            }
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
                seg.path_index = UINT32_MAX;
                seg.resolution_flags =
                    (cell_idx >= 0 &&
                     (size_t)cell_idx < alea_vec_count(&sys->cells) &&
                     alea_cell_entry_is_container(&sys->cells.data[cell_idx]))
                        ? ALEA_RESOLVE_UNDEFINED_FILL : 0;

                if (add_segment(result, &seg) != 0) return -1;
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
static double lattice_next_boundary(const alea_ray_t* ray,
                                    const alea_cell_entry_t* lat_cell,
                                    double t_min,
                                    double t_max);

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

static int lattice_rect_bounds(const alea_cell_entry_t* cell,
                               alea_bbox_t* out_bounds) {
    if (!cell || !out_bounds || cell->lat_pitch[0] <= 0.0 ||
        cell->lat_pitch[1] <= 0.0 || cell->lat_pitch[2] <= 0.0) {
        return -1;
    }
    const int ni = cell->lat_fill_dims[1] - cell->lat_fill_dims[0] + 1;
    const int nj = cell->lat_fill_dims[3] - cell->lat_fill_dims[2] + 1;
    const int nk = cell->lat_fill_dims[5] - cell->lat_fill_dims[4] + 1;
    if (ni <= 0 || nj <= 0 || nk <= 0) return -1;
    *out_bounds = (alea_bbox_t){
        .min_x = cell->lat_lower_left[0],
        .max_x = cell->lat_lower_left[0] + ni * cell->lat_pitch[0],
        .min_y = cell->lat_lower_left[1],
        .max_y = cell->lat_lower_left[1] + nj * cell->lat_pitch[1],
        .min_z = cell->lat_lower_left[2],
        .max_z = cell->lat_lower_left[2] + nk * cell->lat_pitch[2]
    };
    return 0;
}

static int lattice_hex_bounds(const alea_cell_entry_t* cell,
                              alea_bbox_t* out_bounds) {
    if (!cell || !out_bounds || cell->lat_pitch[0] <= 0.0) return -1;
    const int imin = cell->lat_fill_dims[0], imax = cell->lat_fill_dims[1];
    const int jmin = cell->lat_fill_dims[2], jmax = cell->lat_fill_dims[3];
    const int nk = cell->lat_fill_dims[5] - cell->lat_fill_dims[4] + 1;
    if (imax < imin || jmax < jmin || nk <= 0) return -1;

    double min_x = DBL_MAX, max_x = -DBL_MAX;
    double min_y = DBL_MAX, max_y = -DBL_MAX;
    const double pitch = cell->lat_pitch[0];
    const double radius = pitch / M_SQRT3;
    for (int i = imin; i <= imax; i++) {
        for (int j = jmin; j <= jmax; j++) {
            const double x = i * pitch + j * pitch * 0.5;
            const double y = j * pitch * M_SQRT3 * 0.5;
            if (x - radius < min_x) min_x = x - radius;
            if (x + radius > max_x) max_x = x + radius;
            if (y - pitch * 0.5 < min_y) min_y = y - pitch * 0.5;
            if (y + pitch * 0.5 > max_y) max_y = y + pitch * 0.5;
        }
    }
    *out_bounds = (alea_bbox_t){
        .min_x = min_x, .max_x = max_x,
        .min_y = min_y, .max_y = max_y,
        .min_z = nk > 1 ? cell->lat_lower_left[2] : -1e30,
        .max_z = nk > 1 ? cell->lat_lower_left[2] + nk * cell->lat_pitch[2]
                         : 1e30
    };
    return 0;
}

/* The geometry-specific DDA walkers only decide the next interval and its
 * canonical lattice location.  Publishing the synthetic transition and
 * tracing the selected element are shared so full-trace and future early-stop
 * steppers cannot acquire different element-local ray conventions. */
static int raycast_lattice_element_step(
    alea_system_t* sys,
    const alea_ray_t* ray,
    const alea_lattice_location_t* location,
    double t_enter,
    double t_exit,
    bool emit_synthetic_transition,
    alea_raycast_result_t* result) {
    if (!sys || !ray || !result) return -1;
    if (emit_synthetic_transition) {
        alea_ray_hit_t boundary = {
            .t = t_enter,
            .surface_id = 0,
            .primitive_id = ALEA_PRIMITIVE_ID_INVALID
        };
        if (add_hit(result, &boundary) != 0) {
            ALEA_LOG_WARN("add_hit failed (out of memory) - lattice raycast incomplete");
            return -1;
        }
    }
    if (!location || t_exit <= t_enter) return 0;

    alea_ray_t local_ray = *ray;
    local_ray.ox -= location->ox;
    local_ray.oy -= location->oy;
    local_ray.oz -= location->oz;
    raycast_universe_surfaces(sys, &local_ray, location->fill_universe,
                              t_enter, t_exit, result);
    return 0;
}

static bool raycast_has_hit_at(const alea_raycast_result_t* result, double t) {
    if (!result) return false;
    for (size_t i = 0; i < result->hits.count; i++) {
        /* Parsed lattice pitch/bounds can differ slightly from the source CSG
         * plane used to infer them.  Treat a nearby physical hit as the same
         * ownership boundary so the DDA does not publish a micro-segment. */
        if (result->hits.data[i].surface_id > 0 &&
            fabs(result->hits.data[i].t - t) <= SURFACE_SAMPLE_OFFSET)
            return true;
    }
    return false;
}

static int lattice_raycast_interval(const alea_ray_t* ray,
                                    const alea_cell_entry_t* cell,
                                    double t_min, double t_max,
                                    double* out_enter, double* out_exit) {
    if (!ray || !cell || !out_enter || !out_exit) return -1;
    *out_enter = t_min;
    *out_exit = t_max;

    alea_bbox_t bounds;
    if (cell->lat_type == 1) {
        if (lattice_rect_bounds(cell, &bounds) != 0) return -1;
        if (cell->lat_fill_repeating) return 0;
    } else if (cell->lat_type == 2) {
        if (lattice_hex_bounds(cell, &bounds) != 0) return -1;
    } else {
        return -1;
    }
    return ray_bbox_slab_enter_exit(ray, &bounds, t_min, t_max,
                                    out_enter, out_exit) ? 0 : 1;
}

static int lattice_raycast_step_limit(const alea_ray_t* ray,
                                      const alea_cell_entry_t* cell,
                                      double t_enter, double t_exit) {
    const int ni = cell->lat_fill_dims[1] - cell->lat_fill_dims[0] + 1;
    const int nj = cell->lat_fill_dims[3] - cell->lat_fill_dims[2] + 1;
    const int nk = cell->lat_fill_dims[5] - cell->lat_fill_dims[4] + 1;
    if (ni <= 0 || nj <= 0 || nk <= 0) return 0;
    if (cell->lat_type == 1 && cell->lat_fill_repeating) {
        const double span = t_exit - t_enter;
        const double crossings = fabs(ray->dx) * span / cell->lat_pitch[0]
            + fabs(ray->dy) * span / cell->lat_pitch[1]
            + fabs(ray->dz) * span / cell->lat_pitch[2];
        return crossings < (double)INT_MAX - 8.0
            ? (int)ceil(crossings) + 8 : INT_MAX;
    }
    return cell->lat_type == 1 ? 2 * (ni + nj + nk) + 4
                               : 3 * (ni + nj + nk) + 10;
}

/* Full global tracing uses the same boundary query as the cell-aware
 * first-visible walker.  Geometry-specific rectangular/hex math is confined
 * to lattice_next_boundary(); this loop owns interval construction, canonical
 * location lookup, synthetic event publication, and element-local tracing. */
static void raycast_lattice_walk(alea_system_t* sys,
                                 const alea_ray_t* ray,
                                 const alea_cell_entry_t* cell,
                                 double t_min, double t_max,
                                 alea_raycast_result_t* result) {
    double t_enter, t_exit;
    const int interval = lattice_raycast_interval(ray, cell, t_min, t_max,
                                                  &t_enter, &t_exit);
    if (interval != 0 || t_exit <= t_enter) return;

    const int max_steps = lattice_raycast_step_limit(ray, cell, t_enter, t_exit);
    alea_lattice_location_t previous;
    int have_previous = 0;
    double t_current = t_enter;
    if (t_enter > t_min + RAY_EPSILON &&
        !raycast_has_hit_at(result, t_enter)) {
        if (raycast_lattice_element_step(sys, ray, NULL, t_enter, t_enter,
                                         true, result) != 0) {
            return;
        }
    }
    for (int step = 0; step < max_steps && t_current < t_exit; step++) {
        double t_next = lattice_next_boundary(ray, cell, t_current, t_exit);
        if (t_next <= t_current + RAY_EPSILON || t_next > t_exit)
            t_next = t_exit;
        const double t_sample = t_current +
            0.5 * (t_next - t_current);
        double x, y, z;
        alea_ray_point_at(ray, t_sample, &x, &y, &z);
        alea_lattice_location_t location;
        const int has_location = alea_lattice_locate_point(
            sys, cell, x, y, z, &location) == 1;

        const bool changed = !have_previous || !has_location ||
            location.i != previous.i || location.j != previous.j ||
            location.k != previous.k ||
            location.fill_universe != previous.fill_universe;
        const bool emit_synthetic = step > 0 &&
            t_current > t_enter + RAY_EPSILON && changed;
        if (raycast_lattice_element_step(
                sys, ray, has_location ? &location : NULL,
                t_current, t_next, emit_synthetic, result) != 0) {
            return;
        }
        if (has_location) {
            previous = location;
            have_previous = 1;
        } else {
            have_previous = 0;
        }
        t_current = t_next;
    }
    if (t_exit < t_max - RAY_EPSILON &&
        !raycast_has_hit_at(result, t_exit)) {
        (void)raycast_lattice_element_step(sys, ray, NULL, t_exit, t_exit,
                                           true, result);
    }
}

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
    raycast_lattice_walk(sys, ray, lat_cell, t_min, t_max, result);
}

static void raycast_lattice_hex(alea_system_t* sys,
                                const alea_ray_t* ray,
                                const alea_cell_entry_t* lat_cell,
                                double t_min, double t_max,
                                alea_raycast_result_t* result) {
    raycast_lattice_walk(sys, ray, lat_cell, t_min, t_max, result);
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

int alea_raycast_global_reuse_nocache(alea_system_t* sys,
                                      const alea_ray_t* ray,
                                      double t_max,
                                      alea_raycast_result_t* result) {
    if (!sys || !ray || !result) return -1;

    /* This is the reusable counterpart to alea_raycast(). It intentionally
     * assumes the caller prepared the raycast caches, as render loops do. */
    alea_raycast_result_clear(result);
    double effective_t_max = (t_max <= 0) ? DBL_MAX : t_max;
    return raycast_global_pipeline(sys, ray, 0, effective_t_max,
                                   system_has_lattice_cells(sys), false,
                                   result);
}

int alea_raycast_global_breakpoints_reuse_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    alea_raycast_result_t* result) {
    if (!sys || !ray || !result) return -1;

    /* Keep this deliberately below raycast_global_pipeline(): that helper
     * materializes precedence-selected segments as a compatibility service.
     * Coverage diagnostics need only the complete geometric partition. */
    alea_raycast_result_clear(result);
    const double effective_t_max = t_max <= 0 ? DBL_MAX : t_max;
    if (raycast_surfaces_impl(sys, ray, 0, effective_t_max, result) != 0)
        return -1;
    if (system_has_fill_cells(sys))
        raycast_add_fill_hits(sys, ray, 0, effective_t_max, result);
    if (system_has_lattice_cells(sys))
        raycast_add_lattice_hits(sys, ray, 0, effective_t_max, result);
    dedup_sorted_hits(result);
    return 0;
}

void alea_ray_boundary_event_result_init(alea_ray_boundary_event_result_t* result) {
    if (!result) return;
    memset(result, 0, sizeof(*result));
}

void alea_ray_boundary_event_result_clear(alea_ray_boundary_event_result_t* result) {
    if (!result) return;
    alea_vec_clear(&result->events);
}

void alea_ray_boundary_event_result_free(alea_ray_boundary_event_result_t* result) {
    if (!result) return;
    alea_vec_free(&result->events);
}

static const alea_ray_hit_t* boundary_event_find_hit(
    const alea_raycast_result_t* trace, double t, int surface_id) {
    for (size_t i = 0; i < trace->hits.count; i++) {
        const alea_ray_hit_t* hit = &trace->hits.data[i];
        if (hit->surface_id == surface_id &&
            fabs(hit->t - t) <= RAY_EPSILON)
            return hit;
    }
    return NULL;
}

static int boundary_event_append(alea_ray_boundary_event_result_t* events,
                                 const alea_ray_segment_t* before,
                                 const alea_ray_segment_t* after,
                                 double t,
                                 alea_ray_boundary_event_kind_t kind,
                                 int surface_id,
                                 const alea_ray_hit_t* hit,
                                 const alea_ray_boundary_event_options_internal_t* options) {
    if (options &&
        ((options->max_events && events->events.count >= options->max_events) ||
         (options->max_output_bytes &&
          events->events.count >=
              options->max_output_bytes / sizeof(*events->events.data)))) {
        alea_set_error_detail(ALEA_ERR_OVERFLOW,
                              "boundary-event query output limit exceeded");
        return -1;
    }
    alea_ray_boundary_event_t event = {
        .t = t,
        .kind = kind,
        .surface_id = surface_id,
        .primitive_id = hit ? hit->primitive_id : UINT32_MAX,
        .cell_before = before->cell_id,
        .cell_after = after->cell_id,
        .material_before = before->material_id,
        .material_after = after->material_id,
        .resolution_flags = before->resolution_flags | after->resolution_flags,
        .nx = hit ? hit->nx : 0,
        .ny = hit ? hit->ny : 0,
        .nz = hit ? hit->nz : 0
    };
    return alea_vec_push(&events->events, event, alea_ray_boundary_event_t);
}

static int boundary_events_append_group(
    const alea_raycast_result_t* trace,
    const alea_ray_segment_t* before, const alea_ray_segment_t* after,
    bool include_all_coincident_physical,
    const alea_ray_boundary_event_options_internal_t* options,
    alea_ray_boundary_event_result_t* events) {
    const double t = before->t_exit;
    int ownership_changed = before->cell_id != after->cell_id ||
                            before->material_id != after->material_id;
    int lowest_physical_id = -1;
    int* physical_ids = NULL;
    size_t physical_count = 0;
    size_t physical_capacity = 0;
    int saw_synthetic = 0;

    for (size_t i = 0; i < trace->hits.count; i++) {
        const alea_ray_hit_t* hit = &trace->hits.data[i];
        if (fabs(hit->t - t) > RAY_EPSILON) continue;
        if (hit->surface_id == 0) {
            saw_synthetic = 1;
            continue;
        }
        if (hit->surface_id <= 0) continue;
        if (!include_all_coincident_physical) {
            if (lowest_physical_id < 0 || hit->surface_id < lowest_physical_id)
                lowest_physical_id = hit->surface_id;
            continue;
        }
        size_t j;
        for (j = 0; j < physical_count; j++)
            if (physical_ids[j] == hit->surface_id) break;
        if (j == physical_count) {
            if (physical_count == physical_capacity) {
                size_t next_capacity = physical_capacity ? physical_capacity * 2 : 16;
                if (next_capacity < physical_capacity ||
                    next_capacity > (size_t)-1 / sizeof(*physical_ids)) {
                    alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                                          "too many coincident boundary surfaces");
                    free(physical_ids);
                    return -1;
                }
                int* next_ids = realloc(physical_ids,
                                        next_capacity * sizeof(*physical_ids));
                if (!next_ids) {
                    alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,
                                          "failed to grow coincident boundary surface list");
                    free(physical_ids);
                    return -1;
                }
                physical_ids = next_ids;
                physical_capacity = next_capacity;
            }
            physical_ids[physical_count++] = hit->surface_id;
        }
    }

    /* A segment can carry a synthetic boundary generated by traversal even
     * when its t differs microscopically from the surface hit group. */
    if (before->exit_surface_id == 0 || after->enter_surface_id == 0)
        saw_synthetic = 1;

    if (ownership_changed &&
        (include_all_coincident_physical ? physical_count != 0 : lowest_physical_id >= 0)) {
        /* Surface ID is the deterministic reporting tie-breaker: the lowest
         * positive MC surface ID wins in the normal path.  Diagnostic mode
         * preserves every physical participant at this crossing. */
        if (include_all_coincident_physical) {
            for (size_t i = 0; i < physical_count; i++) {
                for (size_t j = i + 1; j < physical_count; j++) {
                    if (physical_ids[j] < physical_ids[i]) {
                        int tmp = physical_ids[i];
                        physical_ids[i] = physical_ids[j];
                        physical_ids[j] = tmp;
                    }
                }
            }
        }
        size_t emit_count = include_all_coincident_physical ? physical_count : 1;
        for (size_t i = 0; i < emit_count; i++) {
            int surface_id = include_all_coincident_physical ? physical_ids[i] :
                                                            lowest_physical_id;
            const alea_ray_hit_t* hit = boundary_event_find_hit(trace, t,
                                                                 surface_id);
            if (boundary_event_append(events, before, after, t,
                                      ALEA_RAY_BOUNDARY_EVENT_PHYSICAL,
                                      surface_id, hit, options) != 0) {
                free(physical_ids);
                return -1;
            }
        }
    } else if (ownership_changed && !saw_synthetic) {
        if (boundary_event_append(events, before, after, t,
                                  ALEA_RAY_BOUNDARY_EVENT_UNRESOLVED,
                                  -1, NULL, options) != 0) {
            free(physical_ids);
            return -1;
        }
    }

    /* Physical events are emitted first; the synthetic event follows at the
     * same distance. This matches hit ordering and remains deterministic. */
    if (saw_synthetic &&
        boundary_event_append(events, before, after, t,
                              ALEA_RAY_BOUNDARY_EVENT_SYNTHETIC_LATTICE,
                              0, NULL, options) != 0) {
        free(physical_ids);
        return -1;
    }
    free(physical_ids);
    return 0;
}

static int boundary_events_from_trace(
    const alea_raycast_result_t* trace,
    const alea_ray_boundary_event_options_internal_t* options,
    alea_ray_boundary_event_result_t* events);

int alea_raycast_boundary_events_with_options(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    const alea_ray_boundary_event_options_internal_t* options,
    alea_raycast_result_t* trace, alea_ray_boundary_event_result_t* events) {
    if (!trace || !events)
        return -1;
    alea_ray_boundary_event_result_clear(events);
    /* Boundary events need the same ownership resolver as public first-hit
     * queries.  The legacy global path eventually falls back to a recursive
     * point lookup, which is both slower and cannot resolve filled/lattice
     * models reliably once the hierarchical cache is the active backend. */
    alea_raycast_result_clear(trace);
    double effective_t_max = (t_max <= 0) ? DBL_MAX : t_max;
    if (raycast_global_pipeline(sys, ray, 0, effective_t_max,
                                system_has_lattice_cells(sys), true,
                                trace) != 0) {
        alea_set_error_detail(ALEA_ERR_INVALID_STATE,
                              "boundary-event trace failed while resolving ray intervals");
        return -1;
    }
    return boundary_events_from_trace(trace, options, events);
}

/* Derive events without retracing.  This is the bridge used by the query
 * dispatcher while batch/event materialization remains a later phase. */
static int boundary_events_from_trace(
    const alea_raycast_result_t* trace,
    const alea_ray_boundary_event_options_internal_t* options,
    alea_ray_boundary_event_result_t* events) {
    if (!trace || !events)
        return -1;
    alea_ray_boundary_event_result_clear(events);

    for (size_t i = 0; i + 1 < trace->segments.count; i++) {
        const alea_ray_segment_t* before = &trace->segments.data[i];
        const alea_ray_segment_t* after = &trace->segments.data[i + 1];
        if (fabs(before->t_exit - after->t_enter) > RAY_EPSILON)
            continue;  /* Incomplete trace: do not invent an event. */

        if (boundary_events_append_group(
                trace, before, after,
                options && options->include_all_coincident_physical,
                options,
                events) != 0) {
            alea_ray_boundary_event_result_clear(events);
            return -1;
        }
    }
    return 0;
}

int alea_raycast_boundary_events_reuse_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    alea_raycast_result_t* trace, alea_ray_boundary_event_result_t* events) {
    return alea_raycast_boundary_events_with_options(sys, ray, t_max, NULL,
                                                      trace, events);
}

static void ray_query_output_clear(alea_ray_query_output_t* output) {
    if (!output) return;
    memset(output, 0, sizeof(*output));
    output->first_cell_id = -1;
    output->first_visible.cell_id = -1;
    output->first_visible.material_id = 0;
    output->first_visible.surface_id = -1;
    output->first_visible.primitive_id = UINT32_MAX;
}

static bool ray_plan_accepts_material(const alea_ray_plan_t* plan,
                                       const alea_ray_segment_t* seg,
                                       bool require_nonvoid) {
    if (seg->cell_id < 0) return false;
    if (require_nonvoid && seg->material_id == 0) return false;
    return plan->material_filter < 0 ||
           seg->material_id == plan->material_filter;
}

static const alea_ray_hit_t* ray_query_enter_hit(
    const alea_raycast_result_t* trace, const alea_ray_segment_t* seg) {
    if (seg->enter_hit_index >= 0 &&
        (size_t)seg->enter_hit_index < trace->hits.count)
        return &trace->hits.data[seg->enter_hit_index];
    return boundary_event_find_hit(trace, seg->t_enter, seg->enter_surface_id);
}

/* Keep scalar SEGMENTS queries consistent with the compact executor: an
 * interval that begins before the requested range becomes a cross-section,
 * not a surface entry.  This mutates only the reusable logical contents and
 * therefore retains all caller-owned capacity. */
static void ray_query_clip_segments(alea_raycast_result_t* trace,
                                    double t_min) {
    if (t_min <= 0.0) return;
    size_t write = 0;
    for (size_t i = 0; i < trace->segments.count; i++) {
        alea_ray_segment_t segment = trace->segments.data[i];
        if (segment.t_exit <= t_min + RAY_EPSILON) continue;
        if (segment.t_enter < t_min) {
            segment.t_enter = t_min;
            segment.enter_surface_id = -1;
            segment.enter_hit_index = -1;
        }
        trace->segments.data[write++] = segment;
    }
    trace->segments.count = write;
}

static int ray_query_fast_forward_trace(alea_system_t* sys,
                                        const alea_ray_t* ray,
                                        double t_max, bool emit_hits,
                                        alea_raycast_result_t* trace) {
    alea_raycast_result_clear(trace);
    trace->ray = *ray;
    return raycast_cell_aware_impl(sys, ray, t_max, true, emit_hits,
                                   NULL, NULL, NULL, 0, -1, 0, -1, false,
                                   trace);
}

static int ray_query_fast_reverse_trace(alea_system_t* sys,
                                        const alea_ray_t* ray,
                                        double t_max, bool emit_hits,
                                        alea_raycast_result_t* trace) {
    if (!isfinite(t_max) || t_max <= 0 || t_max == DBL_MAX) {
        alea_set_error_detail(ALEA_ERR_INVALID_ARG,
                              "FAST_REVERSE requires a finite positive t_max");
        return -1;
    }
    double ex, ey, ez;
    alea_ray_point_at(ray, t_max, &ex, &ey, &ez);
    alea_ray_t reverse_ray;
    alea_ray_init_normalized(&reverse_ray, ex, ey, ez,
                             -ray->dx, -ray->dy, -ray->dz);

    alea_raycast_result_t reverse;
    alea_raycast_result_init(&reverse);
    int rc = ray_query_fast_forward_trace(sys, &reverse_ray, t_max,
                                           emit_hits, &reverse);
    if (rc != 0) {
        alea_raycast_result_free(&reverse);
        return -1;
    }
    alea_raycast_result_clear(trace);
    trace->ray = *ray;
    if (alea_vec_reserve(&trace->segments, reverse.segments.count,
                         alea_ray_segment_t) != 0) {
        alea_raycast_result_free(&reverse);
        return -1;
    }
    for (size_t i = 0; i < reverse.segments.count; i++) {
        const alea_ray_segment_t* src =
            &reverse.segments.data[reverse.segments.count - 1 - i];
        alea_ray_segment_t dst = *src;
        dst.t_enter = t_max - src->t_exit;
        dst.t_exit = t_max - src->t_enter;
        dst.enter_surface_id = src->exit_surface_id;
        dst.exit_surface_id = src->enter_surface_id;
        dst.enter_hit_index = -1;
        trace->segments.data[trace->segments.count++] = dst;
    }
    /* Reverse-hit normalization is intentionally deferred.  Fast boundary
     * provenance requests retain the verified forward trace in bidirectional
     * mode, and FAST_REVERSE exposes segment ownership only for now. */
    alea_raycast_result_free(&reverse);
    return 0;
}

static int ray_query_segments_match(const alea_raycast_result_t* forward,
                                    const alea_raycast_result_t* reverse) {
    if (forward->segments.count != reverse->segments.count)
        return 0;
    for (size_t i = 0; i < forward->segments.count; i++) {
        const alea_ray_segment_t* a = &forward->segments.data[i];
        const alea_ray_segment_t* b = &reverse->segments.data[i];
        if (fabs(a->t_enter - b->t_enter) > RAY_EPSILON ||
            fabs(a->t_exit - b->t_exit) > RAY_EPSILON ||
            a->cell_id != b->cell_id || a->material_id != b->material_id)
            return 0;
    }
    return 1;
}

int alea_raycast_query_reuse_nocache(
    alea_system_t* sys, const alea_ray_t* ray,
    const alea_ray_query_t* query, alea_raycast_result_t* trace,
    alea_ray_boundary_event_result_t* events,
    alea_ray_query_output_t* output) {
    /* Clear every valid caller-owned output before validation too: malformed
     * descriptors must not leave a previous ray's answer observable. */
    ray_query_output_clear(output);
    if (events) alea_ray_boundary_event_result_clear(events);
    if (trace) alea_raycast_result_clear(trace);
    if (!sys || !ray || !query || !trace) return -1;
    alea_ray_plan_t plan;
    if (alea_ray_query_lower(query, &plan) != 0)
        return -1;
    if (plan.product == ALEA_RAY_QUERY_BOUNDARY_EVENTS && !events)
        return -1;

    const double t_max = plan.t_max;
    if (plan.backend == ALEA_RAY_QUERY_BACKEND_AUTO &&
        plan.product == ALEA_RAY_QUERY_FIRST_CELL && output) {
        if (alea_raycast_hier_first_cell_nocache(
                sys, ray, plan.t_min, t_max, plan.material_filter,
                trace, &output->first_cell_id,
                &output->first_cell_t) != 0)
            goto fail;
        return 0;
    }
    if (plan.backend == ALEA_RAY_QUERY_BACKEND_AUTO &&
        plan.product == ALEA_RAY_QUERY_FIRST_VISIBLE) {
        /* The hierarchical stepper validates each interval before returning,
         * then stops before allocating its segment/hit representation. */
        if (!output || alea_raycast_hier_first_visible_nocache(
                            sys, ray, plan.t_min, t_max,
                            plan.material_filter,
                            plan.requirements.need_normal != 0,
                            trace,
                            &output->first_visible) != 0)
            goto fail;
        return 0;
    }
    if (plan.backend == ALEA_RAY_QUERY_BACKEND_AUTO &&
        plan.product == ALEA_RAY_QUERY_ANY_HIT) {
        int any_hit = 0;
        if (alea_raycast_hier_any_hit_nocache(
                sys, ray, plan.t_min, t_max, plan.material_filter,
                trace, &any_hit) != 0)
            goto fail;
        if (output) output->any_hit = any_hit != 0;
        return 0;
    }
    if (plan.backend == ALEA_RAY_QUERY_BACKEND_AUTO &&
        plan.product == ALEA_RAY_QUERY_SEGMENTS) {
        /* Selected intervals always use the coherent hierarchical walker.
         * Lattices used to fall through to the global tracer here, which
         * made AUTO depend on model topology rather than query semantics. */
        if (ray_query_fast_forward_trace(sys, ray, t_max, false, trace) != 0)
            goto fail;
    } else if (plan.backend == ALEA_RAY_QUERY_BACKEND_FAST_FORWARD ||
        plan.backend == ALEA_RAY_QUERY_BACKEND_FAST_REVERSE ||
        plan.backend == ALEA_RAY_QUERY_BACKEND_FAST_FORWARD_REVERSE) {
        const bool emit_hits = plan.product == ALEA_RAY_QUERY_BOUNDARY_EVENTS ||
            (plan.product == ALEA_RAY_QUERY_FIRST_VISIBLE &&
             plan.requirements.need_surface_identity);
        if (plan.backend == ALEA_RAY_QUERY_BACKEND_FAST_FORWARD) {
            if (ray_query_fast_forward_trace(sys, ray, t_max, emit_hits, trace) != 0)
                goto fail;
        } else if (plan.backend == ALEA_RAY_QUERY_BACKEND_FAST_REVERSE) {
            if (ray_query_fast_reverse_trace(sys, ray, t_max, emit_hits, trace) != 0)
                goto fail;
        } else {
            alea_raycast_result_t reverse;
            alea_raycast_result_init(&reverse);
            int rc = ray_query_fast_forward_trace(sys, ray, t_max, emit_hits, trace);
            if (rc == 0)
                rc = ray_query_fast_reverse_trace(sys, ray, t_max, false, &reverse);
            if (rc == 0 && !ray_query_segments_match(trace, &reverse) && output)
                output->directional_mismatch = true;
            alea_raycast_result_free(&reverse);
            if (rc != 0) goto fail;
        }
    } else if (alea_raycast_global_reuse_nocache(sys, ray, t_max, trace) != 0) {
        goto fail;
    }

    if (plan.product == ALEA_RAY_QUERY_BOUNDARY_EVENTS) {
        const alea_ray_boundary_event_options_internal_t event_options = {
            .max_events = plan.max_events,
            .max_output_bytes = plan.max_output_bytes
        };
        if (boundary_events_from_trace(trace, &event_options, events) != 0)
            goto fail;
        size_t write = 0;
        for (size_t i = 0; i < events->events.count; i++) {
            const alea_ray_boundary_event_t event = events->events.data[i];
            if (event.t + RAY_EPSILON < plan.t_min ||
                event.t > t_max + RAY_EPSILON)
                continue;
            events->events.data[write++] = event;
        }
        events->events.count = write;
        if ((plan.max_events && write > plan.max_events) ||
            (plan.max_output_bytes &&
             write > plan.max_output_bytes / sizeof(*events->events.data)))
            goto fail;
        return 0;
    }

    if (plan.product == ALEA_RAY_QUERY_SEGMENTS) {
        ray_query_clip_segments(trace, plan.t_min);
        const size_t count = trace->segments.count;
        if ((plan.max_events && count > plan.max_events) ||
            (plan.max_output_bytes &&
             count > plan.max_output_bytes / sizeof(*trace->segments.data)))
            goto fail;
        return 0;
    }

    for (size_t i = 0; i < trace->segments.count; i++) {
        const alea_ray_segment_t* seg = &trace->segments.data[i];
        if (seg->t_exit <= plan.t_min + RAY_EPSILON)
            continue;

        if (plan.product == ALEA_RAY_QUERY_FIRST_CELL) {
            if (!ray_plan_accepts_material(&plan, seg, false)) continue;
            if (output) {
                output->first_cell_id = seg->cell_id;
                output->first_cell_t = fmax(seg->t_enter, plan.t_min);
            }
            return 0;
        }

        if (!ray_plan_accepts_material(&plan, seg, true)) continue;
        if (plan.product == ALEA_RAY_QUERY_ANY_HIT) {
            if (output) output->any_hit = true;
            return 0;
        }

        /* FIRST_VISIBLE: an interval clipped by t_min has no reportable
         * geometry boundary, matching the renderer's cross-section rule. */
        if (output) {
            alea_ray_first_visible_result_t* hit = &output->first_visible;
            const bool crossed_boundary =
                seg->t_enter >= plan.t_min - RAY_EPSILON;
            const alea_ray_hit_t* enter = crossed_boundary
                ? ray_query_enter_hit(trace, seg) : NULL;
            hit->found = true;
            hit->t = fmax(seg->t_enter, plan.t_min);
            hit->cell_id = seg->cell_id;
            hit->material_id = seg->material_id;
            hit->density = seg->density;
            hit->resolution_flags = seg->resolution_flags;
            if (enter) {
                hit->surface_id = enter->surface_id;
                hit->primitive_id = enter->primitive_id;
                if (plan.requirements.need_normal) {
                    hit->nx = enter->nx;
                    hit->ny = enter->ny;
                    hit->nz = enter->nz;
                }
            }
        }
        return 0;
    }
    return 0;

fail:
    alea_raycast_result_clear(trace);
    if (events) alea_ray_boundary_event_result_clear(events);
    ray_query_output_clear(output);
    return -1;
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
                                   NULL, NULL, NULL, 0, -1, 0, -1, false, result);
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
    if (!sys) return -1;
    alea_ray_t ray;
    if (alea_ray_init(&ray, ox, oy, oz, dx, dy, dz) != 0) return -1;
    if (alea_raycast_ensure_hier_caches(sys) != 0) return -1;
    static ALEA_THREAD_LOCAL alea_raycast_result_t hier_tls_result;
    static ALEA_THREAD_LOCAL int hier_tls_init = 0;
    if (!hier_tls_init) {
        alea_raycast_result_init(&hier_tls_result);
        hier_tls_init = 1;
    }
    int first_cell = -1;
    if (alea_raycast_hier_first_cell_nocache(
            sys, &ray, 0, t_max, -1, &hier_tls_result, &first_cell, out_t) != 0)
        return -1;
    return first_cell;
}

int alea_ray_is_occluded(alea_system_t* sys,
                        double ox, double oy, double oz,
                        double dx, double dy, double dz,
                        double t_max) {
    if (!sys) return 0;

    if (raycast_prefers_hier_mode(sys)) {
        if (alea_raycast_ensure_hier_caches(sys) != 0) return 0;
        static ALEA_THREAD_LOCAL alea_raycast_result_t hier_tls_result;
        static ALEA_THREAD_LOCAL int hier_tls_init = 0;
        if (!hier_tls_init) {
            alea_raycast_result_init(&hier_tls_result);
            hier_tls_init = 1;
        }
        alea_ray_t ray;
        if (alea_ray_init(&ray, ox, oy, oz, dx, dy, dz) != 0) return 0;
        int occluded = 0;
        if (alea_raycast_hier_any_hit_nocache(
                sys, &ray, 0, t_max, -1, &hier_tls_result, &occluded) != 0)
            return 0;
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

/* Persistent portion of one hierarchical ray walk.  Keeping this explicit
 * lets the scalar adapter run to completion while a future packet executor
 * advances the same walk one verified interval at a time. */
typedef struct {
    double t_current;
    int prev_cell_idx;
    int prev_surface_id;
    int pending_enter_hit_index;
    alea_raycast_boundary_event_t enter_event;
    alea_hier_ray_path_t current_path;
    double pending_lattice_entry_sample;
    int iterations_remaining;
} alea_ray_walk_t;

/* One verified open interval produced by the selected-owner walk.  This is
 * deliberately private while compatibility APIs still publish their legacy
 * segment/hit vectors. Keeping it separate from those vectors makes the
 * geometric step independent of path capture and hit materialization. */
typedef struct {
    double t_enter;
    double t_exit;
    int cell_index;
    int cell_id;
    int material_id;
    double density;
    int enter_surface_id;
    int exit_surface_id;
    alea_raycast_boundary_event_t enter_event;
    alea_raycast_boundary_event_t exit_event;
    const alea_hier_ray_path_t* path;
    uint8_t resolution_flags;
} alea_ray_selected_interval_t;

static void alea_ray_walk_init(alea_ray_walk_t* state) {
    memset(state, 0, sizeof(*state));
    state->prev_cell_idx = -2;
    state->prev_surface_id = -1;
    state->pending_enter_hit_index = -1;
    state->enter_event.surface_id = -1;
    state->enter_event.primitive_id = ALEA_PRIMITIVE_ID_INVALID;
    alea_matrix_identity(&state->enter_event.transform);
    state->iterations_remaining = 10000;
}

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

/* Populate the one normal needed by FIRST_VISIBLE without allocating or
 * appending a hit record. */
static void boundary_event_first_visible_normal(
    alea_system_t* sys, const alea_ray_t* ray,
    const alea_raycast_boundary_event_t* ev,
    bool include_normal, alea_ray_first_visible_result_t* visible) {
    if (!ev->has_physical_surface ||
        ev->primitive_id == ALEA_PRIMITIVE_ID_INVALID ||
        ev->primitive_id >= alea_vec_count(&sys->primitives))
        return;
    visible->surface_id = ev->surface_id;
    visible->primitive_id = ev->primitive_id;
    if (!include_normal) return;
    const alea_primitive_entry_t* prim = &sys->primitives.data[ev->primitive_id];
    alea_primitive_data_t prim_data;
    if (!raycast_primitive_copy_payload(sys, ev->primitive_id,
                                        prim->type, &prim_data))
        return;
    double wx, wy, wz;
    alea_ray_point_at(ray, ev->t, &wx, &wy, &wz);
    alea_matrix_transform_point_inverse(&ev->transform, &wx, &wy, &wz);
    double nlx, nly, nlz;
    primitive_normal_at(prim->type, &prim_data, wx, wy, wz,
                        &nlx, &nly, &nlz);
    boundary_event_world_normal(&ev->transform, nlx, nly, nlz,
                                &visible->nx, &visible->ny, &visible->nz);
}

/* Compatibility publication adapter for a verified selected interval. The
 * walker does not need public paths, hit indices, or segment storage to
 * establish ownership and endpoints. */
static int ray_selected_interval_publish(
    alea_system_t* sys, const alea_ray_t* ray,
    const alea_ray_selected_interval_t* interval, bool use_hier_lookup,
    bool emit_hits, alea_raycast_result_t* result, int* io_prev_cell_idx,
    int* io_pending_enter_hit_index) {
    int exit_hit_index = -1;
    if (emit_hits) {
        exit_hit_index = boundary_event_emit_hit(
            sys, ray, &interval->exit_event, result);
    }

    uint32_t path_index = UINT32_MAX;
    if (use_hier_lookup && capture_hier_path(result, interval->path,
                                              &path_index) != 0) {
        return -1;
    }

    if (interval->cell_index == *io_prev_cell_idx &&
        result->segments.count > 0 &&
        result->segments.data[result->segments.count - 1].path_index == path_index) {
        alea_ray_segment_t* previous =
            &result->segments.data[result->segments.count - 1];
        previous->t_exit = interval->t_exit;
        previous->exit_surface_id = interval->exit_surface_id;
    } else {
        const alea_ray_segment_t segment = {
            .t_enter = interval->t_enter,
            .t_exit = interval->t_exit,
            .cell_id = interval->cell_id,
            .material_id = interval->material_id,
            .density = interval->density,
            .enter_surface_id = interval->enter_surface_id,
            .exit_surface_id = interval->exit_surface_id,
            .enter_hit_index = *io_pending_enter_hit_index,
            .resolution_flags = interval->resolution_flags,
            .path_index = path_index
        };
        if (add_segment(result, &segment) != 0) return -1;
    }

    *io_prev_cell_idx = interval->cell_index;
    if (emit_hits) *io_pending_enter_hit_index = exit_hit_index;
    return 0;
}

static bool ray_selected_interval_first_visible(
    alea_system_t* sys, const alea_ray_t* ray,
    const alea_ray_selected_interval_t* interval, double t_min,
    int material_filter, bool include_normal,
    alea_ray_first_visible_result_t* out_visible) {
    if (interval->cell_id < 0 || interval->material_id == 0 ||
        (material_filter >= 0 && interval->material_id != material_filter) ||
        interval->t_exit <= t_min + RAY_EPSILON) {
        return false;
    }
    out_visible->found = true;
    out_visible->t = fmax(interval->t_enter, t_min);
    out_visible->cell_id = interval->cell_id;
    out_visible->material_id = interval->material_id;
    out_visible->density = interval->density;
    out_visible->resolution_flags = interval->resolution_flags;
    if (interval->t_enter >= t_min - RAY_EPSILON) {
        out_visible->surface_id = interval->enter_event.surface_id;
        boundary_event_first_visible_normal(
            sys, ray, &interval->enter_event, include_normal, out_visible);
    }
    return true;
}

/**
 * Test ray against a specific cell's surfaces only.
 * Returns closest hit distance, or DBL_MAX if no hit.
 *
 * When `out_primitive_id` is non-NULL it receives the canonical primitive id
 * of the closest surface (ALEA_PRIMITIVE_ID_INVALID if none), so the caller
 * can build a boundary event without re-scanning the cell.
 */
static double raycast_tree_closest_surface(alea_system_t* sys,
                                           const alea_ray_t* ray,
                                           alea_node_id_t node_id,
                                           double t_min,
                                           double t_max,
                                           alea_raycast_result_t* result,
                                           int* out_surface_id,
                                           uint32_t* out_primitive_id) {
    if (node_id >= alea_vec_count(&sys->nodes)) return DBL_MAX;
    const alea_node_t* node = &sys->nodes.data[node_id];
    alea_operation_t op = ALEA_GET_OPERATION(node);
    if (op == ALEA_OP_PRIMITIVE) {
        if (node->primitive.primitive_id >= alea_vec_count(&sys->primitives))
            return DBL_MAX;
        const alea_primitive_entry_t* prim =
            &sys->primitives.data[node->primitive.primitive_id];
        alea_primitive_data_t prim_data;
        if (!raycast_primitive_copy_payload(sys, node->primitive.primitive_id,
                                            prim->type, &prim_data)) {
            return DBL_MAX;
        }
        if (result) {
            result->surfaces_tested++;
            if ((unsigned)prim->type < ALEA_RAYCAST_PRIM_TYPE_BINS)
                result->prim_type_tests[prim->type]++;
        }
        double closest = DBL_MAX;
        double intersections[4];
        int count = ray_intersect_primitive(ray, prim->type, &prim_data,
                                            intersections);
        for (int i = 0; i < count; i++) {
            if (intersections[i] > t_min && intersections[i] < t_max &&
                intersections[i] < closest) {
                closest = intersections[i];
            }
        }
        if (closest < DBL_MAX) {
            *out_surface_id = node->primitive.mc_surface_id;
            if (out_primitive_id) *out_primitive_id = node->primitive.primitive_id;
        }
        return closest;
    }

    double left = raycast_tree_closest_surface(
        sys, ray, node->operation.left, t_min, t_max, result,
        out_surface_id, out_primitive_id);
    if (op == ALEA_OP_COMPLEMENT) return left;
    int right_surface_id = -1;
    uint32_t right_primitive_id = ALEA_PRIMITIVE_ID_INVALID;
    double right = raycast_tree_closest_surface(
        sys, ray, node->operation.right, t_min, t_max, result,
        &right_surface_id, &right_primitive_id);
    if (right < left) {
        *out_surface_id = right_surface_id;
        if (out_primitive_id) *out_primitive_id = right_primitive_id;
        return right;
    }
    return left;
}

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

    /* A few imported/container cells deliberately have no usable surface
     * index.  Their CSG is still exact occurrence support, so retain a
     * bounded tree fallback for ancestor-event queries and traversal. */
    if (!cell->surface_indices || cell->surface_index_count == 0) {
        return raycast_tree_closest_surface(
            sys, ray, cell->root_node_id, t_min, t_max, result,
            out_surface_id, out_primitive_id);
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

static alea_raycast_ancestor_cell_stat_t*
raycast_ancestor_cell_stat(alea_raycast_result_t* result,
                           uint32_t cell_index, int cell_id) {
    for (size_t i = 0; i < ALEA_RAYCAST_ANCESTOR_HOT_CELLS; ++i) {
        alea_raycast_ancestor_cell_stat_t* stat = &result->ancestor_hot_cells[i];
        if (stat->queries != 0 && stat->cell_index == cell_index)
            return stat;
    }
    for (size_t i = 0; i < ALEA_RAYCAST_ANCESTOR_HOT_CELLS; ++i) {
        alea_raycast_ancestor_cell_stat_t* stat = &result->ancestor_hot_cells[i];
        if (stat->queries == 0) {
            stat->cell_index = cell_index;
            stat->cell_id = cell_id;
            return stat;
        }
    }
    return NULL;
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
    if (result) {
        uint32_t depth = (uint32_t)(path->count - 1);
        result->ancestor_surface_queries++;
        result->ancestor_path_entries_examined += depth;
        if (depth > result->ancestor_max_path_depth)
            result->ancestor_max_path_depth = depth;
    }

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
        int tested_before = result ? result->surfaces_tested : 0;
        double t = raycast_cell_surfaces(sys, &local_ray, cell, t_min, t_max,
                                         result, &surface_id, &prim_id);
        int tested = result ? result->surfaces_tested - tested_before : 0;
        alea_raycast_ancestor_cell_stat_t* stat = result
            ? raycast_ancestor_cell_stat(result, entry->cell_index, cell->mc_cell_id)
            : NULL;
        if (stat) {
            stat->queries++;
            stat->surface_tests += (uint64_t)tested;
        } else if (result) {
            result->ancestor_unattributed_queries++;
            result->ancestor_unattributed_surface_tests += (uint64_t)tested;
        }
        if (t < closest_t) {
            closest_t = t;
            *out_surface_id = surface_id;
            if (out_primitive_id) *out_primitive_id = prim_id;
            if (out_transform) *out_transform = entry->transform;
            if (stat) stat->winning_events++;
            else if (result) result->ancestor_unattributed_winning_events++;
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

    alea_bbox_t lat_bbox;
    if (lattice_rect_bounds(lat_cell, &lat_bbox) != 0) return t_max;

    double t_enter = 0.0, t_exit = t_max;
    if (!lat_cell->lat_fill_repeating &&
        !ray_bbox_slab_enter_exit(ray, &lat_bbox, 0.0, t_max,
                                  &t_enter, &t_exit)) return t_max;

    double sample_t = t_min + RAY_EPSILON;
    double sx = ray->ox + sample_t * ray->dx;
    double sy = ray->oy + sample_t * ray->dy;
    double sz = ray->oz + sample_t * ray->dz;

    int i = (int)floor((sx - ll[0]) / px);
    int j = (lat_cell->lat_fill_repeating || nj > 1)
        ? (int)floor((sy - ll[1]) / py) : 0;
    int k = (lat_cell->lat_fill_repeating || nk > 1)
        ? (int)floor((sz - ll[2]) / pz) : 0;

    if (!lat_cell->lat_fill_repeating) {
        if (i < 0) i = 0;
        if (i >= ni) i = ni - 1;
        if (j < 0) j = 0;
        if (j >= nj) j = nj - 1;
        if (k < 0) k = 0;
        if (k >= nk) k = nk - 1;
    }

    (void)t_enter;
    (void)t_exit;
    double t_next = t_max;
    if (fabs(ray->dx) > RAY_EPSILON) {
        double boundary = ll[0] + ((ray->dx > 0.0) ? (i + 1) : i) * px;
        double t = (boundary - ray->ox) / ray->dx;
        if (t > t_min + RAY_EPSILON && t < t_next) t_next = t;
    }
    if ((lat_cell->lat_fill_repeating || nj > 1) &&
        fabs(ray->dy) > RAY_EPSILON) {
        double boundary = ll[1] + ((ray->dy > 0.0) ? (j + 1) : j) * py;
        double t = (boundary - ray->oy) / ray->dy;
        if (t > t_min + RAY_EPSILON && t < t_next) t_next = t;
    }
    if ((lat_cell->lat_fill_repeating || nk > 1) &&
        fabs(ray->dz) > RAY_EPSILON) {
        double boundary = ll[2] + ((ray->dz > 0.0) ? (k + 1) : k) * pz;
        double t = (boundary - ray->oz) / ray->dz;
        if (t > t_min + RAY_EPSILON && t < t_next) t_next = t;
    }

    return (t_next < t_max) ? t_next : t_max;
}

static double lattice_hex_next_boundary(const alea_ray_t* ray,
                                        const alea_cell_entry_t* lat_cell,
                                        double t_min,
                                        double t_max) {
    double p = lat_cell->lat_pitch[0];
    if (p <= 0.0) return t_max;

    alea_bbox_t lat_bbox;
    if (lattice_hex_bounds(lat_cell, &lat_bbox) != 0) return t_max;

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

typedef enum {
    LATTICE_ENTRY_NO_ENTRY,
    LATTICE_ENTRY_FUTURE_ENTRY,
    LATTICE_ENTRY_ALREADY_INSIDE,
    LATTICE_ENTRY_ERROR
} lattice_entry_kind_t;

typedef struct {
    lattice_entry_kind_t kind;
    double t;
    double sample;
} lattice_entry_query_t;

typedef struct {
    const alea_system_t* sys;
    const alea_ray_t* ray;
    double t_min;
    double t_max;
    double closest;
    alea_raycast_result_t* result;
} lattice_ancestor_event_ctx_t;

static int visit_lattice_ancestor_surface(uint32_t cell_index,
                                          const alea_matrix_t* transform,
                                          void* userdata) {
    lattice_ancestor_event_ctx_t* ctx = userdata;
    if (!ctx || !transform || cell_index >= alea_vec_count(&ctx->sys->cells)) {
        return -1;
    }
    const alea_cell_entry_t* cell = &ctx->sys->cells.data[cell_index];
    alea_ray_t local_ray;
    transform_ray_inverse(transform, ctx->ray, &local_ray);
    int surface_id = -1;
    uint32_t primitive_id = ALEA_PRIMITIVE_ID_INVALID;
    int tested_before = ctx->result ? ctx->result->surfaces_tested : 0;
    double t = raycast_cell_surfaces((alea_system_t*)ctx->sys, &local_ray,
                                     cell, ctx->t_min, ctx->t_max,
                                     ctx->result, &surface_id, &primitive_id);
    if (ctx->result) {
        ctx->result->lattice_entry_ancestor_surface_tests +=
            ctx->result->surfaces_tested - tested_before;
    }
    if (t < ctx->closest) ctx->closest = t;
    return 0;
}

/* Return the next exact enclosing-cell surface for this occurrence.  The
 * support AABB only found the candidate; these surfaces determine when a
 * repeating lattice can actually become active. */
static int lattice_next_ancestor_event(const alea_system_t* sys,
                                       uint32_t placement_index,
                                       const alea_ray_t* ray,
                                       double t_min,
                                       double t_max,
                                       alea_raycast_result_t* result,
                                       double* out_t) {
    if (!out_t) return -1;
    *out_t = DBL_MAX;
    lattice_ancestor_event_ctx_t ctx = {
        .sys = sys,
        .ray = ray,
        .t_min = t_min,
        .t_max = t_max,
        .closest = DBL_MAX,
        .result = result
    };
    int rc = alea_hier_spatial_visit_lattice_placement_ancestors(
        (alea_system_t*)sys, placement_index,
        visit_lattice_ancestor_surface, &ctx);
    if (rc != 0) return -1;
    *out_t = ctx.closest;
    if (result && ctx.closest < DBL_MAX) {
        result->lattice_entry_ancestor_events++;
    }
    return 0;
}

static int lattice_placement_ancestors_active(const alea_system_t* sys,
                                               uint32_t placement_index,
                                               const alea_ray_t* ray,
                                               double t,
                                               double t_max) {
    double sample = t + RAY_EPSILON;
    if (sample >= t_max) return 0;
    double x, y, z;
    alea_ray_point_at(ray, sample, &x, &y, &z);
    return alea_hier_spatial_check_lattice_placement_ancestors(
        (alea_system_t*)sys, placement_index, x, y, z);
}

/* Find the first interval after t_min where both the lattice element and its
 * enclosing occurrence are active.  A repeating lattice is locally valid in
 * every pitch, so while the exact ancestor CSG is false we jump directly to
 * its next surface rather than walking pitches that cannot contribute. */
static lattice_entry_query_t lattice_first_valid_entry(
    const alea_system_t* sys,
    const alea_ray_t* local_ray,
    const alea_ray_t* world_ray,
    const alea_cell_entry_t* cell,
    uint32_t placement_index,
    double t_min,
    double t_max,
    alea_raycast_result_t* result) {
    lattice_entry_query_t out = {
        .kind = LATTICE_ENTRY_NO_ENTRY,
        .t = DBL_MAX,
        .sample = -1.0
    };
    double t_enter, t_exit;
    if (lattice_raycast_interval(local_ray, cell, t_min, t_max,
                                 &t_enter, &t_exit) != 0 ||
        t_exit <= t_enter + RAY_EPSILON) {
        if (result) result->lattice_entry_no_entry_results++;
        return out;
    }

    const int max_steps = lattice_raycast_step_limit(local_ray, cell,
                                                      t_enter, t_exit);
    double t_current = t_enter;
    int ancestor_active = lattice_placement_ancestors_active(
        sys, placement_index, world_ray, t_current, t_exit);
    if (ancestor_active < 0) {
        out.kind = LATTICE_ENTRY_ERROR;
        return out;
    }
    double t_ancestor = DBL_MAX;
    if (lattice_next_ancestor_event(sys, placement_index, world_ray,
                                    t_current + RAY_EPSILON, t_exit,
                                    result, &t_ancestor) != 0) {
        out.kind = LATTICE_ENTRY_ERROR;
        return out;
    }

    for (int step = 0; step < max_steps && t_current < t_exit; step++) {
        /* Outside the occurrence support, lattice periodicity is irrelevant.
         * Move to the next exact ancestor event in one step. */
        if (!ancestor_active) {
            if (t_ancestor >= t_exit - RAY_EPSILON) break;
            t_current = t_ancestor;
            ancestor_active = lattice_placement_ancestors_active(
                sys, placement_index, world_ray, t_current, t_exit);
            if (ancestor_active < 0) {
                out.kind = LATTICE_ENTRY_ERROR;
                return out;
            }
            if (lattice_next_ancestor_event(
                    sys, placement_index, world_ray,
                    t_current + RAY_EPSILON, t_exit, result,
                    &t_ancestor) != 0) {
                out.kind = LATTICE_ENTRY_ERROR;
                return out;
            }
            continue;
        }

        if (result) result->lattice_entry_dda_steps++;
        double t_lattice = lattice_next_boundary(local_ray, cell,
                                                 t_current, t_exit);
        double t_next = fmin(t_lattice, t_ancestor);
        if (t_next <= t_current + RAY_EPSILON || t_next > t_exit)
            t_next = t_exit;
        double x, y, z;
        alea_ray_point_at(local_ray, t_current + 0.5 * (t_next - t_current),
                          &x, &y, &z);
        alea_lattice_location_t location;
        const int valid = alea_lattice_locate_point(sys, cell, x, y, z,
                                                     &location) == 1;
        if (valid) {
            /* CSG containment alone is insufficient in overlapping decks:
             * publishing a synthetic transition for a shadowed fill would
             * make the next resolver sample ambiguous.  Leave the ordinary
             * global surface walker in charge until ownership changes. */
            int canonical =
                alea_hier_spatial_check_lattice_placement_canonical_ancestors(
                    (alea_system_t*)sys, placement_index,
                    world_ray->ox + (t_current + RAY_EPSILON) * world_ray->dx,
                    world_ray->oy + (t_current + RAY_EPSILON) * world_ray->dy,
                    world_ray->oz + (t_current + RAY_EPSILON) * world_ray->dz);
            if (canonical < 0) {
                out.kind = LATTICE_ENTRY_ERROR;
                return out;
            }
            if (!canonical) {
                if (result) result->lattice_entry_canonical_rejections++;
                return out;
            }
            out.t = t_current;
            out.sample = t_current +
                fmin(DDA_SAMPLE_OFFSET, 0.25 * (t_next - t_current));
            if (t_current <= t_min + RAY_EPSILON) {
                out.kind = LATTICE_ENTRY_ALREADY_INSIDE;
                if (result) result->lattice_entry_already_inside_results++;
            } else {
                out.kind = LATTICE_ENTRY_FUTURE_ENTRY;
                if (result) result->lattice_entry_future_entry_results++;
            }
            return out;
        }

        if (t_ancestor <= t_lattice + RAY_EPSILON) {
            ancestor_active = lattice_placement_ancestors_active(
                sys, placement_index, world_ray, t_next, t_exit);
            if (ancestor_active < 0) {
                out.kind = LATTICE_ENTRY_ERROR;
                return out;
            }
            if (lattice_next_ancestor_event(
                    sys, placement_index, world_ray,
                    t_next + RAY_EPSILON, t_exit, result,
                    &t_ancestor) != 0) {
                out.kind = LATTICE_ENTRY_ERROR;
                return out;
            }
        }
        t_current = t_next;
    }
    if (result) result->lattice_entry_no_entry_results++;
    return out;
}

typedef struct {
    const alea_system_t* sys;
    const alea_ray_t* ray;
    double t_min;
    double t_max;
    lattice_entry_query_t closest;
    alea_raycast_result_t* result;
} lattice_entry_visit_ctx_t;

static int visit_lattice_entry(uint32_t placement_index,
                               uint32_t lattice_cell_index,
                               const alea_matrix_t* transform,
                               double placement_enter,
                               double placement_exit,
                               void* userdata) {
    (void)placement_index;
    lattice_entry_visit_ctx_t* ctx = userdata;
    if (!ctx || !transform ||
        lattice_cell_index >= alea_vec_count(&ctx->sys->cells)) {
        return -1;
    }
    const alea_cell_entry_t* cell =
        &ctx->sys->cells.data[lattice_cell_index];
    if (ctx->result) ctx->result->lattice_entry_candidates++;
    alea_ray_t local_ray;
    transform_ray_inverse(transform, ctx->ray, &local_ray);
    lattice_entry_query_t entry = lattice_first_valid_entry(
        ctx->sys, &local_ray, ctx->ray, cell, placement_index,
        fmax(ctx->t_min, placement_enter),
        fmin(ctx->t_max, placement_exit), ctx->result);
    if (entry.kind == LATTICE_ENTRY_ERROR) return -1;
    if (entry.kind != LATTICE_ENTRY_FUTURE_ENTRY ||
        entry.t >= ctx->closest.t) {
        return 0;
    }
    {
        double x, y, z;
        alea_ray_point_at(ctx->ray, entry.sample, &x, &y, &z);
        int valid = alea_hier_spatial_check_lattice_placement_canonical_ancestors(
            (alea_system_t*)ctx->sys, placement_index, x, y, z);
        if (valid < 0) return -1;
        if (!valid) return 0;
        ctx->closest = entry;
    }
    return 0;
}

static int system_first_lattice_entry(alea_system_t* sys,
                                      const alea_ray_t* ray,
                                      double t_min,
                                      double t_max,
                                      alea_raycast_result_t* result,
                                      lattice_entry_query_t* out_entry) {
    if (!out_entry) return -1;
    out_entry->kind = LATTICE_ENTRY_NO_ENTRY;
    out_entry->t = DBL_MAX;
    out_entry->sample = -1.0;
    if (result) result->lattice_entry_calls++;
    lattice_entry_visit_ctx_t ctx = {
        .sys = sys,
        .ray = ray,
        .t_min = t_min,
        .t_max = t_max,
        .closest = {
            .kind = LATTICE_ENTRY_NO_ENTRY,
            .t = DBL_MAX,
            .sample = -1.0
        },
        .result = result
    };
    if (alea_hier_spatial_visit_lattice_placements_ray(
            sys, ray->ox, ray->oy, ray->oz,
            ray->inv_dx, ray->inv_dy, ray->inv_dz,
            t_min, t_max, visit_lattice_entry, &ctx) != 0) {
        return -1;
    }
    *out_entry = ctx.closest;
    return 0;
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
        /* Deepest cell along the first containing chain */
        const alea_cell_hit_t* hit =
            &hits[deepest_first_chain_index(hits, num_hits)];
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
        const alea_hier_ray_path_entry_t* parent_entry =
            &current_path->entries[parent];
        int found;
        if (parent_entry->is_lattice &&
            (size_t)parent_entry->cell_index < alea_vec_count(&sys->cells)) {
            const alea_cell_entry_t* lattice =
                &sys->cells.data[parent_entry->cell_index];
            double lx = px, ly = py, lz = pz;
            alea_matrix_transform_point_inverse(&parent_entry->transform,
                                                &lx, &ly, &lz);
            alea_lattice_location_t location;
            if (alea_lattice_locate_point(sys, lattice, lx, ly, lz,
                                          &location) != 1) {
                continue;
            }
            found = alea_hier_path_enter_lattice_location(
                sys, current_path, parent, px, py, pz, &location,
                &hit_with_transform, &candidate_path);
        } else {
            /* A terminal-child boundary does not change a still-containing
             * ordinary fill owner.  Preserve that tracked owner rather than
             * issuing a deck-order universe query for every ancestor restart.
             * Lattice transitions above remain canonical because their
             * selected element can change at the boundary. */
            found = alea_hier_spatial_find_path_from_parent_coherent(
                sys, current_path, parent, px, py, pz,
                &hit_with_transform, &candidate_path);
        }
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
    ent->lat_i = 0;
    ent->lat_j = 0;
    ent->lat_k = 0;
    ent->lat_ox = 0.0;
    ent->lat_oy = 0.0;
    ent->lat_oz = 0.0;
    alea_matrix_identity(&ent->transform);

    if (cell->lat_type != 0 && cell->lat_fill) {
        alea_lattice_location_t location;
        if (alea_lattice_locate_point(sys, cell, px, py, pz, &location) != 1)
            return -1;
        ent->fill_universe = location.fill_universe;
        ent->is_lattice = 1;
        ent->lat_fill_universe = location.fill_universe;
        ent->lat_i = location.i;
        ent->lat_j = location.j;
        ent->lat_k = location.k;
        ent->lat_ox = location.ox;
        ent->lat_oy = location.oy;
        ent->lat_oz = location.oz;
    }

    if (cell->fill_universe > 0 || ent->is_lattice) {
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

static int raycast_cell_aware_resume(alea_system_t* sys,
                                     const alea_ray_t* ray,
                                     double effective_t_max,
                                     bool use_hier_lookup,
                                     bool emit_hits,
                                     alea_ray_first_visible_result_t* first_visible,
                                     int* first_cell_id,
                                     double* first_cell_t,
                                     double first_cell_t_min,
                                     int first_cell_material_filter,
                                     double visible_t_min,
                                     int visible_material_filter,
                                     bool visible_wants_normal,
                                     alea_raycast_result_t* result,
                                     alea_ray_walk_t* state,
                                     int step_budget,
                                     alea_ray_selected_interval_t* out_selected) {
    if (!state || step_budget <= 0) return -1;
    /* Current position along ray */
    double t_current = state->t_current;
    int prev_cell_idx = state->prev_cell_idx;
    int prev_surface_id = state->prev_surface_id;
    /* Hit index of the surface crossed at t_current (the next segment's
     * enter boundary), or -1. Only maintained when emit_hits is set. */
    int pending_enter_hit_index = state->pending_enter_hit_index;
    /* The physical or synthetic event that entered the current interval.
     * FIRST_VISIBLE reports an entering surface, so it must not reuse the
     * current iteration's exit event. */
    alea_raycast_boundary_event_t enter_event = state->enter_event;
    alea_hier_ray_path_t current_path = state->current_path;
    double pending_lattice_entry_sample = state->pending_lattice_entry_sample;
    const bool need_boundary_event = emit_hits || first_visible || out_selected;

    while (t_current < effective_t_max && step_budget-- > 0 &&
           state->iterations_remaining-- > 0) {
        result->step_iterations++;
        /* Cell resolution samples just past the crossing at t_current, where
         * quadric sign evaluation is numerically noisy. Midpoint verification
         * after the crossing search jumps back here once with a sample point
         * in the segment interior when the resolved cell fails to contain it. */
        double t_sample = pending_lattice_entry_sample > t_current
            ? pending_lattice_entry_sample : t_current + RAY_EPSILON;
        pending_lattice_entry_sample = -1.0;
        int resolve_attempt = 0;
        /* Snapshot of the attempt-0 outcome. The retry is accepted only when
         * it verifies strictly better; otherwise this state is restored, so a
         * failed retry can never degrade the pre-verification answer. */
        int saved_valid = 0;
        int saved_cell_idx = -1;
        int saved_cell_id = -1;
        int saved_material_id = 0;
        double saved_density = 0;
        double saved_t_next = 0;
        int saved_hit_surface_id = -1;
        int saved_next_enter_surface_id = -1;
        alea_raycast_boundary_event_t saved_bevent;
        alea_hier_ray_path_t saved_path;
        memset(&saved_bevent, 0, sizeof(saved_bevent));
        saved_path.count = 0;
resolve_cell:;
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
                double px, py, pz;
                alea_ray_point_at(ray, t_sample, &px, &py, &pz);
                found_via_neighbor = raycast_find_neighbor(sys, prev_cell_idx,
                                                           prev_surface_id,
                                                           px, py, pz,
                                                           &cell_id, &cell_idx,
                                                           &material_id,
                                                           &density);
                /* In hier mode a filled container is not a terminal answer:
                 * reject it so the full lookup descends into the fill. */
                if (found_via_neighbor && use_hier_lookup) {
                    const alea_cell_entry_t* nb = &sys->cells.data[cell_idx];
                    if (nb->fill_universe > 0 ||
                        (nb->lat_type != 0 && nb->lat_fill)) {
                        found_via_neighbor = 0;
                        cell_idx = -1;
                        cell_id = -1;
                        material_id = 0;
                        density = 0;
                    } else if (nb->universe_id == 0) {
                        /* Neighbor resolution at the root bypasses the full
                         * hierarchy lookup.  Refresh its path explicitly:
                         * retaining the previous root entry would attach the
                         * wrong owner occurrence to this segment. */
                        current_path.count = 1;
                        alea_hier_ray_path_entry_t* entry =
                            &current_path.entries[0];
                        entry->cell_index = (uint32_t)cell_idx;
                        entry->cell_id = nb->mc_cell_id;
                        entry->material_id = nb->material_id;
                        entry->universe_id = nb->universe_id;
                        entry->fill_universe = nb->fill_universe;
                        entry->depth = 0;
                        entry->is_lattice = 0;
                        entry->lat_fill_universe = 0;
                        entry->lat_i = 0;
                        entry->lat_j = 0;
                        entry->lat_k = 0;
                        entry->lat_ox = 0.0;
                        entry->lat_oy = 0.0;
                        entry->lat_oz = 0.0;
                        alea_matrix_identity(&entry->transform);
                    }
                }
            }
        }

        /* In hierarchical mode, avoid a full root-to-deepest lookup when the
         * sampled point remains inside the cached placement path. */
        if (!found_via_neighbor && use_hier_lookup && prev_cell_idx >= 0 &&
            (size_t)prev_cell_idx < alea_vec_count(&sys->cells)) {
            const alea_cell_entry_t* prev_cell = &sys->cells.data[prev_cell_idx];
            double px, py, pz;
            alea_ray_point_at(ray, t_sample, &px, &py, &pz);
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
                    sys, &current_path, current_path.count - 1, px, py, pz,
                    &cell_transform, &lattice_cell_index, &lattice_transform);
                if (in_cell < 0) in_cell = 0;
            }
            if (in_cell && !alea_cell_entry_is_container(prev_cell)) {
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
            alea_ray_point_at(ray, t_sample, &px, &py, &pz);
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
        if (need_boundary_event) {
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
            if (need_boundary_event) {
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
                    if (fabs(t_lattice_surface - t_lattice_next) <=
                        SURFACE_SAMPLE_OFFSET) {
                        const double post_t =
                            fmax(t_lattice_surface, t_lattice_next) +
                            10.0 * RAY_EPSILON;
                        double post_x, post_y, post_z;
                        alea_ray_point_at(&lattice_ray, post_t,
                                          &post_x, &post_y, &post_z);
                        alea_lattice_location_t post_location;
                        if (alea_lattice_locate_point(
                                sys, lattice_cell, post_x, post_y, post_z,
                                &post_location) != 1) {
                            pending_lattice_entry_sample = post_t;
                        }
                    }
                    if (need_boundary_event) {
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
                        if (need_boundary_event) {
                            bevent.t = t_lattice_next;
                            bevent.surface_id = lattice_surface_id;
                            bevent.primitive_id = lattice_prim_id;
                            bevent.has_physical_surface = (lattice_surface_id > 0);
                            bevent.is_synthetic_lattice_boundary = true;
                            bevent.transform = lattice_transform;
                        }
                    } else {
                        hit_surface_id = 0;
                        if (need_boundary_event) {
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
                    if (need_boundary_event) {
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
                    if (need_boundary_event) {
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
            if (need_boundary_event) {
                /* Void-region global search runs in the world frame. */
                bevent.t = t_next;
                bevent.surface_id = hit_surface_id;
                bevent.primitive_id = void_prim_id;
                bevent.has_physical_surface = (hit_surface_id > 0);
                bevent.is_synthetic_lattice_boundary = false;
                alea_matrix_identity(&bevent.transform);
            }
        }

        if (use_hier_lookup && lattice_cell_index < 0 &&
            (cell_idx < 0 ||
             ((size_t)cell_idx < alea_vec_count(&sys->cells) &&
              alea_cell_entry_is_container(&sys->cells.data[cell_idx])))) {
            lattice_entry_query_t lattice_entry;
            if (system_first_lattice_entry(
                sys, ray, t_current + RAY_EPSILON, effective_t_max,
                result, &lattice_entry) != 0) {
                return -1;
            }
            if (lattice_entry.kind == LATTICE_ENTRY_FUTURE_ENTRY &&
                lattice_entry.t < t_next - SURFACE_SAMPLE_OFFSET) {
                t_next = lattice_entry.t;
                hit_surface_id = 0;
                next_enter_surface_id = 0;
                pending_lattice_entry_sample = lattice_entry.sample;
                if (need_boundary_event) {
                    bevent.t = lattice_entry.t;
                    bevent.surface_id = 0;
                    bevent.primitive_id = ALEA_PRIMITIVE_ID_INVALID;
                    bevent.has_physical_surface = false;
                    bevent.is_synthetic_lattice_boundary = true;
                }
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

        /* Verify the resolved cell actually owns this segment. No surface of
         * the chosen cell lies inside the open interval, so containment is
         * constant across it: probing an interior point — well away from the
         * numerically noisy boundary sample — catches both bogus adjacency
         * transitions (non-convex cells re-entered across their own surface)
         * and on-boundary sign noise picking a coincident sibling cell. On
         * failure, redo the resolution once, sampling at the probe point. */
        if (cell_idx >= 0 && (size_t)cell_idx < alea_vec_count(&sys->cells) &&
            t_next - t_current > 1e-6) {
            /* Probe at an irrational fraction of the interval, not the exact
             * midpoint: a segment spanning lattice elements has periodic
             * internal boundaries, and its midpoint can land exactly on one,
             * where containment is numerically undefined. */
            double t_probe = t_current +
                0.381966011250105 * (t_next - t_current);
            double mx, my, mz;
            alea_ray_point_at(ray, t_probe, &mx, &my, &mz);
            const alea_cell_entry_t* seg_cell = &sys->cells.data[cell_idx];
            int probe_ok;
            if (!use_hier_lookup && seg_cell->universe_id != 0) {
                /* Flat mode resolves nested cells through the flat spatial
                 * query; their CSG lives in a local frame we don't track
                 * here, so containment against world coordinates would be
                 * meaningless. Trust the resolution. */
                probe_ok = 1;
            } else if (!use_hier_lookup || seg_cell->universe_id == 0) {
                probe_ok = seg_cell->root_node_id != ALEA_NODE_ID_INVALID &&
                    alea_contains_point(sys, seg_cell->root_node_id,
                                        mx, my, mz);
            } else {
                probe_ok = alea_hier_spatial_check_path_containment(
                    sys, &current_path, current_path.count - 1, mx, my, mz,
                    NULL, NULL, NULL) > 0;
            }
            /* A retry must not climb the hierarchy: a shallower answer than
             * attempt-0's is a "found nothing deeper" fallback (e.g. a fill
             * container over an unsupported lattice region), not a
             * correction. Equal-depth siblings (overlaps, re-crossed
             * quartics) and deeper resolutions are legitimate. */
            int depth_ok = !saved_valid || !use_hier_lookup ||
                           current_path.count >= saved_path.count;
            if (!probe_ok || !depth_ok) {
                if (resolve_attempt < 2 && probe_ok == 0 && depth_ok) {
                    /* Re-resolve, sampling at this interval's probe point.
                     * A second retry handles the case where the region just
                     * past t_current belongs to yet another cell than the
                     * one found at the first retry's (stale) probe. */
                    if (!saved_valid) {
                        saved_valid = 1;
                        saved_cell_idx = cell_idx;
                        saved_cell_id = cell_id;
                        saved_material_id = material_id;
                        saved_density = density;
                        saved_t_next = t_next;
                        saved_hit_surface_id = hit_surface_id;
                        saved_next_enter_surface_id = next_enter_surface_id;
                        saved_bevent = bevent;
                        saved_path = current_path;
                    }
                    resolve_attempt++;
                    t_sample = t_probe;
                    goto resolve_cell;
                }
                if (saved_valid && !depth_ok) {
                    /* A retry that climbed the hierarchy is a "found nothing
                     * deeper" fallback, never a correction: keep attempt 0.
                     * A depth-consistent retry that merely fails its probe is
                     * kept as-is — its error is bounded by its own short
                     * interval, whereas restoring attempt-0 would reinstate
                     * the full overrun that triggered the retry. */
                    cell_idx = saved_cell_idx;
                    cell_id = saved_cell_id;
                    material_id = saved_material_id;
                    density = saved_density;
                    t_next = saved_t_next;
                    hit_surface_id = saved_hit_surface_id;
                    next_enter_surface_id = saved_next_enter_surface_id;
                    bevent = saved_bevent;
                    current_path = saved_path;
                }
            }
        } else if (resolve_attempt >= 1 && saved_valid) {
            /* Retry resolved to void or a degenerate interval: keep the
             * attempt-0 answer. */
            cell_idx = saved_cell_idx;
            cell_id = saved_cell_id;
            material_id = saved_material_id;
            density = saved_density;
            t_next = saved_t_next;
            hit_surface_id = saved_hit_surface_id;
            next_enter_surface_id = saved_next_enter_surface_id;
            bevent = saved_bevent;
            current_path = saved_path;
        }

        bevent.t = t_next;
        const alea_ray_selected_interval_t selected = {
            .t_enter = t_current,
            .t_exit = t_next,
            .cell_index = cell_idx,
            .cell_id = cell_id,
            .material_id = material_id,
            .density = density,
            .enter_surface_id = prev_surface_id,
            .exit_surface_id = hit_surface_id,
            .enter_event = enter_event,
            .exit_event = bevent,
            .path = &current_path,
            .resolution_flags =
                (cell_idx >= 0 &&
                 (size_t)cell_idx < alea_vec_count(&sys->cells) &&
                 alea_cell_entry_is_container(&sys->cells.data[cell_idx]))
                    ? ALEA_RESOLVE_UNDEFINED_FILL : 0
        };

        /* The selected-interval facade consumes the geometric result before
         * any compatibility publication or query-specific stop policy. */
        if (out_selected) {
            prev_cell_idx = cell_idx;
            if (next_enter_surface_id < 0)
                next_enter_surface_id = hit_surface_id;
            prev_surface_id = next_enter_surface_id;
            enter_event = bevent;
            t_current = t_next;
            state->t_current = t_current;
            state->prev_cell_idx = prev_cell_idx;
            state->prev_surface_id = prev_surface_id;
            state->pending_enter_hit_index = pending_enter_hit_index;
            state->enter_event = enter_event;
            state->current_path = current_path;
            state->pending_lattice_entry_sample = pending_lattice_entry_sample;
            *out_selected = selected;
            out_selected->path = &state->current_path;
            return 0;
        }

        /* FIRST_CELL and FIRST_VISIBLE consume a verified interval before
         * hit/segment materialization, so work behind the answer is avoided. */
        if (first_cell_id && selected.cell_id >= 0 &&
            (first_cell_material_filter < 0 ||
             selected.material_id == first_cell_material_filter) &&
            selected.t_exit > first_cell_t_min + RAY_EPSILON) {
            *first_cell_id = selected.cell_id;
            if (first_cell_t) *first_cell_t = fmax(selected.t_enter, first_cell_t_min);
            return 0;
        }

        /* FIRST_VISIBLE consumes this completed interval immediately.  It is
         * deliberately after ownership verification but before hit/segment
         * materialization, so work behind the winning interval is avoided. */
        if (first_visible && selected.cell_id >= 0 && selected.material_id != 0 &&
            (visible_material_filter < 0 ||
             selected.material_id == visible_material_filter) &&
            selected.t_exit > visible_t_min + RAY_EPSILON) {
            memset(first_visible, 0, sizeof(*first_visible));
            first_visible->found = true;
            first_visible->t = fmax(selected.t_enter, visible_t_min);
            first_visible->cell_id = selected.cell_id;
            first_visible->material_id = selected.material_id;
            first_visible->density = selected.density;
            /* Synthetic lattice entries are observable boundaries in the
             * canonical trace (surface id 0).  Preserve that identity only
             * when visibility begins at the interval entry; a clipped query
             * begins in the interval interior and has no entering boundary. */
            first_visible->surface_id =
                selected.t_enter >= visible_t_min - RAY_EPSILON
                    ? selected.enter_event.surface_id : -1;
            first_visible->primitive_id = UINT32_MAX;
            first_visible->resolution_flags = selected.resolution_flags;
            if (selected.t_enter >= visible_t_min - RAY_EPSILON)
                boundary_event_first_visible_normal(sys, ray, &selected.enter_event,
                                                    visible_wants_normal,
                                                    first_visible);
            return 0;
        }

        if (!first_visible && !first_cell_id) {
            if (ray_selected_interval_publish(
                    sys, ray, &selected, use_hier_lookup, emit_hits, result,
                    &prev_cell_idx, &pending_enter_hit_index) != 0)
                return -1;
        } else {
            /* Preserve the stepper's neighbor/path locality without retaining
             * the interval that FIRST_VISIBLE has already classified. */
            prev_cell_idx = cell_idx;
        }

        if (next_enter_surface_id < 0)
            next_enter_surface_id = hit_surface_id;
        prev_surface_id = next_enter_surface_id;
        enter_event = bevent;

        /* Move past the intersection */
        t_current = t_next;

        /* Persist only the state needed to resume at the next verified
         * interval. Per-step candidate and retry variables remain local. */
        state->t_current = t_current;
        state->prev_cell_idx = prev_cell_idx;
        state->prev_surface_id = prev_surface_id;
        state->pending_enter_hit_index = pending_enter_hit_index;
        state->enter_event = enter_event;
        state->current_path = current_path;
        state->pending_lattice_entry_sample = pending_lattice_entry_sample;

        /* If we hit nothing, we're done */
        if (t_next >= effective_t_max - RAY_EPSILON) {
            break;
        }
    }

    if (t_current < effective_t_max && state->iterations_remaining > 0 &&
        step_budget <= 0) {
        return 1;  /* yielded after the requested number of intervals */
    }
    return 0;
}

/* The compatibility adapter runs the selected walker to completion.  Packet
 * and fixed-output consumers use the one-step façade below instead. */
static int raycast_cell_aware_impl(alea_system_t* sys,
                                   const alea_ray_t* ray,
                                   double effective_t_max,
                                   bool use_hier_lookup,
                                   bool emit_hits,
                                   alea_ray_first_visible_result_t* first_visible,
                                   int* first_cell_id,
                                   double* first_cell_t,
                                   double first_cell_t_min,
                                   int first_cell_material_filter,
                                   double visible_t_min,
                                   int visible_material_filter,
                                   bool visible_wants_normal,
                                   alea_raycast_result_t* result) {
    alea_ray_walk_t state;
    alea_ray_walk_init(&state);
    return raycast_cell_aware_resume(
        sys, ray, effective_t_max, use_hier_lookup, emit_hits,
        first_visible, first_cell_id, first_cell_t, first_cell_t_min,
        first_cell_material_filter, visible_t_min, visible_material_filter,
        visible_wants_normal, result, &state, INT_MAX, NULL);
}

/* Advance the coherent hierarchical walk by one verified interval. The
 * returned interval is scratch-backed through `walk` and must be consumed
 * before the next call. Returns 1 when another interval may follow, 0 for
 * the terminal interval, and 2 when the walk was already exhausted. */
static int alea_ray_walk_next_selected(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    alea_raycast_result_t* scratch, alea_ray_walk_t* walk,
    alea_ray_selected_interval_t* out_interval) {
    if (!out_interval || walk->t_current >= t_max - RAY_EPSILON)
        return 2;
    const int rc = raycast_cell_aware_resume(
        sys, ray, t_max, true, false, NULL, NULL, NULL, 0.0, -1,
        0.0, -1, false, scratch, walk, 1, out_interval);
    if (rc != 0) return rc;
    return out_interval->t_exit >= t_max - RAY_EPSILON ? 0 : 1;
}

static int alea_ray_walk_next_first_visible(
    alea_system_t* sys, const alea_ray_t* ray, double t_max,
    double t_min, int material_filter, bool include_normal,
    alea_raycast_result_t* scratch, alea_ray_walk_t* walk,
    alea_ray_first_visible_result_t* visible) {
    alea_ray_selected_interval_t interval;
    const int rc = alea_ray_walk_next_selected(
        sys, ray, t_max, scratch, walk, &interval);
    if (rc < 0) return rc;
    if (rc == 2) return 0;
    if (ray_selected_interval_first_visible(
            sys, ray, &interval, t_min, material_filter, include_normal,
            visible))
        return 0;
    return rc == 0 ? 0 : 1;
}

#ifndef RAYCAST_FIRST_VISIBLE_PACKET_WIDTH
#define RAYCAST_FIRST_VISIBLE_PACKET_WIDTH 4
#endif

typedef struct {
    alea_ray_t ray;
    alea_ray_walk_t state;
    alea_raycast_result_t scratch;
    alea_ray_first_visible_result_t visible;
    size_t ray_index;
    uint8_t active;
} raycast_first_visible_lane_t;

static void raycast_first_visible_lane_init(
    raycast_first_visible_lane_t* lane,
    size_t ray_index,
    const double* origins_xyz,
    const double* directions_xyz) {
    memset(lane, 0, sizeof(*lane));
    lane->ray_index = ray_index;
    const double* o = &origins_xyz[ray_index * 3];
    const double* d = &directions_xyz[ray_index * 3];
    if (alea_ray_init(&lane->ray, o[0], o[1], o[2],
                      d[0], d[1], d[2]) != 0) {
        return;
    }
    alea_ray_walk_init(&lane->state);
    alea_raycast_result_init(&lane->scratch);
    lane->scratch.ray = lane->ray;
    lane->visible.cell_id = -1;
    lane->visible.surface_id = -1;
    lane->visible.primitive_id = UINT32_MAX;
    lane->active = 1;
}

static void raycast_first_visible_lane_store(
    const raycast_first_visible_lane_t* lane,
    alea_ray_first_visible_batch_result_t* result) {
    const size_t i = lane->ray_index;
    const alea_ray_first_visible_result_t* visible = &lane->visible;
    result->found[i] = visible->found ? 1 : 0;
    result->t[i] = visible->t;
    result->cell_ids[i] = visible->cell_id;
    result->material_ids[i] = visible->material_id;
    if (result->densities) result->densities[i] = visible->density;
    if (result->surface_ids) result->surface_ids[i] = visible->surface_id;
    if (result->primitive_ids) result->primitive_ids[i] = visible->primitive_id;
    if (result->resolution_flags)
        result->resolution_flags[i] = visible->resolution_flags;
    if (result->normals_xyz) {
        result->normals_xyz[i * 3] = visible->nx;
        result->normals_xyz[i * 3 + 1] = visible->ny;
        result->normals_xyz[i * 3 + 2] = visible->nz;
    }
}

int alea_raycast_hier_first_visible_batch_execute_nocache(
    alea_system_t* sys, const double* origins_xyz, const double* directions_xyz,
    size_t ray_count, const alea_ray_batch_query_t* query,
    alea_ray_first_visible_batch_result_t* result, int* statuses) {
    if (!sys || !origins_xyz || !directions_xyz || !query || !result ||
        !statuses) {
        return -1;
    }

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (size_t packet_start = 0; packet_start < ray_count;
         packet_start += RAYCAST_FIRST_VISIBLE_PACKET_WIDTH) {
        raycast_first_visible_lane_t lanes[RAYCAST_FIRST_VISIBLE_PACKET_WIDTH];
        const size_t remaining = ray_count - packet_start;
        const size_t packet_end = packet_start +
            (remaining < RAYCAST_FIRST_VISIBLE_PACKET_WIDTH
                 ? remaining : RAYCAST_FIRST_VISIBLE_PACKET_WIDTH);
        size_t active_count = 0;
        for (size_t i = packet_start; i < packet_end; i++) {
            raycast_first_visible_lane_init(&lanes[i - packet_start], i,
                                            origins_xyz, directions_xyz);
            if (!lanes[i - packet_start].active) {
                statuses[i] = -1;
                continue;
            }
            active_count++;
        }

        while (active_count > 0 && !alea_interrupted()) {
            for (size_t slot = 0; slot < packet_end - packet_start; slot++) {
                raycast_first_visible_lane_t* lane = &lanes[slot];
                if (!lane->active) continue;
                const size_t i = lane->ray_index;
                const double t_min = query->t_mins ? query->t_mins[i] : 0.0;
                const double t_max = query->t_maxs ? query->t_maxs[i] : 0.0;
                const double effective_t_max = t_max <= 0.0 ? DBL_MAX : t_max;
                const int rc = alea_ray_walk_next_first_visible(
                    sys, &lane->ray, effective_t_max, t_min,
                    query->material_filter,
                    (query->fields & ALEA_RAY_QUERY_FIELD_SURFACE_NORMAL) != 0,
                    &lane->scratch, &lane->state, &lane->visible);
                if (rc < 0) {
                    statuses[i] = -1;
                    lane->active = 0;
                    active_count--;
                    continue;
                }
                if (rc == 0) {
                    raycast_first_visible_lane_store(lane, result);
                    lane->active = 0;
                    active_count--;
                }
            }
        }
        for (size_t slot = 0; slot < packet_end - packet_start; slot++) {
            raycast_first_visible_lane_t* lane = &lanes[slot];
            if (lane->active) {
                statuses[lane->ray_index] = -1;
                lane->active = 0;
            }
            alea_raycast_result_free(&lane->scratch);
        }
    }
    return alea_interrupted() ? -1 : 0;
}

int alea_raycast_hier_any_hit_batch_execute_nocache(
    alea_system_t* sys, const double* origins_xyz, const double* directions_xyz,
    size_t ray_count, const alea_ray_batch_query_t* query,
    alea_ray_any_hit_batch_result_t* result, int* statuses) {
    if (!sys || !origins_xyz || !directions_xyz || !query || !result ||
        !statuses) {
        return -1;
    }

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (size_t packet_start = 0; packet_start < ray_count;
         packet_start += RAYCAST_FIRST_VISIBLE_PACKET_WIDTH) {
        raycast_first_visible_lane_t lanes[RAYCAST_FIRST_VISIBLE_PACKET_WIDTH];
        const size_t remaining = ray_count - packet_start;
        const size_t packet_end = packet_start +
            (remaining < RAYCAST_FIRST_VISIBLE_PACKET_WIDTH
                 ? remaining : RAYCAST_FIRST_VISIBLE_PACKET_WIDTH);
        size_t active_count = 0;
        for (size_t i = packet_start; i < packet_end; i++) {
            raycast_first_visible_lane_init(&lanes[i - packet_start], i,
                                            origins_xyz, directions_xyz);
            if (!lanes[i - packet_start].active) {
                statuses[i] = -1;
                continue;
            }
            active_count++;
        }

        while (active_count > 0 && !alea_interrupted()) {
            for (size_t slot = 0; slot < packet_end - packet_start; slot++) {
                raycast_first_visible_lane_t* lane = &lanes[slot];
                if (!lane->active) continue;
                const size_t i = lane->ray_index;
                const double t_min = query->t_mins ? query->t_mins[i] : 0.0;
                const double t_max = query->t_maxs ? query->t_maxs[i] : 0.0;
                const double effective_t_max = t_max <= 0.0 ? DBL_MAX : t_max;
                const int rc = alea_ray_walk_next_first_visible(
                    sys, &lane->ray, effective_t_max, t_min,
                    query->material_filter, false,
                    &lane->scratch, &lane->state, &lane->visible);
                if (rc < 0) {
                    statuses[i] = -1;
                    lane->active = 0;
                    active_count--;
                    continue;
                }
                if (rc == 0) {
                    result->hits[i] = lane->visible.found ? 1 : 0;
                    lane->active = 0;
                    active_count--;
                }
            }
        }
        for (size_t slot = 0; slot < packet_end - packet_start; slot++) {
            raycast_first_visible_lane_t* lane = &lanes[slot];
            if (lane->active) {
                statuses[lane->ray_index] = -1;
                lane->active = 0;
            }
            alea_raycast_result_free(&lane->scratch);
        }
    }
    return alea_interrupted() ? -1 : 0;
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
                                   NULL, NULL, NULL, 0, -1, 0, -1, false, result);
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
                                   NULL, NULL, NULL, 0, -1, 0, -1, false, result);
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
                                   NULL, NULL, NULL, 0, -1, 0, -1, false, result);
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
                                   NULL, NULL, NULL, 0, -1, 0, -1, false, result);
}

int alea_raycast_hier_first_visible_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_min, double t_max,
    int material_filter, int include_normal, alea_raycast_result_t* scratch,
    alea_ray_first_visible_result_t* out_visible) {
    if (!sys || !ray || !scratch || !out_visible || t_min < 0 ||
        (t_max > 0 && t_min > t_max))
        return -1;
    alea_raycast_result_clear(scratch);
    memset(out_visible, 0, sizeof(*out_visible));
    out_visible->cell_id = -1;
    out_visible->surface_id = -1;
    out_visible->primitive_id = UINT32_MAX;
    const double effective_t_max = t_max <= 0 ? DBL_MAX : t_max;
    scratch->ray = *ray;
    alea_ray_walk_t walk;
    alea_ray_walk_init(&walk);
    for (;;) {
        alea_ray_selected_interval_t interval;
        const int rc = alea_ray_walk_next_selected(
            sys, ray, effective_t_max, scratch, &walk, &interval);
        if (rc < 0) return -1;
        if (rc == 2) return 0;
        if (ray_selected_interval_first_visible(
                sys, ray, &interval, t_min, material_filter,
                include_normal != 0, out_visible))
            return 0;
        if (rc == 0) return 0;
    }
}

int alea_raycast_hier_first_cell_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_min, double t_max,
    int material_filter, alea_raycast_result_t* scratch,
    int* out_cell_id, double* out_t) {
    if (!sys || !ray || !scratch || !out_cell_id || t_min < 0 ||
        (t_max > 0 && t_min > t_max))
        return -1;
    *out_cell_id = -1;
    alea_raycast_result_clear(scratch);
    scratch->ray = *ray;
    const double effective_t_max = t_max <= 0 ? DBL_MAX : t_max;
    alea_ray_walk_t walk;
    alea_ray_walk_init(&walk);
    for (;;) {
        alea_ray_selected_interval_t interval;
        const int rc = alea_ray_walk_next_selected(
            sys, ray, effective_t_max, scratch, &walk, &interval);
        if (rc < 0) return -1;
        if (rc == 2) return 0;
        if (interval.cell_id >= 0 &&
            (material_filter < 0 || interval.material_id == material_filter) &&
            interval.t_exit > t_min + RAY_EPSILON) {
            *out_cell_id = interval.cell_id;
            if (out_t) *out_t = fmax(interval.t_enter, t_min);
            return 0;
        }
        if (rc == 0) return 0;
    }
}

int alea_raycast_hier_any_hit_nocache(
    alea_system_t* sys, const alea_ray_t* ray, double t_min, double t_max,
    int material_filter, alea_raycast_result_t* scratch, int* out_hit) {
    if (!sys || !ray || !scratch || !out_hit || t_min < 0 ||
        (t_max > 0 && t_min > t_max))
        return -1;
    *out_hit = 0;
    alea_raycast_result_clear(scratch);
    scratch->ray = *ray;
    const double effective_t_max = t_max <= 0 ? DBL_MAX : t_max;
    alea_ray_walk_t walk;
    alea_ray_walk_init(&walk);
    for (;;) {
        alea_ray_selected_interval_t interval;
        const int rc = alea_ray_walk_next_selected(
            sys, ray, effective_t_max, scratch, &walk, &interval);
        if (rc < 0) return -1;
        if (rc == 2) return 0;
        if (interval.cell_id >= 0 && interval.material_id != 0 &&
            (material_filter < 0 || interval.material_id == material_filter) &&
            interval.t_exit > t_min + RAY_EPSILON) {
            *out_hit = 1;
            return 0;
        }
        if (rc == 0) return 0;
    }
}

int alea_raycast_hier_segments_nocache(alea_system_t* sys,
                                       const alea_ray_t* ray,
                                       double t_max,
                                       alea_raycast_result_t* result) {
    if (!sys || !ray || !result) return -1;
    double effective_t_max = (t_max <= 0) ? DBL_MAX : t_max;
    result->ray = *ray;
    return raycast_cell_aware_impl(sys, ray, effective_t_max, true, false,
                                   NULL, NULL, NULL, 0, -1, 0, -1, false, result);
}
